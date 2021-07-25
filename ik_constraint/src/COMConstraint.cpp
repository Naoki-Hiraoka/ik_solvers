#include <ik_constraint/COMConstraint.h>
#include <cnoid/Jacobian>

namespace IK{
  bool COMConstraint::checkConvergence () {
    cnoid::Vector3 error = this->targetPos_ - this->robot_->centerOfMass();

    // 収束判定と、ついでにcalc_errorの返り値の計算
    if(this->error_.rows()!=(this->weight_.array() > 0.0).count()) this->error_ = Eigen::VectorXd((this->weight_.array() > 0.0).count());
    bool converged = true;
    int idx=0;
    for(size_t i=0; i<3; i++){
      if(this->weight_[i]>0.0) {
        if(std::fabs(error[i]) > this->precision_[i]) converged = false;
        this->error_[idx] = std::min(std::max(error[i],-this->maxError_[i]),this->maxError_[i]) * this->weight_[i];
        idx++;
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "COM" << std::endl;
      std::cerr << this->robot_->centerOfMass() << std::endl;
      std::cerr << "targetPos" << std::endl;
      std::cerr << this->targetPos_ << std::endl;
    }

    return converged;
  }

  const Eigen::VectorXd& COMConstraint::calc_error (){
    return this->error_;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& COMConstraint::calc_jacobian (const std::vector<cnoid::BodyPtr>& bodies) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_bodies_same(bodies,this->jacobian_bodies_)
       || this->robot_ != this->jacobian_robot_
       || (this->weight_.array() > 0.0).count() != this->jacobian_.rows()){
      this->jacobian_bodies_ = bodies;
      this->jacobian_robot_ = this->robot_;

      int rows = (this->weight_.array() > 0.0).count();
      int cols = 0;
      for(size_t i=0; i < bodies.size(); i++){
        cols += 6 + bodies[i]->numJoints();
      }
      this->jacobian_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(rows,cols);

      std::vector<Eigen::Triplet<double> > tripletList;
      tripletList.reserve(100);//適当
      if(this->robot_){
        int idx = 0;
        for(size_t b=0;b<bodies.size();b++){
          if(bodies[b] == this->robot_){
            for(size_t row=0;row<rows;row++){
              for(size_t col=0;col<6+bodies[b]->numJoints();col++){
                tripletList.push_back(Eigen::Triplet<double>(row,idx+col,1));
              }
            }
            break;
          }
          idx += 6 + bodies[b]->numJoints();
        }
      }

      this->jacobian_.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    // calc CM jacobian
    Eigen::MatrixXd _CMJ;
    cnoid::calcCMJacobian(this->robot_,nullptr,_CMJ); // [joint root]の順
    Eigen::MatrixXd CMJ(_CMJ.rows(),_CMJ.cols()); // [root joint]の順
    CMJ.leftCols(6) = _CMJ.rightCols(6);
    CMJ.rightCols(CMJ.cols()-6) = _CMJ.leftCols(_CMJ.cols()-6);

    int col_idx = 0;
    for(size_t b=0;b<bodies.size();b++){
      if(bodies[b] == this->robot_){
        int row_idx=0;
        for(size_t i=0;i<3;i++){
          if(this->weight_[i]>0.0){
            for(size_t j=0;j<CMJ.cols();j++){
              this->jacobian_.coeffRef(row_idx,col_idx+j) = CMJ(i,j);
            }
            row_idx++;
          }
        }
        break;
      }
      col_idx += 6 + bodies[b]->numJoints();
    }

    return this->jacobian_;
  }

  std::vector<cnoid::SgNodePtr>& COMConstraint::getDrawOnObjects() {

    if(!this->COMMarker_){
      this->COMMarker_ = new cnoid::CrossMarker(0.2/*size*/,cnoid::Vector3f(0.5,1.0,0.0)/*color*/,10/*width*/);
    }
    if(!this->targetMarker_){
      this->targetMarker_ = new cnoid::CrossMarker(0.2/*size*/,cnoid::Vector3f(0.5,1.0,0.0)/*color*/,10/*width*/);
    }

    // update position
    this->COMMarker_->setTranslation(this->robot_->centerOfMass().cast<Eigen::Vector3f::Scalar>());
    this->targetMarker_->setTranslation(this->targetPos_.cast<Eigen::Vector3f::Scalar>());

    this->drawOnObjects_ = std::vector<cnoid::SgNodePtr>{this->COMMarker_,this->targetMarker_};

    return this->drawOnObjects_;
  }

}
