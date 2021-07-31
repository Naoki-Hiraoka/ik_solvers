#include <ik_constraint/COMConstraint.h>
#include <cnoid/Jacobian>

namespace IK{
  bool COMConstraint::checkConvergence () {
    cnoid::Vector3 error = this->robot_->centerOfMass() - this->targetPos_;

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
      std::cerr << this->robot_->centerOfMass().transpose() << std::endl;
      std::cerr << "targetPos" << std::endl;
      std::cerr << this->targetPos_.transpose() << std::endl;
    }

    return converged;
  }

  const Eigen::VectorXd& COMConstraint::calc_error (){
    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "error" << std::endl;
      std::cerr << this->error_.transpose() << std::endl;
    }

    return this->error_;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& COMConstraint::calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_joints_same(joints,this->jacobian_joints_)
       || this->robot_ != this->jacobian_robot_
       || (this->weight_.array() > 0.0).count() != this->jacobian_.rows()){
      this->jacobian_joints_ = joints;
      this->jacobian_robot_ = this->robot_;

      int rows = (this->weight_.array() > 0.0).count();
      this->jacobianColMap_.clear();
      int cols = 0;
      for(size_t i=0;i<this->jacobian_joints_.size();i++){
        this->jacobianColMap_[this->jacobian_joints_[i]] = cols;
        cols += this->getJointDOF(this->jacobian_joints_[i]);
      }
      this->jacobian_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(rows,cols);

      std::vector<Eigen::Triplet<double> > tripletList;
      tripletList.reserve(100);//適当
      if(this->jacobian_robot_){
        if(this->jacobianColMap_.find(this->jacobian_robot_->rootLink()) != this->jacobianColMap_.end()){
          int idx = this->jacobianColMap_[this->jacobian_robot_->rootLink()];
          for(size_t row=0;row<rows;row++){
            for(size_t i=0;i<this->getJointDOF(this->jacobian_robot_->rootLink());i++){
              tripletList.push_back(Eigen::Triplet<double>(row,idx+i,1));
            }
          }
        }
        for(size_t j=0;j<this->jacobian_robot_->numJoints();j++){
          if(this->jacobianColMap_.find(this->jacobian_robot_->joint(j)) != this->jacobianColMap_.end()){
            int idx = this->jacobianColMap_[this->jacobian_robot_->joint(j)];
            for(size_t row=0;row<rows;row++){
              for(size_t i=0;i<this->getJointDOF(this->jacobian_robot_->joint(j));i++){
                tripletList.push_back(Eigen::Triplet<double>(row,idx+i,1));
              }
            }
          }
        }

      }

      this->jacobian_.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    if(this->jacobian_robot_){
      Eigen::MatrixXd CMJ;
      cnoid::calcCMJacobian(this->jacobian_robot_,nullptr,CMJ); // [joint root]の順

      if(this->jacobianColMap_.find(this->jacobian_robot_->rootLink()) != this->jacobianColMap_.end()){
        int col_idx = this->jacobianColMap_[this->jacobian_robot_->rootLink()];
        int row_idx=0;
        for(size_t i=0;i<3;i++){
          if(this->weight_[i]>0.0){
            for(size_t j=0;j<this->getJointDOF(this->jacobian_robot_->rootLink());j++){
              this->jacobian_.coeffRef(row_idx,col_idx+j) = this->weight_[i] * CMJ(i,this->jacobian_robot_->numJoints()+j);
            }
            row_idx++;
          }
        }
      }
      for(size_t j=0;j<this->jacobian_robot_->numJoints();j++){
        if(this->jacobianColMap_.find(this->jacobian_robot_->joint(j)) != this->jacobianColMap_.end()){
          int col_idx = this->jacobianColMap_[this->jacobian_robot_->joint(j)];
          int row_idx=0;
          for(size_t i=0;i<3;i++){
            if(this->weight_[i]>0.0){
              for(size_t d=0;d<this->getJointDOF(this->jacobian_robot_->joint(j));d++){
                this->jacobian_.coeffRef(row_idx,col_idx+d) = this->weight_[i] * CMJ(i,j+d);
              }
              row_idx++;
            }
          }
        }
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
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
