#include <ik_constraint/COMConstraint.h>
#include <cnoid/Jacobian>

namespace IK{
  bool COMConstraint::checkConvergence () {
    // A - B
    cnoid::Vector3 A_p = this->A_robot_ ? this->A_robot_->centerOfMass() + this->A_localp() : this->A_localp();
    cnoid::Vector3 B_p = this->B_robot_ ? this->B_robot_->centerOfMass() + this->B_localp() : this->B_localp();
    cnoid::Vector3 error = A_p - B_p;
    error = (this->eval_R_.transpose() * error).eval();

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
      std::cerr << "A COM" << std::endl;
      std::cerr << A_p.transpose() << std::endl;
      std::cerr << "B COM" << std::endl;
      std::cerr << B_p.transpose() << std::endl;
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
       || this->A_robot_ != this->jacobian_A_robot_
       || this->B_robot_ != this->jacobian_B_robot_
       || (this->weight_.array() > 0.0).count() != this->jacobian_.rows()){
      this->jacobian_joints_ = joints;
      this->jacobian_A_robot_ = this->A_robot_;
      this->jacobian_B_robot_ = this->B_robot_;

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

      if(this->jacobian_A_robot_ != this->jacobian_B_robot_){
        for(int i=0;i<2;i++){
          cnoid::BodyPtr robot = (i==0) ? this->jacobian_A_robot_ : this->jacobian_B_robot_;
          if(!robot) continue;
          if(this->jacobianColMap_.find(robot->rootLink()) != this->jacobianColMap_.end()){
            int idx = this->jacobianColMap_[robot->rootLink()];
            for(size_t row=0;row<rows;row++){
              for(size_t j=0;j<this->getJointDOF(robot->rootLink());j++){
                tripletList.push_back(Eigen::Triplet<double>(row,idx+j,1));
              }
            }
          }
          for(size_t j=0;j<robot->numJoints();j++){
            if(this->jacobianColMap_.find(robot->joint(j)) != this->jacobianColMap_.end()){
              int idx = this->jacobianColMap_[robot->joint(j)];
              for(size_t row=0;row<rows;row++){
                for(size_t k=0;k<this->getJointDOF(robot->joint(j));k++){
                  tripletList.push_back(Eigen::Triplet<double>(row,idx+k,1));
                }
              }
            }
          }
        }
      }

      this->jacobian_.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    if(this->jacobian_A_robot_ != this->jacobian_B_robot_){
      for(size_t i=0;i<2; i++){
        cnoid::BodyPtr robot = (i==0) ? this->jacobian_A_robot_ : this->jacobian_B_robot_;
        if(!robot) continue;
        int sign = (i==0) ? 1 : -1;

        Eigen::MatrixXd CMJ;
        cnoid::calcCMJacobian(robot,nullptr,CMJ); // [joint root]の順
        CMJ = (this->eval_R_.transpose() * CMJ).eval();

        if(this->jacobianColMap_.find(robot->rootLink()) != this->jacobianColMap_.end()){
          int col_idx = this->jacobianColMap_[robot->rootLink()];
          int row_idx=0;
          for(size_t j=0;j<3;j++){
            if(this->weight_[j]>0.0){
              for(size_t k=0;k<this->getJointDOF(robot->rootLink());k++){
                this->jacobian_.coeffRef(row_idx,col_idx+k) = this->weight_[j] * sign * CMJ(j,robot->numJoints()+k);
              }
              row_idx++;
            }
          }
        }
        for(size_t j=0;j<robot->numJoints();j++){
          if(this->jacobianColMap_.find(robot->joint(j)) != this->jacobianColMap_.end()){
            int col_idx = this->jacobianColMap_[robot->joint(j)];
            int row_idx=0;
            for(size_t k=0;k<3;k++){
              if(this->weight_[k]>0.0){
                for(size_t d=0;d<this->getJointDOF(robot->joint(j));d++){
                  this->jacobian_.coeffRef(row_idx,col_idx+d) = this->weight_[k] * sign * CMJ(k,j+d);
                }
                row_idx++;
              }
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

}
