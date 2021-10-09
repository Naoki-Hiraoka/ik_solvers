#include <ik_constraint/JointAngleConstraint.h>
#include <iostream>

namespace IK{
  bool JointAngleConstraint::checkConvergence () {
    if(!this->joint_) return true;

    double error = std::min(std::max(this->joint_->q() - targetq_,-this->maxError_), this->maxError_);

    if(this->error_.rows() != 1) this->error_ = Eigen::VectorXd(1);
    this->error_[0] = this->weight_ * error;

    if(this->debuglevel_>=1){
      std::cerr << "JointAngleConstraint" << std::endl;
      std::cerr << "q" << std::endl;
      std::cerr << this->joint_->q() << std::endl;
      std::cerr << "targetq" << std::endl;
      std::cerr << this->targetq_ << std::endl;
    }

    return std::fabs(error) < this->precision_;
  }

  const Eigen::VectorXd& JointAngleConstraint::calc_error () {
    if(this->debuglevel_>=1){
      std::cerr << "JointAngleConstraint" << std::endl;
      std::cerr << "error" << std::endl;
      std::cerr << this->error_ << std::endl;
    }

    return this->error_;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& JointAngleConstraint::calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) {
    if(!this->is_joints_same(joints,this->jacobian_joints_) ||
       this->joint_ != this->jacobian_joint_){
      this->jacobian_joints_ = joints;
      this->jacobian_joint_ = this->joint_;
      this->jacobianColMap_.clear();
      int cols = 0;
      for(size_t i=0; i < this->jacobian_joints_.size(); i++){
        this->jacobianColMap_[this->jacobian_joints_[i]] = cols;
        cols += this->getJointDOF(this->jacobian_joints_[i]);
      }

      this->jacobian_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(1,cols);

      if(this->jacobianColMap_.find(this->jacobian_joint_) != this->jacobianColMap_.end()){
        if(this->jacobian_joint_->isRotationalJoint() || this->jacobian_joint_->isPrismaticJoint()){
          this->jacobian_.insert(0,this->jacobianColMap_[this->jacobian_joint_]) = 1;
        }
      }

    }

    if(this->jacobianColMap_.find(this->jacobian_joint_) != this->jacobianColMap_.end()){
      if(this->jacobian_joint_->isRotationalJoint() || this->jacobian_joint_->isPrismaticJoint()){
        this->jacobian_.coeffRef(0,this->jacobianColMap_[this->jacobian_joint_]) = this->weight_;
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "JointAngleConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
    }
    return this->jacobian_;
  }
}
