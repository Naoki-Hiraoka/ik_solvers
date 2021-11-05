#include <ik_constraint/JointLimitConstraint.h>
#include <iostream>

namespace IK{
  bool JointLimitConstraint::checkConvergence () {
    if(!this->joint_ || !(this->joint_->isRotationalJoint() || this->joint_->isPrismaticJoint())) {
      if(this->minineq_.rows() != 1) this->minineq_ = Eigen::VectorXd::Zero(1);
      if(this->maxineq_.rows() != 1) this->maxineq_ = Eigen::VectorXd::Zero(1);
      return true;
    }

    double lower = this->joint_->q_lower() - this->joint_->q();
    double upper = this->joint_->q_upper() - this->joint_->q();

    if(this->minineq_.rows() != 1) this->minineq_ = Eigen::VectorXd(1);
    this->minineq_[0] = std::min(lower, this->maxError_) * this->weight_;
    if(this->maxineq_.rows() != 1) this->maxineq_ = Eigen::VectorXd(1);
    this->maxineq_[0] = std::max(upper, -this->maxError_) * this->weight_;

    if(this->debuglevel_>=1){
      std::cerr << "JointLimitConstraint" << std::endl;
      std::cerr << "q" << std::endl;
      std::cerr << this->joint_->q() << std::endl;
      std::cerr << "q_upper" << std::endl;
      std::cerr << this->joint_->q_upper() << std::endl;
      std::cerr << "q_lower" << std::endl;
      std::cerr << this->joint_->q_lower() << std::endl;
    }

    return lower<this->precision_ && upper>-this->precision_;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& JointLimitConstraint::calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) {
    if(!this->is_joints_same(joints,this->jacobianineq_joints_) ||
       this->joint_ != this->jacobianineq_joint_){
      this->jacobianineq_joints_ = joints;
      this->jacobianineq_joint_ = this->joint_;
      this->jacobianineqColMap_.clear();
      int cols = 0;
      for(size_t i=0; i < this->jacobianineq_joints_.size(); i++){
        this->jacobianineqColMap_[this->jacobianineq_joints_[i]] = cols;
        cols += this->getJointDOF(this->jacobianineq_joints_[i]);
      }

      this->jacobianineq_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(1,cols);

      if(this->jacobianineqColMap_.find(this->jacobianineq_joint_) != this->jacobianineqColMap_.end()){
        if(this->jacobianineq_joint_->isRotationalJoint() || this->jacobianineq_joint_->isPrismaticJoint()){
          this->jacobianineq_.insert(0,this->jacobianineqColMap_[this->jacobianineq_joint_]) = 1;
        }
      }

    }

    if(this->jacobianineqColMap_.find(this->jacobianineq_joint_) != this->jacobianineqColMap_.end()){
      if(this->jacobianineq_joint_->isRotationalJoint() || this->jacobianineq_joint_->isPrismaticJoint()){
        this->jacobianineq_.coeffRef(0,this->jacobianineqColMap_[this->jacobianineq_joint_]) = this->weight_;
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "JointLimitConstraint" << std::endl;
      std::cerr << "jacobianineq" << std::endl;
      std::cerr << this->jacobianineq_ << std::endl;
    }
    return this->jacobianineq_;
  }

  const Eigen::VectorXd& JointLimitConstraint::calc_minineq () {
    if(this->debuglevel_>=1){
      std::cerr << "JointLimitConstraint" << std::endl;
      std::cerr << "minineq" << std::endl;
      std::cerr << this->minineq_ << std::endl;
    }

    return this->minineq_;
  }

  const Eigen::VectorXd& JointLimitConstraint::calc_maxineq () {
    if(this->debuglevel_>=1){
      std::cerr << "JointLimitConstraint" << std::endl;
      std::cerr << "maxineq" << std::endl;
      std::cerr << this->maxineq_ << std::endl;
    }

    return this->maxineq_;
  }

}
