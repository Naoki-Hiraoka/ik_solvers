#include <ik_constraint/JointVelocityConstraint.h>
#include <iostream>

namespace IK{
  bool JointVelocityConstraint::checkConvergence () {
    if(this->joint_) {
      if (this->joint_->isRotationalJoint() || this->joint_->isPrismaticJoint()) {

        double lower = (this->joint_->dq_lower() - this->joint_->dq()) * dt_;
        double upper = (this->joint_->dq_upper() - this->joint_->dq()) * dt_;

        if(this->minineq_.rows() != 1) this->minineq_ = Eigen::VectorXd(1);
        this->minineq_[0] = std::min(lower, this->maxError_) * this->weight_;
        if(this->maxineq_.rows() != 1) this->maxineq_ = Eigen::VectorXd(1);
        this->maxineq_[0] = std::max(upper, -this->maxError_) * this->weight_;

        if(this->debuglevel_>=1){
          std::cerr << "JointVelocityConstraint" << std::endl;
          std::cerr << "dq" << std::endl;
          std::cerr << this->joint_->dq() << std::endl;
          std::cerr << "dq_upper" << std::endl;
          std::cerr << this->joint_->dq_upper() << std::endl;
          std::cerr << "dq_lower" << std::endl;
          std::cerr << this->joint_->dq_lower() << std::endl;
        }

        return lower<this->precision_ && upper>-this->precision_;
      }else if (this->joint_->isFreeJoint()){
        cnoid::Vector6 lower, upper;
        lower.head<3>() = (cnoid::Vector3::Ones() * this->joint_->dq_lower() - this->joint_->v()) * dt_;
        lower.tail<3>() = (cnoid::Vector3::Ones() * this->joint_->dq_lower() - this->joint_->w()) * dt_;
        upper.head<3>() = (cnoid::Vector3::Ones() * this->joint_->dq_upper() - this->joint_->v()) * dt_;
        upper.tail<3>() = (cnoid::Vector3::Ones() * this->joint_->dq_upper() - this->joint_->w()) * dt_;

        if(this->minineq_.rows() != 6) this->minineq_ = Eigen::VectorXd(6);
        for(size_t i=0;i<6;i++) this->minineq_[i] = std::min(lower[i], this->maxError_) * this->weight_;
        if(this->maxineq_.rows() != 6) this->maxineq_ = Eigen::VectorXd(6);
        for(size_t i=0;i<6;i++) this->maxineq_[i] = std::max(upper[i], -this->maxError_) * this->weight_;

        if(this->debuglevel_>=1){
          std::cerr << "JointVelocityConstraint" << std::endl;
          std::cerr << "v" << std::endl;
          std::cerr << this->joint_->v() << std::endl;
          std::cerr << "w" << std::endl;
          std::cerr << this->joint_->w() << std::endl;
          std::cerr << "dq_upper" << std::endl;
          std::cerr << this->joint_->dq_upper() << std::endl;
          std::cerr << "dq_lower" << std::endl;
          std::cerr << this->joint_->dq_lower() << std::endl;
        }

        return (lower.array()<this->precision_).count()==6 && (upper.array()>-this->precision_).count()==6;
      }
    }
    if(this->error_.rows() != 0) this->error_ = Eigen::VectorXd::Zero(0);
    return true;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& JointVelocityConstraint::calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) {
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

      int rows;
      if (this->joint_->isRotationalJoint() || this->joint_->isPrismaticJoint()) rows=1;
      else if (this->joint_->isFreeJoint()) rows = 6;
      else rows = 0;

      this->jacobianineq_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(rows,cols);

      if(this->jacobianineqColMap_.find(this->jacobianineq_joint_) != this->jacobianineqColMap_.end()){
        if(this->jacobianineq_joint_->isRotationalJoint() || this->jacobianineq_joint_->isPrismaticJoint()){
          for(size_t i=0;i<rows;i++){
            this->jacobianineq_.insert(rows,this->jacobianineqColMap_[this->jacobianineq_joint_]+i) = 1;
          }
        }
      }

    }

    if(this->jacobianineqColMap_.find(this->jacobianineq_joint_) != this->jacobianineqColMap_.end()){
      int rows;
      if (this->jacobianineq_joint_->isRotationalJoint() || this->jacobianineq_joint_->isPrismaticJoint()) rows=1;
      else if (this->jacobianineq_joint_->isFreeJoint()) rows = 6;
      else rows = 0;

      for(size_t i=0;i<rows;i++){
        this->jacobianineq_.coeffRef(i,this->jacobianineqColMap_[this->jacobianineq_joint_]+i) = this->weight_;
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "JointVelocityConstraint" << std::endl;
      std::cerr << "jacobianineq" << std::endl;
      std::cerr << this->jacobianineq_ << std::endl;
    }
    return this->jacobianineq_;
  }

  const Eigen::VectorXd& JointVelocityConstraint::calc_minineq () {
    if(this->debuglevel_>=1){
      std::cerr << "JointVelocityConstraint" << std::endl;
      std::cerr << "minineq" << std::endl;
      std::cerr << this->minineq_ << std::endl;
    }

    return this->minineq_;
  }

  const Eigen::VectorXd& JointVelocityConstraint::calc_maxineq () {
    if(this->debuglevel_>=1){
      std::cerr << "JointVelocityConstraint" << std::endl;
      std::cerr << "maxineq" << std::endl;
      std::cerr << this->maxineq_ << std::endl;
    }

    return this->maxineq_;
  }

}
