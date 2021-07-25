#include <ik_constraint/JointAngleConstraint.h>
#include <iostream>

namespace IK{
  bool JointAngleConstraint::checkConvergence () {
    if(!this->joint_) return true;

    double error = this->joint_->q() - targetq_;

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

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& JointAngleConstraint::calc_jacobian (const std::vector<cnoid::BodyPtr>& bodies) {
    if(!this->is_bodies_same(bodies,this->jacobian_bodies_) ||
       this->joint_ != this->jacobian_joint_){
      this->jacobian_bodies_ = bodies;
      this->jacobian_joint_ = this->joint_;

      int cols = 0;
      for(size_t i=0; i < bodies.size(); i++){
        cols += 6 + bodies[i]->numJoints();
      }

      this->jacobian_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(1,cols);

      int idx = 0;
      for(size_t b=0;b<bodies.size();b++){
        if(bodies[b] == this->joint_->body()){
          this->jacobian_.insert(0,idx+6+this->joint_->jointId()) = 1;
          this->jacobian_col_ = idx+6+this->joint_->jointId();
          break;
        }
        idx += 6 + bodies[b]->numJoints();
      }
    }

    this->jacobian_.coeffRef(0,this->jacobian_col_) = this->weight_;

    if(this->debuglevel_>=1){
      std::cerr << "JointAngleConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
    }
    return this->jacobian_;
  }
}
