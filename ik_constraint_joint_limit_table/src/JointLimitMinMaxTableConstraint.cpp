#include <ik_constraint_joint_limit_table/JointLimitMinMaxTableConstraint.h>
#include <iostream>

namespace ik_constraint_joint_limit_table{
  bool JointLimitMinMaxTableConstraint::checkConvergence () {
    if(this->minineq_.rows() != 1) this->minineq_ = Eigen::VectorXd::Zero(1);
      if(this->maxineq_.rows() != 1) this->maxineq_ = Eigen::VectorXd::Zero(1);

    if(!this->joint_ || !(this->joint_->isRotationalJoint() || this->joint_->isPrismaticJoint())) {
      return true;
    }

    double lower = this->joint_->q_lower();
    double upper = this->joint_->q_upper();
    for(size_t i=0;i<this->jointLimitTables_.size();i++){
      if(this->jointLimitTables_[i]->getSelfJoint() == this->joint_){
        lower = std::max(lower, this->jointLimitTables_[i]->getLlimit());
        upper = std::min(upper, this->jointLimitTables_[i]->getUlimit());
      }
    }
    lower -= this->joint_->q();
    upper -= this->joint_->q();

    this->minineq_[0] = std::min(this->weight_ * lower, this->maxError_);
    this->maxineq_[0] = std::max(this->weight_ * upper, -this->maxError_);

    if(this->debuglevel_>=1){
      std::cerr << "JointLimitMinMaxTableConstraint " << this->joint_->name() << std::endl;
      std::cerr << "q: " << this->joint_->q() << std::endl;
      std::cerr << "upper: " << upper << std::endl;
      std::cerr << "lower: " << lower << std::endl;
      std::cerr << "tables: " << this->jointLimitTables_.size()<<std::endl;
    }

    return lower<this->precision_ && upper>-this->precision_;
  }
}
