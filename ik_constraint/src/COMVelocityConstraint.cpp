#include <ik_constraint/COMVelocityConstraint.h>
#include <ik_constraint/Jacobian.h>
#include <cnoid/Jacobian>

namespace IK{
  bool COMVelocityConstraint::checkConvergence () {
    if(this->robot_) {
      Eigen::MatrixXd CMJ;
      cnoid::calcCMJacobian(this->robot_,nullptr,CMJ); // [joint root]の順
      cnoid::VectorX dq(this->robot_->numJoints()+6);
      for(int i=0;i<this->robot_->numJoints();i++) dq[i] = this->robot_->joint(i)->dq();
      dq.segment<3>(this->robot_->numJoints()) = this->robot_->rootLink()->v();
      dq.tail<3>() = this->robot_->rootLink()->w();
      cnoid::Vector3 dCM = this->eval_R_.transpose() * CMJ * dq;
      cnoid::Vector3 upper = (this->maxVel_ - dCM) * this->dt_;
      cnoid::Vector3 lower = (this->minVel_ - dCM) * this->dt_;

      if(this->minineq_.rows() != 3) this->minineq_ = Eigen::VectorXd(3);
      for(int i=0;i<3;i++) this->minineq_[i] = std::min(lower[i], this->maxError_[i]) * this->weight_[i];
      if(this->maxineq_.rows() != 3) this->maxineq_ = Eigen::VectorXd(3);
      for(int i=0;i<3;i++) this->maxineq_[i] = std::max(upper[i], -this->maxError_[i]) * this->weight_[i];

      bool converged = true;
      for(size_t i=0; i<3; i++){
        if(this->weight_[i]>0.0) {
          if(upper[i] < -this->precision_[i]) converged = false;
          if(lower[i] >  this->precision_[i]) converged = false;
        }
      }

      if(this->debuglevel_>=1){
        std::cerr << "COMVelocityConstraint" << std::endl;
        std::cerr << "dCM" << std::endl;
        std::cerr << dCM << std::endl;
        std::cerr << "upper" << std::endl;
        std::cerr << upper << std::endl;
        std::cerr << "lower" << std::endl;
        std::cerr << lower << std::endl;
      }

      return converged;
    }

    if(this->minineq_.rows() != 0) this->minineq_ = Eigen::VectorXd(0);
    if(this->maxineq_.rows() != 0) this->maxineq_ = Eigen::VectorXd(0);
    return true;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& COMVelocityConstraint::calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_joints_same(joints,this->jacobianineq_joints_)
       || this->robot_ != this->jacobianineq_robot_){
      this->jacobianineq_joints_ = joints;
      this->jacobianineq_robot_ = this->robot_;

      IK::calcCMJacobianShape(this->jacobianineq_joints_,
                              this->jacobianineq_robot_,
                              nullptr,
                              this->jacobianineq_full_,
                              this->jacobianineqColMap_);
    }

    IK::calcCMJacobianCoef(this->jacobianineq_joints_,
                           this->jacobianineq_robot_,
                           nullptr,
                           this->jacobianineqColMap_,
                           this->jacobianineq_full_);

    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R.insert(i,j) = this->eval_R_(i,j);

    this->jacobianineq_ = eval_R.transpose() * this->jacobianineq_full_;
    for(size_t i=0;i<3;i++) this->jacobianineq_.row(i) *= this->weight_[i];

    if(this->debuglevel_>=1){
      std::cerr << "COMVelocityConstraint" << std::endl;
      std::cerr << "jacobiaineq" << std::endl;
      std::cerr << this->jacobianineq_ << std::endl;
    }

    return this->jacobianineq_;
  }

  const Eigen::VectorXd& COMVelocityConstraint::calc_minineq () {
    if(this->debuglevel_>=1){
      std::cerr << "COMVelocityConstraint" << std::endl;
      std::cerr << "minineq" << std::endl;
      std::cerr << this->minineq_ << std::endl;
    }

    return this->minineq_;
  }

  const Eigen::VectorXd& COMVelocityConstraint::calc_maxineq () {
    if(this->debuglevel_>=1){
      std::cerr << "COMVelocityConstraint" << std::endl;
      std::cerr << "maxineq" << std::endl;
      std::cerr << this->maxineq_ << std::endl;
    }

    return this->maxineq_;
  }

}
