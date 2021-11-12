#include <iostream>

#include <ik_constraint/CollisionConstraint.h>
#include <ik_constraint/Jacobian.h>
#include <cnoid/EigenUtil>

namespace IK{

  CollisionConstraint::CollisionConstraint()
  {
  }

  bool CollisionConstraint::checkConvergence () {
    // 収束判定と、ついでにcalc_minineq/maxineqの返り値の計算

    if(this->minineq_.rows()!=1) this->minineq_ = Eigen::VectorXd::Zero(1);
    if(this->maxineq_.rows()!=1) this->maxineq_ = Eigen::VectorXd::Zero(1);

    double distance;
    bool ret;
    if(!this->computeDistance(this->A_link_, this->B_link_,
                              distance, this->currentDirection_, this->A_currentLocalp_,this->B_currentLocalp_)){
      this->minineq_[0] = -1e10;
      this->maxineq_[0] = 1e10;
      ret = true;
    }else{
      this->minineq_[0] = std::min((this->tolerance_ - distance) / this->velocityDamper_, this->maxError_) * this->weight_;
      this->maxineq_[0] = 1e10;

      ret =  distance - this->tolerance_ > - this->precision_;
    }

    if(this->debuglevel_>=1){
      std::cerr << "CollisionConstraint " << this->A_link_->name() << " - " << this->B_link_->name() << std::endl;
      std::cerr << "distance: " << distance << std::endl;
    }

    return ret;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& CollisionConstraint::calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_joints_same(joints,this->jacobianineq_joints_)
       || this->A_link_ != this->jacobianineq_A_link_
       || this->B_link_ != this->jacobianineq_B_link_){
      this->jacobianineq_joints_ = joints;
      this->jacobianineq_A_link_ = this->A_link_;
      this->jacobianineq_B_link_ = this->B_link_;

      IK::calc6DofJacobianShape(this->jacobianineq_joints_,//input
                                this->jacobianineq_A_link_,//input
                                this->jacobianineq_B_link_,//input
                                this->jacobianineq_full_,
                                this->jacobianineqColMap_,
                                this->path_A_joints_,
                                this->path_B_joints_,
                                this->path_BA_joints_,
                                this->path_BA_joints_numUpwardConnections_
                                );
    }

    cnoid::Position A_localpos = cnoid::Position::Identity();
    A_localpos.translation() = this->A_currentLocalp_;
    cnoid::Position B_localpos = cnoid::Position::Identity();
    B_localpos.translation() = this->B_currentLocalp_;
    IK::calc6DofJacobianCoef(this->jacobianineq_joints_,//input
                             this->jacobianineq_A_link_,//input
                             A_localpos,//input
                             this->jacobianineq_B_link_,//input
                             B_localpos,//input
                             this->jacobianineqColMap_,//input
                             this->path_A_joints_,//input
                             this->path_B_joints_,//input
                             this->path_BA_joints_,//input
                             this->path_BA_joints_numUpwardConnections_,//input
                             this->jacobianineq_full_
                             );

    Eigen::SparseMatrix<double,Eigen::RowMajor> dir(3,1);
    for(int i=0;i<3;i++) dir.insert(i,0) = this->currentDirection_[i];
    this->jacobianineq_ = dir.transpose() * this->jacobianineq_full_.topRows<3>() * this->weight_;

    if(this->debuglevel_>=1){
      std::cerr << "CollisionConstraint " << (this->A_link_ ? this->A_link_->name() : "world") << " - " << (this->B_link_ ? this->B_link_->name() : "world") << std::endl;
      std::cerr << "direction" << std::endl;
      std::cerr << dir << std::endl;
      std::cerr << "jacobianineq" << std::endl;
      std::cerr << this->jacobianineq_ << std::endl;
    }
    return this->jacobianineq_;
  }


  const Eigen::VectorXd& CollisionConstraint::calc_minineq () {
    if(this->debuglevel_>=1){
      std::cerr << "CollisionConstraint" << std::endl;
      std::cerr << "minineq" << std::endl;
      std::cerr << this->minineq_ << std::endl;
    }

    return this->minineq_;
  }

  const Eigen::VectorXd& CollisionConstraint::calc_maxineq () {
    if(this->debuglevel_>=1){
      std::cerr << "CollisionConstraint" << std::endl;
      std::cerr << "maxineq" << std::endl;
      std::cerr << this->maxineq_ << std::endl;
    }

    return this->maxineq_;
  }

}
