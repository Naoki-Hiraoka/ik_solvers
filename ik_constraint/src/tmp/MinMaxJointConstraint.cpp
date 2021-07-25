#include "MinMaxJointConstraint.h"

namespace IK{
  MinMaxJointConstraint::MinMaxJointConstraint(const cnoid::Link* _joint, std::shared_ptr<RobotConfig::JointLimitTable> _minmaxtable):
    joint(_joint),
    minmaxtable(_minmaxtable),
    maxvel(0.1),
    llimit(-1e30),
    ulimit(1e30)
  {
    return;
  }

  const Eigen::VectorXd& MinMaxJointConstraint::calc_minineq () {
    if(this->minineq.rows()!=1) this->minineq = Eigen::VectorXd(1);
    double limit = std::max(this->joint->q_lower(),this->llimit);
    if (this->minmaxtable) limit = std::max(limit,this->minmaxtable->getLlimit());
    this->minineq[0] = std::max(limit - this->joint->q(), -maxvel);
    return this->minineq;
  }

  const Eigen::VectorXd& MinMaxJointConstraint::calc_maxineq () {
    if(this->maxineq.rows()!=1) this->maxineq = Eigen::VectorXd(1);
    double limit = std::min(this->joint->q_upper(),this->ulimit);
    if (this->minmaxtable) limit = std::min(ulimit,this->minmaxtable->getUlimit());
    this->maxineq[0] = std::min(limit - this->joint->q(), maxvel);
    return this->maxineq;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& MinMaxJointConstraint::calc_jacobianineq (const std::vector<cnoid::Body*>& bodies) {
    if(!this->is_bodies_same(bodies,this->jacobianineq_bodies) || this->jacobianineq.rows()!=1){
      this->jacobianineq_bodies = bodies;

      int dim = 0;
      for(size_t i=0; i < bodies.size(); i++){
        dim += 6 + bodies[i]->numJoints();
      }

      this->jacobianineq = Eigen::SparseMatrix<double,Eigen::RowMajor>(1,dim);

      int idx = 0;
      for(size_t b=0;b<bodies.size();b++){
        if(bodies[b] == this->joint->body()){
          this->jacobianineq.insert(0,idx+6+this->joint->jointId()) = 1;
          break;
        }
        idx += 6 + bodies[b]->numJoints();
      }
    }

    return this->jacobianineq;
  }
}
