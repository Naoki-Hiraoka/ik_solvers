#ifndef MINMAXJOINTCONSTRAINT_H
#define MINMAXJOINTCONSTRAINT_H

#include "IKConstraint.h"
#include "../../RobotConfig/RobotConfig.h"

namespace IK{
  class MinMaxJointConstraint : public IKConstraint
  {
  public:
    MinMaxJointConstraint(const cnoid::Link* joint, std::shared_ptr<RobotConfig::JointLimitTable> minmaxtable);

    const Eigen::VectorXd& calc_minineq () override;
    const Eigen::VectorXd& calc_maxineq () override;
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobianineq (const std::vector<cnoid::Body*>& bodies) override;

    void set_maxvel (double _maxvel) {maxvel = _maxvel;}
    void set_llimit (double _llimit) {llimit = _llimit;}
    void set_ulimit (double _ulimit) {ulimit = _ulimit;}
  private:
    const cnoid::Link* joint;
    std::shared_ptr<RobotConfig::JointLimitTable> minmaxtable;
    double maxvel;
    double llimit;
    double ulimit;
  };
}

#endif
