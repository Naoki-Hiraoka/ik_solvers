#ifndef IKCONSTRAINT_JOINTLIMITMINMAXTABLECONSTRAINT_H
#define IKCONSTRAINT_JOINTLIMITMINMAXTABLECONSTRAINT_H

#include <ik_constraint/JointLimitConstraint.h>
#include <joint_limit_table/JointLimitTable.h>

namespace IK{
  class JointLimitMinMaxTableConstraint : public JointLimitConstraint
  {
  public:
    //jointのqをq_upperとq_lowerの間かつmin-max-tableの間にさせる.
    //  maxError: エラーの頭打ち
    //  precision: 収束判定の閾値
    //  weight: コスト関数の重み. error * weight^2 * error.

    const cnoid::LinkPtr& joint() const { return joint_;}
    cnoid::LinkPtr& joint() { return joint_;}
    const std::vector<std::shared_ptr<joint_limit_table::JointLimitTable> >& jointLimitTables() const { return jointLimitTables_;}
    std::vector<std::shared_ptr<joint_limit_table::JointLimitTable> >& jointLimitTables() { return jointLimitTables_;}
    const double& maxError() const { return maxError_;}
    double& maxError() { return maxError_;}
    const double& precision() const { return precision_;}
    double& precision() { return precision_;}
    const double& weight() const { return weight_;}
    double& weight() { return weight_;}

    virtual bool checkConvergence () override;

  protected:
    std::vector<std::shared_ptr<joint_limit_table::JointLimitTable> > jointLimitTables_;
  };
}

#endif
