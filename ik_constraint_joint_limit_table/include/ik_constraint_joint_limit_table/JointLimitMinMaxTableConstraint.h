#ifndef IKCONSTRAINT_JOINTLIMITMINMAXTABLECONSTRAINT_H
#define IKCONSTRAINT_JOINTLIMITMINMAXTABLECONSTRAINT_H

#include <ik_constraint/JointLimitConstraint.h>
#include <joint_limit_table/JointLimitTable.h>

namespace ik_constraint_joint_limit_table{
  class JointLimitMinMaxTableConstraint : public IK::JointLimitConstraint {
  public:
    //jointのqをq_upperとq_lowerの間かつmin-max-tableの間にさせる.

    const std::vector<std::shared_ptr<joint_limit_table::JointLimitTable> >& jointLimitTables() const { return jointLimitTables_;}
    std::vector<std::shared_ptr<joint_limit_table::JointLimitTable> >& jointLimitTables() { return jointLimitTables_;}

    virtual bool checkConvergence () override;

  protected:
    std::vector<std::shared_ptr<joint_limit_table::JointLimitTable> > jointLimitTables_;
  };
}

#endif
