#ifndef IKCONSTRAINT_JOINTANGLECONSTRAINT_H
#define IKCONSTRAINT_JOINTANGLECONSTRAINT_H

#include <ik_constraint/IKConstraint.h>
#include <cnoid/EigenUtil>

namespace IK{
  class JointAngleConstraint : public IKConstraint
  {
  public:
    //jointのqとtargetqを一致させる.
    //  maxError: エラーの頭打ち
    //  precision: 収束判定の閾値
    //  weight: コスト関数の重み. error * weight^2 * error.

    const cnoid::LinkPtr& joint() const { return joint_;}
    cnoid::LinkPtr& joint() { return joint_;}
    const double& targetq() const { return targetq_;}
    double& targetq() { return targetq_;}
    const double& maxError() const { return maxError_;}
    double& maxError() { return maxError_;}
    const double& precision() const { return precision_;}
    double& precision() { return precision_;}
    const double& weight() const { return weight_;}
    double& weight() { return weight_;}

    bool checkConvergence () override;
    const Eigen::VectorXd& calc_error () override;
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) override;

  private:
    cnoid::LinkPtr joint_ = nullptr;
    double targetq_ = 0.0;
    double precision_ = 1e10;
    double maxError_ = 1e-2;
    double weight_ = 1.0;

    cnoid::LinkPtr jacobian_joint_ = nullptr; //前回jacobian_を計算した時のjoint
  };
}

#endif
