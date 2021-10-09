#ifndef IKCONSTRAINT_JOINTVELOCITYCONSTRAINT_H
#define IKCONSTRAINT_JOINTVELOCITYCONSTRAINT_H

#include <ik_constraint/IKConstraint.h>
#include <cnoid/EigenUtil>

namespace IK{
  class JointVelocityConstraint : public IKConstraint
  {
  public:
    //jointのdqを上下限以下にする. radの次元で評価する
    //  dt: [s]
    //  maxError: エラーの頭打ち[rad]
    //  precision: 収束判定の閾値[rad]
    //  weight: コスト関数の重み. error * weight^2 * error.

    const cnoid::LinkPtr& joint() const { return joint_;}
    cnoid::LinkPtr& joint() { return joint_;}
    const double& dt() const { return dt_;}
    double& dt() { return dt_;}

    const double& maxError() const { return maxError_;}
    double& maxError() { return maxError_;}
    const double& precision() const { return precision_;}
    double& precision() { return precision_;}
    const double& weight() const { return weight_;}
    double& weight() { return weight_;}

    bool checkConvergence () override;
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) override;
    const Eigen::VectorXd& calc_minineq () override;
    const Eigen::VectorXd& calc_maxineq () override;

  private:
    cnoid::LinkPtr joint_ = nullptr;
    double dt_ = 0.1;
    double precision_ = 1e10;
    double maxError_ = 1e-2;
    double weight_ = 1.0;

    cnoid::LinkPtr jacobianineq_joint_ = nullptr; //前回jacobian_を計算した時のjoint
  };
}

#endif
