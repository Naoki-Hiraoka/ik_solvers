#ifndef COMCONSTRAINT_H
#define COMCONSTRAINT_H

#include <ik_constraint/IKConstraint.h>
#include <cnoid/EigenUtil>
#include <cnoid/SceneMarkers>
#include <iostream>

namespace IK{
  class COMConstraint : public IKConstraint
  {
  public:
    COMConstraint()
      :
      C_(Eigen::SparseMatrix<double,Eigen::RowMajor>(0,3)) // 下の定義の箇所で初期化するとUbuntu16でエラーになった
    {}
    //robotの重心をworld座標系のtargetPosに位置させる.
    //  maxError: エラーの頭打ち
    //  precision: 収束判定の閾値
    //  weight: コスト関数の重み. error * weight^2 * error. 0の成分はjacobianやerrorに含まれない
    //状態が更新される度に, 手動でcalcForwardKinematics()とcalcCenterOfMass()を呼ぶ必要が有る.
    const cnoid::BodyPtr& A_robot() const { return A_robot_;}
    cnoid::BodyPtr& A_robot() { return A_robot_;}
    const cnoid::Vector3& A_localp() const { return A_localp_;}
    cnoid::Vector3& A_localp() { return A_localp_;}
    const cnoid::BodyPtr& B_robot() const { return B_robot_;}
    cnoid::BodyPtr& B_robot() { return B_robot_;}
    const cnoid::Vector3& B_localp() const { return B_localp_;}
    cnoid::Vector3& B_localp() { return B_localp_;}
    const cnoid::Matrix3d& eval_R() const { return eval_R_;}
    cnoid::Matrix3d& eval_R() { return eval_R_;}

    // for equality
    const cnoid::Vector3& maxError() const { return maxError_;}
    cnoid::Vector3& maxError() { return maxError_;}
    const cnoid::Vector3& precision() const { return precision_;}
    cnoid::Vector3& precision() { return precision_;}
    const cnoid::Vector3& weight() const { return weight_;}
    cnoid::Vector3& weight() { return weight_;}

    // for inequality. (c * 重心位置)をdlとduの範囲内にする.
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& C() const { return C_;}
    Eigen::SparseMatrix<double,Eigen::RowMajor>& C() { return C_;}
    const cnoid::VectorX& dl() const { return dl_;}
    cnoid::VectorX& dl() { return dl_;}
    const cnoid::VectorX& du() const { return du_;}
    cnoid::VectorX& du() { return du_;}
    const cnoid::VectorX& maxCError() const { return maxCError_;}
    cnoid::VectorX& maxCError() { return maxCError_;}
    const cnoid::VectorX& CPrecision() const { return CPrecision_;}
    cnoid::VectorX& CPrecision() { return CPrecision_;}

    // 収束判定
    bool checkConvergence () override;

    const Eigen::VectorXd& calc_error () override;
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) override;
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) override;
    const Eigen::VectorXd& calc_minineq () override;
    const Eigen::VectorXd& calc_maxineq () override;
  protected:
    cnoid::BodyPtr A_robot_ = nullptr;
    cnoid::Vector3 A_localp_ = cnoid::Vector3::Zero();
    cnoid::BodyPtr B_robot_ = nullptr;
    cnoid::Vector3 B_localp_ = cnoid::Vector3::Zero();
    cnoid::Matrix3d eval_R_ = cnoid::Matrix3d::Identity();

    cnoid::Vector3 maxError_ = 0.1 * cnoid::Vector3::Ones();
    cnoid::Vector3 precision_ = 1e-4 * cnoid::Vector3::Ones();
    cnoid::Vector3 weight_ = cnoid::Vector3::Ones();

    Eigen::SparseMatrix<double,Eigen::RowMajor> C_;
    cnoid::VectorX dl_;
    cnoid::VectorX du_;
    cnoid::VectorX maxCError_;
    cnoid::VectorX CPrecision_;

    cnoid::BodyPtr jacobian_A_robot_ = nullptr;// 前回のjacobian計算時のrobot
    cnoid::BodyPtr jacobian_B_robot_ = nullptr;// 前回のjacobian計算時のrobot
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_full_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_full_local_;
    cnoid::BodyPtr jacobianineq_A_robot_ = nullptr;// 前回のjacobian計算時のrobot
    cnoid::BodyPtr jacobianineq_B_robot_ = nullptr;// 前回のjacobian計算時のrobot
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobianineq_full_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobianineq_full_local_;

  };
}

#endif
