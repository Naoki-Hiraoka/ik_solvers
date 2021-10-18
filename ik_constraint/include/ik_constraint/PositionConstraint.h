#ifndef POSITIONCONSTRAINT_H
#define POSITIONCONSTRAINT_H

#include <ik_constraint/IKConstraint.h>
#include <cnoid/EigenUtil>
#include <cnoid/LinkPath>
#include <iostream>

namespace IK{
  class PositionConstraint : public IKConstraint
  {
  public:
    PositionConstraint();

    //A_link中のA_localposの部位とB_link中のB_localposの部位を一致させる.
    //リンクがnullptrならworld座標系を意味する
    //  maxError: エラーの頭打ち eval系
    //  precision: 収束判定の閾値 eval系
    //  weight: コスト関数の重み. error * weight^2 * error. 0の成分はjacobianやerrorに含まれない. eval系
    //状態が更新される度に, 手動でcalcForwardKinematics()を呼ぶ必要が有る.
    const cnoid::LinkPtr& A_link() const { return A_link_;}
    cnoid::LinkPtr& A_link() { return A_link_;}
    const cnoid::Position& A_localpos() const { return A_localpos_;}
    cnoid::Position& A_localpos() { return A_localpos_;}
    const cnoid::LinkPtr& B_link() const { return B_link_;}
    cnoid::LinkPtr& B_link() { return B_link_;}
    const cnoid::Position& B_localpos() const { return B_localpos_;}
    cnoid::Position& B_localpos() { return B_localpos_;}
    const cnoid::Vector6& maxError() const { return maxError_;}
    cnoid::Vector6& maxError() { return maxError_;}
    const cnoid::Vector6& precision() const { return precision_;}
    cnoid::Vector6& precision() { return precision_;}
    const cnoid::Vector6& weight() const { return weight_;}
    cnoid::Vector6& weight() { return weight_;}
    const cnoid::LinkPtr& eval_link() const { return eval_link_;}
    cnoid::LinkPtr& eval_link() { return eval_link_;}
    const cnoid::Matrix3d& eval_localR() const { return eval_localR_;}
    cnoid::Matrix3d& eval_localR() { return eval_localR_;}

    // 収束判定
    bool checkConvergence () override;

    // for debug view
    std::vector<cnoid::SgNodePtr>& getDrawOnObjects() override;

    // エラーを返す. A-B. world系. QPで用いる
    const Eigen::VectorXd& calc_error () override;

    // ヤコビアンを返す. bodyのroot6dof+全関節が変数
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) override;

    // コスト(エラーの二乗和)を返す. 非線形最適化で用いる
    // TODO

    // gradient(-ヤコビアン^T*エラー)を返す
    // TODO

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  private:
    cnoid::LinkPtr A_link_ = nullptr;
    cnoid::Position A_localpos_ = cnoid::Position::Identity();
    cnoid::LinkPtr B_link_ = nullptr;
    cnoid::Position B_localpos_ = cnoid::Position::Identity();
    cnoid::Vector6 maxError_;
    cnoid::Vector6 precision_;
    cnoid::Vector6 weight_ = cnoid::Vector6::Ones();
    cnoid::LinkPtr eval_link_ = nullptr;
    cnoid::Matrix3d eval_localR_ = cnoid::Matrix3d::Identity();

    cnoid::SgLineSetPtr lines_;

    std::vector<cnoid::LinkPtr> path_A_joints_;
    std::vector<cnoid::LinkPtr> path_B_joints_;
    std::vector<cnoid::LinkPtr> path_BA_joints_;
    int path_BA_joints_numUpwardConnections_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_full_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_full_local_;
    cnoid::LinkPtr jacobian_A_link_ = nullptr;// 前回のjacobian計算時のA_link
    cnoid::LinkPtr jacobian_B_link_ = nullptr;// 前回のjacobian計算時のB_link
  };
}

#endif
