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
    //robotの重心をworld座標系のtargetPosに位置させる.
    //  maxError: エラーの頭打ち
    //  precision: 収束判定の閾値
    //  weight: コスト関数の重み. error * weight^2 * error. 0の成分はjacobianやerrorに含まれない
    //状態が更新される度に, 手動でcalcForwardKinematics()とcalcCenterOfMass()を呼ぶ必要が有る.
    const cnoid::BodyPtr& robot() const { return robot_;}
    cnoid::BodyPtr& robot() { return robot_;}
    const cnoid::Vector3& targetPos() const { return targetPos_;}
    cnoid::Vector3& targetPos() { return targetPos_;}
    const cnoid::Vector3& maxError() const { return maxError_;}
    cnoid::Vector3& maxError() { return maxError_;}
    const cnoid::Vector3& precision() const { return precision_;}
    cnoid::Vector3& precision() { return precision_;}
    const cnoid::Vector3& weight() const { return weight_;}
    cnoid::Vector3& weight() { return weight_;}

    // 収束判定
    bool checkConvergence () override;

    // for debug view
    std::vector<cnoid::SgNodePtr>& getDrawOnObjects() override;

    const Eigen::VectorXd& calc_error () override;
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobian (const std::vector<cnoid::BodyPtr>& bodies) override;
  protected:
    cnoid::BodyPtr robot_ = nullptr;
    cnoid::Vector3 targetPos_ = cnoid::Vector3::Zero();
    cnoid::Vector3 maxError_ = 0.1 * cnoid::Vector3::Ones();
    cnoid::Vector3 precision_ = 1e-4 * cnoid::Vector3::Ones();
    cnoid::Vector3 weight_ = cnoid::Vector3::Ones();

    cnoid::BodyPtr jacobian_robot_ = nullptr;// 前回のjacobian計算時のrobot

    cnoid::CrossMarkerPtr COMMarker_;
    cnoid::CrossMarkerPtr targetMarker_;
  };
}

#endif
