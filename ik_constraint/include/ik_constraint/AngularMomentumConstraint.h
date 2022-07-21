#ifndef ANGULARMOMENTUMCONSTRAINT_H
#define ANGULARMOMENTUMCONSTRAINT_H

#include <ik_constraint/IKConstraint.h>
#include <cnoid/EigenUtil>
#include <cnoid/SceneMarkers>
#include <iostream>

namespace IK{
  class AngularMomentumConstraint : public IKConstraint
  {
  public:
    //robotの重心周りの角運動量を目標の値[kg m^2/s]に一致させる
    //  内部の処理では角運動量を重心周りのイナーシャで割って、[rad/s]の次元で扱う
    //  dt: [s]
    //  maxError: エラーの頭打ち[rad]
    //  precision: 収束判定の閾値[rad]
    //  weight: コスト関数の重み. error * weight^2 * error.
    //  targetAngularMomentum: 重心周り. ワールド座標系. [kg m^2/s]
    //状態が更新される度に, 手動でcalcForwardKinematics()とcalcCenterOfMass()を呼ぶ必要が有る.
    const cnoid::BodyPtr& robot() const { return robot_;}
    cnoid::BodyPtr& robot() { return robot_;}
    const cnoid::Matrix3d& eval_R() const { return eval_R_;}
    cnoid::Matrix3d& eval_R() { return eval_R_;}

    const cnoid::Vector3& targetAngularMomentum() const { return targetAngularMomentum_;}
    cnoid::Vector3& targetAngularMomentum() { return targetAngularMomentum_;}
    const double& dt() const { return dt_;}
    double& dt() { return dt_;}

    const cnoid::Vector3& maxError() const { return maxError_;}
    cnoid::Vector3& maxError() { return maxError_;}
    const cnoid::Vector3& precision() const { return precision_;}
    cnoid::Vector3& precision() { return precision_;}
    const cnoid::Vector3& weight() const { return weight_;}
    cnoid::Vector3& weight() { return weight_;}

    // 収束判定
    bool checkConvergence () override;

    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) override;
    const Eigen::VectorXd& calc_error () override;
  protected:
    cnoid::BodyPtr robot_ = nullptr;
    cnoid::Matrix3d eval_R_ = cnoid::Matrix3d::Identity();

    cnoid::Vector3 targetAngularMomentum_ = cnoid::Vector3::Ones();
    double dt_;

    cnoid::Vector3 maxError_ = 0.1 * cnoid::Vector3::Ones();
    cnoid::Vector3 precision_ = 1e-4 * cnoid::Vector3::Ones();
    cnoid::Vector3 weight_ = cnoid::Vector3::Ones();

    cnoid::BodyPtr jacobian_robot_ = nullptr;// 前回のjacobian計算時のrobot
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_full_;

    static void calcAngularMomentumJacobianShape(const std::vector<cnoid::LinkPtr>& joints,//input
                                                 const cnoid::BodyPtr& A_robot,//input
                                                 const cnoid::BodyPtr& B_robot,//input
                                                 Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian,//output
                                                 std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap //output
                                                 );
    static void calcAngularMomentumJacobianCoef(const std::vector<cnoid::LinkPtr>& joints,//input
                                                const cnoid::BodyPtr& A_robot,//input
                                                const cnoid::BodyPtr& B_robot,//input
                                                std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap, //input
                                                Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian//output
                                                );
  };

  namespace cnoid18 {
    // choreonoidのrelease1.7の calcAngularMomentumJacobianにはバグがあり、開発版ではhttps://github.com/s-nakaoka/choreonoid/pull/234 で修正されている. 修正された版の関数(https://github.com/choreonoid/choreonoid/blob/master/src/Body/Jacobian.cpp )を使う
    void calcAngularMomentumJacobian(cnoid::Body* body, cnoid::Link* base, Eigen::MatrixXd& H);
  }
}

#endif
