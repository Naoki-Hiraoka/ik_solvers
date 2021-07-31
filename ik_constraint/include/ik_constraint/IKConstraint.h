#ifndef IKCONSTRAINT_H
#define IKCONSTRAINT_H

#include <cnoid/Body>
#include <cnoid/SceneDrawables>
#include <Eigen/Sparse>
#include <unordered_map>

namespace IK{
  class IKConstraint
  {
  public:

    // 必ず, 状態更新 -> checkConvergence -> (getDrawOnObjects) -> calc_error -> calc_jacobian -> calc_jacobianineq -> calc_min/maxineq -> 状態更新 の順で呼ぶので，同じ処理を何度も行うのではなく最初に呼ばれる関数で1回だけ行って以降はキャッシュを使ってよい
    // 収束判定
    virtual bool checkConvergence ();
    // for debug view
    virtual std::vector<cnoid::SgNodePtr>& getDrawOnObjects();
    // 等式制約のエラーを返す.
    virtual const Eigen::VectorXd& calc_error ();
    // 等式制約のヤコビアンを返す. bodyのroot6dof+全関節が変数
    virtual const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobian (const std::vector<cnoid::LinkPtr>& joints);
    // 不等式制約のヤコビアンを返す. bodyのroot6dof+全関節が変数
    virtual const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints);
    // 不等式制約のmin値を返す
    virtual const Eigen::VectorXd& calc_minineq ();
    // 不等式制約のmax値を返す
    virtual const Eigen::VectorXd& calc_maxineq ();

    // コスト(エラーの二乗和)を返す. 非線形最適化で用いる
    // TODO

    // gradient(-ヤコビアン^T*エラー)を返す
    // TODO

    const int& debuglevel() const { return debuglevel_;}
    int& debuglevel() { return debuglevel_;}
  protected:
    bool is_joints_same(const std::vector<cnoid::LinkPtr>& joints1,const std::vector<cnoid::LinkPtr>& joints2);
    size_t getJointDOF(const cnoid::LinkPtr& joint);

    int debuglevel_ = 0;

    Eigen::VectorXd error_;
    std::vector<cnoid::LinkPtr> jacobian_joints_; // 前回のjacobian計算時のjoints
    std::unordered_map<cnoid::LinkPtr,int> jacobianColMap_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobian_;
    Eigen::VectorXd minineq_;
    Eigen::VectorXd maxineq_;
    std::vector<cnoid::LinkPtr> jacobianineq_joints_; // 前回のjacobianineq計算時のjoints
    std::unordered_map<cnoid::LinkPtr,int> jacobianineqColMap_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> jacobianineq_;
    std::vector<cnoid::SgNodePtr> drawOnObjects_;
  };
}

#endif
