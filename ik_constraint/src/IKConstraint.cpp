#include <ik_constraint/IKConstraint.h>

namespace IK{
  // 収束判定
  bool IKConstraint::checkConvergence (){
    return true;
  }

  // for debug view
  std::vector<cnoid::SgNodePtr>& IKConstraint::getDrawOnObjects(){
    return this->drawOnObjects_;
  }

  // 等式制約のエラーを返す.
  const Eigen::VectorXd& IKConstraint::calc_error (){
    return this->error_;
  }

  // 等式制約のヤコビアンを返す. bodyのroot6dof+全関節が変数
  const Eigen::SparseMatrix<double,Eigen::RowMajor>& IKConstraint::calc_jacobian (const std::vector<cnoid::LinkPtr>& joints){
    if(!this->is_joints_same(joints,this->jacobian_joints_)){
      this->jacobian_joints_ = joints;
      this->jacobianColMap_.clear();
      int num_variables = 0;
      for(size_t i=0;i<this->jacobian_joints_.size();i++){
        jacobianColMap_[this->jacobian_joints_[i]] = num_variables;
        num_variables += this->getJointDOF(this->jacobian_joints_[i]);
      }

      this->jacobian_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,num_variables);
    }
    return this->jacobian_;
  }

  // 不等式制約のヤコビアンを返す bodyのroot6dof+全関節が変数
  const Eigen::SparseMatrix<double,Eigen::RowMajor>& IKConstraint::calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints){
    if(!this->is_joints_same(joints,this->jacobianineq_joints_)){
      this->jacobianineq_joints_ = joints;
      this->jacobianineqColMap_.clear();
      int num_variables = 0;
      for(size_t i=0;i<this->jacobianineq_joints_.size();i++){
        jacobianineqColMap_[this->jacobian_joints_[i]] = num_variables;
        num_variables += this->getJointDOF(this->jacobianineq_joints_[i]);
      }
      this->jacobianineq_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,num_variables);
    }
    return this->jacobianineq_;
  }

  // 不等式制約のmin値を返す
  const Eigen::VectorXd& IKConstraint::calc_minineq (){
    return this->minineq_;
  }

  // 不等式制約のmax値を返す
  const Eigen::VectorXd& IKConstraint::calc_maxineq (){
    return this->maxineq_;
  }

  // コスト(エラーの二乗和)を返す. 非線形最適化で用いる
  // TODO

  // gradient(-ヤコビアン^T*エラー)を返す
  // TODO

  bool IKConstraint::is_joints_same(const std::vector<cnoid::LinkPtr>& joints1,const std::vector<cnoid::LinkPtr>& joints2){
    if (joints1.size() != joints2.size() ) return false;
    for(size_t i=0;i<joints1.size();i++){
      if (joints1[i] != joints2[i] ) return false;
    }
    return true;
  }

  size_t IKConstraint::getJointDOF(const cnoid::LinkPtr& joint) {
    if(joint->isRevoluteJoint() || joint->isPrismaticJoint()) return 1;
    else if(joint->isFreeJoint()) return 6;
    else return 0;
  }

}
