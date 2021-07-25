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
  const Eigen::SparseMatrix<double,Eigen::RowMajor>& IKConstraint::calc_jacobian (const std::vector<cnoid::BodyPtr>& bodies){
    if(!this->is_bodies_same(bodies,this->jacobian_bodies_)){
      this->jacobian_bodies_ = bodies;

      int num_variables = 0;
      for(size_t i=0;i<this->jacobian_bodies_.size();i++){
        num_variables += 6 + this->jacobian_bodies_[i]->numJoints();
      }
      this->jacobian_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,num_variables);
    }
    return this->jacobian_;
  }

  // 不等式制約のヤコビアンを返す bodyのroot6dof+全関節が変数
  const Eigen::SparseMatrix<double,Eigen::RowMajor>& IKConstraint::calc_jacobianineq (const std::vector<cnoid::BodyPtr>& bodies){
    if(!this->is_bodies_same(bodies,this->jacobianineq_bodies_)){
      this->jacobianineq_bodies_ = bodies;

      int num_variables = 0;
      for(size_t i=0;i<this->jacobianineq_bodies_.size();i++){
        num_variables += 6 + this->jacobianineq_bodies_[i]->numJoints();
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

  bool IKConstraint::is_bodies_same(const std::vector<cnoid::BodyPtr>& bodies1,const std::vector<cnoid::BodyPtr>& bodies2){
    if (bodies1.size() != bodies2.size() ) return false;
    for(size_t i=0;i<bodies1.size();i++){
      if (bodies1[i] != bodies2[i] ) return false;
    }
    return true;
  }
}
