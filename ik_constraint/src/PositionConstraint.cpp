#include <ik_constraint/PositionConstraint.h>
#include <ik_constraint/Jacobian.h>

namespace IK{
  PositionConstraint::PositionConstraint () {
    this->maxError_ << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    this->precision_ << 1e-4, 1e-4, 1e-4, 0.001745, 0.001745, 0.001745;
  }

  // 収束判定
  bool PositionConstraint::checkConvergence () {
    const cnoid::Position& A_pos = (this->A_link_) ? this->A_link_->T() * this->A_localpos_ : this->A_localpos_;
    const cnoid::Position& B_pos = (this->B_link_) ? this->B_link_->T() * this->B_localpos_ : this->B_localpos_;

    cnoid::Vector6 error; // world frame
    const cnoid::Vector3 pos_error = A_pos.translation() - B_pos.translation();
    cnoid::Vector3 rot_error = cnoid::Vector3::Zero();
    if((this->weight_.tail<3>().array() > 0.0).count() == 2 &&
       ((this->eval_link_ == this->A_link_) || (this->eval_link_ == this->B_link_)) ) {
      cnoid::Vector3 axis; // evalR local
      if(this->weight_[3] == 0.0) axis = cnoid::Vector3::UnitX();
      else if(this->weight_[4] == 0.0) axis = cnoid::Vector3::UnitY();
      else if(this->weight_[5] == 0.0) axis = cnoid::Vector3::UnitZ();
      cnoid::Matrix3d eval_R = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;
      // A_axisとB_axisを一致させる
      cnoid::Vector3 A_axis; // world frame
      cnoid::Vector3 B_axis; // world frame
      if(this->eval_link_ == this->A_link_){
        A_axis = eval_R * axis;
        B_axis = B_pos.linear() * A_pos.linear().transpose() * A_axis;
      }else{ // this->eval_link_ == this->B_link_
        B_axis = eval_R * axis;
        A_axis = A_pos.linear() * B_pos.linear().transpose() * B_axis;
      }
      Eigen::Vector3d cross = B_axis.cross(A_axis);
      double dot = std::min(1.0,std::max(-1.0,B_axis.dot(A_axis))); // acosは定義域外のときnanを返す
      if(cross.norm()==0){ // 0度 or 180度.
        if(dot == -1){ // 180度
          if(this->weight_[3] == 0.0) rot_error = eval_R * M_PI * cnoid::Vector3::UnitY();
          else if(this->weight_[4] == 0.0) rot_error = eval_R * M_PI * cnoid::Vector3::UnitZ();
          else if(this->weight_[5] == 0.0) rot_error = eval_R * M_PI * cnoid::Vector3::UnitX();
        }else{
          // rot_error.setZero();
        }
      }else{
        double angle = std::acos(dot); // 0~pi
        Eigen::Vector3d axis_ = cross.normalized(); // include sign
        rot_error = angle * axis_;
      }
    }else{
      const cnoid::AngleAxis angleAxis = cnoid::AngleAxis(A_pos.linear() * B_pos.linear().transpose());
      rot_error = angleAxis.angle()*angleAxis.axis();
    }
    error << pos_error , rot_error;

    cnoid::Matrix3d eval_R = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;
    error.head<3>() = (eval_R.transpose() * error.head<3>()).eval();
    error.tail<3>() = (eval_R.transpose() * error.tail<3>()).eval();

    // 収束判定と、ついでにcalc_errorの返り値の計算
    if(this->error_.rows()!=(this->weight_.array() > 0.0).count()) this->error_ = Eigen::VectorXd((this->weight_.array() > 0.0).count());
    bool converged = true;
    int idx=0;
    for(size_t i=0; i<6; i++){
      if(this->weight_[i]>0.0) {
        if(std::fabs(error[i]) > this->precision_[i]) converged = false;
        this->error_[idx] = std::min(std::max(error[i],-this->maxError_[i]),this->maxError_[i]) * this->weight_[i];
        idx++;
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "PositionConstraint" << std::endl;
      std::cerr << "A_pos" << std::endl;
      std::cerr << A_pos.translation().transpose() << std::endl;
      std::cerr << A_pos.linear() << std::endl;
      std::cerr << "B_pos" << std::endl;
      std::cerr << B_pos.translation().transpose() << std::endl;
      std::cerr << B_pos.linear() << std::endl;
    }

    return converged;
  }

  // エラーを返す. A-B. world系. QPで用いる
  const Eigen::VectorXd& PositionConstraint::calc_error () {
    if(this->debuglevel_>=1){
      std::cerr << "PositionConstraint" << std::endl;
      std::cerr << "error" << std::endl;
      std::cerr << this->error_.transpose() << std::endl;
    }
    return this->error_;
  }

  // ヤコビアンを返す. bodyのroot6dof+全関節が変数
  const Eigen::SparseMatrix<double,Eigen::RowMajor>& PositionConstraint::calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_joints_same(joints,this->jacobian_joints_)
       || this->A_link_ != this->jacobian_A_link_
       || this->B_link_ != this->jacobian_B_link_){
      this->jacobian_joints_ = joints;
      this->jacobian_A_link_ = this->A_link_;
      this->jacobian_B_link_ = this->B_link_;

      IK::calc6DofJacobianShape(this->jacobian_joints_,//input
                                this->jacobian_A_link_,//input
                                this->jacobian_B_link_,//input
                                this->jacobian_full_,
                                this->jacobianColMap_,
                                this->path_A_joints_,
                                this->path_B_joints_,
                                this->path_BA_joints_,
                                this->path_BA_joints_numUpwardConnections_
                                );
    }

    IK::calc6DofJacobianCoef(this->jacobian_joints_,//input
                             this->jacobian_A_link_,//input
                             this->A_localpos_,//input
                             this->jacobian_B_link_,//input
                             this->B_localpos_,//input
                             this->jacobianColMap_,//input
                             this->path_A_joints_,//input
                             this->path_B_joints_,//input
                             this->path_BA_joints_,//input
                             this->path_BA_joints_numUpwardConnections_,//input
                             this->jacobian_full_
                             );

    cnoid::Matrix3d eval_R_dense = (this->eval_link_) ? this->eval_link_->R() * this->eval_localR_ : this->eval_localR_;
    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R.insert(i,j) = eval_R_dense(i,j);
    this->jacobian_full_local_.resize(this->jacobian_full_.rows(), this->jacobian_full_.cols());
    this->jacobian_full_local_.topRows<3>() = eval_R.transpose() * this->jacobian_full_.topRows<3>();
    this->jacobian_full_local_.bottomRows<3>() = eval_R.transpose() * this->jacobian_full_.bottomRows<3>();

    this->jacobian_.resize((this->weight_.array() > 0.0).count(),this->jacobian_full_local_.cols());
    for(size_t i=0, idx=0;i<6;i++){
      if(this->weight_[i]>0.0) {
        this->jacobian_.row(idx) = this->weight_[i] * this->jacobian_full_local_.row(i);
        idx++;
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "PositionConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
    }
    return this->jacobian_;
  }

  std::vector<cnoid::SgNodePtr>& PositionConstraint::getDrawOnObjects(){
    if(!this->lines_){
      this->lines_ = new cnoid::SgLineSet;
      this->lines_->setLineWidth(1.0);
      this->lines_->getOrCreateColors()->resize(4);
      this->lines_->getOrCreateColors()->at(0) = cnoid::Vector3f(1.0,1.0,1.0);
      this->lines_->getOrCreateColors()->at(1) = cnoid::Vector3f(1.0,0.0,0.0);
      this->lines_->getOrCreateColors()->at(2) = cnoid::Vector3f(0.0,1.0,0.0);
      this->lines_->getOrCreateColors()->at(3) = cnoid::Vector3f(0.0,0.0,1.0);
      // A, A_x, A_y, A_z, B, B_x, B_y, B_z
      this->lines_->getOrCreateVertices()->resize(8);
      this->lines_->colorIndices().resize(0);
      this->lines_->addLine(0,1); this->lines_->colorIndices().push_back(1); this->lines_->colorIndices().push_back(1);
      this->lines_->addLine(0,2); this->lines_->colorIndices().push_back(2); this->lines_->colorIndices().push_back(2);
      this->lines_->addLine(0,3); this->lines_->colorIndices().push_back(3); this->lines_->colorIndices().push_back(3);
      this->lines_->addLine(4,5); this->lines_->colorIndices().push_back(1); this->lines_->colorIndices().push_back(1);
      this->lines_->addLine(4,6); this->lines_->colorIndices().push_back(2); this->lines_->colorIndices().push_back(2);
      this->lines_->addLine(4,7); this->lines_->colorIndices().push_back(3); this->lines_->colorIndices().push_back(3);
      this->lines_->addLine(0,4); this->lines_->colorIndices().push_back(0); this->lines_->colorIndices().push_back(0);

      this->drawOnObjects_ = std::vector<cnoid::SgNodePtr>{this->lines_};
    }

    const cnoid::Position& A_pos = (this->A_link_) ? this->A_link_->T() * this->A_localpos_ : this->A_localpos_;
    const cnoid::Position& B_pos = (this->B_link_) ? this->B_link_->T() * this->B_localpos_ : this->B_localpos_;

    this->lines_->getOrCreateVertices()->at(0) = A_pos.translation().cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(1) = (A_pos * (0.05 * cnoid::Vector3::UnitX())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(2) = (A_pos * (0.05 * cnoid::Vector3::UnitY())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(3) = (A_pos * (0.05 * cnoid::Vector3::UnitZ())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(4) = B_pos.translation().cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(5) = (B_pos * (0.05 * cnoid::Vector3::UnitX())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(6) = (B_pos * (0.05 * cnoid::Vector3::UnitY())).cast<cnoid::Vector3f::Scalar>();
    this->lines_->getOrCreateVertices()->at(7) = (B_pos * (0.05 * cnoid::Vector3::UnitZ())).cast<cnoid::Vector3f::Scalar>();

    return this->drawOnObjects_;
  }
}
