#include <ik_constraint/PositionConstraint.h>

namespace IK{
  void pushBackTripletList(std::vector<Eigen::Triplet<double> >& tripletList, const cnoid::LinkPtr& joint, int idx){
    if(joint->isFreeJoint()){
      for(size_t d=0;d<6;d++){
        tripletList.push_back(Eigen::Triplet<double>(d,idx+d,1));
      }
      //  0     p[2] -p[1]
      // -p[2]  0     p[0]
      //  p[1] -p[0]  0
      tripletList.push_back(Eigen::Triplet<double>(0,idx+4,1));
      tripletList.push_back(Eigen::Triplet<double>(0,idx+5,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx+3,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx+5,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx+3,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx+4,1));

    } else if(joint->isRotationalJoint()){
      tripletList.push_back(Eigen::Triplet<double>(0,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(3,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(4,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(5,idx,1));

    } else if(joint->isPrismaticJoint()){
      tripletList.push_back(Eigen::Triplet<double>(0,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx,1));

    }
  }

  void fillJacobian(Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian, const cnoid::Vector3& target_p, const cnoid::LinkPtr& joint, int idx, int sign){
    if(joint->isFreeJoint()){
      //root 6dof
      for(size_t j=0;j<6;j++){
        jacobian.coeffRef(j,idx+j) = sign;
      }
      cnoid::Vector3 dp = target_p - joint->p();
      //  0     p[2] -p[1]
      // -p[2]  0     p[0]
      //  p[1] -p[0]  0
      jacobian.coeffRef(0,idx+4)=sign*dp[2];
      jacobian.coeffRef(0,idx+5)=-sign*dp[1];
      jacobian.coeffRef(1,idx+3)=-sign*dp[2];
      jacobian.coeffRef(1,idx+5)=sign*dp[0];
      jacobian.coeffRef(2,idx+3)=sign*dp[1];
      jacobian.coeffRef(2,idx+4)=-sign*dp[0];

    } else if(joint->isRotationalJoint()){
      cnoid::Vector3 omega = joint->R() * joint->a();
      cnoid::Vector3 dp = omega.cross(target_p - joint->p());
      jacobian.coeffRef(0,idx)=sign*dp[0];
      jacobian.coeffRef(1,idx)=sign*dp[1];
      jacobian.coeffRef(2,idx)=sign*dp[2];
      jacobian.coeffRef(3,idx)=sign*omega[0];
      jacobian.coeffRef(4,idx)=sign*omega[1];
      jacobian.coeffRef(5,idx)=sign*omega[2];

    } else if(joint->isPrismaticJoint()){
      cnoid::Vector3 dp = joint->R() * joint->d();
      jacobian.coeffRef(0,idx)=sign*dp[0];
      jacobian.coeffRef(1,idx)=sign*dp[1];
      jacobian.coeffRef(2,idx)=sign*dp[2];
    }
  }

  PositionConstraint::PositionConstraint () {
    this->maxError_ << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    this->precision_ << 1e-4, 1e-4, 1e-4, 0.001745, 0.001745, 0.001745;
  }

  // 収束判定
  bool PositionConstraint::checkConvergence () {
    const cnoid::Position& A_pos = (this->A_link_) ? this->A_link_->T() * this->A_localpos_ : this->A_localpos_;
    const cnoid::Position& B_pos = (this->B_link_) ? this->B_link_->T() * this->B_localpos_ : this->B_localpos_;
    cnoid::Vector6 error;
    error << A_pos.translation() - B_pos.translation() , cnoid::omegaFromRot(A_pos.linear() * B_pos.linear().transpose());

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
      this->jacobianColMap_.clear();
      int num_variables = 0;
      for(size_t i=0;i<this->jacobian_joints_.size();i++){
        this->jacobianColMap_[this->jacobian_joints_[i]] = num_variables;
        num_variables += this->getJointDOF(this->jacobian_joints_[i]);
      }
      this->jacobian_A_link_ = this->A_link_;
      this->jacobian_B_link_ = this->B_link_;

      std::vector<Eigen::Triplet<double> > tripletList;
      tripletList.reserve(100);//適当

      if(!this->A_link_ || !this->B_link_ || !(this->A_link_->body() == this->B_link_->body())){
        // A, Bが関節を共有しない. 別々に処理すれば良い
        for(size_t i=0;i<2;i++){//0:A_link, 1:B_link
          cnoid::LinkPtr target_link = i ? this->B_link_ : this->A_link_;
          std::vector<cnoid::LinkPtr>& path_joints = 1 ? this->path_B_joints_ : this->path_A_joints_;

          if(!target_link) continue;//world固定なので飛ばす

          path_joints.clear();
          cnoid::LinkPath path(target_link);
          for(size_t j=0;j<path.size();j++){
            path_joints.push_back(path[j]);
          }

          for(size_t j=0;j<path_joints.size();j++){
            cnoid::LinkPtr joint = path_joints[j];
            if(this->jacobianColMap_.find(joint)==this->jacobianColMap_.end()) continue;
            int idx = this->jacobianColMap_[joint];
            pushBackTripletList(tripletList,joint,idx);
          }
        }
      } else { //if(!A_link || !B_link || !(A_link->body() == B_link->body()))
        //A,Bが関節を共有する. 一つのpathで考える
        this->path_BA_joints_.clear();
        {
          cnoid::LinkPath path(this->B_link_,this->A_link_);
          size_t j=0;
          for(;!path.isDownward(j);j++) this->path_BA_joints_.push_back(path[j]);
          this->path_BA_joints_numUpwardConnections_ = j;
          j++;
          for(;j<path.size();j++) this->path_BA_joints_.push_back(path[j]);
        }

        for(size_t j=0;j<this->path_BA_joints_.size();j++){
          cnoid::LinkPtr joint = this->path_BA_joints_[j];
          if(this->jacobianColMap_.find(joint)==this->jacobianColMap_.end()) continue;
          int idx = this->jacobianColMap_[joint];
          pushBackTripletList(tripletList,joint,idx);
        }
      }

      this->jacobian_full_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(6,num_variables);
      this->jacobian_full_.setFromTriplets(tripletList.begin(), tripletList.end());

    }

    if(!this->A_link_ || !this->B_link_ || !(this->A_link_->body() == this->B_link_->body())){
      // A, Bが関節を共有しない. 別々に処理すれば良い
      for(size_t i=0;i<2;i++){//0:A_link, 1:B_link
        int sign = i ? -1 : 1;
        cnoid::LinkPtr target_link = i ? this->B_link_ : this->A_link_;
        const cnoid::Position& target_localpos = i ? this->B_localpos_ : this->A_localpos_;
        std::vector<cnoid::LinkPtr>& path_joints = i ? this->path_B_joints_ : this->path_B_joints_;

        if(!target_link) continue;//world固定なので飛ばす

        const cnoid::Position target_position = target_link->T() * target_localpos;
        const cnoid::Vector3 target_p = target_position.translation();

        for(size_t j=0;j<path_joints.size();j++){
          cnoid::LinkPtr joint = path_joints[j];
          if(this->jacobianColMap_.find(joint)==this->jacobianColMap_.end()) continue;
          int idx = this->jacobianColMap_[joint];
          fillJacobian(this->jacobian_full_,target_p,joint,idx,sign);
        }
      }
    } else { //if(!A_link || !B_link || !(A_link->body() == B_link->body()))
      //A,Bが関節を共有する. 一つのpathで考える
      const cnoid::Vector3& target_p = this->A_link_->T() * this->A_localpos_.translation();
      for(size_t j=0;j<this->path_BA_joints_.size();j++){
          cnoid::LinkPtr joint = this->path_BA_joints_[j];
          if(this->jacobianColMap_.find(joint)==this->jacobianColMap_.end()) continue;
          int idx = this->jacobianColMap_[joint];
          if(j<this->path_BA_joints_numUpwardConnections_){
            fillJacobian(this->jacobian_full_,target_p,joint,idx,-1);
          } else {
            fillJacobian(this->jacobian_full_,target_p,joint,idx,1);
          }
      }
    }

    this->jacobian_.resize((this->weight_.array() > 0.0).count(),this->jacobian_full_.cols());
    for(size_t i=0, idx=0;i<6;i++){
      if(this->weight_[i]>0.0) {
        this->jacobian_.row(idx) = this->weight_[i] * this->jacobian_full_.row(i);
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
