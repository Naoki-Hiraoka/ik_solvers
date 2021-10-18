#include <ik_constraint/COMConstraint.h>
#include <cnoid/Jacobian>

namespace IK{
  bool COMConstraint::checkConvergence () {
    bool converged = true;

    // A - B
    cnoid::Vector3 A_p = this->A_robot_ ? this->A_robot_->centerOfMass() + this->A_localp() : this->A_localp();
    cnoid::Vector3 B_p = this->B_robot_ ? this->B_robot_->centerOfMass() + this->B_localp() : this->B_localp();
    cnoid::Vector3 error = A_p - B_p;
    error = (this->eval_R_.transpose() * error).eval();

    // equalityの収束判定と、ついでにcalc_errorの返り値の計算
    if(this->error_.rows()!=(this->weight_.array() > 0.0).count()) this->error_ = Eigen::VectorXd((this->weight_.array() > 0.0).count());
    int idx=0;
    for(size_t i=0; i<3; i++){
      if(this->weight_[i]>0.0) {
        if(std::fabs(error[i]) > this->precision_[i]) converged = false;
        this->error_[idx] = std::min(std::max(error[i],-this->maxError_[i]),this->maxError_[i]) * this->weight_[i];
        idx++;
      }
    }

    // inequalityの収束判定と、ついでにcalc_max/minineqの返り値の計算
    if(this->C_.cols() != 3 ||
       this->C_.rows() != this->dl_.size() ||
       this->du_.size() != this->dl_.size()){
      std::cerr << "\x1b[31m" << "[COMConstraint::checkConvergence] dimension mismatch" << "\x1b[39m" << std::endl;
      this->C_.resize(0,3);
      this->du_.resize(0);
      this->dl_.resize(0);
    }
    if(this->C_.rows() != this->CPrecision_.size()){
      this->CPrecision_ = cnoid::VectorX::Ones(this->C_.rows()) * 1e-4;
    }
    if(this->C_.rows() != this->maxCError_.size() ){
      this->maxCError_ = cnoid::VectorX::Ones(this->C_.rows()) * 0.1;
    }
    if(this->minineq_.rows() != this->C_.rows()) this->minineq_ = Eigen::VectorXd(this->C_.rows());
    if(this->maxineq_.rows() != this->C_.rows()) this->maxineq_ = Eigen::VectorXd(this->C_.rows());

    cnoid::VectorX Ce = this->C_ * error;
    cnoid::VectorX u = this->du_ - Ce;
    cnoid::VectorX l = this->dl_ - Ce;
    for(size_t i=0; i<u.size(); i++){
      if(u[i] < -this->CPrecision_[i]) converged = false;
      this->maxineq_[i] = std::max(u[i],-this->maxCError_[i]);
      if(l[i] > this->CPrecision_[i]) converged = false;
      this->minineq_[i] = std::min(l[i],this->maxCError_[i]);
    }

    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "A COM "<<A_p.transpose() << std::endl;
      std::cerr << "B COM "<<B_p.transpose() << std::endl;
    }

    return converged;
  }

  const Eigen::VectorXd& COMConstraint::calc_error (){
    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "error" << std::endl;
      std::cerr << this->error_.transpose() << std::endl;
    }

    return this->error_;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& COMConstraint::calc_jacobian (const std::vector<cnoid::LinkPtr>& joints) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_joints_same(joints,this->jacobian_joints_)
       || this->A_robot_ != this->jacobian_A_robot_
       || this->B_robot_ != this->jacobian_B_robot_){
      this->jacobian_joints_ = joints;
      this->jacobian_A_robot_ = this->A_robot_;
      this->jacobian_B_robot_ = this->B_robot_;

      COMConstraint::calcCMJacobianShape(this->jacobian_joints_,
                                         this->jacobian_A_robot_,
                                         this->jacobian_B_robot_,
                                         this->jacobian_full_,
                                         this->jacobianColMap_);
    }

    COMConstraint::calcCMJacobianCoef(this->jacobian_joints_,
                                      this->jacobian_A_robot_,
                                      this->jacobian_B_robot_,
                                      this->jacobianColMap_,
                                      this->jacobian_full_);

    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R.insert(i,j) = this->eval_R_(i,j);
    this->jacobian_full_local_= eval_R.transpose() * this->jacobian_full_;

    this->jacobian_.resize((this->weight_.array() > 0.0).count(),this->jacobian_full_local_.cols());
    for(size_t i=0, idx=0;i<3;i++){
      if(this->weight_[i]>0.0) {
        this->jacobian_.row(idx) = this->weight_[i] * this->jacobian_full_local_.row(i);
        idx++;
      }
    }

    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "jacobian" << std::endl;
      std::cerr << this->jacobian_ << std::endl;
    }

    return this->jacobian_;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& COMConstraint::calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_joints_same(joints,this->jacobianineq_joints_)
       || this->A_robot_ != this->jacobianineq_A_robot_
       || this->B_robot_ != this->jacobianineq_B_robot_){
      this->jacobianineq_joints_ = joints;
      this->jacobianineq_A_robot_ = this->A_robot_;
      this->jacobianineq_B_robot_ = this->B_robot_;

      COMConstraint::calcCMJacobianShape(this->jacobianineq_joints_,
                                         this->jacobianineq_A_robot_,
                                         this->jacobianineq_B_robot_,
                                         this->jacobianineq_full_,
                                         this->jacobianineqColMap_);
    }

    COMConstraint::calcCMJacobianCoef(this->jacobianineq_joints_,
                                      this->jacobianineq_A_robot_,
                                      this->jacobianineq_B_robot_,
                                      this->jacobianineqColMap_,
                                      this->jacobianineq_full_);

    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R.insert(i,j) = this->eval_R_(i,j);

    this->jacobianineq_ = this->C_ * eval_R.transpose() * this->jacobianineq_full_;

    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "jacobiaineq" << std::endl;
      std::cerr << this->jacobianineq_ << std::endl;
    }

    return this->jacobianineq_;
  }

  const Eigen::VectorXd& COMConstraint::calc_minineq () {
    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "minineq" << std::endl;
      std::cerr << this->minineq_ << std::endl;
    }

    return this->minineq_;
  }

  const Eigen::VectorXd& COMConstraint::calc_maxineq () {
    if(this->debuglevel_>=1){
      std::cerr << "COMConstraint" << std::endl;
      std::cerr << "maxineq" << std::endl;
      std::cerr << this->maxineq_ << std::endl;
    }

    return this->maxineq_;
  }

  void COMConstraint::calcCMJacobianShape(const std::vector<cnoid::LinkPtr>& joints, const cnoid::BodyPtr& A_robot, const cnoid::BodyPtr& B_robot, Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian, std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap){
    jacobianColMap.clear();
    int cols = 0;
    for(size_t i=0;i<joints.size();i++){
      jacobianColMap[joints[i]] = cols;
      cols += COMConstraint::getJointDOF(joints[i]);
    }
    jacobian = Eigen::SparseMatrix<double,Eigen::RowMajor>(3,cols);

    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.reserve(100);//適当

    if(A_robot != B_robot){
      for(int i=0;i<2;i++){
        cnoid::BodyPtr robot = (i==0) ? A_robot : B_robot;
        if(!robot) continue;
        if(jacobianColMap.find(robot->rootLink()) != jacobianColMap.end()){
          int idx = jacobianColMap[robot->rootLink()];
          for(size_t row=0;row<3;row++){
            for(size_t j=0;j<COMConstraint::getJointDOF(robot->rootLink());j++){
              tripletList.push_back(Eigen::Triplet<double>(row,idx+j,1));
            }
          }
        }
        for(size_t j=0;j<robot->numJoints();j++){
          if(jacobianColMap.find(robot->joint(j)) != jacobianColMap.end()){
            int idx = jacobianColMap[robot->joint(j)];
            for(size_t row=0;row<3;row++){
              for(size_t k=0;k<COMConstraint::getJointDOF(robot->joint(j));k++){
                tripletList.push_back(Eigen::Triplet<double>(row,idx+k,1));
              }
            }
          }
        }
      }
    }
    jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  void COMConstraint::calcCMJacobianCoef(const std::vector<cnoid::LinkPtr>& joints,//input
                                         const cnoid::BodyPtr& A_robot,//input
                                         const cnoid::BodyPtr& B_robot,//input
                                         std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap, //input
                                         Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian//output
                                         ) {
    if(A_robot != B_robot){
      for(size_t i=0;i<2; i++){
        cnoid::BodyPtr robot = (i==0) ? A_robot : B_robot;
        if(!robot) continue;
        int sign = (i==0) ? 1 : -1;

        Eigen::MatrixXd CMJ;
        cnoid::calcCMJacobian(robot,nullptr,CMJ); // [joint root]の順

        if(jacobianColMap.find(robot->rootLink()) != jacobianColMap.end()){
          int col_idx = jacobianColMap[robot->rootLink()];
          for(size_t j=0;j<3;j++){
            for(size_t k=0;k<COMConstraint::getJointDOF(robot->rootLink());k++){
              jacobian.coeffRef(j,col_idx+k) = sign * CMJ(j,robot->numJoints()+k);
            }
          }
        }
        for(size_t j=0;j<robot->numJoints();j++){
          if(jacobianColMap.find(robot->joint(j)) != jacobianColMap.end()){
            int col_idx = jacobianColMap[robot->joint(j)];
            for(size_t k=0;k<3;k++){
              for(size_t d=0;d<COMConstraint::getJointDOF(robot->joint(j));d++){
                jacobian.coeffRef(k,col_idx+d) = sign * CMJ(k,j+d);
              }
            }
          }
        }
      }
    }

  }
}
