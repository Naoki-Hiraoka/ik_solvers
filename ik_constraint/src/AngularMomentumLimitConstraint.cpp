#include <ik_constraint/AngularMomentumLimitConstraint.h>
#include <ik_constraint/AngularMomentumConstraint.h> // for common functions
#include <cnoid/Jacobian>

namespace IK{
  bool AngularMomentumLimitConstraint::checkConvergence () {
    if(this->robot_) {
      Eigen::MatrixXd AMJ;
      cnoid18::calcAngularMomentumJacobian(this->robot_,nullptr,AMJ); // [joint root]の順. comまわり
      cnoid::VectorX dq(this->robot_->numJoints()+6);
      for(int i=0;i<this->robot_->numJoints();i++) dq[i] = this->robot_->joint(i)->dq();
      dq.segment<3>(this->robot_->numJoints()) = this->robot_->rootLink()->v();
      dq.tail<3>() = this->robot_->rootLink()->w();
      cnoid::Vector3 dAM = this->eval_R_.transpose() * AMJ * dq; //eval_R系
      cnoid::Vector3 upper = (this->maxAngularMomentum_ - dAM) * this->dt_; //eval_R系  [kg m^2]
      cnoid::Vector3 lower = (this->minAngularMomentum_ - dAM) * this->dt_; //eval_R系  [kg m^2]

      cnoid::Matrix3 I = AMJ.block<3,3>(0,this->robot_->numJoints()+3); //world系
      cnoid::Matrix3 I_evalR = this->eval_R_.transpose() * I * this->eval_R_;
      cnoid::Matrix3 I_evalR_inv;
      I_evalR_inv <<
        1.0/I_evalR(0,0), 0.0, 0.0,
        0.0, 1.0/I_evalR(1,1), 0.0,
        0.0, 0.0, 1.0/I_evalR(2,2);
      cnoid::Vector3 upper_scaled = I_evalR_inv * upper; //eval_R系  [rad]
      cnoid::Vector3 lower_scaled = I_evalR_inv * lower; //eval_R系  [rad]

      if(this->minineq_.rows() != 3) this->minineq_ = Eigen::VectorXd(3);
      for(int i=0;i<3;i++) this->minineq_[i] = std::min(lower_scaled[i], this->maxError_[i]) * this->weight_[i];
      if(this->maxineq_.rows() != 3) this->maxineq_ = Eigen::VectorXd(3);
      for(int i=0;i<3;i++) this->maxineq_[i] = std::max(upper_scaled[i], -this->maxError_[i]) * this->weight_[i];

      bool converged = true;
      for(size_t i=0; i<3; i++){
        if(this->weight_[i]>0.0) {
          if(upper[i] < -this->precision_[i]) converged = false;
          if(lower[i] >  this->precision_[i]) converged = false;
        }
      }

      if(this->debuglevel_>=1){
        std::cerr << "AngularMomentumLimitConstraint" << std::endl;
        std::cerr << "dAM" << std::endl;
        std::cerr << dAM << std::endl;
        std::cerr << "I_evalR" << std::endl;
        std::cerr << I_evalR << std::endl;
        std::cerr << "upper" << std::endl;
        std::cerr << upper << std::endl;
        std::cerr << "lower" << std::endl;
        std::cerr << lower << std::endl;
        std::cerr << "upper_scaled" << std::endl;
        std::cerr << upper_scaled << std::endl;
        std::cerr << "lower_scaled" << std::endl;
        std::cerr << lower_scaled << std::endl;
      }

      return converged;
    }

    return true;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& AngularMomentumLimitConstraint::calc_jacobianineq (const std::vector<cnoid::LinkPtr>& joints) {
    // 行列の初期化. 前回とcol形状が変わっていないなら再利用
    if(!this->is_joints_same(joints,this->jacobianineq_joints_)
       || this->robot_ != this->jacobianineq_robot_){
      this->jacobianineq_joints_ = joints;
      this->jacobianineq_robot_ = this->robot_;

      AngularMomentumLimitConstraint::calcAngularMomentumJacobianShape(this->jacobianineq_joints_,
                                                                  this->jacobianineq_robot_,
                                                                  nullptr,
                                                                  this->jacobianineq_full_,
                                                                  this->jacobianineqColMap_);
    }

    AngularMomentumLimitConstraint::calcAngularMomentumJacobianCoef(this->jacobianineq_joints_,
                                                               this->jacobianineq_robot_,
                                                               nullptr,
                                                               this->jacobianineqColMap_,
                                                               this->jacobianineq_full_);

    Eigen::MatrixXd AMJ;
    cnoid18::calcAngularMomentumJacobian(this->robot_,nullptr,AMJ); // [joint root]の順. comまわり
    cnoid::Matrix3 I = AMJ.block<3,3>(0,this->robot_->numJoints()+3); //world系
    Eigen::SparseMatrix<double,Eigen::RowMajor> I_inv(3,3);
    for(int i=0;i<3;i++) I_inv.insert(i,i) = 1.0/I(i,i);

    Eigen::SparseMatrix<double,Eigen::RowMajor> eval_R(3,3);
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) eval_R.insert(i,j) = this->eval_R_(i,j);

    this->jacobianineq_ = eval_R.transpose() * I_inv * this->jacobianineq_full_;
    for(size_t i=0;i<3;i++) this->jacobianineq_.row(i) *= this->weight_[i];

    if(this->debuglevel_>=1){
      std::cerr << "AngularMomentumLimitConstraint" << std::endl;
      std::cerr << "jacobiaineq" << std::endl;
      std::cerr << this->jacobianineq_ << std::endl;
    }

    return this->jacobianineq_;
  }

  const Eigen::VectorXd& AngularMomentumLimitConstraint::calc_minineq () {
    if(this->debuglevel_>=1){
      std::cerr << "AngularMomentumLimitConstraint" << std::endl;
      std::cerr << "minineq" << std::endl;
      std::cerr << this->minineq_ << std::endl;
    }

    return this->minineq_;
  }

  const Eigen::VectorXd& AngularMomentumLimitConstraint::calc_maxineq () {
    if(this->debuglevel_>=1){
      std::cerr << "AngularMomentumLimitConstraint" << std::endl;
      std::cerr << "maxineq" << std::endl;
      std::cerr << this->maxineq_ << std::endl;
    }

    return this->maxineq_;
  }

  void AngularMomentumLimitConstraint::calcAngularMomentumJacobianShape(const std::vector<cnoid::LinkPtr>& joints,
                                                                   const cnoid::BodyPtr& A_robot,
                                                                   const cnoid::BodyPtr& B_robot,
                                                                   Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian,
                                                                   std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap){
    jacobianColMap.clear();
    int cols = 0;
    for(size_t i=0;i<joints.size();i++){
      jacobianColMap[joints[i]] = cols;
      cols += getJointDOF(joints[i]);
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
            if(robot->rootLink()->isFreeJoint()){
              for(size_t j=3;j<getJointDOF(robot->rootLink());j++){//並進は無視
                tripletList.push_back(Eigen::Triplet<double>(row,idx+j,1));
              }
            }else{
              for(size_t j=0;j<getJointDOF(robot->rootLink());j++){
                tripletList.push_back(Eigen::Triplet<double>(row,idx+j,1));
              }
            }
          }
        }
        for(size_t j=0;j<robot->numJoints();j++){
          if(jacobianColMap.find(robot->joint(j)) != jacobianColMap.end()){
            int idx = jacobianColMap[robot->joint(j)];
            for(size_t row=0;row<3;row++){
              for(size_t k=0;k<getJointDOF(robot->joint(j));k++){
                tripletList.push_back(Eigen::Triplet<double>(row,idx+k,1));
              }
            }
          }
        }
      }
    }
    jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  void AngularMomentumLimitConstraint::calcAngularMomentumJacobianCoef(const std::vector<cnoid::LinkPtr>& joints,//input
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

        Eigen::MatrixXd AMJ;
        cnoid18::calcAngularMomentumJacobian(robot,nullptr,AMJ); // [joint root]の順. comまわり

        if(jacobianColMap.find(robot->rootLink()) != jacobianColMap.end()){
          int col_idx = jacobianColMap[robot->rootLink()];
          for(size_t j=0;j<3;j++){
            if(robot->rootLink()->isFreeJoint()){
              for(size_t k=3;k<getJointDOF(robot->rootLink());k++){ // 並進は無視
                jacobian.coeffRef(j,col_idx+k) = sign * AMJ(j,robot->numJoints()+k);
              }
            }else{
              for(size_t k=0;k<getJointDOF(robot->rootLink());k++){
                jacobian.coeffRef(j,col_idx+k) = sign * AMJ(j,robot->numJoints()+k);
              }
            }
          }
        }
        for(size_t j=0;j<robot->numJoints();j++){
          if(jacobianColMap.find(robot->joint(j)) != jacobianColMap.end()){
            int col_idx = jacobianColMap[robot->joint(j)];
            for(size_t k=0;k<3;k++){
              for(size_t d=0;d<getJointDOF(robot->joint(j));d++){
                jacobian.coeffRef(k,col_idx+d) = sign * AMJ(k,j+d);
              }
            }
          }
        }
      }
    }
  }
}
