#include <fullbody_inverse_kinematics_solver/FullbodyInverseKinematicsSolverFast.h>
#include <Eigen/Sparse>
#include <iostream>
#include <iomanip>
#include <Eigen/SparseCholesky>

namespace fik {
  inline Eigen::SparseMatrix<double,Eigen::RowMajor> toSelectionMatrix(const Eigen::VectorXd& in)
  {
    Eigen::SparseMatrix<double,Eigen::RowMajor> ret((in.array() > 0.0).count(), in.size());
    if (ret.rows() != 0 && ret.cols() != 0) {
      for (int row=0, col=0; col< in.size(); ++col) {
        if(in(col) > 0.0){
          ret.insert(row, col) = 1;
          row++;
        }
      }
    }
    return ret;
  }

  inline Eigen::SparseMatrix<double,Eigen::RowMajor> toDiagonalMatrix(const Eigen::VectorXd& in)
  {
    Eigen::SparseMatrix<double,Eigen::RowMajor> ret(in.size(), in.size());
    for (int i=0; i< in.size(); ++i) {
      ret.insert(i, i) = in[i];
    }
    return ret;
  }

  bool checkIKConvergence(const std::vector<std::shared_ptr<IK::IKConstraint> >& ikc_list) {
    bool converged = true;
    for ( int i=0; i<ikc_list.size(); i++ ) {
      // checkConvergence()の結果をキャッシュして後に利用するものがあるので、全IKConstraintに対してcheckConvergence()を呼んでおく必要が有る.
      if (!ikc_list[i]->checkConvergence()) converged = false;
    }
    return converged;
  }

  void solveFullbodyIKOnceFast (const cnoid::BodyPtr& robot, const std::vector<std::shared_ptr<IK::IKConstraint> >& ikc_list, const cnoid::VectorX& dq_weight_all, cnoid::VectorX& jlim_avoid_weight_old, double wn, int debugLevel, double we, bool enable_jlim_avoid) {
    // Solvability-unconcerned Inverse Kinematics by Levenberg-Marquardt Method [sugihara:RSJ2009]
    // q = q + S_q^T * H^(-1) * g // S_q: selection matrix
    // g = J^T * We * e
    // H = J^T * We * J + Wn
    // Wn = (e^T * We * e * \bar{we} + \bar{wn}) * Wq // Wq: modify to insert dq weight
    // Weは既にIKConstraintクラスのJ,eに含まれている
    const Eigen::SparseMatrix<double,Eigen::RowMajor> S_q = toSelectionMatrix(dq_weight_all);
    const size_t VALID_Q_NUM = S_q.rows();
    if (VALID_Q_NUM == 0) return;

    Eigen::SparseMatrix<double,Eigen::RowMajor> H(VALID_Q_NUM,VALID_Q_NUM);
    Eigen::VectorXd g = Eigen::VectorXd::Zero(VALID_Q_NUM);
    double sumError = 0.0;

    std::vector<cnoid::LinkPtr> joints;
    joints.push_back(robot->rootLink());
    for(int j=0;j<robot->numJoints();j++){
      joints.push_back(robot->joint(j));
    }
    for(size_t i=0;i<ikc_list.size();i++){
      const Eigen::VectorXd& error = ikc_list[i]->calc_error();
      const Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian = ikc_list[i]->calc_jacobian(joints) * S_q.transpose();

      H += Eigen::SparseMatrix<double,Eigen::RowMajor>(jacobian.transpose() * jacobian);
      g -= jacobian.transpose() * error; // minus
      sumError += error.squaredNorm();
    }

    // joint limit avoidance
    // T. F. Chang and R.-V. Dubey: "A weighted least-norm solution based
    // scheme for avoiding joint limits for redundant manipulators", in IEEE
    // Trans. On Robotics and Automation, 11((2):286-292, April 1995.
    if(jlim_avoid_weight_old.rows() != 6+robot->numJoints()) jlim_avoid_weight_old = Eigen::VectorXd::Zero(6+robot->numJoints());
    Eigen::VectorXd dq_weight_all_jlim = Eigen::VectorXd::Ones(6+robot->numJoints());
    if(enable_jlim_avoid){
      for ( int j = 0; j < robot->numJoints() ; j++ ) {
        double jang = robot->joint(j)->q();
        double jmax = robot->joint(j)->q_upper();
        double jmin = robot->joint(j)->q_lower();
        const double e = 0.017453;
        double jlim_avoid_weight;
        if ( jmax - jmin > 2 * e) {
          if ( jang > jmax - e ) jang = jmax - e;
          if ( jang < jmin + e ) jang = jmin + e;
          jlim_avoid_weight = std::fabs( (std::pow((jmax - jmin),2) * (( 2 * jang) - jmax - jmin)) / (4 * std::pow((jmax - jang),2) * std::pow((jang - jmin),2)) );
          if (std::isnan(jlim_avoid_weight)) jlim_avoid_weight = 0;
        } else {
          jlim_avoid_weight = std::numeric_limits<double>::max();
        }
        if (( jlim_avoid_weight - jlim_avoid_weight_old(6+j) ) >= 0 ) { // add weight only if q approaching to the limit
          dq_weight_all_jlim(6+j) += jlim_avoid_weight;
        }
        jlim_avoid_weight_old(6+j) = jlim_avoid_weight;
      }
    }
    const Eigen::VectorXd dq_weight_all_final = S_q * static_cast<Eigen::VectorXd>(dq_weight_all.array() * dq_weight_all_jlim.array());
    Eigen::SparseMatrix<double,Eigen::RowMajor> Wn = (sumError * we + wn) * toDiagonalMatrix(dq_weight_all_final);
    H += Wn;

    // solve
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    Eigen::VectorXd dq_all = solver.compute(Eigen::SparseMatrix<double>(H)).solve(g); // dq_all = H.inverse() * g; is slow

#define dbg(var)  std::cout << #var"= " << (var) << std::endl
#define dbgn(var) std::cout << #var"= " << std::endl <<(var) << std::endl
#define dbgv(var) std::cout << #var"= " << (var.transpose()) << std::endl
    if (debugLevel >= 1) {
      dbgn(H);
      dbg(H.rows());
      dbg(H.cols());
      dbgn(Wn);
      dbg(Wn.rows());
      dbg(Wn.cols());
      dbgv(g);
      dbg(g.rows());
      dbg(g.cols());
      dbg(sumError);
      dbgv(dq_all);
      std::cout<<std::endl;
    }

    if (!dq_all.allFinite()) {
      std::cerr <<"[FullbodyIK] ERROR nan/inf is found" << std::endl;
      return;
    }
    dq_all = (S_q.transpose() * dq_all).eval();

    // update rootlink pos rot
    robot->rootLink()->p() += dq_all.head<3>();
    if(dq_all.head<6>().tail<3>().norm() != 0){
      const cnoid::AngleAxis dR = cnoid::AngleAxis(dq_all.head<6>().tail<3>().norm(), cnoid::Vector3(dq_all.head<6>().tail<3>().normalized()));
      robot->rootLink()->R() = cnoid::Matrix3(dR * cnoid::AngleAxis(robot->rootLink()->R()));
      // 単純に3x3行列の空間でRを操作していると、だんだん数値誤差によって回転行列でなくなってしまう
      // robot->rootLink()->R() = (dR * robot->rootLink()->R()).eval();
    }
    if (!robot->rootLink()->R().isUnitary()) {
      std::cerr <<"[FullbodyIK] WARN robot->rootLink()->R is not Unitary, something wrong !" << std::endl;
    }
    // update joint angles
    for (size_t i = 0; i < robot->numJoints(); ++i) {
      if(dq_weight_all[6+i] > 0.0) {// dq_weight_all == 0.0の関節は更新しない
        robot->joint(i)->q() += dq_all(6+i);
        if(robot->joint(i)->q() > robot->joint(i)->q_upper()) robot->joint(i)->q() = robot->joint(i)->q_upper();
        if(robot->joint(i)->q() < robot->joint(i)->q_lower()) robot->joint(i)->q() = robot->joint(i)->q_lower();
      }
    }
  }

  class InitialJointState {
  public:
    InitialJointState() {}
    InitialJointState(const cnoid::Position& T_): T(T_) {}
    InitialJointState(double q_): q(q_) {}
    cnoid::Position T;
    double q;
  };

  int solveFullbodyIKLoopFast (const cnoid::BodyPtr& robot, const std::vector<std::shared_ptr<IK::IKConstraint> >& ikc_list, cnoid::VectorX& jlim_avoid_weight_old, const cnoid::VectorX& dq_weight_all, const size_t max_iteration, double wn, int debugLevel, double dt, double we, bool enable_jlim_avoid) {

    // 開始時の関節角度を記憶. 後で速度の計算に使う
    std::vector<InitialJointState> initialJointState(robot->numJoints()+1);
    for(size_t i=0;i<robot->numJoints();i++){
      initialJointState[i].q = robot->joint(i)->q();
    }
    initialJointState[robot->numJoints()].T = robot->rootLink()->T();

    size_t loop = 0;
    while (loop < max_iteration) {
      // 速度を更新. dq_weight_all == 0.0の関節は更新しない
      for(int i=0;i<3;i++) if(dq_weight_all[i] > 0.0) robot->rootLink()->v()[i] = (robot->rootLink()->p()[i] - initialJointState[robot->numJoints()].T.translation()[i]) / dt;
      if((dq_weight_all.head<6>().tail<3>().array() > 0.0).count() > 0.0) {
        cnoid::AngleAxis angleAxis = cnoid::AngleAxis(robot->rootLink()->R() * initialJointState[robot->numJoints()].T.linear().transpose());
        robot->rootLink()->w() = angleAxis.angle()*angleAxis.axis() / dt;
      }
      for(size_t i=0;i<robot->numJoints();i++) {
        if(dq_weight_all[6+i] > 0.0) {
          robot->joint(i)->dq() = (robot->joint(i)->q() - initialJointState[i].q) / dt;
          // check dq limit
          if(robot->joint(i)->dq() > robot->joint(i)->dq_upper()) {
            robot->joint(i)->q() -= (robot->joint(i)->dq() - robot->joint(i)->dq_upper()) * dt;
            robot->joint(i)->dq() = robot->joint(i)->dq_upper();
          }else if(robot->joint(i)->dq() < robot->joint(i)->dq_lower()) {
            robot->joint(i)->q() += (robot->joint(i)->dq_lower() - robot->joint(i)->dq()) * dt;
            robot->joint(i)->dq() = robot->joint(i)->dq_lower();
          }
        }
      }

      robot->calcForwardKinematics(true);
      robot->calcCenterOfMass();
      if (checkIKConvergence(ikc_list)) return loop;
      solveFullbodyIKOnceFast(robot, ikc_list, dq_weight_all, jlim_avoid_weight_old, wn, debugLevel, we, enable_jlim_avoid);
      ++loop;
    }

    // 速度を更新. dq_weight_all == 0.0の関節は更新しない
    for(int i=0;i<3;i++) if(dq_weight_all[i] > 0.0) robot->rootLink()->v()[i] = (robot->rootLink()->p()[i] - initialJointState[robot->numJoints()].T.translation()[i]) / dt;
    if((dq_weight_all.head<6>().tail<3>().array() > 0.0).count() > 0.0) {
      cnoid::AngleAxis angleAxis = cnoid::AngleAxis(robot->rootLink()->R() * initialJointState[robot->numJoints()].T.linear().transpose());
      robot->rootLink()->w() = angleAxis.angle()*angleAxis.axis() / dt;
    }
    for(size_t i=0;i<robot->numJoints();i++) {
      if(dq_weight_all[6+i] > 0.0)
        robot->joint(i)->dq() = (robot->joint(i)->q() - initialJointState[i].q) / dt;
      // check dq limit
      if(robot->joint(i)->dq() > robot->joint(i)->dq_upper()) {
        robot->joint(i)->q() -= (robot->joint(i)->dq() - robot->joint(i)->dq_upper()) * dt;
        robot->joint(i)->dq() = robot->joint(i)->dq_upper();
      }else if(robot->joint(i)->dq() < robot->joint(i)->dq_lower()) {
        robot->joint(i)->q() += (robot->joint(i)->dq_lower() - robot->joint(i)->dq()) * dt;
        robot->joint(i)->dq() = robot->joint(i)->dq_lower();
      }
    }

    robot->calcForwardKinematics(true);
    robot->calcCenterOfMass();

    return loop;
  }

}
