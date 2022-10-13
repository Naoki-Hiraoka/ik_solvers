#include <prioritized_inverse_kinematics_solver/PrioritizedInverseKinematicsSolver.h>
#include <Eigen/Sparse>
#include <iostream>
#include <iomanip>
#include <set>
#include <unordered_map>
#include <cnoid/TimeMeasure>

namespace prioritized_inverse_kinematics_solver {
  bool checkIKConvergence(const std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > >& ikc_list) {
    bool converged = true;
    for ( int i=0; i<ikc_list.size(); i++ ) {
      for(size_t j=0;j<ikc_list[i].size(); j++){
        // checkConvergence()の結果をキャッシュして後に利用するものがあるので、全IKConstraintに対してcheckConvergence()を呼んでおく必要が有る.
        if (!ikc_list[i][j]->checkConvergence()) converged = false;
      }
    }
    return converged;
  }

  void solveIKOnce (const std::vector<cnoid::LinkPtr>& variables,
                    const std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > >& ikc_list,
                    std::vector<std::shared_ptr<prioritized_qp_base::Task> >& prevTasks,
                    const IKParam& param,
                    std::function<void(std::shared_ptr<prioritized_qp_base::Task>&,int)> taskGeneratorFunc) {
    // Solvability-unconcerned Inverse Kinematics by Levenberg-Marquardt Method [sugihara:RSJ2009]
    // H = J^T * We * J + Wn
    // Wn = (e^T * We * e + \bar{wn}) * Wq // Wq: modify to insert dq weight
    // Weは既にIKConstraintクラスのJ,eに含まれている

    cnoid::TimeMeasure timer;
    if(param.debugLevel>0) timer.begin();

    double dim = 0;
    for(size_t i=0;i<variables.size();i++) dim+=IK::IKConstraint::getJointDOF(variables[i]);

    if(prevTasks.size() != ikc_list.size()) {
      prevTasks.clear();
      prevTasks.resize(ikc_list.size(),nullptr);
    }
    for(size_t i=0;i<ikc_list.size();i++){
      taskGeneratorFunc(prevTasks[i],param.debugLevel);

      if(i!=0) prevTasks[i]->toSolve() = true;
      else prevTasks[i]->toSolve() = false;

      int num_eqs = 0;
      int num_ineqs = 0;
      std::vector<std::reference_wrapper<const Eigen::VectorXd> > errors;errors.reserve(ikc_list[i].size());
      std::vector<std::reference_wrapper<const Eigen::SparseMatrix<double,Eigen::RowMajor> > > jacobians;jacobians.reserve(ikc_list[i].size());
      std::vector<std::reference_wrapper <const Eigen::VectorXd> > minineqs;minineqs.reserve(ikc_list[i].size());
      std::vector<std::reference_wrapper<const Eigen::VectorXd> > maxineqs;maxineqs.reserve(ikc_list[i].size());
      std::vector<std::reference_wrapper<const Eigen::SparseMatrix<double,Eigen::RowMajor> > > jacobianineqs;jacobianineqs.reserve(ikc_list[i].size());

      for(size_t j=0; j<ikc_list[i].size(); j++){
        errors.emplace_back(ikc_list[i][j]->calc_error());
        jacobians.emplace_back(ikc_list[i][j]->calc_jacobian(variables));
        jacobianineqs.emplace_back(ikc_list[i][j]->calc_jacobianineq(variables));
        minineqs.emplace_back(ikc_list[i][j]->calc_minineq());
        maxineqs.emplace_back(ikc_list[i][j]->calc_maxineq());

        num_eqs += errors[j].get().rows();
        num_ineqs += minineqs[j].get().rows();
      }

      prevTasks[i]->A().resize(num_eqs, dim);
      prevTasks[i]->b().resize(num_eqs);
      prevTasks[i]->C().resize(num_ineqs, dim);
      prevTasks[i]->dl().resize(num_ineqs);
      prevTasks[i]->du().resize(num_ineqs);
      prevTasks[i]->wa() = cnoid::VectorXd::Ones(num_eqs);
      prevTasks[i]->wc() = cnoid::VectorXd::Ones(num_ineqs);

      int idx_eq = 0;
      int idx_ineq = 0;
      for(size_t j=0;j<ikc_list[i].size(); j++){
        prevTasks[i]->A().middleRows(idx_eq,errors[j].get().rows()) = jacobians[j].get();
        prevTasks[i]->b().segment(idx_eq,errors[j].get().rows()) = - errors[j].get();
        idx_eq += errors[j].get().rows();

        prevTasks[i]->C().middleRows(idx_ineq,minineqs[j].get().rows()) = jacobianineqs[j].get();
        prevTasks[i]->dl().segment(idx_ineq,minineqs[j].get().rows()) = minineqs[j].get();
        prevTasks[i]->du().segment(idx_ineq,minineqs[j].get().rows()) = maxineqs[j].get();
        idx_ineq += minineqs[j].get().rows();
      }

      double sumError = 0;
      sumError += prevTasks[i]->b().squaredNorm();
      for(size_t j=0;j<prevTasks[i]->dl().size(); j++) {
        if(prevTasks[i]->dl()[j]>0) sumError += std::pow(prevTasks[i]->dl()[j],2);
        if(prevTasks[i]->du()[j]<0) sumError += std::pow(prevTasks[i]->du()[j],2);
      }
      prevTasks[i]->w() = cnoid::VectorXd::Ones(dim) * (sumError * ((param.weVec.size()==ikc_list.size())?param.weVec[i]:param.we)+ ((param.wnVec.size()==ikc_list.size())?param.wnVec[i]:param.wn));

      if(param.debugLevel>0) prevTasks[i]->name() = std::string("Task") + std::to_string(i);
    }

    // solve
    cnoid::VectorX result;
    if(!prioritized_qp_base::solve(prevTasks, result, param.debugLevel)){
      std::cerr <<"[PrioritizedIK] prioritized_qp_base::solve failed" << std::endl;
      return;
    }

    if (!result.allFinite()) {
      std::cerr <<"[PrioritizedIK] ERROR nan/inf is found" << std::endl;
      return;
    }

    size_t idx = 0;
    for(size_t i=0;i<variables.size();i++){
      if(variables[i]->isRevoluteJoint() || variables[i]->isPrismaticJoint()){
        // update joint angles
        variables[i]->q() += result[idx];
        if(variables[i]->q() > variables[i]->q_upper()) variables[i]->q() = variables[i]->q_upper();
        if(variables[i]->q() < variables[i]->q_lower()) variables[i]->q() = variables[i]->q_lower();
      }else if(variables[i]->isFreeJoint()) {
        // update rootlink pos rot
        variables[i]->p() += result.segment<3>(idx);
        if(result.segment<3>(idx+3).norm() != 0){
          variables[i]->R() = cnoid::Matrix3(cnoid::AngleAxis(result.segment<3>(idx+3).norm(), cnoid::Vector3(result.segment<3>(idx+3).normalized())) * cnoid::AngleAxis(variables[i]->R()));
          // 単純に3x3行列の空間でRを操作していると、だんだん数値誤差によって回転行列でなくなってしまう
          //const cnoid::Matrix3 dR = cnoid::Matrix3(cnoid::AngleAxis(result.segment<3>(idx+3).norm(), cnoid::Vector3(result.segment<3>(idx+3).normalized())));
          //variables[i]->R() = (dR * variables[i]->R()).eval();
        }
        if (!variables[i]->R().isUnitary()) {
          std::cerr <<"[PrioritizedIK] WARN robot->rootLink()->R is not Unitary, something wrong !" << std::endl;
        }
      }

      idx += IK::IKConstraint::getJointDOF(variables[i]);
    }

    if(param.debugLevel>0) {
      double time = timer.measure();
      std::cerr << "[PrioritizedIK] solveIKOnce time: " << time << "[s]" << std::endl;
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

  int solveIKLoop (const std::vector<cnoid::LinkPtr>& variables,
                   const std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > >& ikc_list,
                   std::vector<std::shared_ptr<prioritized_qp_base::Task> >& prevTasks,
                   const IKParam& param,
                   std::function<void(std::shared_ptr<prioritized_qp_base::Task>&,int)> taskGeneratorFunc) {
    std::set<cnoid::BodyPtr> bodies;
    for(size_t i=0;i<variables.size();i++){
      if(variables[i]->body()) bodies.insert(variables[i]->body());
    }

    std::unordered_map<cnoid::LinkPtr, InitialJointState> initialJointStateMap;
    for(size_t i=0;i<variables.size();i++){
      if(variables[i]->isFreeJoint()) initialJointStateMap[variables[i]] = InitialJointState(variables[i]->T());
      else if(variables[i]->isRotationalJoint() || variables[i]->isPrismaticJoint()) initialJointStateMap[variables[i]] = InitialJointState(variables[i]->q());
      else initialJointStateMap[variables[i]] = InitialJointState();
    }

    size_t loop = 0;
    while (loop < param.maxIteration) {
      for(size_t i=0;i<variables.size();i++){
        if(variables[i]->isFreeJoint()) {
          cnoid::Position& initialT = initialJointStateMap[variables[i]].T;
          variables[i]->v() = (variables[i]->p() - initialT.translation()) / param.dt;
          cnoid::AngleAxis angleAxis = cnoid::AngleAxis(variables[i]->R() * initialT.linear().transpose());
          variables[i]->w() = angleAxis.angle()*angleAxis.axis() / param.dt;
        }
        else if(variables[i]->isRotationalJoint() || variables[i]->isPrismaticJoint()) {
          double initialq = initialJointStateMap[variables[i]].q;
          variables[i]->dq() = (variables[i]->q() - initialq) / param.dt;
        }
      }
      for(std::set<cnoid::BodyPtr>::iterator it=bodies.begin(); it != bodies.end(); it++){
        (*it)->calcForwardKinematics(true);
        (*it)->calcCenterOfMass();
      }

      if (checkIKConvergence(ikc_list)) return loop;
      solveIKOnce(variables, ikc_list, prevTasks, param, taskGeneratorFunc);
      ++loop;
    }

    for(size_t i=0;i<variables.size();i++){
      if(variables[i]->isFreeJoint()) {
        cnoid::Position& initialT = initialJointStateMap[variables[i]].T;
        variables[i]->v() = (variables[i]->p() - initialT.translation()) / param.dt;
        cnoid::AngleAxis angleAxis = cnoid::AngleAxis(variables[i]->R() * initialT.linear().transpose());
        variables[i]->w() = angleAxis.angle()*angleAxis.axis() / param.dt;
      }
      else if(variables[i]->isRotationalJoint() || variables[i]->isPrismaticJoint()) {
        double initialq = initialJointStateMap[variables[i]].q;
        variables[i]->dq() = (variables[i]->q() - initialq) / param.dt;
      }
    }
    for(std::set<cnoid::BodyPtr>::iterator it=bodies.begin(); it != bodies.end(); it++){
      (*it)->calcForwardKinematics(true);
      (*it)->calcCenterOfMass();
    }

    return loop;
  }

}
