#ifndef PRIORITIZED_INVERSE_KINEMATICS_SOLVER_H
#define PRIORITIZED_INVERSE_KINEMATICS_SOLVER_H

#include <cnoid/Body>
#include <ik_constraint/IKConstraint.h>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>
#include <prioritized_qp_osqp/prioritized_qp_osqp.h>

namespace prioritized_inverse_kinematics_solver {
  /*
    variables: 動かして良いjoint (free jointは6DOF扱い)
    ikc_list: タスクたち. vectorの前の要素の方が高優先度. 0番目の要素は必ず満たすと仮定しQPを解かない
    prevTasks: 前回のtasksを入れる. 自動的に更新される.
   */
  class IKParam {
  public:
    size_t maxIteration = 1;
    double wn = 1e-6;
    std::vector<double> wnVec; // wnVec.size() == ikc_list.size()の場合、wnの代わりにこっちを使う
    double we = 1e0;
    std::vector<double> weVec; // weVec.size() == ikc_list.size()の場合、weの代わりにこっちを使う
    int debugLevel = 0;
    double dt = 0.1;
  };
  int solveIKLoop (const std::vector<cnoid::LinkPtr>& variables,
                   const std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > >& ikc_list,
                   std::vector<std::shared_ptr<prioritized_qp_base::Task> >& prevTasks,
                   const IKParam& param = IKParam(),
                   std::function<void(std::shared_ptr<prioritized_qp_base::Task>&,int)> taskGeneratorFunc = [](std::shared_ptr<prioritized_qp_base::Task>& task, int debugLevel){
                     std::shared_ptr<prioritized_qp_osqp::Task> taskOSQP = std::dynamic_pointer_cast<prioritized_qp_osqp::Task>(task);
                     if(!taskOSQP){
                       task = std::make_shared<prioritized_qp_osqp::Task>();
                       taskOSQP = std::dynamic_pointer_cast<prioritized_qp_osqp::Task>(task);
                     }
                     taskOSQP->settings().verbose = debugLevel;
                     taskOSQP->settings().max_iter = 4000;
                     taskOSQP->settings().eps_abs = 1e-3;// 大きい方が速いが，不正確. 1e-5はかなり小さい. 1e-4は普通
                     taskOSQP->settings().eps_rel = 1e-3;// 大きい方が速いが，不正確. 1e-5はかなり小さい. 1e-4は普通
                     taskOSQP->settings().scaled_termination = true;// avoid too severe termination check
                   }
                   );

}

#endif
