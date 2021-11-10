#ifndef PRIORITIZED_INVERSE_KINEMATICS_SOLVER_H
#define PRIORITIZED_INVERSE_KINEMATICS_SOLVER_H

#include <cnoid/Body>
#include <ik_constraint/IKConstraint.h>
#include <prioritized_qp_base/PrioritizedQPBaseSolver.h>
#include <prioritized_qp/PrioritizedQPSolver.h>

namespace prioritized_inverse_kinematics_solver {
  /*
    variables: 動かして良いjoint (free jointは6DOF扱い)
    ikc_list: タスクたち. vectorの前の要素の方が高優先度. 0番目の要素は必ず満たすと仮定しQPを解かない
    prevTasks: 前回のtasksを入れる. 自動的に更新される.
   */
  int solveIKLoop (const std::vector<cnoid::LinkPtr>& variables,
                   const std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > >& ikc_list,
                   std::vector<std::shared_ptr<prioritized_qp_base::Task> >& prevTasks,
                   size_t max_iteration = 1,
                   double wn = 1e-6,
                   int debugLevel = 0,
                   double dt = 0.1,
                   std::function<void(std::shared_ptr<prioritized_qp_base::Task>&,int)> taskGeneratorFunc = [](std::shared_ptr<prioritized_qp_base::Task>& task, int debugLevel){
                     std::shared_ptr<prioritized_qp::Task> taskOSQP = std::dynamic_pointer_cast<prioritized_qp::Task>(task);
                     if(!taskOSQP){
                       task = std::make_shared<prioritized_qp::Task>();
                       taskOSQP = std::dynamic_pointer_cast<prioritized_qp::Task>(task);
                     }
                     taskOSQP->solver().settings()->setVerbosity(debugLevel);
                     taskOSQP->solver().settings()->setMaxIteration(4000);
                     taskOSQP->solver().settings()->setAbsoluteTolerance(1e-4);// 大きい方が速いが，不正確
                     taskOSQP->solver().settings()->setRelativeTolerance(1e-4);// 大きい方が速いが，不正確
                     taskOSQP->solver().settings()->setScaledTerimination(true);// avoid too severe termination check
                                                       }
                   );

}

#endif
