#ifndef PRIORITIZED_INVERSE_KINEMATICS_SOLVER_H
#define PRIORITIZED_INVERSE_KINEMATICS_SOLVER_H

#include <cnoid/Body>
#include <ik_constraint/IKConstraint.h>
#include <prioritized_qp/PrioritizedQPSolver.h>

namespace prioritized_inverse_kinematics_solver {
  /*
    variables: 動かして良いjoint (free jointは6DOF扱い)
    ikc_list: タスクたち. vectorの前の要素の方が高優先度.
    prevTasks: 前回のtasksを入れる. 自動的に更新される.
   */
  int solveIKLoop (const std::vector<cnoid::LinkPtr>& variables,
                   const std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > >& ikc_list,
                   std::vector<std::shared_ptr<prioritized_qp::Task> >& prevTasks,
                   size_t max_iteration = 1,
                   double wn = 1e-6,
                   int debugLevel = 0,
                   double dt = 0.1);

}

#endif
