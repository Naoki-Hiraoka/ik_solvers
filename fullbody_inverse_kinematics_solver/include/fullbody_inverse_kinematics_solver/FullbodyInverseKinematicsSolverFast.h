#ifndef FULLBODY_INVERSE_KINEMATICS_SOLVER_FAST_H
#define FULLBODY_INVERSE_KINEMATICS_SOLVER_FAST_H

#include <iomanip>
#include <cnoid/JointPath>
#include <cnoid/Body>
#include <ik_constraint/IKConstraint.h>

namespace fik {
  /*
    jlim_avoid_weight_old: 変更される. 6+jointのサイズ,順番. 最初は0を入れよ
    dq_weight_all: 6+jointのサイズ,順番. 迷ったら1を入れよ. 0だとその関節を動かさない
   */
  int solveFullbodyIKLoopFast (const cnoid::BodyPtr& robot,
                               const std::vector<std::shared_ptr<IK::IKConstraint> >& ikc_list,
                               cnoid::VectorX& jlim_avoid_weight_old,
                               const cnoid::VectorX& dq_weight_all,
                               const size_t max_iteration = 1,
                               double wn = 1e-6,
                               int debugLevel = 0);

}

#endif
