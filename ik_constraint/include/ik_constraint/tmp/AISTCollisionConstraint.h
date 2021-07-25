#ifndef AISTCOLLISIONCONSTRAINT_H
#define AISTCOLLISIONCONSTRAINT_H

#include "CollisionConstraint.h"
#include <cnoid/AISTCollisionDetector>

namespace IK{
  // VClipの方が速い
  class AISTCollisionConstraint : public CollisionConstraint
  {
  public:
    AISTCollisionConstraint(cnoid::Link* A_link, cnoid::Link* B_link);

    double detectDistance(cnoid::Vector3& A_v, cnoid::Vector3& B_v) override;
  private:
    cnoid::AISTCollisionDetectorPtr _collisionDetector;
    cnoid::CollisionDetector::GeometryHandle A_handle;
    cnoid::CollisionDetector::GeometryHandle B_handle;
  };
}


#endif
