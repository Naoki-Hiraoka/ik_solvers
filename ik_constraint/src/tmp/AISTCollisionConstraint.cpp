#include "AISTCollisionConstraint.h"

namespace IK{

  AISTCollisionConstraint::AISTCollisionConstraint(cnoid::Link* _A_link, cnoid::Link* _B_link):
    CollisionConstraint(_A_link,_B_link),
    _collisionDetector(new cnoid::AISTCollisionDetector)
  {
    A_handle = *_collisionDetector->addGeometry(A_link->collisionShape());
    B_handle = *_collisionDetector->addGeometry(B_link->collisionShape());
    _collisionDetector->makeReady();
  }

  double AISTCollisionConstraint::detectDistance(cnoid::Vector3& A_v, cnoid::Vector3& B_v){
    _collisionDetector->updatePosition(A_handle,A_link->position());
    _collisionDetector->updatePosition(B_handle,B_link->position());
    return _collisionDetector->detectDistance(A_handle,B_handle,A_v,B_v);
  }

}
