#include <ik_constraint/ClientCollisionConstraint.h>

namespace IK {
  ClientCollisionConstraint::ClientCollisionConstraint()
  {
  }


  bool ClientCollisionConstraint::computeDistance(const cnoid::LinkPtr A_link, const cnoid::LinkPtr B_link, double& distance, cnoid::Vector3& direction/*B->A*/, cnoid::Vector3& A_v, cnoid::Vector3& B_v) {
    if (this->direction_.norm() == 0) return false;
    direction = this->direction_.normalized();
    A_v = this->A_localp_;
    B_v = this->B_localp_;

    cnoid::Vector3 A_v_global = A_link ? A_link->T()*A_v : A_v;
    cnoid::Vector3 B_v_global = B_link ? B_link->T()*B_v : B_v;

    distance = direction_.dot(A_v_global - B_v_global);

    return true;
  }
}
