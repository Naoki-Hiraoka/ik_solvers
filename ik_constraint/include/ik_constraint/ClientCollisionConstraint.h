#ifndef IK_CONSTRAINT_CLIENTCOLLISIONCONSTRAINT_H
#define IK_CONSTRAINT_CLIENTCOLLISIONCONSTRAINT_H

#include <ik_constraint/CollisionConstraint.h>

namespace IK{
  class ClientCollisionConstraint : public CollisionConstraint
  {
  public:
    ClientCollisionConstraint();

    // A_linkとB_linkの干渉を回避する
    //  最近傍点と、方向ベクトルはgiven. 距離のみ計算
    //  direction: B->A

    const cnoid::Vector3& A_localp() const { return A_localp_;}
    cnoid::Vector3& A_localp() { return A_localp_;}
    const cnoid::Vector3& B_localp() const { return B_localp_;}
    cnoid::Vector3& B_localp() { return B_localp_;}
    const cnoid::Vector3& direction() const { return direction_;}
    cnoid::Vector3& direction() { return direction_;}

  protected:
    //A_v, B_vはlocal系
    virtual bool computeDistance(const cnoid::LinkPtr A_link, const cnoid::LinkPtr B_link, double& distance, cnoid::Vector3& direction/*B->A*/, cnoid::Vector3& A_v, cnoid::Vector3& B_v) override;

    cnoid::Vector3 A_localp_ = cnoid::Vector3::Zero();
    cnoid::Vector3 B_localp_ = cnoid::Vector3::Zero();
    cnoid::Vector3 direction_ = cnoid::Vector3::UnitX();
  };


}

#endif
