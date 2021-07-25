#ifndef VCLIPCOLLISIONCONSTRAINT_H
#define VCLIPCOLLISIONCONSTRAINT_H

#include "CollisionConstraint.h"
#include "VclipLinkPair.h"

namespace IK{
  class VClipCollisionConstraint : public CollisionConstraint
  {
  public:
    VClipCollisionConstraint(cnoid::Link* A_link, cnoid::Link* B_link, Vclip::Polyhedron* A_vclip_model, Vclip::Polyhedron* B_vclip_model);

    double detectDistance(cnoid::Vector3& A_v, cnoid::Vector3& B_v) override;
  private:
    std::shared_ptr<VclipLinkPair> vcliplinkpair;
    double a_v[3],b_v[3];
  };
}


#endif
