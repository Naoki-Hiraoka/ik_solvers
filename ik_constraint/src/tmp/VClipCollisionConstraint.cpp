#include "VClipCollisionConstraint.h"

namespace IK{

  VClipCollisionConstraint::VClipCollisionConstraint(cnoid::Link* _A_link, cnoid::Link* _B_link, Vclip::Polyhedron* A_vclip_model, Vclip::Polyhedron* B_vclip_model):
    CollisionConstraint(_A_link,_B_link),
    vcliplinkpair(std::make_shared<VclipLinkPair>(A_link,A_vclip_model,B_link,B_vclip_model))
  {
  }

  double VClipCollisionConstraint::detectDistance(cnoid::Vector3& A_v, cnoid::Vector3& B_v){
    double d = this->vcliplinkpair->computeDistance(this->a_v,this->b_v);
    A_v[0]=a_v[0];A_v[1]=a_v[1];A_v[2]=a_v[2];
    B_v[0]=b_v[0];B_v[1]=b_v[1];B_v[2]=b_v[2];
    return d;
  }

}
