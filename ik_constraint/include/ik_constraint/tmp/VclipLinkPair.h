#include <cnoid/Body>
#include "vclip.h"

namespace IK{
  class VclipLinkPair {
  public:
    VclipLinkPair(const cnoid::Link* link0, Vclip::Polyhedron* pqp_model0, const cnoid::Link* link1, Vclip::Polyhedron* pqp_model1, double tolerance=0);
    ~VclipLinkPair();
    bool checkCollision();
    double computeDistance(double *q1, double *q2);//q1,q2はworld系
    const cnoid::Link* link(int index) { return links_[index]; }
    double getTolerance() { return tolerance_; }
    void setTolerance(double t) { tolerance_ = t; }

  private:
    const cnoid::Link *links_[2];
    Vclip::Polyhedron *Vclip_Model1, *Vclip_Model2;
    Vclip::FeaturePair Feature_Pair;
    double tolerance_;
  };

};

