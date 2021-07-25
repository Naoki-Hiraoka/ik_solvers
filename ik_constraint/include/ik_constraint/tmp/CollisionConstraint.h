#ifndef COLLISIONCONSTRAINT_H
#define COLLISIONCONSTRAINT_H

#include "IKConstraint.h"
#include "../../RobotConfig/RobotConfig.h"
#include <cnoid/JointPath>

namespace IK{
  class CollisionConstraint : public IKConstraint
  {
  public:
    CollisionConstraint(cnoid::Link* A_link, cnoid::Link* B_link);

    const Eigen::VectorXd& calc_minineq () override;
    const Eigen::VectorXd& calc_maxineq () override;
    const Eigen::SparseMatrix<double,Eigen::RowMajor>& calc_jacobianineq (const std::vector<cnoid::Body*>& bodies) override;
    std::vector<cnoid::SgNodePtr>& getDrawOnObjects() override;

    void setTolerance(double _tolerance) {tolerance=tolerance;}
    void set_maxvel(double _maxvel) {maxvel=_maxvel;}

    //A_v, B_vはworld系
    virtual double detectDistance(cnoid::Vector3& A_v, cnoid::Vector3& B_v)=0;
  protected:
    cnoid::Link* A_link;
    cnoid::Link* B_link;
    double tolerance;
    double maxvel;

    cnoid::SgLineSetPtr lines;

    double current_distance;
    cnoid::Vector3 A_current_localp;
    cnoid::Vector3 B_current_localp;
    cnoid::Vector3 prev_BA;

    cnoid::JointPath path_A;
    cnoid::JointPath path_B;
    cnoid::JointPath path_BA;
  };


}

#endif
