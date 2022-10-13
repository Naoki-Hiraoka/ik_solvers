#include <prioritized_inverse_kinematics_solver/PrioritizedInverseKinematicsSolver.h>
#include <prioritized_qp_osqp/prioritized_qp_osqp.h>
#include <ik_constraint/PositionConstraint.h>
#include <ik_constraint/COMConstraint.h>
#include <ik_constraint/JointAngleConstraint.h>

#include <cnoid/BodyLoader>
#include <ros/package.h>

int main(void){
  // load robot
  std::string modelfile = ros::package::getPath("choreonoid") + "/share/model/SR1/SR1.body";
  cnoid::BodyLoader bodyLoader;
  cnoid::BodyPtr robot = bodyLoader.load(modelfile);

  // reset manip pose
  std::vector<double> reset_manip_pose{
    0.0, -0.349066, 0.0, 0.820305, -0.471239, 0.0,// rleg
      0.523599, 0.0, 0.0, -1.74533, 0.15708, -0.113446, 0.637045,// rarm
      0.0, -0.349066, 0.0, 0.820305, -0.471239, 0.0,// lleg
      0.523599, 0.0, 0.0, -1.74533, -0.15708, -0.113446, -0.637045,// larm
      0.0, 0.0, 0.0};

  for(int j=0; j < robot->numJoints(); ++j){
    robot->joint(j)->q() = reset_manip_pose[j];
  }
  robot->calcForwardKinematics();
  robot->calcCenterOfMass();

  // setup tasks
  std::vector<std::shared_ptr<IK::IKConstraint> > constraints0;
  std::vector<std::shared_ptr<IK::IKConstraint> > constraints1;
  {
    // task: rleg to target
    std::shared_ptr<IK::PositionConstraint> constraint = std::make_shared<IK::PositionConstraint>();
    constraint->A_link() = robot->link("RLEG_ANKLE_R");
    constraint->A_localpos().translation() = cnoid::Vector3(0.0,0.0,-0.04);
    constraint->B_link() = nullptr;
    constraint->B_localpos().translation() = cnoid::Vector3(0.0,-0.2,-0.0);
    constraints1.push_back(constraint);
  }
  {
    // task: lleg to target
    std::shared_ptr<IK::PositionConstraint> constraint = std::make_shared<IK::PositionConstraint>();
    constraint->A_link() = robot->link("LLEG_ANKLE_R");
    constraint->A_localpos().translation() = cnoid::Vector3(0.0,0.0,-0.04);
    constraint->B_link() = nullptr;
    constraint->B_localpos().translation() = cnoid::Vector3(0.0,0.2,0.0);
    constraints1.push_back(constraint);
  }
  {
    // task: COM to target
    std::shared_ptr<IK::COMConstraint> constraint = std::make_shared<IK::COMConstraint>();
    constraint->A_robot() = robot;
    constraint->B_localp() = cnoid::Vector3(0.0,0.0,0.7);
    constraints1.push_back(constraint);
  }

  std::vector<std::shared_ptr<IK::IKConstraint> > constraints2;
  {
    // task: rarm to target. never reach
    std::shared_ptr<IK::PositionConstraint> constraint = std::make_shared<IK::PositionConstraint>();
    constraint->A_link() = robot->link("RARM_WRIST_R");
    constraint->A_localpos().translation() = cnoid::Vector3(0.0,0.0,-0.02);
    constraint->B_link() = nullptr;
    constraint->B_localpos().translation() = cnoid::Vector3(0.6,-0.2,0.8);
    constraint->B_localpos().linear() = cnoid::Matrix3(cnoid::AngleAxis(-1.5,cnoid::Vector3(0,1,0)));
    constraints2.push_back(constraint);
  }
  {
    // task: larm to target. rotation-axis nil
    std::shared_ptr<IK::PositionConstraint> constraint = std::make_shared<IK::PositionConstraint>();
    constraint->A_link() = robot->link("LARM_WRIST_R");
    constraint->A_localpos().translation() = cnoid::Vector3(0.0,0.0,-0.02);
    constraint->B_link() = nullptr;
    constraint->B_localpos().translation() = cnoid::Vector3(0.3,0.2,0.8);
    constraint->B_localpos().linear() = cnoid::Matrix3(cnoid::AngleAxis(-1.5,cnoid::Vector3(0,1,0)));
    for(size_t i=0;i<3;i++)constraint->weight()[3+i] = 0.0;
    constraints2.push_back(constraint);
  }

  std::vector<std::shared_ptr<IK::IKConstraint> > constraints3;
  {
    // task: joint angle to target
    std::shared_ptr<IK::JointAngleConstraint> constraint = std::make_shared<IK::JointAngleConstraint>();
    constraint->joint() = robot->link("CHEST");
    constraint->targetq() = 0.1;
    constraints3.push_back(constraint);
  }


  std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks;
  std::vector<cnoid::LinkPtr> variables;
  variables.push_back(robot->rootLink());
  for(size_t i=0;i<robot->numJoints();i++){
    variables.push_back(robot->joint(i));
  }
  std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > > constraints{constraints0,constraints1,constraints2,constraints3};
  for(size_t i=0;i<constraints.size();i++){
    for(size_t j=0;j<constraints[i].size();j++){
      constraints[i][j]->debuglevel() = 1;//debug
    }
  }
  prioritized_inverse_kinematics_solver::IKParam param;
  param.debugLevel = 1;
  param.maxIteration = 40;
  int loop = prioritized_inverse_kinematics_solver::solveIKLoop(variables,
                                                                constraints,
                                                                tasks,
                                                                param);

  std::cerr << "loop: " << loop << std::endl;

  for(size_t i=0;i<constraints.size();i++){
    for(size_t j=0;j<constraints[i].size();j++){
      constraints[i][j]->debuglevel() = 0;//not debug
      if(constraints[i][j]->checkConvergence()) std::cerr << "constraint " << i << " " << j << ": converged"<< std::endl;
      else std::cerr << "constraint " << i << ": NOT converged"<< std::endl;
    }
  }
}
