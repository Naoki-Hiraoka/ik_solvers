## aim to do the same procedure as SampleSR1.cpp
## exec(open('sample.py').read())

import cnoid.Util
from cnoid_pyutil import *
from cnoid import FullbodyIK
import numpy as np

fname = cnoid.Util.getShareDirectory() + '/model/SR1/SR1.body'
robot = loadRobot(fname)

reset_pose_av = np.array([ 0.0, -0.349066, 0.0, 0.820305, -0.471239, 0.0, ## rleg
                           0.523599, 0.0, 0.0, -1.74533, 0.15708, -0.113446, 0.637045, ## rarm
                           0.0, -0.349066, 0.0, 0.820305, -0.471239, 0.0, ## lleg
                           0.523599, 0.0, 0.0, -1.74533, -0.15708, -0.113446, -0.637045, ## larm
                           0.0, 0.0, 0.0 ]);

for idx in range(len(reset_pose_av)):
    robot.joint(idx).q = reset_pose_av[idx]

robot.calcForwardKinematics()
robot.calcCenterOfMass()

flushRobotView('SR1')

constraints = FullbodyIK.Constraints()

## task: rarm to target
constraint = FullbodyIK.PositionConstraint()
constraint.A_link = robot.link('RARM_WRIST_R');
constraint.A_localpos = cnoidPosition(translation=np.array([0.0, 0.0, -0.02]))
#constraint.B_link() = nullptr;
constraint.B_localpos = cnoidPosition(translation=np.array([0.3, -0.2, 0.8]), rotation=cnoid.Util.angleAxis(-1.5, np.array([0,1,0])))
#cnoid.Util.AngleAxis(-1.5, np.array([0,1,0])).toRotationMatrix()
## constraint.B_localpos().linear() = cnoid::Matrix3(cnoid::AngleAxis(-1.5,cnoid::Vector3(0,1,0)));
constraints.push_back(constraint)

## task: larm to target. rotation-axis nil
constraint = FullbodyIK.PositionConstraint()
constraint.A_link = robot.link('LARM_WRIST_R')
constraint.A_localpos = cnoidPosition(translation=np.array([0.0, 0.0,-0.02]))
#constraint->B_link() = nullptr;
constraint.B_localpos = cnoidPosition(translation=np.array([0.3, 0.2, 0.8]), rotation=cnoid.Util.angleAxis(-1.5, np.array([0,1,0])))
#constraint->B_localpos().linear() = cnoid::Matrix3(cnoid::AngleAxis(-1.5,cnoid::Vector3(0,1,0)));
weight = constraint.weight
for idx in range(3):
    weight[idx+3] = 0.0

constraint.weight = weight ## not require??
constraints.push_back(constraint)

## task: rleg to target
constraint = FullbodyIK.PositionConstraint()
constraint.A_link = robot.link('RLEG_ANKLE_R')
constraint.A_localpos = cnoidPosition(translation=np.array([0.0, 0.0, -0.04]))
#constraint->B_link() = nullptr;
constraint.B_localpos = cnoidPosition(translation=np.array([0.0, -0.2, 0.0]))
constraints.push_back(constraint)

## task: lleg to target
constraint = FullbodyIK.PositionConstraint()
constraint.A_link = robot.link('LLEG_ANKLE_R')
constraint.A_localpos = cnoidPosition(translation=np.array([0.0, 0.0, -0.04]))
#constraint->B_link() = nullptr;
constraint.B_localpos = cnoidPosition(translation=np.array([0.0, 0.2, 0.0]))
constraints.push_back(constraint)

## task: COM to target
constraint = FullbodyIK.COMConstraint()
constraint.A_robot = robot
constraint.B_localp = np.array([0.0, 0.0, 0.7])
constraints.push_back(constraint)


## task: joint angle to target
constraint = FullbodyIK.JointAngleConstraint()
constraint.joint = robot.link('CHEST')
constraint.targetq = 0.1
constraints.push_back(constraint)


for const in constraints:
    const.debuglevel = 1


jlim_avoid_weight_old = np.zeros(6 + robot.getNumJoints())
dq_weight_all = np.ones(6 + robot.getNumJoints())

loop = FullbodyIK.solveFullbodyIKLoopFast(robot,
                                          constraints,
                                          jlim_avoid_weight_old,
                                          dq_weight_all,
                                          20,
                                          1e-6,
                                          1)

flushRobotView('SR1')

cntr = 0
for const in constraints:
    const.debuglevel = 1
    if const.checkConvergence():
        print('constraint %d (%s) : converged'%(cntr, const))
    else:
        print('constraint %d (%s) : NOT converged'%(cntr, const))
    cntr = cntr + 1

