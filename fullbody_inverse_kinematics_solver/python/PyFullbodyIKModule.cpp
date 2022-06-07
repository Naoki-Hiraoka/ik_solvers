/**
   @author YoheiKakiuchi
*/

#include <cnoid/PyUtil>
// #include "../BodyUtil.h"
#include "fullbody_inverse_kinematics_solver/FullbodyInverseKinematicsSolverFast.h"

#include <ik_constraint/AngularMomentumConstraint.h>
#include <ik_constraint/COMConstraint.h>
#include <ik_constraint/COMVelocityConstraint.h>
#include <ik_constraint/ClientCollisionConstraint.h>
//#include <ik_constraint/CollisionConstraint.h>
#include <ik_constraint/JointAngleConstraint.h>
#include <ik_constraint/JointLimitConstraint.h>
#include <ik_constraint/JointVelocityConstraint.h>
#include <ik_constraint/PositionConstraint.h>
#include <ik_constraint/IKConstraint.h>

#include <memory>
#include <pybind11/pybind11.h>

using namespace cnoid;
namespace py = pybind11;

//IKConstraint
//AngularMomentum
//COM
//COMVelocity
//ClientCollision
//JointAngle
//JointLimit
//JointVelocity
//Position

typedef std::shared_ptr<IK::IKConstraint >              IKConstraintPtr;
typedef std::shared_ptr<IK::AngularMomentumConstraint > AngularMomentumConstraintPtr;
typedef std::shared_ptr<IK::COMConstraint >             COMConstraintPtr;
typedef std::shared_ptr<IK::COMVelocityConstraint >     COMVelocityConstraintPtr;
typedef std::shared_ptr<IK::ClientCollisionConstraint > ClientCollisionConstraintPtr;
typedef std::shared_ptr<IK::JointAngleConstraint >      JointAngleConstraintPtr;
typedef std::shared_ptr<IK::JointLimitConstraint >      JointLimitConstraintPtr;
typedef std::shared_ptr<IK::JointVelocityConstraint >   JointVelocityConstraintPtr;
typedef std::shared_ptr<IK::PositionConstraint >        PositionConstraintPtr;

using Matrix4RM = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
using Matrix3RM = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;

class Constraints
{
public:
  std::vector< IKConstraintPtr > ikc_list;

public:
  Constraints() { };

public:
  void push_back(IKConstraintPtr &ptr) {
    ikc_list.push_back(ptr);
  }
  size_t size() { return ikc_list.size(); }
  IKConstraintPtr &at(int i) { return ikc_list.at(i); }
};

int solveFullbodyIKLoopFast (const cnoid::BodyPtr& robot,
                             //const std::vector<std::shared_ptr<IK::IKConstraint> >& ikc_list,
                             Constraints &const_,
                             cnoid::VectorX& jlim_avoid_weight_old,
                             const cnoid::VectorX& dq_weight_all,
                             const size_t max_iteration,
                             double wn,
                             int debugLevel)
                             //const size_t max_iteration = 1,
                             //double wn = 1e-6,
                             //int debugLevel = 0)
{
  return fik::solveFullbodyIKLoopFast(robot,
                                      const_.ikc_list,
                                      jlim_avoid_weight_old,
                                      dq_weight_all,
                                      max_iteration,
                                      wn,
                                      debugLevel);
}

class pyIKConstraint : public IK::IKConstraint
{
public:
  using IK::IKConstraint::IKConstraint;

  // trampoline (one for each virtual function)
  virtual bool checkConvergence () override {
    PYBIND11_OVERLOAD_PURE(
      bool, /* Return type */
      IK::IKConstraint,      /* Parent class */
      checkConvergence        /* Name of function in C++ (must match Python name) */
    );
  }
};


PYBIND11_MODULE(FullbodyIK, m)
{
    m.doc() = "fullbody inverse kinematics module";

    py::module::import("cnoid.Util");

    py::class_< Constraints > (m, "Constraints")
      .def(py::init<>())
      .def("size", &Constraints::size)
      .def("at", &Constraints::at)
      .def("__iter__", [](const Constraints &s) { return py::make_iterator(s.ikc_list.begin(), s.ikc_list.end()); },
           py::keep_alive<0, 1>())
      .def("push_back", (void (Constraints::*)(IKConstraintPtr &)) &Constraints::push_back)
      ;

    py::class_< IK::IKConstraint, IKConstraintPtr, pyIKConstraint > (m, "IKConstraint")
      .def(py::init<>())
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
      .def("checkConvergence", &IK::IKConstraint::checkConvergence)
      ;

    py::class_<IK::AngularMomentumConstraint, AngularMomentumConstraintPtr, IK::IKConstraint > (m, "AngularMomentumConstraint")
      .def(py::init<>());

    py::class_<IK::COMConstraint, COMConstraintPtr, IK::IKConstraint > (m, "COMConstraint")
      .def(py::init<>())
      .def_property("A_robot",
                    (cnoid::BodyPtr & (IK::COMConstraint::*)()) &IK::COMConstraint::A_robot,
                    &IK::COMConstraint::set_A_robot)
      .def_property("B_robot",
                    (cnoid::BodyPtr & (IK::COMConstraint::*)()) &IK::COMConstraint::B_robot,
                    &IK::COMConstraint::set_B_robot)
      .def_property("A_localp",
                    (cnoid::Vector3 & (IK::COMConstraint::*)()) &IK::COMConstraint::A_localp,
                    &IK::COMConstraint::set_A_localp)
      .def_property("B_localp",
                    (cnoid::Vector3 & (IK::COMConstraint::*)()) &IK::COMConstraint::B_localp,
                    &IK::COMConstraint::set_B_localp)
      .def_property("eval_R",
                    (cnoid::Matrix3d & (IK::COMConstraint::*)()) &IK::COMConstraint::eval_R,
                    &IK::COMConstraint::set_eval_R)
      .def_property("maxError",
                    (cnoid::Vector3 & (IK::COMConstraint::*)()) &IK::COMConstraint::maxError,
                    &IK::COMConstraint::set_maxError)
      .def_property("precision",
                    (cnoid::Vector3 & (IK::COMConstraint::*)()) &IK::COMConstraint::precision,
                    &IK::COMConstraint::set_precision)
      .def_property("weight",
                    (cnoid::Vector3 & (IK::COMConstraint::*)()) &IK::COMConstraint::weight,
                    &IK::COMConstraint::set_weight)
      ;

    py::class_<IK::COMVelocityConstraint, COMVelocityConstraintPtr, IK::IKConstraint > (m, "COMVelocityConstraint")
      .def(py::init<>());

    py::class_<IK::ClientCollisionConstraint, ClientCollisionConstraintPtr, IK::IKConstraint > (m, "ClientCollisionConstraint")
      .def(py::init<>());

    //py::class_<IK::CollisionConstraint, CollisionConstraintPtr > (m, "CollisionConstraint")
    //  .def(py::init<>());
    py::class_<IK::JointAngleConstraint, JointAngleConstraintPtr, IK::IKConstraint > (m, "JointAngleConstraint")
      .def(py::init<>())
      .def_property("joint",
                    (cnoid::LinkPtr & (IK::JointAngleConstraint::*)()) &IK::JointAngleConstraint::joint,
                    &IK::JointAngleConstraint::set_joint)
      .def_property("targetq",
                    (double & (IK::JointAngleConstraint::*)()) &IK::JointAngleConstraint::targetq,
                    &IK::JointAngleConstraint::set_targetq)
      .def_property("maxError",
                    (double & (IK::JointAngleConstraint::*)()) &IK::JointAngleConstraint::maxError,
                    &IK::JointAngleConstraint::set_maxError)
      .def_property("precision",
                    (double & (IK::JointAngleConstraint::*)()) &IK::JointAngleConstraint::precision,
                    &IK::JointAngleConstraint::set_precision)
      .def_property("weight",
                    (double & (IK::JointAngleConstraint::*)()) &IK::JointAngleConstraint::weight,
                    &IK::JointAngleConstraint::set_weight)
      ;

    py::class_<IK::JointLimitConstraint, JointLimitConstraintPtr, IK::IKConstraint > (m, "JointLimitConstraint")
      .def(py::init<>());

    py::class_<IK::JointVelocityConstraint, JointVelocityConstraintPtr, IK::IKConstraint > (m, "JointVelocityConstraint")
      .def(py::init<>());

    py::class_<IK::PositionConstraint, PositionConstraintPtr, IK::IKConstraint > (m, "PositionConstraint")
      .def(py::init<>())
      .def_property("A_link",
                    (cnoid::LinkPtr & (IK::PositionConstraint::*)()) &IK::PositionConstraint::A_link,
                    &IK::PositionConstraint::set_A_link)
      .def_property("B_link",
                    (cnoid::LinkPtr & (IK::PositionConstraint::*)()) &IK::PositionConstraint::B_link,
                    &IK::PositionConstraint::set_B_link)
      .def_property("eval_link",
                    (cnoid::LinkPtr & (IK::PositionConstraint::*)()) &IK::PositionConstraint::eval_link,
                    &IK::PositionConstraint::set_eval_link)
      .def_property("A_localpos",
                    //(cnoid::Position & (IK::PositionConstraint::*)()) &IK::PositionConstraint::A_localpos,
                    //&IK::PositionConstraint::set_A_localpos)
                    [](IK::PositionConstraint& self) -> Isometry3::MatrixType& { return self.A_localpos().matrix(); },
                    [](IK::PositionConstraint& self, Eigen::Ref<const Matrix4RM> in_pos) {
                      Position p(in_pos); self.set_A_localpos(p); })
      .def_property("B_localpos",
                    //(cnoid::Position & (IK::PositionConstraint::*)()) &IK::PositionConstraint::B_localpos,
                    //&IK::PositionConstraint::set_B_localpos)
                    [](IK::PositionConstraint& self) -> Isometry3::MatrixType& { return self.B_localpos().matrix(); },
                    [](IK::PositionConstraint& self, Eigen::Ref<const Matrix4RM> in_pos) {
                      Position p(in_pos); self.set_B_localpos(p); })
      .def_property("maxError",
                    (cnoid::Vector6 & (IK::PositionConstraint::*)()) &IK::PositionConstraint::maxError,
                    &IK::PositionConstraint::set_maxError)
      .def_property("precision",
                    (cnoid::Vector6 & (IK::PositionConstraint::*)()) &IK::PositionConstraint::precision,
                    &IK::PositionConstraint::set_precision)
      .def_property("weight",
                    (cnoid::Vector6 & (IK::PositionConstraint::*)()) &IK::PositionConstraint::weight,
                    &IK::PositionConstraint::set_weight)
      .def_property("eval_localR",
                    (cnoid::Matrix3d & (IK::PositionConstraint::*)()) &IK::PositionConstraint::eval_localR,
                    &IK::PositionConstraint::set_eval_localR)
      ;

    m.def("solveFullbodyIKLoopFast", &solveFullbodyIKLoopFast,
          py::arg("robot"),
          py::arg("const"),
          py::arg("jlim_avoid_weight_old"),
          py::arg("dq_weight_all"),
          py::arg("max_iteration") = 1,
          py::arg("double wn") = 1e-6,
          py::arg("int debugLevel") = 0);
}
