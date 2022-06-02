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
//COM
//COMVelocity
//ClientCollision
//JointAngle
//JointLimit
//JointVelocity
//Position

class Constraints
{
public:
  std::vector< std::shared_ptr <IK::IKConstraint > > ikc_list;

public:
  Constraints() { };

public:
  void push_back(std::shared_ptr <IK::IKConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  void push_back(std::shared_ptr <IK::COMConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  void push_back(std::shared_ptr <IK::COMVelocityConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  void push_back(std::shared_ptr <IK::ClientCollisionConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  void push_back(std::shared_ptr <IK::JointAngleConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  void push_back(std::shared_ptr <IK::JointLimitConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  void push_back(std::shared_ptr <IK::JointVelocityConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  void push_back(std::shared_ptr <IK::PositionConstraint > &ptr) {
    ikc_list.push_back(ptr);
  }
  size_t size() { return ikc_list.size(); }
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

#if 0
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
#endif

PYBIND11_MODULE(FullbodyIK, m)
{
    m.doc() = "fullbody inverse kinematics module";

    py::module::import("cnoid.Util");

#if 0
    py::class_< IK::IKConstraint, pyIKConstraint > (m, "IKConstraint")
      .def(py::init<>())
      .def("checkConvergence", &IK::IKConstraint::checkConvergence);
#endif

    py::class_< Constraints > (m, "Constraints")
      .def(py::init<>())
      .def("size", &Constraints::size)
      .def("__iter__", [](const Constraints &s) { return py::make_iterator(s.ikc_list.begin(), s.ikc_list.end()); },
           py::keep_alive<0, 1>())
      .def("push_back", (void (Constraints::*)(std::shared_ptr< IK::COMConstraint> &)) &Constraints::push_back)
      .def("push_back", (void (Constraints::*)(std::shared_ptr< IK::COMVelocityConstraint> &)) &Constraints::push_back)
      .def("push_back", (void (Constraints::*)(std::shared_ptr< IK::ClientCollisionConstraint> &)) &Constraints::push_back)
      .def("push_back", (void (Constraints::*)(std::shared_ptr< IK::JointAngleConstraint> &)) &Constraints::push_back)
      .def("push_back", (void (Constraints::*)(std::shared_ptr< IK::JointLimitConstraint> &)) &Constraints::push_back)
      .def("push_back", (void (Constraints::*)(std::shared_ptr< IK::JointVelocityConstraint> &)) &Constraints::push_back)
      .def("push_back", (void (Constraints::*)(std::shared_ptr< IK::PositionConstraint> &)) &Constraints::push_back)
      ;

    py::class_<IK::AngularMomentumConstraint, std::shared_ptr< IK::AngularMomentumConstraint> > (m, "AngularMomentumConstraint")
      .def(py::init<>());

    py::class_<IK::COMConstraint, std::shared_ptr< IK::COMConstraint> > (m, "COMConstraint")
      .def(py::init<>())
      .def("checkConvergence", &IK::COMConstraint::checkConvergence)
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
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
      ;

    py::class_<IK::COMVelocityConstraint, std::shared_ptr< IK::COMVelocityConstraint> > (m, "COMVelocityConstraint")
      .def(py::init<>())
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
      .def("checkConvergence", &IK::COMVelocityConstraint::checkConvergence)
      ;

    py::class_<IK::ClientCollisionConstraint, std::shared_ptr< IK::ClientCollisionConstraint> > (m, "ClientCollisionConstraint")
      .def(py::init<>())
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
      .def("checkConvergence", &IK::ClientCollisionConstraint::checkConvergence)
      ;
    //py::class_<IK::CollisionConstraint, std::shared_ptr< IK::CollisionConstraint> > (m, "CollisionConstraint")
    //  .def(py::init<>());
    py::class_<IK::JointAngleConstraint, std::shared_ptr< IK::JointAngleConstraint> > (m, "JointAngleConstraint")
      .def(py::init<>())
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
      .def("checkConvergence", &IK::JointAngleConstraint::checkConvergence)
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

    py::class_<IK::JointLimitConstraint, std::shared_ptr< IK::JointLimitConstraint> > (m, "JointLimitConstraint")
      .def(py::init<>())
      .def("checkConvergence", &IK::JointLimitConstraint::checkConvergence)
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
      ;
    py::class_<IK::JointVelocityConstraint, std::shared_ptr< IK::JointVelocityConstraint> > (m, "JointVelocityConstraint")
      .def("checkConvergence", &IK::JointVelocityConstraint::checkConvergence)
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
      .def(py::init<>())
      ;

    py::class_<IK::PositionConstraint, std::shared_ptr< IK::PositionConstraint> > (m, "PositionConstraint")
      .def(py::init<>())
      .def("checkConvergence", &IK::PositionConstraint::checkConvergence)
      .def_property("debuglevel", (int & (IK::IKConstraint::*)())&IK::IKConstraint::debuglevel, &IK::IKConstraint::set_debuglevel)
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
                    (cnoid::Position & (IK::PositionConstraint::*)()) &IK::PositionConstraint::A_localpos,
                    &IK::PositionConstraint::set_A_localpos)
      .def_property("B_localpos",
                    (cnoid::Position & (IK::PositionConstraint::*)()) &IK::PositionConstraint::B_localpos,
                    &IK::PositionConstraint::set_B_localpos)
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
