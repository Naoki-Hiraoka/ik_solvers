#ifndef LPSCFRCONSTRAINT_H
#define LPSCFRCONSTRAINT_H

#include "SCFRConstraint.h"
#include <clpeigen/clpeigen.h>

namespace IK{
  // Testing Static Equilibuium for Legged Robots (2008)の方法

  class LPSCFRConstraint : public SCFRConstraint
  {
  public:
    LPSCFRConstraint(cnoid::Body* robot, const std::vector<std::shared_ptr<RobotConfig::EndEffector> >& endeffectors);

    void seteqs(double _eps){eps = _eps;}
    void setmaxiter(double _maxiter){maxiter = _maxiter;};
  protected:
    //SCFR_M, SCFR_u, SCFR_lをセットする
    void calcProjection(const Eigen::SparseMatrix<double,Eigen::RowMajor>& A, const Eigen::VectorXd& b, const Eigen::SparseMatrix<double,Eigen::RowMajor>& C, const Eigen::VectorXd& d) override;
    void updateSCFRlines() override;

    double calcArea(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c);
    double calcMidArea(const std::tuple<Eigen::Vector2d,Eigen::Vector2d,double>& v1, const std::tuple<Eigen::Vector2d,Eigen::Vector2d,double>& v2);

    clpeigen::solver solver;

    double eps;
    int maxiter;

    std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> > Y; //半時計回り. first: innervertex, second: outervector, third: この点と次の点で構成させるOuterとInnerの間の三角系の面積

  };
}

#endif
