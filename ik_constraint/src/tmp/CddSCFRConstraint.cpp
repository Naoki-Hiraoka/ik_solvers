#include "CddSCFRConstraint.h"
#include <cddeigen/cddeigen.h>
#include <qhulleigen/qhulleigen.h>

namespace IK{
  CddSCFRConstraint::CddSCFRConstraint(cnoid::Body* robot, const std::vector<std::shared_ptr<RobotConfig::EndEffector> >& endeffectors):
    SCFRConstraint(robot,endeffectors)
  {
  }

  // Ax = b, Cx >= d
  void CddSCFRConstraint::calcProjection(const Eigen::SparseMatrix<double,Eigen::RowMajor>& A, const Eigen::VectorXd& b, const Eigen::SparseMatrix<double,Eigen::RowMajor>& C, const Eigen::VectorXd& d){
    /*
      INPUT:
        A_eq   x + b_eq    = 0
        A_ineq x + b_ineq >= 0
      OUTPUT:
        x = V y + R_nonneg z + R_free w (sum y = 1, y >= 0, z >= 0)
    */
    Eigen::MatrixXd V, R_nonneg, R_free;
    cddeigen::HtoVgmp(A,-b,C,-d,V,R_nonneg,R_free,this->debuglevel);

    // get p_x, p_y component
    this->V2 = V.topRows<2>();
    this->R_nonneg2 = R_nonneg.topRows<2>();
    this->R_free2 = R_free.topRows<2>();

    // hullをとることで点を減らす
    if (this->V2.cols() > 3) qhulleigen::convexhull(this->V2,this->V2);
    if (this->R_nonneg2.cols() > 3) qhulleigen::convexhull(this->R_nonneg2,this->R_nonneg2);
    if (this->R_free2.cols() > 3) qhulleigen::convexhull(this->R_free2,this->R_free2);

    // 近い点を除去することで点を減らす
    std::vector<Eigen::Vector2d> V2_filterd;
    for(size_t i=0;i<this->V2.cols();i++){
      Eigen::Vector2d v = V2.col(i);
      if(std::find_if(V2_filterd.begin(), V2_filterd.end(), [&](Eigen::Vector2d x) { return (x-v).norm() <= 1e-3; }) == V2_filterd.end()){
        V2_filterd.push_back(v);
      }
    }
    this->V2.resize(2,V2_filterd.size());
    for(size_t i=0;i<V2_filterd.size();i++){
      V2.col(i) = V2_filterd[i];
    }
    std::vector<Eigen::Vector2d> R_nonneg2_filterd;
    for(size_t i=0;i<this->R_nonneg2.cols();i++){
      Eigen::Vector2d v = R_nonneg2.col(i);
      if(std::find_if(R_nonneg2_filterd.begin(), R_nonneg2_filterd.end(), [&](Eigen::Vector2d x) { return (x-v).norm() <= 1e-3; }) == R_nonneg2_filterd.end()){
        R_nonneg2_filterd.push_back(v);
      }
    }
    this->R_nonneg2.resize(2,R_nonneg2_filterd.size());
    for(size_t i=0;i<R_nonneg2_filterd.size();i++){
      R_nonneg2.col(i) = R_nonneg2_filterd[i];
    }
    std::vector<Eigen::Vector2d> R_free2_filterd;
    for(size_t i=0;i<this->R_free2.cols();i++){
      Eigen::Vector2d v = R_free2.col(i);
      if(std::find_if(R_free2_filterd.begin(), R_free2_filterd.end(), [&](Eigen::Vector2d x) { return (x-v).norm() <= 1e-3; }) == R_free2_filterd.end()){
        R_free2_filterd.push_back(v);
      }
    }
    this->R_free2.resize(2,R_free2_filterd.size());
    for(size_t i=0;i<R_free2_filterd.size();i++){
      R_free2.col(i) = R_free2_filterd[i];
    }

    /*
      INPUT:
        x = V y + R_nonneg z + R_free w (sum y = 1, y >= 0, z >= 0)
      OUTPUT:
        A_eq   x + b_eq    = 0
        A_ineq x + b_ineq >= 0
    */
    Eigen::MatrixXd A2, C2;
    Eigen::VectorXd b2, d2;

    cddeigen::VtoHgmp(this->V2,this->R_nonneg2,this->R_free2,A2,b2,C2,d2,this->debuglevel);

    // create SCFR
    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.reserve(10);//適当
    for(size_t i=0;i<A2.rows();i++){
      for(size_t j=0;j<A2.cols();j++){
        if(A2(i,j)!=0) tripletList.push_back(Eigen::Triplet<double>(i,j,A2(i,j)));
      }
    }
    for(size_t i=0;i<C2.rows();i++){
      for(size_t j=0;j<C2.cols();j++){
        if(C2(i,j)!=0) tripletList.push_back(Eigen::Triplet<double>(A2.rows()+i,j,C2(i,j)));
      }
    }
    this->SCFR_M = Eigen::SparseMatrix<double,Eigen::RowMajor>(A2.rows()+C2.rows(),2);
    this->SCFR_M.setFromTriplets(tripletList.begin(), tripletList.end());
    this->SCFR_u = Eigen::VectorXd(b2.rows()+d2.rows());
    this->SCFR_u.head(b2.rows()) = -b2;
    for(size_t i=b2.rows();i<this->SCFR_u.rows();i++) this->SCFR_u[i] = 1e30;
    this->SCFR_l = Eigen::VectorXd(b2.rows()+d2.rows());
    this->SCFR_l.head(b2.rows()) = -b2;
    this->SCFR_l.tail(d2.rows()) = -d2;

    return;
  }

  void CddSCFRConstraint::updateSCFRlines(){
    if(this->initial_p){
      this->calcSCFR();
      this->initial_p = false;
    }

    if(!this->SCFRlines){
      this->SCFRlines = new cnoid::SgLineSet;
      this->SCFRlines->setLineWidth(3.0);
      this->SCFRlines->getOrCreateColors()->resize(1);
      this->SCFRlines->getOrCreateColors()->at(0) = cnoid::Vector3f(0.0,1.0,0.0);

      this->SCFRlines->getOrCreateVertices()->resize(0);
      this->SCFRlines->colorIndices().resize(0);

      // calculate V hull
      std::vector<cnoid::Vector2> V_hull;
      std::vector<std::vector<int> > hull_index;
      V_hull.reserve(10);//適当
      hull_index.reserve(10);//適当
      if(V2.cols() <= 2){
        for(size_t i=0;i<V2.cols();i++){
          V_hull.push_back(cnoid::Vector2(V2(0,i),V2(1,i)));
        }
      }else{
        // convex hull by qhull.
        Eigen::MatrixXd _V_hull;
        qhulleigen::convexhull(this->V2,_V_hull,hull_index);
        for(size_t i=0;i<_V_hull.cols();i++){
          V_hull.push_back(_V_hull.col(i));
        }
      }

      this->SCFRlines->getOrCreateVertices()->resize(V_hull.size());
      for(size_t i=0;i<V_hull.size();i++){
        this->SCFRlines->vertices()->at(i) = cnoid::Vector3f(V_hull[i][0],V_hull[i][1],0);
      }
      for(size_t i=0;i<hull_index.size();i++){
        this->SCFRlines->addLine(hull_index[i][0],hull_index[i][1]);
        this->SCFRlines->colorIndices().push_back(0);
        this->SCFRlines->colorIndices().push_back(0);
      }
    }

    // z座標は毎回更新する
    double z = this->robot->centerOfMass()[2];
    for(size_t i=0;i<this->SCFRlines->vertices()->size();i++){
      this->SCFRlines->vertices()->at(i)[2] = z;
    }
  }

}
