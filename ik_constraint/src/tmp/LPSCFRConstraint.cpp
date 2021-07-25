#include "LPSCFRConstraint.h"

namespace IK{
  LPSCFRConstraint::LPSCFRConstraint(cnoid::Body* robot, const std::vector<std::shared_ptr<RobotConfig::EndEffector> >& endeffectors):
    SCFRConstraint(robot,endeffectors),
    eps(0.05),
    maxiter(30)
  {
  }

  // Ax = b, Cx >= d
  void LPSCFRConstraint::calcProjection(const Eigen::SparseMatrix<double,Eigen::RowMajor>& A, const Eigen::VectorXd& b, const Eigen::SparseMatrix<double,Eigen::RowMajor>& C, const Eigen::VectorXd& d){

    // initialize solver
    // lbM <= Mx <= ubM
    // lb <= x <= ub
    Eigen::SparseMatrix<double,Eigen::RowMajor> M(A.rows()+C.rows(),A.cols());
    M.topRows(A.rows()) = A;
    M.bottomRows(C.rows()) = C;
    Eigen::VectorXd lbM(A.rows()+C.rows());
    lbM.head(b.rows()) = b;
    lbM.tail(d.rows()) = d;
    Eigen::VectorXd ubM(A.rows()+C.rows());
    ubM.head(b.rows()) = b;
    for(size_t i=b.rows();i<ubM.rows();i++) ubM[i] = std::numeric_limits<double>::max();
    Eigen::VectorXd lb(A.cols());
    for(size_t i=0;i<lb.rows();i++) lb[i] = -std::numeric_limits<double>::max();
    Eigen::VectorXd ub(A.cols());
    for(size_t i=0;i<ub.rows();i++) ub[i] = std::numeric_limits<double>::max();

    Eigen::VectorXd o = Eigen::VectorXd::Zero(M.cols());
    this->solver.initialize(o,M,lbM,ubM,lb,ub,this->debuglevel);

    Eigen::VectorXd solution;
    this->Y.clear();

    // first 4 solve
    std::vector<Eigen::Vector2d> outer(4);
    outer[0] = Eigen::Vector2d(1,0);
    outer[1] = Eigen::Vector2d(0,1);
    outer[2] = Eigen::Vector2d(-1,0);
    outer[3] = Eigen::Vector2d(0,-1);
    for(size_t i=0;i<outer.size();i++){
      o.head<2>() = outer[i];
      this->solver.updateObjective(o);
      this->solver.solve();
      this->solver.getSolution(solution);
      Y.push_back(std::tuple<Eigen::Vector2d,Eigen::Vector2d,double>(solution.head<2>(),outer[i],0));
    }

    //calc initial area_Y_inner area_Y_mid(=area_Y_outer - area_Y_inner)
    double area_Y_inner = 0;
    double area_Y_mid = 0;
    area_Y_inner += this->calcArea(std::get<0>(*std::next(Y.begin(),0)),
                                   std::get<0>(*std::next(Y.begin(),1)),
                                   std::get<0>(*std::next(Y.begin(),2)));
    area_Y_inner += this->calcArea(std::get<0>(*std::next(Y.begin(),2)),
                                   std::get<0>(*std::next(Y.begin(),3)),
                                   std::get<0>(*std::next(Y.begin(),0)));
    for(std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> >::iterator it=Y.begin();it!=Y.end();it++){
      std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> >::iterator nextit = std::next(it);
      if(nextit==Y.end()) nextit = Y.begin();
      double area = this->calcMidArea(*it,*nextit);
      std::get<2>(*it) = area;
      area_Y_mid += area;

    }

    for(size_t i=0;i<maxiter;i++){
      if(area_Y_inner / (area_Y_inner+area_Y_mid) > 1 - this->eps) break; //近似が概ね収束

      // 最大値を探す
      double maxvalue = -1;
      std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> >::iterator maxit;
      for(std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> >::iterator it=Y.begin();it!=Y.end();it++){
        if(std::get<2>(*it) > maxvalue){
          maxit = it;
          maxvalue = std::get<2>(*it);
        }
      }
      std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> >::iterator nextit = std::next(maxit);
      if (nextit == Y.end()) nextit = Y.begin();

      Eigen::Vector2d p1 = std::get<0>(*maxit);
      Eigen::Vector2d p2 = std::get<0>(*nextit);
      Eigen::Vector2d n = Eigen::Vector2d((p2-p1)[1],-(p2-p1)[0]).normalized();

      // solve LP
      o.head<2>() = n;
      this->solver.updateObjective(o);
      this->solver.solve();
      this->solver.getSolution(solution);
      Eigen::Vector2d pX = solution.head<2>();

      area_Y_inner += this->calcArea(p1,pX,p2);
      area_Y_mid -= std::get<2>(*maxit);
      std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> vX(pX,n,0);
      std::get<2>(*maxit) = this->calcMidArea(*maxit,vX);
      area_Y_mid += std::get<2>(*maxit);
      std::get<2>(vX) = this->calcMidArea(vX,*nextit);
      area_Y_mid += std::get<2>(vX);
      Y.insert(std::next(maxit),vX);
    }

    this->SCFR_M = Eigen::SparseMatrix<double,Eigen::RowMajor>(Y.size(),2);
    this->SCFR_u = Eigen::VectorXd(Y.size());
    this->SCFR_l = Eigen::VectorXd(Y.size());
    size_t i=0;
    for(std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> >::iterator it=Y.begin();it!=Y.end();it++){
      this->SCFR_M.insert(i,0) = std::get<1>(*it)[0];
      this->SCFR_M.insert(i,1) = std::get<1>(*it)[1];
      this->SCFR_u[i] = std::get<1>(*it).dot(std::get<0>(*it));
      this->SCFR_l[i] = -1e30;
      i++;
    }
  }

  void LPSCFRConstraint::updateSCFRlines(){
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

      this->SCFRlines->getOrCreateVertices()->resize(Y.size());
      {
        size_t i=0;
        for(std::list<std::tuple<Eigen::Vector2d,Eigen::Vector2d,double> >::iterator it=Y.begin();it!=Y.end();it++){
          this->SCFRlines->vertices()->at(i) = cnoid::Vector3f(std::get<0>(*it)[0],std::get<0>(*it)[1],0);
          i++;
        }
      }
      for(size_t i=0;i<Y.size();i++){
        int next = (i+1 != Y.size())? i+1 : 0;
        this->SCFRlines->addLine(i,next);
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

  double LPSCFRConstraint::calcArea(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c){
      Eigen::Vector2d v1 = b - a;
      Eigen::Vector2d v2 = c - a;
      return (v1[0]*v2[1]-v1[1]*v2[0])/2;
  }

  double LPSCFRConstraint::calcMidArea(const std::tuple<Eigen::Vector2d,Eigen::Vector2d,double>& v1, const std::tuple<Eigen::Vector2d,Eigen::Vector2d,double>& v2){
    const Eigen::Vector2d p1 = std::get<0>(v1);
    const Eigen::Vector2d p2 = std::get<0>(v2);

    // 2点が一致している
    if((p1 - p2).norm() < 1e-6) return 0;

    // 2点が無限遠にある
    if(p1[0] > std::numeric_limits<double>::max()/2 && p2[0] > std::numeric_limits<double>::max()/2) return 0;
    if(p1[0] < -std::numeric_limits<double>::max()/2 && p2[0] < -std::numeric_limits<double>::max()/2) return 0;
    if(p1[1] > std::numeric_limits<double>::max()/2 && p2[1] > std::numeric_limits<double>::max()/2) return 0;
    if(p1[1] < -std::numeric_limits<double>::max()/2 && p2[1] < -std::numeric_limits<double>::max()/2) return 0;

    const Eigen::Vector2d n1 = std::get<1>(v1);
    const Eigen::Vector2d n2 = std::get<1>(v2);

    // pX: outerの交点を求める
    Eigen::Matrix<double,2,2> tmp;
    tmp.row(0) = n1.transpose(); tmp.row(1) = n2.transpose();
    Eigen::Vector2d pX = tmp.inverse() * Eigen::Vector2d(p1.dot(n1),p2.dot(n2));

    return this->calcArea(p1,pX,p2);
  }

}
