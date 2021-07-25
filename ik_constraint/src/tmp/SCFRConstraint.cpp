#include "SCFRConstraint.h"
#include <cnoid/Jacobian>

namespace IK{
  SCFRConstraint::SCFRConstraint(cnoid::Body* _robot, const std::vector<std::shared_ptr<RobotConfig::EndEffector> >& _endeffectors):
    robot(_robot),
    endeffectors(_endeffectors),
    initial_p(true),
    g(9.80665)
  {
  }

  const Eigen::VectorXd& SCFRConstraint::calc_minineq (){
    Eigen::Vector2d cm(this->robot->centerOfMass()[0],this->robot->centerOfMass()[1]);
    this->minineq = this->SCFR_l - this->SCFR_M * cm;
    return this->minineq;
  }

  const Eigen::VectorXd& SCFRConstraint::calc_maxineq (){
    Eigen::Vector2d cm(this->robot->centerOfMass()[0],this->robot->centerOfMass()[1]);
    this->maxineq = this->SCFR_u - this->SCFR_M * cm;
    return this->maxineq;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& SCFRConstraint::calc_jacobianineq (const std::vector<cnoid::Body*>& bodies) {
    int dim = 0;
    for(size_t i=0;i<bodies.size();i++){
      dim += 6 + bodies[i]->numJoints();
    }

    if(this->initial_p){
      this->calcSCFR();
      this->initial_p = false;
    }

    // calc CM jacobian
    Eigen::MatrixXd _CMJ;
    cnoid::calcCMJacobian(this->robot,nullptr,_CMJ); // [joint root]の順
    Eigen::MatrixXd CMJ(_CMJ.rows(),_CMJ.cols()); // [root joint]の順
    CMJ.leftCols(6) = _CMJ.rightCols(6);
    CMJ.rightCols(CMJ.cols()-6) = _CMJ.leftCols(_CMJ.cols()-6);

    Eigen::SparseMatrix<double,Eigen::RowMajor> CMJ_sparse(2,dim);
    int idx = 0;
    for(size_t b=0;b<bodies.size();b++){
      if(bodies[b] == this->robot){
        for(size_t i=0;i<2;i++){
          for(size_t j=0;j<CMJ.cols();j++){
            CMJ_sparse.insert(i,idx+j) = CMJ(i,j);
          }
        }

        break;
      }
      idx += 6 + bodies[b]->numJoints();
    }

    this->jacobianineq = this->SCFR_M * CMJ_sparse;
    return this->jacobianineq;
  }

  std::vector<cnoid::SgNodePtr>& SCFRConstraint::getDrawOnObjects() {

    if(!this->COMmarker){
      this->COMmarker = new cnoid::CrossMarker(0.2/*size*/,cnoid::Vector3f(0.5,1.0,0.0)/*color*/,10/*width*/);
    }

    // update position
    this->updateSCFRlines();
    this->COMmarker->setTranslation(this->robot->centerOfMass().cast<Eigen::Vector3f::Scalar>());

    this->drawOnObjects = std::vector<cnoid::SgNodePtr>{COMmarker,SCFRlines};

    for(size_t i=0;i<this->endeffectors.size();i++){
      std::vector<cnoid::SgNodePtr> contactlines = this->endeffectors[i]->getcontact()->getDrawOnObjects(this->endeffectors[i]->getlink()->T()*this->endeffectors[i]->getlocalpos());
      std::copy(contactlines.begin(),contactlines.end(),std::back_inserter(this->drawOnObjects));
    }

    return this->drawOnObjects;
  }

  void SCFRConstraint::calcSCFR(){
    Eigen::SparseMatrix<double,Eigen::RowMajor> A;
    Eigen::VectorXd b;
    Eigen::SparseMatrix<double,Eigen::RowMajor> C;
    Eigen::VectorXd d;

    // 制約多面体を計算
    this->calcPolyhedra(A,b,C,d);

    // 制約多面体を射影しSCFRを計算
    this->calcProjection(A,b,C,d);

    // 各要素を正規化
    for(size_t i=0;i<this->SCFR_M.rows();i++){
      double norm = this->SCFR_M.row(i).norm();
      this->SCFR_M.coeffRef(i,0) /= norm;
      this->SCFR_M.coeffRef(i,1) /= norm;
      this->SCFR_l[i] /= norm;
      this->SCFR_u[i] /= norm;
    }

  }

  void SCFRConstraint::calcPolyhedra(Eigen::SparseMatrix<double,Eigen::RowMajor>& A, Eigen::VectorXd& b, Eigen::SparseMatrix<double,Eigen::RowMajor>& C, Eigen::VectorXd& d){
    // x = [px py w1 w2 ...]^T
    // wはエンドエフェクタ座標系．エンドエフェクタまわり

    // Grasp Matrix Gx = h
    Eigen::SparseMatrix<double,Eigen::RowMajor> G(6,2+6*this->endeffectors.size());
    Eigen::VectorXd h=Eigen::VectorXd::Zero(6);
    std::vector<Eigen::Triplet<double> > G_tripletList;
    G_tripletList.reserve(6*6*this->endeffectors.size());
    G_tripletList.push_back(Eigen::Triplet<double>(3,1,-this->robot->mass()*this->g));
    G_tripletList.push_back(Eigen::Triplet<double>(4,0,this->robot->mass()*this->g));
    h[2] = this->robot->mass()*this->g;
    for (size_t i=0;i<this->endeffectors.size();i++){
      const cnoid::Position pos = this->endeffectors[i]->getlink()->T() * this->endeffectors[i]->getlocalpos();
      const cnoid::Matrix3& R = pos.linear();
      const cnoid::Matrix3& p_x_R = cnoid::hat(pos.translation()) * R;

      for(size_t row=0;row<3;row++){
        for(size_t col=0;col<3;col++){
          G_tripletList.push_back(Eigen::Triplet<double>(row,2+6*i+col,R(row,col)));
        }
      }
      for(size_t row=0;row<3;row++){
        for(size_t col=0;col<3;col++){
          G_tripletList.push_back(Eigen::Triplet<double>(3+row,2+6*i+col,p_x_R(row,col)));
        }
      }
      for(size_t row=0;row<3;row++){
        for(size_t col=0;col<3;col++){
          G_tripletList.push_back(Eigen::Triplet<double>(3+row,2+6*i+3+col,R(row,col)));
        }
      }
    }
    G.setFromTriplets(G_tripletList.begin(), G_tripletList.end());

    //接触力制約．Aw = b, Cw >= d
    std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor> > As(this->endeffectors.size());
    std::vector<Eigen::VectorXd> bs(this->endeffectors.size());
    std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor> > Cs(this->endeffectors.size());
    std::vector<Eigen::VectorXd> ds(this->endeffectors.size());
    int num_eqs = 0;
    int num_ineqs = 0;
    for (size_t i=0;i<this->endeffectors.size();i++){
      this->endeffectors[i]->getcontact()->getContactConstraint(As[i],bs[i],Cs[i],ds[i]);
      num_eqs += bs[i].rows();
      num_ineqs += ds[i].rows();
    }
    Eigen::SparseMatrix<double,Eigen::RowMajor> A_appended(num_eqs,2+6*this->endeffectors.size());
    Eigen::VectorXd b_appended(num_eqs);
    Eigen::SparseMatrix<double,Eigen::RowMajor> C_appended(num_ineqs,2+6*this->endeffectors.size());
    Eigen::VectorXd d_appended(num_ineqs);
    std::vector<Eigen::Triplet<double> > A_tripletList;
    std::vector<Eigen::Triplet<double> > C_tripletList;
    A_tripletList.reserve(num_eqs*2);//適当
    C_tripletList.reserve(num_ineqs*2);//適当
    int idx_eq=0;
    int idx_ineq=0;
    for(size_t i=0;i<this->endeffectors.size();i++){
      for (int k=0; k < As[i].outerSize(); ++k){
        for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(As[i],k); it; ++it){
          A_tripletList.push_back(Eigen::Triplet<double>(idx_eq+it.row(),2+i*6+it.col(),it.value()));
        }
      }
      b_appended.segment(idx_eq,bs[i].rows()) = bs[i];
      idx_eq += bs[i].rows();

      for (int k=0; k < Cs[i].outerSize(); ++k){
        for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(Cs[i],k); it; ++it){
          C_tripletList.push_back(Eigen::Triplet<double>(idx_ineq+it.row(),2+i*6+it.col(),it.value()));
        }
      }
      d_appended.segment(idx_ineq,ds[i].rows()) = ds[i];
      idx_ineq += ds[i].rows();
    }
    A_appended.setFromTriplets(A_tripletList.begin(), A_tripletList.end());
    C_appended.setFromTriplets(C_tripletList.begin(), C_tripletList.end());

    // all. Ax = b, Cx >= d
    A = Eigen::SparseMatrix<double,Eigen::RowMajor>(G.rows()+A_appended.rows(),2+6*this->endeffectors.size());
    A.middleRows(0,G.rows()) = G;
    A.middleRows(G.rows(),A_appended.rows()) = A_appended;
    b = Eigen::VectorXd(h.rows()+b_appended.rows());
    b.segment(0,h.rows()) = h;
    b.segment(h.rows(),b_appended.rows()) = b_appended;
    C = C_appended;
    d = d_appended;

  }
}
