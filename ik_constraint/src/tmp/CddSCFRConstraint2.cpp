#include "CddSCFRConstraint2.h"
#include <cddeigen/cddeigen.h>
#include <qhulleigen/qhulleigen.h>

namespace IK{
  CddSCFRConstraint2::CddSCFRConstraint2(cnoid::Body* robot, const std::vector<std::shared_ptr<RobotConfig::EndEffector> >& endeffectors):
    CddSCFRConstraint(robot,endeffectors)
  {
  }

  void CddSCFRConstraint2::calcPolyhedra(Eigen::SparseMatrix<double,Eigen::RowMajor>& A, Eigen::VectorXd& b, Eigen::SparseMatrix<double,Eigen::RowMajor>& C, Eigen::VectorXd& d) {
    // x = [px py w]^T
    // wはワールド座標系．原点まわり
    this->debuglevel=1;
    // Grasp Matrix Gx = h
    Eigen::SparseMatrix<double,Eigen::RowMajor> G(6,2+6);
    Eigen::VectorXd h=Eigen::VectorXd::Zero(6);
    std::vector<Eigen::Triplet<double> > G_tripletList;
    G_tripletList.reserve(8);
    G_tripletList.push_back(Eigen::Triplet<double>(3,1,-this->robot->mass()*this->g));
    G_tripletList.push_back(Eigen::Triplet<double>(4,0,this->robot->mass()*this->g));
    h[2] = this->robot->mass()*this->g;
    for (size_t i=0;i<6;i++){
      G_tripletList.push_back(Eigen::Triplet<double>(i,2+i,1));
    }
    G.setFromTriplets(G_tripletList.begin(), G_tripletList.end());

    //接触力制約．Aw = b, Cw >= d
    std::vector<Eigen::MatrixXd> Vs(this->endeffectors.size());
    std::vector<Eigen::MatrixXd> R_nonnegs(this->endeffectors.size());
    std::vector<Eigen::MatrixXd> R_frees(this->endeffectors.size());
    int num_V_col = 0;
    int num_R_nonneg_col = 0;
    int num_R_free_col = 0;
    for (size_t i=0;i<this->endeffectors.size();i++){
      Eigen::SparseMatrix<double,Eigen::RowMajor> Ai;
      Eigen::VectorXd bi;
      Eigen::SparseMatrix<double,Eigen::RowMajor> Ci;
      Eigen::VectorXd di;
      // エンドエフェクタ系，エンドエフェクタまわりの制約を取得
      this->endeffectors[i]->getcontact()->getContactConstraint(Ai,bi,Ci,di);

      // SPAN表現へ
      // Ax - b =0, Cx - d >= 0
      cddeigen::HtoVgmp(Ai,-bi,Ci,-di,Vs[i],R_nonnegs[i],R_frees[i],this->debuglevel);

      // ワールド系，原点まわりへ
      Eigen::Matrix<double,6,6> T=Eigen::Matrix<double,6,6>::Zero();
      const cnoid::Position pos = this->endeffectors[i]->getlink()->T() * this->endeffectors[i]->getlocalpos();
      const cnoid::Matrix3& R = pos.linear();
      const cnoid::Matrix3& p_x_R = cnoid::hat(pos.translation()) * R;
      T.topLeftCorner<3,3>() = R;
      T.bottomLeftCorner<3,3>() = p_x_R;
      T.bottomRightCorner<3,3>() = R;
      std::cerr <<"Vs[i]"<<std::endl;
      std::cerr <<Vs[i]<<std::endl;
      Vs[i] = (T * Vs[i]).eval();
      R_nonnegs[i] = (T * R_nonnegs[i]).eval();
      R_frees[i] = (T * R_frees[i]).eval();

      // 列数を足す
      num_V_col += Vs[i].cols();
      num_R_nonneg_col += R_nonnegs[i].cols();
      num_R_free_col += R_frees[i].cols();
    }

    // 全エンドエフェクタを合わせる
    Eigen::MatrixXd V_appended(6,num_V_col);
    Eigen::MatrixXd R_nonneg_appended(6,num_R_nonneg_col);
    Eigen::MatrixXd R_free_appended(6,num_R_free_col);
    int idx_V=0;
    int idx_R_nonneg=0;
    int idx_R_free=0;
    for(size_t i=0;i<this->endeffectors.size();i++){
      std::cerr <<Vs[i]<<std::endl;
      V_appended.block(0,idx_V,6,Vs[i].cols()) = Vs[i];
      idx_V += Vs[i].cols();
      R_nonneg_appended.block(0,idx_R_nonneg,6,R_nonnegs[i].cols()) = R_nonnegs[i];
      idx_R_nonneg += R_nonnegs[i].cols();
      R_free_appended.block(0,idx_R_free,6,R_frees[i].cols()) = R_frees[i];
      idx_R_free += R_frees[i].cols();
    }
    std::cerr << V_appended<<std::endl;
    // FACE表現へ Aw + b =0, Cw + d >= 0
    Eigen::MatrixXd A_w;
    Eigen::VectorXd b_w;
    Eigen::MatrixXd C_w;
    Eigen::VectorXd d_w;
    cddeigen::VtoHgmp(V_appended,R_nonneg_appended,R_free_appended,A_w,b_w,C_w,d_w,this->debuglevel);

    // all. Ax = b, Cx >= d
    A = Eigen::SparseMatrix<double,Eigen::RowMajor>(G.rows()+A_w.rows(),2+6);
    A.middleRows(0,G.rows()) = G;
    for(size_t i=0;i<A_w.rows();i++){
      for(size_t j=0;j<A_w.cols();j++){
        if(A_w(i,j)!=0) A.insert(G.rows()+i,2+j) = A_w(i,j);
      }
    }
    b = Eigen::VectorXd(G.rows()+A_w.rows());
    b.head(h.rows()) = h;
    b.tail(b_w.rows()) = -b_w;
    C = Eigen::SparseMatrix<double,Eigen::RowMajor>(C_w.rows(),2+C_w.cols());
    for(size_t i=0;i<C_w.rows();i++){
      for(size_t j=0;j<C_w.cols();j++){
        if(C_w(i,j)!=0) C.insert(i,2+j) = C_w(i,j);
      }
    }
    d = -d_w;
  }
}
