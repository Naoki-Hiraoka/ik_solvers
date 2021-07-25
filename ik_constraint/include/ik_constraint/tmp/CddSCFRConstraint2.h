#ifndef CDDSCFRCONSTRAINT2_H
#define CDDSCFRCONSTRAINT2_H

#include "CddSCFRConstraint.h"

namespace IK{
  // CddSCFRConstraintの，Cdd計算をする空間の次元をエンドエフェクタ数に依存しない(8次元)ようにしたもの
  // 制約が増えるためむしろ遅くなる．

  class CddSCFRConstraint2 : public CddSCFRConstraint
  {
  public:
    CddSCFRConstraint2(cnoid::Body* robot, const std::vector<std::shared_ptr<RobotConfig::EndEffector> >& endeffectors);
  protected:
    // calcSCFRの中で呼ばれる．制約多面体を求める
    void calcPolyhedra(Eigen::SparseMatrix<double,Eigen::RowMajor>& A, Eigen::VectorXd& b, Eigen::SparseMatrix<double,Eigen::RowMajor>& C, Eigen::VectorXd& d) override;

  };
}

#endif
