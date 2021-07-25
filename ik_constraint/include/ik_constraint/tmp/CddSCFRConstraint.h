#ifndef CDDSCFRCONSTRAINT_H
#define CDDSCFRCONSTRAINT_H

#include "SCFRConstraint.h"

namespace IK{
  // GMPを使用しないと高頻度で計算が失敗する
  // GMPを使用すると非常に遅い
  class CddSCFRConstraint : public SCFRConstraint
  {
  public:
    CddSCFRConstraint(cnoid::Body* robot, const std::vector<std::shared_ptr<RobotConfig::EndEffector> >& endeffectors);
  protected:
    //SCFR_M, SCFR_u, SCFR_lをセットする
    void calcProjection(const Eigen::SparseMatrix<double,Eigen::RowMajor>& A, const Eigen::VectorXd& b, const Eigen::SparseMatrix<double,Eigen::RowMajor>& C, const Eigen::VectorXd& d) override;
    void updateSCFRlines() override;

    //SCFRのV-表現
    Eigen::MatrixXd V2;
    Eigen::MatrixXd R_nonneg2;
    Eigen::MatrixXd R_free2;
  };
}

#endif
