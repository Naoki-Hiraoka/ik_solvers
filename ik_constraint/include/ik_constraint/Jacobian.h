#ifndef IK_CONSTAINT_JACOBIAN_H
#define IK_CONSTAINT_JACOBIAN_H

#include <vector>
#include <unordered_map>

#include <cnoid/Body>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace IK {
  // world座標系で見た、A - B のヤコビアン. linkがnullptrの場合、localposはworld座標を意味する.
  //   jacobianを新たにコンストラクトし、非ゼロ要素に1を入れる.
  void calc6DofJacobianShape(const std::vector<cnoid::LinkPtr>& joints, //input
                             cnoid::LinkPtr& A_link, //input
                             cnoid::LinkPtr& B_link, //input
                             Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian, //output
                             std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap, //output
                             std::vector<cnoid::LinkPtr>& path_A_joints, //output
                             std::vector<cnoid::LinkPtr>& path_B_joints, //output
                             std::vector<cnoid::LinkPtr>& path_BA_joints, //output
                             int& path_BA_joints_numUpwardConnections //output
                             );
  // world座標系で見た、A - B のヤコビアン. linkがnullptrの場合、localposはworld座標を意味する.
  //   jacobianの形状は上の関数で既に整えられている前提.
  void calc6DofJacobianCoef(const std::vector<cnoid::LinkPtr>& joints, //input
                            const cnoid::LinkPtr& A_link, //input
                            const cnoid::Position& A_localpos, //input
                            const cnoid::LinkPtr& B_link, //input
                            const cnoid::Position& B_localpos, //input
                            std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap, //input
                            const std::vector<cnoid::LinkPtr>& path_A_joints, //input
                            const std::vector<cnoid::LinkPtr>& path_B_joints, //input
                            const std::vector<cnoid::LinkPtr>& path_BA_joints, //input
                            const int& path_BA_joints_numUpwardConnections, //input
                            Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian //output
                            );

}

#endif
