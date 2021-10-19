#include <ik_constraint/Jacobian.h>

#include <cnoid/LinkPath>
#include <cnoid/Jacobian>

namespace IK {
  size_t getJointDOF(const cnoid::LinkPtr& joint) {
    if(joint->isRevoluteJoint() || joint->isPrismaticJoint()) return 1;
    else if(joint->isFreeJoint()) return 6;
    else return 0;
  }

  void pushBackTripletList(std::vector<Eigen::Triplet<double> >& tripletList, const cnoid::LinkPtr& joint, int idx){
    if(joint->isFreeJoint()){
      for(size_t d=0;d<6;d++){
        tripletList.push_back(Eigen::Triplet<double>(d,idx+d,1));
      }
      //  0     p[2] -p[1]
      // -p[2]  0     p[0]
      //  p[1] -p[0]  0
      tripletList.push_back(Eigen::Triplet<double>(0,idx+4,1));
      tripletList.push_back(Eigen::Triplet<double>(0,idx+5,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx+3,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx+5,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx+3,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx+4,1));

    } else if(joint->isRotationalJoint()){
      tripletList.push_back(Eigen::Triplet<double>(0,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(3,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(4,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(5,idx,1));

    } else if(joint->isPrismaticJoint()){
      tripletList.push_back(Eigen::Triplet<double>(0,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(1,idx,1));
      tripletList.push_back(Eigen::Triplet<double>(2,idx,1));

    }
  }

  void fillJacobian(Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian, const cnoid::Vector3& target_p, const cnoid::LinkPtr& joint, int idx, int sign){
    if(joint->isFreeJoint()){
      //root 6dof
      for(size_t j=0;j<6;j++){
        jacobian.coeffRef(j,idx+j) = sign;
      }
      cnoid::Vector3 dp = target_p - joint->p();
      //  0     p[2] -p[1]
      // -p[2]  0     p[0]
      //  p[1] -p[0]  0
      jacobian.coeffRef(0,idx+4)=sign*dp[2];
      jacobian.coeffRef(0,idx+5)=-sign*dp[1];
      jacobian.coeffRef(1,idx+3)=-sign*dp[2];
      jacobian.coeffRef(1,idx+5)=sign*dp[0];
      jacobian.coeffRef(2,idx+3)=sign*dp[1];
      jacobian.coeffRef(2,idx+4)=-sign*dp[0];

    } else if(joint->isRotationalJoint()){
      cnoid::Vector3 omega = joint->R() * joint->a();
      cnoid::Vector3 dp = omega.cross(target_p - joint->p());
      jacobian.coeffRef(0,idx)=sign*dp[0];
      jacobian.coeffRef(1,idx)=sign*dp[1];
      jacobian.coeffRef(2,idx)=sign*dp[2];
      jacobian.coeffRef(3,idx)=sign*omega[0];
      jacobian.coeffRef(4,idx)=sign*omega[1];
      jacobian.coeffRef(5,idx)=sign*omega[2];

    } else if(joint->isPrismaticJoint()){
      cnoid::Vector3 dp = joint->R() * joint->d();
      jacobian.coeffRef(0,idx)=sign*dp[0];
      jacobian.coeffRef(1,idx)=sign*dp[1];
      jacobian.coeffRef(2,idx)=sign*dp[2];
    }
  }

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
                             ){
    jacobianColMap.clear();
    int num_variables = 0;
    for(size_t i=0;i<joints.size();i++){
      jacobianColMap[joints[i]] = num_variables;
      num_variables += getJointDOF(joints[i]);
    }

    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.reserve(100);//適当

    if(!A_link || !B_link || !(A_link->body() == B_link->body())){
      // A, Bが関節を共有しない. 別々に処理すれば良い
      for(size_t i=0;i<2;i++){//0:A_link, 1:B_link
        cnoid::LinkPtr target_link = i ? B_link : A_link;
        std::vector<cnoid::LinkPtr>& path_joints = i ? path_B_joints : path_A_joints;

        if(!target_link) continue;//world固定なので飛ばす

        path_joints.clear();
        cnoid::LinkPath path(target_link);
        for(size_t j=0;j<path.size();j++){
          path_joints.push_back(path[j]);
        }

        for(size_t j=0;j<path_joints.size();j++){
          cnoid::LinkPtr joint = path_joints[j];
          if(jacobianColMap.find(joint)==jacobianColMap.end()) continue;
          int idx = jacobianColMap[joint];
          pushBackTripletList(tripletList,joint,idx);
        }
      }
    } else { //if(!A_link || !B_link || !(A_link->body() == B_link->body()))
      //A,Bが関節を共有する. 一つのpathで考える
      path_BA_joints.clear();
      {
        cnoid::LinkPath path(B_link,A_link);
        size_t j=0;
        for(;!path.isDownward(j);j++) path_BA_joints.push_back(path[j]);
        path_BA_joints_numUpwardConnections = j;
        j++;
        for(;j<path.size();j++) path_BA_joints.push_back(path[j]);
      }

      for(size_t j=0;j<path_BA_joints.size();j++){
        cnoid::LinkPtr joint = path_BA_joints[j];
        if(jacobianColMap.find(joint)==jacobianColMap.end()) continue;
        int idx = jacobianColMap[joint];
        pushBackTripletList(tripletList,joint,idx);
      }
    }

    jacobian = Eigen::SparseMatrix<double,Eigen::RowMajor>(6,num_variables);
    jacobian.setFromTriplets(tripletList.begin(), tripletList.end());

  }

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
                            ) {
    if(!A_link || !B_link || !(A_link->body() == B_link->body())){
      // A, Bが関節を共有しない. 別々に処理すれば良い
      for(size_t i=0;i<2;i++){//0:A_link, 1:B_link
        int sign = i ? -1 : 1;
        cnoid::LinkPtr target_link = i ? B_link : A_link;
        const cnoid::Position& target_localpos = i ? B_localpos : A_localpos;
        const std::vector<cnoid::LinkPtr>& path_joints = i ? path_B_joints : path_A_joints;

        if(!target_link) continue;//world固定なので飛ばす

        const cnoid::Position target_position = target_link->T() * target_localpos;
        const cnoid::Vector3 target_p = target_position.translation();

        for(size_t j=0;j<path_joints.size();j++){
          cnoid::LinkPtr joint = path_joints[j];
          if(jacobianColMap.find(joint)==jacobianColMap.end()) continue;
          int idx = jacobianColMap[joint];
          fillJacobian(jacobian,target_p,joint,idx,sign);
        }
      }
    } else { //if(!A_link || !B_link || !(A_link->body() == B_link->body()))
      //A,Bが関節を共有する. 一つのpathで考える
      const cnoid::Vector3& target_p = A_link->T() * A_localpos.translation();
      for(size_t j=0;j<path_BA_joints.size();j++){
          cnoid::LinkPtr joint = path_BA_joints[j];
          if(jacobianColMap.find(joint)==jacobianColMap.end()) continue;
          int idx = jacobianColMap[joint];
          if(j<path_BA_joints_numUpwardConnections){
            fillJacobian(jacobian,target_p,joint,idx,-1);
          } else {
            fillJacobian(jacobian,target_p,joint,idx,1);
          }
      }
    }
  }

  void calcCMJacobianShape(const std::vector<cnoid::LinkPtr>& joints,
                           const cnoid::BodyPtr& A_robot,
                           const cnoid::BodyPtr& B_robot,
                           Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian,
                           std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap){
    jacobianColMap.clear();
    int cols = 0;
    for(size_t i=0;i<joints.size();i++){
      jacobianColMap[joints[i]] = cols;
      cols += getJointDOF(joints[i]);
    }
    jacobian = Eigen::SparseMatrix<double,Eigen::RowMajor>(3,cols);

    std::vector<Eigen::Triplet<double> > tripletList;
    tripletList.reserve(100);//適当

    if(A_robot != B_robot){
      for(int i=0;i<2;i++){
        cnoid::BodyPtr robot = (i==0) ? A_robot : B_robot;
        if(!robot) continue;
        if(jacobianColMap.find(robot->rootLink()) != jacobianColMap.end()){
          int idx = jacobianColMap[robot->rootLink()];
          for(size_t row=0;row<3;row++){
            for(size_t j=0;j<getJointDOF(robot->rootLink());j++){
              tripletList.push_back(Eigen::Triplet<double>(row,idx+j,1));
            }
          }
        }
        for(size_t j=0;j<robot->numJoints();j++){
          if(jacobianColMap.find(robot->joint(j)) != jacobianColMap.end()){
            int idx = jacobianColMap[robot->joint(j)];
            for(size_t row=0;row<3;row++){
              for(size_t k=0;k<getJointDOF(robot->joint(j));k++){
                tripletList.push_back(Eigen::Triplet<double>(row,idx+k,1));
              }
            }
          }
        }
      }
    }
    jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  void calcCMJacobianCoef(const std::vector<cnoid::LinkPtr>& joints,//input
                          const cnoid::BodyPtr& A_robot,//input
                          const cnoid::BodyPtr& B_robot,//input
                          std::unordered_map<cnoid::LinkPtr,int>& jacobianColMap, //input
                          Eigen::SparseMatrix<double,Eigen::RowMajor>& jacobian//output
                          ) {
    if(A_robot != B_robot){
      for(size_t i=0;i<2; i++){
        cnoid::BodyPtr robot = (i==0) ? A_robot : B_robot;
        if(!robot) continue;
        int sign = (i==0) ? 1 : -1;

        Eigen::MatrixXd CMJ;
        cnoid::calcCMJacobian(robot,nullptr,CMJ); // [joint root]の順

        if(jacobianColMap.find(robot->rootLink()) != jacobianColMap.end()){
          int col_idx = jacobianColMap[robot->rootLink()];
          for(size_t j=0;j<3;j++){
            for(size_t k=0;k<getJointDOF(robot->rootLink());k++){
              jacobian.coeffRef(j,col_idx+k) = sign * CMJ(j,robot->numJoints()+k);
            }
          }
        }
        for(size_t j=0;j<robot->numJoints();j++){
          if(jacobianColMap.find(robot->joint(j)) != jacobianColMap.end()){
            int col_idx = jacobianColMap[robot->joint(j)];
            for(size_t k=0;k<3;k++){
              for(size_t d=0;d<getJointDOF(robot->joint(j));d++){
                jacobian.coeffRef(k,col_idx+d) = sign * CMJ(k,j+d);
              }
            }
          }
        }
      }
    }

  }
}




