#include "CollisionConstraint.h"
#include <cnoid/EigenUtil>
#include <cnoid/TimeMeasure>

namespace IK{

  CollisionConstraint::CollisionConstraint(cnoid::Link* _A_link, cnoid::Link* _B_link):
    A_link(_A_link),
    B_link(_B_link),
    tolerance(0.01),
    maxvel(0.1),
    A_current_localp(cnoid::Vector3::Zero()),
    B_current_localp(cnoid::Vector3::Zero()),
    current_distance(1.0),
    prev_BA(cnoid::Vector3::UnitX())
  {
  }

  //先にcalc_jacobianineqが呼ばれている前提
  const Eigen::VectorXd& CollisionConstraint::calc_minineq () {
    if(this->minineq.rows()!=1) this->minineq = Eigen::VectorXd(1);
    this->minineq[0] = std::min(tolerance - this->current_distance, maxvel);
    return this->minineq;
  }

  //先にcalc_jacobianineqが呼ばれている前提
  const Eigen::VectorXd& CollisionConstraint::calc_maxineq () {
    if(this->maxineq.rows()!=1){
      this->maxineq = Eigen::VectorXd(1);
      this->maxineq[0] = 1e30;
    }
    return this->maxineq;
  }

  const Eigen::SparseMatrix<double,Eigen::RowMajor>& CollisionConstraint::calc_jacobianineq (const std::vector<cnoid::Body*>& bodies) {
    if(!this->is_bodies_same(bodies,this->jacobianineq_bodies) || this->jacobianineq.rows()!=1){
      this->jacobianineq_bodies = bodies;

      cnoid::TimeMeasure timer;
      if(this->debuglevel>0){
        timer.begin();
      }

      std::vector<Eigen::Triplet<double> > tripletList;
      tripletList.reserve(20);//適当

      if(! (A_link->body() == B_link->body())){
        for(size_t i=0;i<2;i++){//0:A_link, 1:B_link
          cnoid::Link* target_link = i ? B_link : A_link;
          cnoid::JointPath& path = i? this->path_B : this->path_A;

          int idx = 0;
          for(size_t b=0;b<bodies.size();b++){
            if(bodies[b] == target_link->body()){
              //root 6dof
              for(size_t j=0;j<3;j++){
                tripletList.push_back(Eigen::Triplet<double>(0,idx+j,1));
              }
              tripletList.push_back(Eigen::Triplet<double>(0,idx+3,1));
              tripletList.push_back(Eigen::Triplet<double>(0,idx+4,1));
              tripletList.push_back(Eigen::Triplet<double>(0,idx+5,1));

              //joints
              path.setPath(target_link);
              for(size_t j=0;j<path.numJoints();j++){
                int col = idx+6+path.joint(j)->jointId();
                tripletList.push_back(Eigen::Triplet<double>(0,col,1));
              }
              break;
            }
            idx += 6 + bodies[b]->numJoints();
          }
        }
      }else{// if(! (A_link->body() == B_link->body()))
        int idx = 0;
        for(size_t b=0;b<bodies.size();b++){
          if(bodies[b] == A_link->body()){
            //joints
            this->path_BA.setPath(B_link,A_link);

            for(size_t j=0;j<this->path_BA.numJoints();j++){
              int col = idx+6+this->path_BA.joint(j)->jointId();
              tripletList.push_back(Eigen::Triplet<double>(0,col,1));
            }
            break;
          }
          idx += 6 + bodies[b]->numJoints();
        }
      }

      int dim = 0;
      for(size_t i=0; i < bodies.size(); i++){
        dim += 6 + bodies[i]->numJoints();
      }

      this->jacobianineq = Eigen::SparseMatrix<double,Eigen::RowMajor>(1,dim);
      this->jacobianineq.setFromTriplets(tripletList.begin(), tripletList.end());

      if(this->debuglevel>0){
        double time = timer.measure();
        std::cerr << " CollisionConstraint initialize jacobianineq time: " << time
                  << std::endl;
      }
    }

    cnoid::TimeMeasure timer;
    if(this->debuglevel>0){
      timer.begin();
    }

    cnoid::Vector3 A_v, B_v;
    this->current_distance = this->detectDistance(A_v,B_v);
    this->A_current_localp = this->A_link->T().inverse() * A_v;
    this->B_current_localp = this->B_link->T().inverse() * B_v;

    //jacobian A-B
    cnoid::Vector3 BA;
    if(std::abs(this->current_distance)>0.001){
      BA = (A_v - B_v).normalized();
      prev_BA = BA;
    }else{
    //干渉が発生していると，正しく離れる方向を指し示さないことが多い
      BA = prev_BA.normalized();
      //(A_link->wc() - B_link->wc()).normalized()
    }

    if(! (A_link->body() == B_link->body())){
      for(size_t i=0;i<2;i++){//0:A_link, 1:B_link
        int sign = i ? -1 : 1;
        cnoid::Link* target_link = i ? B_link : A_link;
        cnoid::Vector3& target_p = i ? B_v : A_v;
        cnoid::JointPath& path = i? this->path_B : this->path_A;

        int idx = 0;
        for(size_t b=0;b<bodies.size();b++){
          if(bodies[b] == target_link->body()){
            //root 6dof
            for(size_t j=0;j<3;j++){
              this->jacobianineq.coeffRef(0,idx+j) = sign*BA[j];
            }
            cnoid::Vector3 dp = target_p - bodies[b]->rootLink()->p();
            Eigen::Matrix<double, 1, 3> BA_minusdphat = BA.transpose() * - cnoid::hat(dp);
            this->jacobianineq.coeffRef(0,idx+3)=sign*BA_minusdphat[0];
            this->jacobianineq.coeffRef(0,idx+4)=sign*BA_minusdphat[1];
            this->jacobianineq.coeffRef(0,idx+5)=sign*BA_minusdphat[2];

            //joints
            for(size_t j=0;j<path.numJoints();j++){
              int col = idx+6+path.joint(j)->jointId();
              cnoid::Vector3 omega = path.joint(j)->R() * path.joint(j)->a();
              if(!path.isJointDownward(j)) omega = -omega;
              cnoid::Vector3 dp = omega.cross(target_p - path.joint(j)->p());
              this->jacobianineq.coeffRef(0,col)=sign*BA.dot(dp);
            }
            break;
          }
          idx += 6 + bodies[b]->numJoints();
        }
      }
    }else{// if(! (A_link->body() == B_link->body()))
      int idx = 0;
      for(size_t b=0;b<bodies.size();b++){
        if(bodies[b] == A_link->body()){
          //joints
          this->path_BA.setPath(B_link,A_link);

          const cnoid::Vector3& target_p = A_v;

          for(size_t j=0;j<this->path_BA.numJoints();j++){
            int col = idx+6+this->path_BA.joint(j)->jointId();
            cnoid::Vector3 omega = this->path_BA.joint(j)->R() * this->path_BA.joint(j)->a();
            if(!this->path_BA.isJointDownward(j)) omega = -omega;
            cnoid::Vector3 dp = omega.cross(target_p - this->path_BA.joint(j)->p());
            this->jacobianineq.coeffRef(0,col)=BA.dot(dp);
          }
          break;
        }
        idx += 6 + bodies[b]->numJoints();
      }
    }

    if(this->debuglevel>0){
      double time = timer.measure();
      std::cerr << " CollisionConstraint calc jacobianineq time: " << time
                << std::endl;
    }

    return this->jacobianineq;
  }

  std::vector<cnoid::SgNodePtr>& CollisionConstraint::getDrawOnObjects(){
    if(!this->lines){
      this->lines = new cnoid::SgLineSet;
      this->lines->setLineWidth(1.0);
      this->lines->getOrCreateColors()->resize(3);
      this->lines->getOrCreateColors()->at(0) = cnoid::Vector3f(0.3,0.0,0.0);
      this->lines->getOrCreateColors()->at(1) = cnoid::Vector3f(0.6,0.0,0.0);
      this->lines->getOrCreateColors()->at(2) = cnoid::Vector3f(1.0,0.0,0.0);
      // A, B
      this->lines->getOrCreateVertices()->resize(2);
      this->lines->colorIndices().resize(0);
      this->lines->addLine(0,1); this->lines->colorIndices().push_back(0); this->lines->colorIndices().push_back(0);

      this->drawOnObjects = std::vector<cnoid::SgNodePtr>{this->lines};
    }

    cnoid::Vector3 A_v = this->A_link->T() * this->A_current_localp;
    cnoid::Vector3 B_v = this->B_link->T() * this->B_current_localp;
    double d = this->current_distance;

    this->lines->vertices()->at(0) = A_v.cast<cnoid::Vector3f::Scalar>();
    this->lines->vertices()->at(1) = B_v.cast<cnoid::Vector3f::Scalar>();
    if (d < tolerance) {
      this->lines->setLineWidth(3.0);
      this->lines->colorIndices().at(0) = 2;
      this->lines->colorIndices().at(1) = 2;
    } else if (d < tolerance * 2) {
      this->lines->setLineWidth(2.0);
      this->lines->colorIndices().at(0) = 1;
      this->lines->colorIndices().at(1) = 1;
    } else {
      this->lines->setLineWidth(1.0);
      this->lines->colorIndices().at(0) = 0;
      this->lines->colorIndices().at(1) = 0;
    }

    return this->drawOnObjects;
  }

}
