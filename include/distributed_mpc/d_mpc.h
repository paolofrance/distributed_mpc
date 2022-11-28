#pragma once

#include <ros/ros.h>
#include <eigen3/Eigen/Dense>

class DistMPC
{
public:

  DistMPC(const int& n_dofs, const int& N, const double& dt);
  
  bool c2d (const Eigen::MatrixXd& A,const Eigen::MatrixXd& B, const double& dt,Eigen::MatrixXd& Ad,Eigen::MatrixXd& Bd);
  
  void setSysParams(const Eigen::MatrixXd& A,
                    const Eigen::MatrixXd& B,
                    const Eigen::MatrixXd& C);

  void setCostsParams(const Eigen::MatrixXd& Q1,
                      const Eigen::MatrixXd& R1,
                      const Eigen::MatrixXd& Q2,
                      const Eigen::MatrixXd& R2);
  
  void setC2DSysParams(const Eigen::MatrixXd& A,
                       const Eigen::MatrixXd& B,
                       const Eigen::MatrixXd& C);

  
  Eigen::MatrixXd distCoopMPCGain();
  
  Eigen::MatrixXd distMPCGain(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& B,
                              const Eigen::MatrixXd& C,
                              const Eigen::MatrixXd& Q1,
                              const Eigen::MatrixXd& R1,
                              const Eigen::MatrixXd& Q2,
                              const Eigen::MatrixXd& R2,
                              const int& N);
  
  
  void setAlpha(const double& alpha);
  void setInitialState(const Eigen::VectorXd& x);
  void setCurrentState(const Eigen::VectorXd& x);

  bool setReference(const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  void getControlInputs(Eigen::VectorXd& u1, Eigen::VectorXd& u2);
  
  Eigen::MatrixXd blkdiag(const Eigen::MatrixXd& a, int count);
  
  Eigen::VectorXd controlInputs(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  Eigen::VectorXd step(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  Eigen::VectorXd step(const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  Eigen::VectorXd getCurrentState();

protected:
  
  Eigen::MatrixXd Ylowtriangular(const Eigen::MatrixXd& A,const Eigen::MatrixXd& B,const Eigen::MatrixXd& C, const int& N);
  Eigen::MatrixXd YColMat(const Eigen::MatrixXd& A,const Eigen::MatrixXd& C, const int& N);
  
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;
  Eigen::MatrixXd C_;
  Eigen::MatrixXd Q1_;
  Eigen::MatrixXd R1_;
  Eigen::MatrixXd Q2_;
  Eigen::MatrixXd R2_;
  
  Eigen::MatrixXd K_mpc_;
  
  Eigen::VectorXd reference_;
  Eigen::VectorXd X_;

  bool state_ok_;
  
  int n_dofs_;
  int N_;
  double dt_;
  double alpha_;
};

