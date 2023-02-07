#pragma once

#include <ros/ros.h>
#include <Eigen/Dense>

class DistMPC
{
public:

  DistMPC(const int& n_dofs, const int& N);
  
  bool c2d (const Eigen::MatrixXd& A,const Eigen::MatrixXd& B, const double& dt,Eigen::MatrixXd& Ad,Eigen::MatrixXd& Bd);
  
  void setSysParams(const Eigen::MatrixXd& A,
                    const Eigen::MatrixXd& B,
                    const Eigen::MatrixXd& C);
  
  bool getSysParams(Eigen::MatrixXd& A,
                    Eigen::MatrixXd& B,
                    Eigen::MatrixXd& C);

  void setCostsParams(const Eigen::MatrixXd& Q1,
                      const Eigen::MatrixXd& R1,
                      const Eigen::MatrixXd& Q2,
                      const Eigen::MatrixXd& R2);
  
  void setC2DSysParams(const Eigen::MatrixXd& A,
                       const Eigen::MatrixXd& B,
                       const Eigen::MatrixXd& C,
                       const double& dt);

  
  Eigen::MatrixXd distCoopMPCGain();
  
  Eigen::MatrixXd distMPCGain(const Eigen::MatrixXd& Q1,
                              const Eigen::MatrixXd& R1,
                              const Eigen::MatrixXd& Q2,
                              const Eigen::MatrixXd& R2,
                              const int& N);
  
  
  bool setAlpha(const double& alpha);
  bool setCurrentState(const Eigen::VectorXd& x);

  bool setReference(const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  
  Eigen::MatrixXd blkdiag(const Eigen::MatrixXd& a, int count);
  
  Eigen::VectorXd controlInputs(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  Eigen::VectorXd step(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  Eigen::VectorXd step(const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r);
  Eigen::VectorXd getCurrentState();
  bool getControlInputs(Eigen::VectorXd& u1, Eigen::VectorXd& u2);

protected:
  
  bool computeControlInputs(Eigen::VectorXd& u1, Eigen::VectorXd& u2);
  Eigen::MatrixXd Ylowtriangular(const Eigen::MatrixXd& A,const Eigen::MatrixXd& B,const Eigen::MatrixXd& C, const int& N);
  Eigen::MatrixXd YColMat(const Eigen::MatrixXd& A,const Eigen::MatrixXd& C, const int& N);
  Eigen::MatrixXd pow(const Eigen::MatrixXd& mat, const int& pow);

  
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;
  Eigen::MatrixXd C_;
  Eigen::VectorXd X_;
  
  Eigen::MatrixXd theta_;
  Eigen::MatrixXd psi_; 
  
  Eigen::MatrixXd Q1_;
  Eigen::MatrixXd R1_;
  Eigen::MatrixXd Q2_;
  Eigen::MatrixXd R2_;
  
  Eigen::MatrixXd K_mpc_;
  Eigen::VectorXd U_;
  
  Eigen::VectorXd reference_;

  bool state_ok_;
  bool sys_params_set_;
  bool cost_params_set_;
  bool gains_set_;
  bool alpha_set_;
  bool control_set_;
  
  int n_dofs_;
  int N_;
  double alpha_;
};

