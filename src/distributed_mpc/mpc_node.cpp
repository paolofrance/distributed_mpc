#include "ros/ros.h"
#include <distributed_mpc/d_mpc.h>

std::vector<double> range(double min, double max, double dt) {
    std::vector<double> range;
    for(int i=0; i<max/dt; i++) {
        range.push_back(min + i*dt);
    }
    return range;
}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "d_mpc");
  ros::NodeHandle n;
  
  int N=50;
  double dt = 0.008;
  
  std::vector<double> time = range(0.0,10.0,dt);
  
  int n_dofs=2;
  
  
  Eigen::VectorXd ref_h; 
  Eigen::VectorXd ref_r; 
  
  if (n_dofs == 1)
  {
    ref_h.resize(n_dofs*time.size());
    ref_r.resize(n_dofs*time.size());

    for (int i = 0;i<ref_h.size();i++)
    {
      ref_h(i) = std::sin(time[i]);
      ref_r(i) = 0.5*std::sin(time[i]);
    }
  }
  
  else if (n_dofs == 2)
  {
    ref_h.resize(n_dofs*time.size());
    ref_r.resize(n_dofs*time.size());
    
    for (int i = 0;i<ref_h.size()/n_dofs;i++)
    {
      ref_h(i*n_dofs) = std::sin(time[i]);
      ref_h(i*n_dofs+1) = std::cos(time[i]);
      ref_r(i*n_dofs) = 0.5*std::sin(time[i]);
      ref_r(i*n_dofs+1) = 0.5*std::cos(time[i]);
    }
  }
  
    
  double m,c,k;
  
  m=10;
  k=0.1;
  c=50;
  
  Eigen::MatrixXd Ac(2 * n_dofs, 2 * n_dofs);
  Ac << Eigen::MatrixXd::Zero(n_dofs, n_dofs), Eigen::MatrixXd::Identity(n_dofs,n_dofs),
      -k / m * Eigen::MatrixXd::Identity(n_dofs,n_dofs), -c / m * Eigen::MatrixXd::Identity(n_dofs,n_dofs);

  Eigen::MatrixXd Bc(2 * n_dofs, n_dofs);
  Bc << Eigen::MatrixXd::Zero(n_dofs, n_dofs),
      1 / m * Eigen::MatrixXd::Identity(n_dofs,n_dofs);

  Eigen::MatrixXd Cc(n_dofs, 2 * n_dofs);
  Cc << Eigen::MatrixXd::Identity(n_dofs,n_dofs), Eigen::MatrixXd::Zero(n_dofs,n_dofs);

  Eigen::MatrixXd Dc = Eigen::MatrixXd::Zero(n_dofs, n_dofs);
  
  DistMPC dMPC(n_dofs,N);
  dMPC.setC2DSysParams(Ac,Bc,Cc,dt);
  
  Eigen::VectorXd x=Eigen::VectorXd::Zero(2*n_dofs);
  dMPC.setCurrentState(x);
  
  Eigen::MatrixXd Qh = Eigen::MatrixXd::Zero(2*n_dofs, 2*n_dofs); Qh.setZero();
  Qh.topLeftCorner(n_dofs, n_dofs) = Eigen::MatrixXd::Identity(n_dofs, n_dofs);
  
  Eigen::MatrixXd Qr = Eigen::MatrixXd::Zero(2*n_dofs, 2*n_dofs); Qr.setZero();
  Qr.bottomRightCorner(n_dofs, n_dofs) = Eigen::MatrixXd::Identity(n_dofs, n_dofs);

  Eigen::MatrixXd Rh = 0.0005 * Eigen::MatrixXd::Identity(n_dofs, n_dofs);
  Eigen::MatrixXd Rr = 0.0001 * Eigen::MatrixXd::Identity(n_dofs, n_dofs);

  double alpha = 0.8;

  dMPC.setCostsParams(Qh,Rh,Qr,Rr);

  dMPC.setAlpha(alpha);  
  
  Eigen::MatrixXd K_mpc = dMPC.distCoopMPCGain();
  dMPC.setReference(ref_h.segment(0,N*n_dofs),ref_r.segment(0,N*n_dofs));  
  
  Eigen::VectorXd X(2*n_dofs); X.setZero();
  
  Eigen::MatrixXd A,B,C;
  dMPC.getSysParams(A,B,C);
  
  Eigen::VectorXd Xa; Xa.resize(4*n_dofs);
  Xa << X,X;
  
  
  
  for (int i = 0;i<10;i++)  
  {
    auto mid = std::chrono::steady_clock::now();
    
    Eigen::MatrixXd K_mpc = dMPC.distCoopMPCGain();
//     dMPC.step(ref_h.segment(i*n_dofs,N*n_dofs), ref_r.segment(i*n_dofs,N*n_dofs));
    

    Eigen::VectorXd control = dMPC.controlInputs(X,ref_h.segment(i*n_dofs,N*n_dofs), ref_r.segment(i*n_dofs,N*n_dofs));

    Eigen::VectorXd uh,ur;
    dMPC.getControlInputs(uh,ur);

    X = A*X+B*uh+B*ur;
    
    ROS_INFO_STREAM("state: "<<X.transpose());
    ROS_INFO_STREAM("uh: "<<uh.transpose());
    ROS_INFO_STREAM("ur: "<<ur.transpose());
    
    auto end = std::chrono::steady_clock::now();
    ROS_INFO_STREAM("time to compute: "<<std::chrono::duration_cast<std::chrono::microseconds>(end- mid).count());
    
  }
  
  
  
  
  
  
  
  return 0;
}





