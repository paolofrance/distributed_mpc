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
  
  int N=2;
  int n_dofs=1;
  double dt = 0.01;
  DistMPC dMPC(n_dofs,N,dt);
  
  std::vector<double> time = range(0.0,10.0,dt);
  
  
  Eigen::VectorXd ref_h; ref_h.resize(time.size());
  Eigen::VectorXd ref_r; ref_r.resize(time.size());
  
  for (int i = 0;i<time.size();i++)
  {
    ref_h(i) = std::sin(time[i]);
    ref_r(i) = 0.5*std::sin(time[i]);
  }
  
  
  Eigen::MatrixXd Ac;
  Eigen::MatrixXd Bc;
  Eigen::MatrixXd C;
  Ac.resize(2,2);
  Bc.resize(2,1);
  C.resize(1,2);
  
  double m,c,k;
  
  m=10;
  k=0.1;
  c=50;

  
  Ac << 0, 1,
        -k/m, -c/m;
        
  Bc << 0,
        1/m;
       
  C << 1, 0;
  
  dMPC.setC2DSysParams(Ac,Bc,C);
  
  Eigen::VectorXd x=Eigen::VectorXd::Zero(4);
  dMPC.setInitialState(x);

  Eigen::MatrixXd Qh; Qh.resize(2,2); 
  Qh <<1,0,
       0,0;
  Eigen::MatrixXd Qr; Qr.resize(2,2);
  Qr <<0,0,
       0,1;
  Eigen::MatrixXd Rh; Rh.resize(1,1); Rh<< .0005;
  Eigen::MatrixXd Rr; Rr.resize(1,1); Rr<< .0001;
  double alpha = 0.8;

  dMPC.setCostsParams(Qh,Rh,Qr,Rr);

  dMPC.setAlpha(alpha);  
  Eigen::MatrixXd K_mpc = dMPC.distCoopMPCGain();
  dMPC.setReference(ref_h.segment(0,N),ref_r.segment(0,N));
  Eigen::VectorXd uh,ur;
  dMPC.getControlInputs(uh,ur);
  
  for (int i = 0;i<100;i++)  
  {
    dMPC.step(ref_h.segment(i,i+N-1), ref_r.segment(i,i+N-1));
    
    ROS_INFO_STREAM(dMPC.getCurrentState());
    
  }
  
  
  
  
  
  
  
  return 0;
}





