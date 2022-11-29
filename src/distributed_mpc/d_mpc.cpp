#include <distributed_mpc/d_mpc.h>


DistMPC::DistMPC(const int& n_dofs, const int& N, const double& dt): n_dofs_(n_dofs), N_(N),dt_(dt)
{
  A_.resize(2*n_dofs_,2*n_dofs_);
  B_.resize(2*n_dofs_,1);
  C_.resize(n_dofs_,2*n_dofs_);
  X_.resize(2*n_dofs_);
  
  reference_.resize(2*(2*X_.size() + 2*N ),1);
  
  reference_.setZero();
  state_ok_       = false;
  sys_params_set_ = false;
  cost_params_set_= false;
  gains_set_      = false;
}

bool DistMPC::c2d (const Eigen::MatrixXd& A,const Eigen::MatrixXd& B, const double& dt, Eigen::MatrixXd& Ad,Eigen::MatrixXd& Bd)
{
  Eigen::MatrixXd I  = Eigen::MatrixXd::Identity(A.rows(), A.rows());  
  Ad = (Eigen::MatrixXd::Identity(A.rows(), A.rows()) + 0.5*A*dt) * (Eigen::MatrixXd::Identity(A.rows(), A.rows()) - 0.5*A*dt).inverse();
  Bd = A.inverse() * ( Ad-Eigen::MatrixXd::Identity(A.rows(), A.rows() ) )*B;
  return true;
}

void DistMPC::setSysParams( const Eigen::MatrixXd& A,
                            const Eigen::MatrixXd& B,
                            const Eigen::MatrixXd& C)
{
  A_ = A;
  B_ = B;
  C_ = C;
  
  Eigen::MatrixXd Aa = blkdiag(A,2);
  Eigen::MatrixXd Ba;Ba.resize(2*B_.rows(),B_.cols()); 
  Ba << B,
        B;
  Eigen::MatrixXd Ca = blkdiag(C,2);  
  
  theta_ = Ylowtriangular(Aa,Ba,Ca,N_);
  psi_ = YColMat(Aa,Ca,N_); 
  
  sys_params_set_ = true;
  
}

bool DistMPC::getSysParams(Eigen::MatrixXd& A,
                           Eigen::MatrixXd& B,
                           Eigen::MatrixXd& C)
{
  if(!sys_params_set_)
  {
    ROS_ERROR("system params not yet set. return");
    return false;
  }
  
  A = A_;
  B = B_;
  C = C_;
  
  return true;
}

void DistMPC::setC2DSysParams(const Eigen::MatrixXd& A,
                              const Eigen::MatrixXd& B,
                              const Eigen::MatrixXd& C)
{
  Eigen::MatrixXd Ad, Bd;
  c2d (A,B,dt_,Ad,Bd);
  
  setSysParams(Ad,Bd,C);
}


void DistMPC::setCostsParams( const Eigen::MatrixXd& Q1,
                              const Eigen::MatrixXd& R1,
                              const Eigen::MatrixXd& Q2,
                              const Eigen::MatrixXd& R2)
{
  Q1_ = Q1;
  R1_ = R1;
  Q2_ = Q2;
  R2_ = R2;
  cost_params_set_ = true;
}

Eigen::MatrixXd DistMPC::distCoopMPCGain()
{
  if (!sys_params_set_)
  {
    ROS_ERROR_STREAM("System parameters not set. use setC2DSysParams or setSysParams to set the parameters before !");
  }
  if (!cost_params_set_)
  {
    ROS_ERROR_STREAM("Cost parameters not set. use setCostsParams to set the parameters before !");
  }

  Eigen::MatrixXd Q = alpha_*Q1_ + (1-alpha_)*Q2_;
  Eigen::MatrixXd R1 = alpha_*R1_;
  Eigen::MatrixXd R2 = (1-alpha_)*R2_;

  return distMPCGain(Q,R1,Q,R2,N_);
}


Eigen::MatrixXd DistMPC::distMPCGain( const Eigen::MatrixXd& Q1,
                                      const Eigen::MatrixXd& R1,
                                      const Eigen::MatrixXd& Q2,
                                      const Eigen::MatrixXd& R2,
                                      const int& N)
{
  Eigen::MatrixXd Qh_bar = blkdiag(Q1,N);
  Eigen::MatrixXd Rh_bar = blkdiag(R1,N);
  Eigen::MatrixXd Qr_bar = blkdiag(Q2,N);
  Eigen::MatrixXd Rr_bar = blkdiag(R2,N);
  
  Eigen::MatrixXd L1 = ( theta_.transpose() * Qh_bar * theta_ + Rh_bar ).inverse() * theta_.transpose() * Qh_bar;
  Eigen::MatrixXd L2 = ( theta_.transpose() * Qr_bar * theta_ + Rr_bar ).inverse() * theta_.transpose() * Qr_bar;
    
  Eigen::MatrixXd gamma_1(L1.rows(), L1.cols() + 4); gamma_1 << -L1*psi_, L1;
  Eigen::MatrixXd gamma_2(L2.rows(), L2.cols() + 4); gamma_2 << -L2*psi_, L2;
  
  Eigen::MatrixXd lambda_1 = L1*theta_;
  Eigen::MatrixXd lambda_2 = L2*theta_;

  Eigen::MatrixXd K1(2*N,2*N); K1.setZero();
  K1 << Eigen::MatrixXd::Identity(lambda_1.rows(), lambda_1.cols()), -lambda_1,
        -lambda_2, Eigen::MatrixXd::Identity(lambda_2.rows(), lambda_2.cols());        
  Eigen::MatrixXd K2(2*N,gamma_1.cols()+gamma_2.cols());K2.setZero();
  K2 << gamma_1, Eigen::MatrixXd::Zero(gamma_1.rows(), gamma_1.cols()),
        Eigen::MatrixXd::Zero(gamma_2.rows(), gamma_2.cols()), gamma_2;

  K_mpc_ = K1.inverse() *K2;
  
  gains_set_ = true;
  
  return K_mpc_;
}


void DistMPC::setAlpha(const double& alpha)
{
  if(alpha>1 || alpha <0)
  {
    ROS_ERROR_STREAM("weight alpha must be 0 < alpha < 1 . Current value of alpha: "<<alpha);
  }
  alpha_ = alpha;
}

bool DistMPC::setCurrentState(const Eigen::VectorXd& x)
{  
  if (x.size() != 2*n_dofs_)
  {
    ROS_ERROR_STREAM("State size is not correct. got: "<< x.size()<<", required: "<< 2*n_dofs_);
    return false;
  }

  X_ = x;
  state_ok_=true;
  return true;
}

bool DistMPC::setReference(const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r)
{
  if (!state_ok_)
  {
    ROS_ERROR("State not set properly! ");
    return false;
  }
  if(ref_h.size()<N_ || ref_r.size()<N_)
  {
    ROS_ERROR_STREAM("Reference length is not correct! required: " << N_ <<", got ref_h.size(): "<<ref_h.size()<<" and ref_r.size(): "<< ref_r.size() );
    return false;    
  }
  
  Eigen::VectorXd ref_vec; ref_vec.resize(2*ref_h.size());

  for(int i=0; i<ref_h.size(); i++ )
  {
    ref_vec(2*i) = ref_h(i);
    ref_vec(2*i+1) = ref_r(i);
  }

  reference_.segment(0,2*X_.size())<< X_,X_;  
  reference_.segment(2*X_.size(),ref_vec.size()) = ref_vec;
  reference_.segment(2*X_.size()+ref_vec.size(),(2*X_.size()+ref_vec.size())) = reference_.segment(0, N_+ref_vec.size());

  state_ok_ = false;

  return true;
}

bool DistMPC::getControlInputs(Eigen::VectorXd& u1, Eigen::VectorXd& u2)
{
  if(!gains_set_)
  {
      ROS_ERROR_STREAM("Gains are not yet computed. use distMPCGain or distCoopMPCGain functions before");
      return false;
  }

  Eigen::VectorXd U = K_mpc_*reference_;
    
  if(U.size()< 2*N_)
  {
    ROS_ERROR_STREAM("Control vector length is not correct! required: " << 2*N_ <<", got U.size(): "<<U.size() );
    return false;    
  }
  u1 = U.segment(0,n_dofs_);
  u2 = U.segment(N_,n_dofs_);
  return true;
}

Eigen::VectorXd DistMPC::controlInputs(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r)
{
  if (x.size() != 2*n_dofs_)
  {
    ROS_ERROR_STREAM("State size is not correct. got: "<< x.size()<<", required: "<< 2*n_dofs_);
  }
  setCurrentState(x);
  setReference(ref_h.segment(0,N_),ref_r.segment(0,N_));
  
  Eigen::VectorXd u1,u2;
  if (!getControlInputs(u1,u2))
    ROS_ERROR("something wrong in control inputs computation");
  
  Eigen::VectorXd ret; ret.resize(2*u1.size());
  ret << u1,
         u2;
  return ret;
}

Eigen::VectorXd DistMPC::step(const Eigen::VectorXd& x, const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r)
{
  
  if (x.size() != 2*n_dofs_)
  {
    ROS_ERROR_STREAM("State size is not correct. got: "<< x.size()<<", required: "<< 2*n_dofs_);
  }
  if (!sys_params_set_)
  {
    ROS_ERROR_STREAM("System parameters not set. use setC2DSysParams or setSysParams to set the parameters before !");
  }
  
  Eigen::VectorXd u = controlInputs(x,ref_h,ref_r);
  
  X_ = A_*X_ +B_*u.segment(0,n_dofs_)+B_*u.segment(n_dofs_,n_dofs_);

  ROS_INFO_STREAM("X: "<<X_.transpose());
  
  return X_;
}

Eigen::VectorXd DistMPC::step(const Eigen::VectorXd& ref_h, const Eigen::VectorXd& ref_r)
{
  return step(X_,ref_h,ref_r);
}

Eigen::VectorXd DistMPC::getCurrentState(){return X_;};


Eigen::MatrixXd DistMPC::blkdiag(const Eigen::MatrixXd& a, int count)
{
    Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(a.rows() * count, a.cols() * count);
    for (int i = 0; i < count; ++i)
    {
        bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
    }
    return bdm;
}


Eigen::MatrixXd DistMPC::Ylowtriangular(const Eigen::MatrixXd& A,const Eigen::MatrixXd& B,const Eigen::MatrixXd& C, const int& N)
{
  Eigen::MatrixXd ret; ret.resize(2*N,N); ret.setZero();
  Eigen::MatrixXd tmp;
  Eigen::VectorXd b;
  b.resize(N);

  for (size_t i=0;i<N;i++)
  {
    tmp=C;
    for(int j=0;j<i;j++)
      tmp = tmp * A;
    tmp = tmp * B;
    b(i)=tmp(0);
  }  
  
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<=i;j++)
    {
      ret.block(i*2,j,2,1)<< b(i-j),
                             b(i-j);
    }
  }
  
  return ret;
}


Eigen::MatrixXd DistMPC::YColMat(const Eigen::MatrixXd& A,const Eigen::MatrixXd& C, const int& N)
{
  
  Eigen::MatrixXd T; T.resize(2*N,4); T.setZero();

  Eigen::MatrixXd tmp;
  
  for(int i=0;i<N;i++)
  {
    tmp=C*A;
    for(int j=0;j<i;j++)
      tmp = tmp * A;
    
    T.block(2*i,0,2,4) = tmp;
  }
  return T;
}





