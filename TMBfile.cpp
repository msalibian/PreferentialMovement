#include <TMB.hpp>
#include<cmath>

template<class Type>
struct spde_t_smooth2{
  Eigen::SparseMatrix<Type> M0;        // G0 eqn (10) in Lindgren
  Eigen::SparseMatrix<Type> M1;        // G1 eqn (10) in Lindgren
  Eigen::SparseMatrix<Type> M2;        // G2 eqn (10) in Lindgren
  Eigen::SparseMatrix<Type> M3;
  spde_t_smooth2(SEXP x){  /* x = List passed from R */
    M0 = tmbutils::asSparseMatrix<Type>(getListElement(x,"M0"));
    M1 = tmbutils::asSparseMatrix<Type>(getListElement(x,"M1"));
    M2 = tmbutils::asSparseMatrix<Type>(getListElement(x,"M2"));
    M3 = tmbutils::asSparseMatrix<Type>(getListElement(x,"M3"));
  }
};

template<class Type>
Eigen::SparseMatrix <Type> Q_spde_smooth2(spde_t_smooth2<Type> spde, Type kappa){
  Type kappa_pow2 = kappa*kappa;
  Type kappa_pow4 = kappa_pow2*kappa_pow2;
  Type kappa_pow6 = kappa_pow4*kappa_pow2;
  return kappa_pow6*spde.M0 + Type(3.0)*kappa_pow4*spde.M1 + Type(3.0)*kappa_pow2*spde.M2 + spde.M3;
}
// define negative joint log-liklihood ([X,Y,S,beta]=[Y|S,X,beta][X|S,beta][S][beta])
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(tsim);     // Time points where X is simulated
  DATA_VECTOR(Y1);        // Observations taken in longitude. Must have same length as iobs.
  DATA_VECTOR(Y2);        // Observations taken in latitude. Must have same length as iobs.
  DATA_VECTOR(Y);         // Sampled temperatures
  DATA_IVECTOR(trackId);  // Vector which points to when new tracks start/end
  DATA_IVECTOR(meshidxloc); // INLA pointer
  DATA_VECTOR(Ind);         // Used to create mean vector
  
  
  
  PARAMETER_VECTOR(S);        // Random field
  PARAMETER_VECTOR(beta);     // Latent behaviour states
  
  DATA_STRUCT(spde,spde_t_smooth2);
  // Field parameters
  PARAMETER(mu);
  PARAMETER(log_papertau);
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  PARAMETER(alpha);
  PARAMETER(log_d);
  PARAMETER(log_sdbehav);
  // Transform log parameters
  Type kappa = exp(log_kappa);
  Type tau = exp(log_tau);
  Type d = exp(log_d);
  Type sdbehav = exp(log_sdbehav);
  // Create sigma (marginal variance)
  Type sigma = sqrt(1 / (8 * M_PI * exp(2*log_papertau) * exp(4*log_kappa)));
  Type range = sqrt(8*2) / kappa;
  
  // initiate i and j
  int i,j;
  // ans will be the resulting likelihood
  Type ans=0;
  // create mean vector
  vector<Type> muvec(S.size());
  // create sparse matrix Q for SPDE
  SparseMatrix<Type> Q = Q_spde_smooth2(spde,kappa);
  muvec=Ind*mu;
  // evaluate [S]
  // Negative log likelihood
  ans += SCALE(GMRF(Q), 1/exp(log_papertau))(S - mu);      
  // create objects
  matrix<Type> covcond(Y.size(),Y.size());
  vector<Type> muveccond(Y.size());
  vector<Type> dt(tsim.size()-1);
  vector<Type> diffx(tsim.size()-1);
  vector<Type> diffy(tsim.size()-1);
  vector<Type> gradx(Y.size());
  vector<Type> grady(Y.size());
  // calculate gradients and covariance/mean of [Y|S,X]
  for (i=0;i<Y.size();i++)
  {
    covcond(i,i)=tau*tau;
    gradx(i) =  ((S(meshidxloc(Y.size() + 2*i + 0)) - S(meshidxloc(i)))/(Y1(Y.size() + 2*i + 0) - Y1(i)));
    grady(i) =  ((S(meshidxloc(Y.size() + 2*i + 1)) - S(meshidxloc(i)))/(Y2(Y.size() + 2*i + 1) - Y2(i)));
    muveccond(i) = S(meshidxloc(i));
    for ( j=0;j<i;j++)
    {
      covcond(i,j)=0;
      covcond(j,i)=covcond(i,j);
    }
  }
  // evaluate [Y|S,X]
  density::MVNORM_t<Type> neg_log_density_cond(covcond);
  ans += neg_log_density_cond(Y-muveccond);
  // calculate time differences/scaled movement
  for(int i=0;i<dt.size();i++){
    dt(i) = tsim(i+1)-tsim(i);
    diffx(i) = (Y1(i+1) - Y1(i))/dt(i);
    diffy(i) = (Y2(i+1) - Y2(i))/dt(i);
  }
  // calculate [X|S]
  for(int k=0; k<(trackId.size()-1); k++){
    for(int i=(trackId(k)+1);i<(trackId(k+1)-1);i++){
      ans -= dnorm(Y1(i+1),(Y1(i) + (((1-beta(i+1))*diffx(i-1))-(beta(i+1)*(alpha*gradx(i)*S(meshidxloc(i)))))*dt(i)), d*sqrt(dt(i)),1);
      ans -= dnorm(Y2(i+1),(Y2(i) + (((1-beta(i+1))*diffy(i-1))-(beta(i+1)*(alpha*grady(i)*S(meshidxloc(i)))))*dt(i)), d*sqrt(dt(i)),1);
    }
    // Calculate [beta]
    ans -= dnorm(beta(trackId(k)), Type(0.5), Type(0.000001), 1);
    for(int i=(trackId(k));i<(trackId(k+1)-1);i++){
      ans -= dnorm(beta(i+1), beta(i), dt(i)*sdbehav, 1);
    }
  }
  // record sigma
  REPORT(sigma);
  REPORT(range);
  // return final likelihood
  return ans;
}