#include <Rcpp.h>
using namespace Rcpp;

// Project adapter: native fixed-phi BayesPrism Gibbs conditional, corrected
// by a Metropolis accept/reject for a typed soft observation likelihood.
// Never convert a lower bound to a point or post-process a finished theta.
double observation_loss(const NumericVector& theta,const NumericVector& sizes,
                        const NumericMatrix& groups,const NumericVector& lo,
                        const NumericVector& hi,const NumericVector& scales,double strength){
  double total=0;int K=theta.size();for(int k=0;k<K;k++)total+=theta[k]/sizes[k];
  double value=0;
  for(int h=0;h<groups.nrow();h++){
    double q=0;for(int k=0;k<K;k++)q+=groups(h,k)*theta[k]/sizes[k]/total;
    double z=(std::max(lo[h]-q,0.0)+std::max(q-hi[h],0.0))/scales[h];
    value+=z<=1 ? .5*z*z : z-.5;
  }return strength*value;
}

double pair_step_sd(const NumericVector& theta,const NumericVector& sizes,
                   const NumericMatrix& groups,const NumericVector& scales,
                   int i,int j,double ai,double aj,double strength){
 double total=0;int K=theta.size();for(int k=0;k<K;k++)total+=theta[k]/sizes[k];
 double s=theta[i]+theta[j],u=theta[i]/s,dt=s*u*(1-u),dtotal=dt*(1/sizes[i]-1/sizes[j]);
 double curvature=ai*aj/(ai+aj)+1;
 for(int h=0;h<groups.nrow();h++){
  double num=0;for(int k=0;k<K;k++)num+=groups(h,k)*theta[k]/sizes[k];
  double dn=dt*(groups(h,i)/sizes[i]-groups(h,j)/sizes[j]);
  double dq=(dn*total-num*dtotal)/(total*total);
  curvature+=strength*dq*dq/(scales[h]*scales[h]);
 }
 return std::min(1.0,2/std::sqrt(curvature));
}

// [[Rcpp::export]]
List tbm_flow_sample(NumericVector X,NumericMatrix phi,double alpha,
 IntegerVector keep,NumericVector sizes,NumericMatrix groups,NumericVector lo,
 NumericVector hi,NumericVector scales,double strength){
 RNGScope scope;int K=phi.nrow(),G=phi.ncol(),L=keep.size(),next=0,accepted=0,local_accepted=0,local_proposed=0;
 if(sizes.size()!=K||groups.ncol()!=K||groups.nrow()!=lo.size()||lo.size()!=hi.size()||hi.size()!=scales.size())stop("Observation dimension mismatch");
 for(int k=0;k<K;k++)if(!R_finite(sizes[k])||sizes[k]<=0)stop("Invalid gene-matched cell size");
 NumericVector theta(K,1.0/K),nk(K),mean(K),sq(K),proposal(K);
 NumericMatrix Z(G,K),trace(L,K);std::vector<double> prob(K);std::vector<int> draws(K);
 double loss=observation_loss(theta,sizes,groups,lo,hi,scales,strength);
 for(int it=1;it<=keep[L-1];it++){
  std::fill(nk.begin(),nk.end(),0.0);bool save=next<L&&it==keep[next];
  for(int g=0;g<G;g++){
   if(X[g]<0||X[g]>2147483647.0||std::abs(X[g]-std::round(X[g]))>1e-6)stop("Noninteger mixture counts");
   double ps=0;for(int k=0;k<K;k++){prob[k]=phi(k,g)*theta[k];ps+=prob[k];}
   if(ps<=0)stop("Zero assignment probability");
   for(int k=0;k<K;k++)prob[k]/=ps;
   R::rmultinom((int)std::round(X[g]),prob.data(),K,draws.data());
   for(int k=0;k<K;k++){nk[k]+=draws[k];if(save)Z(g,k)+=draws[k];}
  }
  long double total=0;for(int k=0;k<K;k++){proposal[k]=R::rgamma(nk[k]+alpha,1);total+=proposal[k];}
  if(total<=0)stop("Invalid Dirichlet draw");
  for(int k=0;k<K;k++)proposal[k]/=(double)total;
  double ploss=observation_loss(proposal,sizes,groups,lo,hi,scales,strength);
  // Independence proposal is the native full conditional, so the posterior
  // density/proposal terms cancel; only observation likelihood remains.
  bool positive=true;for(int k=0;k<K;k++)if(proposal[k]<=0)positive=false;
  bool accept=(strength==0||groups.nrow()==0||positive)&&(ploss<=loss||std::log(R::runif(0,1))<loss-ploss);
  if(accept){std::copy(proposal.begin(),proposal.end(),theta.begin());loss=ploss;accepted++;}
  // Exact pairwise log-ratio Metropolis sweeps prevent an RNA-only Dirichlet
  // proposal from freezing under a strong protein likelihood. This changes
  // proposal efficiency, not the target. The state-dependent proposal's
  // reverse density and the logit Jacobian are both included explicitly.
  if(strength>0&&groups.nrow()>0){
   for(int i=0;i<K;i++){
    int j=(int)std::floor(R::runif(0,K-1));if(j>=i)j++;
    double s=theta[i]+theta[j];if(s<=0||theta[i]<=0||theta[j]<=0)continue;
    double u=theta[i]/s;if(u<=0||u>=1)continue;
    double v=std::log(u)-std::log1p(-u),ai=nk[i]+alpha,aj=nk[j]+alpha;
    double sd=pair_step_sd(theta,sizes,groups,scales,i,j,ai,aj,strength);
    double vn=v+R::rnorm(0,sd);local_proposed++;
    if(vn< -690||vn>34)continue;
    double un=vn>=0?1/(1+std::exp(-vn)):std::exp(vn)/(1+std::exp(vn));
    if(un<=0||un>=1)continue;
    NumericVector candidate=clone(theta);candidate[i]=s*un;candidate[j]=s*(1-un);
    double nl=observation_loss(candidate,sizes,groups,lo,hi,scales,strength);
    double sn=pair_step_sd(candidate,sizes,groups,scales,i,j,ai,aj,strength),dv=vn-v;
    double logr=ai*(std::log(un)-std::log(u))+aj*(std::log1p(-un)-std::log1p(-u))-(nl-loss)
                +std::log(sd/sn)+.5*dv*dv*(1/(sd*sd)-1/(sn*sn));
    if(logr>=0||std::log(R::runif(0,1))<logr){std::copy(candidate.begin(),candidate.end(),theta.begin());loss=nl;local_accepted++;}
   }
  }
  if(save){for(int k=0;k<K;k++){mean[k]+=theta[k];sq[k]+=theta[k]*theta[k];trace(next,k)=theta[k];}next++;}
  if(it%25==0)checkUserInterrupt();
 }
 NumericVector cv(K);double maxerr=0;
 for(int k=0;k<K;k++){mean[k]/=L;cv[k]=mean[k]>0?std::sqrt(std::max(0.0,sq[k]/L-mean[k]*mean[k]))/mean[k]:NA_REAL;for(int g=0;g<G;g++)Z(g,k)/=L;}
 for(int g=0;g<G;g++){double total=0;for(int k=0;k<K;k++)total+=Z(g,k);maxerr=std::max(maxerr,std::abs(total-X[g]));}
 return List::create(_["Z_n"]=Z,_["theta_n"]=mean,_["theta.cv_n"]=cv,_["theta_trace"]=trace,
 _["acceptance_rate"]=(double)accepted/keep[L-1],_["max_gene_read_conservation_error"]=maxerr,
 _["final_observation_loss"]=loss,_["n_observations"]=groups.nrow(),
 _["local_proposal_acceptance_rate"]=local_proposed>0?(double)local_accepted/local_proposed:NA_REAL,
 _["local_proposals"]=local_proposed);
}
