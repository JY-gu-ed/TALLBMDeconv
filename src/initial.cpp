#include <Rcpp.h>
using namespace Rcpp;
// Native symmetric-Dirichlet BayesPrism Gibbs kernel, C++ execution only.
// [[Rcpp::export]]
List tbm_bp_initial(NumericVector X,NumericMatrix phi,double alpha,IntegerVector keep){
 RNGScope scope;int K=phi.nrow(),G=phi.ncol(),L=keep.size(),next=0;
 NumericVector theta(K,1.0/K),nk(K),mean(K),sq(K);
 NumericMatrix Z(G,K),trace(L,K);
 std::vector<double> prob(K),gamma(K);std::vector<int> draws(K);
 for(int it=1;it<=keep[L-1];it++){
  std::fill(nk.begin(),nk.end(),0.0);bool save=next<L&&it==keep[next];
  for(int g=0;g<G;g++){
   double ps=0;for(int k=0;k<K;k++){prob[k]=phi(k,g)*theta[k];ps+=prob[k];}
   if(ps<=0 && X[g]==0)continue;
   if(ps<=0)stop("Zero assignment probability for expressed gene");
   for(int k=0;k<K;k++)prob[k]/=ps;
   R::rmultinom((int)std::round(X[g]),prob.data(),K,draws.data());
   for(int k=0;k<K;k++){nk[k]+=draws[k];if(save)Z(g,k)+=draws[k];}
  }
  long double total=0;for(int k=0;k<K;k++){gamma[k]=R::rgamma(nk[k]+alpha,1);total+=gamma[k];}
  double denom=(double)total;
  for(int k=0;k<K;k++)theta[k]=gamma[k]/denom;
  if(save){for(int k=0;k<K;k++){mean[k]+=theta[k];sq[k]+=theta[k]*theta[k];trace(next,k)=theta[k];}next++;}
  if(it%25==0)checkUserInterrupt();
 }
 NumericVector cv(K);
 for(int k=0;k<K;k++){mean[k]/=L;cv[k]=std::sqrt(std::max(0.0,sq[k]/L-mean[k]*mean[k]))/mean[k];for(int g=0;g<G;g++)Z(g,k)/=L;}
 return List::create(_["Z_n"]=Z,_["theta_n"]=mean,_["theta.cv_n"]=cv,_["theta_trace"]=trace);
}
