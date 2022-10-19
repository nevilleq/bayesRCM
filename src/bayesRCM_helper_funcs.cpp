#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]
arma::mat matinv(arma::mat x) {
  return( inv(x) ) ;
}

// [[Rcpp::export()]]
arma::mat matsolve(arma::mat x, arma::mat y) {
  return( solve(x,y) ) ;
}

// [[Rcpp::export()]]
arma::mat matprod(arma::mat x, arma::mat y) {
  return( x * y ) ;
}

// [[Rcpp::export()]]
arma::mat matABA(arma::mat x, arma::mat y) {
  return( x * y * x.t() ) ;
}

// [[Rcpp::export()]]
arma::mat matABinvA(arma::mat x, arma::mat y) {
  return( x * solve(y , x.t()) ) ;
}

// [[Rcpp::export()]]
double mattr(arma::mat x) {
  return( trace(x) ) ;
}

// [[Rcpp::export()]]
double matdet(arma::mat x) {
  return( det(x) ) ;
}

// [[Rcpp::export()]]
arma::mat matchol(arma::mat x) {
  return( chol(x) ) ;
}

// [[Rcpp::export()]]
double log_multi_gamma(int p, double n) {
  //  multivariate gamma function at value n,p
  //   \Gamma_p(n/2)= \pi^{p(p-1)/4}\Pi_{j=1}^p \Gamma\left[ (n+1-j)/2\right].
  double f;
  const double pi= std::atan(1.0)*4;
  f = (p*(p-1)/4)*log(pi);
  for( int j = 1; j <= p; j ++ ) {
    f = f + lgamma(n+(1-j)/2);
  }
  return(f);
}

// [[Rcpp::export]]
double log_iwishart_InvA_const(double df, arma::mat S){
  // Generates the normalizing constant "cons."  for a Inv-Wishart(df,S).
  // Nonsingular pdf is p(K) = cons. |K|^{-(df+2|p|)/2} exp(-trace(inv(K) S)/2)
  int p = S.n_rows;
  double iwc;
  iwc = (df+p-1)/2*(log(det(S))-p*log(2))-log_multi_gamma(p,(df+p-1)/2);
  return(iwc);
}

// [[Rcpp::export]]
double log_J(double h, arma::mat B, double a11){
  double y; arma::mat B1;
  const double pi=std::atan(1.0)*4;
  B1 = B(1,1);
  y =  log(2*pi/B(1,1))/2-log_iwishart_InvA_const(h,B1)+
    (h-1)/2*log(a11)-(B(0,0)-pow(B(0,1),2)/B(1,1))*a11/2;
  return(y);
}

// [[Rcpp::export]]
double log_H(double nu, arma::mat V, arma::mat Omega, int i, int j){
  // (i,j) = 0
  arma::mat Omega0, Ome12, Ome22, Omega0_ij;
  arma::mat c;
  Omega0 = Omega;
  Omega0(i-1,j-1) = 0;  Omega0(j-1,i-1) = 0;
  Ome12 = Omega0.row(j-1); Ome12.shed_cols(j-1,j-1);
  Ome22 = Omega0;
  Ome22.shed_rows(j-1,j-1); Ome22.shed_cols(j-1,j-1);
  c = Ome12*solve(Ome22,Ome12.t());
  Omega0_ij << Omega(i-1,i-1) << 0 << arma::endr << 0 << c(0,0) << arma::endr;
  
  // (i,j) = 1, note j>i
  arma::mat Omega1_ij, Ome11, A;
  Ome12 = Omega;
  Ome12.swap_rows(i,j-1);  Ome12 = Ome12.rows(i-1,i);
  Ome12.shed_cols(j-1,j-1); Ome12.shed_cols(i-1,i-1);
  
  Ome22 = Omega;
  Ome22.shed_rows(j-1,j-1); Ome22.shed_rows(i-1,i-1);
  Ome22.shed_cols(j-1,j-1); Ome22.shed_cols(i-1,i-1);
  Omega1_ij = Ome12*solve(Ome22,Ome12.t());
  
  Ome11 << Omega(i-1,i-1) << Omega(i-1,j-1)  << arma::endr << Omega(i-1,j-1)  << Omega(j-1,j-1)  << arma::endr;
  A = Ome11-Omega1_ij;
  
  double a11, f; arma::mat V_ij, V_jj;
  a11 = A(0,0); V_jj = V(j-1,j-1);
  V_ij << V(i-1,i-1) << V(i-1,j-1)  << arma::endr << V(i-1,j-1)  << V(j-1,j-1)  << arma::endr;
  f = -log_iwishart_InvA_const(nu,V_jj)-log_J(nu,V_ij,a11) + (nu-2)/2*(log(a11))  - trace(V_ij*(Omega0_ij-Omega1_ij))/2;
  
  return(f);
}

// [[Rcpp::export]]
double log_dWish(arma::mat Omega, double df, arma::mat S){
  return( (df-2)/2*log(det(Omega))-trace(S*Omega)/2+log_iwishart_InvA_const(df,S) );
}

// [[Rcpp::export]]
double log_GWish_NOij_pdf(double b, arma::mat D, arma::mat Omega, int i, int j, int edgeij){
  // Compute log p(Omega\omega(i,j) ) upto the normalizing constant of G-Wishart
  
  double f;
  
  if(edgeij ==0) {
    
    Omega(i-1,j-1) = 0; Omega(j-1,i-1)=0;
    
    arma::mat Ome12, Ome22, Omega_new, D1;
    Ome12 = Omega.row(j-1); Ome12.shed_cols(j-1,j-1);
    Ome22 = Omega;
    Ome22.shed_rows(j-1,j-1); Ome22.shed_cols(j-1,j-1);
    
    arma::mat c;
    c = Ome12*solve(Ome22,Ome12.t());
    Omega_new = Omega; Omega_new(j-1,j-1)=c(0,0);
    
    D1 = D(j-1,j-1);
    f = -log_iwishart_InvA_const(b,D1)+(b-2)/2*log(det(Ome22))-trace(D*Omega_new)/2;
    
  } else {
    
    arma::mat Ome11,Ome12, Ome22, A;
    Ome12 = Omega;
    Ome12.swap_rows(i,j-1);  Ome12 = Ome12.rows(i-1,i);
    Ome12.shed_cols(j-1,j-1); Ome12.shed_cols(i-1,i-1);
    
    Ome22 = Omega;
    Ome22.shed_rows(j-1,j-1); Ome22.shed_rows(i-1,i-1);
    Ome22.shed_cols(j-1,j-1); Ome22.shed_cols(i-1,i-1);
    
    Ome11 << Omega(i-1,i-1) << Omega(i-1,j-1)  << arma::endr << Omega(i-1,j-1)  << Omega(j-1,j-1)  << arma::endr;
    
    A = Ome11-Ome12*solve(Ome22,Ome12.t());
    
    double log_Joint, logK2by2, logKii;
    log_Joint = (b-2)/2*log(det(Omega))-trace(D*Omega)/2;
    
    arma::mat D_ij, V, Dii, Aii;
    D_ij  << D(i-1,i-1) << D(i-1,j-1) << arma::endr << D(i-1,j-1) << D(j-1,j-1) << arma::endr;
    logK2by2 = log_dWish(A,b,D_ij);
    
    V = inv(D_ij); Dii = 1/V(1,1); Aii=A(0,0);
    logKii = log_dWish(Aii,b+1,Dii);
    
    f = log_Joint+logKii-logK2by2;
  }
  return(f);
}
