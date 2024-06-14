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
  double y, a_temp; arma::mat B1;
  const double pi=std::atan(1.0)*4;
  B1 = B(1,1);
  if (a11 < 0) {
    a_temp = 1;
  } else {
    a_temp = a11;
  }
  y = log(2*pi/B(1,1))/2-log_iwishart_InvA_const(h,B1)+
    (h-1)/2*log(a_temp)-(B(0,0)-pow(B(0,1),2)/B(1,1))*a11/2;
  //when a11 is negative this evaluate to nan
  return(y);
  
  //print test
  //Rcpp::Rcout << "a11: " << a11 << std::endl;
  //Rcpp::Rcout << "y:"  << y << std::endl;
  //Rcpp::Rcout << "a_temp:" << a_temp << std::endl;
  
}

// [[Rcpp::export]]
double log_H(double nu, arma::mat V, arma::mat Omega, int i, int j){
  // (i,j) = 0
  arma::mat Omega0, Ome12, Ome22, Omega0_ij;
  arma::mat c;
  Omega0 = Omega;
  Omega0(i-1,j-1) = 0;  Omega0(j-1,i-1) = 0;
  Ome12 = Omega0.row(j-1); 
  //Rcpp::Rcout << "Ome12: " << Ome12 << std::endl;
  Ome12.shed_cols(j-1,j-1);
  //Rcpp::Rcout << "Ome12: " << Ome12 << std::endl;
  Ome22 = Omega0;
  Ome22.shed_rows(j-1,j-1); Ome22.shed_cols(j-1,j-1);
  //Rcpp::Rcout << "Ome22: " << Ome22 << std::endl;
  c = Ome12*solve(Ome22,Ome12.t()); // check
  if(c(0,0) == 0) {
    c(0,0) = 0.0000001;
  }
  //Rcpp::Rcout << "c: " << c << std::endl;
  Omega0_ij << Omega(i-1,i-1) << 0 << arma::endr << 0 << c(0,0) << arma::endr;
  //Rcpp::Rcout << "Omega0_ij: " << Omega0_ij << std::endl;
  
  // (i,j) = 1, note j>i
  arma::mat Omega1_ij, Ome11, A, test;
  Ome12 = Omega;
  Ome12.swap_rows(i,j-1);  Ome12 = Ome12.rows(i-1,i);
  //Rcpp::Rcout << "Ome12: " << Ome12 << std::endl;
  Ome12.shed_cols(j-1,j-1); Ome12.shed_cols(i-1,i-1);
  //Rcpp::Rcout << "Ome12: " << Ome12 << std::endl;
  
  Ome22 = Omega;
  Ome22.shed_rows(j-1,j-1); Ome22.shed_rows(i-1,i-1);
  //Rcpp::Rcout << "Ome22: " << Ome22 << std::endl;
  Ome22.shed_cols(j-1,j-1); Ome22.shed_cols(i-1,i-1);
  //Rcpp::Rcout << "Ome22: " << Ome22 << std::endl;
  Omega1_ij = Ome12*solve(Ome22,Ome12.t()); //check
  if(Omega1_ij(1,1) == 0) {
    Omega1_ij(1,1) = 0.0000001;
  } 
  if(Omega1_ij(0,0) == 0) {
    Omega1_ij(0,0) = 0.0000001;
  } 
  
  Ome11 << Omega(i-1,i-1) << Omega(i-1,j-1)  << arma::endr << Omega(i-1,j-1)  << Omega(j-1,j-1)  << arma::endr;
  //Rcpp::Rcout << "Ome11: " << Ome11 << std::endl;
  A = Ome11-Omega1_ij;
  //Rcpp::Rcout << "A: " << A << std::endl;
  //test << 0 << 0 << arma::endr << 0 << 1 << arma::endr;
  
  //Need to fix log(det(Omega))
  double a11, f, det_0, det_1; arma::mat V_ij, V_jj;
  a11 = A(0,0); V_jj = V(j-1,j-1);
  V_ij << V(i-1,i-1) << V(i-1,j-1)  << arma::endr << V(i-1,j-1)  << V(j-1,j-1)  << arma::endr;
  // Rcpp::Rcout << "V_ij: " << V_ij << std::endl;
  // Rcpp::Rcout << "V_jj: " << V_jj << std::endl;
  det_0 = matdet(Omega0_ij);
  det_1 = matdet(Omega1_ij);
  //Handle case when det may be < 0 due to imputing small number in Omega0_ij, Omega1_ij;
  if(det_0 < 0) {
    det_0 = 0.001;
  }
  if(det_1 < 0) {
    det_1 = 0.001;
  }
  //check determinants
  // Rcpp::Rcout << "Omega0_ij: " << Omega0_ij << std::endl;
  // Rcpp::Rcout << "det_0 " << det_0 << std::endl;
  // Rcpp::Rcout << "Omega1_ij: " << Omega1_ij << std::endl;
  // Rcpp::Rcout << "det_1 " << det_1 << std::endl;
  f = -log_iwishart_InvA_const(nu,V_jj) - log_J(nu,V_ij,a11) + (nu-2)/2*(log(det_0) - log(det_1)) - trace(V_ij*(Omega0_ij-Omega1_ij))/2;
  
  // //print test
  // Rcpp::Rcout << "a11:"  << a11 << std::endl;
  // Rcpp::Rcout << "f:"  << f << std::endl;
  // Rcpp::Rcout << "Constant" << log_iwishart_InvA_const(nu,V_jj) << std::endl;
  // Rcpp::Rcout << "Log_J" << log_J(nu,V_ij,a11) << std::endl;
  // Rcpp::Rcout << "Trace" << trace(V_ij*(Omega0_ij-Omega1_ij))/2 << std::endl;
  
  
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
    Rcpp::Rcout << "c:"  << c(0, 0) << std::endl;
    Rcpp::Rcout << "b:"  << b << std::endl;
    Rcpp::Rcout << "D1:"  << D1 << std::endl;
    Rcpp::Rcout << "Ome22:"  << Ome22 << std::endl;
    Rcpp::Rcout << "Omega_new:"  << Omega_new << std::endl;
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
    
    Rcpp::Rcout << "A:"  << A << std::endl;
    
    double log_Joint, logK2by2, logKii;
    log_Joint = (b-2)/2*log(det(Omega))-trace(D*Omega)/2;
    
    arma::mat D_ij, V, Dii, Aii;
    D_ij  << D(i-1,i-1) << D(i-1,j-1) << arma::endr << D(i-1,j-1) << D(j-1,j-1) << arma::endr;
    logK2by2 = log_dWish(A,b,D_ij);
    
    Rcpp::Rcout << "D_ij:"  << D_ij << std::endl;
    V = inv(D_ij); Dii = 1/V(1,1); Aii=A(0,0); //error D_ij not invertible
    Rcpp::Rcout << "V"  << D_ij << std::endl;
    Rcpp::Rcout << "Dii:"  << Dii << std::endl;
    Rcpp::Rcout << "Aii:"  << Aii << std::endl;
    logKii = log_dWish(Aii,b+1,Dii);
    
    f = log_Joint+logKii-logK2by2;
  }
  return(f);
}
