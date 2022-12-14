// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// matinv
arma::mat matinv(arma::mat x);
RcppExport SEXP _bayesRCM_matinv(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matinv(x));
    return rcpp_result_gen;
END_RCPP
}
// matsolve
arma::mat matsolve(arma::mat x, arma::mat y);
RcppExport SEXP _bayesRCM_matsolve(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matsolve(x, y));
    return rcpp_result_gen;
END_RCPP
}
// matprod
arma::mat matprod(arma::mat x, arma::mat y);
RcppExport SEXP _bayesRCM_matprod(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matprod(x, y));
    return rcpp_result_gen;
END_RCPP
}
// matABA
arma::mat matABA(arma::mat x, arma::mat y);
RcppExport SEXP _bayesRCM_matABA(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matABA(x, y));
    return rcpp_result_gen;
END_RCPP
}
// matABinvA
arma::mat matABinvA(arma::mat x, arma::mat y);
RcppExport SEXP _bayesRCM_matABinvA(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(matABinvA(x, y));
    return rcpp_result_gen;
END_RCPP
}
// mattr
double mattr(arma::mat x);
RcppExport SEXP _bayesRCM_mattr(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mattr(x));
    return rcpp_result_gen;
END_RCPP
}
// matdet
double matdet(arma::mat x);
RcppExport SEXP _bayesRCM_matdet(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matdet(x));
    return rcpp_result_gen;
END_RCPP
}
// matchol
arma::mat matchol(arma::mat x);
RcppExport SEXP _bayesRCM_matchol(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matchol(x));
    return rcpp_result_gen;
END_RCPP
}
// log_multi_gamma
double log_multi_gamma(int p, double n);
RcppExport SEXP _bayesRCM_log_multi_gamma(SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(log_multi_gamma(p, n));
    return rcpp_result_gen;
END_RCPP
}
// log_iwishart_InvA_const
double log_iwishart_InvA_const(double df, arma::mat S);
RcppExport SEXP _bayesRCM_log_iwishart_InvA_const(SEXP dfSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(log_iwishart_InvA_const(df, S));
    return rcpp_result_gen;
END_RCPP
}
// log_J
double log_J(double h, arma::mat B, double a11);
RcppExport SEXP _bayesRCM_log_J(SEXP hSEXP, SEXP BSEXP, SEXP a11SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type a11(a11SEXP);
    rcpp_result_gen = Rcpp::wrap(log_J(h, B, a11));
    return rcpp_result_gen;
END_RCPP
}
// log_H
double log_H(double nu, arma::mat V, arma::mat Omega, int i, int j);
RcppExport SEXP _bayesRCM_log_H(SEXP nuSEXP, SEXP VSEXP, SEXP OmegaSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(log_H(nu, V, Omega, i, j));
    return rcpp_result_gen;
END_RCPP
}
// log_dWish
double log_dWish(arma::mat Omega, double df, arma::mat S);
RcppExport SEXP _bayesRCM_log_dWish(SEXP OmegaSEXP, SEXP dfSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(log_dWish(Omega, df, S));
    return rcpp_result_gen;
END_RCPP
}
// log_GWish_NOij_pdf
double log_GWish_NOij_pdf(double b, arma::mat D, arma::mat Omega, int i, int j, int edgeij);
RcppExport SEXP _bayesRCM_log_GWish_NOij_pdf(SEXP bSEXP, SEXP DSEXP, SEXP OmegaSEXP, SEXP iSEXP, SEXP jSEXP, SEXP edgeijSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type edgeij(edgeijSEXP);
    rcpp_result_gen = Rcpp::wrap(log_GWish_NOij_pdf(b, D, Omega, i, j, edgeij));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayesRCM_matinv", (DL_FUNC) &_bayesRCM_matinv, 1},
    {"_bayesRCM_matsolve", (DL_FUNC) &_bayesRCM_matsolve, 2},
    {"_bayesRCM_matprod", (DL_FUNC) &_bayesRCM_matprod, 2},
    {"_bayesRCM_matABA", (DL_FUNC) &_bayesRCM_matABA, 2},
    {"_bayesRCM_matABinvA", (DL_FUNC) &_bayesRCM_matABinvA, 2},
    {"_bayesRCM_mattr", (DL_FUNC) &_bayesRCM_mattr, 1},
    {"_bayesRCM_matdet", (DL_FUNC) &_bayesRCM_matdet, 1},
    {"_bayesRCM_matchol", (DL_FUNC) &_bayesRCM_matchol, 1},
    {"_bayesRCM_log_multi_gamma", (DL_FUNC) &_bayesRCM_log_multi_gamma, 2},
    {"_bayesRCM_log_iwishart_InvA_const", (DL_FUNC) &_bayesRCM_log_iwishart_InvA_const, 2},
    {"_bayesRCM_log_J", (DL_FUNC) &_bayesRCM_log_J, 3},
    {"_bayesRCM_log_H", (DL_FUNC) &_bayesRCM_log_H, 5},
    {"_bayesRCM_log_dWish", (DL_FUNC) &_bayesRCM_log_dWish, 3},
    {"_bayesRCM_log_GWish_NOij_pdf", (DL_FUNC) &_bayesRCM_log_GWish_NOij_pdf, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesRCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
