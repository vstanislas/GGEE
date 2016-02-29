// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// SVD
List SVD(arma::mat x);
RcppExport SEXP GGEE_SVD(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    __result = Rcpp::wrap(SVD(x));
    return __result;
END_RCPP
}
// cbind1
NumericMatrix cbind1(NumericVector x, NumericVector y);
RcppExport SEXP GGEE_cbind1(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    __result = Rcpp::wrap(cbind1(x, y));
    return __result;
END_RCPP
}
// cbind2
NumericMatrix cbind2(NumericMatrix a, NumericVector b);
RcppExport SEXP GGEE_cbind2(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    __result = Rcpp::wrap(cbind2(a, b));
    return __result;
END_RCPP
}
// IntProd
List IntProd(NumericMatrix g1, NumericMatrix g2, int g1num, int g2num);
RcppExport SEXP GGEE_IntProd(SEXP g1SEXP, SEXP g2SEXP, SEXP g1numSEXP, SEXP g2numSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type g1(g1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type g2(g2SEXP);
    Rcpp::traits::input_parameter< int >::type g1num(g1numSEXP);
    Rcpp::traits::input_parameter< int >::type g2num(g2numSEXP);
    __result = Rcpp::wrap(IntProd(g1, g2, g1num, g2num));
    return __result;
END_RCPP
}
