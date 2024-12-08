// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// exchCpp
NumericMatrix exchCpp(int t, double alpha);
RcppExport SEXP _SA24204141_exchCpp(SEXP tSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(exchCpp(t, alpha));
    return rcpp_result_gen;
END_RCPP
}
// ar1Cpp
NumericMatrix ar1Cpp(int t, double alpha);
RcppExport SEXP _SA24204141_ar1Cpp(SEXP tSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(ar1Cpp(t, alpha));
    return rcpp_result_gen;
END_RCPP
}
// statCpp
NumericMatrix statCpp(int t, NumericVector a);
RcppExport SEXP _SA24204141_statCpp(SEXP tSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(statCpp(t, a));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204141_exchCpp", (DL_FUNC) &_SA24204141_exchCpp, 2},
    {"_SA24204141_ar1Cpp", (DL_FUNC) &_SA24204141_ar1Cpp, 2},
    {"_SA24204141_statCpp", (DL_FUNC) &_SA24204141_statCpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204141(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
