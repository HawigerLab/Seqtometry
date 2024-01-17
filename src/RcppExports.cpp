// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// wks
double wks(const NumericVector& gene_ranks, const NumericVector& sig_idx, const NumericVector& aux);
RcppExport SEXP _Seqtometry_wks(SEXP gene_ranksSEXP, SEXP sig_idxSEXP, SEXP auxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type gene_ranks(gene_ranksSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sig_idx(sig_idxSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type aux(auxSEXP);
    rcpp_result_gen = Rcpp::wrap(wks(gene_ranks, sig_idx, aux));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Seqtometry_wks", (DL_FUNC) &_Seqtometry_wks, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Seqtometry(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
