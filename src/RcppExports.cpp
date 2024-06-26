// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gradN
Rcpp::NumericVector gradN(Function Model, List Data, NumericVector par, double h, int order);
RcppExport SEXP _YABS_gradN(SEXP ModelSEXP, SEXP DataSEXP, SEXP parSEXP, SEXP hSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type Model(ModelSEXP);
    Rcpp::traits::input_parameter< List >::type Data(DataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(gradN(Model, Data, par, h, order));
    return rcpp_result_gen;
END_RCPP
}
// harmwg
SEXP harmwg(Function Model, List Data, int Iterations, int Status, int Thinning, double ACC, NumericMatrix DevianceMat, int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples, NumericMatrix PPD, int Adapt, NumericMatrix Sigma);
RcppExport SEXP _YABS_harmwg(SEXP ModelSEXP, SEXP DataSEXP, SEXP IterationsSEXP, SEXP StatusSEXP, SEXP ThinningSEXP, SEXP ACCSEXP, SEXP DevianceMatSEXP, SEXP LIVSEXP, SEXP MonitorSEXP, SEXP Mo0SEXP, SEXP samplesSEXP, SEXP PPDSEXP, SEXP AdaptSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type Model(ModelSEXP);
    Rcpp::traits::input_parameter< List >::type Data(DataSEXP);
    Rcpp::traits::input_parameter< int >::type Iterations(IterationsSEXP);
    Rcpp::traits::input_parameter< int >::type Status(StatusSEXP);
    Rcpp::traits::input_parameter< int >::type Thinning(ThinningSEXP);
    Rcpp::traits::input_parameter< double >::type ACC(ACCSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DevianceMat(DevianceMatSEXP);
    Rcpp::traits::input_parameter< int >::type LIV(LIVSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Monitor(MonitorSEXP);
    Rcpp::traits::input_parameter< List >::type Mo0(Mo0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PPD(PPDSEXP);
    Rcpp::traits::input_parameter< int >::type Adapt(AdaptSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(harmwg(Model, Data, Iterations, Status, Thinning, ACC, DevianceMat, LIV, Monitor, Mo0, samples, PPD, Adapt, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// harm
SEXP harm(Function Model, List Data, int Iterations, int Status, int Thinning, double ACC, NumericMatrix DevianceMat, int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples, NumericMatrix PPD, int Adapt, NumericMatrix Sigma);
RcppExport SEXP _YABS_harm(SEXP ModelSEXP, SEXP DataSEXP, SEXP IterationsSEXP, SEXP StatusSEXP, SEXP ThinningSEXP, SEXP ACCSEXP, SEXP DevianceMatSEXP, SEXP LIVSEXP, SEXP MonitorSEXP, SEXP Mo0SEXP, SEXP samplesSEXP, SEXP PPDSEXP, SEXP AdaptSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type Model(ModelSEXP);
    Rcpp::traits::input_parameter< List >::type Data(DataSEXP);
    Rcpp::traits::input_parameter< int >::type Iterations(IterationsSEXP);
    Rcpp::traits::input_parameter< int >::type Status(StatusSEXP);
    Rcpp::traits::input_parameter< int >::type Thinning(ThinningSEXP);
    Rcpp::traits::input_parameter< double >::type ACC(ACCSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DevianceMat(DevianceMatSEXP);
    Rcpp::traits::input_parameter< int >::type LIV(LIVSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Monitor(MonitorSEXP);
    Rcpp::traits::input_parameter< List >::type Mo0(Mo0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PPD(PPDSEXP);
    Rcpp::traits::input_parameter< int >::type Adapt(AdaptSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(harm(Model, Data, Iterations, Status, Thinning, ACC, DevianceMat, LIV, Monitor, Mo0, samples, PPD, Adapt, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// gcharm
SEXP gcharm(Function Model, List Data, int Iterations, int Status, int Thinning, double ACC, NumericMatrix DevianceMat, double h, int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples, NumericMatrix PPD, int Adapt, NumericMatrix Sigma);
RcppExport SEXP _YABS_gcharm(SEXP ModelSEXP, SEXP DataSEXP, SEXP IterationsSEXP, SEXP StatusSEXP, SEXP ThinningSEXP, SEXP ACCSEXP, SEXP DevianceMatSEXP, SEXP hSEXP, SEXP LIVSEXP, SEXP MonitorSEXP, SEXP Mo0SEXP, SEXP samplesSEXP, SEXP PPDSEXP, SEXP AdaptSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type Model(ModelSEXP);
    Rcpp::traits::input_parameter< List >::type Data(DataSEXP);
    Rcpp::traits::input_parameter< int >::type Iterations(IterationsSEXP);
    Rcpp::traits::input_parameter< int >::type Status(StatusSEXP);
    Rcpp::traits::input_parameter< int >::type Thinning(ThinningSEXP);
    Rcpp::traits::input_parameter< double >::type ACC(ACCSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DevianceMat(DevianceMatSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type LIV(LIVSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Monitor(MonitorSEXP);
    Rcpp::traits::input_parameter< List >::type Mo0(Mo0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PPD(PPDSEXP);
    Rcpp::traits::input_parameter< int >::type Adapt(AdaptSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(gcharm(Model, Data, Iterations, Status, Thinning, ACC, DevianceMat, h, LIV, Monitor, Mo0, samples, PPD, Adapt, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// ohss
SEXP ohss(Function Model, List Data, int Iterations, int Status, int Thinning, double ACC, NumericMatrix DevianceMat, int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples, NumericMatrix PPD, int Adapt);
RcppExport SEXP _YABS_ohss(SEXP ModelSEXP, SEXP DataSEXP, SEXP IterationsSEXP, SEXP StatusSEXP, SEXP ThinningSEXP, SEXP ACCSEXP, SEXP DevianceMatSEXP, SEXP LIVSEXP, SEXP MonitorSEXP, SEXP Mo0SEXP, SEXP samplesSEXP, SEXP PPDSEXP, SEXP AdaptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type Model(ModelSEXP);
    Rcpp::traits::input_parameter< List >::type Data(DataSEXP);
    Rcpp::traits::input_parameter< int >::type Iterations(IterationsSEXP);
    Rcpp::traits::input_parameter< int >::type Status(StatusSEXP);
    Rcpp::traits::input_parameter< int >::type Thinning(ThinningSEXP);
    Rcpp::traits::input_parameter< double >::type ACC(ACCSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type DevianceMat(DevianceMatSEXP);
    Rcpp::traits::input_parameter< int >::type LIV(LIVSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Monitor(MonitorSEXP);
    Rcpp::traits::input_parameter< List >::type Mo0(Mo0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PPD(PPDSEXP);
    Rcpp::traits::input_parameter< int >::type Adapt(AdaptSEXP);
    rcpp_result_gen = Rcpp::wrap(ohss(Model, Data, Iterations, Status, Thinning, ACC, DevianceMat, LIV, Monitor, Mo0, samples, PPD, Adapt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_YABS_gradN", (DL_FUNC) &_YABS_gradN, 5},
    {"_YABS_harmwg", (DL_FUNC) &_YABS_harmwg, 14},
    {"_YABS_harm", (DL_FUNC) &_YABS_harm, 14},
    {"_YABS_gcharm", (DL_FUNC) &_YABS_gcharm, 15},
    {"_YABS_ohss", (DL_FUNC) &_YABS_ohss, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_YABS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
