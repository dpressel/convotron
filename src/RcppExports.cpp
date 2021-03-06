// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ctron_xcorr1
NumericVector ctron_xcorr1(const NumericVector& x, const NumericVector& y);
RcppExport SEXP convotron_ctron_xcorr1(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_xcorr1(x, y));
    return __result;
END_RCPP
}
// ctron_conv1
NumericVector ctron_conv1(const NumericVector& x, const NumericVector& y);
RcppExport SEXP convotron_ctron_conv1(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_conv1(x, y));
    return __result;
END_RCPP
}
// ctron_fftfilt1
NumericVector ctron_fftfilt1(const NumericVector& x, const NumericVector& y, bool corr);
RcppExport SEXP convotron_ctron_fftfilt1(SEXP xSEXP, SEXP ySEXP, SEXP corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type corr(corrSEXP);
    __result = Rcpp::wrap(ctron_fftfilt1(x, y, corr));
    return __result;
END_RCPP
}
// ctron_fftfilt1mx
NumericMatrix ctron_fftfilt1mx(const NumericMatrix& x, const NumericMatrix& y, bool corr);
RcppExport SEXP convotron_ctron_fftfilt1mx(SEXP xSEXP, SEXP ySEXP, SEXP corrSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type corr(corrSEXP);
    __result = Rcpp::wrap(ctron_fftfilt1mx(x, y, corr));
    return __result;
END_RCPP
}
// ctron_xcorr1mx
NumericMatrix ctron_xcorr1mx(const NumericMatrix& x, const NumericMatrix& y);
RcppExport SEXP convotron_ctron_xcorr1mx(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_xcorr1mx(x, y));
    return __result;
END_RCPP
}
// ctron_conv1mx
NumericMatrix ctron_conv1mx(const NumericMatrix& x, const NumericMatrix& y);
RcppExport SEXP convotron_ctron_conv1mx(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_conv1mx(x, y));
    return __result;
END_RCPP
}
// ctron_xcorr1mxv
NumericMatrix ctron_xcorr1mxv(const NumericMatrix& x, const NumericVector& y);
RcppExport SEXP convotron_ctron_xcorr1mxv(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_xcorr1mxv(x, y));
    return __result;
END_RCPP
}
// ctron_conv1mxv
NumericMatrix ctron_conv1mxv(const NumericMatrix& x, const NumericVector& y);
RcppExport SEXP convotron_ctron_conv1mxv(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_conv1mxv(x, y));
    return __result;
END_RCPP
}
// ctron_xcorr2mx
NumericMatrix ctron_xcorr2mx(const NumericMatrix& x, const NumericMatrix& y);
RcppExport SEXP convotron_ctron_xcorr2mx(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_xcorr2mx(x, y));
    return __result;
END_RCPP
}
// ctron_conv2mx
NumericMatrix ctron_conv2mx(const NumericMatrix& x, const NumericMatrix& y);
RcppExport SEXP convotron_ctron_conv2mx(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type y(ySEXP);
    __result = Rcpp::wrap(ctron_conv2mx(x, y));
    return __result;
END_RCPP
}
