# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Valid region cross-correlation of a vector with another vector.
#' This is equivalent to convolve(x, y, type='filter') 
#'
#' @param x first input signal (length(x) >= length(y))
#' @param y second input signal (length(y) < length(x))
#' @return a vector that is of width length(x) - length(y) + 1
ctron_xcorr1 <- function(x, y) {
    .Call('convotron_ctron_xcorr1', PACKAGE = 'convotron', x, y)
}

#' Valid region convolution of a vector with another vector
#' This is equivalent to convolve(x, rev(y), type='filter')
#'
#' @param x first input signal (length(x) >= length(y))
#' @param y second input signal (length(y) < length(x))
#' @return a vector that is of width length(x) - length(y) + 1
ctron_conv1 <- function(x, y) {
    .Call('convotron_ctron_conv1', PACKAGE = 'convotron', x, y)
}

#' Valid region cross-correlation of a vector with another vector using fftw.
#' This is equivalent to convolve(x, y, type='filter') when corr=TRUE
#'
#' @param x first input signal (length(x) >= length(y))
#' @param y second input signal (length(y) < length(x))
#' @return a vector that is of width length(x) - length(y) + 1
ctron_fftfilt1 <- function(x, y, corr) {
    .Call('convotron_ctron_fftfilt1', PACKAGE = 'convotron', x, y, corr)
}

#' Valid region convolution/cross-correlation of a matrix with another
#' matrix using fftw.  Each matrix must have the same number of rows:
#' one for each convolution.
#'
#' @param x each row is an input signal
#' @param y each row is an input signal
#' @return a matrix that is of width length(x[i,]) - length(y[i,]) + 1
ctron_fftfilt1mx <- function(x, y, corr) {
    .Call('convotron_ctron_fftfilt1mx', PACKAGE = 'convotron', x, y, corr)
}

#' Valid region cross-correlation of a matrix with another matrix.
#' Each matrix must have the same number of rows: one for each cross-corr.
#'
#' @param x each row is an input signal
#' @param y each row is an input signal
#' @return a matrix that is of width length(x[i,]) - length(y[i,]) + 1
ctron_xcorr1mx <- function(x, y) {
    .Call('convotron_ctron_xcorr1mx', PACKAGE = 'convotron', x, y)
}

#' Valid region convolution of a matrix with another matrix.
#' Each matrix must have the same number of rows: one for each convolution.
#'
#' @param x each row is an input signal
#' @param y each row is an input signal
#' @return a matrix that is of width length(x[i,]) - length(y[i,]) + 1
ctron_conv1mx <- function(x, y) {
    .Call('convotron_ctron_conv1mx', PACKAGE = 'convotron', x, y)
}

#' Valid region convolution of a matrix with a single filter.
#'
#' @param x each row is an input signal
#' @param y is an input signal
#' @return a matrix that is of width length(x[i,]) - length(y) + 1
ctron_xcorr1mxv <- function(x, y) {
    .Call('convotron_ctron_xcorr1mxv', PACKAGE = 'convotron', x, y)
}

#' Valid region convolution of a matrix with a single filter.
#'
#' @param x each row is an input signal
#' @param y is an input signal
#' @return a matrix that is of width length(x[i,]) - length(y) + 1
ctron_conv1mxv <- function(x, y) {
    .Call('convotron_ctron_conv1mxv', PACKAGE = 'convotron', x, y)
}

#' Valid region cross-correlation of a matrix with a filter matrix
#'
#' @param x each row is an input signal of dims (iR, iC)
#' @param y is an input signal of dims (fR, fC), fR <= iR, fC <= iC
#' @return a matrix that is of dims (iR - fR + 1, iC - fC + 1)
ctron_xcorr2mx <- function(x, y) {
    .Call('convotron_ctron_xcorr2mx', PACKAGE = 'convotron', x, y)
}

#' Valid region convolution of a matrix with a filter matrix
#'
#' @param x each row is an input signal of dims (iR, iC)
#' @param y is an input signal of dims (fR, fC), fR <= iR, fC <= iC
#' @return a matrix that is of dims (iR - fR + 1, iC - fC + 1)
ctron_conv2mx <- function(x, y) {
    .Call('convotron_ctron_conv2mx', PACKAGE = 'convotron', x, y)
}

