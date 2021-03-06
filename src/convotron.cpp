#include <Rcpp.h>
#include <fftw3.h>
using namespace Rcpp;

// Fast power of 2 function
int nextPowerOf2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}

// 1D FFT filtering
void fftfilt1(int wide, fftw_plan& p, 
	      fftw_complex* xwide, fftw_complex* ywide, bool corr)
{

    // FFTs
    p = fftw_plan_dft_1d(wide, xwide, xwide, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_execute_dft(p, ywide, ywide);
    
    // conj, followed by complex multiply
    for (int i = 0; i < wide; i++)
    {
        ywide[i][1] = -ywide[i][1];
        double xwr = xwide[i][0];
        double xwi = xwide[i][1];
        xwide[i][0] = xwr * ywide[i][0] - xwi * ywide[i][1];
        xwide[i][1] = xwr * ywide[i][1] + xwi * ywide[i][0];
    }

    // IFFT
    p = fftw_plan_dft_1d(wide, xwide, xwide, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

}

//' Valid region cross-correlation of a vector with another vector.
//' This is equivalent to convolve(x, y, type='filter') 
//'
//' @param x first input signal (length(x) >= length(y))
//' @param y second input signal (length(y) < length(x))
//' @return a vector that is of width length(x) - length(y) + 1
// [[Rcpp::export]]
NumericVector ctron_xcorr1(const NumericVector& x, const NumericVector& y)
{

    int iT = x.size();
    int fT = y.size();
    
    // valid-region cross-correlation output
    NumericVector z(iT - fT + 1, 0.);
    int oT = z.size();

    for (int i = 0; i < oT; ++i)
    {
    	for (int j = 0; j < fT; ++j)
    	{
    	    z[i] += x[i + j] * y[j];
    	}
    } 
    return z;
}

//' Valid region convolution of a vector with another vector
//' This is equivalent to convolve(x, rev(y), type='filter')
//'
//' @param x first input signal (length(x) >= length(y))
//' @param y second input signal (length(y) < length(x))
//' @return a vector that is of width length(x) - length(y) + 1
// [[Rcpp::export]]
NumericVector ctron_conv1(const NumericVector& x, const NumericVector& y)
{

    int iT = x.size();
    int fT = y.size();
    
    // Valid region convolution output
    NumericVector z(iT - fT + 1, 0.);
    int oT = z.size();

    // Same as xcorr by y is time-reversed
    for (int i = 0; i < oT; ++i)
    {
    	for (int j = 0; j < fT; ++j)
    	{
    	    z[i] += x[i + j] * y[fT - 1 - j];
    	}
    } 
    return z;
}



//' Valid region cross-correlation of a vector with another vector using fftw.
//' This is equivalent to convolve(x, y, type='filter') when corr=TRUE
//'
//' @param x first input signal (length(x) >= length(y))
//' @param y second input signal (length(y) < length(x))
//' @return a vector that is of width length(x) - length(y) + 1
// [[Rcpp::export]]
NumericVector ctron_fftfilt1(const NumericVector& x, const NumericVector& y, bool corr)
{

    int xsz = x.size();
    int ysz = y.size();

    // fft size
    int wide = nextPowerOf2(xsz + ysz - 1);

    // valid-region convolution/cross-correlation size
    int narrow = xsz - ysz + 1;
    NumericVector z(narrow);

    fftw_plan p;

    // Buffers for FFTs
    fftw_complex* xwide = 
	(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wide);
    fftw_complex* ywide =
	(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wide);

    int mn = std::min<int>(xsz, ysz);

    // Zero memory upfront
    for (int i = mn; i < wide; ++i)
    {
        xwide[i][0] = xwide[i][1] = ywide[i][0] = ywide[i][1] = 0.;
    }

    // copy x into a complex array
    for (int i = 0; i < xsz; ++i)
    {
        xwide[i][0] = x[i];
        xwide[i][1] = 0;
    }

    // copy y into a complex array
    if (corr)
    {
        for (int i = 0; i < ysz; ++i)
        {
            ywide[i][0] = y[i];
            ywide[i][1] = 0;
        }
    }
    else
    {
        for (int i = 0; i < ysz; ++i)
        {
            ywide[i][0] = y[ysz - i - 1];
            ywide[i][1] = 0;
        }
    }

    fftfilt1(wide, p, xwide, ywide, corr);

    // Copy to output
    for (int i = 0; i < narrow; ++i)
    {
        double re = xwide[i][0] / wide;
        z[i] = re;
    }

    // Cleanup
    fftw_destroy_plan(p);
    fftw_free(xwide);
    fftw_free(ywide);

    // Go home
    return z;
}

//' Valid region convolution/cross-correlation of a matrix with another
//' matrix using fftw.  Each matrix must have the same number of rows:
//' one for each convolution.
//'
//' @param x each row is an input signal
//' @param y each row is an input signal
//' @return a matrix that is of width length(x[i,]) - length(y[i,]) + 1
// [[Rcpp::export]]
NumericMatrix ctron_fftfilt1mx(const NumericMatrix& x, const NumericMatrix& y, bool corr)
{

    int xsz = x.ncol();
    int ysz = y.ncol();
    int rows = x.rows();
    
    // fft size
    int wide = nextPowerOf2(xsz + ysz - 1);

    // valid-region convolution/cross-correlation size
    int narrow = xsz - ysz + 1;
    NumericMatrix z(rows, narrow);

    fftw_plan p, pinv;

    // Buffers for FFTs
    int doubleWide = 2 * sizeof(fftw_complex) * wide;
    fftw_complex* xy = (fftw_complex*) fftw_malloc(doubleWide);

    int rank = 1;
    int n[] = {wide};
    int howmany = 2;
    int istride = 1;
    int ostride = 1;
    int idist = wide;
    int odist = wide;
    int *inembed = n;
    int *onembed = n;
    
    p = fftw_plan_many_dft(rank, n, howmany,
			   xy, inembed,
			   istride, idist,
			   xy, onembed,
			   ostride, odist,
			   FFTW_FORWARD, FFTW_ESTIMATE);

    howmany = 1;
    
    pinv = fftw_plan_many_dft(rank, n, howmany,
			      xy, inembed,
			      istride, idist,
			      xy, onembed,
			      ostride, odist,
			      FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int k = 0; k < rows; ++k)
    {
	memset(xy, 0, doubleWide);
      
	// copy x into a complex array
	for (int i = 0; i < xsz; ++i)
	{
	    xy[i][0] = x(k, i);
	}

	// copy y into a complex array
	if (corr)
	{
	    for (int i = 0, j = wide; i < ysz; ++i, ++j)
	    {
		xy[j][0] = y(k, i);
	    }
	}
	else
	{
	    for (int i = 0, j = wide; i < ysz; ++i, ++j)
	    {
		xy[j][0] = y(k, ysz - i - 1);
	    }
	}

	fftw_execute(p);
	

	// conj, followed by complex multiply
	for (int i = 0, j = wide; i < wide; i++, j++)
	{
	    
	    xy[j][1] = -xy[j][1];
	    double xwr = xy[i][0];
	    double xwi = xy[i][1];
	    xy[i][0] = xwr * xy[j][0] - xwi * xy[j][1];
	    xy[i][1] = xwr * xy[j][1] + xwi * xy[j][0];
	}

	// IFFT
	//p = fftw_plan_dft_1d(wide, xy, xy, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(pinv);

	// Copy to output
	for (int i = 0; i < narrow; ++i)
	{
	    double re = xy[i][0] / wide;
	    z(k, i) = re;
	}
    }
    // Cleanup
    fftw_destroy_plan(p);
    fftw_destroy_plan(pinv);
    fftw_free(xy);

    // Go home
    return z;
}


//' Valid region cross-correlation of a matrix with another matrix.
//' Each matrix must have the same number of rows: one for each cross-corr.
//'
//' @param x each row is an input signal
//' @param y each row is an input signal
//' @return a matrix that is of width length(x[i,]) - length(y[i,]) + 1
// [[Rcpp::export]]
NumericMatrix ctron_xcorr1mx(const NumericMatrix& x, const NumericMatrix& y)
{
    int iT = x.ncol();
    int fT = y.ncol();
    int rows = x.nrow();

    // Should we perform filtering using the same kernel for all rows
    // or should we apply a different kernel
    int oT = iT - fT + 1;
    NumericMatrix z(rows, oT);

    for (int k = 0; k < rows; ++k)
    {
	for (int i = 0; i < oT; ++i)
	{
	    for (int j = 0; j < fT; ++j)
	    {
		z(k, i) += x(k, i + j) * y(k, j);
	    }
	}
    }

    return z;
}

//' Valid region convolution of a matrix with another matrix.
//' Each matrix must have the same number of rows: one for each convolution.
//'
//' @param x each row is an input signal
//' @param y each row is an input signal
//' @return a matrix that is of width length(x[i,]) - length(y[i,]) + 1
// [[Rcpp::export]]
NumericMatrix ctron_conv1mx(const NumericMatrix& x, const NumericMatrix& y)
{
    int iT = x.ncol();
    int fT = y.ncol();
    int rows = x.nrow();

    // Should we perform filtering using the same kernel for all rows
    // or should we apply a different kernel
    int oT = iT - fT + 1;
    NumericMatrix z(rows, oT);

    for (int k = 0; k < rows; ++k)
    {
	for (int i = 0; i < oT; ++i)
	{
	    for (int j = 0; j < fT; ++j)
	    {
		z(k, i) += x(k, i + j) * y(k, fT - 1 - j);
	    }
	}
    }
    return z;
}

//' Valid region convolution of a matrix with a single filter.
//'
//' @param x each row is an input signal
//' @param y is an input signal
//' @return a matrix that is of width length(x[i,]) - length(y) + 1
// [[Rcpp::export]]
NumericMatrix ctron_xcorr1mxv(const NumericMatrix& x, const NumericVector& y)
{
    int iT = x.ncol();
    int fT = y.size();
    int rows = x.nrow();

    // Should we perform filtering using the same kernel for all rows
    // or should we apply a different kernel
    int oT = iT - fT + 1;
    NumericMatrix z(rows, oT);

    for (int k = 0; k < rows; ++k)
    {
	for (int i = 0; i < oT; ++i)
	{
	    for (int j = 0; j < fT; ++j)
	    {
		z(k, i) += x(k, i + j) * y[j];
	    }
	}
    }
    return z;
}

//' Valid region convolution of a matrix with a single filter.
//'
//' @param x each row is an input signal
//' @param y is an input signal
//' @return a matrix that is of width length(x[i,]) - length(y) + 1
// [[Rcpp::export]]
NumericMatrix ctron_conv1mxv(const NumericMatrix& x, const NumericVector& y)
{
    int iT = x.ncol();
    int fT = y.size();
    int rows = x.nrow();

    // Should we perform filtering using the same kernel for all rows
    // or should we apply a different kernel
    int oT = iT - fT + 1;
    NumericMatrix z(rows, oT);

    for (int k = 0; k < rows; ++k)
    {
	for (int i = 0; i < oT; ++i)
	{
	    for (int j = 0; j < fT; ++j)
	    {
		z(k, i) += x(k, i + j) * y[fT - 1 - j];
	    }
	}
    }
    return z;
}

//' Valid region cross-correlation of a matrix with a filter matrix
//'
//' @param x each row is an input signal of dims (iR, iC)
//' @param y is an input signal of dims (fR, fC), fR <= iR, fC <= iC
//' @return a matrix that is of dims (iR - fR + 1, iC - fC + 1)
// [[Rcpp::export]]
NumericMatrix ctron_xcorr2mx(const NumericMatrix& x, const NumericMatrix& y)
{
    int iR = x.nrow();
    int iC = x.ncol();
    int fR = y.nrow();
    int fC = y.ncol();

    int oR = iR - fR + 1;
    int oC = iC - fC + 1;

    NumericMatrix z(oR, oC);

    for (int i = 0; i < oR; ++i)
    {
	for (int j = 0; j < oC; ++j)
	{

	    double acc = 0.;
	    for (int m = 0; m < fR; ++m)
	    {
		for (int n = 0; n < fC; ++n)
		{
		    acc += x(i + m, j + n) * y(m, n);
		}
	    }
	    z(i, j) = acc;
	}
    }
    return z;
}

//' Valid region convolution of a matrix with a filter matrix
//'
//' @param x each row is an input signal of dims (iR, iC)
//' @param y is an input signal of dims (fR, fC), fR <= iR, fC <= iC
//' @return a matrix that is of dims (iR - fR + 1, iC - fC + 1)
// [[Rcpp::export]]
NumericMatrix ctron_conv2mx(const NumericMatrix& x, const NumericMatrix& y)
{
    int iR = x.nrow();
    int iC = x.ncol();
    int fR = y.nrow();
    int fC = y.ncol();

    int oR = iR - fR + 1;
    int oC = iC - fC + 1;

    NumericMatrix z(oR, oC);

    for (int i = 0; i < oR; ++i)
    {
	for (int j = 0; j < oC; ++j)
	{

	    double acc = 0.;
	    for (int m = 0; m < fR; ++m)
	    {
		for (int n = 0; n < fC; ++n)
		{
		    acc += x(i + m, j + n) * y(fR - m - 1, fC - n - 1);
		}
	    }
	    z(i, j) = acc;
	}
    }
    return z;
}

