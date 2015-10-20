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


void fftfilt1(int wide, fftw_plan& p, fftw_complex* xwide, fftw_complex* ywide, bool corr)
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

    fftw_plan p;

    // Buffers for FFTs
    fftw_complex* xwide = 
	(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wide);
    fftw_complex* ywide =
	(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wide);

    int mn = std::min<int>(xsz, ysz);


    for (int k = 0; k < rows; ++k)
    {

	// Zero memory upfront
	for (int i = mn; i < wide; ++i)
	{
	    xwide[i][0] = xwide[i][1] = ywide[i][0] = ywide[i][1] = 0.;
	}

	// copy x into a complex array
	for (int i = 0; i < xsz; ++i)
	{
	    xwide[i][0] = x(k, i);
	    xwide[i][1] = 0;
	}

	// copy y into a complex array
	if (corr)
	{
	    for (int i = 0; i < ysz; ++i)
	    {
		ywide[i][0] = y(k, i);
		ywide[i][1] = 0;
	    }
	}
	else
	{
	    for (int i = 0; i < ysz; ++i)
	    {
		ywide[i][0] = y(k, ysz - i - 1);
	    }
	}

	fftfilt1(wide, p, xwide, ywide, corr);

	// Copy to output
	for (int i = 0; i < narrow; ++i)
	{
	    double re = xwide[i][0] / wide;
	    z(k, i) = re;
	}
    }
    // Cleanup
    fftw_destroy_plan(p);
    fftw_free(xwide);
    fftw_free(ywide);

    // Go home
    return z;
}


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

    if (y.nrow() == rows)
    {
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
    }
    else
    {
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
    }
    return z;
}

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

    if (y.nrow() == rows)
    {
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
    }
    else
    {
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
    }
    return z;
}
