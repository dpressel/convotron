
test_that('ctron_xcorr1() , convolve(), and ctron_fftfilt1() are same', {
    A = runif(50) * 2 - 1
    B = runif(5) * 2 - 1
    convolveOut = convolve(A, B, type='filter')
    xcorr1Out = ctron_xcorr1(A, B)
    fft1Out = ctron_fftfilt1(A, B, corr=TRUE)
    delta = sum(abs(convolveOut - xcorr1Out))
    expect_equal(0, delta, 1e-6)
    delta = sum(abs(convolveOut - fft1Out))
    expect_equal(0, delta, 1e-6)
})

test_that('ctron_conv1() , convolve(), and ctron_fftfilt1(corr=FALSE) are same', {
    A = runif(50) * 2 - 1
    B = runif(5) * 2 - 1
    convolveOut = convolve(A, rev(B), type='filter')
    conv1Out = ctron_conv1(A, B)
    fft1Out = ctron_fftfilt1(A, B, corr=FALSE)
    delta = sum(abs(convolveOut - conv1Out))
    expect_equal(0, delta, 1e-6)
    delta = sum(abs(convolveOut - fft1Out))
    expect_equal(0, delta, 1e-6)
})

# Do a 2-dimensional convolution in R, based on
# https://github.com/cran/WaveletCo/blob/master/R/conv2.R
# Modified to be over the valid region and return a real value
conv2 <- function(X, A) {
  X.pad <- matrix(0, ncol = NCOL(X) + NCOL(A)-1, nrow = NROW(X)+NROW(A)-1);
  X.pad[1:NROW(X), 1:NCOL(X)] <- X
  A.pad <- matrix(0, ncol = NCOL(X) + NCOL(A)-1, nrow = NROW(X)+NROW(A)-1);
  A.pad[1:NROW(A), 1:NCOL(A)] <- A
  
  X.fft <- fft(X.pad);
  A.fft <- fft(A.pad);
  M <- Re(fft(X.fft * A.fft, inverse = TRUE))/length(X.fft)
  

  N.row <- (NROW(A) + (0:(NROW(X)-NROW(A))))
  N.col <- (NCOL(A) + (0:(NCOL(X)-NCOL(A))))
  XC <- M[N.row, N.col]
 
  return(XC);

}

# Do xcorr on a matrix in a loop using R's builtin convolve function
loopconvR <- function(Dx, Kx) {
  
  height = dim(Dx)[1]
  dlen = dim(Dx)[2]
  klen = dim(Kx)[2] 
  narrow = dlen - klen + 1
  mx = matrix(nrow=height, ncol=narrow)
  for (i in 1:height) {
    a = convolve(Dx[i,], Kx[i,], type='filter')
    mx[i,] = a
  }
  return(mx)
}

test_that('ctron_xcorr1mx(), loopconvR() and ctron_fftfilt1mx() are same', {
    A = matrix(runif(300)*2-1, 6, 50)
    B = matrix(runif(30)*2-1, 6, 5)
    convolveOut = loopconvR(A, B)
    expect_equal(dim(convolveOut)[1], 6)
    expect_equal(dim(convolveOut)[2], 46)

    xcorr1Out = ctron_xcorr1mx(A, B)
    expect_equal(dim(xcorr1Out)[1], 6)
    expect_equal(dim(xcorr1Out)[2], 46)

    delta = sum(sum(abs(convolveOut - xcorr1Out)))
    expect_equal(0, delta, 1e-6)

    fft1Out = ctron_fftfilt1mx(A, B, corr=TRUE)
    expect_equal(dim(fft1Out)[1], 6)
    expect_equal(dim(fft1Out)[2], 46)

    delta = sum(sum(abs(convolveOut - fft1Out)))
    expect_equal(0, delta, 1e-6)
			     
})

test_that('ctron_conv2mx(), conv2() produce the same answer', {
    A = matrix(runif(3000)*2-1, 60, 50)
    B = matrix(runif(30)*2-1, 5, 6)

    conv2Out = ctron_conv2mx(A, B)
    conv2ROut = conv2(A, B)
    delta = sum(sum(abs(conv2ROut - conv2Out)))
    expect_equal(0, delta, 1e-6)		     


})