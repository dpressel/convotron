
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