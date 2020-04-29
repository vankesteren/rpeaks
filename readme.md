<p align="center">
  <img src="./img/logo.png" width="200px"></img>
</p>

Fast implementation of Pan & Tompkins (1985). It is programmed relatively efficiently, using `Rcpp` and `RcppArmadillo` to process long ecg data much faster than alternatives. Default processing parameters are taken directly from the original paper.

This package gratefully borrows parts of the processing code from `rsleep::detect_rpeaks()`.

## Speed comparison

For an ecg of 100 seconds:
```r
long_ecg <- rep(rsleep::example_ecg_200hz, 10)
microbenchmark::microbenchmark(
  rsleep = rsleep::detect_rpeaks(long_ecg, 200, 5, 15, 1),
  rpeaks = rpeaks::rpeaks_pan_tompkins(long_ecg, 200)
)
```
```
#> Unit: milliseconds
#>   expr        min         lq       mean     median         uq        max neval
#> rsleep 326.356301 339.053551 352.451999 351.001601 362.584601 417.046201   100
#> rpeaks   2.424701   2.663251   3.236545   2.810151   2.982601   7.274501   100
```
`rpeaks` is more than 100 times as fast on average.

## Details
This algorithm uses a butterworth filter of order 1 for the band-pass step, and a 3rd-order length-5 Savitzky-Golay smoothing filter to compute the derivative of the band-passed signal. Peak detection on the preprocessed signal works in a simplified way: we take the first value above the lower bound (3 * the mean signal value) which is higher than its  neighbours, and not within the refractory period after the previous R peak.


### Reference
Pan, J., & Tompkins, W. J. (1985). A real-time QRS detection algorithm. _IEEE transactions on biomedical engineering, (3)_, 230-236.