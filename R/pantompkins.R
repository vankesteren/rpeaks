## usethis namespace: start
#' @useDynLib rpeaks, .registration = TRUE
#' @importFrom Rcpp evalCpp
## usethis namespace: end
NULL

#' Detect R peaks in ecg data using the Pan-Tompkins algorithm
#' 
#' Fast implementation of the Pan-Tompkins algorithm. The default parameters are
#' taken from the original pan and tompkins paper.
#' 
#' @param ecg raw ecg data vector
#' @param sample_rate sampling rate in Hz of the ecg
#' @param integration_window size of the integration window in seconds
#' @param refractory refractory period in seconds (minimum time between peaks)
#' @param band_low lower bound of the band-pass filter in Hz
#' @param band_high upper bound of the band-pass filter in Hz
#' 
#' @details This algorithm uses a butterworth filter of order 1 for the 
#' band-pass step, and a 3rd-order length-5 Savitzky-Golay smoothing filter to 
#' compute the derivative of the band-passed signal. Peak detection on the 
#' preprocessed signal works in a simplified way: we take the first value above
#' the lower bound (3 * the mean signal value) which is higher than its 
#' neighbours, and not within the refractory period after the previous R peak.
#' 
#' @references Pan, J., & Tompkins, W. J. (1985). A real-time QRS detection 
#' algorithm. IEEE transactions on biomedical engineering, (3), 230-236.
#' 
#' @export
rpeaks_pan_tompkins <- function(ecg, sample_rate, integration_window = 0.15,
                                refractory = 0.2, band_low = 5, 
                                band_high = 15) {
  # Nyquist frequency
  nyq <- sample_rate / 2
  n <- length(ecg)
  
  # Expand ecg with 1 second at the start and the end
  ecg <- c(ecg[sample_rate:1], ecg, ecg[n:(n - sample_rate)])
  
  # band-pass filter
  bandpass <- signal::butter(
    n    = 1, 
    W    = c(band_low / nyq, band_high / nyq), 
    type = "pass"
  )
  ecg <- signal::filter(bandpass, ecg)
  
  # smooth derivative filter
  deriv <- signal::sgolay(p = 3, n = 5, m = 1)
  ecg <- signal::filter(deriv, ecg)
  
  # square the derivatives and remove the pre and post second
  ecg <- (ecg*ecg)[(sample_rate + 1):(sample_rate + n)]
  
  # fast integration of the signal using a boxcar (C++)
  integrator <- rep(1, round(integration_window*sample_rate))
  ecg <- fast_conv(ecg, integrator)
  
  # fast peak detection (C++)
  rpeak_idx <- detect_peaks(
    signal      = ecg, 
    lower_bound = mean(ecg)*3, 
    refractory  = round(sample_rate * refractory)
  )
  
  # return time in seconds -- peaks lag behind by integration window / 2
  (rpeak_idx - round(length(integrator) / 2)) / sample_rate
}