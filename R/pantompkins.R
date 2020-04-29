## usethis namespace: start
#' @useDynLib rpeaks, .registration = TRUE
#' @importFrom Rcpp evalCpp
## usethis namespace: end
NULL

#' Detect rpeaks using the pan-tompkins algorithm
#' 
#' Fast implementation of the pan-tompkins algorithm. The default parameters are
#' taken from the original pan and tompkins paper
#' 
#' @param ecg raw ecg data vector
#' @param sample_rate sampling rate in Hz of the ecg
#' @param integration_window size of the integration window in seconds
#' @param refractory refractory period in seconds (minimum time between peaks)
#' 
#' @references Pan, J., & Tompkins, W. J. (1985). A real-time QRS detection algorithm. IEEE transactions on biomedical engineering, (3), 230-236.
#' 
#' @export
rpeaks_pan_tompkins <- function(ecg, sample_rate, integration_window = 0.15,
                                refractory = 0.2) {
  # Nyquist frequency
  nyq <- sample_rate / 2
  
  # expand ecg with 1 second at the start and the end
  ecg_expand <- 
    c(ecg[sample_rate:1], ecg, ecg[length(ecg):(length(ecg) - sample_rate)])
  
  # band-pass filter
  band <- signal::butter(n = 1, W = c(5 / nyq, 15 / nyq), type = "pass")
  ecg_filt <- signal::filter(band, ecg_expand)
  
  # smooth derivative filter
  deriv <- signal::sgolay(p = 3, n = 5, m = 1)
  ecg_deriv <- signal::filter(deriv, ecg_filt)
  
  # square the derivatives
  ecg_deriv_2 <- (ecg_deriv*ecg_deriv)[(sample_rate + 1):(sample_rate + length(ecg))]
  
  # fast integration (boxcar filter)
  integrator <- rep(1, round(integration_window*sample_rate))
  ecg_int <- fast_conv(ecg_deriv_2, integrator)
  
  # fast peak detection
  rpeak_idx <- detect_peaks(
    signal      = ecg_int, 
    lower_bound = mean(ecg_int)*3, 
    refractory  = round(sample_rate * refractory)# 200ms refractory period)
  )
  
  # return time in seconds peaks lag behind by integration window / 2
  (rpeak_idx - round(length(integrator) / 2)) * sample_rate
}