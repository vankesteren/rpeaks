library(rpeaks)
library(rsleep)
library(tidyverse)
library(firatheme)
library(patchwork)

# test timing
seconds <- seq(1, 2000, length.out = 40)
rpeaks_t  <- numeric(40)
rsleep_t  <- numeric(40)
for (i in 21:40) {
  cat("Testing for seconds:", seconds[i], "\n")
  ecg <- rep(rsleep::example_ecg_200hz, seconds[i])
  t1 <- Sys.time()
  peaks <- rpeaks_pan_tompkins(ecg = ecg, sample_rate = 200)
  t2 <- Sys.time()
  peaks <- detect_rpeaks(ecg, 200, 5, 15, 1)
  t3 <- Sys.time()
  rpeaks_t[i] <- t2 - t1
  rsleep_t[i] <- t3 - t2
}

# create plots
p1 <- 
  tibble(secs = seconds, rsleep = rsleep_t, rpeaks = rpeaks_t) %>%
  pivot_longer(cols = c("rsleep", "rpeaks"), names_to = "method") %>% 
  ggplot(aes(x = secs, y = value*1000, color = method)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "ECG length (seconds)", y = "Processing time (ms)",
       colour = "") +
  theme_fira() +
  scale_colour_fira() + theme(legend.position = "top")

p2 <- 
  tibble(secs = seconds, better = rsleep_t / rpeaks_t) %>% 
  ggplot(aes(x = secs, y = better)) +
  geom_point(colour = firaCols[3]) +
  geom_smooth(method = "lm", se = FALSE, color = firaCols[3]) +
  labs(x = "ECG length (seconds)", y = "Speed improvement") +
  theme_fira() +
  scale_colour_fira()

# save plot
firaSave(
  plot     = p1 / p2,
  filename = "img/speedcomparison.png", 
  device   = "png", 
  width    = 7, 
  height   = 7, 
  dpi      = 300
)