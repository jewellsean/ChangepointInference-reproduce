required_packages <- c("microbenchmark", "doParallel", "tidyverse", "gridExtra")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

require(ChangepointInference)
require(microbenchmark)
library(tidyverse)
library(gridExtra)
library(doParallel)
cl <- makeCluster(12)
registerDoParallel(cl)

set.seed(1)
sd <- 1
delta <- 1.5
window_size <- 50
dd <- NULL
for (length_series in c(100, 200, 300, 500, 1000, 2000)) {
  
  n_changepts <- 10 * round(log10(length_series))
  changepts <- unique(sort(sample(2:length_series, size = n_changepts, replace = FALSE)))
  n_changepts <- length(changepts)
  changepts_aug <- c(0, changepts, length_series)

  mean_i <- 0
  underlying_mean <- 0 * numeric(length_series)
    for (i in 1:n_changepts) {
      underlying_mean[(changepts_aug[i] + 1):changepts_aug[i + 1]] <- mean_i
      if (mean_i == 0) {
        mean_i <- delta
      } else {
        mean_i <- 0
      }
    }
    
  
  
  lam <-  log(length_series)
  
  df <- foreach (seed = 1:50, .combine = rbind, .packages=c('ChangepointInference', 'microbenchmark')) %dopar% {
    set.seed(seed)  
    y <- underlying_mean + rnorm(length_series, sd = sd)
      
    df_t <- NULL
    ff <- changepoint_estimates(y, "L0", lam)
    Kp <- max(c(length(ff$change_pts), 1))
      out <- microbenchmark("l0_est" = {changepoint_estimates(y, "L0", lam)},
                            "l0_inference" = {changepoint_inference(y, "L0-fixed", lam, window_size)},
                            "bs_est" = {changepoint_estimates(y, "BS", Kp)},
                            "bs_inference" = {changepoint_inference(y, "BS-adaptive-M", Kp)},
                            times = 1L)
      df_tmp <- data.frame(procedure = out$expr, time = out$time, seed = seed)
      df_t <- rbind(df_t, df_tmp)
      df_t
  }
  
  df$length <- length_series
  dd <- rbind(dd, df)
}
  

stopCluster(cl)


p <- dd %>% mutate(time = time / 1e9, length = as.factor(length),
procedure = factor(procedure, levels = c("bs_est", "l0_est", "bs_inference", "l0_inference"),
labels = c("Estimate changepoints with binary segmentation (Approaches 1-3)",
"Estimate changepoints with L0 segmentation (Approach 4)",
"Inference with binary segmentation (Approaches 1-3)",
"Inference with L0 segmentation (Approach 4)"
))) %>%
  ggplot(aes(x = length, y = time, color = procedure)) +
  geom_boxplot() +
  scale_y_log10() +
  xlab("Length of series T") +
  ylab("Time (s)") +
  labs(color="") +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_text(size=0)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

grid.arrange(p, nrow = 1)