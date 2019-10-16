## This script can be used to reproduce the results in d_est.csv and d_vs_true.csv
## To replicate the exact results the following parameters need to be changed
## n = 2000
## n_changepts = 50 
## n_replicates = 100
## max_n_search_path_l0 = 500
## Uncomment the last two lines to save the output CSVs


## Load required packages
required_packages <- c("tidyverse")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")
require(ChangepointInference)
require(tidyverse)

## utility functions
nearest_changepoint <- function(pt, cps) {
  M <- length(pt)
  out <- numeric(M)
  for (i in 1:M) {
    out[i] <- cps[which.min(abs(pt[i] - cps))]
  }
  return(out)
}

true_changepts_df <- function(df, true_changepts) {
  est_cps <- df$change_pts
  d_tmp <- data.frame(true_changepts = true_changepts, seed = df$seed[1], delta = df$delta[1], window_size = df$window_size[1], type = df$type[1]) %>%
    mutate(nearest_est_cp = nearest_changepoint(true_changepts, est_cps),
    dist_nearest_est = abs(true_changepts - nearest_est_cp)) %>%
    left_join(df, by = c("seed", "delta", "window_size", "type", "nearest_est_cp" = "change_pts"))
  return(d_tmp)
}

### Simulation parameter setup
deltas <- c(0.00000000000000000001, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
sd <- 1
n_replicates <- 5
if (n_replicates != 100) {
  warning("Set n_replicates = 100 to replicate paper results!")
}
# underlying_mean skelton
set.seed(1112)
n <- 50 
if (n != 2000) {
  warning("Set n = 2000 to replicate paper results!")
}
n_changepts <- 2
if (n_changepts != 50) {
  warning("Set n_changepts = 50 to replicate paper results!")
}
changepts <- unique(sort(sample(2:n, size = n_changepts, replace = FALSE)))
n_changepts <- length(changepts)
changepts_aug <- c(0, changepts, n)

## initial guess for tuning parameter
lam <- (sd ^ 2) * log(n) 
max_n_search_path_l0 <- 10
if (max_n_search_path_l0 != 500) {
  warning("Set max_n_search_path_l0 = 500 to replicate paper results!")
}
source("utils_sample_splitting.R")
## window sizes
window_sizes <- c(1, 30, 50)

# simulations with the same underyling mean, different noisy observations of this mean
d_out <- NULL
# agg. by true changepoints
d_true <- NULL
for (seed in 1:n_replicates) {
  for(delta in deltas) {
    print(paste0("At delta: ", delta))
    mean_i <- 0
    underlying_mean <- 0 * numeric(n)
    for (i in 1:n_changepts) {
      underlying_mean[(changepts_aug[i] + 1):changepts_aug[i + 1]] <- mean_i
      if (mean_i == 0) {
        mean_i <- delta
      } else {
        mean_i <- 0
      }
    }
    
    n <- length(underlying_mean)
    true_changepts <- which(underlying_mean[2:n] != underlying_mean[1:(n - 1)])
    k <- length(true_changepts)
    
    set.seed(seed)
    y <- underlying_mean + rnorm(n, sd = sd)
    
    ## attempt to properly calibrate tuning parameter for L0 segmentation to the correct number of changepoints
    fit <- changepoint_estimates(y, "L0", lam)
    if (length(fit$change_pts) != k) {
      paths <- estimate_paths(y, lambda_min = 1e-2, lambda_max = 1e1, max_iters = max_n_search_path_l0)
      distances <- abs(paths$path_stats$changepoints_n - k)
      distances[paths$path_stats$changepoints_n == 0] <- Inf
      lam_it <- paths$path_stats$lambda[which.min(distances)]
      lam <- lam_it
    }
    
    for (window_size in window_sizes) {
      ## l0 inference
      out <- changepoint_inference(y, "L0-fixed", lam, window_size, sd^2)
      if (length(out$change_pts) > 0) {
        df <- data.frame(change_pts = out$change_pts, pvals = out$pvals, approx = out$approximation_error, window_size = window_size, delta = delta, type = "L0-fixed", seed = seed)
        df_true_tmp <- true_changepts_df(df, true_changepts)
      } else { 
        df <- data.frame(change_pts = NA, pvals = NA, approx = NA, window_size = window_size, delta = delta, type = "L0-fixed", seed = seed)
        df_true_tmp <- true_changepts_df(df, true_changepts)
      }
      
      d_out <- rbind(d_out, df)
      d_true <- rbind(d_true, df_true_tmp)
      
      ## BS approaches
      out <- changepoint_inference(y, "BS-fixed", k, window_size = window_size, sig = sd ^ 2)
      if (length(out$change_pts) > 0) {
        df <- data.frame(change_pts = out$change_pts, pvals = out$pvals, approx = out$approximation_error, window_size = window_size, delta = delta,
                         type = "BS-fixed", seed = seed)
        df_true_tmp <- true_changepts_df(df, true_changepts)
      } else {
        df <- data.frame(change_pts = NA, pvals = NA, approx = NA, window_size = window_size, delta = delta, type = "BS-fixed", seed = seed)
        df_true_tmp <- true_changepts_df(df, true_changepts)
      }
      
      d_true <- rbind(d_true, df_true_tmp)
      d_out <- rbind(d_out, df)
    }
    
    
    for (type_i in c("BS-adaptive-M", "BS-adaptive-M-O-D")) {
      out <- changepoint_inference(y, type_i, k, sig = sd ^ 2)
      if (length(out$change_pts) > 0) {
        df <- data.frame(change_pts = out$change_pts, pvals = out$pvals, approx = out$approximation_error, window_size = NA, delta = delta, type = type_i, seed = seed)
        df_true_tmp <- true_changepts_df(df, true_changepts)
      } else {
        df <- data.frame(change_pts = NA, pvals = NA, approx = NA, window_size = NA, delta = delta, type = type_i, seed = seed)
        df_true_tmp <- true_changepts_df(df, true_changepts)
      }
      d_true <- rbind(d_true, df_true_tmp)
      d_out <- rbind(d_out, df)
    }
    
    
    ## sample splitting
    df_sample_split <- sample_split_l0_bs(y, k)
    df_sample_split$window_size <- NA
    df_sample_split$delta <- delta
    df_sample_split$seed <- seed
    d_out <- rbind(d_out, df_sample_split)
    
    if (!is.na(df_sample_split[df_sample_split$type == 'l0-sample-split', ]$change_pts)) {
      df_true_tmp1 <- true_changepts_df(df_sample_split[df_sample_split$type == 'l0-sample-split', ], true_changepts)  
    } else {
      df_true_tmp1 <- NULL
    }
    
    df_true_tmp2 <- true_changepts_df(df_sample_split[df_sample_split$type == 'bs-sample-split', ], true_changepts)
    df_true_tmp <- rbind(df_true_tmp1, df_true_tmp2)
    d_true <- rbind(d_true, df_true_tmp)
    
  }
}

# write.csv(d_out, file = "d_est.csv", row.names = FALSE)
# write.csv(d_true, file = "d_vs_true.csv", row.names = FALSE)





