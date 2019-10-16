## sample split p-values 
require(ChangepointInference)

sample_split_p_values <- function(est_change_pts, y_even, n ) {
  n_length_even <- length(y_even)
  n_estimated_changepts <- length(est_change_pts)
  pvals <- NA * numeric(n_estimated_changepts)
  aug_estimated_changepts <- c(0, est_change_pts, n_length_even)
  for (chg_i in 2:(n_estimated_changepts + 1)) {
    y_sub_1 <- y_even[(aug_estimated_changepts[chg_i - 1] + 1):aug_estimated_changepts[chg_i]]
    y_sub_2 <- y_even[(aug_estimated_changepts[chg_i] + 1):aug_estimated_changepts[chg_i + 1]]
    if ((length(y_sub_1) + length(y_sub_2)) > 2) {
      pval <- t.test(y_sub_1, y_sub_2, var.equal = TRUE)$p.value  
    } else { 
      pval <- NA
    }
    pvals[chg_i - 1] <- pval
  }
  return(data.frame(change_pts = seq(2, n, by = 2)[est_change_pts], pvals = pvals))
}


sample_split_l0_bs <- function(y, k) {
  y_odd <- y[seq(1, n, by = 2)]
  y_even <- y[seq(2, n, by = 2)]
  
  ### L0 paths 
  
  fit_odd_path <- estimate_paths(y_odd, max_iters = max_n_search_path_l0)
  distances <- abs(paths$path_stats$changepoints_n - k)
  distances[paths$path_stats$changepoints_n == 0] <- Inf
  lam_it <- paths$path_stats$lambda[which.min(distances)]
  fit_odd <- changepoint_estimates(y_odd, "L0", lam_it)
  if (length(fit_odd$change_pts) > 0) {
    l0_sample_split <- sample_split_p_values(fit_odd$change_pts, y_even, n)
    df1 <- data.frame(change_pts = l0_sample_split$change_pts, pvals = l0_sample_split$pvals, approx = 0, type = "l0-sample-split")  
  } else { 
    df1 <- data.frame(change_pts = NA, pvals = NA, approx = 0, type = "l0-sample-split")  
  }
  
  
  ## BS paths 
  fit_odd <- changepoint_estimates(y_odd, "BS", k)
  bs_sample_split <- sample_split_p_values(fit_odd$change_pts, y_even, n)
  df2 <- data.frame(change_pts = bs_sample_split$change_pts, pvals = bs_sample_split$pvals, approx = 0, type = "bs-sample-split")
  
  df <- rbind(df1, df2)
  return(df)
}