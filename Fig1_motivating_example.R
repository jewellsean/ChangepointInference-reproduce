## motivating example
library(ChangepointInference)
library(tidyverse)
library(gridExtra)

set.seed(1)
n <- 200
n_changepoints <- 20

true_changepts <- sort(unique(floor(runif(n_changepoints, 1, n))))
n_changepts <- length(true_changepts)
changepts_aug <- c(0, true_changepts, n)
underlying_mean <- 0 * numeric(n)
for (i in 1:n_changepts) {
  underlying_mean[(changepts_aug[i] + 1):changepts_aug[i + 1]] <- rnorm(1, mean = 0, sd = 4)
}
true_changepts <- which(underlying_mean[2:n] != underlying_mean[1:(n - 1)])
k <- length(true_changepts)

set.seed(1)
y <- underlying_mean + rnorm(n)

out1 <- changepoint_inference(y, type = "BS-adaptive-M-O-D", tuning_parameter = k, sig = 1)
out2 <- changepoint_inference(y, type = "BS-adaptive-M", tuning_parameter = k, sig = 1)

df <- rbind(data.frame(change_pts = out1$change_pts, pvals = out1$pvals, type = "M-O-D"), 
            data.frame(change_pts = out2$change_pts, pvals = out2$pvals, type = "M"))

df_underlying <- NULL
last_chg_pt <- 0
for (i in 1:k) {
  df_underlying <- rbind(df_underlying, 
                         data.frame(x0 = last_chg_pt + 1, xend = true_changepts[i], y0 = underlying_mean[true_changepts[i]], yend = underlying_mean[true_changepts[i]]))
  last_chg_pt <- true_changepts[i]
} 
df_underlying <- rbind(df_underlying, 
                       data.frame(x0 = last_chg_pt + 1, xend = n, y0 = underlying_mean[n], yend = underlying_mean[n]))

plot_data_changepts_rejection <- function(out, title) {
  # plot estimated changepoints
  chg_cols <- data.frame(x = out$change_pts, pvals = out$pvals) %>%
    mutate(reject = if_else(pvals < 0.05, "rejection", "fail_to_reject"))
  
  df <- data.frame(x = 1:length(y), y = y) 
  
  q <- df %>% 
    ggplot() + 
    # geom_step(data = data.frame(x = 1:length(y), mu = underlying_mean), aes(x, mu), color = "darkblue", lwd = 1) + 
    geom_segment(data = df_underlying, aes(x=x0, y=y0, xend=xend, yend=yend), color = "darkblue", lwd = 1) +
    geom_vline(data = chg_cols, aes(xintercept = x, color = reject), lwd = 0.5) + 
    geom_point(aes(x, y), alpha = 0.5) +
    ylab("") + 
    xlab("") + 
    theme_bw() + 
    ggtitle(title) +
    theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  return(q)
}

q1 <- plot_data_changepts_rejection(out1, "a)")
q2 <- plot_data_changepts_rejection(out2, "b)")

grid.arrange(q1, q2)
