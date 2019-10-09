required_packages <- c("microbenchmark", "doParallel", "tidyverse", "magrittr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")

require(ChangepointInference)
require(microbenchmark)
library(doParallel)
library(tidyverse)
library(magrittr)

cl <- makeCluster(12)
registerDoParallel(cl)


set.seed(1)
sd <- 0.5
u <- rep(c(3, 5), each = 1000)
thj <- 1000
n <- length(u)
lam <- 2 * log(n)
window_sizes <- 10 ^ seq(1, 3, length.out = 10)

df <- foreach (seed = 1:50, .combine = rbind, .packages=c('ChangepointInference', 'microbenchmark')) %dopar% {
    y <- u + rnorm(n)
    df_t <- NULL
    for (window_size in window_sizes) {
        out <- microbenchmark("window" = {changepoint_inference(y, "L0-fixed", lam, window_size)}, times = 1L)
        df_tmp <- data.frame(window_size = window_size, time = out$time, seed = seed)
        df_t <- rbind(df_t, df_tmp)
    }
    df_t
}

stopCluster(cl)

df %>% mutate(time = time / 1e9) %>%
    group_by(window_size) %>%
    summarize(mean_time = mean(time), sd_time = sd(time)) %>%
    ggplot(aes(x = window_size, y = mean_time)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    geom_smooth(method = "lm", se = TRUE) +
    xlab("h") +
    ylab("Time (s)") +
    theme_bw()

out <- df %>% mutate(time = time / 1e9) %>%
    group_by(window_size) %>%
    summarize(mean_time = mean(time), sd_time = sd(time), n = n()) %$%
coef(lm(log10(mean_time) - 2 * log10(window_size) ~ 1))

p <- df %>% mutate(time = time / 1e9) %>%
    group_by(window_size) %>%
    summarize(mean_time = mean(time), sd_time = sd(time)) %>%
    ggplot(aes(x = window_size, y = mean_time)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    geom_abline(slope = 2, intercept = out[1], color = "red") +
    xlab("h") +
    ylab("Time (s)") +
    theme_bw()