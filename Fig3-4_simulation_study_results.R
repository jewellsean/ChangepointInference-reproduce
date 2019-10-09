## All experiment figures  
library(tidyverse)
library(magrittr)
library(latex2exp)
library(gridExtra)

d_vs_true <- read_csv("d_vs_true.csv")
d_est <- read_csv("d_est.csv")

set.seed(1112)
n <- 2000
n_changepts <- 50
delta <- 3

changepts <- unique(sort(sample(2:n, size = n_changepts, replace = FALSE)))
n_changepts <- length(changepts)

changepts_aug <- c(0, changepts, n)

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

DETECTION_THRESHOLD <- 2
DETECTION_THRESHOLD_SAMPLE_SPLIT <- DETECTION_THRESHOLD
REJECTION_THRESHOLD <- 0.05

## Figre 3:

p_error <- d_est %>% filter(delta == 1e-20) %>% 
  mutate(window_size = if_else(is.na(window_size), -1, window_size)) %>% 
  filter(type %in% c("L0-fixed", "BS-fixed", "BS-adaptive-M", "BS-adaptive-M-O-D")) %>% 
  mutate(type = factor(type, levels = c("BS-adaptive-M-O-D", "BS-adaptive-M", "BS-fixed", "L0-fixed"), labels = c("Approach 1", "Approach 2", "Approach 3", "Approach 4"))) %>% 
  filter(window_size %in% c("-1", "50"), !is.na(type)) %>% 
  ggplot(aes(sample = pvals)) +
  stat_qq(distribution = qunif, size = 0.01, alpha = 0.2) +
  stat_qq_line(distribution = qunif) +
  facet_grid(cols = vars(type)) + 
  xlab("Theoretical Unif[0, 1] quantiles") + 
  ylab(TeX("Observed $p$-value quantiles")) +
  ggtitle("b)") +
  theme_bw() +
  theme(legend.position = "none")


set.seed(1)
y <- underlying_mean + rnorm(n)
q <- data.frame(x = 1:length(y), y = y) %>% 
  ggplot() + 
  geom_line(data = data.frame(x = 1:length(y), mu = underlying_mean), aes(x, mu), color = "darkblue", lwd = 1) + 
  geom_vline(data = data.frame(x = true_changepts), aes(xintercept = x),  color = "grey", lwd = 0.5) + 
  geom_point(aes(x, y), alpha = 0.5) +
  ylab("") + 
  xlab("") + 
  ggtitle("a)") +
  theme_bw() + 
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



p <- d_vs_true %>% filter(window_size %in% c(50, NA)) %>% 
  mutate(correct_detection = if_else(dist_nearest_est <= DETECTION_THRESHOLD, 1, 0), 
         rejected = if_else(pvals < REJECTION_THRESHOLD, 1, 0)) %>% 
  group_by(delta, window_size, type) %>% 
  summarize(power = mean(correct_detection * rejected)) %>% 
  filter(type %in% c("BS-fixed", "BS-adaptive-M-O-D", "BS-adaptive-M")) %>% 
  mutate(type2 = factor(type, levels = c("BS-adaptive-M-O-D", "BS-adaptive-M", "BS-fixed"), labels = 
                          c("Conditioning on the estimated changepoints, order, and signs with BS (Approach 1)", 
                              "Conditioning on the estimated changepoints with BS (Approach 2)", 
                              "Conditioning on the jth estimated changepoint with BS (Approach 3)"
                                                                                         ))) %>% 
  ggplot(aes(delta, power, color = type2)) + 
  geom_line() +
  geom_point() + 
  xlab(TeX("$\\delta$")) + 
  ylab("Power") + 
  ggtitle("c)") +
  theme_bw() + 
  theme(legend.position="bottom", legend.title = element_text(size=0)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

grid.arrange(q,p_error,  p)

## max improvements 
max_improvements <- d_vs_true %>% filter(window_size %in% c(50, NA), delta > 0.1) %>%
  mutate(correct_detection = if_else(dist_nearest_est <= DETECTION_THRESHOLD, 1, 0), 
         rejected = if_else(pvals < REJECTION_THRESHOLD, 1, 0)) %>% 
  group_by(delta, type) %>% 
  summarize(power = mean(correct_detection * rejected)) %>% 
  filter(type %in% c("BS-fixed", "BS-adaptive-M-O-D", "BS-adaptive-M"), !is.na(power)) %>% 
  spread(type, power) %>% 
  mutate(diff_0_3 = `BS-fixed` - `BS-adaptive-M`, 
         diff_3_1 = `BS-adaptive-M` - `BS-adaptive-M-O-D`) 

max(max_improvements$diff_0_3, na.rm = TRUE)
max(max_improvements$diff_3_1, na.rm = TRUE)


## Figure 4

p1 <- d_vs_true %>% 
    filter(window_size %in% c(50, NA), delta > 0.1) %>% 
    mutate(pvals = if_else(is.na(pvals), 1, pvals),
           correct_detection = if_else(type %in% c("l0-sample-split", "bs-sample-split"), if_else(dist_nearest_est <= DETECTION_THRESHOLD_SAMPLE_SPLIT, 1, 0),
                                       if_else(dist_nearest_est <= DETECTION_THRESHOLD, 1, 0)), 
         rejected = if_else(pvals < REJECTION_THRESHOLD, 1, 0)) %>% 
    group_by(delta, window_size, type) %>% 
    summarize(power = mean(correct_detection * rejected), detection_prob = mean(correct_detection)) %>% 
    filter(type %in% c("L0-fixed", "BS-fixed", "BS-adaptive-M-O-D", "BS-adaptive-M", "l0-sample-split", "bs-sample-split")) %>% 
    mutate(type2 = factor(type, levels = c("BS-adaptive-M-O-D", "BS-adaptive-M", "BS-fixed", "L0-fixed", "bs-sample-split", "l0-sample-split"), 
                        labels = c("Conditioning on the estimated changepoints, order, and signs with BS (Approach 1)", 
                                   "Conditioning on the estimated changepoints with BS (Approach 2)", 
                                   "Conditioning on the jth estimated changepoint with BS (Approach 3)", 
                                   "Conditioning on the jth estimated changepoint with L0 segmentation (Approach 4)", 
                                   "Sample splitting with binary segmentation (Approach 5)",
                                   "Sample splitting with L0 segmentation (Approach 6)" 
                                   ))) %>% 
  ggplot(aes(delta, power, color = type2)) + 
  geom_line() + 
  xlab(TeX("$\\delta$")) + 
  ylab("Power") + 
  labs(color="") +  
  theme_bw() +
  ggtitle("a)") +
  theme(legend.position="bottom", legend.title = element_text(size=0)) +
  guides(color=guide_legend(nrow=6,byrow=TRUE))


p2 <- d_vs_true %>% 
  mutate(pvals = if_else(is.na(pvals), 1, pvals),
         correct_detection = if_else(type %in% c("l0-sample-split", "bs-sample-split"), if_else(dist_nearest_est <= DETECTION_THRESHOLD_SAMPLE_SPLIT, 1, 0),
         if_else(dist_nearest_est <= DETECTION_THRESHOLD, 1, 0)), 
         rejected = if_else(pvals < REJECTION_THRESHOLD, 1, 0)) %>% 
  filter(type %in% c("L0-fixed", "BS-fixed", "l0-sample-split", "bs-sample-split"), delta > 0.1) %>% 
  group_by(delta, window_size, type) %>% 
  summarize(detection_prob = mean(correct_detection)) %>% 
  mutate(type2 = factor(type, levels = c("BS-fixed", "L0-fixed", "bs-sample-split", "l0-sample-split"), 
                        labels = c("Binary segmentation (Approaches 1-3)", 
                                   "L0 segmentation (Approach 4)", 
                                   "Sample splitting with binary segmentation (Approach 5)",
                                   "Sample splitting with L0 segmentation (Approach 6)"
                                   ))) %>% 
  ggplot(aes(delta, detection_prob, linetype = type2)) + 
  geom_line() + 
  xlab(TeX("$\\delta$")) + 
  ylab("Detection probability") + 
  labs(color="") +  
  theme_bw() +
  ggtitle("b)") +
  theme(legend.position="bottom", legend.title = element_text(size=0)) +
  guides(linetype=guide_legend(nrow=4,byrow=TRUE))


p <- d_vs_true %>% 
  mutate(pvals = if_else(is.na(pvals), 1, pvals),
         correct_detection = if_else(type %in% c("l0-sample-split", "bs-sample-split"), if_else(dist_nearest_est <= DETECTION_THRESHOLD_SAMPLE_SPLIT, 1, 0),
                                     if_else(dist_nearest_est <= DETECTION_THRESHOLD, 1, 0)), 
         rejected = if_else(pvals < REJECTION_THRESHOLD, 1, 0)) %>%
  filter(window_size %in% c(1, 30, 50)) %>% 
  group_by(delta, window_size, type) %>% 
  summarize(power = mean(correct_detection * rejected)) %>% 
  ungroup() %>% 
  filter(type %in% c("L0-fixed", "BS-fixed")) %>% 
  mutate(type2 = factor(type, levels = c("BS-fixed", "L0-fixed"), 
                        labels = c("Conditioning on the jth estimated changepoint with BS (Approach 3)", 
                                   "Conditioning on the jth estimated changepoint with L0 segmentation (Approach 4)")), 
         window_size = factor(window_size, levels = c(1, 30, 50), labels =  c("Window size: 1", "Window size: 30", "Window size: 50"))) %>% 
  ggplot(aes(delta, power, color = type2)) + 
  geom_line() + 
  xlab(TeX("$\\delta$")) + 
  ylab("Power") + 
  labs(color="") +  
  theme_bw() +
  facet_wrap(~window_size) +
  ggtitle("c)") +
  theme(legend.position="bottom", legend.title = element_text(size=0)) +
  guides(color=guide_legend(nrow=1,byrow=TRUE))

lay <- rbind(c(1, 2), c(4, 4))
grid.arrange(p1, p2, p, layout_matrix = lay)