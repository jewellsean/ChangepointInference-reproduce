library(ChangepointInference)
library(changepoint) ## for dataset 
library(tidyverse)
library(gridExtra)
data(HC1, package = "changepoint")

n <- 2000
y <- HC1[1:n] / sd(HC1[1:n])
## plot scaled data
plot(y, cex = 0.1)

## L0 inference 
fit <- changepoint_inference(y, "L0-fixed", log(length(y)), window_size = 50)


REJECTION_THRESHOLD <- 0.05

K <- length(fit$change_pts)

fit_bs_thj <- changepoint_inference(y, "BS-fixed", K, window_size = 50)
fit_bs_order_signs <- changepoint_inference(y, "BS-adaptive-M-O-D", K)
fit_bs_cps <- changepoint_inference(y, "BS-adaptive-M", K)

d1 <- data.frame(cps = fit_bs_thj$change_pts, pvals = fit_bs_thj$pvals, type = "bs_thj")
d2 <- data.frame(cps = fit_bs_order_signs$change_pts, pvals = fit_bs_order_signs$pvals, type = "bs_order_signs")
d3 <- data.frame(cps = fit_bs_cps$change_pts, pvals = fit_bs_cps$pvals, type = "bs_cps")
d4 <- data.frame(cps = fit$change_pts, pvals = fit$pvals, type = "l0_thj")

dd <- rbind(d1, d2, d3, d4)

dres <- dd %>% mutate(rejected = as.factor(if_else(pvals < REJECTION_THRESHOLD , 1, 0)), 
                      type = factor(type, levels = c("bs_order_signs", "bs_cps", "bs_thj", "l0_thj"), 
                                    labels = c("Conditioning on the estimated changepoints, order, and signs with BS (Approach 1)", 
                                               "Conditioning on the estimated changepoints with BS (Approach 2)", 
                                               "Conditioning on the jth estimated changepoint with BS (Approach 3)", 
                                               "Conditioning on the jth estimated changepoint with L0 segmentation (Approach 4)"
                                    )))

dat <- data.frame(position = 1:length(y), y = y)

p <- dat %>% ggplot(aes(x = position, y = y)) + 
  geom_point(size = 0.1, alpha = 0.8) +
  xlab("Position") + 
  ylab("G-C content") + 
  geom_vline(data = dres, aes(xintercept = cps, color = rejected)) + 
  facet_wrap(~type) + 
  labs(color="") +  
  theme_bw() +
  theme(legend.position="none", 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
grid.arrange(p)
