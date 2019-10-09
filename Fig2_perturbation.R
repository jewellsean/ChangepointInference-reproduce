library(ChangepointInference)
library(tidyverse)
library(latex2exp)
library(gridExtra)

## plotting parameters
col_red <- "#d95f02"
col_blue <- "#1f78b4"

## generate sample data with one changepoint 
set.seed(1)
mu <- rep(c(1, 2), each = 100)
n <- length(mu)
y <- mu + rnorm(n, sd = 0.1)

## BS-adaptive-M inference
k <- 1
fit <- changepoint_inference(y, "BS-adaptive-M", 1, return_conditioning_sets = TRUE)

v <- construct_v_tL_tR(1, fit$change_pts, n, n)
vTy <- round(sum(v * y))

ymax <- 3
p0 <- data.frame(x = 1:n, y = y, mu = mu) %>% 
  ggplot() + 
  geom_point(aes(x, y), color = "gray", size = 0.5) +
  geom_line(aes(x, mu), color = "darkblue") + 
  ggtitle(TeX(paste0("a) Original data ($\\phi =", vTy, ")$"))) + 
  ylab(TeX("$y'(\\phi)$")) +
  coord_cartesian(ylim = c(0, ymax)) +
  theme_bw()

phi <- 0
yphi <- construct_perturbed_data_tL_tR(y, 1, fit$change_pts, n, phi)

p1 <- data.frame(x = 1:n, y = yphi, mu = mu) %>% 
  ggplot() + 
  geom_point(aes(x, y), color = col_red, size = 0.5) +
  ggtitle(TeX(paste0("b) Perturbed data ($\\phi = $", phi, ")"))) +
  geom_line(aes(x, mu), color = "darkblue") + 
  ylab('') +
  coord_cartesian(ylim = c(0, ymax)) +
  theme_bw()

phi <- -2
yphi <- construct_perturbed_data_tL_tR(y, 1, fit$change_pts, n, phi)
p2 <- data.frame(x = 1:n, y = yphi, mu = mu) %>% 
  ggplot() + 
  geom_point(aes(x, y), color = col_blue, size = 0.5) +
  ylab('') +
  ggtitle(TeX(paste0("c) Perturbed data ($\\phi = $", phi, ")"))) +
  geom_line(aes(x, mu), color = "darkblue") + 
  coord_cartesian(ylim = c(0, ymax)) +
  theme_bw()

q <- fit$conditioning_sets[[1]] %>% 
  mutate(y = 1, contained = factor(contained, levels = c(0, 1), labels = c("Not in model", "In model"))) %>% 
  ggplot() + 
  geom_rect(aes(xmin = min_mean, xmax = max_mean, ymin = -10, ymax = 10, fill = contained)) +
  coord_cartesian(xlim = c(-2, 2)) + 
  ggtitle(TeX(paste0("d) The set of values of $\\phi$ such that $M(y'(\\phi)) = M(y)$"))) + 
  xlab(TeX("$\\phi$")) + 
  ylab('') + 
  scale_fill_manual(values=c(col_red, col_blue)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position="none", legend.title = element_blank())

lay <- rbind(c(1, 2, 3), c(4, 4, 4))
grid.arrange(p0, p1, p2, q, layout_matrix = lay)