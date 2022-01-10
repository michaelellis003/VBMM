library(VBMM)
library(tidyverse)

mix <- simulate_mixture(
    N = 10000,
    B = 2,
    K = 3,
    mu = c(-20, 10, 35),
    sigmasq = c(0.5, 0.5, 0.5),
    phi = c(1/3, 2/3),
    p = list(c(1/3, 1/3, 1/3), c(1/6, 1/6, 2/3))
)

df <- data.frame(y = mix$y, w = mix$w, z=mix$z)

ggplot(df) +
    geom_histogram(aes(x = y, fill = factor(z))) +
    facet_wrap(~w)

y <- mix$y

iter <- 2000
output <- VBMIM(
    y,
    B = 2,
    K = 20,
    max_iter = iter
)

round(output$mu_q_mu, 2)
round(output$sigmasq_q_mu, 3)

sigma_sq <- output$B_q_sigmasq/output$A_q_sigmasq
round(output$B_q_sigmasq/output$A_q_sigmasq, 2)

round(colMeans(output$phi_q_phi), 2)
round(colMeans(output$pi_q_pi[1, , ]), 2)
round(colMeans(output$pi_q_pi[2, , ]), 2)

x <- seq(-40, 60, length.out = 1000)
dat <- matrix(0, nrow = length(x), ncol = 20)
for (k in 1:K) {
    dat[, k] <- colMeans(output$pi_q_pi[1, , ])[k]*dnorm(x, output$mu_q_mu[k], sqrt(sigma_sq[k]))
}
matplot(x, dat, type = "l", lty = 1)

x <- seq(-40, 60, length.out = 1000)
dat_true <- matrix(0, nrow = length(x), ncol = 3)
for (k in 1:3) {
    dat_true[, k] <- c(1/3, 1/3, 1/3)[k]*dnorm(x, c(-20, 10, 35)[k], sqrt(c(0.5, 0.5, 0.5)[k]))
}
matplot(x, dat_true, type = "l", add = T, lty = 2)
