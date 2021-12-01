library(VBMM)
library(tidyverse)

mix <- simulate_mixture(
    N = 500,
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

iter <- 100
output <- vbmm(
        y,
        B = 2,
        K = 20,
        max_iter = iter
)

output_old <- vbmm_old(
        y,
        B = 2,
        K = 20,
        mu0 = 0,
        sigmasq0 = 10^8,
        A0 = 0.01,
        B0 = 0.01,
        alpha_phi0 = 1,
        beta_V0 = 0.5,
        max_iter = iter)

## plot elbos
df <- data.frame(
    iteration = 1:iter,
    elbo = output$elbos
)

ggplot(df, aes(x = iteration, y = elbo)) +
    geom_point() +
    geom_line()

output$variational_parameters$mu_mu
output$variational_parameters$sigmasq_mu

round(colMeans(output$variational_parameters$phi_w), 2)
round(colMeans(output$variational_parameters$p_z[, 1, ]), 2)
round(colMeans(output$variational_parameters$p_z[, 2, ]), 2)

