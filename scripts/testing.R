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
w <- mix$w

B <- 2
K <- 20
iter <- 500
output <- VBMIM(
    y,
    w,
    B,
    K,
    max_iter = iter
)

x <- seq(-40, 60, length.out = 1000)
sigma_sq <- output$B_q_sigmasq/output$A_q_sigmasq
for(b in 1:B) {

    preds <- matrix(0, nrow = length(x), ncol = 20)
    for (k in 1:20) {
        preds[, k] <- colMeans(output$pi_q_pi[b, , ])[k]*dnorm(x, output$mu_q_mu[k], sqrt(sigma_sq[k]))
    }
    matplot(x, preds, type = "l", lty = 1)
    hist(y[w == b], freq = FALSE, add = TRUE, breaks = 500)
}

round(output$mu_q_mu, 2)
round(output$sigmasq_q_mu, 3)

round(output$B_q_sigmasq/output$A_q_sigmasq, 2)

round(output$alpha_q_phi/sum(output$alpha_q_phi), 2)
round(colMeans(output$pi_q_pi[1, , ]), 2)
round(colMeans(output$pi_q_pi[2, , ]), 2)

elbo_df <- data.frame(
    iteration = 1:iter,
    elbos = output$elbos
)

ggplot(elbo_df, aes(x = iteration, y = elbos)) +
    geom_point() +
    geom_line()




## mcmc ------------------------------------------------------------------------
n_mcmc <- 5000
burnin <- 2500
output_mcmc <- mcmc(
    y,
    w,
    B = 2,
    K = 20,
    n_mcmc,
    burnin,
    n_message = 500,
    prior_parameters = NULL
)
