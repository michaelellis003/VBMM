library(VBMM)
library(tidyverse)
library(MCMCpack)

## B = 1 -----------------------------------------------------------------------
N <- 1000
K <- 3
mu <- c(-20, 5, 25)
sigmasq <- c(0.5, 0.5, 0.5)
probs <- c(1/2, 1/4, 1/4)

sample_clusters <- sample(1:K, size = N, replace = TRUE)
y <- rnorm(N, mean = mu[sample_clusters], sd = sqrt(sigmasq[sample_clusters]))

B <- 1
K <- 20
max_iter <- 50
output <- vbmm(
    y = y,
    B = B,
    K = K,
    max_iter,
    w = NULL
)

plot_q_mu(output, type = "boxplot")
round(output$mu_q_mu, 2)
plot_q_sigmasq(output, type = "boxplot")
plot_elbo(output)

# B = 2 ------------------------------------------------------------------------
## simulate data
mix <- simulate_mixture(
    N = 1000,
    B = 3,
    K = 3,
    mu = c(-1/2, 0, 1/2),
    sigmasq = c(0.5, 0.5, 0.5),
    phi = c(2/3, 1/6, 1/6),
    p = list(c(1/3, 1/3, 1/3), c(1/2, 1/4, 1/4), c(1/8, 1/4, 5/8))
)

y <- mix$y
w <- mix$w
B <- mix$B
K <- 20
max_iter <- 50
output <- vbmm(
    y = y,
    w = w,
    B = B,
    K = K,
    max_iter = max_iter
)

library(ggtern)
alpha_q_phi <- output$alpha_q_phi
phi_sim <- MCMCpack::rdirichlet(1000, alpha_q_phi)
dat <- data.frame(
    x = phi_sim[, 1],
    y = phi_sim[, 2],
    z = phi_sim[, 3]
)

dat_true <- data.frame(
    x = mix$phi[1],
    y = mix$phi[2],
    z = mix$phi[3]
)

ggtern(data = dat, aes(x=x, y=y, z=z)) +
    geom_point(alpha = 0.01) +
    stat_density_tern(geom = "polygon",
                      n = 200,
                      aes(fill = ..level..)) +
    geom_point(data = dat_true, color = "orange", size = 2, shape = 17)

plot_q_mu(output, type = "boxplot")
round(output$mu_q_mu, 2)

plot_q_sigmasq(output, type = "boxplot")
round(output$B_q_sigmasq/output$A_q_sigmasq, 2)
plot_elbo(output)
