library(VBMM)
library(tidyverse)

N <- 500
K <- 3
mu <- c(-20, 5, 25)
sigmasq <- c(0.5, 0.5, 0.5)
probs <- c(1/2, 1/4, 1/4)

sample_clusters <- sample(1:K, size = N, replace = TRUE)
y <- rnorm(N, mean = mu[sample_clusters], sd = sqrt(sigmasq[sample_clusters]))

B <- 1
K <- 20
max_iter <- 2000
output <- vbmm(
    y = y,
    B = B,
    K = K,
    max_iter,
    w = NULL
)

plot_q_mu(output, type = "boxplot")
plot_q_sigmasq(output, type = "boxplot")
plot_elbo(output)


round(output$mu_q_mu, 2)
