test_that("stick_breaking_expectation function works", {

    B <- 2
    K <- 2
    alpha_q_V <- matrix()
    beta_q_V <- rep

    alpha_q_V <- matrix(alpha_q_V, nrow = K, ncol = 1)
    beta_q_V <-  matrix(beta_q_V, nrow = K, ncol = 1)

    E_pi <- matrix(E_pi, nrow = K, ncol = 1)
    tmp <- stick_breaking_expectation(B=1, K, alpha_q_V, beta_q_V)
})
