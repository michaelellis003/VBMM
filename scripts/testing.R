
mix <- simulate_mixture(
    N = 500,
    B = 2,
    K = 3,
    mu = c(-20, 5, 25),
    sigmasq = c(0.5, 0.5, 0.5),
    phi = c(1/2, 1/2),
    p = list(c(1/3, 1/3, 1/3), c(1/3, 1/3, 1/3))
)

y <- mix$y

output <- vbmm(
        y,
        B = 2,
        K = 20,
        mu0 = 0,
        sigmasq0 = 10^8,
        A0 = 0.01,
        B0 = 0.01,
        alpha_phi0 = 1,
        beta_V0 = 0.5,
        max_iter = 5000
)

output$mu_q_mu
output$A_q_sigmasq/output$B_q_sigmasq
output$alpha_q_phi/sum(output$alpha_q_phi)
