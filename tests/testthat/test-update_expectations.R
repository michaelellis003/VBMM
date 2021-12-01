test_that("update_expectations function works", {

    mod <- list(
        y = matrix(0, nrow = 10, ncol = 2),
        N = 10,
        B = 2,
        K = 2
    )

    N <- mod$N
    B <- mod$B
    K <- mod$K

    variational_parameters <- list(
        shape1_V = matrix(1, nrow = K-1, ncol = B),
        shape2_V = matrix(1, nrow = K-1, ncol = B),
        concentration_phi = c(1, 1),
        mu_mu = c(0, 0),
        sigmasq_mu = c(1, 1),
        shape_sigmasq = c(1, 1),
        scale_sigmasq = c(1, 1),
        phi_w = matrix(1/2, nrow = N, ncol = B),
        p_z = array(1/2, dim = c(N, B, K))
    )

    expected_values <- update_expectations(mod = mod,
                        variational_parameters = variational_parameters,
                        expectations = NULL
                        )

    expect_equal(
        expected_values,
        list(
            E_n = matrix(N/(K*B), nrow = K, ncol = B),
            E_w_n = rep(N/B, B),
            E_z_n = rep(N/K, K),
            E_log_phi = digamma(variational_parameters$concentration_phi) +
                digamma(sum(variational_parameters$concentration_phi)),
            E_sigma_sq = variational_parameters$scale_sigmasq/variational_parameters$shape_sigmasq,
            E_log_sigma_sq = log(variational_parameters$scale_sigmasq) -
                digamma(variational_parameters$shape_sigmasq),
            E_SS = matrix(1, nrow = N, ncol = K),
            E_log_p = matrix(-1, nrow = K, ncol = B)
        )
    )


})
