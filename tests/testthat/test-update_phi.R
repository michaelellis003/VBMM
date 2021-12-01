test_that("update_phi function works", {

    E_log_phi <- log(c(1/4, 3/4))
    E_log_p <- log(c(1/6, 1/6, 2/3))

    E_sigma_sq <- c(0.5, 0.5, 0.5)
    E_log_sigma_sq <- log(c(0.5, 0.5, 0.5))

    mu_mu <- c(-10, 0, 10)
    sigmasq_mu <- c(0, 0, 0)

    y <- c(-1, 0)

    dens <- matrix(0, nrow = 2, ncol = 3)
    for(k in 1:3) {
        dens[, k] <- E_log_p[k] - (1/2)*log(2*pi) - (1/2)*E_log_sigma_sq[k] -
            (1/(2*E_sigma_sq[k]))*((y - mu_mu[k])^2 + sigmasq_mu[k])
    }
    log_phi_tilde <- E_log_phi + rowSums(dens)
    log_phi_tilde <- matrix(log_phi_tilde, nrow = 2, ncol = 2, byrow = TRUE)
    multinomial_logit(log_phi_tilde)

})
