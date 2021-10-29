update_p <- function(
    y,
    N,
    B,
    K,
    E_log_phi,
    E_log_p,
    E_sigma_sq,
    E_log_sigma_sq,
    E_SS
) {

    p_tilde <- array(NA, dim = c(N, B, K))
    log_norm_dens <- (1/2)*(log(2*pi) + t(apply(E_SS, 1,
                                                function(x) E_log_sigma_sq + (1/E_sigma_sq)*x)))

    log_mix_dens <- array(0, dim = c(N, B, K))
    p_q_p <- list()

    for(b in 1:B) {
        for(i in 1:N) {
            p_tilde[i, b, ] <- E_log_p[, b] - log_norm_dens[i, ] + sum(E_log_phi)
        }
        p_q_p[[b]] <- multinomial_logit(p_tilde[, b, ])
    }

    return(p_q_p)
}
