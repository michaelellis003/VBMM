#' Title
#'
#' @param mod
#' @param expectations
#' @param variational_parameters
#'
#' @return
#' @export
#'
#' @examples
update_phi <- function(
    mod,
    expectations,
    variational_parameters
) {

    # unpack model
    y <- mod$y[,1]
    N <- mod$N
    B <- mod$B
    K <- mod$K

    # unpack expectations
    E_log_phi <- expectations$E_log_phi
    E_log_p <- expectations$E_log_p
    E_sigma_sq <- expectations$E_sigma_sq
    E_log_sigma_sq <- expectations$E_log_sigma_sq
    E_SS <- expectations$E_SS

    log_mix_dens <- array(0, dim = c(N, B, K))
    log_phi_tilde <- matrix(0, nrow = N, ncol = B)

    lnd <- log_normal_dens(E_SS,
                           E_log_sigma_sq,
                           E_sigma_sq)

    for(b in 1:B) {
        for(k in 1:K) {
            log_mix_dens[, b, k] <- E_log_p[k, b] - lnd[, k]
        }
        log_phi_tilde[, b] <- E_log_phi[b] + rowSums(log_mix_dens[, b, ])
    }

    phi_w <- multinomial_logit(log_phi_tilde)

    return(phi_w)
}
