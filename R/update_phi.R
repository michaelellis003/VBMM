#' Title
#'
#' @param y
#' @param N
#' @param B
#' @param E_log_phi
#' @param E_log_p
#' @param E_sigma_sq
#' @param E_log_sigma_sq
#' @param E_SS
#'
#' @return
#' @export
#'
#' @examples
update_phi <- function(
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

    phi_tilde <- matrix(0, nrow = N, ncol = B)
    log_norm_dens <- (1/2)*(log(2*pi) + t(apply(E_SS, 1,
                                                function(x) E_log_sigma_sq + (1/E_sigma_sq)*x)))

    log_mix_dens <- array(0, dim = c(N, B, K))

    for(b in 1:B) {
        for(i in 1:N) {
            log_mix_dens[i, b, ] <- E_log_p[, b] - log_norm_dens[i, ]
        }
        phi_tilde[, b] <- E_log_phi[b] + rowSums(log_mix_dens[, b, ])
    }

    phi_q_phi <- multinomial_logit(phi_tilde)

    return(phi_q_phi)
}
