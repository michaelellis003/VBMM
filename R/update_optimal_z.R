#' Title
#'
#' @param y
#' @param N
#' @param K
#' @param E_log_pi
#' @param E_sigma_sq
#' @param E_log_sigma_sq
#' @param E_SS
#'
#' @return
#' @export
#'
#' @examples
update_optimal_z <- function(
    y,
    N,
    K,
    E_log_pi,
    E_sigma_sq,
    E_log_sigma_sq,
    E_SS
) {

    tilde_pi_q_pi <- matrix(0, nrow = N, ncol = K)
    log_norm_dens <- (1/2)*(log(2*pi) + t(apply(E_SS, 1,
                                       function(x) E_log_sigma_sq + (1/E_sigma_sq)*x)))

    for(k in 1:K) {
        tilde_pi_q_pi[, k] <- E_log_pi[k] - log_norm_dens[, k]
    }

    pi_q_pi <- normalize_pi(tilde_pi_q_pi)

    return(pi_q_pi)
}
