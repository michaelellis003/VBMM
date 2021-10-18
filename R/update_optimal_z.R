#' Title
#'
#' @param y
#' @param N
#' @param K
#' @param E_sigma_sq
#' @param E_log_sigma_sq
#' @param E_log_V
#' @param E_log_1_V
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
    E_sigma_sq,
    E_log_sigma_sq,
    E_log_V,
    E_log_1_V,
    E_SS
) {

    tilde_pi_q_pi <- matrix(0, nrow = N, ncol = K)
    log_norm_dens <- (1/2)*(log(2*pi) + t(apply(E_SS, 1,
                                       function(x) E_log_sigma_sq + (1/E_sigma_sq)*x)))

    # k = 1
    tilde_pi_q_pi[, 1] <- E_log_V[1] - log_norm_dens[, 1]

    for(k in 2:K) {
        ### stick-breaking representation
        if (k < K) {
            tilde_pi_q_pi[, k] <- E_log_V[k] + sum(E_log_1_V[1:(k-1)]) - log_norm_dens[, k]
        } else { # if k == K then E[log V_K] = E[log 1] = 0
            tilde_pi_q_pi[, k] <- sum(E_log_1_V[1:(k-1)]) - log_norm_dens[, k]
        }
    }

    pi_q_pi <- calculate_pi(tilde_pi_q_pi)

    return(pi_q_pi)
}
