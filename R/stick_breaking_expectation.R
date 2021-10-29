#' Title
#'
#' @param B
#' @param K
#' @param alpha_q_V
#' @param beta_q_V
#'
#' @return
#' @export
#'
#' @examples
stick_breaking_expectation <- function(B, K, alpha_q_V, beta_q_V) {

    E_log_V <- digamma(alpha_q_V) - digamma(alpha_q_V + beta_q_V) # E[log V_bk]
    E_log_1_V <- digamma(beta_q_V) - digamma(alpha_q_V + beta_q_V) # E[log 1-V_bk]

    E_log_p <- matrix(NA, nrow = K, ncol = B)

    # k = 1
    E_log_p[1, ] <- E_log_V[1, ]

    # k = 2
    E_log_p[2, ] <- E_log_V[2, ] + E_log_1_V[1, ]

    for(k in 3:K) {
        ### stick-breaking representation
        if (k < K) {
            E_log_p[k, ] <- E_log_V[k, ] + colSums(E_log_1_V[1:(k-1), ])
        } else { # if k == K then E[log V_K] = E[log 1] = 0
            E_log_p[k, ] <- colSums(E_log_1_V[1:(k-1), ])
        }
    }

    return(E_log_p)

}
