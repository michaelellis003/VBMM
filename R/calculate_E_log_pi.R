#' Title
#'
#' @param K
#' @param E_log_V
#' @param E_log_1_V
#'
#' @return
#' @export
#'
#' @examples
calculate_E_log_pi <- function(K, E_log_V, E_log_1_V) {

    E_log_pi <- rep(NA, K)

    # k = 1
    E_log_pi[1] <- E_log_V[1]

    for(k in 2:K) {
        ### stick-breaking representation
        if (k < K) {
            E_log_pi[k] <- E_log_V[k] + sum(E_log_1_V[1:(k-1)])
        } else { # if k == K then E[log V_K] = E[log 1] = 0
            E_log_pi[k] <- sum(E_log_1_V[1:(k-1)])
        }
    }

    return(E_log_pi)

}
