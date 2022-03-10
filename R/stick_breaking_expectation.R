#' Title
#'
#' @param B
#' @param K
#' @param E_log_V
#' @param E_log_1_V
#'
#' @return
#' @export
#'
#' @examples
stick_breaking_expectation <- function(B, K, E_log_V, E_log_1_V) {

    if (B > 1) {
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
    } else {
        E_log_p <- matrix(NA, nrow = K, ncol = B)

        # k = 1
        E_log_p[1] <- E_log_V[1]

        # k = 2
        E_log_p[2] <- E_log_V[2] + E_log_1_V[1]

        for(k in 3:K) {
            ### stick-breaking representation
            if (k < K) {
                E_log_p[k] <- E_log_V[k] + sum(E_log_1_V[1:(k-1)])
            } else { # if k == K then E[log V_K] = E[log 1] = 0
                E_log_p[k] <- sum(E_log_1_V[1:(k-1)])
            }
        }
    }

    return(E_log_p)

}
