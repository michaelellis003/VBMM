#' Title
#'
#' @param log_prob_tilde
#'
#' @return
#' @export
#'
#' @examples
multinomial_logit <- function(log_prob_tilde) {

    N <- nrow(log_prob_tilde)
    K <- ncol(log_prob_tilde)

    phi <- matrix(0, nrow = N, ncol = K)
    for(i in 1:N) {
        a <- max(log_prob_tilde[i, ])
        phi[i, ] <- 1/exp(a + log(sum(exp(log_prob_tilde[i, ] - a))) - log_prob_tilde[i, ])
        # for(k in 1:K) {
        #     a <- max(log_prob_tilde[i, ])
        #     phi[i, ] <- 1/exp(a + log(sum(exp(log_prob_tilde[i, ] - a))) - log_prob_tilde[i, ])
        #     if(all(phi[i, ] == rep(1,K))) {
        #         phi[i, ] <- rep(1/K, K)
        #     }
        # }
    }

    # get maximum values in each row
    # A <- apply(phi_tilde, 1, function(x) max(x))
    #
    # tmp <- apply(phi_tilde, 1, function(x) log(sum(exp(x-max(x)))))
    #
    # phi_q_phi <- 1/exp(apply(phi_tilde, 2, function(x) A + tmp - x))

    return(phi)
}
