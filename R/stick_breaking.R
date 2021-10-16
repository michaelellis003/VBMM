#' Stick Breaking Representation for Dirichlet Process
#'
#' @param V
#' @param K
#'
#' @return
#' @export
#'
#' @examples
stick_breaking <- function(V, K) {

    pi <- rep(NA, K)

    for(k in 1:K) {
        if(k == 1){
            pi[k] <- V[k]
        } else {
            stick_prob <- V[k]
            for(j in 1:(k-1)){
                stick_prob <- stick_prob*(1 - V[j])
            }
            pi[k] <- stick_prob
        }
    }

    return(pi)
}
