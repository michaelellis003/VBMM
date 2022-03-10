#' Title
#'
#' @param log_prob_tilde
#'
#' @return
#' @export
#'
#' @examples
multinomial_logit <- function(log_prob_tilde) {

    # get maximum values in each row
    A <- apply(log_prob_tilde, 1, function(x) max(x))
    tmp <- apply(log_prob_tilde, 1, function(x) log(sum(exp(x-max(x)))))
    probs <- 1/exp(apply(log_prob_tilde, 2, function(x) A + tmp - x))

    return(probs)
}
