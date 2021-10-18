#' Title
#'
#' @param tilde_pi_q_pi
#'
#' @return
#' @export
#'
#' @examples
calculate_pi <- function(tilde_pi_q_pi) {

    # get maximum values in each row
    A <- apply(tilde_pi_q_pi, 1, function(x) max(x))

    tmp <- apply(tilde_pi_q_pi, 1, function(x) log(sum(exp(x-max(x)))))

    pi_q_pi <- 1/exp(apply(tilde_pi_q_pi, 2, function(x) A + tmp - x))

    return(pi_q_pi)
}
