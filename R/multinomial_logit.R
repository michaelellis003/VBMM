#' Title
#'
#' @param phi_tilde
#'
#' @return
#' @export
#'
#' @examples
multinomial_logit <- function(phi_tilde) {

    # get maximum values in each row
    A <- apply(phi_tilde, 1, function(x) max(x))

    tmp <- apply(phi_tilde, 1, function(x) log(sum(exp(x-max(x)))))

    phi_q_phi <- 1/exp(apply(phi_tilde, 2, function(x) A + tmp - x))

    return(phi_q_phi)
}
