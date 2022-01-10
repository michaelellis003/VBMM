#' Title
#'
#' @return
#' @export
#'
#' @examples
default_priors <- function() {

    return(
        list(
            mu_prior_mu = 0,
            sigmasq_prior_mu = 10^2,
            shape_prior_sigmasq = 0.01,
            scale_prior_sigmasq = 0.01,
            shape1_prior_V = 1,
            shape2_prior_V = 0.5,
            concentration_prior_phi = 2
        )
    )

}
