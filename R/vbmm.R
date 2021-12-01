#' Title
#'
#' @param y
#' @param B
#' @param K
#' @param max_iter
#' @param prior_parameters
#' @param initial_values
#'
#' @return
#' @export
#'
#' @examples
vbmm <- function(
    y,
    B = 2,
    K = 20,
    max_iter = 30,
    prior_parameters = NULL,
    initial_values = NULL
) {

    N <- length(y) # sample size
    y_mat <- matrix(rep(y, K), nrow = N, ncol = K)
    mod <- list(
        y = y_mat,
        N = N,
        B = B,
        K = K
    )
    elbos <- rep(-Inf, max_iter)

    # get prior parameters
    if(is.null(prior_parameters)) {
        prior_parameters <- default_priors()
    }

    # initialize variational parameters
    variational_parameters <-
        update_variational_parameters(mod,
                                      variational_parameters = NULL,
                                      expectations = NULL,
                                      prior_parameters = prior_parameters,
                                      initialize = TRUE,
                                      initial_values = initial_values)

    # calculate initial expectations
    expected_values <- update_expectations(mod,
                                           variational_parameters = variational_parameters)

    for(i in 1:max_iter) {
        if(i %% 100 == 0) {
            message("Iteration ", i, " of ", max_iter)
        }

        ## update variational parameters
        variational_parameters <-
            update_variational_parameters(mod,
                                          variational_parameters = variational_parameters,
                                          expectations = expected_values,
                                          prior_parameters = prior_parameters,
                                          initialize = FALSE,
                                          initial_values = NULL)

        print(round(colMeans(variational_parameters$phi_w), 2))
        print(round(colMeans(variational_parameters$p_z[, 1, ]), 2))
        print(round(colMeans(variational_parameters$p_z[, 2, ]), 2))

        if(all(round(colMeans(variational_parameters$phi_w), 2) == c(1,1))) {

            return(
                list(
                    mod = mod,
                    variational_parameters = variational_parameters,
                    expected_values = expected_values,
                    elbos = elbos
                )
            )

        }

        ## update expected values
        expected_values <- update_expectations(mod,
                                               variational_parameters = variational_parameters)

        ## update elbo
        elbos[i] <- elbo(mod,
                         variational_parameters,
                         expectations = expected_values,
                         prior_parameters)

    }

    return(
        list(
            mod = mod,
            variational_parameters = variational_parameters,
            expected_values = expected_values,
            elbos = elbos
        )
    )
}
