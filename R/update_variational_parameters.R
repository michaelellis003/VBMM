#' Title
#'
#' @param mod
#' @param variational_parameters
#' @param expectations
#' @param prior_parameters
#' @param initialize
#' @param initial_values
#'
#' @return
#' @export
#'
#' @examples
update_variational_parameters <- function(
    mod,
    variational_parameters = NULL,
    expectations = NULL,
    prior_parameters = NULL,
    initialize = FALSE,
    initial_values = NULL
) {

    # unpack model
    y <- mod$y[,1]
    N <- mod$N
    B <- mod$B
    K <- mod$K

    if(initialize) {
        variational_parameters <- list()

        if(is.null(initial_values)) {

            variational_parameters$shape1_V <- matrix(1, nrow = (K-1), ncol = B)
            variational_parameters$shape2_V <- matrix(runif((K-1)*B, 0, 1), nrow = (K-1), ncol = B)

            variational_parameters$concentration_phi <- runif(B, 0, 1)

            variational_parameters$mu_mu <- sort(rnorm(K, mean = mean(y), sd = sd(y)))
            variational_parameters$sigmasq_mu <- 1/rgamma(K, 0.1, 0.1)

            variational_parameters$shape_sigmasq <- runif(K, 0, 1)
            variational_parameters$scale_sigmasq <- runif(K, 0, 1)

            variational_parameters$phi_w <- matrix(0, nrow = N, ncol = B)
            variational_parameters$phi_w[, 1] <- runif(N, 0.5, 1)
            variational_parameters$phi_w[, 2] <- 1 - variational_parameters$phi_w[, 1]

            variational_parameters$p_z <- array(0, dim = c(N, B, K))
            variational_parameters$p_z[, 1, 1] <- runif(N, 0.5, 1)
            variational_parameters$p_z[, 1, 2] <- 1 - variational_parameters$p_z[, 1, 1]

            variational_parameters$p_z[, 2, 1] <- runif(N, 0.5, 1)
            variational_parameters$p_z[, 2, 2] <- 1 - variational_parameters$p_z[, 2, 1]

            # V <- matrix(NA, nrow = K-1, ncol = B)
            # for(b in 1:B) {
            #     for(k in 1:(K-1)) {
            #         V[k, b] <- rbeta(1, variational_parameters$shape1_V[k, b],
            #                          variational_parameters$shape2_V[k, b])
            #     }
            # }
            #
            # for (b in 1:B) {
            #     variational_parameters$p_z[, b, ] <- matrix(stick_breaking(V[, b], K), N, K, byrow = TRUE)
            # }

        }
    } else {

        E_n <- expectations$E_n
        E_w_n <- expectations$E_w_n
        E_z_n <- expectations$E_z_n
        E_log_phi <- expectations$E_log_phi
        E_sigma_sq <- expectations$E_sigma_sq
        E_log_sigma_sq <- expectations$E_log_sigma_sq
        E_SS <- expectations$E_SS
        E_log_p <- expectations$E_log_p

        ## optimal density for w_ib
        variational_parameters$phi_w <-
            update_phi(mod = mod,
                       expectations = expectations,
                       variational_parameters = variational_parameters)


        ## optimal density for z_ibk
        variational_parameters$p_z <-
            update_p(mod = mod,
                     expectations = expectations,
                     variational_parameters = variational_parameters)

        SS_wt_tmp <- array(NA, dim = c(N, B, K))
        y_wt_tmp <- array(NA, dim = c(N, B, K))
        for(b in 1:B) {
            for(k in 1:K) {
                y_wt_tmp[ , b, k] <- y * variational_parameters$phi_w[ , b] *
                    variational_parameters$p_z[, b, k]

                SS_wt_tmp[, b, k] <- ((y - variational_parameters$mu_mu[k])^2 +
                    variational_parameters$sigmasq_mu[k]) *
                    variational_parameters$phi_w[ , b] *
                    variational_parameters$p_z[, b, k]
            }
        }

        SS_wt <- SS_wt_tmp[, 1, ]
        y_wt <- y_wt_tmp[, 1, ]
        for(b in 2:B) {
            y_wt <- y_wt + y_wt_tmp[, b, ]
            SS_wt <- SS_wt + SS_wt_tmp[, b, ]
        }

        ## optimal density for \phi
        concentration_prior_phi <- prior_parameters$concentration_prior_phi
        for(b in 1:B) {
            variational_parameters$concentration_phi[b] <-
                concentration_prior_phi +
                sum(variational_parameters$phi_w[, b])
        }

        for(k in 1:K){

            for(b in 1:B) {

                ## optimal density for V_bk
                if(k < K) {
                    variational_parameters$shape1_V[k, b] <-
                        prior_parameters$shape1_prior_V + E_n[k, b]

                    variational_parameters$shape2_V[k, b] <-
                        prior_parameters$shape2_prior_V + sum(E_n[(k+1):K, b])
                }
            }

            ## optimal density for mu_k
            sigmasq_prior_mu <- prior_parameters$sigmasq_prior_mu
            mu_prior_mu <- prior_parameters$mu_prior_mu

            variational_parameters$sigmasq_mu[k] <-
                1/(1/sigmasq_prior_mu + (1/E_sigma_sq[k])*E_z_n[k])

            variational_parameters$mu_mu[k] <- variational_parameters$sigmasq_mu[k] *
                (mu_prior_mu/sigmasq_prior_mu + (1/E_sigma_sq[k])*sum(y_wt[, k]))


            ## optimal density for sigma^2_k
            shape_prior_sigmasq <- prior_parameters$shape_prior_sigmasq
            scale_prior_sigmasq <- prior_parameters$scale_prior_sigmasq

            variational_parameters$shape_sigmasq[k] <- shape_prior_sigmasq + E_z_n[k]/2
            variational_parameters$scale_sigmasq[k] <- scale_prior_sigmasq +
                (1/2)*(sum(SS_wt[, k]))
        }

    }

    return(variational_parameters)
}
