#' Title
#'
#' @param mod
#' @param variational_parameters
#' @param expectations
#'
#' @return
#' @export
#'
#' @examples
update_expectations <- function(
    mod,
    variational_parameters = NULL,
    expectations = NULL
) {

    # unpack model
    y <- mod$y[,1]
    y_mat <- mod$y
    N <- mod$N
    B <- mod$B
    K <- mod$K

    # unpack variational parameters
    shape1_V <- variational_parameters$shape1_V
    shape2_V <- variational_parameters$shape2_V
    concentration_phi <- variational_parameters$concentration_phi
    mu_mu <- variational_parameters$mu_mu
    sigmasq_mu <- variational_parameters$sigmasq_mu
    shape_sigmasq <- variational_parameters$shape_sigmasq
    scale_sigmasq <- variational_parameters$scale_sigmasq
    phi_w <- variational_parameters$phi_w
    p_z <- variational_parameters$p_z

    # Calculate expectations

    # expected count of observations in each w and z cluster
    E_n <- matrix(NA, nrow = K, ncol = B)
    for(b in 1:B) {
        for(k in 1:K) {
            E_n[k, b] <- sum(p_z[, b, k] * phi_w[, b])
        }
    }

    # tmp <- matrix(NA, nrow = K, ncol = B)
    # for(b in 1:B) {
    #     tmp[, b] <- t(phi_w[, b] %*% as.matrix(p_z[, b, ]))
    # }

    E_w_n <- colSums(E_n) # expected count of observations in each w cluster
    E_z_n <- rowSums(E_n) # expected count of observations in each z cluster

    E_log_phi = digamma(concentration_phi) + digamma(sum(concentration_phi)) # E[log \phi_b]

    E_sigma_sq <- scale_sigmasq/shape_sigmasq # E[sigma^2_k]
    E_log_sigma_sq <- log(scale_sigmasq) - digamma(shape_sigmasq) # E[log sigma^2_k]

    E_SS <- t(apply(y_mat, 1, function(x) (x-mu_mu)^2 + sigmasq_mu)) # E[(y - mu_k)^2]
    E_log_p <- stick_breaking_expectation(B, K,
                                          alpha_q_V = shape1_V,
                                          beta_q_V = shape2_V) #E[log p_bk]

    expectations <- list(
        E_n = E_n,
        E_w_n = E_w_n,
        E_z_n = E_z_n,
        E_log_phi = E_log_phi,
        E_sigma_sq = E_sigma_sq,
        E_log_sigma_sq = E_log_sigma_sq,
        E_SS = E_SS,
        E_log_p = E_log_p
    )

    return(expectations)

}
