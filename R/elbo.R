#' Title
#'
#' @param mod
#' @param variational_parameters
#' @param expectations
#' @param prior_parameters
#'
#' @return
#' @export
#'
#' @examples
elbo <- function(
    mod,
    variational_parameters,
    expectations,
    prior_parameters
) {

    # unpack model
    y <- mod$y[,1]
    N <- mod$N
    B <- mod$B
    K <- mod$K

    # unpack expected values
    E_n <- expectations$E_n
    E_w_n <- expectations$E_w_n
    E_z_n <- expectations$E_z_n
    E_log_phi <- expectations$E_log_phi
    E_sigma_sq <- expectations$E_sigma_sq
    E_log_sigma_sq <- expectations$E_log_sigma_sq
    E_SS <- expectations$E_SS
    E_log_p <- expectations$E_log_p

    # unpack variational parameters
    shape1_V <- variational_parameters$shape1_V
    shape2_v <- variational_parameters$shape2_V
    concentration_phi <- variational_parameters$concentration_phi
    mu_mu <- variational_parameters$mu_mu
    sigmasq_mu <- variational_parameters$sigmasq_mu
    shape_sigmasq <- variational_parameters$shape_sigmasq
    scale_sigmasq <- variational_parameters$scale_sigmasq
    phi_w <- variational_parameters$phi_w
    p_z <- variational_parameters$p_z

    # unpack prior parameters
    concentration_prior_phi <- prior_parameters$concentration_prior_phi
    shape1_prior_V <- prior_parameters$shape1_prior_V
    shape2_prior_V <- prior_parameters$shape2_prior_V
    sigmasq_prior_mu <- prior_parameters$sigmasq_prior_mu
    mu_prior_mu <- prior_parameters$mu_prior_mu
    shape_prior_sigmasq <- prior_parameters$shape_prior_sigmasq
    scale_prior_sigmasq <- prior_parameters$scale_prior_sigmasq

    elbo <- 0
    post <- 0

    ## posterior - p
    loglikelihood <- 0
    z_cond <- 0
    w_cond <- 0
    lnd <- log_normal_dens(E_SS, E_log_sigma_sq, E_sigma_sq)

    for(i in 1:N) {
        for(b in 1:B) {
            loglikelihood <- loglikelihood + phi_w[i, b]*sum((-1)*lnd[i, ] * p_z[i, b, ])
            z_cond <- z_cond + sum(p_z[i, b, ]*E_log_p[, b])
            w_cond <- w_cond + phi_w[i, b]*E_log_phi[b]
        }
    }

    E_SS0 <- (mu_mu - mu_prior_mu)^2  + sigmasq_mu # E[(mu_k - mu0)^2]
    E_log_p0 <- stick_breaking_expectation(B, K, matrix(shape1_prior_V, nrow = (K-1), ncol = B),
                                           matrix(shape2_prior_V, nrow = (K-1), ncol = B))

    mu_prior <- sum((-1/2)*log(2*pi) - (1/2)*E_log_sigma_sq - (1/(2*E_sigma_sq))*E_SS0)

    sigmasq_prior <- sum(shape_prior_sigmasq*log(scale_prior_sigmasq) - lgamma(shape_prior_sigmasq) +
        (shape_prior_sigmasq + 1)*(1/E_log_sigma_sq) - scale_prior_sigmasq/E_sigma_sq)

    phi_prior <- lgamma(B*concentration_prior_phi) - B*lgamma(concentration_prior_phi) +
        sum((concentration_prior_phi-1) * E_log_phi)

    p_prior <- sum(rowSums(E_log_p0))

    post <- loglikelihood + z_cond + w_cond + mu_prior + sigmasq_prior + phi_prior + p_prior

    ## variational - q
    var <- 0
    var_z <- 0
    var_w <- 0
    for(b in 1:B) {
        for(k in 1:K) {
            tmp <- as.numeric(t(p_z[, b, k]) %*% log(p_z[, b, k]))
            if(is.nan(tmp)) {
                var_z <- var_z
            } else {
                var_z <- var_z + tmp
            }
        }
        tmp1 <- as.numeric(t(phi_w[, b]) %*% log(phi_w[, b]))
        if(is.nan(tmp1)) {
            var_w <- var_w
        } else {
            var_w <- var_w + tmp1
        }
    }

    var_mu <- sum(rowSums((-1)*lnd))

    var_sigma_sq <- sum(shape_sigmasq*log(scale_sigmasq) - lgamma(shape_sigmasq) -
                            (shape_sigmasq + 1) * E_log_sigma_sq - scale_sigmasq/E_sigma_sq)

    var_p <- sum(rowSums(E_log_p))
    var_phi <- lgamma(sum(concentration_phi)) - sum(lgamma(concentration_prior_phi)) +
        sum((concentration_phi-1) * E_log_phi)

    var <- var_z + var_w + var_mu + var_sigma_sq + var_p + var_phi

    elbo <- post - var

    return(elbo)
}
