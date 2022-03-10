#' Title
#'
#' @param mixture_model
#'
#' @return
#' @export
#'
#' @examples
elbo <- function(
    mixture_model
) {

    # unpack model
    y <- mixture_model$y
    N <- length(y)
    w <- mixture_model$w
    B <- mixture_model$B
    K <- mixture_model$K
    prior_parameters <- mixture_model$prior_parameters
    variational_parameters <- mixture_model$variational_parameters
    expected_values <- mixture_model$expected_values

    # unpack prior parameters
    mu_mu_k <- prior_parameters$mu_mu_k
    sigmasq_mu_k <- prior_parameters$sigmasq_mu_k
    alpha_sigmasq_k <- prior_parameters$alpha_sigmasq_k
    beta_sigmasq_k <- prior_parameters$beta_sigmasq_k
    a_lambda <- prior_parameters$a_lambda
    b_lambda <- prior_parameters$b_lambda

    # unpack variational parameters
    alpha_q_V <- variational_parameters$alpha_q_V
    beta_q_V <- variational_parameters$beta_q_V
    mu_q_mu <- variational_parameters$mu_q_mu
    sigmasq_q_mu <- variational_parameters$sigmasq_q_mu
    A_q_sigmasq <- variational_parameters$A_q_sigmasq
    B_q_sigmasq <- variational_parameters$B_q_sigmasq
    p_q_p <- variational_parameters$p_q_p
    log_tilde_p_q_p <- variational_parameters$log_tilde_p_q_p
    alpha_q_lambda <- variational_parameters$alpha_q_lambda
    beta_q_lambda <- variational_parameters$beta_q_lambda

    # unpack expected values
    E_sigma_sq <- expected_values$E_sigma_sq
    E_log_sigma_sq <- expected_values$E_log_sigma_sq
    E_log_V <- expected_values$E_log_V
    E_log_1_V <- expected_values$E_log_1_V
    E_log_p <- expected_values$E_log_p
    E_lambda <- expected_values$E_lambda

    # calculate other expectations needed
    E_log_lambda <- digamma(alpha_q_lambda) - digamma(beta_q_lambda)

    if (B == 1) {
        Nw <- N
        w_mat <- matrix(1, nrow = N, ncol = B)
    } else {
        w_mat <- matrix(0, nrow = N, ncol = B)
        w_mat[cbind(1:N, w)] <- 1
        Nw <- colSums(w_mat)
    }

    # log joint density --------------------------------------------------------
    log_normal <- matrix(0, N, K)
    log_joint_density <- matrix(0, N, B)
    for (b in 1:B) {
        for (k in 1:K) {
            E_RSS <- (y-mu_q_mu[k])^2 + sigmasq_q_mu[k]
            log_likelihood <- -log(2*pi)/2 - E_log_sigma_sq[k]/2 -
                E_RSS/(2*E_sigma_sq[k])

            log_normal[, K] <- p_q_p[b, , k]*(E_log_p[k] + log_likelihood)
        }
        log_joint_density[, b] <- w_mat[, b]*rowSums(log_normal)
    }

    log_joint_density <- sum(colSums(log_joint_density)) +
        sum(-log(2*pi)/2 - E_log_sigma_sq/2 - ((mu_q_mu - mu_mu_k)^2)/(2*E_sigma_sq)) +
        sum(alpha_sigmasq_k*log(beta_sigmasq_k) - lgamma(alpha_sigmasq_k) -
                (alpha_sigmasq_k + 1)*E_log_sigma_sq - beta_sigmasq_k/E_sigma_sq) +
        sum(a_lambda*log(b_lambda) - lgamma(a_lambda) + (a_lambda-1)*E_log_lambda - b_lambda*E_lambda)

        if(B == 1) {
            log_joint_density <- log_joint_density +
                sum(K*(2*E_log_lambda + (E_lambda-1)*E_log_1_V))
        } else {
            log_joint_density <- log_joint_density +
                sum(rowSums(apply(E_log_1_V, 1, function(x) K*(2*E_log_lambda + (E_lambda-1)*x))))
        }

    # # optimal density ----------------------------------------------------------
    optimal_density <- matrix(0, nrow = N, ncol = B)
    log_cat_desity <- matrix(0, nrow = N, ncol = K)
    for (b in 1:B) {
        log_sum_tilde_p <- log(rowSums(exp(log_tilde_p_q_p[b, , ])))
        for(k in 1:K) {
            log_cat_desity[, k] <- p_q_p[b, , k]*(log_tilde_p_q_p[b, , k] + log_sum_tilde_p)
        }
        optimal_density[, b] <- w_mat[, b]*rowSums(log_cat_desity)
    }

    optimal_density <- colSums(optimal_density)

    optimal_density <- sum(optimal_density) +
        sum(-log(2*pi)/2 - log(sigmasq_mu_k)/2 - sigmasq_mu_k/(2*E_log_sigma_sq)) +
        sum(A_q_sigmasq*log(B_q_sigmasq) - lgamma(A_q_sigmasq) - (A_q_sigmasq+1)*E_log_sigma_sq -
                B_q_sigmasq/E_sigma_sq) +
        sum(alpha_q_lambda*log(beta_q_lambda) -
                lgamma(alpha_q_lambda) + (alpha_q_lambda-1)*E_log_lambda - beta_q_lambda*E_lambda) +
        sum(colSums(lgamma(alpha_q_V + beta_q_V) - lgamma(alpha_q_V) - lgamma(beta_q_V) +
                        (alpha_q_V - 1)*E_log_V + (beta_q_V - 1)*E_log_1_V))

    elbo <- log_joint_density - optimal_density

    return(elbo)
}
