#' Title
#'
#' @param mixture_model
#'
#' @return
#' @export
#'
#' @examples
elbo_dpgmm <- function(
    mixture_model
) {

    # unpack model
    y <- mixture_model$y
    N <- length(y)
    B <- mixture_model$B
    K <- mixture_model$K
    prior_parameters <- mixture_model$prior_parameters
    variational_parameters <- mixture_model$variational_parameters

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
    tilde_p_q_p <- variational_parameters$tilde_p_q_p
    alpha_q_lambda <- variational_parameters$alpha_q_lambda
    beta_q_lambda <- variational_parameters$beta_q_lambda

    # calculate expectations needed
    E_sigmasq <- B_q_sigmasq/A_q_sigmasq
    E_log_sigmasq <- log(B_q_sigmasq) - digamma(A_q_sigmasq)
    E_log_p <- stick_breaking_expectation(B, K, alpha_q_V, beta_q_V)

    # log joint density --------------------------------------------------------
    log_joint_density <- matrix(0, nrow = N, ncol = K)
    for(i in 1:N) {
        for(k in 1:K) {
            log_joint_density[i, k] <- (-1*p_q_p[i, k]/2)*((log(2*pi) + E_log_sigmasq[k]) +
                                                               E_sigmasq[k]*((y[i] - mu_q_mu[k])^2 + sigmasq_q_mu[k])) +
                p_q_p[i, k]*E_log_p[k]
        }
    }
    log_joint_density <- colSums(log_joint_density)

    log_joint_density <- sum(log_joint_density - log(2*pi)/2 - log(sigmasq_mu_k)/2 -
                                 (1/(2*sigmasq_mu_k))*((mu_q_mu - mu_mu_k)^2 + sigmasq_q_mu) +
                                 (log(beta_sigmasq_k) - lgamma(alpha_sigmasq_k) -
                                      (alpha_sigmasq_k + 1)*E_log_sigmasq + beta_sigmasq_k/E_sigmasq) +
                                 K*(2*(digamma(alpha_q_lambda) - log(beta_q_lambda)) +
                                        (alpha_q_lambda/(alpha_q_lambda+beta_q_lambda)-1)*(digamma(beta_q_V) - digamma(alpha_q_V + beta_q_V))))

    log_joint_density <- log_joint_density +
        a_lambda*log(b_lambda) - lgamma(a_lambda) +
        (a_lambda-1)*(digamma(alpha_q_lambda) - log(beta_q_lambda)) -
        beta_q_lambda*(alpha_q_lambda/(alpha_q_lambda+beta_q_lambda))


    # optimal density ----------------------------------------------------------
    optimal_density <- matrix(0, nrow = N, ncol = K)
    for(i in 1:N) {
        for(k in 1:K) {
            optimal_density[i, k] <- p_q_p[i, k]*(tilde_p_q_p[i, k] + sum(tilde_p_q_p[i, ]))
        }
    }
    optimal_density <- colSums(optimal_density)

    optimal_density <- sum(optimal_density -
                               log(2*pi)/2 - log(sigmasq_mu_k)/2 - (1/(2*E_log_sigmasq))*(sigmasq_mu_k) +
                               A_q_sigmasq*log(B_q_sigmasq) - lgamma(A_q_sigmasq) - (A_q_sigmasq+1)*E_log_sigmasq -
                               B_q_sigmasq/E_sigmasq +
                               lgamma(alpha_q_V + beta_q_V) - lgamma(alpha_q_V) - lgamma(beta_q_V) +
                               (alpha_q_V - 1)*(digamma(alpha_q_V) - digamma(alpha_q_V + beta_q_V)) +
                               (beta_q_V  - 1)*(digamma(beta_q_V) - digamma(alpha_q_V + beta_q_V)))

    optimal_density <- optimal_density + alpha_q_lambda*log(beta_q_lambda) -
        lgamma(alpha_q_lambda) + (alpha_q_lambda-1)*(digamma(alpha_q_lambda) - digamma(alpha_q_lambda + beta_q_lambda)) -
        beta_q_lambda*(alpha_q_lambda/(alpha_q_lambda + beta_q_lambda))

    elbo <- log_joint_density - optimal_density

    return(elbo)
}
