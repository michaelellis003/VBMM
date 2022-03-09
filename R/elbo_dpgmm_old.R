elbo <- function(
    mixture_model
) {

    # unpack model
    y <- mixture_model$y
    N <- length(y)
    B <- mixture_model$B
    K <- mixture_model$K
    max_iter <- mixture_model$max_iter
    prior_parameters <- mixture_model$prior_parameters
    variational_parameters <- mixture_model$variational_parameters

    # unpack prior parameters
    mu_mu_k <- prior_parameters$mu_mu_k
    sigmasq_mu_k <- prior_parameters$sigmasq_mu_k
    alpha_sigmasq_k <- prior_parameters$alpha_sigmasq_k
    beta_sigmasq_k <- prior_parameters$beta_sigmasq_k
    lambda <- prior_parameters$lambda
    a_lambda <- prior_parameters$a_lambda
    b_lambda <- prior_parameters$b_lambda

    # unpack variational parameters
    alpha_q_V <- variational_parameters$alpha_q_V
    beta_q_V <- variational_parameters$beta_q_V
    alpha_q_phi <- variational_parameters$alpha_q_phi
    mu_q_mu <- variational_parameters$mu_q_mu
    sigmasq_q_mu <- variational_parameters$sigmasq_q_mu
    A_q_sigmasq <- variational_parameters$A_q_sigmasq
    B_q_sigmasq <- variational_parameters$B_q_sigmasq
    p_q_p <- variational_parameters$p_q_p
    tilde_p_q_p <- variational_parameters$tilde_p_q_p

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
             (alpha_sigmasq_k - 1)*E_log_sigmasq + beta_sigmasq_k/E_sigmasq) +
        K*(lgamma((1+lambda)) - lgamma(lambda) + (lambda-1)*(digamma(beta_q_V) -
                                                                 digamma(alpha_q_V + beta_q_V))))


    # optimal density ----------------------------------------------------------
    optimal_density <- matrix(0, nrow = N, ncol = K)
    for(i in 1:N) {
        for(k in 1:K) {
            optimal_density[i, k] <- p_q_p[i, k]*(tilde_p_q_p[i, k] + sum(tilde_p_q_p[i, ]))
        }
    }
    optimal_density <- colSums(optimal_density)

    for(k in 1:K) {
        optimal_density <- sum(optimal_density -
            log(2*pi)/2 - log(sigmasq_mu_k)/2 - (1/(2*E_log_sigmasq))*(sigmasq_mu_k) +
            A_q_sigmasq*log(B_q_sigmasq) - lgamma(A_q_sigmasq) - (A_q_sigmasq+1)*E_log_sigmasq -
            B_q_sigmasq/E_sigmasq +
            lgamma(alpha_q_V + beta_q_V) - lgamma(alpha_q_V) - lgamma(beta_q_V) +
            (alpha_q_V + 1)*(alpha_q_V/(alpha_q_V + beta_q_V)) +
            (beta_q_V+1)*(beta_q_V/(alpha_q_V + beta_q_V)))

    }

    elbo <- log_joint_density - optimal_density

    # optimal variational density ----------------------------------------------


    # ## posterior - p
    # loglikelihood <- 0
    # z_cond <- 0
    # w_cond <- 0
    #
    # ## calculate expectations
    # E_sigma_sq <- B_q_sigmasq/A_q_sigmasq
    # E_log_sigma_sq <- log(B_q_sigmasq) - digamma(A_q_sigmasq)
    # E_log_pi <- stick_breaking_expectation(B, K, alpha_q_V, beta_q_V)
    # E_log_phi <- digamma(alpha_q_phi) + digamma(sum(alpha_q_phi))
    #
    # for(i in 1:N) {
    #     for(b in 1:B) {
    #         loglikelihood <- loglikelihood + phi_q_phi[i, b]*
    #             sum((-1/2)*log(2*pi) - (1/2)*E_log_sigma_sq -
    #             (1/2)*(1/E_sigma_sq)*((y[i] - mu_q_mu)^2 - sigmasq_q_mu)*
    #                 pi_q_pi[b, i,  ])
    #         z_cond <- z_cond + sum(pi_q_pi[b, i, ]*E_log_pi[, b])
    #         w_cond <- w_cond + phi_q_phi[i, b]*E_log_phi[b]
    #     }
    # }
    #
    # E_SS0 <- (mu_q_mu - mu0)^2  + sigmasq_q_mu # E[(mu_k - mu0)^2]
    # E_log_p0 <- stick_breaking_expectation(B, K, matrix(1, nrow = (K-1), ncol = B),
    #                                        matrix(beta, nrow = (K-1), ncol = B))
    #
    # mu_prior <- sum((-1/2)*log(2*pi) - (1/2)*E_log_sigma_sq - (1/(2*E_sigma_sq))*E_SS0)
    #
    # sigmasq_prior <- sum(A0*log(B0) - lgamma(A0) +
    #     (A0 + 1)*(1/E_log_sigma_sq) - B0/E_sigma_sq)
    #
    # phi_prior <- lgamma(B) - B*lgamma(1) +
    #     sum((1-1) * E_log_phi)
    #
    # p_prior <- sum(rowSums(E_log_p0))
    #
    # post <- loglikelihood + z_cond + w_cond + mu_prior + sigmasq_prior + phi_prior + p_prior
    #
    # ## variational - q
    # var <- 0
    # var_z <- 0
    # var_w <- 0
    # for(b in 1:B) {
    #     for(k in 1:K) {
    #         tmp <- as.numeric(t(pi_q_pi[b, , k]) %*% log(pi_q_pi[b, , k]))
    #         if(is.nan(tmp)) {
    #             var_z <- var_z
    #         } else {
    #             var_z <- var_z + tmp
    #         }
    #     }
    #     tmp1 <- as.numeric(t(phi_q_phi[, b]) %*% log(phi_q_phi[, b]))
    #     if(is.nan(tmp1)) {
    #         var_w <- var_w
    #     } else {
    #         var_w <- var_w + tmp1
    #     }
    # }
    #
    # lnd <- matrix(0, nrow = N, ncol = K)
    # for(k in 1:K){
    #     lnd[, k] <- (-1/2)*log(2*pi) - (1/2)*E_log_sigma_sq[k] -
    #         (1/2)*(1/E_sigma_sq[k])*((y - mu_q_mu[k])^2 - sigmasq_q_mu[k])
    # }
    # var_mu <- sum(rowSums(lnd))
    #
    # var_sigma_sq <- sum(A_q_sigmasq*log(B_q_sigmasq) - lgamma(A_q_sigmasq) -
    #                         (A_q_sigmasq + 1) * E_log_sigma_sq - B_q_sigmasq/E_sigma_sq)
    #
    # var_p <- sum(rowSums(E_log_pi))
    # var_phi <- lgamma(sum(1)) - sum(lgamma(1)) +
    #     sum((1-1) * E_log_phi)
    #
    # var <- var_z + var_w + var_mu + var_sigma_sq + var_p + var_phi
    #
    # elbo <- post - var

    return(elbo)
}
