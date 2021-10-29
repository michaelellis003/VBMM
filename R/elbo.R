elbo <- function(
    y,
    K,
    mu0,
    sigmasq0,
    A,
    B,
    beta0,
    mu_q_mu,
    sigmasq_q_mu,
    E_sigma_sq,
    E_log_sigma_sq,
    pi_q_pi,
    E_n,
    E_log_pi,
    E_log_V,
    E_log_1_V,
    E_SS
) {
    N <- length(y)
    likeli <- rep(NA, N)
    for(i in 1:N) {
        likeli[i] <- sum(pi_q_pi[i, ]*(-1/2)*(log(2*pi) + E_log_sigma_sq + (1/E_log_sigma_sq)*E_SS[i, ])) +
            sum(pi_q_pi[i, ]*E_log_pi)
    }
    likeli <- sum(likeli)

    V_beta_dens <- K*(lgamma(1+beta0) - lgamma(beta0)) + (beta0 - 1)*sum(E_log_1_V)

    mu_norm_dens <- (-K/2)*(log(2*pi) + log(sigmasq0)) -
        (1/(2*sigmasq0))*sum((mu_q_mu - mu0)^2 + sigmasq_q_mu)

    sigmasq_ig_dens <- K*(A*log(B) - lgamma(A)) - (A+1)*sum(E_log_sigma_sq) -
        B*sum(1/E_sigma_sq)


}
