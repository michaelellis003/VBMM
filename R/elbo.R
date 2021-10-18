elbo <- function(
    y,
    K,
    mu0,
    sigmasq0,
    A,
    B,
    beta,
    E_n,
    pi_q_pi,
    alpha_q_V,
    beta_q_V,
    mu_q_mu,
    sigmasq_q_mu,
    A_q_sigmasq,
    B_q_sigmasq
) {

    elbo <- digamma(alpha_q_V) - digamma(alpha_q_V + beta_q_V) + sum(E_n[(k+1):K])


}
