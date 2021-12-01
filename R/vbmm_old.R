#' Title
#'
#' @param y
#' @param B
#' @param K
#' @param mu0
#' @param sigmasq0
#' @param A0
#' @param B0
#' @param alpha_phi0
#' @param beta_V0
#' @param max_iter
#'
#' @return
#' @export
#'
#' @examples
vbmm_old <- function(
    y,
    B = 2,
    K = 20,
    mu0 = 0,
    sigmasq0 = 10^8,
    A0 = 0.01,
    B0 = 0.01,
    alpha_phi0 = 1,
    beta_V0 = 0.5,
    max_iter = 30
) {

    N <- length(y) # sample size
    elbos <- rep(-Inf, max_iter)
    y_mat <- matrix(rep(y, K), nrow = N, ncol = K)

    ## storing parameters
    alpha_q_phi <- rep(alpha_phi0, B)

    alpha_q_V <- matrix(1, nrow = (K-1), ncol = B)
    beta_q_V <- matrix(runif(K-1, 0, 1), nrow = (K-1), ncol = B)

    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- 1/rgamma(K, 1, 1)

    A_q_sigmasq <- runif(K, 0, 1)
    B_q_sigmasq <- runif(K, 0, 1)

    phi_q_phi <- matrix(0, nrow = N, ncol = B)


    ## Calculate initial expectations
    E_n <- matrix(NA, nrow = K, ncol = B)
    E_log_phi = digamma(alpha_q_phi) + digamma(sum(alpha_q_phi)) # E[log \phi_b]
    E_sigma_sq <- B_q_sigmasq/A_q_sigmasq # E[sigma^2_k]
    E_log_sigma_sq <- log(B_q_sigmasq) - digamma(A_q_sigmasq) # E[log sigma^2_k]
    E_SS <- t(apply(y_mat, 1, function(x) (x-mu_q_mu)^2 + sigmasq_q_mu)) # E[(y - mu_k)^2]
    E_log_p <- stick_breaking_expectation(B, K, alpha_q_V, beta_q_V) #E[log p_bk]

    for(m in 1:max_iter) {
        message("Iteration ", m, " of ", max_iter)

        ## optimal density for w_ib
        phi_q_phi <- update_phi_old(y, N, B, K, E_log_phi, E_log_p, E_sigma_sq,
                              E_log_sigma_sq, E_SS)


        ## optimal density for z_ibk
        p_q_p <- update_p_old(y, N, B, K, E_log_phi, E_log_p, E_sigma_sq, E_log_sigma_sq, E_SS)

        # expected count of observations in each w and z cluster
        for(b in 1:B) {
            E_n[, b] <- t(phi_q_phi[, b] %*% p_q_p[[b]])
        }

        E_w_n <- colSums(E_n) # expected count of observations in each w cluster
        E_z_n <- rowSums(E_n) # expected count of observations in each z cluster

        ## optimal density for \phi
        for(b in 1:B) {
            alpha_q_phi[b] <- alpha_phi0 + sum(phi_q_phi[, b])
        }

        for(k in 1:K){

            ## optimal density for V_bk
            for(b in 1:B) {
                if(k < K) {
                    alpha_q_V[k, b] <- 1 + E_n[k, b]
                    beta_q_V[k, b] <- beta_V0 + sum(E_n[(k+1):K, b])
                }
            }

            ## optimal density for mu_k
            wts <- phi_q_phi[, 1] * p_q_p[[1]][, k]
            for(b in 1:B) {
                wts <- wts + phi_q_phi[, b] * p_q_p[[b]][, k]
            }

            sigmasq_q_mu[k] <- 1/(1/sigmasq0 + (1/E_sigma_sq[k])*E_z_n[k])
            mu_q_mu[k] <- sigmasq_q_mu[k]*(mu0/sigmasq0 + (1/E_sigma_sq[k])*sum(y*wts))

            ## optimal density for sigma^2_k
            A_q_sigmasq[k] <- A0 + E_z_n[k]/2
            B_q_sigmasq[k] <- B0 + (1/2)*(sum(((y - mu_q_mu[k])^2 + sigmasq_q_mu[k])*wts))
        }

        ## update expectations
        E_log_phi = digamma(alpha_q_phi) + digamma(sum(alpha_q_phi)) # E[log \phi_b]
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq # E[sigma^2_k]
        E_log_sigma_sq <- log(B_q_sigmasq) - digamma(A_q_sigmasq) # E[log sigma^2_k]
        E_SS <- t(apply(y_mat, 1, function(x) (x-mu_q_mu)^2 + sigmasq_q_mu)) # E[(y - mu_k)^2]
        E_log_p <- stick_breaking_expectation(B, K, alpha_q_V, beta_q_V) #E[log p_bk]

        ## elbo
        # elbos[m] <- elbo(y, y_mat, N, B, K, mu0, sigmasq0, A0, B0, alpha_phi0, beta_V0,
        #                 phi_q_phi, p_q_p, mu_q_mu, sigmasq_q_mu, A_q_sigmasq, B_q_sigmasq,
        #                 E_n, E_log_phi, E_sigma_sq, E_log_sigma_sq, E_SS, E_log_p)
    }

    return(
        list(
            phi_q_phi = phi_q_phi,
            p_q_p = p_q_p,
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            alpha_q_phi = alpha_q_phi,
            elbos = elbos
        )
    )

}
