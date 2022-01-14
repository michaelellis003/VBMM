#' Title
#'
#' @param y
#' @param B
#' @param K
#' @param mu0
#' @param sigmasq0
#' @param A0
#' @param B0
#' @param beta
#' @param max_iter
#'
#' @return
#' @export
#'
#' @examples
VBMIM_unknown_w <- function(
    y,
    B = 2,
    K = 20,
    mu0 = 0,
    sigmasq0 = 10^2,
    A0 = 0.01,
    B0 = 0.01,
    beta = 0.5,
    max_iter = 30
) {

    N <- length(y) # sample size
    sum_y <- sum(y)
    elbos <- rep(-Inf, max_iter)

    ## initialize parameters
    tilde_pi_q_pi <- array(0, dim = c(B, N, K))
    pi_q_pi <- array(0, dim = c(B, N, K))

    alpha_q_V <- matrix(1, nrow = K-1, ncol = B)
    beta_q_V <- matrix(rbeta((K-1)*B, 1, 1), nrow = K-1, ncol = B)

    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- 1/rgamma(K, 1, 1)

    A_q_sigmasq <- rep(1, K)
    B_q_sigmasq <- rep(1, K)

    dens_phi_q_phi <-  array(0, dim = c(B, N, K))
    tilde_phi_q_phi <- matrix(0, nrow = N, ncol = B)
    phi_q_phi <- matrix(0, nrow = N, ncol = B)
    phi_q_phi[, 1] <- runif(N, min = 0.5, max = 1)
    phi_q_phi[, 2] <- 1 - phi_q_phi[, 1]

    alpha_q_phi <- rep(0.5, B)

    for(m in 1:max_iter) {

        if(m %% 100 == 0) {
            message("Iteration ", m, " of ", max_iter)
        }

        ## expectation of sigma^2
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq

        ## optimal density for z_ik
        for(b in 1:B) {
            for(k in 1:K) {
                ### stick-breaking representation
                E_pi_k <- 0
                if(k == 1) {
                    E_pi_k <- digamma(alpha_q_V[k, b]) - digamma(alpha_q_V[k, b] + beta_q_V[k, b])
                } else if (k > 1 & k < K) {
                    E_pi_k <- digamma(alpha_q_V[k, b]) - digamma(alpha_q_V[k, b] + beta_q_V[k, b]) +
                        sum(digamma(beta_q_V[1:(k-1), b]) -
                                digamma(alpha_q_V[1:(k-1), b] + beta_q_V[1:(k-1), b]))
                } else { # if k == K then E[log V_K] = E[log 1] = 0
                    E_pi_k <- sum(digamma(beta_q_V[1:(k-1), b]) -
                                      digamma(alpha_q_V[1:(k-1), b] + beta_q_V[1:(k-1), b]))
                }

                tilde_pi_q_pi[b, , k] <- E_pi_k -
                    (1/2)*(digamma(A_q_sigmasq[k]) + log(B_q_sigmasq[k]) +
                               (1/E_sigma_sq)*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k])) +
                    digamma(alpha_q_phi[b]) + digamma(sum(alpha_q_phi))

            }

            pi_q_pi[b, , ] <- multinomial_logit(tilde_pi_q_pi[b, , ])
                #t(apply(tilde_pi_q_pi[b, , ], 1, function(x)exp(x)/sum(exp(x))))
        }

        ## optimal density for w_ik
        for (b in 1:B) {
            for (k in 1:K) {
                dens_phi_q_phi[b, , k] <- E_pi_k -
                    (1/2)*(digamma(A_q_sigmasq[k]) + log(B_q_sigmasq[k]) +
                               (1/E_sigma_sq)*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k]))
            }
            tilde_phi_q_phi[, b] <- digamma(alpha_q_phi[b]) + digamma(sum(alpha_q_phi)) +
                rowSums(dens_phi_q_phi[b, , ])
        }

        phi_q_phi <- multinomial_logit(tilde_phi_q_phi)

        # expected count of observations in each w and z cluster
        E_n <- matrix(NA, nrow = K, ncol = B)
        for(b in 1:B) {
            for(k in 1:K) {
                E_n[k, b] <- sum(pi_q_pi[b, , k] * phi_q_phi[, b])
            }
        }

        E_w_n <- colSums(E_n) # expected count of observations in each w cluster
        E_z_n <- rowSums(E_n) # expected count of observations in each z cluster

        ## optimal density for \phi
        for(b in 1:B) {
            alpha_q_phi[b] <- 1 + sum(phi_q_phi[, b])
        }

        wts <- matrix(0, nrow = N, ncol = K)
        for (b in 1:B) {
            wts <- wts + pi_q_pi[b, , ] * phi_q_phi[, b]
        }

        for(k in 1:K){

            ## optimal density for V_k
            for (b in 1:B) {
                if(k < K) {
                    alpha_q_V[k, b] <- 1 + E_n[k, b]
                    beta_q_V[k, b] <- beta + sum(E_n[(k+1):K, b])
                }
            }

            ## optimal density for mu_k
            sigmasq_q_mu[k] <- 1/(1/sigmasq0 + (1/E_sigma_sq[k])*E_z_n[k])
            mu_q_mu[k] <- sigmasq_q_mu[k]*(mu0/sigmasq0 + (1/E_sigma_sq[k])*sum(y*wts[, k]))

            ## optimal density for sigma^2_k
            A_q_sigmasq[k] <- A0 + E_z_n[k]/2
            B_q_sigmasq[k] <- B0 + (1/2)*(sum(((y - mu_q_mu[k])^2 + sigmasq_q_mu[k])*wts[, k]))
        }

    }

    return(
        list(
            pi_q_pi = pi_q_pi,
            phi_q_phi = phi_q_phi,
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
