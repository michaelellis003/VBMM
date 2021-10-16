#' Variational Inference for Infinite Mixture Model
#'
#' @param y
#' @param K
#' @param mu0
#' @param sigmasq0
#' @param A
#' @param B
#' @param beta
#' @param max_iter
#'
#' @return
#' @export
#'
#' @examples
vi <- function(
    y,
    K = 20,
    mu0 = 0,
    sigmasq0 = 10^8,
    A = 0.01,
    B = 0.01,
    beta = 0.5,
    max_iter = 30
) {


    N <- length(y) # sample size
    sum_y <- sum(y)
    elbos <- rep(-Inf, max_iter)

    ## initialize parameters
    tilde_pi_q_pi <- matrix(0, nrow = N, ncol = K)
    pi_q_pi <- matrix(0, nrow = N, ncol = K)

    alpha_q_V <- rep(1, K-1)
    beta_q_V <- rbeta(K-1, 1, 1)

    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- 1/rgamma(K, 1, 1)

    A_q_sigmasq <- rep(1, K)
    B_q_sigmasq <- rep(1, K)

    for(m in 1:max_iter) {

        ## expectation of sigma^2
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq

        ## optimal density for z_ik
        for(k in 1:K) {
            ### stick-breaking representation
            E_pi_k <- 0
            if(k == 1) {
                E_pi_k <- digamma(alpha_q_V[k]) - digamma(alpha_q_V[k] + beta_q_V[k])
            } else if (k > 1 & k < K) {
                E_pi_k <- digamma(alpha_q_V[k]) - digamma(alpha_q_V[k] + beta_q_V[k]) +
                    sum(digamma(beta_q_V[1:(k-1)]) -
                            digamma(alpha_q_V[1:(k-1)] + beta_q_V[1:(k-1)]))
            } else { # if k == K then E[log V_K] = E[log 1] = 0
                E_pi_k <- sum(digamma(beta_q_V[1:(k-1)]) -
                                  digamma(alpha_q_V[1:(k-1)] + beta_q_V[1:(k-1)]))
            }

            tilde_pi_q_pi[, k] <- E_pi_k -
                (1/2)*(digamma(A_q_sigmasq[k]) + log(B_q_sigmasq[k]) +
                           (1/E_sigma_sq)*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k]))
        }

        pi_q_pi <- t(apply(tilde_pi_q_pi, 1, function(x)exp(x)/sum(exp(x))))

        E_n <- colSums(pi_q_pi) # a count of observations in each group

        for(k in 1:K){

            ## optimal density for V_k
            if(k < K) {
                alpha_q_V[k] <- 1 + E_n[k]
                beta_q_V[k] <- beta + sum(E_n[(k+1):K])
            }

            ## optimal density for mu_k
            sigmasq_q_mu[k] <- 1/(1/sigmasq0 + (1/E_sigma_sq[k])*E_n[k])
            mu_q_mu[k] <- sigmasq_q_mu[k]*(mu0/sigmasq0 + (1/E_sigma_sq[k])*sum(y*pi_q_pi[, k]))

            ## optimal density for sigma^2_k
            A_q_sigmasq[k] <- A + E_n[k]/2
            B_q_sigmasq[k] <- B + (1/2)*(sum(((y - mu_q_mu[k])^2 + sigmasq_q_mu[k])*pi_q_pi[, k]))
        }

    }

    return(
        list(
            pi_q_pi = pi_q_pi,
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            elbos = elbos
        )
    )
}
