#' Title
#'
#' @param y
#' @param w
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
VBMIM <- function(
    y,
    w,
    B = 2,
    K = 20,
    prior_parameters = NULL,
    max_iter = 30
) {

    N <- length(y) # sample size
    sum_y <- sum(y)
    elbos <- rep(-Inf, max_iter)

    if(is.null(prior_parameters)) {
        mu0 <- 0
        sigmasq0 <- 10^2
        A0 <- 0.01
        B0 <- 0.01
        beta <- 0.5

        prior_parameters <- list(
            mu0 = mu0,
            sigmasq0 = sigmasq0,
            A0 = A0,
            B0 = B0,
            beta = beta
        )
    }

    ## initialize parameters
    tilde_pi_q_pi <- array(0, dim = c(B, N, K))
    pi_q_pi <- array(0, dim = c(B, N, K))

    alpha_q_V <- matrix(1, nrow = K-1, ncol = B)
    beta_q_V <- matrix(rbeta((K-1)*B, 1, 1), nrow = K-1, ncol = B)

    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- runif(K, 0, 1)

    A_q_sigmasq <- runif(K, 0, 1)
    B_q_sigmasq <- runif(K, 0, 1)

    alpha_q_phi <- runif(B, 0, 1)

    phi_q_phi <- matrix(0, nrow = N, ncol = B)
    for(i in 1:N) {
        phi_q_phi[i, w[i]] <- 1
    }

    for(m in 1:max_iter) {

        if(m %% 100 == 0) {
            message("Iteration ", m, " of ", max_iter)
        }

        ## calculate expectations
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq
        E_log_sigma_sq <- log(B_q_sigmasq) - digamma(A_q_sigmasq)
        E_log_pi <- stick_breaking_expectation(B, K, alpha_q_V, beta_q_V)
        E_log_phi <- digamma(alpha_q_phi) + digamma(sum(alpha_q_phi))

        ## optimal density for z_ik
        for(b in 1:B) {
            for(k in 1:K) {
                # tilde_pi_q_pi[b, , k] <- E_pi[k, b] - (1/2)^log(2*pi) -
                #     (1/2)*(log(B_q_sigmasq[k]) - digamma(A_q_sigmasq[k])) -
                #     (1/2)*(1/E_sigma_sq[k])*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k]) +
                #     digamma(alpha_q_phi[b]) + digamma(sum(alpha_q_phi))

                tilde_pi_q_pi[b, , k] <- E_log_pi[k, b] - (1/2)^log(2*pi) -
                    (1/2)*(E_log_sigma_sq[k]) -
                    (1/2)*(1/E_sigma_sq[k])*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k]) +
                    E_log_phi[b]
            }
            pi_q_pi[b, , ] <- multinomial_logit(tilde_pi_q_pi[b, , ])
        }

        E_n <- matrix(NA, nrow = K, ncol = B) # expected count of observations in w and z cluster
        wts <- matrix(0, nrow = N, ncol = K)
        for(b in 1:B) {
            E_n[, b] <- t(phi_q_phi[, b]) %*% pi_q_pi[b, , ]

            ## optimal density for \phi
            alpha_q_phi[b] <- 1 + sum(phi_q_phi[, b])

            wts <- wts + pi_q_pi[b, , ] * phi_q_phi[, b]
        }

        E_w_n <- colSums(E_n) # expected count of observations in each w cluster
        E_z_n <- rowSums(E_n) # expected count of observations in each z cluster

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

        variational_parameters <- list(
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            alpha_q_phi = alpha_q_phi,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            pi_q_pi = pi_q_pi,
            phi_q_phi = phi_q_phi
        )

        elbos[m] <- elbo(y, w, B, K,
                         variational_parameters,
                         prior_parameters)

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
