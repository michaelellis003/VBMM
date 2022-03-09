#' Title
#'
#' @param mixture_model
#'
#' @return
#' @export
#'
#' @examples
vi_mdpgmm <- function(
    mixture_model
) {

    # unpack model
    y <- mixture_model$y
    w <- mixture_model$w
    B <- mixture_model$B
    K <- mixture_model$K
    max_iter <- mixture_model$max_iter
    prior_parameters <- mixture_model$prior_parameters

    if (is.null(prior_parameters)) {

        mu_mu_k <- 0
        sigmasq_mu_k <- 100
        alpha_sigmasq_k <- 0.1
        beta_sigmasq_k <- 0.1
        a_lambda <- 0.1
        b_lambda <- 0.1
        alpha_phi <- 0.5


        prior_parameters <- list(
            mu_mu_k = mu_mu_k,
            sigmasq_mu_k = sigmasq_mu_k,
            alpha_sigmasq_k = alpha_sigmasq_k,
            beta_sigmasq_k = beta_sigmasq_k,
            a_lambda = a_lambda,
            b_lambda = b_lambda,
            alpha_phi = alpha_phi
        )

        mixture_model[["prior_parameters"]] <- prior_parameters
    }

    N <- length(y) # sample size
    elbos <- rep(-Inf, max_iter)

    ## initialize parameters
    tilde_p_q_p <- array(0, dim = c(B, N, K))
    p_q_p <- array(0, dim = c(B, N, K))

    alpha_q_V <- matrix(1, nrow = K-1, ncol = B)
    beta_q_V <- matrix(rbeta((K-1)*B, 1, 1), nrow = K-1, ncol = B)

    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- runif(K, 0, 1)

    A_q_sigmasq <- runif(K, 0, 1)
    B_q_sigmasq <- runif(K, 0, 1)

    alpha_q_phi <- runif(B, 0, 1)

    phi_q_phi <- matrix(0, nrow = N, ncol = B)
    phi_q_phi[cbind(1:N, w)] <- 1

    alpha_q_lambda <- runif(B)
    beta_q_lambda <- runif(B)

    for(m in 1:max_iter) {

        if(m %% 100 == 0) {
            message("Iteration ", m, " of ", max_iter)
        }

        ## calculate expectations
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq
        E_log_sigma_sq <- log(B_q_sigmasq) - digamma(A_q_sigmasq)
        E_log_p <- stick_breaking_expectation(B, K, alpha_q_V, beta_q_V)
        E_log_phi <- digamma(alpha_q_phi) + digamma(sum(alpha_q_phi))

        ## optimal density for z_ik
        for(b in 1:B) {
            for(k in 1:K) {
                tilde_p_q_p[b, , k] <- (E_log_p[k, b] - (1/2)^log(2*pi) -
                    (1/2)*(E_log_sigma_sq[k]) -
                    (1/2)*(1/E_sigma_sq[k])*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k]))*(w == b)

            }
            p_q_p[b, , ] <- multinomial_logit(tilde_p_q_p[b, , ])
        }

        E_n <- matrix(NA, nrow = K, ncol = B) # expected count of observations in w and z cluster
        wts <- matrix(0, nrow = N, ncol = K)
        for(b in 1:B) {
            E_n[, b] <- t(phi_q_phi[, b]) %*% p_q_p[b, , ]
            wts <- wts + p_q_p[b, , ] * phi_q_phi[, b]

            ## optimal density for \phi
            alpha_q_phi[b] <- 1 + sum(phi_q_phi[, b])

            ## optimal density for \lambda
            Nb = sum(w == b)
            alpha_q_lambda[b] <- a_lambda + Nb - 1
            beta_q_lambda[b] <- b_lambda - sum(digamma(beta_q_V[1:(K-1), b]) -
                                                digamma(alpha_q_V[1:(K-1), b] + beta_q_V[1:(K-1), b]))
        }

        E_w_n <- colSums(E_n) # expected count of observations in each w cluster
        E_z_n <- rowSums(E_n) # expected count of observations in each z cluster

        for(k in 1:K){

            ## optimal density for V_k
            for (b in 1:B) {
                if(k < K) {
                    alpha_q_V[k, b] <- 1 + E_n[k, b]
                    beta_q_V[k, b] <- alpha_q_lambda[b]/beta_q_lambda[b] + sum(E_n[(k+1):K, b])
                }
            }

            ## optimal density for mu_k
            sigmasq_q_mu[k] <- 1/(1/sigmasq_mu_k + (1/E_sigma_sq[k])*sum(E_n[k, ]))
            mu_q_mu[k] <- sigmasq_q_mu[k]*(mu_mu_k/sigmasq_mu_k + (1/E_sigma_sq[k])*sum(y*wts[, k]))

            ## optimal density for sigma^2_k
            A_q_sigmasq[k] <- alpha_sigmasq_k + sum(E_n[k, ])/2
            B_q_sigmasq[k] <- beta_sigmasq_k + (1/2)*(sum(((y - mu_q_mu[k])^2 + sigmasq_q_mu[k])*wts[, k]))
        }

        variational_parameters <- list(
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            alpha_q_phi = alpha_q_phi,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            p_q_p = p_q_p,
            phi_q_phi = phi_q_phi
        )

        mixture_model[["variational_parameters"]] <- variational_parameters

        # elbos[m] <- elbo(y, w, B, K,
        #                  variational_parameters,
        #                  prior_parameters)

    }

    return(
        list(
            p_q_p = p_q_p,
            phi_q_phi = phi_q_phi,
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            alpha_q_phi = alpha_q_phi,
            elbos = elbos,
            mixture_model = mixture_model
        )
    )
}
