#' Title
#'
#' @param mixture_model
#'
#' @return
#' @export
#'
#' @examples
vi_dpgmm <- function(
    mixture_model
) {

    # unpack model
    y <- mixture_model$y
    B <- mixture_model$B
    K <- mixture_model$K
    max_iter <- mixture_model$max_iter
    prior_parameters <- mixture_model$prior_parameters

    if (is.null(prior_parameters)) {

        mu_mu_k <- 0
        sigmasq_mu_k <- 100
        alpha_sigmasq_k <- 0.01
        beta_sigmasq_k <- 0.01
        a_lambda <- 0.01
        b_lambda <- 0.01

        prior_parameters <- list(
            mu_mu_k = mu_mu_k,
            sigmasq_mu_k = sigmasq_mu_k,
            alpha_sigmasq_k = alpha_sigmasq_k,
            beta_sigmasq_k = beta_sigmasq_k,
            a_lambda = a_lambda,
            b_lambda = b_lambda
        )

        mixture_model[["prior_parameters"]] <- prior_parameters
    }

    N <- length(y) # sample size
    elbos <- rep(-Inf, max_iter)

    ## initialize parameters
    tilde_p_q_p <- matrix(0, nrow = N, ncol = K)
    p_q_p <- matrix(0, nrow = N, ncol = K)

    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- runif(K)

    A_q_sigmasq <- runif(K)
    B_q_sigmasq <- runif(K)

    alpha_q_V <- runif(K)
    beta_q_V <- runif(K)
    alpha_q_V[K] <- 1
    beta_q_V[K] <- 1

    alpha_q_lambda <- runif(1)
    beta_q_lambda <- runif(1)

    for(m in 1:max_iter) {
        if(m %% 100 == 0) {
            message("Iteration ", m, " of ", max_iter)
        }

        ## calculate some expectations
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq
        E_log_sigmasq <- log(B_q_sigmasq) - digamma(A_q_sigmasq)
        E_log_p <- stick_breaking_expectation(B, K, alpha_q_V, beta_q_V)

        ## optimal density for z_ik
        for(k in 1:K) {
            tilde_p_q_p[ , k] <- (E_log_p[k] - (1/2)^log(2*pi) -
                                        (1/2)*(E_log_sigmasq[k]) -
                                        (1/2)*(1/E_sigma_sq[k])*((y-mu_q_mu[k])^2 + sigmasq_q_mu[k]))

        }
        p_q_p <- multinomial_logit(tilde_p_q_p[ , ])
        E_n <- colSums(p_q_p) # a count of observations in each group

        for(k in 1:K){

            ## optimal density for V_k
            if(k < K) {
                alpha_q_V[k] <- 1 + E_n[k]
                beta_q_V[k] <- alpha_q_lambda/beta_q_lambda + sum(E_n[(k+1):K])
            }

            ## optimal density for mu_k
            sigmasq_q_mu[k] <- 1/(1/sigmasq_mu_k + (1/E_sigma_sq[k])*E_n[k])
            mu_q_mu[k] <- sigmasq_q_mu[k]*(mu_mu_k/sigmasq_mu_k + (1/E_sigma_sq[k])*sum(y*p_q_p[, k]))

            ## optimal density for sigma^2_k
            A_q_sigmasq[k] <- alpha_sigmasq_k + E_n[k]/2
            B_q_sigmasq[k] <- beta_sigmasq_k + (1/2)*(sum(((y - mu_q_mu[k])^2 + sigmasq_q_mu[k])*p_q_p[, k]))
        }

        alpha_q_lambda <- a_lambda + N - 1
        beta_q_lambda <- b_lambda - sum(digamma(beta_q_V[1:(K-1)]) -
                                            digamma(alpha_q_V[1:(K-1)] + beta_q_V[1:(K-1)]))


        variational_parameters <- list(
            tilde_p_q_p = tilde_p_q_p,
            p_q_p = p_q_p,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            alpha_q_lambda = alpha_q_lambda,
            beta_q_lambda = beta_q_lambda
        )

        mixture_model[["variational_parameters"]] <- variational_parameters

        elbos[m] <- elbo_dpgmm(mixture_model)
    }

    return(
        list(
            p_q_p = p_q_p,
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            elbos = elbos,
            mixture_model = mixture_model
        )
    )
}
