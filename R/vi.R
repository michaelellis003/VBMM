#' Title
#'
#' @param mixture_model
#'
#' @return
#' @export
#'
#' @examples
vi <- function(
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

    if (B == 1) {
        Nw <- N
        w_mat <- matrix(1, nrow = N, ncol = B)
    } else {
        w_mat <- matrix(0, nrow = N, ncol = B)
        w_mat[cbind(1:N, w)] <- 1
        Nw <- colSums(w_mat)
    }

    ## initialize parameters
    alpha_q_V <- matrix(1, nrow = K-1, ncol = B)
    beta_q_V <- matrix(rbeta((K-1)*B, 1, 1), nrow = K-1, ncol = B)

    mu_q_mu <- rnorm(K, mean = mean(y), sd = sd(y))
    sigmasq_q_mu <- runif(K, 0, 1)

    A_q_sigmasq <- runif(K, 0, 1)
    B_q_sigmasq <- runif(K, 0, 1)

    ## optimal density for \lambda
    alpha_q_lambda <- a_lambda + Nw - 1
    beta_q_lambda <- runif(B, 0, 1)

    ## optimal density for \phi
    alpha_q_phi <- 1 + Nw

    ## calculate expectations
    E_sigma_sq <- B_q_sigmasq/A_q_sigmasq
    E_log_sigma_sq <- log(B_q_sigmasq) - digamma(A_q_sigmasq)
    E_log_V <- digamma(alpha_q_V) - digamma(alpha_q_V + beta_q_V) # E[log V_bk]
    E_log_1_V <- digamma(beta_q_V) - digamma(alpha_q_V + beta_q_V) # E[log 1-V_bk]
    E_log_p <- stick_breaking_expectation(B, K, E_log_V, E_log_1_V)

    for(m in 1:max_iter) {

        if(m %% 100 == 0) {
            message("Iteration ", m, " of ", max_iter)
        }

        log_tilde_p_q_p<- array(0, dim = c(B, N, K))
        p_q_p <- array(0, dim = c(B, N, K))
        E_n <- matrix(0, nrow = K, ncol = B) # expected count of observations in w and z cluster
        E_RSS <- matrix(0, nrow = N, ncol = K)
        for(b in 1:B) {

            for(k in 1:K) {
                ## optimal density for z_ik
                E_RSS[, k] <- (y-mu_q_mu[k])^2 + sigmasq_q_mu[k]
                log_likelihood <- -log(2*pi)/2 - E_log_sigma_sq[k]/2 -
                    E_RSS[, k]/(2*E_sigma_sq[k])

                log_tilde_p_q_p[b, , k] <- (E_log_p[k, b] + log_likelihood)*w_mat[, b]

            }
            p_q_p[b, , ] <- multinomial_logit(log_tilde_p_q_p[b, , ])

            # update expectation
            E_n[, b] <- t(w_mat[, b]) %*% p_q_p[b, , ]

            ## optimal density for \lambda
            beta_q_lambda[b] <- b_lambda - sum(E_log_1_V[, b])
        }

        ## lambda expectation
        E_lambda <- alpha_q_lambda/beta_q_lambda

        for(k in 1:K){

            ## optimal density for V_k
            for (b in 1:B) {
                if(k < K) {
                    alpha_q_V[k, b] <- 1 + E_n[k, b]
                    beta_q_V[k, b] <- E_lambda[b] + sum(E_n[(k+1):K, b])
                }
            }

            if ( B == 1) {
                wts <-p_q_p[, , k]
            } else {
                wts <-rowSums(t(p_q_p[, , k]) * w_mat)
            }

            ## optimal density for mu_k
            sigmasq_q_mu[k] <- 1/(1/sigmasq_mu_k + (1/E_sigma_sq[k])*sum(E_n[k, ]))
            mu_q_mu[k] <- sigmasq_q_mu[k]*(mu_mu_k/sigmasq_mu_k + (1/E_sigma_sq[k])*sum(y*wts))

            ## optimal density for sigma^2_k
            A_q_sigmasq[k] <- alpha_sigmasq_k + sum(E_n[k, ])/2
            B_q_sigmasq[k] <- beta_sigmasq_k + (1/2)*(sum(((y - mu_q_mu[k])^2 + sigmasq_q_mu[k])*wts))
        }

        ## update expectations
        E_sigma_sq <- B_q_sigmasq/A_q_sigmasq
        E_log_sigma_sq <- log(B_q_sigmasq) - digamma(A_q_sigmasq)
        E_log_V <- digamma(alpha_q_V) - digamma(alpha_q_V + beta_q_V) # E[log V_bk]
        E_log_1_V <- digamma(beta_q_V) - digamma(alpha_q_V + beta_q_V) # E[log 1-V_bk]
        E_log_p <- stick_breaking_expectation(B, K, E_log_V, E_log_1_V)

        expected_values <- list(
            E_sigma_sq = E_sigma_sq,
            E_log_sigma_sq = E_log_sigma_sq,
            E_log_V = E_log_V,
            E_log_1_V = E_log_1_V,
            E_log_p = E_log_p,
            E_lambda = E_lambda
        )

        variational_parameters <- list(
            alpha_q_V = alpha_q_V,
            beta_q_V = beta_q_V,
            mu_q_mu = mu_q_mu,
            sigmasq_q_mu = sigmasq_q_mu,
            A_q_sigmasq = A_q_sigmasq,
            B_q_sigmasq = B_q_sigmasq,
            p_q_p = p_q_p,
            log_tilde_p_q_p = log_tilde_p_q_p,
            alpha_q_lambda = alpha_q_lambda,
            beta_q_lambda = alpha_q_lambda,
            alpha_q_phi = alpha_q_phi
        )

        mixture_model[["variational_parameters"]] <- variational_parameters
        mixture_model[["expected_values"]] <- expected_values

        elbos[m] <- elbo(mixture_model)

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
            alpha_q_phi = alpha_q_phi,
            elbos = elbos,
            mixture_model = mixture_model
        )
    )
}
