#' Title
#'
#' @param y
#' @param w
#' @param B
#' @param K
#' @param n_mcmc
#' @param burnin
#' @param n_message
#' @param prior_parameters
#'
#' @return
#' @export
#'
#' @examples
mcmc <- function(
    y,
    w,
    B = 2,
    K = 20,
    n_mcmc = 5000,
    burnin = 2500,
    n_message = 500,
    prior_parameters = NULL
) {

    N <- length(y) # sample size
    n_save <- n_mcmc-burnin

    ## prior parameters
    if(is.null(prior_parameters)) {
        mu0 <- 0
        sigmasq0 = 10^2
        A0 <- 0.1
        B0 <- 0.1
        alpha0 <- 1
        a0_beta <- 0.1
        b0_beta <- 0.1

        prior_parameters <- list(
            mu0 = mu0,
            sigmasq0 = sigmasq0,
            A0 = A0,
            B0 = B0,
            alpha0 = alpha0,
            a0_beta = a0_beta,
            b0_beta = b0_beta
        )
    }


    ## initial parameters
    mu <- rep(mu0, K)
    sigmasq <- rep(sigmasq0, K)

    phi <- runif(B, 0, 1)

    V <- matrix(0, nrow = K, ncol = B)
    beta <- runif(B, 0, 1)
    pi <- matrix(0, nrow = K, ncol = B)
    z <- array(0, dim = c(N, K, B))
    for(i in 1:N) {
        for(b in 1:B) {
            V[, b] <- c(rbeta(K-1, 1, beta[1]), 1)
            pi[, b] <- stick_breaking(V, K)

            ind <- sample(c(1:K), size = 1, replace = TRUE, prob = pi[, b])
            z[i, ind, b] <- 1
        }
    }

    ## encode w
    tmp <- matrix(0, nrow = N, ncol = B)
    for(i in 1:N) {
        tmp[i, w[i]] <- 1
    }
    w <- tmp

    ## save parameters
    mu_save <- matrix(NA, nrow = n_save, ncol = K)
    sigmasq_save <- matrix(NA, nrow = n_save, ncol = K)
    phi_save <- matrix(NA, nrow = n_save, ncol = B)
    beta_save <- matrix(NA, nrow = n_save, ncol = B)
    pi_save <- array(NA, dim = c(n_save, K, B))
    z_save <- array(NA, dim = c(n_save, N, K, B))

    for(i in 1:n_mcmc) {
        if(i %% n_message == 0) {
            message("Iteration ", i , " out of ", n_mcmc)
        }

        counts <- matrix(0, nrow = K, ncol = B)
        for(b in 1:B) {
            counts[, b] <- t(w[, b]) %*% z[, , b]
        }
        n_k <- rowSums(counts)

        for(k in 1:K) {

            wts <- rowSums(w * z[, k, ])

            ## update parameters for mu
            a <- 1/(1/sigmasq0 + n_k[k]/sigmasq[k])
            b <- mu0/sigmasq0 + sum(y*wts)/sigmasq[k]

            ## sample mu
            mu[k] <- rnorm(1, a*b, sqrt(a))

            ## updated parameters for sigma squared
            weighted_SS <- ((y - mu[k])^2)*wts

            updated_A <- A0 + n_k/2
            updated_B <- B0 + weighted_SS/2

            ## sample sigma squared
            sigmasq[k] <- 1/rgamma(1, updated_A, updated_B)

            ## updated parameters for pi tilde
            updated_shape <- 1 + counts[k, ]
            if( k != K) {
                updated_beta <- beta + colSums(counts[k:K, ])

            } else {
                updated_beta <- beta + counts[K, ]
            }

            ## sample V
            if(k < K) {
                V[k, ] <- mapply(rbeta, n = 1, shape1 = updated_shape,
                                 shape2 = updated_beta)
            }
        }

        ## update beta
        updated_a0_beta <- a0_beta + K - 1
        for(b in 1:B) {
            updated_b0_beta <- b0_beta - sum(log(1 - V[1:(K-1), b]))
            beta[b] <- rgamma(1, updated_a0_beta, updated_b0_beta)
        }

        ## update alpha for dirichlet
        updated_alpha <- alpha0 + colSums(w)
        phi <- MCMCpack::rdirichlet(1, updated_alpha)

        z_k_dens <- matrix(NA, nrow = N, ncol = K)
        z_probs <- array(NA, dim = c(N, K, B))
        for(b in 1:B) {
            ## update pi
            pi[, b] <- stick_breaking(V[, b], K)
            for(k in 1:K) {
                z_k_dens[, k] <- log(phi[b]) + log(pi[k, b]) +
                    dnorm(y, mu[k], sqrt(sigmasq[k]), log = TRUE)
            }

            z_probs[, , b] <- multinomial_logit(z_k_dens)
        }

        ## sample z
        z <- array(0, dim = c(N, K, B))
        for(j in 1:N) {
            for(b in 1:B) {
                ind <- sample(c(1:K), size = 1, replace = TRUE, prob = z_probs[j, , b])
                z[j, ind, b] <- 1
            }
        }

        if(i > burnin) {
            index <- i-burnin
            mu_save[index, ] <- mu
            sigmasq_save[index, ] <- sigmasq
            beta_save[index, ] <- beta
            phi_save[index, ] <- phi
            pi_save[index, ,] <- pi
            z_save[index, , ,] <- z
        }
    }

    return(
        list(
            mu_save = mu_save,
            sigmasq_save = sigmasq_save,
            beta_save = beta_save,
            pi_save = pi_save,
            phi_save = phi_save,
            z_save = z_save
        )
    )

}
