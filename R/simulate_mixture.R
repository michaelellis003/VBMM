#' Title
#'
#' @param N
#' @param B
#' @param K
#' @param mu
#' @param sigmasq
#' @param phi
#' @param p
#'
#' @return
#' @export
#'
#' @examples
simulate_mixture <- function(
    N = 500,
    B = 2,
    K = 3,
    mu = c(-20, 5, 25),
    sigmasq = c(0.5, 0.5, 0.5),
    phi = c(1/2, 1/2),
    p = list(c(1/3, 1/3, 1/3), c(1/2, 1/4, 1/4))
) {

    z <- rep(NA, N)
    w <- sample(1:B, size = N, replace = TRUE, prob = phi)
    y <- rep(0, N)

    for(i in 1:N) {
        z[i] <- sample(1:K, size = 1, replace = TRUE, prob = p[[w[i]]])
        y[i] <- rnorm(1, mean = mu[z[i]], sd = sqrt(sigmasq[z[i]]))
    }

    return(
        list(
            y = y,
            N = N,
            B = B,
            K = K,
            mu = mu,
            sigmasq = sigmasq,
            phi = phi,
            p = p,
            w = w,
            z = z
        )
    )
}
