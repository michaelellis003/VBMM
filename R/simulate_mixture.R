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

    sample_parents <- sample(1:B, size = N, replace = TRUE, prob = phi)
    y <- list()

    for(b in 1:B) {
        Nb <- sum(sample_parents == b)
        sample_clusters <- sample(1:K, size = Nb, replace = TRUE, prob = p[[b]])
        y[[b]] <- rnorm(Nb, mean = mu[sample_clusters], sd = sqrt(sigmasq[sample_clusters]))
    }

    y <- unlist(y)

    return(
        list(
            y = y,
            N = N,
            B = B,
            K = K,
            mu = mu,
            sigmasq = sigmasq,
            phi = phi,
            p = p
        )
    )
}
