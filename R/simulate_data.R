#' Simulate Data From Mixture Model
#'
#' @param N
#' @param K
#' @param mu
#' @param sigmasq
#' @param probs
#' @param plot
#'
#' @return
#' @export
#'
#' @examples
simulate_data <- function(
    N = 100,    # sample size
    K = 3,      # Number of mixtures
    mu = c(-20, 5, 25),
    sigmasq = c(0.5, 0.5, 0.5),
    probs = c(1/3, 1/3, 1/3),
    plot = TRUE
) {


    sample_clusters <- sample(1:K, size = N, replace = TRUE)
    y <- rnorm(N, mean = mu[sample_clusters], sd = sqrt(sigmasq[sample_clusters]))

    if(plot) {
        df <- data.frame(
            Sample = 1:N,
            y = y
        )

        ggplot2::ggplot(df, aes(x = y)) +
            geom_density(kernel = "gaussian")
    }

    return(
        list(
            y = y,
            N = N,
            K = K,
            mu = mu,
            sigmasq = sigmasq,
            probs = probs
        )
    )
}
