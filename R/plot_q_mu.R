#' Title
#'
#' @param vi_output
#' @param type
#'
#' @return
#' @export
#'
#' @examples
plot_q_mu <- function(vi_output,
                      type = "boxplot") {

    mu <- vi_output$mu_q_mu
    sigmasq <- vi_output$sigmasq_q_mu
    K <- length(mu)
    B <- vi_output$mixture_model$B

    df <- data.frame()
    for (k in 1:K) {
        df <- rbind(
            df,
            data.frame(
                K = factor(rep(k, 1000)),
                mu = rnorm(1000, mean = mu[k], sd = sqrt(sigmasq[k]))
            )
        )
    }

    if (type == "boxplot") {
        p <- ggplot(df, aes(x=K, y=mu)) +
            geom_boxplot()
    } else if (type == "histogram") {
        p <- ggplot(df) +
            geom_histogram(aes(x=mu, color = K, fill = K)) +
            facet_grid(~K)
    } else {
        stop("Type must be boxplot or histogram")
    }

    return(p)
}
