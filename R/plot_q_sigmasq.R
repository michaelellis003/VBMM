#' Title
#'
#' @param vi_output
#' @param type
#'
#' @return
#' @export
#'
#' @examples
plot_q_sigmasq <- function(vi_output,
                           type = "boxplot") {

    A <- vi_output$A_q_sigmasq
    B <- vi_output$B_q_sigmasq
    K <- length(A)
    B <- vi_output$mixture_model$B

    df <- data.frame()
    for (k in 1:K) {
        df <- rbind(
            df,
            data.frame(
                K = factor(rep(k, 1000)),
                sigmasq = 1/rgamma(1000, A[k], B[k])
            )
        )
    }

    if (type == "boxplot") {
        p <- ggplot(df, aes(x=K, y=sigmasq)) +
            geom_boxplot()
    } else if (type == "histogram") {
        p <- ggplot(df) +
            geom_histogram(aes(x=sigmasq, color = K, fill = K)) +
            facet_grid(~K)
    } else {
        stop("Type must be boxplot or histogram")
    }

    return(p)
}
