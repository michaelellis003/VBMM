#' Title
#'
#' @param vi_output
#'
#' @return
#' @export
#'
#' @examples
plot_elbo <- function(vi_output) {

    elbos <- vi_output$elbos
    num_iter <- length(elbos)

    elbos <- data.frame(
        Iteration = 1:num_iter,
        ELBO = elbos
    )

    p <- ggplot(elbos, aes(x = Iteration, y = ELBO)) +
        geom_line() +
        geom_point()

    return(p)
}
