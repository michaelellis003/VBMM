#' Title
#'
#' @param y
#' @param B
#' @param K
#' @param max_iter
#' @param w
#' @param prior_parameters
#' @param initial_values
#'
#' @return
#' @export
#'
#' @examples
vbmm <- function(
    y,
    B = 2,
    K = 20,
    max_iter = 30,
    w = NULL,
    prior_parameters = NULL,
    initial_values = NULL
) {

    mixture_model <- list()

    if (B%%1 != 0  | B < 1 ) {
        stop("B must be an integer greater than 0")
    } else {

        if (B == 1) {

            N <- length(y)
            mixture_model["Type"] <- "DPGMM"

            if (!is.null(w)) {
                stop("if B is 1 then w must be NULL")
            } else {
                mixture_model[["w"]] <- rep(1, N)
            }

        } else {

            mixture_model[["w"]] <- w
            mixture_model["Type"] <- "MDPGMM"

        }

        mixture_model["B"] <- B
    }

    if (K%%1 != 0  | K < 2 ) {
        stop("K must be an integer greater than 1")
    } else {
        mixture_model["K"] <- K
    }


    mixture_model[["prior_parameters"]] <- NULL
    mixture_model[["initial_values"]] <- NULL
    mixture_model[["y"]] <- y
    mixture_model[["max_iter"]] <- max_iter
    mixture_model[["variational_parameters"]] <- NULL

    output <- vi(mixture_model)

    return(output)
}
