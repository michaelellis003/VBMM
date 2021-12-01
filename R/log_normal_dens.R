#' Title
#'
#' @param E_SS
#' @param E_log_sigma_sq
#' @param E_sigma_sq
#'
#' @return
#' @export
#'
#' @examples
log_normal_dens <- function(
    E_SS,
    E_log_sigma_sq,
    E_sigma_sq
) {

    return(
        (1/2)*(log(2*pi) + t(apply(E_SS, 1,
                    function(x) E_log_sigma_sq + (1/E_sigma_sq)*x)))
    )

}
