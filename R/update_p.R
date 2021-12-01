#' Title
#'
#' @param mod
#' @param expectations
#' @param variational_parameters
#'
#' @return
#' @export
#'
#' @examples
update_p <- function(
    mod,
    expectations,
    variational_parameters
) {

    # unpack model
    y <- mod$y[,1]
    N <- mod$N
    B <- mod$B
    K <- mod$K

    # unpack expectations
    E_log_phi <- expectations$E_log_phi
    E_log_p <- expectations$E_log_p
    E_sigma_sq <- expectations$E_sigma_sq
    E_log_sigma_sq <- expectations$E_log_sigma_sq
    E_SS <- expectations$E_SS

    p_tilde <- array(NA, dim = c(N, B, K))
    log_mix_dens <- array(0, dim = c(N, B, K))
    p_z <- array(1/NA, dim = c(N, B, K))

    lnd <- log_normal_dens(E_SS,
                           E_log_sigma_sq,
                           E_sigma_sq)

    for(b in 1:B) {
        for(i in 1:N) {
            p_tilde[i, b, ] <- E_log_p[, b] - lnd[i, ] + E_log_phi[b]
        }
        p_z[, b, ] <- multinomial_logit(p_tilde[, b, ])
    }

    return(p_z)
}
