mcmc_plot <- function(output_mcmc, type = "hist") {
    # params <- names(output_mcmc)
    #
    # ## clean mu data
    # mu_data <- as.data.frame(output_mcmc$mu_save)
    # gathercols <- colnames(mu_data) <- paste0("mu_", 1:ncol(mu_data))
    # mu_data["Sample"] <- 1:nrow(mu_data)
    # mu_data <- tidyr::gather(mu_data, Parameter, Sample, gathercols, factor_key=TRUE)
    #
    # ## clean sigma_sq data
    # sigmasq_data <- as.data.frame(output_mcmc$sigmasq_save)
    # gathercols <- colnames(sigmasq_data) <- paste0("sigmasq_", 1:ncol(sigmasq_data))
    # sigmasq_data["Sample"] <- 1:nrow(sigmasq_data)
    # sigmasq_data <- tidyr::gather(sigmasq_data, Parameter, gathercols, factor_key=TRUE)
    #
    # ## clean beta data
    # beta_data <- as.data.frame(output_mcmc$beta_save)
    # gathercols <- colnames(beta_data) <- paste0("beta_", 1:ncol(beta_data))
    # beta_data["Sample"] <- 1:nrow(beta_data)
    # beta_data <- tidyr::gather(beta_data, Parameter, Sample, gathercols, factor_key=TRUE)
    #
    # ## clean phi data
    # phi_data <- as.data.frame(output_mcmc$phi_save)
    # gathercols <- colnames(phi_data) <- paste0("phi_", 1:ncol(phi_data))
    # phi_data["Sample"] <- 1:nrow(phi_data)
    # phi_data <- tidyr::gather(phi_data, Parameter, Sample, gathercols, factor_key=TRUE)
    #
    # ## clean pi data
    # B <- dim(output_mcmc$pi_save)[3]
    # for(b in 1:B) {
    #     pi_data <- as.data.frame(output_mcmc$pi_save[, , b])
    #     gathercols <- colnames(pi_data) <- paste0("pi_", 1:ncol(pi_data))
    #     pi_data["Sample"] <- 1:nrow(pi_data)
    #     pi_data <- tidyr::gather(pi_data, Parameter, Sample, gathercols, factor_key=TRUE)
    #
    #     if(type == 'hist') {
    #         pi_data <- ggplot(pi_data, aes(x = Sample)) +
    #             geom_histogram() +
    #             facet_wrap(~Parameter, scales = "free")
    #
    #         print(pi_data)
    #     }
    # }
    #
    # if(type == 'hist') {
    #     mu_plot <- ggplot(mu_data, aes(x = Sample)) +
    #         geom_histogram() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(mu_plot)
    #
    #     sigmasq_plot <- ggplot(sigmasq_data, aes(x = Sample)) +
    #         geom_histogram() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(sigmasq_plot)
    #
    #     beta_plot <- ggplot(beta_data, aes(x = Sample)) +
    #         geom_histogram() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(beta_plot)
    #
    #     phi_plot <- ggplot(phi_data, aes(x = Sample)) +
    #         geom_histogram() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(phi_plot)
    # } elseif(type = "trace") {
    #
    #     mu_plot <- ggplot(mu_data, aes(x = 1:)) +
    #         geom_line() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(mu_plot)
    #
    #     sigmasq_plot <- ggplot(sigmasq_data, aes(x = Sample)) +
    #         geom_histogram() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(sigmasq_plot)
    #
    #     beta_plot <- ggplot(beta_data, aes(x = Sample)) +
    #         geom_histogram() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(beta_plot)
    #
    #     phi_plot <- ggplot(phi_data, aes(x = Sample)) +
    #         geom_histogram() +
    #         facet_wrap(~Parameter, scales = "free")
    #
    #     print(phi_plot)
    # }

}
