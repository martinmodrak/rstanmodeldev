#' @importFrom magrittr %>%
#' @export
sbc <- function(model, generator, N_steps, ...) {
    true_list <- list()
    observed_list <- list()
    for (i in 1:N_steps) {
        data <- generator()
        true_list[[i]] <- data$true
        observed_list[[i]] <- data$observed
    }

    fits <- sampling_multi(model, observed_list, ...)

    param_stats <- fits %>% purrr::imap_dfr(function(fit, data_id) {
        samples <- rstan::extract(fit, pars = names(true_list[[data_id]]))
        eval <- evaluate_all_params(samples, true_list[[data_id]])
        eval %>% dplyr::mutate(run = data_id)
    })
    diagnostics <- fits %>% purrr::imap_dfr(function(fit, run_id) {
        data.frame(run = run_id, n_divergent = rstan::get_num_divergent(fit), n_treedepth = rstan::get_num_max_treedepth(fit),
            n_chains_low_bfmi = length(rstan::get_low_bfmi_chains(fit)), total_time = sum(rstan::get_elapsed_time(fit)), min_n_eff = min(rstan::summary(fit)$summary[,
                "n_eff"], na.rm = TRUE), max_Rhat = max(rstan::summary(fit)$summary[, "Rhat"], na.rm = TRUE))
    })


    #Try to recover memory
    rm(fits);
    gc()

    return(list(params = param_stats, diagnostics = diagnostics, data = observed_list, true_values = true_list))
}

#' @export
#' @importFrom magrittr %>%
plot_sbc_histogram <- function(params, binwidth = 10, caption = NULL) {
    # CI - taken from https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R

    if (100%%binwidth != 0) {
        stop("binwidth has to divide 100")
    }
    n_simulations <- length(unique(params$run))
    CI <- qbinom(c(0.005, 0.5, 0.995), size = n_simulations, prob = binwidth/100)
    ci_lower <- CI[1]
    ci_mean <- CI[2]
    ci_upper <- CI[3]

    # The visualisation style taken as well from
    # https://github.com/seantalts/simulation-based-calibration/blob/master/Rsbc/generate_plots_sbc_inla.R
    params %>% ggplot(aes(x = order_within)) + geom_segment(aes(x = 0, y = ci_mean, xend = 100, yend = ci_mean), colour = "grey25") +
        geom_polygon(data = data.frame(x = c(-10, 0, -10, 110, 100, 110, -10), y = c(ci_lower, ci_mean, ci_upper, ci_upper,
            ci_mean, ci_lower, ci_lower)), aes(x = x, y = y), fill = "grey45", color = "grey25", alpha = 0.5) + geom_histogram(breaks = seq(1,
        101, by = binwidth), closed = "left", fill = "#A25050", colour = "black") + facet_wrap(~param_name, scales = "free_y") +
        ggtitle("Posterior order within 100 samples")

}

#' @export
#' @importFrom magrittr %>%
plot_sbc_scatter <- function(params, caption = NULL, plot_stat = "median", x_axis_trans = "identity", y_axis_trans = "identity") {
    if (!plot_stat %in% c("median", "mean")) {
        stop("Invalid plot_stat")
    }

    n_simulations <- length(unique(params$run))

    point_alpha <- 1/((n_simulations * 0.03) + 1)
    params %>% ggplot(aes_string(x = "true_value", y = plot_stat)) + geom_point(alpha = point_alpha) + geom_abline(slope = 1,
                                                                                                                         intercept = 0, color = "blue") + facet_wrap(~param_name, scales = "free") + scale_x_continuous(trans = x_axis_trans) +
              scale_y_continuous(trans = y_axis_trans) + ggtitle(paste0(plot_stat, " of marginal posteriors vs. true value - ", caption))
}

summarise_sbc_diagnostics <- function(sbc_results) {
    sbc_results$diagnostics %>% summarise(has_divergence = mean(n_divergent > 0), has_treedepth = mean(n_treedepth > 0), has_low_bfmi = mean(n_chains_low_bfmi >
        0), median_total_time = median(total_time), low_neff = mean(min_n_eff < 100), high_Rhat = mean(max_Rhat > 1.1))
}

sbc_power <- function(sbc_params, effect_threshold = 0.1) {
    sbc_params %>% mutate(sign_determined = sign(q_lower * q_upper) == 1, sign_determined_correct = sign_determined & sign(q_lower) ==
        sign(true_value), effect_above_threshold = sign_determined_correct & (q_lower > effect_threshold | q_upper < -effect_threshold)) %>%
        gather("stat", "value", sign_determined, sign_determined_correct, effect_above_threshold) %>% ggplot(aes(x = abs(true_value),
        y = value)) + geom_jitter(width = 0, height = 0.1, alpha = 0.1) + geom_smooth()
}
