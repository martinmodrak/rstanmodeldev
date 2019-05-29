missing_arg <- function() quote(expr = )

order_within_samples <- function(x, samples, n_samples = 99) {
    thinned_samples <- sample(samples, n_samples)
    sum(thinned_samples < x) + 1
}

evaluate_single_param_indices <- function(
  samples, param_name, indices, true_value,
  q_probs = c(0.025, 0.975)) {
    if (is.null(indices)) {
        param_samples <- samples[[param_name]]
    } else {
        # A magic form to run samples[[param_name]][,indices[1], ... , indices[N]]
        # based on the length of indices
        param_samples <- do.call(`[`,
                                 append(list(samples[[param_name]], missing_arg()), indices))
    }

    if (length(indices) > 0) {
        indices_str <- do.call(paste, append(indices, list(sep = ",")))
        full_name <- paste0(param_name, "[", indices_str, "]")
    } else {
        full_name <- param_name
    }

    return(tibble(param_name = full_name,
                  true_value = true_value,
                  median = median(param_samples),
                  IQR = IQR(param_samples),
                  quantile = ecdf(param_samples)(true_value),
                  order_within = order_within_samples(true_value, param_samples),
                  q_lower = quantile(param_samples,q_probs[1]),
                  q_upper = quantile(param_samples, q_probs[2]))
           )
}

#' @export
evaluate_single_param <- function(samples, param_name, param_values) {
    result <-  list()
    dimensions <- dim(samples[[param_name]])[-1]  #The first dimension is the number of samples
    num_dimensions <- length(dimensions)
    next_element <- 1
    if (num_dimensions == 0) {
        result[[next_element]] <-
          evaluate_single_param_indices(
            samples, param_name, NULL, param_values)
        next_element <- next_element + 1
    } else if (num_dimensions == 1) {
        for (i in 1:dimensions[1]) {
            result[[next_element]] <-
              evaluate_single_param_indices(
                samples, param_name, list(i), param_values[i])
            next_element <- next_element + 1
        }
    } else if (num_dimensions == 2) {
        for (i in 1:dimensions[1]) {
            for (j in 1:dimensions[2]) {
                result[[next_element]] <-
                  evaluate_single_param_indices(
                    samples, param_name, list(i, j), param_values[i, j])
                next_element <- next_element + 1
            }
        }
    } else {
        stop("3+ dimensional parameters not supported yet")
    }
    return(do.call(rbind, result))
}

#' @export
evaluate_all_params <- function(samples, true_params) {
    result <- list()
    next_element <- 1
    for (param_name in names(true_params)) {
        if (!param_name %in% names(samples)) {
            next
        }
        param_values <- get(param_name, true_params)
        result[[next_element]] <-
          evaluate_single_param(samples, param_name, param_values)
        next_element <- next_element + 1
    }
    return(do.call(rbind, result))
}

#' @export
evaluation_summary <- function(fit, true_params, print_params_result = TRUE) {
    samples <- rstan::extract(fit)
    eval_result <- evaluate_all_params(samples, true_params)

    if (print_params_result) {
        # Add convergence diagnostics
        diags <- rstan::summary(fit)$summary[eval_result$param_name, , drop = FALSE]

        eval_result <- eval_result %>%
          mutate(n_eff = diags[, "n_eff"], Rhat = diags[, "Rhat"])

        print(eval_result)
    }
    quantiles <- eval_result$quantile
    within25 <- mean(quantiles >= 0.375 & quantiles <= 0.625)
    within50 <- mean(quantiles >= 0.25 & quantiles <= 0.75)
    within95 <- mean(quantiles >= 0.025 & quantiles <= 0.975)
    cat("\nWithin 25% interval:", within25,
        "\nWithin 50% interval:", within50,
        "\nWithin 95% interval:", within95, "\n")
}

#' @export
average_sampling_time <- function(fits) {
    time_list <- lapply(fits, rstan::get_elapsed_time)
    all__times <- Reduce(rbind, time_list, array(0, c(0, 2)))
    warmup_times <- all__times[, "warmup"]
    sample_times <- all__times[, "sample"]
    return(list(total = mean(warmup_times + sample_times), sample = mean(sample_times)))
}

#' @export
launch_shinystan_nonblocking <- function(fit) {
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package \"future\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("shinystan", quietly = TRUE)) {
    stop("Package \"shinystan\" needed for this function to work. Please install it.",
         call. = FALSE)
  }


  future::plan(future::multisession)
  future::future(shinystan::launch_shinystan(fit))
}



# Check transitions that ended with a divergence
#' @export
check_div <- function(fit) {
    sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    divergent <- do.call(rbind, sampler_params)[, "divergent__"]
    n <- sum(divergent)
    N <- length(divergent)

    print(sprintf("%s of %s iterations ended with a divergence (%s%%)",
                  n, N, 100 * n / N))
    if (n > 0)
        print("  Try running with larger adapt_delta to remove the divergences")
}

# Check transitions that ended prematurely due to maximum tree depth limit
#' @export
check_treedepth <- function(fit, max_depth = 10) {
    sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    treedepths <- do.call(rbind, sampler_params)[, "treedepth__"]
    n <- length(treedepths[sapply(treedepths, function(x) x == max_depth)])
    N <- length(treedepths)

    print(sprintf("%s of %s iterations saturated the maximum tree depth of %s (%s%%)",
                  n, N, max_depth, 100 * n / N))
    if (n > 0)
        print("  Run again with max_depth set to a larger value to avoid saturation")
}

# Checks the energy Bayesian fraction of missing information (E-BFMI)
#' @export
check_energy <- function(fit) {
    sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    no_warning <- TRUE
    for (n in 1:length(sampler_params)) {
        energies <- sampler_params[n][[1]][, "energy__"]
        numer <- sum(diff(energies) ^ 2) / length(energies)
        denom <- var(energies)
        if (numer / denom < 0.2) {
            print(sprintf("Chain %s: E-BFMI = %s", n, numer / denom))
            no_warning <- FALSE
        }
    }
    if (no_warning)
        print("E-BFMI indicated no pathological behavior")
    else
      print("  E-BFMI below 0.2 indicates you may need to reparameterize your model")
}

# Checks the effective sample size per iteration
#' @export
check_n_eff <- function(fit) {
    fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
    N <- dim(fit_summary)[[1]]

    iter <- dim(rstan::extract(fit)[[1]])[[1]]

    no_warning <- TRUE
    for (n in 1:N) {
        ratio <- fit_summary[, "n_eff"][n] / iter
        if (!is.na(ratio) && ratio < 0.02) {
            print(sprintf("n_eff / iter for parameter %s is %s!",
                          rownames(fit_summary)[n], ratio))
            no_warning <- FALSE
        }
    }
    if (no_warning)
        print("n_eff / iter looks reasonable for all parameters")
    else
      print("  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated")
}

# Checks the potential scale reduction factors
#' @export
check_rhat <- function(fit) {
    fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
    N <- dim(fit_summary)[[1]]

    no_warning <- TRUE
    for (n in 1:N) {
        rhat <- fit_summary[, 6][n]
        if (!is.na(rhat) && (rhat > 1.1 || is.infinite(rhat))) {
            print(sprintf("Rhat for parameter %s is %s!",
                          rownames(fit_summary)[n], rhat))
            no_warning <- FALSE
        }
    }
    if (no_warning)
        print("Rhat looks reasonable for all parameters")
    else
      print("  Rhat above 1.1 indicates that the chains very likely have not mixed")
}

#' @export
check_all_diagnostics <- function(fit) {
    check_n_eff(fit)
    check_rhat(fit)
    check_div(fit)
    check_treedepth(fit)
    check_energy(fit)
}
