#' @param models either a single stanmodel object which is used for
#'   all data or a list of models the same length as data
#' @param data a list of datasets to pass to each fit
#' @param map_fun a function to process each single chain fit. It has to take 3 params:
#'   a fit (single chain), data index and chain index. By default it is a noop.
#'   This function is run in parallel.
#'   Note that this function must be runnable in new RStudio sessions. You may
#'   use the `R_session_init` or `map_fun_dependencies` parameters to ensure it is loaded.
#' @param combine_fun a function called to combine the results of applying `map_fun`
#'   to individual chain fits. It should accept a list of values of the `map_fun` calls.
#'   Default value is [sflist2stanfit()] which creates a combined multichain fit.
#'   As of now, this function is run serially after all fits have bin computed.
#' @param ids_to_compute if given the computation is done only for certain items in
#'   the data list (useful for restarting long computations after a failure)
#' @param init_per_item optional list of inits of the same length as data
#' @param control_per_item optional list of control arguments of the same length as data
#' @param map_fun_dependencies a list of package names that need te be loaded for
#'   `map_fun` to run. IMPORTANT: when developing packages, you need to install
#'   the latest version (the packages are loaded from the default library)
#' @param cache_dir if not NULL, fits will be cached in this directory
#' @param ... Other params to be passed the [sampling()].
#'
#' @return A list of length `length(data)` containing the result of applying
#'   `combine_fun` to each set of `chains` result of applying `map_fun`
#'   to each fit.
#' @export
sampling_multi <- function(models, data, map_fun = sampling_multi_noop, combine_fun = rstan::sflist2stanfit, chains = 4, cores = parallel::detectCores(),
    init = NULL, control = NULL, init_per_item = NULL, control_per_item = NULL, map_fun_dependencies = c(), R_session_init_expr = NULL,
    cache_dir = NULL, ids_to_compute = 1:length(data), cluster_options = list(), ...) {

    # TODO: argument validation
    if (!is.null(cache_dir) && !dir.exists(cache_dir)) {
        stop(paste0("Cache dir '", cache_dir, "'  does not exist"))
    }

    if(is.null(cluster_options$spec)) {
        cluster_options$spec <- cores
    }
    if(is.null(cluster_options$useXDR)) {
        cluster_options$useXDR <- FALSE
    }
    cl <- do.call(parallel::makeCluster, cluster_options)
    on.exit(parallel::stopCluster(cl))

    # Prepare argument lists for all calls
    .dotlists_per_item <- list()
    for (data_id in 1:length(data)) {
        if (is.list(models)) {
            model <- models[[data_id]]
        } else {
            model <- models
        }
        if (class(model) != "stanmodel") {
            stop(paste0("Model for data_id ", data_id, " is not of class 'stanmodel'"))
        }

        for (chain_id in 1:chains) {

            .dotlist <- c(list(model, chains = 1L, cores = 0L), list(...))

            if (is.list(init_per_item)) {
                .dotlist$init <- init_per_item[[data_id]]
                # TODO check that init is not set
            } else {
                .dotlist$init <- init
            }

            # Handle per-chain inits
            if (is.list(.dotlist$init))
                .dotlist$init <- .dotlist$init[chain_id]

            if (is.list(control_per_item)) {
                .dotlist$control <- control_per_item[[data_id]]
                # TODO check that control is not set
            } else {
                .dotlist$control <- control
            }

            .dotlist$chain_id <- chain_id

            .dotlist$data <- data[[data_id]]


            .dotlists_per_item[[(data_id - 1) * chains + chain_id]] <- .dotlist
        }
    }


    fit_fun <- function(i) {
        cached <- FALSE
        not_cached_msg <- ""
        arg_list <- .dotlists_per_item[[i]]
        if (!is.null(cache_dir)) {
            filename <- sprintf("%s/%06d.rds", cache_dir, i)
            if (file.exists(filename)) {
                cached_data <- readRDS(filename)

                if (!is.list(cached_data)) {
                  not_cached_msg <- "Unexpected data format"
                  cached <- FALSE
                } else {
                  fit_from_cache <- cached_data$fit

                  model <- arg_list[[1]]

                  if (class(model) != "stanmodel") {
                    stop("Argument 1 is not a stan model")
                  }

                  if (is.null(cached_data$data) || is.null(fit_from_cache)) {
                    not_cached_msg <- "Does not contain fit or data"
                    cached <- FALSE
                  } else if (class(fit_from_cache) != "stanfit") {
                    not_cached_msg <- "Not a stanfit"
                    cached <- FALSE
                  } else if (model@model_code != fit_from_cache@stanmodel@model_code) {
                    not_cached_msg <- "Model code changed"
                    cached <- FALSE
                  } else if (!identical(cached_data$data, arg_list$data)) {
                    not_cached_msg <- "Different data"
                    cached <- FALSE
                  } else {
                    single_fit <- fit_from_cache
                    cached <- TRUE
                  }
                }
            }
        }

        if (!cached) {
            single_fit <- do.call(rstan::sampling, args = arg_list)
            if (!is.null(cache_dir)) {
                saveRDS(list(data = arg_list$data, fit = single_fit), filename)
            }
        }

        # TODO should catch error from map_fun
        list(cached = cached, not_cached_msg = not_cached_msg, result = map_fun(single_fit, data_id, chain_id))
    }



    dependencies <- c("rstan", "Rcpp", map_fun_dependencies)
    .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
        dirname(system.file(package = d))
    })))
    .paths <- .paths[.paths != ""]
    parallel::clusterExport(cl, varlist = c(".paths", "dependencies"), envir = environment())
    parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
    parallel::clusterEvalQ(cl, expr = for (dep in dependencies) {
        suppressPackageStartupMessages(require(dep, quietly = TRUE, character.only = TRUE))
    })

    parallel::clusterExport(cl, varlist = ".dotlists_per_item", envir = environment())

    parallel::clusterExport(cl, varlist = c("map_fun", "R_session_init_expr"), envir = environment())

    parallel::clusterEvalQ(cl, expr = R_session_init_expr)

    ids <- rep(ids_to_compute, each = chains)
    items <- (rep(ids_to_compute, each = chains) - 1) * chains + rep(1:chains, times = length(ids_to_compute))

    results_raw <- parallel::parLapplyLB(cl, X = items, fun = fit_fun, chunk.size = 1)
    results_flat <- purrr::map(results_raw, function(x) {
        x$result
    })

    results <- list()
    if (!is.null(cache_dir)) {
        cache_info <- purrr::map_dfr(results_raw, function(x) {
            data.frame(cached = x$cached, msg = x$not_cached_msg, stringsAsFactors = FALSE)
        })
        cat(paste0(sum(cache_info$cached), " out of ", length(items), " chains read from cache\n"))
        if (any(cache_info$msg != "")) {
            cat("Reasons for ignoring cache:\n")
            message_stats <- aggregate(cached ~ msg, cache_info, FUN = length, subset = cache_info$msg != "")
            names(message_stats) <- c("message", "count")
            print(message_stats)
        }
    }
    for (i in 1:length(data)) {
        results[[i]] <- combine_fun(results_flat[ids == i])
    }
    results
}

#' @export
sampling_multi_noop <- function(fit, data_id, chain_id) {
    fit
}

#' @return Returns a function to be used as `map_fun` in [sampling_multi()]
#'   that stores all fits in RDS files
#' @export
sampling_multi_store_file_generator <- function(base_dir, base_name) {
    force(base_dir)
    force(base_name)
    function(fit, data_id, chain_id) {
        filename <- paste0(base_dir, "/", base_name, "_", data_id, "_", chain_id, ".rds")
        saveRDS(fit, filename)
        filename
    }
}
