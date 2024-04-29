library(sdmTMB)
r2.sdmTMB <- function(x, which_fixef = NULL, method = NULL) {
  if (!inherits(x, "sdmTMB")) {
    cli::cli_abort("'x' must be a model of class sdmTMB.", call. = FALSE)
  }
  if (isTRUE(x$reml)) {
    cli::cli_abort("r2.sdmTMB() does not yet work with REML", call. = FALSE)
  }
  if (length(x$family$family) > 1) {
    cli::cli_abort("r2.sdmTMB() does not work for delta (hurdle) models yet.", call. = FALSE)
  }
  if (x$family$family == "student") {
    cli::cli_inform("Family is student, but the variance does not (yet) account for the degrees of freedom.")
  }
  # if (x$family$family == "binomial" & is.null(method)) {
  #   cli::cli_inform("`method` not specified, using the theoretical approach")
  # }
  # if (x$family$family == "tweedie" & is.null(method)) {
  #   cli::cli_inform("`method` not specified, using the lognormal approximation.")
  # }
  if (!x$family$family %in% c("student", "gaussian", "binomial", "tweedie", "Gamma", "poisson")) {
    cli::cli_abort("r2.sdmTMB() currently only works for Gaussian, binomial, Gamma, Poisson, and Tweedie models.", call. = FALSE)
  }
  if (!is.null(x$spatial_varying)) {
    cli::cli_abort("r2.sdmTMB() currently does not work with spatially varying coefficient models.", call. = FALSE)
  }
  
  lp <- x$tmb_obj$env$last.par.best
  r <- x$tmb_obj$report(lp)
  
  varF <- var(r$eta_fixed_i[, 1L]) # FIXME delta
  
  if (isTRUE(x$smoothers$has_smooths)) {
    varSmooths <- var(r$eta_smooth_i)
  } else {
    varSmooths <- 0
  }
  
  b <- tidy(x, "ran_par")
  sigma <- function(x) {
    .b <- tidy(x, "ran_par")
    .b$estimate[.b$term == "phi"]
  }
  
  varO <- varE <- varV <- varG <- 0
  if (x$tmb_data$include_spatial == 1L && x$tmb_data$no_spatial) {
    varO <- b$estimate[b$term == "sigma_O"]^2 # spatial variance
  }
  if (x$tmb_data$spatial_only == 0L && !x$tmb_data$no_spatial) {
    varE <- b$estimate[b$term == "sigma_E"]^2 # spatiotemporal variance
  }
  if (x$tmb_data$random_walk == 1L) {
    if (!identical(x$time_varying, ~1)) {
      cli::cli_abort("r2.sdmTMB() currently only works with time-varying intercepts.", call. = FALSE)
    }
    varV <- b$estimate[b$term == "sigma_V"]^2 # time-varying variance
  }
  if (x$tmb_data$nobs_RE > 0) {
    varG <- b$estimate[b$term == "sigma_G"]^2 # random effect variance
  }
  
  if (x$family$family %in% c("student", "gaussian")) {
    varR <- sigma(x)^2
    # denominator <- varF_all + varO + varE + varR + varV + varG
  } else if (x$family$family == "binomial" && x$family$link == "logit") {
    # "theoretical" method of Nakagawa supp. row 115
    varR <- pi^2 / 3 # FIXME: CHECK THIS
  } else if (x$family$family %in% c("tweedie", "Gamma", "poisson")) {
    varR <- get_distribution_variance(x)
  } else {
    cli::cli_abort("Family not implemented", call. = FALSE)
  }
  
  denominator <- varF + varSmooths + varO + varE + varR + varV + varG
  varF_all <- varF + varSmooths
  
  marginal <- varF_all / denominator
  
  cond_rf_sp <- cond_rf_spt <- cond_tv <- cond_re <- cond_all <- cond_smooth <- cond_fixed <- NULL
  if (varO != 0) {
    cond_rf_sp <- varO / denominator
  }
  if (varE != 0) {
    cond_rf_spt <- varE / denominator
  }
  if (varV != 0) {
    cond_tv <- varV / denominator
  }
  if (varG != 0) {
    cond_re <- varG / denominator
  }
  if (varSmooths != 0) {
    cond_smooth <- varSmooths / denominator
  }
  if (varF != 0) {
    cond_fixed <- varF / denominator
  }
  cond_all <- (denominator - varR) / denominator
  
  out <- list(
    conditional = cond_all,
    marginal = marginal,
    partial_smoothers = cond_smooth,
    partial_fixed = if (!is.null(cond_smooth)) cond_fixed else NULL,
    partial_spatial = cond_rf_sp,
    partial_spatiotemporal = cond_rf_spt,
    partial_time_varying = cond_tv,
    partial_random_intercepts = cond_re
  )
  out[vapply(out, is.null, logical(1L))] <- NULL
  ret <- t(as.data.frame(lapply(out, `[`, 1L)))
  out <- data.frame(component = row.names(ret), R2 = ret[, 1L, drop = TRUE], stringsAsFactors = FALSE)
  row.names(out) <- NULL
  
  if (requireNamespace("tibble", quietly = TRUE)) {
    out <- tibble::as_tibble(out)
  }
  out
}