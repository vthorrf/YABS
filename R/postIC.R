postIC <- function(x, posterior = NULL) {
  ## Validate pointwise log-likelihood input
  x <- as.matrix(x)

  if (!is.numeric(x) || length(dim(x)) != 2L) {
    stop("'x' must be a numeric matrix of pointwise log-likelihood values.",
         call. = FALSE)
  }
  if (nrow(x) < 2L) {
    stop("'x' must contain at least two posterior draws.", call. = FALSE)
  }
  if (ncol(x) < 1L) {
    stop("'x' must contain at least one observation.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("'x' must contain only finite, non-missing values.", call. = FALSE)
  }

  ## Basic sample properties are needed when validating posterior draws
  S <- nrow(x)
  n <- ncol(x)

  if (!is.null(posterior)) {
    posterior <- as.matrix(posterior)

    if (!is.numeric(posterior) || length(dim(posterior)) != 2L) {
      stop("'posterior' must be a numeric matrix.", call. = FALSE)
    }
    if (nrow(posterior) != S) {
      stop("'posterior' and 'x' must contain the same number of draws.",
           call. = FALSE)
    }
    if (ncol(posterior) < 1L) {
      stop("'posterior' must contain at least one parameter.", call. = FALSE)
    }
    if (anyNA(posterior) || any(!is.finite(posterior))) {
      stop("'posterior' must contain only finite, non-missing values.",
           call. = FALSE)
    }
  }

  ## Numerically stable log(mean(exp(z)))
  log_mean_exp <- function(z) {
    z_max <- max(z)
    z_max + log(mean(exp(z - z_max)))
  }

  ## Basic posterior quantities
  loglik_total   <- rowSums(x)
  deviance_draws <- -2 * loglik_total
  Dbar           <- mean(deviance_draws)

  mean_log_lik_i <- colMeans(x)
  var_log_lik_i  <- apply(x, 2L, stats::var)
  lppd_i         <- apply(x, 2L, log_mean_exp)
  lppd           <- sum(lppd_i)

  ## WBIC estimated by importance reweighting ordinary posterior draws
  wbic <- NA_real_
  wbic_beta <- NA_real_
  wbic_ess <- NA_real_
  wbic_relative_ess <- NA_real_
  wbic_max_weight <- NA_real_
  wbic_perplexity <- NA_real_
  wbic_relative_perplexity <- NA_real_

  if (n > 1L) {
    wbic_beta <- 1 / log(n)
    log_weights <- (wbic_beta - 1) * loglik_total
    log_weights <- log_weights - max(log_weights)
    weights <- exp(log_weights)
    weights <- weights / sum(weights)

    wbic <- -2 * sum(weights * loglik_total)
    wbic_ess <- 1 / sum(weights^2)
    wbic_relative_ess <- wbic_ess / S
    wbic_max_weight <- max(weights)

    positive_weights <- weights[weights > 0]
    wbic_entropy <- -sum(positive_weights * log(positive_weights))
    wbic_perplexity <- exp(wbic_entropy)
    wbic_relative_perplexity <- wbic_perplexity / S
  }

  ## WAIC
  p_waic1 <- 2 * sum(lppd_i - mean_log_lik_i)
  p_waic2 <- sum(var_log_lik_i)
  waic1 <- -2 * (lppd - p_waic1)
  waic2 <- -2 * (lppd - p_waic2)

  ## Conditional predictive ordinates and raw importance-sampling LOO
  log_cpo  <- -apply(-x, 2L, log_mean_exp)
  cpo      <- exp(log_cpo)
  lpml     <- sum(log_cpo)
  looic_is <- -2 * lpml
  p_loo_is <- lppd - lpml

  ## Optional PSIS-LOO
  loo_fit <- NULL
  looic <- NA_real_
  se_looic <- NA_real_
  elpd_loo <- NA_real_
  se_elpd_loo <- NA_real_
  p_loo <- NA_real_
  se_p_loo <- NA_real_
  pareto_k <- NULL

  if (requireNamespace("loo", quietly = TRUE)) {
    loo_attempt <- tryCatch(loo::loo(x, r_eff = 1),
                            error = function(e) e)

    if (inherits(loo_attempt, "error")) {
      warning("PSIS-LOO could not be computed: ",
              conditionMessage(loo_attempt),
              call. = FALSE)
    } else {
      loo_fit       <- loo_attempt
      loo_estimates <- loo_fit$estimates

      elpd_loo    <- unname(loo_estimates["elpd_loo", "Estimate"])
      se_elpd_loo <- unname(loo_estimates["elpd_loo", "SE"])
      p_loo       <- unname(loo_estimates["p_loo", "Estimate"])
      se_p_loo    <- unname(loo_estimates["p_loo", "SE"])
      looic       <- unname(loo_estimates["looic", "Estimate"])
      se_looic    <- unname(loo_estimates["looic", "SE"])
      pareto_k    <- loo_fit$diagnostics$pareto_k
    }
  }

  ## Parameterization-invariant DIC
  p_V <- 0.5 * var(deviance_draws)
  DIC_i <- Dbar + p_V

  ## Posterior correlational ICOMP analogues
  C_R <- NA_real_
  ICOMP <- NA_real_
  ICOMP_PEU <- NA_real_
  ICOMP_PEU_LN <- NA_real_
  CICOMP <- NA_real_

  if (!is.null(posterior)) {
    d <- ncol(posterior)
    posterior_sd <- apply(posterior, 2L, sd)

    if (any(!is.finite(posterior_sd)) || any(posterior_sd <= 0)) {
      warning("At least one posterior parameter has zero or non-finite variance. ",
              "Posterior ICOMP criteria were not computed. Remove fixed or ",
              "deterministic columns from 'posterior'.",
              call. = FALSE)
    } else {
      if (d == 1L) {
        ## The one-dimensional correlation matrix is [1].
        C_R <- 0
      } else {
        R_post <- cor(posterior)
        chol_R <- tryCatch(chol(R_post),
                           error = function(e) NULL)

        if (is.null(chol_R)) {
          ## For a singular correlation matrix, -log|R| diverges.
          C_R <- Inf

          warning("The posterior correlation matrix is singular or not positive ",
                  "definite. Posterior ICOMP criteria were set to Inf.",
                  call. = FALSE)
        } else {
          log_det_R <- 2 * sum(log(diag(chol_R)))
          C_R <- max(0, -0.5 * log_det_R)
        }
      }

      if (!is.na(C_R)) {
        ## Posterior analogue of correlational ICOMP(IFIM)
        ICOMP <- Dbar + 2 * C_R

        ## Posterior expected-utility analogue; preferred ICOMP result
        ICOMP_PEU <- Dbar + d + 2 * C_R

        ## Log-sample-size weighted PEU analogue
        ICOMP_PEU_LN <- Dbar + d + 2 * log(n) * C_R

        ## Consistent ICOMP analogue
        CICOMP <- Dbar + d * (log(n) + 1) + 2 * C_R
      }
    }
  }

  ## Final results
  IC <- list(Dbar = Dbar,
             
             WBIC = wbic,
             WBIC_beta = wbic_beta,
             WBIC_ESS = wbic_ess,
             WBIC_relative_ESS = wbic_relative_ess,
             WBIC_max_weight = wbic_max_weight,
             WBIC_perplexity = wbic_perplexity,
             WBIC_relative_perplexity = wbic_relative_perplexity,
             
             lppd = lppd,
             WAIC = waic2,
             WAIC1 = waic1,
             p_WAIC1 = p_waic1,
             WAIC2 = waic2,
             p_WAIC2 = p_waic2,
             
             CPO = cpo,
             log_CPO = log_cpo,
             LPML = lpml,
             LOOIC_IS = looic_is,
             p_LOO_IS = p_loo_is,
             
             LOOIC = looic,
             SE_LOOIC = se_looic,
             elpd_LOO = elpd_loo,
             SE_elpd_LOO = se_elpd_loo,
             p_LOO = p_loo,
             SE_p_LOO = se_p_loo,
             Pareto_k = pareto_k,
             PSIS_LOO = loo_fit,
             
             DIC_i = DIC_i,
             p_V = p_V,
             
             ICOMP = ICOMP,
             ICOMP_PEU = ICOMP_PEU,
             ICOMP_PEU_LN = ICOMP_PEU_LN,
             CICOMP = CICOMP,
             C_R = C_R)
  return(IC)
}
