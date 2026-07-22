VB <- function(Model, Data, Initial.Values = NULL, iterations = NULL,
               status = NULL, thinning = NULL, adapt = NULL,
               method = c("sagva", "qnsagva", "mcvi", "svgd"),
               n.particles = 50, par.cov = NULL,
               # Retained for backward compatibility with the preliminary SAGVA API.
               # The corrected SAGVA algorithm does not use a line search.
               w.min = 1e-8, w.max = 1.0,
               armijo.c1 = 1e-4, wolfe.c2 = 0.9,
               ls.max.iter = 20,
               # Shared numerical differentiation and convergence controls
               h = 1e-4, stop.tolerance = NULL,
               learning.rate = NULL, min.iter = NULL, patience = 50,
               # QN-SAGVA
               memory = 10, damping = 0.2,
               # MCVI / Hamiltonian VI
               T.mcmc = 3, L.leapfrog = 3, k.rb = 5,
               mcvi.spsa.scale = 1e-3, mcvi.epsilon = NULL,
               # SVGD
               step.size = 0.01, h.scale = 1.0, use.adam = TRUE,
               # Initialization and numerical recovery
               opt.init = TRUE, n.starts = 3, init.jitter = 0.5,
               particle.scale = 1.5, retry.on.failure = TRUE,
               control = list(), seed = NULL) {

  method <- match.arg(method)
  if (!is.function(Model)) stop("'Model' must be a function.")
  if (!is.list(Data)) stop("'Data' must be a list.")
  if (!is.null(seed)) set.seed(seed)
  if (is.null(stop.tolerance)) {
    stop.tolerance <- switch(method, sagva   = 5e-4,
                                     qnsagva = 5e-4,
                                     mcvi    = 1e-3,
                                     svgd    = 5e-3)
  }
  if (is.null(Initial.Values)) {
    if (!is.function(Data$PGF)) {
      stop("Supply 'Initial.Values' or a Data$PGF function.")
    }
    Initial.Values <- Data$PGF(Data)
  }
  Initial.Values <- as.numeric(Initial.Values)
  if (!length(Initial.Values) || any(!is.finite(Initial.Values))) {
    stop("'Initial.Values' must be a finite numeric vector.")
  }

  iterations <- if (is.null(iterations)) 1000L else as.integer(iterations)
  if (!is.finite(iterations) || iterations < 1L) iterations <- 1000L
  thinning <- if (is.null(thinning)) 1L else as.integer(thinning)
  if (!is.finite(thinning) || thinning < 1L) thinning <- 1L
  status <- if (is.null(status)) max(1L, round(iterations / 10)) else as.integer(status)
  if (!is.finite(status) || status < 1L) status <- max(1L, round(iterations / 10))
  adapt <- if (is.null(adapt)) floor(iterations / 2) else as.integer(adapt)
  adapt <- min(max(adapt, 0L), iterations)
  min.iter <- if (is.null(min.iter)) max(50L, adapt) else as.integer(min.iter)
  min.iter <- min(max(min.iter, 1L), iterations)
  patience <- min(max(as.integer(patience), 1L), max(1L, iterations - min.iter + 1L))
  n.particles <- max(as.integer(n.particles), 2L)
  n.starts <- max(as.integer(n.starts), 1L)
  memory <- max(as.integer(memory), 1L)
  T.mcmc <- max(as.integer(T.mcmc), 1L)
  L.leapfrog <- max(as.integer(L.leapfrog), 1L)
  k.rb <- max(as.integer(k.rb), 1L)

  logpost <- function(par) {
    out <- tryCatch(Model(as.numeric(par), Data), error = function(e) NULL)
    if (is.null(out) || is.null(out$LP) || length(out$LP) != 1L ||
        !is.finite(out$LP)) return(-Inf)
    as.numeric(out$LP)
  }
  objective <- function(par) {
    lp <- logpost(par)
    if (!is.finite(lp)) .Machine$double.xmax / 100 else -lp
  }

  stabilize_cov <- function(S, d, eig.floor = 1e-8, eig.ceiling = 1e8) {
    S <- tryCatch(as.matrix(S), error = function(e) NULL)
    if (is.null(S) || !identical(dim(S), c(d, d)) || any(!is.finite(S))) {
      return(diag(0.1^2, d))
    }
    S <- (S + t(S)) / 2
    ee <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
    if (is.null(ee)) return(diag(0.1^2, d))
    values <- pmin(pmax(ee$values, eig.floor), eig.ceiling)
    out <- ee$vectors %*% (values * t(ee$vectors))
    (out + t(out)) / 2
  }

  d <- length(Initial.Values)
  if (!is.null(par.cov)) {
    if (!is.matrix(par.cov) || !identical(dim(par.cov), c(d, d))) {
      stop("'par.cov' must be a square covariance matrix matching Initial.Values.")
    }
    par.cov <- stabilize_cov(par.cov, d)
  }

  startTime <- proc.time()
  map.info <- NULL

  if (isTRUE(opt.init)) {
    cat("Improving initial values with multi-start MAP estimation.\n")
    opt.control <- modifyList(list(maxit = 1000L, reltol = 1e-10), control)

    if (!is.null(par.cov)) {
      jitter.sd <- sqrt(pmax(diag(par.cov), 1e-10)) * init.jitter
    } else {
      jitter.sd <- pmax(abs(Initial.Values), 1) * init.jitter
    }

    starts <- vector("list", n.starts)
    starts[[1L]] <- Initial.Values
    if (n.starts > 1L) {
      for (i in 2:n.starts) starts[[i]] <- Initial.Values + rnorm(d, 0, jitter.sd)
    }

    fits <- lapply(starts, function(st) {
      fit <- tryCatch(
        optim(st, objective, method = "BFGS", hessian = TRUE, control = opt.control),
        error = function(e) NULL
      )
      if (is.null(fit) || !is.finite(fit$value)) {
        nm <- tryCatch(
          optim(st, objective, method = "Nelder-Mead", control = opt.control),
          error = function(e) NULL
        )
        if (!is.null(nm) && is.finite(nm$value)) {
          fit <- tryCatch(
            optim(nm$par, objective, method = "BFGS", hessian = TRUE,
                  control = opt.control),
            error = function(e) nm
          )
        }
      }
      fit
    })

    ok <- vapply(fits, function(x) !is.null(x) && is.finite(x$value) &&
                   all(is.finite(x$par)), logical(1))
    if (any(ok)) {
      candidates <- fits[ok]
      best <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "value"))]]
      Initial.Values <- as.numeric(best$par)
      map.info <- list(par = Initial.Values, LP = -best$value,
                       convergence = best$convergence,
                       n.starts = n.starts)
      if (identical(best$convergence, 0L)) {
        cat("MAP optimization converged successfully.\n")
      } else {
        cat("The best MAP run did not fully converge; its finite solution will be used.\n")
      }
    } else {
      warning("All MAP starts failed; using the supplied initial values.")
    }
  }

  # A Laplace covariance is used only as an initialization.
  if (is.null(par.cov)) {
    H <- tryCatch(optimHess(Initial.Values, objective), error = function(e) NULL)
    if (!is.null(H) && identical(dim(H), c(d, d)) && all(is.finite(H))) {
      H <- (H + t(H)) / 2
      eh <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
      if (!is.null(eh)) {
        precision.values <- pmin(pmax(eh$values, 1e-6), 1e8)
        par.cov <- eh$vectors %*% ((1 / precision.values) * t(eh$vectors))
      }
    }
    if (is.null(par.cov)) par.cov <- diag(0.1^2, d)
  }
  par.cov <- stabilize_cov(par.cov, d)

  MO0 <- tryCatch(Model(Initial.Values, Data), error = function(e) NULL)
  required <- c("LP", "parm", "yhat", "Dev", "Monitor")
  if (is.null(MO0) || !all(required %in% names(MO0))) {
    stop("Model(par, Data) must return LP, parm, yhat, Dev, and Monitor.")
  }

  n.draws <- max(1L, floor(iterations / thinning))
  make_storage <- function() list(
    thinned = matrix(as.numeric(MO0$parm), n.draws, length(MO0$parm), byrow = TRUE),
    postpred = matrix(as.numeric(MO0$yhat), n.draws, length(MO0$yhat), byrow = TRUE),
    Dev = matrix(as.numeric(MO0$Dev), n.draws, length(MO0$Dev), byrow = TRUE),
    Mon = matrix(as.numeric(MO0$Monitor), n.draws, length(MO0$Monitor), byrow = TRUE)
  )

  run_once <- function(covariance = par.cov, retry = FALSE) {
    st <- make_storage()
    lr <- if (is.null(learning.rate)) -1 else as.numeric(learning.rate)

    switch(method,
      sagva = {
        cat("Algorithm: stochastic approximation for Gaussian variational approximation\n\n")
        sagva(Model = Model, Data = Data, Iterations = iterations, Status = status,
              InitialValues = Initial.Values, InitialCov = covariance,
              Thinning = thinning, thinned = st$thinned, postpred = st$postpred,
              Dev = st$Dev, Mon = st$Mon, h = h, learning_rate = lr,
              Stop_Tolerance = stop.tolerance, Min_Iterations = min.iter,
              Patience = patience)
      },
      qnsagva = {
        cat("Algorithm: limited-memory quasi-Newton SAGVA\n\n")
        qnsagva(Model = Model, Data = Data, Iterations = iterations, Status = status,
                InitialValues = Initial.Values, InitialCov = covariance,
                Thinning = thinning, thinned = st$thinned, postpred = st$postpred,
                Dev = st$Dev, Mon = st$Mon, memory = memory, h = h,
                damping = damping, learning_rate = lr,
                Stop_Tolerance = stop.tolerance, Min_Iterations = min.iter,
                Patience = patience)
      },
      mcvi = {
        cat("Algorithm: Hamiltonian Markov-chain variational inference\n\n")
        eps0 <- if (is.null(mcvi.epsilon)) -1 else as.numeric(mcvi.epsilon)
        if (retry && eps0 > 0) eps0 <- eps0 / 2
        mcvi(Model = Model, Data = Data, Iterations = iterations, Status = status,
             thinned = st$thinned, postpred = st$postpred, Dev = st$Dev,
             Mon = st$Mon, InitialValues = Initial.Values,
             InitialCov = covariance, Thinning = thinning,
             T_mcmc = T.mcmc, L_leapfrog = L.leapfrog, K_rb = k.rb,
             h = h, learning_rate = if (lr > 0) lr else 5e-3,
             spsa_scale = mcvi.spsa.scale, initial_epsilon = eps0,
             Stop_Tolerance = stop.tolerance, Min_Iterations = min.iter,
             Patience = patience)
      },
      svgd = {
        cat("Algorithm: Stein variational gradient descent\n\n")
        n.core <- max(1L, ceiling(0.8 * n.particles))
        n.wide <- n.particles - n.core
        scale.now <- if (retry) particle.scale / 2 else particle.scale
        core <- MASS::mvrnorm(n.core, Initial.Values,
                              Sigma = covariance * scale.now^2)
        wide <- if (n.wide > 0L) {
          MASS::mvrnorm(n.wide, Initial.Values,
                        Sigma = covariance * (3 * scale.now)^2)
        } else matrix(numeric(0), 0L, d)
        particles <- rbind(core, wide)
        particles[1L, ] <- Initial.Values
        svgd(Model = Model, Data = Data, Iterations = iterations, Status = status,
             InitialParticles = particles, Thinning = thinning,
             thinned = st$thinned, postpred = st$postpred,
             Dev = st$Dev, Mon = st$Mon,
             step_size = if (retry) step.size / 2 else step.size,
             h_scale = h.scale, use_adam = use.adam, h_grad = h,
             Stop_Tolerance = stop.tolerance, Min_Iterations = min.iter,
             Patience = patience)
      }
    )
  }

  cat("Variational inference method:", method, "\n\n")
  fit <- run_once()

  numerical_failure <- function(x) {
    if (is.null(x) || is.null(x$n_iter) || !is.finite(x$n_iter)) return(TRUE)
    if (method %in% c("sagva", "qnsagva", "mcvi")) {
      if (is.null(x$mean) || any(!is.finite(x$mean)) ||
          is.null(x$VarCov) || any(!is.finite(x$VarCov))) return(TRUE)
    }
    if (method == "svgd" && (is.null(x$particles) || any(!is.finite(x$particles)))) {
      return(TRUE)
    }
    !is.null(x$criterion) && !is.finite(x$criterion)
  }

  if (isTRUE(retry.on.failure) && numerical_failure(fit)) {
    warning("The first VB run was numerically unstable; retrying with a tighter initialization.")
    fit2 <- run_once(stabilize_cov(par.cov * 0.25, d), retry = TRUE)
    if (!numerical_failure(fit2)) fit <- fit2
  }

  if (method == "svgd") {
    posterior <- as.matrix(fit$particles)
    latent <- as.matrix(fit$latent_particles)
    ppred <- as.matrix(fit$particles_yhat)
    Monitor <- as.matrix(fit$particles_Mon)
    dev <- as.matrix(fit$particles_Dev)
    LP <- as.numeric(fit$particles_lp)
  } else {
    posterior <- as.matrix(fit$thinned)
    latent <- as.matrix(fit$latent)
    ppred <- as.matrix(fit$postpred)
    Monitor <- as.matrix(fit$Mon)
    dev <- as.matrix(fit$Dev)
    LP <- as.numeric(fit$LP)
  }

  pnames <- Data[["parm.names"]]
  if (is.null(pnames) || length(pnames) != ncol(posterior)) {
    pnames <- paste0("parm", seq_len(ncol(posterior)))
  }
  colnames(posterior) <- pnames
  if (ncol(latent) == length(Initial.Values)) {
    colnames(latent) <- if (length(Data[["latent.parm.names"]]) == ncol(latent))
      Data[["latent.parm.names"]] else paste0("latent", seq_len(ncol(latent)))
  }

  dev.vec <- if (ncol(dev) == 1L) as.numeric(dev) else rowSums(dev)
  Dbar <- mean(dev.vec, na.rm = TRUE)
  pV <- stats::var(dev.vec, na.rm = TRUE) / 2
  DIC <- list(DIC = Dbar + pV, Dbar = Dbar, pD = pV,
              definition = "variance-based effective number of parameters")

  elbo <- if (method == "svgd") NA_real_ else as.numeric(fit$LowerBound)
  ess <- if (method == "svgd") rep(NA_real_, ncol(posterior)) else
    rep(nrow(posterior), ncol(posterior))
  names(ess) <- pnames

  stopTime <- proc.time()
  elapsedTime <- unname((stopTime - startTime)[3])

  diagnostics <- list(
    converged = isTRUE(fit$converged),
    criterion = if (!is.null(fit$criterion)) as.numeric(fit$criterion) else NA_real_,
    optimization.iterations = if (!is.null(fit$n_iter)) as.integer(fit$n_iter) else iterations,
    invalid.updates = if (!is.null(fit$invalid_updates)) fit$invalid_updates else
      if (!is.null(fit$invalid_trajectories)) fit$invalid_trajectories else
        if (!is.null(fit$invalid_gradients)) fit$invalid_gradients else 0L
  )

  extras <- list(
    latent.posterior = data.frame(latent),
    variational.mean = if (!is.null(fit$mean)) fit$mean else NULL,
    variational.VarCov = if (!is.null(fit$VarCov)) fit$VarCov else NULL,
    MAP = map.info,
    diagnostics = diagnostics
  )

  if (method == "qnsagva") {
    extras$accepted.curvature.pairs <- fit$accepted_curvature_pairs
  }
  if (method == "mcvi") {
    extras$q0.mean <- fit$q0_mean
    extras$q0.VarCov <- fit$q0_VarCov
    extras$inverse.model <- list(A = fit$A_r, C = fit$C_r,
                                 b = fit$b_r, log.sd = fit$log_sd_r)
    extras$hmc.step.size <- fit$epsilon
    extras$ELBO.history <- fit$ELBO_history
  }
  if (method == "svgd") {
    extras$final.particles <- fit$particles
    extras$latent.particles <- fit$latent_particles
    extras$particles.lp <- fit$particles_lp
    extras$kernel.bandwidth <- fit$bandwidth
    extras$Stein.update.history <- fit$Stein_update_history
  }

  vb.info <- list(method = method,
                  optimization.iterations = diagnostics$optimization.iterations,
                  posterior.draws = nrow(posterior),
                  thinning = thinning,
                  min.iter = min.iter,
                  patience = patience,
                  n.particles = if (method == "svgd") n.particles else NA_integer_,
                  ELBO = elbo,
                  elapsed.mins = elapsedTime / 60)

  Result <- c(list(posterior = data.frame(posterior),
                   yhat = ppred,
                   LP = LP,
                   Monitor = if (ncol(Monitor) == 1L) as.numeric(Monitor) else Monitor,
                   DIC = DIC,
                   ESS = ess,
                   PSRF = rep(NA_real_, ncol(posterior)),
                   MPSRF = NA_real_,
                   vb.info = vb.info), extras)

  class(Result) <- "YABS_VB"
  cat("\nDone! Elapsed time:", round(elapsedTime, 2), "seconds.\n")
  Result
}
