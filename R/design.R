# Core model for manuscript v29. No external package is needed for inference.
check_probability <- function(x, name) {
  if (any(!is.finite(x) | x < 0 | x > 1)) stop(name, " must be in [0, 1].")
}

make_spec <- function(y = 60, n_h = 200, cap = 30, quality = 1,
                      p0 = .30, improvement = .15, delta = 0,
                      looks = c(20L, 40L, 60L), method = c("HBOP2-FE", "BOP2-FE"),
                      prior_t = c(.5, .5), prior_c = c(.5, .5)) {
  method <- match.arg(method)
  check_probability(c(quality, p0, delta, p0 + delta + improvement), "Design probabilities")
  stopifnot(length(y) == 1L, length(n_h) == 1L, is.finite(y), is.finite(n_h),
            n_h > 0, y >= 0, y <= n_h, is.finite(cap), cap >= 0,
            length(looks) >= 1L, looks[1] > 0, all(diff(looks) > 0),
            all(looks == as.integer(looks)), length(prior_t) == 2L,
            length(prior_c) == 2L, all(c(prior_t, prior_c) > 0), improvement > 0)
  # Compute effective size first so equivalent cap/discount pairs have identical
  # floating-point representations as well as identical mathematical meaning.
  effective <- quality * min(n_h, cap)
  w <- effective / n_h
  list(y = y, n_h = n_h, cap = cap, quality = quality, p0 = p0,
       improvement = improvement, delta = delta, looks = looks, method = method,
       prior_t = prior_t, prior_c = prior_c, w = w, ess = effective,
       a_c = prior_c[1] + w * y, b_c = prior_c[2] + w * (n_h - y))
}

# Numerically stable upper beta tail; quadrature errors are retained for auditing.
posterior_table <- function(spec, tolerance = 1e-10) {
  do.call(rbind, lapply(seq_along(spec$looks), function(k) {
    n <- spec$looks[k]
    rows <- lapply(0:n, function(x) {
      at <- spec$prior_t[1] + x; bt <- spec$prior_t[2] + n - x
      if (spec$method == "BOP2-FE") {
        value <- pbeta(spec$p0 + spec$delta, at, bt, lower.tail = FALSE)
        err <- 0
      } else if (spec$delta >= 1) {
        value <- 0; err <- 0
      } else {
        ans <- integrate(function(c) pbeta(c + spec$delta, at, bt, lower.tail = FALSE) *
                           dbeta(c, spec$a_c, spec$b_c),
                         0, 1 - spec$delta, rel.tol = tolerance, abs.tol = tolerance,
                         subdivisions = 1500L, stop.on.error = TRUE)
        value <- ans$value; err <- ans$abs.error
      }
      data.frame(look = k, n = n, x = x, pp = value, integration_error = err)
    })
    d <- do.call(rbind, rows)
    if (any(diff(d$pp) < -1e-9)) stop("Posterior probability is not monotone.")
    d
  }))
}

boundary_table <- function(spec, theta, posterior = posterior_table(spec)) {
  stopifnot(all(c("lambda", "gamma", "eta") %in% names(theta)),
            theta$lambda >= 0, theta$lambda <= 1, theta$gamma > 0,
            theta$gamma <= 1, theta$eta > 0, theta$eta <= 3)
  do.call(rbind, lapply(seq_along(spec$looks), function(k) {
    n <- spec$looks[k]; pp <- posterior$pp[posterior$look == k]
    cf <- theta$lambda * (n / max(spec$looks))^theta$gamma
    ce <- 1 - (1 - theta$lambda) * (n / max(spec$looks))^theta$eta
    if (k == length(spec$looks)) ce <- cf <- theta$lambda
    # -1 and n+1 mean an empty futility / efficacy region respectively.
    f <- max(c(-1L, which(pp < cf) - 1L))
    e <- min(c(n + 1L, which(pp >= ce) - 1L))
    stopifnot(f < e)
    data.frame(look = k, n = n, cutoff_futility = cf, cutoff_efficacy = ce,
               futility_max = f, efficacy_min = e)
  }))
}

# Exact forward probability propagation over surviving states. Works for any looks.
exact_oc <- function(p, boundaries, return_stages = FALSE) {
  check_probability(p, "p")
  stopifnot(length(p) == 1L)
  mass <- 1; previous_n <- 0; stages <- list()
  for (k in seq_len(nrow(boundaries))) {
    b <- boundaries[k, ]; inc <- b$n - previous_n
    next_mass <- numeric(b$n + 1L)
    for (x in which(mass > 0) - 1L) {
      index <- x + (0:inc) + 1L
      next_mass[index] <- next_mass[index] + mass[x + 1L] * dbinom(0:inc, inc, p)
    }
    x <- 0:b$n; final <- k == nrow(boundaries)
    fut <- if (final) 0 else sum(next_mass[x <= b$futility_max])
    eff <- if (final) 0 else sum(next_mass[x >= b$efficacy_min])
    fg <- if (final) sum(next_mass[x >= b$efficacy_min]) else 0
    fn <- if (final) sum(next_mass[x < b$efficacy_min]) else 0
    stages[[k]] <- data.frame(look = k, n = b$n, futility = fut, efficacy = eff,
                              final_go = fg, final_no_go = fn,
                              stopped = fut + eff + fg + fn)
    if (!final) next_mass[x <= b$futility_max | x >= b$efficacy_min] <- 0
    mass <- next_mass; previous_n <- b$n
  }
  d <- do.call(rbind, stages)
  if (abs(sum(d$stopped) - 1) > 1e-10) stop("Probability mass is not conserved.")
  out <- data.frame(go = sum(d$efficacy + d$final_go),
                    en = sum(d$n * d$stopped),
                    early_futility = sum(d$futility), early_efficacy = sum(d$efficacy),
                    final_go = sum(d$final_go), final_no_go = sum(d$final_no_go),
                    final_reached = tail(d$stopped, 1),
                    var_n = sum(d$n^2 * d$stopped) - sum(d$n * d$stopped)^2)
  if (return_stages) list(summary = out, stages = d) else out
}

# Independent simulation: patient Bernoulli outcomes and direct PP decisions,
# rather than exact-state propagation or integer-boundary lookups.
simulate_oc <- function(p, spec, theta, B = 100000L, seed = 20261107L,
                        posterior = posterior_table(spec), keep_records = FALSE) {
  check_probability(p, "p"); stopifnot(B > 1, B == as.integer(B))
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  status <- rep("active", B); final_n <- rep(max(spec$looks), B); x <- integer(B)
  previous_n <- 0
  for (k in seq_along(spec$looks)) {
    active <- which(status == "active"); n <- spec$looks[k]
    if (length(active)) {
      inc <- n - previous_n
      x[active] <- x[active] + rowSums(matrix(runif(length(active) * inc) < p,
                                              nrow = length(active), ncol = inc))
      pp <- posterior$pp[posterior$look == k][x[active] + 1L]
      if (k == length(spec$looks)) {
        status[active] <- ifelse(pp >= theta$lambda, "final_go", "final_no_go")
      } else {
        cf <- theta$lambda * (n / max(spec$looks))^theta$gamma
        ce <- 1 - (1 - theta$lambda) * (n / max(spec$looks))^theta$eta
        status[active[pp < cf]] <- "early_futility"
        status[active[pp >= ce]] <- "early_efficacy"
      }
      stopped <- active[status[active] != "active"]
      final_n[stopped] <- n
    }
    previous_n <- n
  }
  indicators <- list(go = status %in% c("early_efficacy", "final_go"), en = final_n,
                     early_futility = status == "early_futility",
                     early_efficacy = status == "early_efficacy",
                     final_go = status == "final_go", final_no_go = status == "final_no_go",
                     final_reached = final_n == max(spec$looks))
  out <- data.frame(metric = names(indicators),
                    estimate = vapply(indicators, mean, numeric(1)),
                    mcse = vapply(indicators, function(z) {
                      if (is.logical(z)) sqrt(mean(z) * (1 - mean(z)) / B) else sd(z) / sqrt(B)
                    }, numeric(1)), B = B, seed = seed)
  if (keep_records) list(summary = out, records = data.frame(status, final_n, x)) else out
}
