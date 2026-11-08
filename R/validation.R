# Independent exhaustive enumeration of all binomial-increment paths.
# This intentionally shares neither DP propagation nor rectangle summation.
enumerate_paths <- function(p, boundaries) {
  increments <- diff(c(0, boundaries$n))
  paths <- expand.grid(lapply(increments, function(n) 0:n))
  cumulative <- t(apply(as.matrix(paths), 1, cumsum))
  probability <- rep(1, nrow(paths))
  for (k in seq_along(increments)) probability <- probability * dbinom(paths[[k]], increments[k], p)
  active <- rep(TRUE, nrow(paths)); go <- ef <- ee <- fg <- numeric(nrow(paths))
  n_stop <- rep(tail(boundaries$n, 1), nrow(paths))
  for (k in seq_len(nrow(boundaries))) {
    final <- k == nrow(boundaries)
    success <- active & cumulative[, k] >= boundaries$efficacy_min[k]
    failure <- active & cumulative[, k] <= boundaries$futility_max[k]
    if (final) {
      fg[success] <- 1
    } else {
      ef[failure] <- 1; ee[success] <- 1
    }
    go[success] <- 1; n_stop[success | failure] <- boundaries$n[k]
    active[success | failure] <- FALSE
  }
  stopifnot(!any(active))
  c(go = sum(go * probability), en = sum(n_stop * probability),
    early_futility = sum(ef * probability), early_efficacy = sum(ee * probability),
    final_go = sum(fg * probability))
}

compare_mc <- function(mc, exact) {
  d <- mc
  d$exact <- as.numeric(exact[1, d$metric])
  d$difference <- d$estimate - d$exact
  # The reference SE uses exact moments, avoiding a zero estimated SE at endpoints.
  d$reference_se <- sqrt(pmax(0, d$exact * (1 - d$exact)) / d$B)
  d$reference_se[d$metric == "en"] <- sqrt(exact$var_n / d$B[1])
  d$z <- ifelse(d$reference_se > 0, d$difference / d$reference_se, 0)
  d
}

finalize_mc_validation <- function(d, family_alpha = .01) {
  # Simultaneous 99% family coverage via Bonferroni. Probability intervals are
  # exact Clopper-Pearson; EN uses a normal interval with known exact variance.
  m <- nrow(d); alpha_each <- family_alpha / m
  d$lower <- d$upper <- NA_real_
  for (i in seq_len(m)) {
    if (d$metric[i] == "en") {
      half <- qnorm(1 - alpha_each / 2) * d$reference_se[i]
      d$lower[i] <- d$estimate[i] - half; d$upper[i] <- d$estimate[i] + half
    } else {
      successes <- as.integer(round(d$estimate[i] * d$B[i]))
      ci <- binom.test(successes, d$B[i], conf.level = 1 - alpha_each)$conf.int
      d$lower[i] <- ci[1]; d$upper[i] <- ci[2]
    }
  }
  d$covered <- d$exact >= d$lower - 1e-12 & d$exact <= d$upper + 1e-12
  d$family_alpha <- family_alpha
  d
}

validate_designs <- function(designs) {
  checks <- list(); add <- function(name, difference, tolerance) {
    checks[[length(checks) + 1L]] <<- data.frame(check = name, difference = difference,
                        tolerance = tolerance, passed = is.finite(difference) && difference <= tolerance)
  }
  for (id in names(designs)) {
    z <- designs[[id]]
    for (p in c(0, .27, .30, .45, .47, 1)) {
      dp <- exact_oc(p, z$boundaries)
      paths <- enumerate_paths(p, z$boundaries)
      add(paste(id, "DP vs all paths", p), max(abs(unlist(dp[names(paths)]) - paths)), 1e-11)
    }
    pp <- z$posterior
    add(paste(id, "integration absolute error"), max(pp$integration_error), 1e-8)
    # Reverse integral over treatment rate is an independent form of Pr(T>C+d).
    if (z$spec$method == "HBOP2-FE") {
      for (k in seq_along(z$spec$looks)) {
        n <- z$spec$looks[k]
        for (x in c(0, floor(n * .3), n)) {
          at <- z$spec$prior_t[1] + x; bt <- z$spec$prior_t[2] + n - x
          reverse <- integrate(function(t) dbeta(t, at, bt) *
                    pbeta(pmax(0, t - z$spec$delta), z$spec$a_c, z$spec$b_c),
                    z$spec$delta, 1, rel.tol = 1e-10, abs.tol = 1e-10,
                    subdivisions = 1500L)$value
          forward <- pp$pp[pp$n == n & pp$x == x]
          add(paste(id, "forward vs reverse integral", n, x), abs(forward - reverse), 1e-8)
        }
      }
    }
  }
  do.call(rbind, checks)
}
