# Vectorized exact calibration for three looks. A joint distribution of cumulative
# X1 and X2 plus the last-increment survival function gives the full stopping rule.
prefix2 <- function(a) {
  z <- apply(a, 2, cumsum)
  t(apply(z, 1, cumsum))
}

rectangle_mass <- function(prefix, lo1, hi1, lo2, hi2) {
  value <- function(i, j) {
    valid <- i >= 0 & j >= 0
    out <- numeric(length(i))
    out[valid] <- prefix[cbind(pmin(i[valid], nrow(prefix) - 1) + 1L,
                               pmin(j[valid], ncol(prefix) - 1) + 1L)]
    out
  }
  ans <- value(hi1, hi2) - value(lo1 - 1, hi2) - value(hi1, lo2 - 1) + value(lo1 - 1, lo2 - 1)
  ans[hi1 < lo1 | hi2 < lo2] <- 0
  pmax(0, ans)
}

grid_policies <- function(spec, grid, posterior) {
  stopifnot(length(spec$looks) == 3L)
  out <- grid
  for (k in 1:3) {
    n <- spec$looks[k]; pp <- posterior$pp[posterior$look == k]
    cf <- grid$lambda * (n / max(spec$looks))^grid$gamma
    ce <- 1 - (1 - grid$lambda) * (n / max(spec$looks))^grid$eta
    if (k == 3L) ce <- cf <- grid$lambda
    # Futility is strict (<), efficacy is inclusive (>=), including exact ties.
    out[[paste0("f", k)]] <- vapply(cf, function(t) sum(pp < t) - 1L, integer(1))
    out[[paste0("e", k)]] <- vapply(ce, function(t) sum(pp < t), integer(1))
  }
  out
}

exact_grid_oc <- function(p, looks, policies) {
  stopifnot(length(looks) == 3L)
  a <- 0:looks[1]; b <- 0:looks[2]
  joint <- outer(a, b, function(x1, x2) dbinom(x1, looks[1], p) *
                   dbinom(x2 - x1, looks[2] - looks[1], p))
  pre <- prefix2(joint)
  d <- policies; lo <- d$f1 + 1L; hi <- d$e1 - 1L
  c1 <- pbinom(hi, looks[1], p) - pbinom(d$f1, looks[1], p)
  early1 <- pbinom(d$e1 - 1, looks[1], p, lower.tail = FALSE)
  early2 <- rectangle_mass(pre, lo, hi, d$e2, rep(looks[2], nrow(d)))
  c2 <- rectangle_mass(pre, lo, hi, d$f2 + 1L, d$e2 - 1L)
  fut2 <- rectangle_mass(pre, lo, hi, rep(0, nrow(d)), d$f2)
  fg <- numeric(nrow(d))
  for (e in unique(d$e3)) {
    ii <- which(d$e3 == e)
    weighted <- sweep(joint, 2, pbinom(e - b - 1, looks[3] - looks[2], p,
                                     lower.tail = FALSE), "*")
    fg[ii] <- rectangle_mass(prefix2(weighted), lo[ii], hi[ii], d$f2[ii] + 1L, d$e2[ii] - 1L)
  }
  data.frame(go = early1 + early2 + fg,
             en = looks[1] + diff(looks)[1] * c1 + diff(looks)[2] * c2,
             early_futility = pbinom(d$f1, looks[1], p) + fut2,
             early_efficacy = early1 + early2, final_go = fg)
}

calibrate <- function(spec, grid, alpha = .10) {
  pp <- posterior_table(spec)
  policies <- grid_policies(spec, grid, pp)
  null <- exact_grid_oc(spec$p0 + spec$delta, spec$looks, policies)
  alt <- exact_grid_oc(spec$p0 + spec$delta + spec$improvement, spec$looks, policies)
  d <- cbind(policies, type1 = null$go, power = alt$go, en_null = null$en,
             null_early_futility = null$early_futility,
             null_early_efficacy = null$early_efficacy, en_alt = alt$en)
  d$feasible <- d$type1 <= alpha
  if (!any(d$feasible)) stop("No feasible design; do not substitute an infeasible candidate.")
  # Stabilize roundoff only at the 12th decimal for mathematically tied policies.
  feasible <- which(d$feasible)
  ord <- with(d[feasible, ], order(-round(power, 12), round(en_null, 12),
                   -round(null_early_futility, 12), round(null_early_efficacy, 12),
                   lambda, gamma, eta))
  chosen <- feasible[ord[1]]
  selected <- d[chosen, , drop = FALSE]
  selected$equivalent_parameter_sets <- sum(apply(d[, paste0(rep(c("f", "e"), 3), rep(1:3, each = 2))],
    1, function(z) all(z == unlist(selected[, paste0(rep(c("f", "e"), 3), rep(1:3, each = 2))]))))
  theta <- as.list(selected[, c("lambda", "gamma", "eta")])
  list(grid = d, selected = selected, theta = theta, posterior = pp,
       boundaries = boundary_table(spec, theta, pp), spec = spec)
}
