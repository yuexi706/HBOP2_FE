#!/usr/bin/env Rscript
source("R/design.R"); source("R/calibration.R"); source("R/validation.R"); source("config/design.R")
tests <- 0L
verify <- function(condition) { stopifnot(condition); tests <<- tests + 1L }
s <- make_spec(); pp <- posterior_table(s)
verify(abs(s$w - .15) < 1e-14 && s$ess == 30 && s$a_c == 9.5 && s$b_c == 21.5)
verify(make_spec(cap = 500)$w == 1)
verify(make_spec(cap = 0)$ess == 0 && make_spec(quality = 0)$ess == 0)
verify(inherits(try(make_spec(y = 201), silent = TRUE), "try-error"))
verify(inherits(try(make_spec(quality = 1.1), silent = TRUE), "try-error"))
verify(inherits(try(make_spec(looks = c(20, 20, 60)), silent = TRUE), "try-error"))
verify(all(vapply(split(pp$pp, pp$n), function(v) all(diff(v) >= -1e-10), logical(1))))
theta <- list(lambda = .77, gamma = .8, eta = 3)
b <- boundary_table(s, theta, pp)
for (p in c(0, .3, .45, 1)) {
  d <- exact_oc(p, b); e <- enumerate_paths(p, b)
  verify(max(abs(unlist(d[names(e)]) - e)) < 1e-11)
  verify(abs(d$go - d$early_efficacy - d$final_go) < 1e-11)
  verify(abs(d$early_futility + d$early_efficacy + d$final_go + d$final_no_go - 1) < 1e-11)
}
# Ties must continue at CF and succeed at CE / final lambda.
tied <- pp; ix <- tied$look == 3 & tied$x == 23
theta_tie <- theta; theta_tie$lambda <- tied$pp[ix]
bt <- boundary_table(s, theta_tie, tied)
verify(bt$efficacy_min[3] == 23 && bt$futility_max[3] == 22)
# Zero borrowing leaves a vague random comparator, which differs from a fixed null.
zero <- posterior_table(make_spec(cap = 0)); fixed <- posterior_table(make_spec(method = "BOP2-FE"))
verify(max(abs(zero$pp - fixed$pp)) > .1)
# Equal effective borrowing means identical statistical designs.
s1 <- make_spec(cap = 20, quality = .75); s2 <- make_spec(cap = 30, quality = .5)
verify(identical(posterior_table(s1), posterior_table(s2)))
z <- calibrate(s, settings$grid)
verify(nrow(z$grid) == 13000 && z$selected$type1 <= .1)
verify(z$selected$power >= max(z$grid$power[z$grid$feasible]) - 1e-12)
verify(max(abs(exact_oc(.30, z$boundaries)$go - z$selected$type1)) < 1e-12)
verify(inherits(try(calibrate(s, settings$grid[1, , drop = FALSE], alpha = 0), silent = TRUE), "try-error"))
# A nonzero margin and a non-default schedule exercise the general engine.
sm <- make_spec(delta = .05, looks = c(8, 16, 24))
bm <- boundary_table(sm, theta)
verify(max(abs(unlist(exact_oc(.35, bm)[names(enumerate_paths(.35, bm))]) - enumerate_paths(.35, bm))) < 1e-11)
ks <- make_spec(172, 637, p0 = .27, improvement = .2)
kb <- boundary_table(ks, settings$case_theta)
verify(identical(kb$futility_max, c(4L, 11L, 20L)) && identical(kb$efficacy_min, c(12L, 18L, 21L)))
verify(round(exact_oc(.27, kb)$go, 3) == .099 && round(exact_oc(.47, kb)$go, 3) == .966)
r1 <- simulate_oc(.3, s, theta, B = 1000, seed = 27, posterior = pp)
r2 <- simulate_oc(.3, s, theta, B = 1000, seed = 27, posterior = pp)
verify(identical(r1, r2))
cat("PASS:", tests, "statistical assertions.\n")
