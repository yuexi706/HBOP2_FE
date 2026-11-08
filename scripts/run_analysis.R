#!/usr/bin/env Rscript
# Run from project root: Rscript --vanilla scripts/run_analysis.R
source("config/design.R")
for (f in c("design", "calibration", "validation")) source(file.path("R", paste0(f, ".R")))
args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
B <- if (quick) 2000L else settings$B
for (d in c("tables", "calibration", "validation", "figures"))
  dir.create(file.path("results", d), recursive = TRUE, showWarnings = FALSE)
write_csv <- function(d, file) write.csv(d, file.path("results", file), row.names = FALSE, na = "")

designs <- list(); cache <- new.env(parent = emptyenv())
register <- function(id, spec) {
  key <- paste(spec$method, spec$y, spec$n_h, spec$ess, spec$p0, spec$improvement, sep = "_")
  if (exists(key, cache, inherits = FALSE)) z <- get(key, cache) else {
    z <- calibrate(spec, settings$grid, settings$alpha); assign(key, z, cache)
  }
  z$spec <- spec
  designs[[id]] <<- z
  write_csv(z$grid, paste0("calibration/", id, "_grid.csv"))
  message(sprintf("%-24s lambda=%.2f gamma=%.2f eta=%.2f alpha=%.6f power=%.6f",
           id, z$theta$lambda, z$theta$gamma, z$theta$eta, z$selected$type1, z$selected$power))
}

register("primary_hb", make_spec())
register("primary_fixed", make_spec(method = "BOP2-FE"))
for (i in seq_len(nrow(settings$sensitivity))) {
  b <- settings$sensitivity[i, ]
  register(sprintf("sensitivity_%02d", i), make_spec(cap = b$cap, quality = b$quality))
}
register("case_hb_optimized", make_spec(172, 637, p0 = .27, improvement = .20))
register("case_fixed", make_spec(172, 637, p0 = .27, improvement = .20, method = "BOP2-FE"))
for (cap in c(10, 50)) register(paste0("case_cap", cap),
                   make_spec(172, 637, cap = cap, p0 = .27, improvement = .20))
register("case_discounted", make_spec(172, 637, quality = .5, p0 = .27, improvement = .20))
case_spec <- make_spec(172, 637, p0 = .27, improvement = .20)
case_pp <- posterior_table(case_spec)
designs$case_published <- list(spec = case_spec, posterior = case_pp,
   theta = settings$case_theta, boundaries = boundary_table(case_spec, settings$case_theta, case_pp),
   selected = NULL)

design_metadata <- function(id, z) {
  s <- z$spec
  data.frame(design = id, method = s$method, y_h = s$y, n_h = s$n_h,
    eca_rate = s$y / s$n_h, p0 = s$p0, delta = s$delta, improvement = s$improvement,
    cap = if (s$method == "HBOP2-FE") s$cap else NA_real_,
    quality = if (s$method == "HBOP2-FE") s$quality else NA_real_,
    w = if (s$method == "HBOP2-FE") s$w else NA_real_,
    ess = if (s$method == "HBOP2-FE") s$ess else NA_real_,
    pmd = if (s$method == "HBOP2-FE") s$a_c / (s$a_c + s$b_c) - s$p0 else NA_real_,
    lambda = z$theta$lambda, gamma = z$theta$gamma, eta = z$theta$eta)
}
metadata <- do.call(rbind, lapply(names(designs), function(id) design_metadata(id, designs[[id]])))
write_csv(metadata, "tables/design_settings.csv")
write_csv(primary_scenarios, "tables/primary_scenarios.csv")
write_csv(case_scenarios, "tables/case_scenarios.csv")
write_csv(settings$sensitivity, "tables/borrowing_configurations.csv")

calibration_summary <- do.call(rbind, lapply(names(designs), function(id) {
  z <- designs[[id]]
  cbind(design_metadata(id, z), exact_type1 = exact_oc(z$spec$p0, z$boundaries)$go,
    exact_power = exact_oc(z$spec$p0 + z$spec$improvement, z$boundaries)$go,
    en_null = exact_oc(z$spec$p0, z$boundaries)$en,
    en_alt = exact_oc(z$spec$p0 + z$spec$improvement, z$boundaries)$en,
    selection = if (id == "case_published") "Manuscript supplied" else "Exact grid optimum",
    grid_candidates = if (id == "case_published") NA_integer_ else nrow(z$grid),
    feasible_candidates = if (id == "case_published") NA_integer_ else sum(z$grid$feasible),
    equivalent_parameter_sets = if (id == "case_published") NA_integer_ else z$selected$equivalent_parameter_sets)
}))
write_csv(calibration_summary, "tables/calibration_summary.csv")
boundary_output <- do.call(rbind, lapply(names(designs), function(id)
                 cbind(design = id, designs[[id]]$boundaries)))
write_csv(boundary_output, "tables/decision_boundaries.csv")
posterior_output <- do.call(rbind, lapply(names(designs), function(id)
                 cbind(design = id, designs[[id]]$posterior)))
write_csv(posterior_output, "validation/posterior_lookup.csv")

exact_rows <- mc_rows <- stopping_rows <- list(); seed_index <- 0L
for (id in names(designs)) {
  z <- designs[[id]]; scenarios <- if (startsWith(id, "case_")) case_scenarios else primary_scenarios
  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, ]; seed_index <- seed_index + 1L
    exact <- exact_oc(sc$p_t, z$boundaries, return_stages = TRUE)
    exact_rows[[length(exact_rows) + 1L]] <- cbind(design_metadata(id, z), sc, exact$summary)
    stopping_rows[[length(stopping_rows) + 1L]] <- cbind(design = id, scenario = sc$scenario, exact$stages)
    mc <- simulate_oc(sc$p_t, z$spec, z$theta, B = B,
                     seed = settings$seed + seed_index, posterior = z$posterior)
    mc_rows[[length(mc_rows) + 1L]] <- cbind(design = id, scenario = sc$scenario,
                   compare_mc(mc, exact$summary))
  }
  message("Exact + Monte Carlo complete: ", id)
}
exact_results <- do.call(rbind, exact_rows)
mc_results <- finalize_mc_validation(do.call(rbind, mc_rows))
write_csv(exact_results, "tables/all_operating_characteristics.csv")
write_csv(subset(exact_results, design %in% c("primary_hb", "primary_fixed")), "tables/primary_operating_characteristics.csv")
write_csv(subset(exact_results, startsWith(design, "sensitivity_")), "tables/borrowing_sensitivity.csv")
write_csv(subset(exact_results, startsWith(design, "case_")), "tables/case_operating_characteristics.csv")
write_csv(do.call(rbind, stopping_rows), "tables/stopping_distribution.csv")
write_csv(mc_results, "validation/monte_carlo_comparison.csv")

checks <- validate_designs(designs)
# Verify vectorized grid engine against DP on deterministic samples spanning the full grid.
for (id in names(designs)) {
  z <- designs[[id]]
  if (is.null(z$grid)) next
  ii <- unique(c(seq(1, nrow(z$grid), length.out = 19) |> as.integer(),
                 which.max(ifelse(z$grid$feasible, z$grid$power, -Inf))))
  discrepancy <- vapply(ii, function(j) {
    g <- z$grid[j, ]; b <- boundary_table(z$spec, as.list(g[c("lambda", "gamma", "eta")]), z$posterior)
    max(abs(c(exact_oc(z$spec$p0, b)$go - g$type1,
              exact_oc(z$spec$p0 + z$spec$improvement, b)$go - g$power)))
  }, numeric(1))
  checks <- rbind(checks, data.frame(check = paste(id, "vectorized grid vs DP"),
               difference = max(discrepancy), tolerance = 1e-11, passed = max(discrepancy) < 1e-11))
}
write_csv(checks, "validation/deterministic_checks.csv")

# Audit the numerical claims of Section 5, rather than force matching in calibration.
kp <- subset(exact_results, design == "case_published")
claim_audit <- data.frame(quantity = c("Type I error", "Null expected N", "Power", "Alternative expected N", "Early efficacy under alternative"),
   manuscript = c(.099, 40.2, .966, 42.7, .662),
   recomputed = c(kp$go[1], kp$en[1], kp$go[2], kp$en[2], kp$early_efficacy[2]),
   rounding_digits = c(3, 1, 3, 1, 3))
claim_audit$matches_rounding <- mapply(function(x, y, d) round(x, d) == y,
                claim_audit$recomputed, claim_audit$manuscript, claim_audit$rounding_digits)
write_csv(claim_audit, "validation/manuscript_case_claims.csv")

curves <- list()
for (id in c("primary_hb", "primary_fixed", "case_published", "case_hb_optimized", "case_fixed",
              sprintf("sensitivity_%02d", seq_len(nrow(settings$sensitivity))))) {
  z <- designs[[id]]
  curves[[id]] <- do.call(rbind, lapply(seq(.10, .70, .005), function(p)
    cbind(design_metadata(id, z), p_t = p, exact_oc(p, z$boundaries))))
}
curve_results <- do.call(rbind, curves)
write_csv(curve_results, "tables/operating_curves.csv")
saveRDS(list(designs = designs, settings = settings, quick = quick, B = B,
             exact = exact_results, mc = mc_results, curves = curve_results,
             checks = checks, claims = claim_audit, calibration = calibration_summary),
        "results/analysis.rds", version = 3)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
if (!all(checks$passed)) stop("Deterministic validation failed; inspect results/validation.")
if (!all(claim_audit$matches_rounding)) stop("Manuscript numerical claims differ at stated precision.")
if (!all(mc_results$covered)) stop("Monte Carlo family validation failed; inspect differences before release.")
message("All validation checks passed. Rendering figures and report.")
source("scripts/make_outputs.R")
source("scripts/write_manifest.R")
