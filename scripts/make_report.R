# Data-driven manuscript tables, figure captions and a self-contained HTML report.
# Sourced by make_outputs.R after all exact and MC checks have passed.
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
escape_html <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s, fixed = TRUE), fixed = TRUE), fixed = TRUE)
escape_tex <- function(s) {
  # All manuscript tables contain plain ASCII text; escape LaTeX metacharacters.
  vapply(as.character(s), function(z) {
    z <- gsub("<=", "\u2264", z, fixed = TRUE)
    z <- gsub(">=", "\u2265", z, fixed = TRUE)
    chars <- strsplit(z, "", fixed = TRUE)[[1]]
    paste0(vapply(chars, function(ch) {
      if (ch %in% c("&", "%", "$", "#", "_", "{", "}")) paste0("\\", ch)
      else if (ch == "\u2264") "\\ensuremath{\\leq}"
      else if (ch == "\u2265") "\\ensuremath{\\geq}"
      else if (ch == "<") "\\ensuremath{<}" else if (ch == ">") "\\ensuremath{>}"
      else ch
    }, character(1)), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}
pretty_method <- function(id) {
  ifelse(id == "case_published", "HBOP2-FE (manuscript)",
    ifelse(id == "case_hb_optimized", "HBOP2-FE (optimized)",
      ifelse(grepl("fixed", id), "BOP2-FE", "HBOP2-FE")))
}
manuscript_table <- function(d, filename, caption, note) {
  csv <- file.path("results", "tables", paste0(filename, ".csv"))
  write.csv(d, csv, row.names = FALSE, na = "")
  tex_rows <- apply(d, 1, function(row) paste0(paste(escape_tex(row), collapse = " & "), " \\\\"))
  latex <- c("\\begin{table*}[t]", "\\centering", paste0("\\caption{", escape_tex(caption), "}"),
    paste0("\\label{tab:", filename, "}"), "\\small", "\\setlength{\\tabcolsep}{4pt}",
    "\\resizebox{\\textwidth}{!}{%", paste0("\\begin{tabular}{", paste(c("l", "l", rep("r", ncol(d) - 2)), collapse = ""), "}"),
    "\\toprule", paste0(paste(escape_tex(names(d)), collapse = " & "), " \\\\"), "\\midrule",
    tex_rows, "\\bottomrule", "\\end{tabular}%", "}", "\\par\\vspace{3pt}",
    "\\begin{minipage}{\\textwidth}\\footnotesize", paste0("\\textit{Note.} ", escape_tex(note)),
    "\\end{minipage}", "\\end{table*}")
  writeLines(latex, sub("\\.csv$", ".tex", csv))
  html_rows <- apply(d, 1, function(row) paste0("<tr><td>", paste(escape_html(row), collapse = "</td><td>"), "</td></tr>"))
  html <- paste0('<table class="publication"><caption>', escape_html(caption), '</caption><thead><tr><th>',
    paste(escape_html(names(d)), collapse = "</th><th>"), '</th></tr></thead><tbody>',
    paste(html_rows, collapse = "\n"), '</tbody></table>')
  paste0(html, '<p class="note">', note, '</p>')
}

main <- subset(ex, design %in% c("primary_hb", "primary_fixed"))
main <- main[order(main$scenario, main$method), ]
get_mc <- function(id, sc, metric, field) mc[mc$design == id & mc$scenario == sc & mc$metric == metric, field]
main_table <- data.frame(Scenario = main$scenario, Method = main$method,
  Target = fmt(main$p_c_target, 2), Treatment = fmt(main$p_t, 2),
  `Exact go` = fmt(main$go),
  `MC go` = mapply(function(id, sc) fmt(get_mc(id, sc, "go", "estimate")), main$design, main$scenario),
  MCSE = mapply(function(id, sc) fmt(get_mc(id, sc, "go", "mcse"), 5), main$design, main$scenario),
  `E(N)` = fmt(main$en, 1), `Early futility` = fmt(main$early_futility),
  `Early efficacy` = fmt(main$early_efficacy), `Final go` = fmt(main$final_go), check.names = FALSE)
t1 <- manuscript_table(main_table, "table_1_primary", "Primary operating characteristics conditional on the fixed ECA (60/200).",
  paste0("N = 60; looks = 20, 40, 60; B = ", format(a$B, big.mark = ","),
    " independent current trials per row. S1, S3 and S4 are target nulls; S2, S5 and S6 are target alternatives. Final go is unconditional, including the probability of reaching the final look. E(N) counts current patients only."))
smallcal <- subset(cal, design %in% c("primary_hb", "primary_fixed", "case_published", "case_fixed"))
t2 <- manuscript_table(data.frame(Setting = ifelse(startsWith(smallcal$design, "case_"), "Case", "Primary"),
  Method = smallcal$method, Lambda = fmt(smallcal$lambda, 2), Gamma = fmt(smallcal$gamma, 2), Eta = fmt(smallcal$eta, 2),
  `Type I error` = fmt(smallcal$exact_type1, 6), Power = fmt(smallcal$exact_power, 6),
  `E(N), null` = fmt(smallcal$en_null, 2), `E(N), alternative` = fmt(smallcal$en_alt, 2), check.names = FALSE),
  "table_2_calibration", "Separately calibrated designs on the prespecified 13,000-point grid.",
  "Maximize compatible power subject to exact type I error <= 0.10. Ties: lower null E(N), higher null futility, lower null early efficacy, then increasing lambda, gamma, eta. Mathematically equal criteria are compared to 12 decimal places. The manuscript case design coincides with the selected optimum.")
ss <- subset(ex, startsWith(design, "sensitivity_"))
sw <- do.call(rbind, lapply(unique(ss$design), function(id) {
  d <- ss[ss$design == id, ]; v <- function(sc, x) d[d$scenario == sc, x]
  data.frame(Cap = d$cap[1], Quality = fmt(d$quality[1], 2), ESS = fmt(d$ess[1], 1),
    PMD = fmt(d$pmd[1], 5), Lambda = fmt(d$lambda[1], 2), Gamma = fmt(d$gamma[1], 2), Eta = fmt(d$eta[1], 2),
    S1 = fmt(v("S1", "go")), S2 = fmt(v("S2", "go")), S3 = fmt(v("S3", "go")),
    S4 = fmt(v("S4", "go")), S5 = fmt(v("S5", "go")), S6 = fmt(v("S6", "go")))
}))
t3 <- manuscript_table(sw, "table_3_borrowing", "Recalibrated borrowing sensitivity using the same ECA (60/200).",
  "S1-S6 columns are exact go probabilities. ESS is the borrowed likelihood size, excluding baseline prior size. PMD is comparator posterior mean minus the planning value 0.30; both are deterministic design diagnostics. Equal ESS pairs are equivalent. Distinct ESS values can also yield identical integer decisions after calibration.")
db <- do.call(rbind, lapply(c("primary_hb", "primary_fixed", "case_published", "case_fixed"), function(id) {
  b <- a$designs[[id]]$boundaries
  data.frame(Setting = ifelse(startsWith(id, "case_"), "Case", "Primary"), Method = pretty_method(id),
    N = b$n, CF = fmt(b$cutoff_futility), CE = fmt(b$cutoff_efficacy),
    `Futility / no-go` = paste0("x <= ", b$futility_max),
    `Efficacy / go` = paste0("x >= ", b$efficacy_min), check.names = FALSE)
}))
t4 <- manuscript_table(db, "table_4_decisions", "Protocol response-count decisions.",
  "Apply rules only if the trial has not stopped previously. Interim futility uses PP < CF, efficacy PP >= CE. At N = 60, use the single go threshold lambda. Counts between thresholds continue at an interim look.")
ks <- subset(ex, design %in% c("case_published", "case_fixed", "case_cap10", "case_cap50", "case_discounted"))
kt <- data.frame(Design = ifelse(ks$design == "case_published", "HBOP2-FE, ESS 30",
     ifelse(ks$design == "case_fixed", "BOP2-FE", paste0("HBOP2-FE, ESS ", ks$ess))),
  Scenario = ks$scenario, Target = fmt(ks$p_c_target, 2), Treatment = fmt(ks$p_t, 2),
  `Exact go` = fmt(ks$go), `E(N)` = fmt(ks$en, 1),
  `Early futility` = fmt(ks$early_futility), `Early efficacy` = fmt(ks$early_efficacy), check.names = FALSE)
t5 <- manuscript_table(kt, "table_5_case", "Illustrative case study and recalibrated borrowing sensitivity.",
  "The same aggregate ECA (172/637) is retained in every row. The compatible planning null is exactly 0.27, not 172/637. K1 and K2 are the compatible null and alternative; K3 and K4 are favorable-ECA target nulls; K5 is the unfavorable-ECA alternative. These are design calculations, not analyses of clinical patient records.")

primary_h <- subset(ex, design == "primary_hb"); primary_f <- subset(ex, design == "primary_fixed")
case_h <- subset(ex, design == "case_published")
text_primary <- sprintf("Conditional on the observed ECA (60 responders among 200 subjects), exact calibration selected (lambda, gamma, eta) = (0.77, 0.80, 3.00) for HBOP2-FE and (0.87, 0.50, 2.60) for fixed-null BOP2-FE. Under the compatible null, exact type I error was %.4f and %.4f, respectively. Under the compatible alternative, power was %.4f versus %.4f. This %.2f percentage-point power increase accompanied larger expected current-trial sample sizes: %.2f versus %.2f under the null and %.2f versus %.2f under the alternative. Early efficacy under the compatible alternative was %.4f versus %.4f. Thus, increased compatible power was accompanied by greater patient use; this comparison does not establish general efficiency superiority.",
 primary_h$go[1], primary_f$go[1], primary_h$go[2], primary_f$go[2],
 100 * (primary_h$go[2] - primary_f$go[2]), primary_h$en[1], primary_f$en[1],
 primary_h$en[2], primary_f$en[2], primary_h$early_efficacy[2], primary_f$early_efficacy[2])
text_stress <- sprintf("With the ECA fixed, target-null rates of 0.35 and 0.40 produced false-positive probabilities of %.4f and %.4f for HBOP2-FE, compared with %.4f and %.4f for BOP2-FE. Neither method controlled type I error under these mismatch configurations. For HBOP2-FE the unfavorable-ECA alternatives S5 and S6 had power %.4f and %.4f. S3 and S6 share the same treatment rate, as do S4 and S5; their exact go probabilities are therefore equal despite different scientific interpretations. Without concurrent comparator data, target comparator truth labels the scenario but does not enter the decision rule.",
 primary_h$go[3], primary_h$go[4], primary_f$go[3], primary_f$go[4], primary_h$go[5], primary_h$go[6])
text_sens <- "Recalibration across the 12 cap/quality configurations produced only two distinct integer stopping policies in this grid: ESS 5 through 30 selected one policy, whereas ESS 37.5 and 50 selected another. Consequently many borrowing configurations had exactly identical operating characteristics. Equal ESS pairs are analytically equivalent, while the additional plateaus at unequal ESS arise from discrete decision thresholds and grid optimization. These results provide no evidence here for a smooth monotone robustness gain from lowering the borrowing cap. The posterior mean displacement is a fixed diagnostic, not an empirical distribution across repeated trials."
text_case <- sprintf("For the illustrative ECA (172/637), the manuscript boundary parameters (0.73, 0.50, 2.80) were recovered by the full grid search. Exact type I error at pT = 0.27 was %.6f with expected N %.4f. At pT = 0.47, power was %.6f, expected N %.4f and early efficacy probability %.6f. The response-count rules were futility at x <= 4 and efficacy at x >= 12 at N = 20; x <= 11 and x >= 18 at N = 40; and final go at x >= 21 at N = 60. All five manuscript numerical claims reproduced at their stated precision. At target nulls 0.32 and 0.37, false-positive probabilities increased to %.4f and %.4f. The low cap, high cap and quality-discount alternatives retained the same decision policy after separate recalibration. No patient-level clinical data were available or synthesized as if observed.",
 case_h$go[1], case_h$en[1], case_h$go[2], case_h$en[2], case_h$early_efficacy[2], case_h$go[3], case_h$go[4])
text_validation <- sprintf("Each of %d scenario-design combinations was independently simulated using %s current trials, with prespecified seeds. Exact probability propagation was checked against exhaustive increment-path enumeration, the vectorized calibration engine against forward propagation, and posterior quadrature against reverse integration. All %d deterministic checks passed. All %d exact values lay inside simultaneous Monte Carlo validation intervals with approximately 99%% family coverage (Bonferroni adjustment; exact Clopper-Pearson intervals for probabilities and normal intervals using exact sample-size variance for mean N). The largest absolute standardized discrepancy was %.3f. These diagnostics assess numerical implementation and simulation precision, not external-source transportability.",
 nrow(mc) / 7, format(a$B, big.mark = ","), nrow(a$checks), nrow(mc), max(abs(mc$z)))
methods <- paste0("We conditioned on an analysis-ready completed ECA and simulated only current treatment outcomes. Independent Jeffreys Beta(0.5,0.5) priors were used for treatment and comparator rates. The ECA likelihood power was w = q min(1, cap/n_h), with effective borrowed size w n_h; comparator data did not update the treatment prior. Posterior superiority probabilities were calculated by one-dimensional adaptive numerical integration using upper beta tails (absolute and relative tolerances 1e-10). Boundaries were calibrated separately for each method and borrowing configuration on lambda = 0.50:0.99 by 0.01, gamma = 0.10:1.00 by 0.10, eta = 0.50:3.00 by 0.10 (13,000 points). Exact binomial calculations implemented the complete sequential stopping rule. Compatible power was maximized subject to type I error <= 0.10; prespecified tie-breakers were null expected sample size, null futility and null early efficacy, then increasing boundary parameters. The fixed-null comparator was the planning value (0.30 primary; 0.27 case). Results are conditional on the ECA and the tested boundary family.")
writeLines(c("# Simulation and case-study results", "", "## Methods", "", methods,
             "", "## Primary simulation", "", text_primary, "", text_stress,
             "", "## Borrowing sensitivity", "", text_sens,
             "", "## Illustrative case study", "", text_case,
             "", "## Numerical validation", "", text_validation), "docs/manuscript_results.md")

captions <- c(
  "Figure 1. Exact operating characteristics for separately calibrated HBOP2-FE and fixed-null BOP2-FE, conditional on ECA (60/200). Panel A: probability of go, including early efficacy and final go. Panel B: expected current-trial sample size. Vertical guides mark compatible null 0.30 and alternative 0.45. Horizontal guide marks the nominal type I error ceiling 0.10. N = 60; looks at 20, 40, 60. HBOP2-FE borrows an effective 30 comparator subjects. Lines are exact calculations; no Monte Carlo uncertainty is attached to them.",
  "Figure 2. Robustness to ECA-target mismatch with ECA (60/200) fixed. Each method/borrowing configuration was separately calibrated at target null 0.30. Panel A: target-null go probability at pT = target comparator truth; values above 0.10 indicate false-positive inflation. Panel B: target-alternative go probability at pT = target comparator truth + 0.15. The vertical guide marks ECA compatibility at 0.30. HBOP2-FE effective sizes 10, 20 and 30 produce identical policies and their curves overlap exactly. Colors and line patterns identify designs; overlap is substantive, not missing data.",
  "Figure 3. Exact sensitivity to effective ECA size after recalibration for each of the 12 cap/quality settings. Duplicate effective sizes are plotted once. Panels show type I error (S1, S3, S4), power (S2, S5, S6), and expected current-trial sample size (S1, S2). All configurations use ECA (60/200). Lines connect the evaluated configurations as visual guides. Plateaus reflect identical integer stopping policies, including some configurations with unequal effective size. Scenario definitions are in Table 1 and the scenario CSV.",
  "Figure 4. Complete exact calibration search for both primary methods with ECA (60/200). Each point is a boundary parameter set on the 13,000-point grid; multiple sets can have identical coordinates. Feasible points satisfy the exact compatible-null ceiling 0.10. Gold diamonds show the maximum-power designs selected by the prespecified tie-breakers. The displays are conditional calibration surfaces, not confidence regions.",
  "Figure 5. Integer response-count decision regions for the illustrative ECA (172/637), comparing separately calibrated HBOP2-FE and fixed-null BOP2-FE. The manuscript HBOP2-FE design and the independently calibrated full-grid optimum coincide at (lambda, gamma, eta) = (0.73, 0.50, 2.80). BOP2-FE uses (0.90, 0.70, 2.60). Blank space beyond each row's sample size represents impossible counts. F denotes maximum futility/no-go count and E minimum efficacy/go count. Apply decisions only to trials still ongoing. ECA effective size is 30 for HBOP2-FE. At the final look there is no continuation region.",
  "Figure 6. Exact illustrative case-study curves with ECA (172/637) fixed. The HBOP2-FE manuscript design is identical to the independently recovered full-grid optimum; one curve represents both. Separately optimized fixed-null BOP2-FE provides the benchmark. Panel A: probability of go. Panel B: expected current-trial sample size. Vertical guides mark the planning null 0.27 and alternative 0.47; the empirical ECA rate 172/637 is not substituted for the planning null.",
  paste0("Figure 7. Independent Monte Carlo validation for all ", nrow(mc) / 7, " scenario-design combinations (B = ", format(a$B, big.mark = ","), " per combination). Histograms show differences from exact values standardized by exact reference SE for go probability and expected N. Reference SE uses exact Bernoulli variance or the exact stopping-sample-size variance, divided by B. Dashed guides at +/-3 aid inspection; they are not the formal simultaneous pass threshold. All seven outcome metrics and family-adjusted intervals are available in monte_carlo_comparison.csv. These comparisons are correlated across metrics; no normality test is implied.")
)
writeLines(paste0(seq_along(captions), ". ", captions), "docs/figure_captions.md")
escape_html <- function(s) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", s, fixed = TRUE), fixed = TRUE), fixed = TRUE)
para <- function(s) paste0("<p>", escape_html(s), "</p>")
figure_files <- c("figure_1_primary_operating_curves", "figure_2_mismatch_robustness", "figure_3_borrowing_sensitivity",
                  "figure_4_exact_calibration", "figure_5_case_decision_regions", "figure_6_case_operating_curves",
                  "figure_7_monte_carlo_validation")
fig_html <- vapply(seq_along(figure_files), function(i) paste0('<figure><img src="figures/', figure_files[i],
  '.svg" alt="Figure ', i, '"><figcaption>', escape_html(captions[i]), '</figcaption></figure>'), character(1))
html <- c('<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
  '<title>HBOP2-FE | Reproducible simulation and case study</title><style>',
  'body{font:16px/1.65 Georgia,serif;color:#24313a;background:#fff;margin:0}main{max-width:1120px;margin:44px auto;padding:0 28px}h1,h2,h3{font-family:Arial,sans-serif;line-height:1.25}h1{font-size:32px;color:#244B73}h2{font-size:23px;margin-top:42px;border-top:1px solid #ccd4da;padding-top:24px}.eyebrow{font:12px Arial,sans-serif;letter-spacing:.14em;color:#65727c;text-transform:uppercase}.summary{background:#eef3f6;padding:18px 22px;border-left:4px solid #244B73}.note,figcaption{font:13px/1.6 Arial,sans-serif;color:#52616b}table{border-collapse:collapse;width:100%;font:12px/1.5 Arial,sans-serif;margin:18px 0}caption{font-weight:bold;text-align:left;color:#24313a;padding:10px 0}thead{border-top:2px solid #344a5b;border-bottom:1px solid #344a5b}th,td{padding:7px 6px;text-align:right;white-space:nowrap}th:first-child,td:first-child{text-align:left}tbody{border-bottom:2px solid #344a5b}tbody tr:nth-child(even){background:#f4f6f7}.table-wrap{overflow-x:auto}figure{margin:28px 0 42px}figure img{width:100%;height:auto}figcaption{margin-top:10px}a{color:#244B73}code{font-size:13px}footer{margin:40px 0;border-top:1px solid #ccd4da;padding-top:18px;font:13px Arial,sans-serif}</style></head><body><main>',
  '<div class="eyebrow">Manuscript v29 / R reproducibility study</div><h1>HBOP2-FE simulation and case study</h1>',
  paste0('<div class="summary"><strong>', if (a$quick) "SMOKE RUN: not for manuscript use." else "Full run: 100,000 independent current trials per scenario-design combination.",
    '</strong><br>Fixed completed ECA; exact full-grid calibration; independent Monte Carlo validation; seven exportable figures.</div>'),
  '<h2>Methods</h2>', para(methods), '<h2>Primary operating characteristics</h2>', para(text_primary),
  paste0('<div class="table-wrap">', t1, '</div>'), fig_html[1],
  '<h2>ECA-target mismatch</h2>', para(text_stress), fig_html[2],
  '<h2>Borrowing sensitivity</h2>', para(text_sens), paste0('<div class="table-wrap">', t3, '</div>'), fig_html[3],
  '<h2>Calibration and protocol boundaries</h2>', paste0('<div class="table-wrap">', t2, '</div>'), fig_html[4],
  paste0('<div class="table-wrap">', t4, '</div>'),
  '<h2>Illustrative case study</h2>', para(text_case), fig_html[5], fig_html[6], paste0('<div class="table-wrap">', t5, '</div>'),
  '<h2>Numerical validation</h2>', para(text_validation), fig_html[7],
  '<h2>Interpretation and scope</h2>',
  para("Calibration controls the tested compatible null only. Severe false-positive inflation under ECA mismatch remains a substantive limitation. Because the two methods have different discrete feasible decision families, differences at their optimized designs cannot be interpreted as a universal benefit of borrowing. This study does not evaluate patient-level propensity-score estimation, randomized concurrent controls, clinical efficacy, delayed outcomes, or unconditional repeated ECA generation. No clinical recommendations follow from these illustrative calculations."),
  '<footer>All displayed numerical results are generated from results/analysis.rds. Source formulas are mapped in docs/methods_audit.md. Raw exact values, simulation estimates, seeds and validation intervals are provided as CSV. Figure PDFs are vector output; SVG, 300 dpi PNG and 600 dpi TIFF are also included. See README.md for reproduction.</footer>',
  '</main></body></html>')
writeLines(html, "results/report.html")
