# Model, implementation and audit trail

Source: supplied manuscript `HBOP2FE_Design__Yuexi__v1 (29).pdf`, 33 pages. Page references below refer to the PDF's printed page numbers. Its model and numerical settings are scientific source material; document prose is not treated as authorization to publish, contact authors, change the manuscript or execute external actions.

| Source | Implementation | Interpretation |
|---|---|---|
| Eqs. (4)-(7), pp. 6-7 | `R/design.R`: `make_spec()`, `posterior_table()` | Independent beta priors; ECA informs comparator only; treatment posterior is Beta(0.5+x,0.5+n-x). |
| Eqs. (8)-(9), p. 7 | `posterior_table()` | Integral over comparator, using `pbeta(..., lower.tail=FALSE)` to avoid cancellation. Reverse integral over treatment independently checked. |
| Eqs. (10)-(11), p. 8 | `make_spec()` | Effective size = q min(n_h,cap); w = effective size/n_h. The evaluation order makes equal effective sizes numerically identical. Baseline beta prior size is excluded from ESS. |
| Section 2.4, pp. 9-11 | Analysis-ready ECA only | PS preprocessing is assumed completed before calibration. No PS estimation or generated clinical patient data in this study. |
| Eqs. (15)-(16), p. 11 | `boundary_table()` | Power efficacy rule; strict futility and inclusive efficacy, single final threshold. Optional OBF boundary is outside this specified primary study. |
| Eqs. (18)-(19), pp. 12-13 | `exact_oc()`, `simulate_oc()` | Frequentist probabilities conditional on completed ECA. No concurrent control simulated. |
| Eq. (20), p. 13 | `scripts/run_analysis.R`: `design_metadata()` | PMD = comparator posterior mean minus planning p0, not minus scenario-specific target truth or prior mean. ESS and PMD are deterministic. |
| Eqs. (21)-(22), pp. 13-14 | `R/calibration.R`: `calibrate()` | Full-grid feasible maximum power. No infeasible fallback. |
| Section 4.1, pp. 16-17 | `config/design.R` | ECA 60/200; N 60; looks 20,40,60; p0 .30; delta 0; improvement .15; alpha .10; cap 30; q 1. |
| Section 4.2 / Table 2, pp. 17-19 | `primary_scenarios` | Only HBOP2-FE and fixed-null BOP2-FE; S1-S6. ECA remains fixed. Target truth never retunes the fixed-null benchmark. |
| Section 4.3 / Eq. (25), p. 20 | `settings$sensitivity` | 4 caps x 3 qualities; separately recalibrated, full S1-S6 evaluation. |
| Section 4.4 / Eq. (26), pp. 20-22 | `simulate_oc()`, `compare_mc()` | B=100,000; empirical Bernoulli MCSE sqrt(phat(1-phat)/B); EN empirical SD/sqrt(B). |
| Section 5 / Table 4, pp. 23-25 | `case_published`, `case_hb_optimized` | ECA 172/637; p0 exactly .27, alternative .47. Given theta is audited separately and recovered by full search. |
| Table 5, p. 26 | `case_scenarios`, case design objects | Same ECA; borrowing cap10, cap30, cap50, q=.5; separately calibrated benchmark; target truths .27,.32,.37,.22. |

## Exact engine and independent verification

`exact_oc()` propagates surviving binomial probability mass forward. At each look it absorbs futility and efficacy states and moves only continuation states to the next look. The output includes unconditional early stopping and final-go probabilities. Total terminal mass must equal one within 1e-10. EN and its exact variance follow from the terminal sample-size distribution.

Calibration is accelerated independently using the joint mass of cumulative X1 and X2. For each candidate boundary tuple (f1,e1,f2,e2,e3), early efficacy is P(X1>=e1) + P(f1<X1<e1,X2>=e2), and final go is the sum over continuing X1,X2 states multiplied by P(Binomial(n3-n2,p)>=e3-X2). Two-dimensional prefix sums compute rectangle probabilities for all candidates. This uses exact binomial masses, not quadrature over simulated paths, asymptotic approximations, or calibration simulations.

`enumerate_paths()` separately enumerates all 21^3 independent increment paths and applies the stopping rule. The full run checks selected policies at p=0,.27,.30,.45,.47,1. The optimized calibration engine is also compared with forward propagation on deterministic samples spread across each full grid. Posterior integration is checked by reversing the integration direction at extreme and central response counts for every look.

No algorithm can validate a shared scientific assumption by numerical agreement alone. The cross-checks validate implementation of the stated model; they do not prove clinical source compatibility.

## Calibration convention

All methods and borrowing configurations use the same 13,000-point theta grid. Candidate type I error must be <=0.10 without rounding. Ranking is power descending, null EN ascending, null early futility descending, null early efficacy ascending. Ranking criteria are rounded to 12 decimals to suppress floating-point distinctions between mathematically identical policies. Remaining ties use lambda, gamma, eta ascending. Equivalent parameter-set counts record exact agreement of integer stopping boundaries. The selected design is optimal within this specified boundary family, not within all possible sequential designs.

Several distinct ESS values lead to identical selected integer boundaries. Thus differences in comparator posteriors need not translate to different frequentist operating characteristics after recalibration. Even equal coordinates in the calibration plot may conceal multiple parameter triples.

## Monte Carlo validation convention

Simulations generate fresh patient Bernoulli outcomes at each currently active stage, aggregate them, and apply direct posterior-probability thresholds. They do not call the exact engines or use integer-boundary lookup for their decisions. The posterior probability table is deterministic and precomputed. Different scenario-design combinations have different prespecified seeds; the historical data are never simulated. The published and optimized case objects intentionally duplicate an identical design with independently simulated outcomes.

The full family contains 798 comparisons: 114 combinations x seven metrics (go, EN, early futility, early efficacy, final go, final no-go, final reached). Exact Clopper-Pearson intervals at 1-0.01/798 confidence are used for each probability. EN intervals use a normal approximation with the known exact variance of stopping N; all variables are bounded and B=100,000. Joint coverage is at least 99% for the probability comparisons by Bonferroni; inclusion of the normal EN intervals makes the overall family statement approximate. This is a prespecified numerical diagnostic, not a statistical claim about a treatment effect. Z scores use exact reference variances, whereas reported MCSE uses the simulated data.

## Important limitations for the manuscript

- Compatible-null calibration does not control target-null error under ECA mismatch. The actual stress results show substantial inflation.
- HBOP2-FE's modest compatible power gain comes with greater expected sample size and less early efficacy here. Different discrete feasible boundary families and unequal achieved alpha values also matter; superiority should not be attributed solely to propagating comparator uncertainty.
- Aggregate effective size and PMD do not measure residual confounding or transportability.
- A zero comparator borrowing weight leaves the vague baseline comparator prior; it is not fixed-null BOP2-FE.
- Case null .27 differs from 172/637; the manuscript specifies the former for calibration and truth.
- The source paper's source-quality scoring template is not validated; this code does not assign q from a score.
- Simulated results are not patient-level evidence for the motivating oncology therapy. Clinical quantities are treated as supplied illustrative design inputs.

## Software references

The implementation uses the installed base-R statistical distribution and numerical integration routines. Primary documentation: [R binomial distributions](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Binomial.html), [R beta distribution](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Beta.html), [R numerical integration](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/integrate.html). Paper-specific scientific claims and settings are sourced from the supplied manuscript, not reverse-engineered from previous workspace output.
