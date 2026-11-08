# Reproduction and Publication Guide

This project implements the supplied manuscript version 29. It is independent of the older workspace folders `HBOP2FE/` and `HBOP2FE V3/`. The analysis conditions on a completed external control arm (ECA), as specified in the current manuscript. The ECA is neither regenerated nor replaced across simulation scenarios.

## Where to Start

1. Open `results/report.html` in a local browser to review the complete report, figures, tables and interpretation.
2. Use `docs/manuscript_results.md` as a starting point for the manuscript Methods and Results sections. These paragraphs are generated from the computed results.
3. Find complete English figure captions in `docs/figure_captions.md`.
4. Integrate `results/tables/table_1_primary.tex` through `table_5_case.tex` into a LaTeX manuscript. Load the `booktabs` and `graphicx` packages.
5. Use the full-precision CSV files for further calculations. They retain numerical results, calibration grids, random seeds, Monte Carlo standard errors (MCSEs) and validation checks. Formatted manuscript tables generally round probabilities to three decimal places and should not be used to recompute results.

## Run the Analysis in R

Set the working directory to the project root (`HBOP2_FE` after cloning the GitHub repository), and run:

```r
source("scripts/install_dependencies.R")
source("tests/run_tests.R")
source("scripts/run_analysis.R")
source("scripts/verify_release.R")
```

The `Rscript --vanilla` commands in the README provide a clean command-line execution path. When using a restored renv environment, follow the README instructions to load that environment explicitly.

The default command performs the full analysis, with 100,000 independent simulated current trials per scenario-design combination. Calibration uses exact probabilities and is therefore unaffected by Monte Carlo error. The optional `--quick` mode produces smoke-test outputs rather than manuscript results; a default full run restores the manuscript outputs.

## Findings for the Manuscript

- In the primary compatible setting, HBOP2-FE has an exact type I error of 0.099412 and power of 0.866070. The separately calibrated fixed-null BOP2-FE benchmark has corresponding values of 0.095683 and 0.839265.
- The approximately 2.68 percentage-point increase in power accompanies larger expected current-trial sample sizes. It should not be described as an unqualified efficiency improvement.
- Both methods show substantial false-positive inflation when the ECA is incompatible with the target population. HBOP2-FE go probabilities under the S3 and S4 target nulls are 0.3260 and 0.6340. These results do not establish type I error control under mismatch.
- Separate recalibration of the 12 cap/quality configurations produces only two distinct integer stopping policies. Overlapping curves and sensitivity plateaus reflect the discrete decisions and optimization results.
- The case-study parameters `(0.73, 0.50, 2.80)` and all five principal numerical claims reproduce successfully. The full grid search selects the same case-study design.
- The compatible case-study null is exactly **0.27**, as specified in the manuscript. It is not replaced by **172/637**. Although the values are close, the distinction matters for exact verification.

## Figure Selection

| Figure | Content |
|---|---|
| 1 | Primary go-probability and expected-sample-size curves. |
| 2 | Target-population mismatch with the ECA fixed, showing false-positive inflation and loss of power. |
| 3 | Borrowing sensitivity after separate recalibration, indexed by effective sample size with equivalent duplicate points removed. |
| 4 | Complete 13,000-point exact calibration searches for both primary methods. |
| 5 | Case-study integer decision regions for HBOP2-FE and fixed-null BOP2-FE. The manuscript HBOP2-FE design coincides with the recovered optimum. |
| 6 | Case-study go-probability and expected-sample-size curves. |
| 7 | Standardized Monte Carlo discrepancies from exact results, suitable for supplementary material. |

Each figure is supplied as vector PDF/SVG, 300 dpi PNG and 600 dpi TIFF. The default width is 7.2 inches. Select the final dimensions and font sizes according to the target journal's requirements. Figures 1-3 and 5-6 can support the main manuscript; calibration and validation figures can be placed in the supplement.

## GitHub and Analysis Scope

The project is hosted at [yuexi706/HBOP2_FE](https://github.com/yuexi706/HBOP2_FE). The repository contains the code, configuration, tests, dependency lockfile, exact results, figures, tables, license and GitHub Actions workflow. The source manuscript PDF is not redistributed; its provenance and checksum are recorded instead.

The case study is an **illustrative design calculation based on the aggregate counts supplied in the manuscript**, not a reanalysis of KEYNOTE-042 participant records. No clinical patient data have been fabricated. The manuscript reserves propensity-score estimation and residual-confounding evaluation for a separate patient-level simulation study; those procedures are outside this primary analysis.
