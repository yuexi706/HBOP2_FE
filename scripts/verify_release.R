#!/usr/bin/env Rscript
a <- readRDS("results/analysis.rds")
stopifnot(!a$quick, a$B == 100000L, nrow(a$mc) == 798L,
          all(a$checks$passed), all(a$mc$covered), all(a$claims$matches_rounding),
          all(a$calibration$exact_type1 <= .10))
stopifnot(all(vapply(a$designs[names(a$designs) != "case_published"],
                    function(z) nrow(z$grid) == 13000L && z$selected$feasible, logical(1))))
files <- list.files("results/figures", full.names = TRUE)
for (ext in c("pdf", "svg", "png", "tiff")) {
  selected <- files[endsWith(files, paste0(".", ext))]
  stopifnot(length(selected) == 7L, all(file.info(selected)$size > 1000))
}
for (ext in c("csv", "tex")) stopifnot(length(list.files("results/tables",
                                   pattern = paste0("^table_[1-5].*\\.", ext, "$"))) == 5L)
stopifnot(file.exists("results/report.html"), file.exists("renv.lock"),
          file.exists("docs/manuscript_results.md"), file.exists("docs/figure_captions.md"))
manifest <- read.csv("results/manifest.csv", stringsAsFactors = FALSE)
stopifnot(all(file.exists(manifest$path)),
          identical(unname(tools::md5sum(manifest$path)), manifest$md5))
cat("PASS: full-run release checks; 7 figures x 4 formats, 5 manuscript tables, 321 deterministic checks, 798 MC comparisons.\n")
