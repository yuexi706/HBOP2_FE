# MD5 checksums detect accidental edits to source and generated files. They are
# integrity checks, not a claim that platform-dependent PDF/font bytes are fixed.
files <- sort(list.files(c("R", "scripts", "config", "tests", "docs", "results", ".github"),
                        full.names = TRUE, recursive = TRUE, all.files = TRUE))
files <- files[!grepl("(^|/)(run|render)\\.log$|(^|/)manifest\\.csv$", files)]
files <- files[!basename(files) %in% c(".DS_Store", ".Rhistory", ".RData", ".Ruserdata")]
files <- sort(c(files, "README.md", "LICENSE", "CITATION.cff", "renv.lock", ".gitignore"))
write.csv(data.frame(path = files, bytes = file.info(files)$size,
                     md5 = unname(tools::md5sum(files))),
          "results/manifest.csv", row.names = FALSE)
