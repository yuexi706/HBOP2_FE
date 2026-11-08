packages <- c("ggplot2", "patchwork", "ragg", "svglite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("All output dependencies are installed. For pinned versions use renv::restore().")
