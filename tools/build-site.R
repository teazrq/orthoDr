# Run from the repository root using shared R packages and the installed source.
options(rgl.useNULL = TRUE)
source("tools/build-getting-started.R", local = new.env(parent = globalenv()))
pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)

# Keep downloads identical to the scripts exercised by the package checks.
destination <- file.path("docs", "examples")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
scripts <- file.path("inst", "examples", c("getting-started.R", "survival-optimizer.R"))
if (!all(file.copy(scripts, destination, overwrite = TRUE)))
  stop("Could not copy the downloadable examples")
