# Reuse the checked guide instead of maintaining duplicate tutorial examples.
guide_section <- function(title) {
  root <- rprojroot::find_package_root_file()
  source_file <- file.path(root, "inst", "examples", "getting-started.R")
  guide <- readLines(file.path(root, "inst", "examples", "getting-started.md"),
                     warn = FALSE)
  signature <- paste0("<!-- Source MD5: ", unname(tools::md5sum(source_file)), " -->")
  if (!signature %in% guide)
    stop("Rebuild the guide with Rscript tools/build-getting-started.R")
  start <- which(guide == paste0("## ", title))
  if (length(start) != 1L)
    stop("Guide section not found: ", title)
  following <- which(grepl("^## ", guide) & seq_along(guide) > start)
  end <- if (length(following)) following[1L] - 1L else length(guide)
  section <- guide[seq.int(start + 1L, end)]
  # The article supplies its own title, so promote the guide's subsections.
  section <- sub("^### ", "## ", section)
  # Article HTML is one directory below docs; pkgdown copies man/figures here.
  gsub("../../man/figures/", "../reference/figures/", section, fixed = TRUE)
}
