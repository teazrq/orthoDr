# Exercise the public, downloadable guide as part of R CMD check.
library(orthoDr)
local({
  e <- new.env(parent = globalenv())
  sys.source(system.file("examples", "getting-started.R", package = "orthoDr"),
             envir = e)
  source_hash <- unname(tools::md5sum(system.file(
    "examples", "getting-started.R", package = "orthoDr")))
  guide <- readLines(system.file("examples", "getting-started.md", package = "orthoDr"))
  stopifnot(
    any(grepl(paste0("Source MD5: ", source_hash), guide, fixed = TRUE)),
    identical(dim(e$reg_fit$B), c(4L, 1L)),
    length(e$reg_prediction$pred) == 5L,
    all(is.finite(e$reg_prediction$pred)),
    ncol(e$surv_prediction$surv) == 5L,
    nrow(e$surv_prediction$surv) == sum(e$event),
    identical(e$surv_prediction$timepoints, sort(e$observed_time[e$event == 1])),
    all(is.finite(e$surv_prediction$surv)),
    all(e$surv_prediction$surv >= 0 & e$surv_prediction$surv <= 1),
    all(apply(e$surv_prediction$surv, 2, function(s) all(diff(s) <= 1e-12))),
    length(e$dose_prediction$pred) == 5L,
    all(is.finite(e$dose_prediction$pred)),
    length(e$pseudo_prediction$pred) == 5L,
    all(is.finite(e$pseudo_prediction$pred)),
    e$pca_distance < 1e-4,
    abs(drop(crossprod(e$optim_fit$B)) - 1) < 1e-8
  )
})

# Keep the original public survival example runnable, including both 3D calls.
check_survival_example <- function() {
  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  on.exit(rgl::close3d(), add = TRUE)
  e <- new.env(parent = globalenv())
  sys.source(system.file("examples", "survival-optimizer.R", package = "orthoDr"),
             envir = e)
  stopifnot(all(is.finite(e$orthoDr.fit$B)), is.finite(e$orthoDr.fit$fn))
}
check_survival_example()
