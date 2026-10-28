# Exercise the public, downloadable guide as part of R CMD check.
library(orthoDr)
check_getting_started <- function() {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  e <- new.env(parent = globalenv())
  sys.source(system.file("examples", "source-md5.R", package = "orthoDr"),
             envir = environment())
  hash_cases <- replicate(3L, tempfile())
  on.exit(unlink(hash_cases), add = TRUE)
  writeBin(charToRaw("x <- 1\n"), hash_cases[1L])
  writeBin(charToRaw("x <- 1\r\n"), hash_cases[2L])
  writeBin(charToRaw("x <- 2\n"), hash_cases[3L])
  hashes <- unname(vapply(hash_cases, canonical_source_md5, character(1)))
  stopifnot(identical(hashes[1L], hashes[2L]),
            !identical(hashes[1L], hashes[3L]))
  sys.source(system.file("examples", "getting-started.R", package = "orthoDr"),
             envir = e)
  source_hash <- canonical_source_md5(system.file(
    "examples", "getting-started.R", package = "orthoDr"))
  guide <- readLines(system.file("examples", "getting-started.md", package = "orthoDr"))
  stopifnot(
    any(grepl(paste0("Source MD5: ", source_hash), guide, fixed = TRUE)),
    identical(dim(e$reg_fit$B), c(4L, 1L)),
    identical(dim(e$reduced_x), c(nrow(e$x_reg), 1L)),
    all(is.finite(e$reduced_x)),
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
    isTRUE(e$dose_fit$converge),
    identical(e$dose_fit$stop.reason, "gradient_tolerance"),
    all(e$dose_fit$gradient.norm <= 1e-6),
    abs(e$dose_fit$fn + e$dose_fit$value - e$dose_fit$penalty) < 1e-12,
    length(e$pseudo_prediction$pred) == 5L,
    all(is.finite(e$pseudo_prediction$pred)),
    identical(dim(e$pseudo_fit$B), c(3L, 2L)),
    e$pca_distance < 1e-4,
    abs(drop(crossprod(e$optim_fit$B)) - 1) < 1e-8
  )
}
check_getting_started()

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
