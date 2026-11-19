test_that("survival entry points reject nonfinite times by original row", {
  set.seed(465)
  x <- matrix(rnorm(48 * 2), 48, 2)
  y <- c(48:25, 24:1)
  event <- rep(c(1, 0), 24)
  bad <- y
  bad[25] <- NA_real_  # A failure at the 24th order statistic.
  expect_error(orthoDr_surv(x, bad, event, ndr = 1), "observation rows 25", fixed = TRUE)
  bad[2] <- Inf         # A censored subject, also part of the risk set.
  expected <- "observation rows 2, 25"
  expect_error(CP_SIR(x, bad, event), expected, fixed = TRUE)
  expect_error(orthoDr_surv(x, bad, event, ndr = 1, B.initial = matrix(c(1, 0))),
               expected, fixed = TRUE)
  for (type in c("2D", "3D")) {
    expect_error(view_dr_surv(x, bad, event, B = c(1, 0), type = type,
                              FUN = function(y) stop("transformation called")),
                 expected, fixed = TRUE)
  }
  saved <- structure(list(x = x, y = bad, censor = event, B = matrix(c(1, 0)),
                           bw = .7, keep.data = TRUE),
                     class = c("orthoDr", "fit", "surv"))
  expect_error(predict(saved, x), "observation rows 2, 25.*Correct the original follow-up times and refit")
})

test_that("survival fitting and prediction preserve supported time objects", {
  x <- matrix(c(-2, -1, .2, .4, 1, 2))
  time <- setNames(c(0L, 2L, -1L, 0L, 1L, 2L), letters[1:6])
  event <- c(1, 0, 0, 1, 0, 1)
  times <- list(time, matrix(time, dimnames = list(names(time), "hours")),
                setNames(as.Date(time, origin = "2020-01-01"), names(time)),
                setNames(as.POSIXct(time, origin = "2020-01-01", tz = "America/Chicago"), names(time)),
                setNames(as.POSIXlt(as.POSIXct(time, origin = "2020-01-01", tz = "UTC")), names(time)),
                setNames(as.difftime(time, units = "hours"), names(time)))
  ref <- orthoDr_surv(x, time, event, keep.data = TRUE, maxitr = 2, ncore = 1)
  for (y in times) {
    expect_equal(CP_SIR(x, y, event), CP_SIR(x, time, event))
    fit <- orthoDr_surv(x, y, event, keep.data = TRUE, maxitr = 2, ncore = 1)
    expect_identical(fit$y, y)
    pred <- predict(fit, x)
    values <- if (is.matrix(y)) y[, 1L] else y
    expect_identical(pred$timepoints, values[c(1L, 4L, 6L)])
    expect_equal(pred$surv, predict(ref, x)$surv)
  }

  # The difference between the distinct times would overflow R's integer type.
  ref$y <- c(-2147483647L, 2147483647L, -2147483647L,
             -2147483647L, 2147483647L, 2147483647L)
  expect_no_warning(out <- predict(ref, x))
  expect_identical(out$timepoints, ref$y[c(1L, 4L, 6L)])
  expect_equal(out$surv[1L, ], out$surv[2L, ])
})

test_that("survival plots validate one transformed coordinate per input row", {
  x <- matrix(1:6)
  t <- setNames(as.difftime(c(-2, -1, 0, 0, 1, 2), units = "hours"), letters[1:6])
  event <- c(1, 0, 1, 0, 1, 1)
  points <- surface <- NULL
  local_mocked_bindings(
    plot.new = function(...) NULL,
    par = function(...) NULL,
    plot = function(x, y, ...) { points <<- y },
    surface3d = function(x, y, z, ...) { surface <<- y },
    box3d = function(...) NULL,
    axis3d = function(...) NULL,
    mtext3d = function(...) NULL,
    .package = "orthoDr"
  )
  for (type in c("2D", "3D")) {
    calls <- 0L
    transform <- function(y) {
      expect_identical(y, t)
      calls <<- calls + 1L
      y
    }
    view_dr_surv(x, t, event, B = 1, FUN = transform, type = type, legend.add = FALSE)
    expect_identical(calls, 1L)
    expect_error(view_dr_surv(x, t, event, B = 1, FUN = function(y) 1, type = type),
                 "FUN(y) must have one value per observation", fixed = TRUE)
    expect_error(view_dr_surv(x, t, event, B = 1, FUN = function(y) matrix(0, 3, 2), type = type),
                 "FUN(y) must be a real numeric or time vector or a one-column matrix", fixed = TRUE)
    expect_error(view_dr_surv(x, t, event, B = 1, type = type,
                              FUN = function(y) c(-Inf, 0, 1, 1, 2, NaN)),
                 "observation rows 1, 6", fixed = TRUE)
  }
  expect_identical(points, t)
  expect_equal(surface[1L, ] * 4 / 1.5, -2:2)
})
