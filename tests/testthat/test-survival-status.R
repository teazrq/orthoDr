test_that("survival entry points reject invalid status in original observation order", {
  x <- matrix(seq_len(7L))
  y <- c(7, 3, 5, 1, 6, 2, 4)
  censor <- c(0, 2, NA, NaN, Inf, -1, .5)
  expected <- paste("Invalid censor at 1-based observation rows",
                    "2 (2), 3 (NA), 4 (NaN), 5 (Inf), 6 (-1), 7 (0.5).",
                    "Use 0 (censored) or 1 (event).")

  expect_error(orthoDr_surv(x, y, censor, B.initial = matrix(1)), expected, fixed = TRUE)
  expect_error(CP_SIR(x, y, matrix(censor)), expected, fixed = TRUE)
  for (type in c("2D", "3D"))
    expect_error(view_dr_surv(x, y, censor, B = 1, type = type), expected, fixed = TRUE)

  saved <- structure(list(x = x, y = y, censor = censor, B = matrix(1),
                          bw = .7, keep.data = TRUE),
                     class = c("orthoDr", "fit", "surv"))
  expect_error(predict(saved, x),
               paste(expected, "Recode the event indicators and refit the model."),
               fixed = TRUE)
})

test_that("survival status validation checks its shape and bounds diagnostics", {
  x <- matrix(1:20)
  expect_error(CP_SIR(x, 1:20, matrix(rep(1, 20), ncol = 2)),
               "numeric or logical vector or a one-column matrix")
  expect_error(CP_SIR(x, 1:20, rep("1", 20)), "numeric or logical vector")
  expect_error(CP_SIR(x, 1:20, rep(1, 19)), "one value per observation")
  expect_error(CP_SIR(x, 1:20, c(rep(2, 18), 0, 1)), "8 (2) and 10 more", fixed = TRUE)
})

test_that("valid status formats preserve retained data and survival predictions", {
  x <- matrix(c(-2, -1, .2, .4, 1, 2))
  y <- c(1, 4, 3, 2, 2, 5)
  status <- c(1L, 0L, 0L, 1L, 0L, 1L)
  reference <- orthoDr_surv(x, y, status, keep.data = TRUE, maxitr = 2, ncore = 1)
  expected <- predict(reference, x)

  named <- setNames(as.logical(status), letters[seq_along(status)])
  column <- matrix(status, dimnames = list(names(named), "censor"))
  for (censor in list(named, column)) {
    fit <- orthoDr_surv(x, y, censor, keep.data = TRUE, maxitr = 2, ncore = 1)
    expect_identical(fit$censor, censor)
    expect_equal(predict(fit, x), expected)
  }
})

test_that("both survival plots accept all-censored observations", {
  x <- matrix(1:4)
  shown <- surface <- NULL
  local_mocked_bindings(
    plot.new = function(...) NULL,
    par = function(...) NULL,
    plot = function(x, y, col, ...) { shown <<- col },
    surface3d = function(x, y, z, ...) { surface <<- z },
    box3d = function(...) NULL,
    axis3d = function(...) NULL,
    mtext3d = function(...) NULL,
    .package = "orthoDr"
  )
  view_dr_surv(x, 1:4, matrix(FALSE, 4L, 1L), B = 1, legend.add = FALSE)
  expect_true(all(shown == "red"))
  view_dr_surv(x, 1:4, integer(4L), B = 1, FUN = identity, type = "3D")
  expect_equal(surface, matrix(1, 100L, 4L))
})
