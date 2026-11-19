# Independent grouped weighted product-limit reference for the plotted steps.
view_surv_reference <- function(time, event, weights, grid) {
  vapply(grid, function(t) {
    probability <- 1
    for (jump in sort(unique(time[event == 1 & time <= t]))) {
      deaths <- sum(weights[time == jump & event == 1])
      if (deaths > 0)
        probability <- probability *
          (1 - deaths / sum(weights[time >= jump]))
    }
    probability
  }, numeric(1))
}

test_that("3D survival uses observed empirical percentiles and every event", {
  n <- 200L
  x <- matrix(seq(-1, 1, length.out = n), ncol = 1)
  y <- (seq_len(n) / n)^4 + .03
  event <- as.numeric(seq_len(n) %in% c(3L, 7L, 151L))
  transformed <- log(y) - mean(log(y))
  expected <- sort(unique(transformed[c(1L, seq(2L, n, 2L), which(event == 1))]))
  calls <- 0L
  shown <- NULL
  local_mocked_bindings(
    surface3d = function(x, y, z, ...) {
      shown <<- list(x = x, y = y, z = z)
    },
    box3d = function(...) NULL,
    axis3d = function(...) NULL,
    mtext3d = function(...) NULL,
    .package = "orthoDr"
  )
  expect_no_warning(view_dr_surv(x, y, event, B = 1, bw = .7,
    FUN = function(t) { calls <<- calls + 1L; log(t) - mean(log(t)) }, type = "3D"))
  expect_identical(calls, 1L)
  expect_equal(dim(shown$z), c(100L, length(expected)))
  expect_equal(shown$y[1, ] * diff(range(transformed)) / 1.5, expected,
               tolerance = 1e-14)
  expect_false(any(transformed[c(5L, 9L)] %in% expected))
  expect_true(all(transformed[event == 1] %in% expected))
  xgrid <- seq(min(x), max(x), length.out = 100)
  for (i in c(1L, 37L, 100L)) {
    weights <- exp(-.5 * ((xgrid[i] - x[, 1]) / .7)^2)
    expect_equal(shown$z[i, ],
      view_surv_reference(transformed, event, weights, expected), tolerance = 1e-12)
  }
})

test_that("surface steps preserve ties and zero-weight tails without recycling", {
  x <- matrix(c(0, 0, 100, 100, 100, 100), ncol = 1)
  y <- c(1, 1, 2, 3, 4, 4)
  event <- c(1, 0, 1, 1, 1, 0)
  shown <- NULL
  local_mocked_bindings(
    surface3d = function(x, y, z, ...) { shown <<- z },
    box3d = function(...) NULL,
    axis3d = function(...) NULL,
    mtext3d = function(...) NULL,
    .package = "orthoDr"
  )
  expect_no_warning(view_dr_surv(x, y, event, B = 1, bw = .001,
                                FUN = identity, type = "3D"))
  expect_equal(dim(shown), c(100L, 4L))
  expect_equal(shown[1, ], c(.5, .375, .25, .125))
  expect_equal(shown[100, ], c(1, .75, .5, .25))
  expect_true(all(is.finite(shown)))
  expect_true(all(apply(shown, 1, function(s) all(diff(s) <= 0))))
  expect_true(all(shown >= 0 & shown <= 1))
  expect_error(view_dr_surv(x, y, event, B = 1, bw = 1e-300,
                           FUN = identity, type = "3D"),
               "Cannot compute finite Gaussian surface weights")
})

test_that("actual non-square rgl surfaces keep predictor and time orientation", {
  old <- options(rgl.useNULL = TRUE)
  on.exit(options(old), add = TRUE)
  invisible(rgl::open3d(useNULL = TRUE))
  on.exit(rgl::close3d(), add = TRUE)
  x <- matrix(c(-1, -.5, .25, .5, 1, 2), ncol = 1)
  y <- c(1, 1, 2, 3, 4, 5)
  event <- c(1, 0, 1, 0, 1, 1)
  expect_no_warning(view_dr_surv(x, y, event, B = 1, bw = .8,
                                FUN = identity, type = "3D"))
  ids <- rgl::ids3d()
  vertices <- rgl::rgl.attrib(ids$id[ids$type == "surface"], "vertices")
  xgrid <- seq(min(x), max(x), length.out = 100)
  timegrid <- sort(unique(y))
  expect_equal(nrow(vertices), 100L * length(timegrid))
  expect_equal(unname(vertices[, 1]), rep(xgrid / diff(range(x)), length(timegrid)),
               tolerance = 1e-7)
  expect_equal(unname(vertices[, 2]), rep(timegrid / diff(range(y)) * 1.5, each = 100),
               tolerance = 1e-7)
  survival <- matrix(vertices[, 3], 100)
  for (i in c(1L, 48L, 100L)) {
    weights <- exp(-.5 * ((xgrid[i] - x[, 1]) / .8)^2)
    expect_equal(survival[i, ], view_surv_reference(y, event, weights, timegrid),
                 tolerance = 1e-7)
  }
})

test_that("one-time 3D survival displays its actual profile as a line", {
  old <- options(rgl.useNULL = TRUE)
  on.exit(options(old), add = TRUE)
  invisible(rgl::open3d(useNULL = TRUE))
  on.exit(rgl::close3d(), add = TRUE)
  x <- matrix(c(-1, -.5, .5, 1), ncol = 1)
  y <- rep(2, 4)
  event <- c(1, 1, 0, 0)
  expect_no_warning(view_dr_surv(x, y, event, B = 1, bw = .6,
                                FUN = identity, type = "3D"))
  ids <- rgl::ids3d()
  expect_false(any(ids$type == "surface"))
  profiles <- lapply(ids$id[ids$type == "linestrip"], function(id)
    rgl::rgl.attrib(id, "vertices"))
  profile <- Filter(function(vertices) nrow(vertices) == 100L, profiles)
  expect_length(profile, 1L)
  vertices <- profile[[1]]
  expect_true(all(is.finite(vertices)))
  expect_equal(unname(vertices[, 2]), rep(3, 100))
  xgrid <- seq(min(x), max(x), length.out = 100)
  expected <- vapply(xgrid, function(t) {
    weights <- exp(-.5 * ((t - x[, 1]) / .6)^2)
    1 - sum(weights[event == 1]) / sum(weights)
  }, numeric(1))
  expect_equal(unname(vertices[, 3]), expected, tolerance = 1e-7)
})

test_that("2D survival keeps its observed points and graphics arguments", {
  shown <- NULL
  local_mocked_bindings(
    plot.new = function(...) stop("Unexpected extra graphics page"),
    par = function(...) NULL,
    plot = function(x, y, ...) { shown <<- list(x = x, y = y, args = list(...)) },
    .package = "orthoDr"
  )
  x <- matrix(c(1, 2, 4), ncol = 1)
  expect_null(view_dr_surv(x, c(1, 2, 5), c(1, 0, 1), B = 2, FUN = log,
                          legend.add = FALSE, xlab = "Index", ylab = "Log time"))
  expect_equal(shown$x, x * 2)
  expect_equal(shown$y, log(c(1, 2, 5)))
  expect_equal(shown$args$xlab, "Index")
  expect_equal(shown$args$ylab, "Log time")
  expect_equal(shown$args$col, c("blue", "red", "blue"))
})

test_that("display-time rounding preserves the last failure", {
  old <- options(rgl.useNULL = TRUE)
  on.exit(options(old), add = TRUE)
  invisible(rgl::open3d(useNULL = TRUE))
  on.exit(rgl::close3d(), add = TRUE)
  view_dr_surv(matrix(0:3), 2^52 + 1:4, c(0,0,0,1), B=1, bw=1, type="3D")
  ids <- rgl::ids3d()
  vertices <- lapply(ids$id[ids$type == "linestrip"], function(id) rgl::rgl.attrib(id, "vertices"))
  profiles <- Filter(function(v) nrow(v) == 100, vertices)
  expect_length(profiles, 1)
  expect_equal(unname(profiles[[1]][,3]), rep(0,100))
})
