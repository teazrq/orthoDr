split_projection_float <- function(value) {
  shifted <- (2^27 + 1) * value
  high <- shifted - (shifted - value)
  list(high = high, low = value - high)
}

compensated_projection <- function(x, origin, basis) {
  delta <- sweep(x, 2L, origin, "-")
  back <- x - delta
  residual <- (x - (delta + back)) + sweep(back, 2L, origin, "-")

  # Retain the residuals of the difference, each product and each sum.
  # Scale only large terms before splitting their significands.
  bound <- .Machine$double.xmax / 2^28
  divisor <- ifelse(abs(delta) >= bound, 2^28, 1)
  small <- delta / divisor
  terms <- split_projection_float(small)
  coeff <- split_projection_float(basis)
  out <- matrix(0, nrow(x), ncol(basis))
  for (k in seq_len(ncol(basis))) {
    total <- error <- numeric(nrow(x))
    for (j in seq_len(ncol(x))) {
      product <- small[, j] * basis[j, k]
      hi <- terms$high[, j]
      lo <- terms$low[, j]
      tail <- ((hi * coeff$high[j, k] - product) + hi * coeff$low[j, k] +
               lo * coeff$high[j, k] + lo * coeff$low[j, k])
      factor <- divisor[, j]
      term <- product * factor
      updated <- total + term
      carried <- updated - total
      error <- error + ((total - (updated - carried)) + (term - carried))
      error <- error + (tail * factor + residual[, j] * basis[j, k])
      total <- updated
    }
    out[, k] <- total + error
  }
  out
}

scaled_prediction_distance <- function(delta, peak, shape, bw, divisor = 1) {
  # Combine opposite scales to keep the standardized distance in range.
  denominator <- peak * bw * shape / divisor
  if (is.finite(denominator) && denominator >= .Machine$double.xmin)
    return(delta / denominator)
  log_scale <- log(peak) + log(bw) + log(shape / divisor)
  exp(log(abs(delta)) - log_scale)
}

prediction_log_weights <- function(object, testx) {
  if (is.null(object$bw))
    stop("This fitted object has no bandwidth. Supply the original bandwidth as object$bw or refit the model.")
  bw <- object$bw
  if (length(bw) != 1L || !is.finite(bw) || bw <= 0)
    stop("The fitted bandwidth must be a positive finite scalar")
  x <- object$x
  basis <- object$B
  if (ncol(testx) != nrow(basis))
    stop("New observations must have the same number of predictors as the training data")
  if (any(!is.finite(basis)))
    stop("Training projection scales must be positive and finite for prediction")
  if (is.complex(testx) || any(!is.finite(testx)))
    stop("testx must be a finite numerical matrix")
  storage.mode(x) <- "double"
  storage.mode(testx) <- "double"

  # Anchor before projection so within-feature differences survive translation.
  # Leave headroom for subtraction and the dot product at extreme scales.
  bound <- 2^ceiling(log2(4 * ncol(x)))
  if (max(abs(x)) > .Machine$double.xmax / bound) {
    x <- x / bound
    testx <- testx / bound
  }
  origin <- x[1L, ]
  centered <- sweep(x, 2L, origin, "-")
  u <- centered %*% basis
  peak <- apply(abs(u), 2L, max)
  magnitude <- drop(apply(abs(centered), 2L, max) %*% abs(basis))
  eps <- .Machine$double.eps
  ill_conditioned <- is.finite(peak) & eps * magnitude > sqrt(eps) * peak
  if (any(ill_conditioned)) {
    u[, ill_conditioned] <- compensated_projection(x, origin, basis[, ill_conditioned, drop = FALSE])
    peak <- apply(abs(u), 2L, max)
  }
  if (any(!is.finite(peak) | peak <= 0))
    stop("Training projection scales must be positive and finite for prediction")
  shape <- apply(sweep(u, 2L, peak, "/"), 2L, sd) * sqrt(2)
  if (any(!is.finite(shape) | shape <= 0))
    stop("Training projection scales must be positive and finite for prediction")

  centered_queries <- sweep(testx, 2L, origin, "-")
  v <- centered_queries %*% basis
  query_magnitude <- sweep(abs(centered_queries) %*% abs(basis), 2L, magnitude, "+")
  direct <- any(ill_conditioned) | rowSums(!is.finite(query_magnitude)) > 0L
  log_weights <- matrix(0, nrow(v), nrow(u))
  for (j in seq_len(ncol(u))) {
    # If projected-coordinate rounding can affect distances at the
    # bandwidth, form each query-to-training difference first.
    rounding <- scaled_prediction_distance(eps * query_magnitude[, j], peak[j], shape[j], bw)
    direct <- direct | rounding > sqrt(eps)
    delta <- outer(v[, j], u[, j], "-")
    scaled <- scaled_prediction_distance(delta, peak[j], shape[j], bw)
    log_weights <- log_weights - scaled^2
  }
  for (i in which(direct)) {
    difference <- compensated_projection(x, testx[i, ], basis)
    divisor <- 1
    if (any(!is.finite(difference))) {
      divisor <- bound
      difference <- compensated_projection(x / divisor, testx[i, ] / divisor, basis)
    }
    log_weights[i, ] <- 0
    for (j in seq_len(ncol(u))) {
      scaled <- scaled_prediction_distance(difference[, j], peak[j], shape[j], bw, divisor)
      log_weights[i, ] <- log_weights[i, ] - scaled^2
    }
  }
  if (anyNA(log_weights)) stop("Gaussian prediction log weights are invalid")
  log_weights
}

gaussian_mean <- function(log_weights, y) {
  shifted <- log_weights - max(log_weights)
  weights <- exp(shifted)
  terms <- weights * y
  # A relative weight can round to zero before multiplication by y.
  tiny <- weights < .Machine$double.xmin
  root <- exp(shifted[tiny] / 2)
  terms[tiny] <- (y[tiny] * root) * root
  total <- sum(weights)
  numerator <- sum(terms)
  if (is.finite(numerator)) return(numerator / total)
  scale <- max(abs(y))
  mean <- sum(terms / scale) / total
  max(-1, min(1, mean)) * scale
}

# Rows of log_weights are query profiles. Return the time-by-profile curve
# alongside the sorted observation indices and their time block assignments.
weighted_km <- function(log_weights, times, event) {
  if (anyNA(log_weights) || any(log_weights > 0))
    stop("Invalid weighted survival log weights")
  ord <- order(times)
  sorted <- times[ord]
  event <- as.logical(event[ord])
  new_time <- c(TRUE, sorted[-1L] != sorted[-length(sorted)])
  starts <- which(new_time)
  ends <- c(starts[-1L] - 1L, length(times))
  blocks <- cumsum(new_time)
  n_blocks <- length(starts)
  n_queries <- nrow(log_weights)
  log_factors <- matrix(0, n_blocks, n_queries)

  if (n_queries == 0L)
    return(list(order = ord, blocks = blocks, surv = log_factors))

  log_weights <- t(log_weights[, ord, drop = FALSE])
  peak <- apply(log_weights, 2L, max)
  if (all(is.finite(peak)) && all(apply(log_weights, 2L, min) - peak >= -512)) {
    # With this log range, every positive mass and every finite ratio
    # remains representable on one common scale, even late in follow-up.
    weights <- exp(sweep(log_weights, 2L, peak, "-"))
    deaths <- rowsum(weights * event, blocks, reorder = FALSE)
    alive <- rowsum(weights * !event, blocks, reorder = FALSE)
    if (n_blocks > 1L) {
      later <- apply(deaths + alive, 2L, function(mass) rev(cumsum(rev(mass))))
      alive[-n_blocks, ] <- alive[-n_blocks, , drop = FALSE] + later[-1L, , drop = FALSE]
    }
    failed <- deaths > 0
    log_factors[failed] <- -log1p(deaths[failed] / alive[failed])
    surv <- matrix(exp(apply(log_factors, 2L, cumsum)), n_blocks, n_queries)
    return(list(order = ord, blocks = blocks, surv = surv))
  }

  block_mass <- function(indices) {
    if (length(indices) == 1L) {
      peak <- log_weights[indices, ]
      return(list(peak = peak, mass = as.numeric(is.finite(peak))))
    }
    logs <- log_weights[indices, , drop = FALSE]
    peak <- apply(logs, 2L, max)
    reference <- ifelse(is.finite(peak), peak, 0)
    list(peak = peak, mass = colSums(exp(sweep(logs, 2L, reference, "-"))))
  }

  merge_mass <- function(a, b) {
    peak <- pmax(a$peak, b$peak)
    reference <- ifelse(is.finite(peak), peak, 0)
    mass <- a$mass * exp(a$peak - reference) + b$mass * exp(b$peak - reference)
    list(peak = peak, mass = mass)
  }

  later <- list(peak = rep(-Inf, n_queries), mass = numeric(n_queries))
  for (i in rev(seq_len(n_blocks))) {
    idx <- seq.int(starts[i], ends[i])
    failed <- event[idx]
    if (any(!failed)) later <- merge_mass(later, block_mass(idx[!failed]))
    if (!any(failed)) next
    died <- block_mass(idx[failed])
    if (i == n_blocks && all(failed)) {
      log_factors[i, ] <- -Inf
    } else {
      if (any(!is.finite(pmax(died$peak, later$peak))))
        stop("Gaussian log distances overflowed in a mixed survival risk set; increase bw")
      # Compute the ratio before adding log masses. Adding to a very
      # negative peak first can erase the count of ties.
      ratio <- died$peak - later$peak + log(died$mass) - log(later$mass)
      log_factors[i, ] <- -pmax(ratio, 0) - log1p(exp(-abs(ratio)))
    }
    later <- merge_mass(later, died)
  }
  surv <- matrix(exp(apply(log_factors, 2L, cumsum)), n_blocks, n_queries)
  list(order = ord, blocks = blocks, surv = surv)
}
