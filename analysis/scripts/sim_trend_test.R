#' Upper and lower exact conditional Cochran-Armitage tail p-values
#'
#' `CATTexact::catt_exact()` returns only the upper-tail p-value
#' P(T >= t_obs), that is, a one-sided test against an increasing
#' trend. The lower tail is obtained from the same function by
#' reflecting the scores (d' = max(d) + min(d) - d) and reversing the
#' group order, so that the reflected scores are again increasing.
#' Because the total event count is fixed under conditioning,
#' sum(d' x) = const - sum(d x), so the upper tail under d' is the
#' lower tail under d. This holds for any increasing score vector and
#' any group sizes.
#'
#' Both score vectors are shifted so that their minimum is 1. The
#' conditional test is invariant to such a shift, but CATTexact 0.1.1
#' returns incorrect p-values for some score vectors that contain 0
#' (e.g. scores (0, 3, 4, 5) disagree with (1, 4, 5, 6) and with
#' direct enumeration; found 2026-11-17).
#'
#' @param scores Strictly increasing group scores.
#' @param n_per_group Group sizes.
#' @param x Event counts per group.
#' @return Named numeric vector with elements `upper` and `lower`.
exact_trend_tails <- function(scores, n_per_group, x) {
  scores <- scores - min(scores) + 1
  scores_rev <- rev(max(scores) + min(scores) - scores)
  c(
    upper = CATTexact::catt_exact(scores, n_per_group, x)$exact.pvalue,
    lower = CATTexact::catt_exact(
      scores_rev, rev(n_per_group), rev(x)
    )$exact.pvalue
  )
}

#' Enumerate every table of a k x 2 design with all p-values
#'
#' Builds the full product-binomial sample space (all event-count
#' vectors x with 0 <= x_i <= n_i) and, for each table, computes the
#' p-values of every method compared in the manuscript. Exact
#' conditional p-values are obtained by direct enumeration of the
#' conditional (multivariate hypergeometric) distribution of
#' T = sum(d x) given the total event count R1. This is an
#' implementation independent of `CATTexact`, used both as a lookup
#' in the simulation and for exact (non-Monte Carlo) size and power.
#'
#' Tables with a degenerate margin (R1 = 0 or R1 = N) have zero
#' Cochran-Armitage variance; the asymptotic and all exact trend
#' p-values are set to NA there, as `catt_exact()` errors on them.
#'
#' @param n_per_group Group sizes.
#' @param scores Strictly increasing group scores.
#' @param tol Tolerance for ties in the trend statistic.
#' @return Data frame with one row per table: counts `x1..xk`, `r1`,
#'   and p-value columns `p_asymp`, `p_asymp_upper` (one-sided,
#'   increasing), `p_fisher`, `p_exact_upper`,
#'   `p_exact_lower`, `p_exact_trend` (equal-tail), `p_exact_midp`,
#'   `p_exact_mass`, and `p_exact_doubled`.
enumerate_trend_tables <- function(n_per_group, scores, tol = 1e-7) {
  k <- length(n_per_group)
  stopifnot(length(scores) == k, all(diff(scores) > 0))
  n_tables <- prod(n_per_group + 1)
  if (n_tables > 2e5) {
    stop(
      "Sample space has ", n_tables, " tables; brute-force ",
      "enumeration is limited to 2e5. Supply `lookup` or use a ",
      "smaller design."
    )
  }
  grid <- as.matrix(expand.grid(
    lapply(n_per_group, function(m) 0:m)
  ))
  colnames(grid) <- paste0("x", seq_len(k))
  n_total <- sum(n_per_group)
  r1 <- rowSums(grid)
  t_obs <- as.vector(grid %*% scores)
  log_w <- rowSums(lchoose(
    matrix(n_per_group, nrow(grid), k, byrow = TRUE), grid
  ))
  s1 <- sum(scores * n_per_group)
  s2 <- sum(scores^2 * n_per_group)
  e_t <- r1 * s1 / n_total
  var_t <- r1 * (n_total - r1) * (n_total * s2 - s1^2) /
    (n_total^2 * (n_total - 1))

  p_up <- p_lo <- p_pt <- p_mass <- rep(NA_real_, nrow(grid))
  for (r in unique(r1)) {
    idx <- which(r1 == r)
    if (r == 0 || r == n_total) next
    pr <- exp(log_w[idx] - lchoose(n_total, r))
    tt <- t_obs[idx]
    dev <- abs(tt - e_t[idx])
    # Base vapply rather than purrr: purrr is not in renv.lock.
    p_up[idx] <- vapply(tt, \(t0) sum(pr[tt >= t0 - tol]), 0)
    p_lo[idx] <- vapply(tt, \(t0) sum(pr[tt <= t0 + tol]), 0)
    p_pt[idx] <- vapply(tt, \(t0) sum(pr[abs(tt - t0) < tol]), 0)
    p_mass[idx] <- vapply(dev, \(d0) sum(pr[dev >= d0 - tol]), 0)
  }

  z <- ifelse(var_t > 0, (t_obs - e_t) / sqrt(var_t), NA_real_)
  p_fisher <- vapply(seq_len(nrow(grid)), \(i) {
    fisher.test(rbind(grid[i, ], n_per_group - grid[i, ]))$p.value
  }, 0)

  data.frame(
    grid,
    r1 = r1,
    p_asymp = 2 * pnorm(-abs(z)),
    p_asymp_upper = pnorm(-z),
    p_fisher = p_fisher,
    p_exact_upper = p_up,
    p_exact_lower = p_lo,
    p_exact_trend = pmin(2 * pmin(p_up, p_lo), 1),
    p_exact_midp = pmin(2 * pmin(p_up, p_lo) - p_pt, 1),
    p_exact_mass = pmin(p_mass, 1),
    p_exact_doubled = pmin(2 * p_up, 1)
  )
}

#' Exact rejection probability of each method under a given DGM
#'
#' Sums the product-binomial probability of every table on which a
#' method rejects. NA p-values (degenerate margins) count as
#' non-rejections. No Monte Carlo error is involved.
#'
#' @param tables Output of `enumerate_trend_tables()`.
#' @param n_per_group Group sizes (as used to build `tables`).
#' @param true_probs Event probability in each group.
#' @param alpha Nominal level.
#' @return Named numeric vector of rejection probabilities.
exact_rejection <- function(tables, n_per_group, true_probs,
                            alpha = 0.05) {
  k <- length(n_per_group)
  x <- as.matrix(tables[, paste0("x", seq_len(k))])
  log_p <- rowSums(sapply(seq_len(k), \(i) {
    dbinom(x[, i], n_per_group[i], true_probs[i], log = TRUE)
  }))
  w <- exp(log_p)
  p_cols <- grep("^p_", names(tables), value = TRUE)
  vapply(
    setNames(p_cols, p_cols),
    \(m) sum(w[!is.na(tables[[m]]) & tables[[m]] < alpha]),
    0
  )
}

sim_trend_test <- function(n_per_group,
                           true_probs,
                           scores = NULL,
                           n_rep = 2000,
                           alpha = 0.05,
                           lookup = NULL) {
  # Morris et al. (2019) §4.1: the RNG seed is set ONCE by the caller;
  # this function does not call set.seed() and must not do so.
  k <- length(true_probs)
  if (is.null(scores)) scores <- seq_len(k)
  stopifnot(length(n_per_group) == k)
  stopifnot(length(scores) == k)
  if (is.null(lookup)) {
    lookup <- enumerate_trend_tables(n_per_group, scores)
  }
  key_lookup <- do.call(
    paste, c(lookup[paste0("x", seq_len(k))], sep = "-")
  )

  results <- vector("list", n_rep)
  rng_states <- vector("list", n_rep)

  for (i in seq_len(n_rep)) {
    rng_states[[i]] <- .Random.seed
    x <- rbinom(k, size = n_per_group, prob = true_probs)
    dat <- rbind(x, n_per_group - x)

    n_total <- sum(dat)
    c_sums <- rowSums(dat)

    t_stat <- sum(scores * x) -
      (c_sums[1] * sum(scores * n_per_group)) / n_total
    var_num <- c_sums[1] * c_sums[2] *
      (n_total * sum(scores^2 * n_per_group) -
         sum(scores * n_per_group)^2)
    var_denom <- n_total^2 * (n_total - 1)
    z_ca <- if (var_num > 0) {
      t_stat / sqrt(var_num / var_denom)
    } else {
      NA_real_
    }

    p_asymp <- if (!is.na(z_ca)) {
      2 * pnorm(-abs(z_ca))
    } else {
      NA_real_
    }

    p_exact <- tryCatch({
      ft <- fisher.test(dat)
      ft$p.value
    }, error = function(e) NA_real_)

    # catt_exact() errors on a degenerate margin (R1 = 0 or N); only
    # that case is mapped to NA. Both tails come from CATTexact.
    tails <- tryCatch(
      exact_trend_tails(scores, n_per_group, x),
      error = function(e) c(upper = NA_real_, lower = NA_real_)
    )
    row_lk <- lookup[match(paste(x, collapse = "-"), key_lookup), ]

    results[[i]] <- data.frame(
      rep = i,
      z_ca = z_ca,
      p_asymp = p_asymp,
      p_fisher = p_exact,
      p_exact_upper = tails[["upper"]],
      p_exact_lower = tails[["lower"]],
      p_exact_trend = min(2 * min(tails), 1),
      p_exact_midp = row_lk$p_exact_midp,
      p_exact_mass = row_lk$p_exact_mass,
      p_exact_doubled = min(2 * tails[["upper"]], 1)
    )
  }

  out <- do.call(rbind, results)
  attr(out, "rng_states") <- rng_states
  out
}

method_labels <- c(
  p_asymp = "CA Asymptotic",
  p_fisher = "Fisher Exact (unordered)",
  p_exact_trend = "CA Exact, equal-tail",
  p_exact_mass = "CA Exact, mass-based",
  p_exact_midp = "CA Exact, mid-p",
  p_exact_doubled = "CA Exact, doubled upper tail"
)

summarize_sim <- function(sim_results, alpha = 0.05) {
  methods <- intersect(names(method_labels), names(sim_results))

  # Monte Carlo SE formulas from Morris, White & Crowther (2019), Table 6.
  summaries <- lapply(methods, function(m) {
    p_col <- sim_results[[m]]
    valid <- !is.na(p_col)
    n_valid <- sum(valid)
    p_valid <- p_col[valid]
    reject <- sum(p_valid < alpha)
    r <- reject / n_valid
    mp <- mean(p_valid)
    data.frame(
      method = method_labels[[m]],
      n_valid = n_valid,
      reject_rate = r,
      mcse_reject = sqrt(r * (1 - r) / n_valid),
      mean_p = mp,
      mcse_mean_p = stats::sd(p_valid) / sqrt(n_valid),
      median_p = stats::median(p_valid)
    )
  })

  do.call(rbind, summaries)
}
