library(tinytest)

# Tests for the two-sided exact Cochran-Armitage p-values used in the
# manuscript. They guard against the directional error found in the
# 2026-11-17 referee report: CATTexact::catt_exact() returns only the
# upper tail, so doubling it is not a two-sided test.
if (!requireNamespace("CATTexact", quietly = TRUE)) {
  exit_file("CATTexact not installed; skipping")
}

find_pkg_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, "DESCRIPTION"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop("Could not locate package root (DESCRIPTION not found)")
    }
    d <- parent
  }
}
sim_script <- file.path(
  find_pkg_root(), "analysis", "scripts", "sim_trend_test.R"
)
if (!file.exists(sim_script)) {
  exit_file("analysis/ not shipped in the built package; skipping")
}
source(sim_script)

n <- rep(10, 4)
d <- 1:4

## -- catt_exact() is one-directional (upper tail only) --
p_inc <- CATTexact::catt_exact(d, n, c(0, 1, 2, 4))$exact.pvalue
p_dec <- CATTexact::catt_exact(d, n, c(4, 2, 1, 0))$exact.pvalue
expect_true(
  p_inc < 0.05 && p_dec > 0.95,
  info = "catt_exact() gives the upper tail; a decreasing table is ~1"
)

## -- exact_trend_tails() recovers the lower tail by reflection --
tl_inc <- exact_trend_tails(d, n, c(0, 1, 2, 4))
tl_dec <- exact_trend_tails(d, n, c(4, 2, 1, 0))
expect_equal(
  tl_dec[["lower"]], tl_inc[["upper"]], tolerance = 1e-12,
  info = "lower tail of the reversed table equals the upper tail"
)
expect_equal(
  min(2 * min(tl_dec), 1), min(2 * min(tl_inc), 1), tolerance = 1e-12,
  info = "equal-tail p-value is symmetric under group reversal"
)

## -- CATTexact tails agree with direct enumeration on every table --
## -- of a small design (unequal n_i, asymmetric unequal scores, --
## -- so the reflected score vector differs from the original) --
n_small <- c(3, 3, 2, 3)
d_small <- c(0, 1, 2, 5)
lk <- enumerate_trend_tables(n_small, d_small)
lk_ok <- lk[!is.na(lk$p_exact_upper), ]
catt_tails <- t(vapply(seq_len(nrow(lk_ok)), \(i) {
  x <- unlist(lk_ok[i, paste0("x", 1:4)])
  exact_trend_tails(d_small, n_small, x)
}, c(upper = 0, lower = 0)))
expect_equal(
  unname(catt_tails[, "upper"]), lk_ok$p_exact_upper,
  tolerance = 1e-10,
  info = "enumerated upper tail matches catt_exact() on all tables"
)
expect_equal(
  unname(catt_tails[, "lower"]), lk_ok$p_exact_lower,
  tolerance = 1e-10,
  info = "enumerated lower tail matches reflected catt_exact()"
)

## -- CATTexact 0.1.1 is not shift-invariant when a score is 0; --
## -- exact_trend_tails() shifts scores to avoid it --
x_z <- c(0, 2, 0, 0)
p_zero <- CATTexact::catt_exact(c(0, 3, 4, 5), n_small, x_z)$exact.pvalue
p_one <- CATTexact::catt_exact(c(1, 4, 5, 6), n_small, x_z)$exact.pvalue
lk_z <- enumerate_trend_tables(n_small, c(0, 3, 4, 5))
p_enum <- lk_z$p_exact_upper[
  lk_z$x1 == 0 & lk_z$x2 == 2 & lk_z$x3 == 0 & lk_z$x4 == 0
]
expect_equal(
  p_one, p_enum, tolerance = 1e-10,
  info = "shifted scores (1, 4, 5, 6) agree with enumeration"
)
expect_equal(
  exact_trend_tails(c(0, 3, 4, 5), n_small, x_z)[["upper"]], p_enum,
  tolerance = 1e-10,
  info = "exact_trend_tails() is protected against the zero-score bug"
)
if (abs(p_zero - p_enum) < 1e-10) {
  message("CATTexact zero-score discrepancy no longer reproduces")
}

## -- tail identities: p_upper + p_lower = 1 + point probability --
pt <- lk_ok$p_exact_upper + lk_ok$p_exact_lower - 1
expect_true(
  all(pt > 0 & pt <= 1 + 1e-12),
  info = "point probability of the observed T is in (0, 1]"
)
expect_true(
  all(lk_ok$p_exact_midp <= lk_ok$p_exact_trend + 1e-12),
  info = "mid-p two-sided p-value never exceeds the equal-tail one"
)

## -- the doubled upper tail never rejects a decreasing trend --
lk10 <- enumerate_trend_tables(n, d)
rej_dec <- exact_rejection(lk10, n, c(0.35, 0.20, 0.10, 0.05))
expect_true(
  rej_dec[["p_exact_doubled"]] < 1e-3,
  info = "doubled upper tail has ~zero power for a decreasing trend"
)
expect_true(
  rej_dec[["p_exact_trend"]] > 0.3,
  info = "equal-tail test detects the decreasing trend"
)

## -- exact size of every two-sided method is at most alpha --
## -- (conditional exact tests); doubled rule is about alpha / 2 --
rej_null <- exact_rejection(lk10, n, rep(0.2, 4))
expect_true(
  all(rej_null[c("p_exact_trend", "p_exact_mass")] <= 0.05),
  info = "conditional exact tests hold size at p = 0.2"
)
expect_true(
  rej_null[["p_exact_doubled"]] <= 0.025,
  info = "doubled upper tail is a one-sided test at level alpha / 2"
)
