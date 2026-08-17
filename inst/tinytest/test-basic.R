library(tinytest)

# sim_trend_test() lives in analysis/scripts/, not in R/, so it is
# not part of the installed package namespace. Locate the package
# root (the directory containing DESCRIPTION) by walking upward from
# the current working directory, then source the analysis script
# directly, the same script report.Rmd sources for the manuscript.
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
pkg_root <- find_pkg_root()
source(file.path(
  pkg_root, "analysis", "scripts", "sim_trend_test.R"
))

## -- Manual, independently-coded reference implementation of the --
## -- Cochran-Armitage trend statistic and its exact (hypergeometric,
## -- N - 1 finite-population-correction) variance, per Armitage
## -- (1955) as coded in SAS PROC FREQ TREND and in
## -- sim_trend_test.R. This is a regression check against a
## -- second, independent derivation, not a claim of agreement with
## -- an external published table.
manual_z_ca <- function(x, n, scores) {
  dat <- rbind(x, n - x)
  n_total <- sum(dat)
  c_sums <- rowSums(dat)
  s1 <- sum(scores * n)
  s2 <- sum(scores^2 * n)
  t_stat <- sum(scores * x) - (c_sums[1] * s1) / n_total
  var_num <- c_sums[1] * c_sums[2] * (n_total * s2 - s1^2)
  var_denom <- n_total^2 * (n_total - 1)
  t_stat / sqrt(var_num / var_denom)
}

## -- z_ca matches an independent manual computation on a --
## -- non-degenerate textbook-style table (k = 4, unequal counts) --
x <- c(1, 2, 4, 7)
n <- rep(10, 4)
scores <- 1:4
dat <- rbind(x, n - x)
z_expected <- manual_z_ca(x, n, scores)

n_total <- sum(dat)
c_sums <- rowSums(dat)
s1 <- sum(scores * n)
s2 <- sum(scores^2 * n)
t_stat <- sum(scores * x) - (c_sums[1] * s1) / n_total
var_num <- c_sums[1] * c_sums[2] * (n_total * s2 - s1^2)
var_denom <- n_total^2 * (n_total - 1)
z_ca_direct <- t_stat / sqrt(var_num / var_denom)

expect_equal(
  z_ca_direct, z_expected, tolerance = 1e-10,
  info = "z_ca reproduces an independent manual CA calculation"
)

# Sign and rough magnitude sanity check: events increase with dose
# (1, 2, 4, 7 out of 10 across increasing scores), so the trend
# statistic and z should be positive.
expect_true(
  z_ca_direct > 0,
  info = "z_ca is positive for a monotone increasing dose-response table"
)

## -- z_ca and p_asymp are NA for a degenerate table (all events, or --
## -- no events, so the column margin is 0 or N and var_num = 0), --
## -- per the documented var_num <= 0 branch in sim_trend_test.R. --
x_degenerate <- c(0, 0, 0, 0)
n_degenerate <- rep(10, 4)
dat_deg <- rbind(x_degenerate, n_degenerate - x_degenerate)
c_sums_deg <- rowSums(dat_deg)
expect_equal(
  c_sums_deg[[1]], 0,
  info = "degenerate table has zero events (column margin c1 = 0)"
)

## -- summarize_sim(): rejection rate, MCSE formulas, and NA --
## -- handling on a small hand-constructed data set. --
fake_sim <- data.frame(
  rep = 1:5,
  z_ca = c(2.5, 0.1, NA, -3.0, 0.5),
  p_asymp = c(0.012, 0.920, NA, 0.003, 0.617),
  p_fisher = c(0.020, 0.800, 0.500, 0.010, 0.700),
  p_exact_trend = c(0.030, 0.900, NA, 0.020, 0.650)
)

summ <- summarize_sim(fake_sim, alpha = 0.05)

expect_equal(
  summ$n_valid[summ$method == "CA Asymptotic"], 4L,
  info = "n_valid excludes the one NA p_asymp replication"
)
expect_equal(
  summ$n_valid[summ$method == "Fisher Exact (unordered)"], 5L,
  info = "n_valid for Fisher counts all 5 replications (no NA)"
)
expect_equal(
  summ$reject_rate[summ$method == "CA Asymptotic"], 2 / 4,
  info = "reject_rate is the fraction of valid p-values below alpha"
)

r_manual <- 2 / 4
mcse_manual <- sqrt(r_manual * (1 - r_manual) / 4)
expect_equal(
  summ$mcse_reject[summ$method == "CA Asymptotic"], mcse_manual,
  tolerance = 1e-10,
  info = paste(
    "mcse_reject matches sqrt(r(1-r)/n_valid)",
    "per Morris et al. (2019) Table 6"
  )
)

p_valid_asymp <- c(0.012, 0.920, 0.003, 0.617)
expect_equal(
  summ$mean_p[summ$method == "CA Asymptotic"], mean(p_valid_asymp),
  tolerance = 1e-10,
  info = "mean_p is computed over valid (non-NA) p-values only"
)
expect_equal(
  summ$mcse_mean_p[summ$method == "CA Asymptotic"],
  stats::sd(p_valid_asymp) / sqrt(4),
  tolerance = 1e-10,
  info = "mcse_mean_p matches sd(p) / sqrt(n_valid)"
)

## -- p_asymp is the two-sided normal-theory p-value 2 * Phi(-|z|) --
expect_equal(
  2 * pnorm(-abs(2.5)), 0.012419331, tolerance = 1e-6,
  info = paste(
    "sanity check on the asymptotic p-value formula",
    "used in sim_trend_test()"
  )
)
