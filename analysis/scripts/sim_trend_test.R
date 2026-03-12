sim_trend_test <- function(n_per_group,
                           true_probs,
                           scores = NULL,
                           n_rep = 2000,
                           alpha = 0.05,
                           seed = 42) {
  set.seed(seed)
  k <- length(true_probs)
  if (is.null(scores)) scores <- seq_len(k)
  stopifnot(length(n_per_group) == k)
  stopifnot(length(scores) == k)

  results <- vector("list", n_rep)

  for (i in seq_len(n_rep)) {
    x <- rbinom(k, size = n_per_group, prob = true_probs)
    dat <- rbind(x, n_per_group - x)

    n_total <- sum(dat)
    r_sums <- colSums(dat)
    c_sums <- rowSums(dat)
    p_hat <- x / n_per_group

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

    p_exact_trend <- tryCatch({
      CATTexact::catt_exact(scores, n_per_group, x)$exact.pvalue
    }, error = function(e) NA_real_)

    results[[i]] <- data.frame(
      rep = i,
      z_ca = z_ca,
      p_asymp = p_asymp,
      p_fisher = p_exact,
      p_exact_trend = p_exact_trend
    )
  }

  out <- do.call(rbind, results)
  out
}

summarize_sim <- function(sim_results, alpha = 0.05) {
  methods <- c("p_asymp", "p_fisher", "p_exact_trend")
  labels <- c(
    "CA Asymptotic",
    "Fisher Exact (unordered)",
    "CA Exact (conditional)"
  )

  summaries <- lapply(seq_along(methods), function(j) {
    p_col <- sim_results[[methods[j]]]
    valid <- !is.na(p_col)
    n_valid <- sum(valid)
    reject <- sum(p_col[valid] < alpha)
    data.frame(
      method = labels[j],
      n_valid = n_valid,
      reject_rate = reject / n_valid,
      mean_p = mean(p_col, na.rm = TRUE),
      median_p = median(p_col, na.rm = TRUE)
    )
  })

  do.call(rbind, summaries)
}
