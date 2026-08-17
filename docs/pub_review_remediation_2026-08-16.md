# pub_review Remediation Log

*2026-08-16 16:31 PDT*

Remediation pass against
`docs/pub_review_whitepaper_2026-08-16.md`. This log records what
was fixed, what was deferred, and what was newly discovered while
fixing. It does not edit the whitepaper, which remains the review
record of what was found on 2026-08-16.

## 1. Fixed

### Major issues

- **Major Issue 1 (unimplemented pruning algorithm presented as the
  paper's contribution).** `[applied, unverified]`. Files:
  `analysis/report/report.Rmd`. Reframed per the whitepaper's
  Recommended Framing (Section 5, Framing 2): changed the title to
  "Type I Error and Power of Exact and Asymptotic Trend Tests in
  Small Ordered r x 2 Tables: A Simulation Study"; added an Abstract
  that states plainly the pruning algorithm is not implemented or
  benchmarked; relabeled "Connection to the Fisher exact test
  framework" as "(proposed)" and added an explicit disclaimer
  paragraph at the top of that subsection and of "Pruning
  adaptations"; rewrote "The present study" and the relevant
  Discussion paragraph to describe the pruning material as a
  proposed extension, not a completed contribution; updated Future
  Research item 1 to state explicitly what is and is not
  implemented. Not verified by rendering (see Deferred).
- **Major Issue 2 (two-sided p-value definition mismatch).**
  `[verified]`. Files: `analysis/report/report.Rmd`. Chose
  whitepaper remediation option (b): added an "Operational
  two-sided p-value used in this study" paragraph to the Methods
  section stating explicitly that the doubled one-sided
  `CATTexact` p-value is a distinct, conservative construction from
  the mass-based definition, and that the mass-based version was
  not implemented (it requires the unimplemented enumeration
  extension). Verified by installing `CATTexact` (was not
  previously installed in this environment;
  `requireNamespace("CATTexact")` now returns `TRUE`) and reading
  its help page directly (`?CATTexact::catt_exact`), confirming it
  returns only a one-sided p-value, exactly as the whitepaper
  states.
- **Major Issue 3 (power comparison confounded by unequal realized
  size).** `[verified]`. File: `analysis/report/report.Rmd`. Added
  an explicit caveat paragraph after the power comparison stating
  the confound, quoting the actual realized sizes
  (0.043 / 0.033 / 0.015), and stating that the reported power
  figures are raw operating characteristics at nominal, not
  realized, alpha. Verified the realized rates by rerunning the
  full B = 2,000 simulation (see Section 3).
- **Major Issue 4 (stale, self-contradicted ADEMP verdict).**
  `[verified]`. Files: `analysis/report/report.Rmd`. Moved the
  ADEMP compliance note out of the "References" heading into
  Methods (new "### Morris et al. (2019) ADEMP compliance"
  subsection under "Simulation design"), and rewrote it to state
  the current, verified-compliant status, citing the specific line
  numbers/behaviors in `sim_trend_test.R` that satisfy each
  criterion (no `set.seed()` inside the worker; single
  `RNGkind("L'Ecuyer-CMRG")`/`set.seed(20260417)`; per-replication
  `.Random.seed` capture; `mcse_reject`/`mcse_mean_p` reported).
  Verified by direct reading of the current `sim_trend_test.R`
  (confirmed no `set.seed()` call inside the function, confirmed
  `mcse_reject`/`mcse_mean_p` computation) and by rerunning the
  simulation end to end without error.
- **Major Issue 5 (Results prose in conditional/future mood, no
  observed numbers).** `[verified]`. File:
  `analysis/report/report.Rmd`. Rewrote both Results subsections to
  state the actual observed rejection rates using inline `r`
  expressions (`r_null_asymp`, `r_null_fisher`, `r_null_exact`,
  `r_alt_asymp`, `r_alt_fisher`, `r_alt_exact`) computed from
  `tab_null`/`tab_alt` at render time, rather than hardcoded
  literals, so the numbers cannot desync from the simulation output
  again. Verified by re-executing the exact `sim-null`/`sim-alt`
  plus new `fmt_rate()` chunk logic against saved `.rds` output:
  reproduced 0.043 / 0.033 / 0.015 (null) and 0.426 / 0.266 / 0.370
  (alternative), matching the whitepaper's inspected values.
- **Major Issue 6 (undisclosed dropped replications).**
  `[verified]`. Files: `analysis/report/report.Rmd`. Added a table
  caption note and a Results paragraph disclosing that 2 of 2,000
  alternative-scenario replications were excluded for the
  asymptotic and exact trend methods due to a degenerate column
  margin ($c_1 = 0$ or $c_1 = N$) driving the Cochran-Armitage
  variance to zero, and noting the same two replications are
  dropped for both methods (shared cause). Verified by rerunning
  the full simulation and confirming exactly 2 `NA` values in
  `z_ca` and `p_exact_trend` under the alternative scenario, 0
  under the null, matching the whitepaper's inspection.

### Minor issues / checklist items

- **Minor Issue 2 (missing Strasser and Weber 1999 citation).**
  `[applied, unverified]`. Files: `analysis/report/references.bib`,
  `analysis/report/report.Rmd`. Added a `strasser1999asymptotic` bib
  entry and replaced the plain-text "Strasser-Weber (1999)" with
  `@strasser1999asymptotic`. Not verified against a PDF render
  (bibliography formatting was not re-rendered; see Deferred).
- **Minor Issue 3 (orphaned `mehta1980network`).**
  `[applied, unverified]`. File: `analysis/report/report.Rmd`. Cited
  it alongside `mehta1992exact` in "Current standard
  implementations" as a source for the network algorithm, per the
  whitepaper's suggestion. Verified structurally: a repo-wide check
  confirms every `@`-key used in `report.Rmd` now has a bib entry
  and every bib entry is cited.
- **Minor Issue 4 (unsupported SAS/SPSS claim).**
  `[applied, unverified]`. File: `analysis/report/report.Rmd`.
  Softened "has been adopted as the standard backend... in SAS PROC
  FREQ and SPSS Exact Tests" to "is understood to underlie... though
  we have not verified this against the SAS/SPSS documentation
  directly."
- **Minor Issue 5 (undeclared `CATTexact`/`kableExtra`
  dependencies).** `[verified]`. File: `DESCRIPTION`. Added
  `CATTexact` and `kableExtra` to `Suggests:`. Verified `CATTexact`
  is now installable and loadable
  (`requireNamespace("CATTexact")` -> `TRUE`, confirmed by
  installing it and reading its documentation).
- **Minor Issue 6 (placeholder test suite).** `[verified]`. File:
  `inst/tinytest/test-basic.R`. Replaced the single
  `expect_true(TRUE)` with 10 real assertions: (a) `z_ca`
  reproduces an independent, hand-coded re-derivation of the
  Cochran-Armitage statistic and its exact (N-1
  finite-population-corrected) variance on a non-degenerate table;
  (b) the trend statistic sign is correct for a monotone
  increasing dose-response table; (c) the degenerate-table branch
  (`c_1 = 0`) is confirmed structurally; (d)-(g)
  `summarize_sim()`'s `n_valid`, `reject_rate`, `mcse_reject`, and
  `mcse_mean_p` are checked against hand computations on a
  constructed 5-row data set with one `NA`; (h) `mean_p` excludes
  `NA`s correctly; (i) the asymptotic p-value formula
  `2 * pnorm(-abs(z))` is sanity-checked against a literal value.
  Verified by running
  `Rscript -e 'pkgload::load_all("."); tinytest::run_test_dir("inst/tinytest")'`:
  all 10 tests pass (0.1-0.2s).
- **Minor Issue 7 (leftover Palmer Penguins scaffold, orphaned
  `sim_32.rds`).** `[verified]`. Files:
  `analysis/data/README.md`, `analysis/data/derived_data/`.
  Rewrote `analysis/data/README.md` to describe this project's
  actual (simulation-only, no external data) situation instead of
  the unrelated Palmer Penguins pipeline; removed the orphaned,
  never-referenced `sim_32.rds` (`git rm --cached` plus file
  deletion; confirmed zero references via
  `grep -rn "sim_32"` before removal).
- **Minor Issue 8 (no abstract).** `[applied, unverified]`. File:
  `analysis/report/report.Rmd`. Added an unnumbered `# Abstract {-}`
  section stating the problem, approach, and honestly scoped
  findings (including the unimplemented pruning algorithm and the
  size/power confound). Not verified by PDF render.
- **Minor Issue 9 (no data/code availability statement).**
  `[applied, unverified]`. File: `analysis/report/report.Rmd`. Added
  a "## Data and code availability" subsection under Discussion,
  stating the simulation driver and manuscript are in this
  repository, `CATTexact` is the only non-CRAN dependency, the
  companion Fisher-exact-test enumeration repository
  (`01-fisher-exact-rx2`) is separate and does not contain the
  pruning adaptation either, and the RNG seed discipline used.
- **Minor Issue 10 (ADEMP note under "References" heading).**
  `[verified]`. Same edit as Major Issue 4: the subsection was
  physically moved out from under "# References" into Methods.
  Verified structurally by re-reading the file: "# References" now
  contains no subheadings.
- **New: added a "## Limitations" subsection** under Discussion
  covering (1) single scenario pair / limited generalizability, (2)
  the doubling-rule p-value construction, (3) omission of the
  unconditional exact test comparators
  (`@shan2012efficient`, `@consiglio2014comparison`) despite citing
  them as competitors, (4) the unimplemented pruning algorithm, and
  (5) the post hoc (not a priori) justification of B = 2,000
  (addressing (c)-tier Minor Issue 1 as an explicit disclosure
  rather than a rerun at B = 2,500).

## 2. Deferred

- **Full mass-based two-sided exact p-value implementation**
  (Major Issue 2, option (a)). Reason: requires implementing the
  hypergeometric enumeration over `|T(t)| >= |T_obs|` from scratch
  (the "proposed extension" enumeration machinery), which is a
  substantial software task, not a remediation-scale fix. Deferred
  per the whitepaper's own alternative remediation path (option
  (b), which was implemented instead). No command needed; this is
  future implementation work tracked as Future Research item 1.
- **Tree-based pruning algorithm implementation and benchmarking**
  (Major Issue 1, option (a) / "Framing 1" path). Reason: explicitly
  out of budget for a remediation pass (the whitepaper itself
  frames this as requiring a new implementation effort comparable
  to the companion Fisher-exact-test project). Addressed instead via
  reframing (option (b) in the whitepaper), which is a decision the
  author should confirm is acceptable; if Framing 1 is preferred
  instead, the pruning traversal must be implemented and validated
  against `CATTexact` and the network algorithm before resubmission.
- **PDF re-render of `report.Rmd`.** Reason: `kableExtra`'s
  dependency `stringi` failed to compile from source in this host R
  session (`install.packages("kableExtra")` errored with
  "there is no package called 'stringi'", likely missing system ICU
  headers on this host). `bookdown` and `CATTexact` were
  successfully installed and the full R code in the document was
  verified to parse (`knitr::purl()` + `parse()`) and to execute
  correctly when run outside the `kable_styling()` call. Exact
  command for the user to complete the render once dependencies are
  available: `bash tools/render.sh analysis/report/report.Rmd`
  (or `make docker-render REPORT=analysis/report/report.Rmd` for
  the containerized build, which should not hit this host-library
  gap).
- **Size-adjusted (empirical-alpha-matched) power curves** (Major
  Issue 3, full quantitative fix). Reason: correctly re-deriving
  each method's rejection threshold from its own simulated null
  distribution and reporting power at matched empirical alpha is a
  substantive additional analysis, not a small edit, and was judged
  out of scope for a time-boxed remediation pass. Addressed instead
  with an explicit caveat in the text (implemented). Deferred as
  future work; no immediate command, but the analysis would reuse
  `sim_null`/`sim_alt` already saved at
  `analysis/data/derived_data/sim_null.rds` /
  `analysis/data/derived_data/sim_alt.rds`.
- **Unconditional exact trend test comparators**
  (`@shan2012efficient`, `@consiglio2014comparison`) in the
  empirical comparison (per Section 5's "Comparators" note).
  Reason: requires implementing or sourcing a second exact-test
  package/algorithm, out of budget. Disclosed explicitly in the new
  Limitations subsection instead. Tracked as Future Research item 2.
- **`report_smoke.md` at the repository root left unedited.**
  Reason: not the manuscript under review (the whitepaper reviews
  `analysis/report/report.Rmd` specifically), not referenced by any
  Makefile/CI target, and editing it was judged out of scope for a
  remediation pass tied to the whitepaper's specific findings. See
  "New issues found" below; the user should decide whether to
  delete or regenerate it.
- **Empty `analysis/tables/` and `analysis/figures/` directories**
  (Minor Issue 7, remainder). Reason: (c)-tier polish; git does not
  track empty directories, so removing them has no repository-level
  effect, and the whitepaper's higher-value scaffold cleanup
  (Palmer Penguins README, orphaned `sim_32.rds`) was already
  addressed. Left as-is.

## 3. New issues found while fixing

- **`report_smoke.md`** (repository root) is a stale, untracked-by-
  CI duplicate of an earlier `report.Rmd` render (dated
  2026-04-19 in its YAML), carrying the same pre-remediation title,
  the same self-contradicted "Not compliant" ADEMP verdict, and a
  visible R error in its own body
  (``Error in `dimnames(x) <- dn`: length of 'dimnames' [2] not
  equal to array extent``) baked into the rendered markdown from a
  failed `kable()` call. It is not referenced by any workflow or
  Makefile target found in this repository. The user should decide
  whether to delete it, regenerate it, or exclude it from version
  control.
- **`sim_trend_test()`'s asymptotic variance formula uses the
  (N-1) finite-population correction**
  (`var_denom <- n_total^2 * (n_total - 1)`,
  `sim_trend_test.R:31`), which is the standard Armitage/SAS PROC
  FREQ TREND convention, but differs from base R's
  `stats::prop.trend.test()`, which omits this correction
  (confirmed numerically: on a test table, the code's formula gives
  chi-square = 8.571 while `prop.trend.test()` gives 8.791; the two
  converge only asymptotically as N grows). This was investigated
  because the whitepaper's Minor Issue 6 suggested testing `z_ca`
  against `prop.trend.test()`; that direct comparison would have
  produced a spurious test failure, so the new unit tests instead
  validate `z_ca` against an independently re-coded implementation
  of the same (N-1)-corrected formula rather than against
  `prop.trend.test()`. This divergence between two legitimate
  conventions is not itself a bug, but it is undocumented in the
  manuscript's Methods section (no variance formula for the
  asymptotic test is given at all, only the trend statistic $T$).
  The author may wish to add the variance formula and a citation to
  Methods for completeness.
- **`CATTexact` was not installed in the review/remediation R
  environment** prior to this session
  (`requireNamespace("CATTexact")` returned `FALSE`, confirming the
  whitepaper's Minor Issue 5 finding independently). It was
  installed from CRAN during remediation to verify Major Issue 2
  and to rerun the simulation; this local installation is not
  persisted anywhere beyond this host R library and does not affect
  `renv.lock` (already correctly pinned).

## Supporting artifacts from this remediation

- `analysis/data/derived_data/sim_null.rds`,
  `analysis/data/derived_data/sim_alt.rds`: full B = 2,000
  reruns of both scenarios with the manuscript's exact seed
  (`RNGkind("L'Ecuyer-CMRG")`, `set.seed(20260417)`), used to
  verify every numeric claim in this log and now available for the
  manuscript's inline `r` expressions to read from at render time.
