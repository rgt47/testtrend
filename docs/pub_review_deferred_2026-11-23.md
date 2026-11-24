# pub_review Deferred-Item Pass: 32-test-for-trend
*2026-11-23 12:58 PST*

Source log: `docs/pub_review_remediation_2026-08-16.md` (answering
`docs/pub_review_whitepaper_2026-08-16.md`). It is the only
`pub_review_remediation_*.md` in the workspace. Its seven Deferred
bullets are the inventory below.

State of the tree at the start of this pass: HEAD is `af32e81`
(2026-09-11). The working tree holds substantial uncommitted work
from a 2026-11-17 referee and voice revision, recorded in
`analysis/report/referee-remediation-log-2026-11-17.md`. That work
is evidence for the status of several items. Where an item is
resolved by it, the fix exists in the working tree only, not in
any commit. This pass made no edits to the manuscript, code, tests
or data. The only file written is this log.

## 1. Inventory

| # | Item (2026-08-16 Deferred) | Class |
|---|---|---|
| 1 | Mass-based two-sided exact p-value (Major 2, option a) | RESOLVED-SINCE |
| 2 | Tree-based pruning algorithm, implementation and benchmark (Major 1, Framing 1) | AUTHOR-DECISION |
| 3 | PDF re-render of `report.Rmd` | RESOLVED-SINCE |
| 4 | Size-adjusted power (Major 3, quantitative fix) | AUTHOR-DECISION |
| 5 | Unconditional exact trend comparators (Shan 2012, Consiglio 2014) | AUTHOR-DECISION |
| 6 | Stale `report_smoke.md` at the repository root | AUTHOR-DECISION |
| 7 | Empty `analysis/tables/` and `analysis/figures/` | NO-LONGER-APPLIES |

Evidence for each classification:

- **Item 1.** `enumerate_trend_tables()` in
  `analysis/scripts/sim_trend_test.R` now computes `p_exact_mass` by
  enumerating the conditional distribution of the centered
  statistic. It is defined in Methods 2.1 of `report.Rmd` and
  reported in the tables. Working tree only (2026-11-17, M1/M6).
  Not cross-checked against `coin` (see Section 5).
- **Item 2.** Sections 'Connection to the Fisher exact test
  framework (proposed)' and 'Pruning adaptations' remain
  explicitly labeled as unimplemented. The 2026-11-17 pass
  completed the sketch (M7) but did not implement it.
- **Item 3.** `report.pdf` (2026-11-17 18:24) is newer than
  `report.Rmd` (18:23). It was rendered through `tools/render.sh`
  and staged as `share/report-2026-11-17-1824-af32e81-wip.pdf`.
  Caveat carried forward: a plain host render fails unless
  CATTexact is on the library path (it lives only in
  `/opt/homebrew/lib/R/4.6/site-library`). The containerized
  render (`make docker-render REPORT=analysis/report/report.Rmd`)
  has never been run.
- **Item 4.** `report.Rmd` lines 932-937 now decline size-adjusted
  power on stated grounds. It would need a randomized or
  recalibrated rule for each method. This is a later, deliberate
  choice, so it is the author's call. Evidence for that call is in
  Section 4.
- **Item 5.** Still absent. The Limitations section and Future
  research item 3 disclose this. See Section 4.
- **Item 6.** `report_smoke.md` is unchanged (dated 2026-04-19) and
  is tracked in git. It still carries the old title, the stale ADEMP
  verdict, and a baked-in `kable()` error.
- **Item 7.** Both directories are still empty. Git does not track
  empty directories, so this has no effect on the repository. No
  action is needed.

## 2. Done now

No deferred item was classified DO-NOW. Items 4 and 5 are now cheap
to compute with the exact-enumeration code added on 2026-11-17. Each
still changes what the paper claims, so the decision belongs to the
author. Two evidence-only runs were made to support those
decisions. Neither changed any workspace file.

- **Test suite re-run.** `[verified]` All 23 assertions pass: 10 in
  `test-basic.R` and 13 in `test-exact-trend.R`, in 5.4 s under
  R 4.6.2. The Homebrew library was appended last to
  `.libPaths()` so that CATTexact loads. Without it,
  `test-exact-trend.R` exits silently and only 10 assertions run.
  Script:
  `/private/tmp/claude-501/-Users-zenn-prj-res-01-fisher-exact-rx2-fisherexacttestrx2-analysis-report/be4b228f-9df5-41a4-9039-7499d4f23926/scratchpad/deferred/tt32_run_tests.R`.
- **Size-calibrated exact power (evidence for item 4).**
  `[verified]` For the main design (k = 4, n_i = 10, scores 1:4),
  `enumerate_trend_tables()` was run, with `exact_rejection()`-style
  sums at the manuscript's increasing alternative
  (0.05, 0.10, 0.20, 0.35). Script:
  `/private/tmp/claude-501/-Users-zenn-prj-res-01-fisher-exact-rx2-fisherexacttestrx2-analysis-report/be4b228f-9df5-41a4-9039-7499d4f23926/scratchpad/deferred/tt32_size_adj.R`.
  The nominal-alpha sizes at p = 0.2 reproduce the 2026-11-17 log:
  asymptotic 0.0434, Fisher 0.0398, equal-tail 0.0297 and mid-p
  0.0434.

  Each method's cutoff was raised to the largest attainable p-value
  whose exact size stays at or below 0.05. The size was taken as
  the supremum over the common null p in {0.02, ..., 0.50}, the
  same grid as the manuscript's size curve. The design is symmetric
  under p to 1 - p, so this grid covers [0, 1]. Results:

  | Method | Sup size, nominal | Cutoff | Sup size, calibrated | Power, calibrated |
  |---|---|---|---|---|
  | Asymptotic | 0.0514 | 0.0460 | 0.0467 | 0.4429 |
  | Fisher (unordered) | 0.0435 | 0.0532 | 0.0485 | 0.2842 |
  | Exact, equal-tail | 0.0396 | 0.0668 | 0.0497 | 0.4429 |
  | Exact, mid-p | 0.0467 | 0.0503 | 0.0497 | 0.4429 |
  | Exact, mass-based | 0.0396 | 0.0668 | 0.0497 | 0.4429 |

  Once each test is calibrated to the same guaranteed size, the
  asymptotic, equal-tail, mid-p and mass-based tests have the same
  power to four decimals. At nominal alpha their exact powers
  differ: 0.4429 for asymptotic and mid-p against 0.3825 for
  equal-tail and mass-based. That gap is therefore entirely an
  artifact of discreteness and conservatism, not a difference in
  the ordering of the sample space. Fisher stays far behind
  (0.284).

  The run also calibrated to the known p = 0.2 only. That is an
  oracle, not a valid test, and is not recommended. It gives
  asymptotic and mid-p 0.474, equal-tail and mass-based 0.454, and
  Fisher 0.307.

## 3. Long runs

None. Every computation relevant to the deferred items uses the
11^4 = 14,641-table enumeration and runs in seconds (enumeration
took 3.2 to 3.9 s here). The full B = 2,000 Monte Carlo rerun is
already done (2026-11-17) and is not needed by any open item.

## 4. Author decisions

- **Item 2: implement the pruning algorithm, or keep Framing 2.**
  Options: (a) implement the tree traversal and benchmark it against
  `CATTexact` and the network algorithm; (b) keep the current
  framing, with the pruning material labeled as a proposed
  extension. The current manuscript is consistently written for
  (b), and section 1.7 says the study adds no new software.
  Recommendation: keep (b) for this paper. If the author wants (a),
  it is a separate paper, built in the `01-fisher-exact-rx2`
  workspace where the rx2 kernels live.
- **Item 4: report calibrated power, or keep declining.** Options:
  (a) keep the current text, which reports power at nominal alpha
  with a caveat; (b) add a calibrated-power row or paragraph. A
  deterministic cutoff calibrated to the supremum size over the
  nuisance grid is valid and needs no randomization. Section 2
  shows it is cheap. The evidence matters for the paper's message:
  under calibration the four trend tests tie (0.4429), so the
  current 'asymptotic has the highest power' statement holds only
  at nominal alpha. Recommendation: (b). Add one sentence and the
  numbers through an inline-R chunk that reuses `lookup`. Then
  replace the 'we do not report size-adjusted power' sentence in
  lines 932-937. Its stated obstacle, a recalibrated rule, is cheap
  here. This needs a re-render.
- **Item 5: add the unconditional exact test.** Options: (a) add it.
  On the same 14,641-table space, take p = sup over p of
  P_p(|Z| >= |z_obs|), with Z the Cochran-Armitage statistic
  (Tang et al. 2000), over a fine p grid. Then compute its exact
  size and power with `exact_rejection()`. The cost is about
  14,641 x grid-size probabilities, seconds of compute. (b) Keep
  the disclosure in Limitations. Against (a): no working
  independent reference implementation is available here to
  validate it (`coin` segfaults; see Section 5), and it adds a
  comparator to a paper whose abstract now limits its scope.
  Recommendation: (a) only if the author wants the
  conditional-versus-unconditional question answered in this
  paper. It is the natural next comparator, and the Introduction
  spends a subsection on it. Otherwise (b).
- **Item 6: `report_smoke.md`.** Options: delete it with
  `git rm report_smoke.md`, or regenerate it. It is tracked,
  stale, and contains an R error. Recommendation: delete. Nothing
  references it. This pass does not delete files.

## 5. New issues found

- **Orphaned derived data.** `analysis/data/derived_data/sim_null.rds`
  and `sim_alt.rds` (tracked, 2026-08-16) are not read by
  `report.Rmd`. The report regenerates the simulations in chunks
  (`grep readRDS` finds nothing). They also predate the 2026-11-17
  columns (`p_exact_lower`, `p_exact_midp`, `p_exact_mass`), and
  there is no decreasing scenario. The 2026-08-16 log's statement
  that they are 'available for the manuscript's inline `r`
  expressions' is no longer accurate. The author should either
  regenerate them from the current chunks and read them in the
  report, or remove them.
- **Tracked render debris at the repository root.** `cache/` (25
  knitr cache files) and `figure/` (3 PNGs), both from April 2026,
  are tracked in git. They are the same kind of artifact that
  commit `af32e81` removed and ignored under `analysis/report/`.
  The root-level copies were missed. Suggested: `git rm -r --cached
  cache figure` and add both to `.gitignore`. Not done, because
  it changes the git index.
- **Still open from 2026-08-16 'New issues' (not a Deferred item).**
  The asymptotic test's variance formula, with its (N - 1)
  finite-population factor, is still not stated in the manuscript.
  Section 1.2 says only 'standardized by its null variance'. The
  code (`var_t` in `enumerate_trend_tables()`, `var_denom` in
  `sim_trend_test()`) uses the factor; `stats::prop.trend.test()`
  does not. Adding the display to Methods 2.1 is a one-paragraph
  edit plus a re-render.
- **Leftovers of the 2026-11-17 pass (context, not counted
  above).** (i) M5(b): other k and unequal score spacing. This is
  cheap with `enumerate_trend_tables()` within its 2e5-table guard.
  (ii) The `coin` cross-check of `p_exact_mass`: re-attempted today
  under R 4.6.2, with the Homebrew library last. `coin::lbl_test()`
  still segfaults in `R_ExpectationCovarianceStatistic`
  (`[verified]`, script `tt32_coin.R` in the scratchpad above). It
  needs a `coin` built against the native R, or the Docker image.
  (iii) m16: the persistent identifier is still a TODO comment in
  'Data and code availability'.
- **The working tree is large and uncommitted.** About 1,900 lines
  of diff, a new test file, two new logs and a new staged PDF
  (2026-11-17). All the RESOLVED-SINCE classifications above depend
  on that uncommitted state. The author should review and commit
  it before any further pass.
