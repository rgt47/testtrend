# Referee Remediation Log: report.Rmd
*2026-11-17 17:57 PST*

Input: `referee-report-2026-11-17.md` (same directory). Changes are
in the working tree only. Nothing was staged, committed or pushed.
The manuscript was re-rendered through `tools/render.sh`. The uncommitted voice-revision edits
(`voice-revision-log-2026-11-17.md`) were kept and built on. New and
rewritten prose follows the 'Short, simple sentences' rule. The
abstract was rewritten because the correctness fix changed its
findings.

Files changed:

- `analysis/report/report.Rmd`
- `analysis/report/references.bib`
- `analysis/scripts/sim_trend_test.R`
- `inst/tinytest/test-exact-trend.R` (new)
- `analysis/report/referee-remediation-log-2026-11-17.md` (this file)
- render outputs: `report.pdf`, `report.tex`, `report_files/`
  figures (including the new `size-curve-1.pdf`),
  `share/MANIFEST.md`, and
  `share/report-2026-11-17-1824-af32e81-wip.pdf`

Evidence labels: **[verified]** means I ran code or tests that
confirm the fix. **[applied, unverified]** means the edit was made
but not executed or rendered.

## 1. Fixed

### Major comments

- **M1, doubled upper tail is not two-sided.** `report.Rmd`,
  `sim_trend_test.R`. Added `exact_trend_tails()`. It gets the lower
  tail from `catt_exact()` by reflecting the scores and reversing
  the group order. The primary exact p-value is now equal-tail,
  2 min(p_up, p_lo). The doubled upper tail is kept only as a
  labeled comparison row. Methods, Results 3.1, Discussion, Figure 3
  and the abstract were rewritten. A decreasing alternative was
  added. [verified] Full B = 2,000 rerun with the original seed and
  RNG stream reproduces the original numbers exactly: null 0.0425,
  0.0330, doubled 0.0150; alternative 0.426, 0.266, 0.370. The new
  equal-tail null rate is 0.0270, which matches the referee's
  replay. Under the decreasing alternative the doubled rule rejects
  0.000 and the equal-tail test 0.378. The direction of
  `catt_exact()` was confirmed: (0,1,2,4) gives 0.012 and (4,2,1,0)
  gives 0.996.
- **M2, uncentered T(t).** `report.Rmd` Methods 2.1. T(t) is now
  centered, with E_0 T = 0. The pruning bounds in 2.3 use the same
  centering. [applied, unverified] (display math; not rendered)
- **M3, 'refinement of the reference set' argument.** Deleted from
  Results 3.1 and the Discussion. It is replaced by the
  granularity explanation the referee asked for. At R1 = 8 the
  attainable equal-tail p-values near alpha are 0.007, 0.020 and
  0.052. As a result the conditional size is 0.020 against Fisher's
  0.048. [verified] (computed by `enumerate_trend_tables()`; the
  values enter the text through inline R)
- **M4, false MCSE claim.** The Simulation size paragraph and
  Limitations item 5 now give the realized MCSE, close to 0.01
  under the alternatives. They no longer claim 'well below 1 pp'.
  They also note that the exact computations do not depend on B.
  B was not increased, because that would shift the RNG stream for
  the later scenarios. [verified] (MCSE values from the rerun:
  0.0099 to 0.0111)
- **M5, design too narrow. PARTIAL.** (a) Exact size over the
  common null probability p = 0.02 to 0.50 is computed by full
  enumeration of the 11^4 sample space (`exact_rejection()`). It is
  shown as a new figure (`size-curve`). (c) A decreasing alternative
  was added, both by Monte Carlo and exactly. A new table
  (`table-exact`) gives exact size and power for all methods. (b)
  is partly done, by exact computation only (`exact-extra` table).
  It adds a threshold shape (0.10, 0.10, 0.10, 0.40), a plateau
  shape (0.05, 0.25, 0.25, 0.25), and an unequal design
  n = (5, 10, 10, 15) under the null and both linear alternatives.
  The new results change the picture in two ways:
  - Under both non-linear shapes, the unordered Fisher test beats
    the equal-tail exact trend test (threshold 0.364 vs 0.331;
    plateau 0.121 vs 0.098).
  - In the unequal design, the mass-based and equal-tail tests
    differ (null size 0.038 vs 0.031; decreasing power 0.400 vs
    0.312). The asymptotic null size is 0.0506.

  Other k and unequal score spacing are still deferred (Section 2).
  [verified] (exact values computed; for example, exact null sizes
  at p = 0.2 are asymptotic 0.0434, Fisher 0.0398, equal-tail
  0.0297, mid-p 0.0434 and doubled 0.0148)
- **M6, missing comparators. PARTIAL.** The mid-p test (equal-tail,
  Lancaster) and the mass-based exact test (|T - E_0 T| ordering,
  by enumeration) were added. The mass-based test uses the same
  ordering as `coin`'s exact linear-by-linear test. I tried to
  cross-check it against `coin`, but `coin`'s compiled code
  segfaulted. The unconditional exact test is deferred. [verified]
  (In the symmetric main design, the mass-based and equal-tail
  p-values agree to 4e-15 on every table, as the manuscript says.
  They differ in the unequal design.)
- **M7, incomplete pruning sketch.** Section 2.3 now states the box
  constraints, the integer constraint, the assumption that scores
  are sorted and the greedy fill order. It adds the complementary
  acceptance case, with a closed-form subtree probability from the
  Vandermonde identity. It remains labeled as unimplemented.
  [applied, unverified] (mathematics checked by hand only)
- **M8, power interpretation.** Results 3.2 and Discussion point 3
  now state that the asymptotic test has the highest power, and
  that its size is close to nominal in this design. Section 1.3 is
  qualified: the case for exact tests is about guaranteed size,
  not an observed failure here. [verified] (numbers come from inline
  R; exact maximum asymptotic size is 0.0514, at p = 0.5)
- **Related: reconciling with Corcoran and Mehta (2002).** I read
  the full PDF. The paper studies this same design (four groups of
  10), but with one-sided tests. The manuscript now says so. An
  inline computation shows that the one-sided asymptotic test in
  this design exceeds alpha for every p >= 0.13 on the grid, with a
  maximum of 0.052. This is consistent with their finding. The idea
  that tail errors partly cancel in the two-sided test is
  labeled as a suggestion, not a result. [verified]

### Minor comments

- **m1** Displays in 1.2 and 2.1 now end with a period before
  'Here ...'. [applied, unverified]
- **m2** The regulatory claim was narrowed to trend tests in rodent
  carcinogenicity analyses, citing the FDA (2001) draft guidance
  (docket and title confirmed by web search). [applied, unverified]
- **m3** 'by analogy, to general r x c tables' deleted. [applied]
- **m4** 'He' changed to 'They' for @agresti1990exact and
  @mehta1992exact. [applied]
- **m5** `@bartholomew1959test` moved to the mention of
  Bartholomew's statistic. [applied]
- **m6** `perm` cited as @fay2010exact (JSS 36(2), DOI confirmed by
  web search). [applied]
- **m7** `CATTexact` wording softened to 'following an idea of'
  Mehta et al. (1992), matching its DESCRIPTION. [verified against
  the installed DESCRIPTION]
- **m8** The unverified SAS/SPSS algorithm claim was removed. SAS
  `EXACT TREND` is stated, since it is documented (confirmed by web
  search). The algorithm SAS uses is stated as unchecked. SPSS was
  dropped. [applied]
- **m9** 'We have shown' changed to describe the companion work as
  unpublished benchmarks in a separate repository. No citable
  source exists. [applied]
- **m10** Groups are now rows of the k x 2 table throughout. Section
  1.5 was changed to 'row totals'. Methods notes that the code
  stores the 2 x k transpose. [applied]
- **m11** 'per Morris Table 6' changed to '@morris2019using,
  Table 6'. [applied]
- **m12** The ADEMP estimands are now the rejection probabilities
  (size and power). [applied]
- **m13** The change-log-style compliance subsection was replaced by
  a short 'Reproducibility of the simulation' paragraph that points
  to `docs/morris-audit-2026-04-17.md`. [applied]
- **m14** The Table 2 caption now says 'the same' replications. A
  `stopifnot()` checks that the NA rows coincide for every trend
  method. [verified in the chunk run; see Section 4]
- **m15** Figure 3 now uses the two-sided equal-tail p-value, log
  scales and alpha reference lines. [applied, unverified] (not
  rendered)

### Missing and questionable references

- `lin2006exact`: authors corrected to Tang, Ng, Guo, Chan and Chan,
  and pages corrected to 847-859 [verified via Crossref]. The key
  was kept, to preserve citation keys. Its misattributed claim
  ('asymptotic test can be anticonservative') was removed from 1.3.
  It is now cited for what its title supports: a comparison of
  exact conditional and unconditional CA tests under several
  models. Its citation for score-function comparisons in Future
  research was removed, because I could not confirm it.
- `williams1988tests`: DOI corrected to 10.2307/2347316 [verified
  via Crossref]. The claim was rewritten to match the paper as
  found by search: continuous approximations are unreliable for
  several small proportions, and exact conditional distributions
  with mid-p are favored. The old claim, 'proposed an exact
  conditional trend alternative', was dropped.
- `agresti1998order`: volume, issue and pages corrected to 28(2):
  139-155 [verified via Crossref].
- `edelmann2025catexact`: year corrected to 2020 [verified from the
  installed DESCRIPTION, Date/Publication 2020-06-25]. The key was
  kept.
- Added `corcoran2002exact` (JMASM 1(1), article 7; DOI
  10.22237/jmasm/1020254880). This was confirmed from the journal
  page and the full PDF. Pages were not found, so none are given.
- `hirji2006exact` is now also cited for mid-p as a remedy for
  discreteness, and for mid-p possibly exceeding alpha. These
  attributions are [unverified]: I did not check the book's text.
- Added `tang2000exact` (Biom J 42(7):795-806). This was confirmed by
  web search.
- Added `lancaster1961significance` (JASA 56(294):223-234). This was
  confirmed by web search from several reference lists.
- Added `fay2010exact` and `fda2001statistical`, as above.

### Suggestions for strengthening

- **S1/S2** Exact size and power by enumeration, and a size curve
  over the nuisance probability. Done (see M5). [verified]
- **S3** Mid-p and mass-based comparators added. `coin` and the
  unconditional test were not added (see Section 2).
- **S5** A direction unit test was added in
  `inst/tinytest/test-exact-trend.R`. [verified] (see Section 4)

## 2. Deferred

- **M5(b), rest of the broader designs** (k in {3, 5, 6}, other
  unequal n_i, unequal score spacing). Deferred for budget. These
  are cheap with the existing code
  (`enumerate_trend_tables(n_per_group, scores)` with
  `exact_rejection()`; guard at 2e5 tables). Listed as Future
  research item 2.
- **M6, unconditional exact test (Tang et al. 2000; Shan et al.
  2012) and the `coin` cross-check.** The unconditional test is a
  substantive new implementation. The only `coin` install is in
  `/opt/homebrew/lib/R/4.6/site-library`. Its compiled code
  (`R_ExpectationCovarianceStatistic`) segfaults under this
  machine's R 4.6.0, even with that library placed last. The
  pure-R packages from that library (CATTexact) load fine. A
  `coin` cross-check of `p_exact_mass` needs a working `coin`
  install, for example in the Docker image.
- **m16, persistent identifier.** No DOI exists. An HTML TODO comment
  was left in Data and code availability. The author should add a
  Zenodo DOI, or the repository URL with the submitted commit.
- **m17, stale report.pdf / report.tex.** Now FIXED: the manuscript
  was re-rendered through `tools/render.sh`. This updated
  `report.pdf`, `report.tex`, the `report_files/` figures and
  `share/MANIFEST.md`, and staged
  `share/report-2026-11-17-1824-af32e81-wip.pdf`. The render used
  the scratch `R_PROFILE_USER` described in Section 3. A
  reproducible container render
  (`make docker-render REPORT=analysis/report/report.Rmd`) was not
  run.
- **S4, implement the tree algorithm.** Out of scope; the section
  stays labeled as unimplemented.
- **Significance (Section 3 of the report), venue and novelty.** The
  framing was narrowed in 1.7 and the abstract. The venue decision
  is the author's.

## 3. New issues found while fixing

- **CATTexact 0.1.1 gives incorrect p-values for some score vectors
  that contain 0.** With n = (3, 3, 2, 3) and x = (0, 2, 0, 0),
  scores (0, 3, 4, 5) give an upper tail of 0.673. Scores
  (1, 4, 5, 6) give 0.509, which matches direct enumeration. The
  conditional test must be invariant to a shift of the scores.
  Scores (0, 1, 2, 5) were not affected, so the trigger is not
  simply 'contains 0'. The manuscript's scores (1:4) are not
  affected, and neither are their reflection. I checked 1,500
  sampled tables of the main design: both CATTexact tails match
  enumeration to 5e-15. `exact_trend_tails()`
  now shifts scores to a minimum of 1. A tinytest records the
  discrepancy. This may be worth reporting upstream to the package
  maintainer.
- **purrr was not used.** purrr is not in `renv.lock` or
  DESCRIPTION. The new code uses base `vapply()`, which avoids adding
  a dependency. This departs from the user's purrr preference.
- **Local R library setup (environment, not manuscript).** The
  native library (`/Library/Frameworks/R.framework/.../library`)
  has knitr, kableExtra and bookdown, but not CATTexact.
  CATTexact is only in the Homebrew library. Putting the Homebrew
  library first makes compiled packages segfault. My early
  segfaults were caused by that ordering, not by broken packages.
  The tests and render here used a scratch `R_PROFILE_USER` that
  sources the project `.Rprofile` and then appends the Homebrew
  library last. A plain host render (`bash tools/render.sh ...`)
  will fail with 'no package CATTexact' until CATTexact is
  installed in the native library.
- **Mid-p and asymptotic tests nearly coincide in this design.**
  Their exact rejection probabilities are equal to four decimals,
  and their Monte Carlo rates are identical. They still disagree on
  124 of the 14,641 tables. The manuscript says only that they are
  'close'.

## 4. Verification record

- `inst/tinytest/`: 23 of 23 assertions pass (10 in `test-basic.R`,
  13 in the new `test-exact-trend.R`). Run them with
  `Rscript -e "tinytest::run_test_dir('inst/tinytest')"` from the
  `testtrend` root. The new file calls `exit_file()` and skips
  silently unless CATTexact is on the library path.
- Full B = 2,000 rerun of all three scenarios: done outside knitr
  with the manuscript's seed and call order. Elapsed time was 10.3
  minutes under CPU contention. An unloaded estimate of about 3
  minutes is extrapolated from a 50-replication run, not measured.
- A full render through `tools/render.sh` completed. Every
  `stopifnot()` guard passed:
  - the same NA rows across the trend methods;
  - mass-based equal to equal-tail in the main design;
  - mass-based differing from equal-tail in the unequal design.

  No undefined citations or cross-references appeared ('??' or raw
  @keys). Inline values were spot-checked against the computed
  values.
- Two render attempts failed before the final one, and both were
  fixed:
  - My own guard caught a false claim, that the one-sided
    asymptotic test exceeds alpha for every p >= 0.13. The prose
    now reports 21 of 49 grid points and calls the difference from
    Corcoran and Mehta unresolved.
  - The knitr cache dropped objects created by `source()`, so
    `sim-setup` is now `cache=FALSE`.
- The size-curve figure was inspected visually. The other three
  figures were rendered but not inspected.
