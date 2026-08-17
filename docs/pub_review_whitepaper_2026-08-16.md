# Referee Review: An Exact Conditional Test for Trend in
# r x 2 Tables via Fisher Enumeration

*Review date: 2026-08-16 16:09 PDT*

Reviewed workspace: `32-test-for-trend/testtrend`. Single manuscript
under review: `analysis/report/report.Rmd`, rendered as
`analysis/report/report.pdf` / `report.tex` (last staged render
`analysis/report/share/report-2026-08-15-1647-eba835e.pdf`).

## 1. Summary of the work under review

The manuscript proposes an exact conditional test for trend in
ordered r x 2 tables, framed as an extension of the authors'
companion tree-based/network-based enumeration machinery for the
Fisher exact test (a separate repository, `01-fisher-exact-rx2`, not
part of this workspace). The Introduction reviews the
Cochran-Armitage trend statistic, the asymptotic-versus-exact
literature, and existing R implementations (`CATTexact`, `coin`,
`perm`). The Methods section (a) restates the exact conditional
p-value as a sum of multivariate hypergeometric probabilities over
tables with `|T(t)| >= |T_obs|`, (b) draws an analogy to the Fisher
exact test's probability-ordered enumeration, and (c) derives
greedy `O(k-j)` bounds `T_min`/`T_max` for pruning a hypothetical
tree traversal of the trend statistic. The empirical section is a
small Monte Carlo study (ADEMP-structured, per Morris, White, and
Crowther 2019) with two scenarios (null and a dose-response
alternative), k = 4 groups of size 10, B = 2,000 replications,
comparing the asymptotic Cochran-Armitage test, the unordered
Fisher exact test, and `CATTexact::catt_exact()` (an existing
package, not a new implementation). The Discussion interprets the
results as confirming type I error control and a power advantage
over the unordered Fisher test. A "Future research" section lists
six extensions, including the C/C++ implementation of the
pruning algorithm that Section 2.3 describes. A final subsection,
nested under the "References" heading, reports a prior Morris et
al. ADEMP audit (`docs/morris-audit-2026-04-17.md`) with verdict
"Not compliant."

## 2. Major issues

1. **The paper's stated algorithmic contribution is never
   implemented, run, or benchmarked anywhere in this repository.**
   Location: `report.Rmd`, Sections "Connection to the Fisher exact
   test framework" and "Pruning adaptations" (lines 261-313); title
   ("...via Fisher Enumeration"). Inspected: a repo-wide search
   (`grep -rn "find_max\|find_min\|T_min\|T_max\|pruning"` over
   `analysis/scripts` and `R/`) finds no matches, and there is no
   `src/` directory or C/C++ code anywhere in the workspace. The
   tree-based pruning bounds `T_min`, `T_max` are derived
   mathematically but the empirical section instead calls
   `CATTexact::catt_exact()`, an existing third-party package
   (Methods chunk `sim-null`/`sim-alt`, `sim_trend_test.R:49-51`).
   "Future research" item 1 (lines 593-599) confirms this directly:
   "The present simulation uses the `CATTexact` package... A direct
   implementation using the tree-based and network-based algorithms
   ... would enable performance comparisons." This means the paper's
   title and its most novel methodological content (the pruning
   adaptation) are entirely unevaluated: no runtime comparison, no
   correctness check against the network algorithm, not even a toy
   worked example. A referee cannot assess a claimed algorithmic
   contribution that exists only as an untested mathematical sketch.
   Remediation: either (a) implement and benchmark the pruning
   algorithm before submission, or (b) reframe the paper entirely
   around what is actually done (a simulation study of an existing
   exact trend test against its asymptotic and unordered
   competitors), moving the pruning derivation to a clearly labeled
   "proposed algorithm, not yet implemented" subsection or removing
   it. See also Section 5 (Recommended framing) below.

2. **The two-sided exact p-value actually computed does not match
   the two-sided exact p-value defined in the Methods section.**
   Location: Methods, "The exact conditional trend test" (lines
   237-259) defines
   `p = sum_{t} Pr(t) * 1[|T(t)| >= |T_obs|]`.
   Inspected: the code that produces every reported p-value instead
   computes `pmin(2 * sim_null$p_exact_trend, 1)` and
   `pmin(2 * sim_alt$p_exact_trend, 1)` (`report.Rmd` lines 369-371,
   381-383), where `sim_null$p_exact_trend` /
   `sim_alt$p_exact_trend` are the one-sided exact p-values from
   `CATTexact::catt_exact()$exact.pvalue`
   (`sim_trend_test.R:49-51`). Doubling a one-sided exact p-value is
   a distinct construction from summing hypergeometric mass over
   `|T(t)| >= |T_obs|`; for a discrete, generally asymmetric null
   distribution of `T`, the two need not agree, and the doubling
   rule is known in the literature to be conservative (it can
   systematically inflate the p-value beyond the mass-based
   two-sided definition). This is very plausibly the mechanism
   behind an unexplained empirical anomaly in Table 1
   (`report.tex:507-509`, inspected in the rendered output): under
   the null, the exact trend test rejects at only 0.015, more
   conservative than the *unordered* Fisher exact test at 0.033,
   even though the ordered test conditions on a reference set that
   is a strict refinement and should, if anything, be closer to
   nominal than an unrelated test. An exact test markedly more
   conservative than a different exact test on the same data,
   with no discussion, is exactly the kind of result a statistical
   referee would ask the authors to explain or fix. Remediation:
   (a) implement the p-value exactly as defined in the Methods
   section (sum over `|T(t)| >= |T_obs|`) and compare it numerically
   against the doubled one-sided `CATTexact` p-value on a
   representative table, or (b) if the doubling rule is retained for
   practical reasons, state this explicitly in the Methods section
   as the operational definition, cite the doubling literature, and
   discuss its known conservatism in the Discussion when
   interpreting the size and power results.

3. **The power comparison at nominal alpha is confounded by unequal
   realized type I error and is not adjusted for or even
   flagged.** Location: Discussion, second paragraph (lines
   556-562); Table 2 (`report.tex:533-535`, inspected). The realized
   null rejection rates (Table 1) are 0.043 (asymptotic), 0.033
   (Fisher), and 0.015 (exact trend) at nominal alpha = 0.05 — a
   roughly 3-fold spread in actual size. The power figures (Table 2)
   are 0.426, 0.266, and 0.370 respectively. The Discussion states
   only that "the exact trend test recovers the power advantage of
   incorporating the ordering, which the unordered Fisher test
   cannot exploit," comparing 0.370 to 0.266 at face value, and
   separately notes the asymptotic test's higher power (0.426) only
   as evidence that it "may exhibit slight inflation" of size,
   without connecting the two facts. Comparing power across tests
   with substantially different realized sizes is not a valid
   like-for-like comparison; a referee would require either
   size-adjusted (empirical-alpha-matched) power curves or an
   explicit acknowledgment that the power ranking partly reflects
   differences in conservatism rather than differences in
   statistical efficiency at a common actual type I error rate.
   Remediation: report power at the empirical (simulation-based)
   critical value for each test, or at minimum discuss the
   confound explicitly and temper the power claim accordingly.

4. **The manuscript contains a self-contradicted, unresolved ADEMP
   compliance verdict.** Location: "References" section, subsection
   "Morris et al. (2019) ADEMP Compliance" (lines 635-650), citing
   `docs/morris-audit-2026-04-17.md`. Inspected: the audit
   (dated 2026-04-17) reports "Not compliant," citing (a)
   `set.seed()` called inside `sim_trend_test()` and reseeded per
   call, (b) two distinct scenario seeds preventing paired
   comparison, and (c) `summarize_sim()` returning no Monte Carlo
   SE. Inspected: the current `sim_trend_test.R` (verified by direct
   reading) contains none of these defects — there is no
   `set.seed()` call inside the function (a comment at line 6-7
   explicitly states "this function does not call `set.seed()` and
   must not do so"), a single `RNGkind("L'Ecuyer-CMRG"); set.seed(
   20260417)` is set once in `report.Rmd` (lines 324-325) before both
   scenarios run in sequence, per-replication `.Random.seed` is
   captured (`sim_trend_test.R:17`), and `summarize_sim()` computes
   `mcse_reject` and `mcse_mean_p` for every method
   (`sim_trend_test.R:88, 90`). The remediation described in the
   audit has evidently been carried out, but the manuscript still
   asserts, in the version submitted for this review, that the study
   is "Not compliant" and repeats the now-stale list of gaps
   verbatim. A referee reading the manuscript as written would
   conclude, correctly based on the text alone, that the simulation
   design is non-compliant with basic reporting standards — an
   incorrect and self-undermining statement given the actual code.
   This also raises a process concern: is the whitepaper/manuscript
   pipeline re-synchronizing prose with code after each remediation?
   Remediation: rewrite this subsection to reflect the current state
   (compliant, with the residual gap that B = 2,000 was not derived
   from a target MCSE a priori — see Minor Issue 1), move it out of
   the "References" heading into Methods or an appendix, and delete
   the stale gap list or explicitly mark it as historical
   ("as of the 2026-04-17 audit; since remediated").

5. **No results are reported in the past tense from the actual
   simulation; the Results section narrates expectations, not
   findings.** Location: "Results," subsections "Type I error
   control" and "Power comparison" (lines 415-421, 451-456).
   Inspected: the prose accompanying both tables is written entirely
   in the conditional/future mood — "The asymptotic Cochran-Armitage
   test *should* reject at approximately 5%," "The exact conditional
   test *is expected* to be slightly conservative," "The exact trend
   test *is expected* to have higher power" — without ever stating
   what the tables in fact show (0.043, 0.033, 0.015 for size;
   0.426, 0.266, 0.370 for power). No numeric value from Tables 1-2
   appears anywhere in the running text. This reads as boilerplate
   written before the simulation was run and never updated to report
   the actual numbers, and it obscures the anomaly in Major Issue 2.
   Remediation: rewrite the Results prose to state the observed
   values explicitly and interpret them (including the anomaly that
   the exact trend test's realized size, 0.015, is well below
   nominal and below the unordered Fisher exact test's size).

6. **Undocumented loss of replications is not disclosed or
   explained.** Location: Table 2, `n_valid` column
   (`report.tex:533-535`, inspected): 1,998 valid replications for
   the asymptotic and exact trend methods versus 2,000 for Fisher,
   under the alternative scenario (B = 2,000 stated in the
   Simulation size subsection, lines 355-360). Inspected: in
   `sim_trend_test.R`, `z_ca` (and hence `p_asymp`) is set to `NA`
   when `var_num <= 0` (lines 32-36), which occurs for degenerate
   tables (e.g., a column margin of 0); `p_exact_trend` presumably
   also returns `NA` on the same or a correlated failure via its own
   `tryCatch`. Neither the Results text nor the Discussion mentions
   that any replications were dropped, why, or whether the 2
   dropped replications share a common cause across the two
   affected methods. Silent listwise deletion of Monte Carlo
   replications, even at this small scale, should be disclosed and,
   per Morris et al. (2019) section 4, ideally addressed by design
   (e.g., excluding degenerate parameter configurations) rather than
   by ad hoc `NA` handling. Remediation: report and briefly explain
   the count and cause of excluded replications for every table
   footnote, and consider whether `n_rep` should be inflated to
   compensate.

## 3. Minor issues

1. **Simulation size not derived from a pre-specified target MCSE.**
   Location: "Simulation size" (lines 355-361). The manuscript
   computes the MCSE *post hoc* for the realized rejection
   proportions (about 0.05 and about 0.8) rather than deriving B
   from a worst-case (p = 0.5) target MCSE as recommended by Morris
   et al. (2019) and as the prior audit's remediation plan step 4
   proposed (B >= 2,500 for MCSE <= 1 pp at p = 0.5). The chosen
   B = 2,000 happens to be adequate for the realized proportions but
   was not justified a priori under the referenced framework.

2. **A citation is used in prose without a corresponding
   bibliography entry.** Location: line 165, "the theory of
   `@mehta1992exact` and Strasser-Weber (1999)." Inspected:
   `references.bib` contains no entry resolvable to Strasser and
   Weber (1999) (the `coin` package's permutation-test framework
   reference); the name is written as plain text rather than a
   `@` citation key, so it does not appear in the reference list and
   is not a resolvable citation. Add a bib entry (Strasser & Weber,
   1999, *Computational Statistics*) and cite it properly.

3. **An unused bibliography entry.** `mehta1980network` is defined
   in `references.bib` but never cited in `report.Rmd` (confirmed by
   cross-referencing all `@`-keys used in the text against the bib
   file). Either cite it (it may be the intended source for the
   "network algorithm" claim currently attributed only to
   `mehta1983network`) or remove it.

4. **Unsupported claims about commercial/other software internals.**
   Location: lines 154-156, "[the network algorithm] underlies the
   commercial software StatXact... and has been adopted as the
   standard backend for exact trend testing in SAS PROC FREQ and
   SPSS Exact Tests." No citation supports the SAS/SPSS claims
   specifically (only `mehta1983network` is cited, which documents
   StatXact). This is plausible but unverified in this review; add
   citations to SAS/SPSS documentation or soften to "is understood
   to underlie."

5. **Dependencies used in the analysis are not declared in
   `DESCRIPTION`.** Inspected: `DESCRIPTION` `Suggests:` lists only
   `here, knitr, rmarkdown, rprojroot, tinytest`. `CATTexact`
   (`sim_trend_test.R:50`) and `kableExtra`
   (`report.Rmd` setup chunk, line 42) are both used but absent from
   `DESCRIPTION`, even though both are pinned in `renv.lock`. This
   is a reproducibility metadata gap: a fresh `renv::restore()`
   driven only by `DESCRIPTION` would not surface these as required
   packages for `R CMD check` or `devtools::install_deps()`.
   `CATTexact` could not be verified as installed in this review's R
   environment (`requireNamespace("CATTexact")` returned `FALSE`),
   underscoring the need for explicit declaration.

6. **No unit tests exist for the analysis code.** Inspected:
   `inst/tinytest/test-basic.R` contains only
   `expect_true(TRUE)`. There is no test of `sim_trend_test()`'s
   trend-statistic computation, no test of `summarize_sim()`'s MCSE
   formulas, and no regression test pinning a known table's exact
   or asymptotic p-value against a hand-computed or externally
   verified value. For a manuscript whose empirical content rests
   entirely on this code, some minimal correctness tests (e.g.,
   comparing `z_ca` against `stats::prop.trend.test()` or a manual
   calculation on a textbook table) would materially strengthen
   confidence in the reported numbers.

7. **Leftover scaffold content unrelated to the manuscript.**
   `analysis/data/README.md` documents a Palmer Penguins data
   pipeline (raw/derived CSVs, `scripts/01_data_preparation.R`) that
   has no connection to this trend-test simulation study and does
   not exist in the repository; it also contains unfilled template
   placeholders (`[Analyst name, email]`, `[YYYY-MM-DD]`).
   `analysis/data/derived_data/sim_32.rds` is present but never
   referenced by any script or the manuscript (confirmed by
   `grep -rn "sim_32"` returning no hits). `analysis/tables/` and
   `analysis/figures/` are empty directories (figures are instead
   written to `figure/` and `analysis/report/report_files/`).
   Remove or repurpose this leftover zzcollab scaffold material
   before submission; an inconsistent `data/` directory undermines
   the paper's own reproducibility narrative.

8. **No abstract.** The manuscript has no Abstract section; it
   begins directly with "# Introduction." Statistical journals in
   the target class (JASA, Biometrics, JCGS, R Journal) require a
   structured or unstructured abstract stating the problem,
   approach, and findings.

9. **No explicit data/code availability statement**, despite the
   paper depending on a companion project (`01-fisher-exact-rx2`)
   not included in this repository. A referee will ask where the
   code for the described (but unimplemented, see Major Issue 1)
   pruning algorithm resides, and whether the simulation code and
   `CATTexact` version used are archived for reproducibility (a DOI
   or tagged release, per journal policy).

10. **Structural placement of the ADEMP compliance note.** The
    subsection discussed in Major Issue 4 is nested under the
    "# References" heading (line 633), which pandoc/citeproc treats
    as the bibliography anchor. Placing prose content there is
    unconventional and will likely render awkwardly (a subsection
    heading immediately followed by an auto-generated reference
    list). Move it to Methods, an appendix, or supplementary
    material.

## 4. What remains to be done

**(a) Required for correctness**

- Resolve the mismatch between the Methods-section p-value
  definition and the doubled one-sided p-value actually computed
  (Major Issue 2); verify numerically on at least one worked table.
- Investigate and explain (or fix) the anomaly that the exact trend
  test is more conservative than the unordered Fisher exact test
  under the null (Major Issues 2-3).
- Rewrite the stale "Not compliant" ADEMP verdict to match the
  current, apparently remediated, simulation code (Major Issue 4).
- Disclose and explain the dropped replications in Table 2 (Major
  Issue 6).
- Rewrite Results prose to report actual observed values rather
  than expected/conditional-mood statements (Major Issue 5).

**(b) Required for acceptance**

- Either implement and benchmark the pruning algorithm described in
  Section 2.3, or reframe the paper so its title and contribution
  claims match what was actually built and evaluated (Major Issue
  1; see Section 5 below).
- Add a size-adjusted or at least explicitly caveated power
  comparison (Major Issue 3).
- Add an abstract, a data/code availability statement, and a
  properly scoped limitations subsection (Minor Issues 8-9).
- Declare `CATTexact` and `kableExtra` in `DESCRIPTION` (Minor
  Issue 5).
- Add the missing Strasser and Weber (1999) citation and remove or
  use the orphaned `mehta1980network` entry (Minor Issues 2-3).
- Add minimal correctness tests for the trend-statistic and MCSE
  computations (Minor Issue 6).

**(c) Desirable polish**

- Justify B = 2,000 a priori from a target MCSE rather than post
  hoc (Minor Issue 1).
- Remove or repurpose the unrelated Palmer Penguins scaffold content
  and the orphaned `sim_32.rds` file (Minor Issue 7).
- Move the ADEMP compliance note out from under the "References"
  heading (Minor Issue 10).
- Support or soften the SAS/SPSS internals claim (Minor Issue 4).

## 5. Recommended framing

**Plausible framings.**

1. *New algorithm paper* (tree-based pruning for the exact trend
   test, analogous to the companion Fisher exact test work) — the
   framing the current title and Methods sections imply.
2. *Simulation/computational comparison paper* — an ADEMP-compliant
   Monte Carlo comparison of existing exact and asymptotic trend
   tests (asymptotic Cochran-Armitage, unordered Fisher, and the
   `CATTexact` conditional exact test) in the small-sample regime,
   contributing evidence rather than new computational machinery.
3. *Methodological note / tutorial* connecting the exact trend test
   to the r x 2 Fisher exact test enumeration framework conceptually,
   aimed at practitioners choosing among existing R tools.

**What the literature already covers.** Inspected via the citation
list and Introduction: the exact conditional trend test itself is
decades old (Mehta, Patel, and Senchaudhuri 1992; Agresti 1990,
1992, 2001), the network algorithm is the established computational
standard (underlying StatXact and, per the manuscript, SAS/SPSS),
and at least three R packages already provide exact or near-exact
trend tests (`CATTexact`, `coin`, `perm`). The manuscript's own
Introduction (lines 174-193) correctly identifies the one place
where a literature gap plausibly exists: no existing implementation
uses tree-based traversal with memoized pruning for the trend
statistic specifically (as opposed to the network/state-space
algorithm). That gap is real and citable, but this manuscript, as
submitted, does not fill it — it only asserts the gap and sketches
bounds. Framing 1 is therefore not currently earned by the content
of this repository.

**Recommendation.** Reframe as Framing 2 (simulation/computational
comparison paper), with Section 2.3's pruning derivation retained
only as a short "proposed extension" note (or moved to the "Future
research" section where it already partially lives) rather than as
a load-bearing contribution. This is the framing that matches what
is actually implemented, run, and verifiable in this repository: a
small but honest ADEMP-structured comparison of three existing trend
procedures. Reasoning: (i) the pruning algorithm is unimplemented
and unbenchmarked (Major Issue 1), so no referee at JCGS or the R
Journal — venues that expect either working software or
demonstrated algorithmic performance — would accept it as the
paper's headline contribution in its current state; (ii) the
Monte Carlo comparison, once the p-value definition mismatch (Major
Issue 2) and the size-power confound (Major Issue 3) are fixed, is a
legitimate, modestly novel contribution because it directly compares
type I error and power of the asymptotic, unordered-exact, and
ordered-exact tests side by side in a small, sparse-table regime
that is exactly where practitioners need such guidance and where
existing comparisons (Williams 1988; Lin and Chin 2006) are older or
narrower in scope.

**Implications under the recommended framing.**

- *Title*: something like "Type I Error and Power of Exact and
  Asymptotic Trend Tests in Small Ordered r x 2 Tables: A Simulation
  Study," dropping "via Fisher Enumeration" until an implementation
  exists.
- *Abstract/Introduction*: lead with the practical question (which
  trend test to use with small, sparse dose-response data) rather
  than with the algorithmic connection to the Fisher exact test
  enumeration machinery; keep the enumeration connection as
  motivation for future work, not as the stated objective.
- *Comparators*: the current three-way comparison is reasonable but
  thin; a referee is likely to ask for more DGM configurations
  (varying k, unequal n_i, unequal spacing of scores, more effect
  sizes) and possibly the unconditional exact test of Shan and Ma
  (2012) or Consiglio, Ricci, and colleagues (already cited),
  which the current design omits from the empirical comparison
  despite citing it as a competitor with claimed higher power.
- *Emphasize*: the ADEMP-structured design, the honest reporting of
  MCSE, and (once fixed) the explanation of the size-power tradeoff
  across the three methods.
- *De-emphasize or move to supplementary/future work*: the detailed
  `T_min`/`T_max` pruning derivation and the framing sentences that
  imply the tree-based algorithm was built and evaluated here.
- *Target journal*: under this framing, a methods-in-practice
  outlet such as *Statistics in Medicine* (already the CSL style
  used) or *Communications in Statistics -- Simulation and
  Computation* is a better fit than JCGS or the R Journal, both of
  which expect either theoretical/algorithmic novelty or an
  installable software artifact as the primary contribution.

If the authors instead prefer to keep Framing 1 (the algorithm
paper), the path to a defensible submission requires implementing
the pruning traversal, validating it against `CATTexact` and the
network algorithm on a battery of tables, and reporting wall-clock
or node-count comparisons — effectively completing Future Research
item 1 and item 6 before this manuscript, not after it.

## 6. Assessment

**Verdict: major revision** (bordering on reject-and-resubmit if
Framing 1 is retained without an implementation).

Justification: the manuscript's literature review, mathematical
setup, and simulation design (ADEMP structure, MCSE reporting,
single-seed RNG discipline) are competently done and largely
verified by direct inspection of the code that produces them. But
three defects are serious enough that a referee at any of the
target-class journals would not accept without substantial rework:
(1) the paper's title and core methodological narrative describe an
algorithm that is never implemented or evaluated anywhere in the
repository (Major Issue 1); (2) the actual computed p-value does not
match its own stated definition, and this mismatch is a plausible
direct cause of a striking, undiscussed empirical anomaly (Major
Issues 2-3); and (3) the manuscript's own embedded audit section
asserts the simulation is non-compliant with reporting standards
when the underlying code, on inspection, appears to have already
been remediated (Major Issue 4) — a manuscript that does not
accurately describe its own code cannot be taken at face value by a
referee for any of its other empirical claims either. None of these
are fatal to the underlying research program, but all three require
author action, not just a copyedit, before this paper is
publication-ready.

## 7. Revision history

- 2026-08-16: Initial pub_review whitepaper. No prior
  `pub_review_whitepaper_*.md` existed in `docs/`. Identified six
  major issues (unimplemented core algorithm; p-value definition
  mismatch between Methods and code; unadjusted power comparison
  confounded by unequal realized size; stale self-contradicted
  Morris ADEMP compliance verdict; expectation-toned rather than
  results-toned Results prose; undisclosed dropped replications)
  and ten minor issues (simulation-size justification, missing and
  orphaned bibliography entries, unsupported SAS/SPSS claim,
  undeclared `CATTexact`/`kableExtra` dependencies, absent unit
  tests, leftover Palmer Penguins scaffold content, missing
  abstract, missing data/code availability statement, and
  misplaced ADEMP note). Recommended reframing the paper as a
  simulation/computational comparison study rather than an
  algorithm paper, given that the tree-based pruning algorithm
  described in Methods is not implemented anywhere in this
  workspace. Overall verdict: major revision.
