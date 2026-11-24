# Referee Report: Type I Error and Power of Exact and Asymptotic Trend Tests in Small Ordered r x 2 Tables
*2026-11-17 17:17 PST*

Manuscript: `analysis/report/report.Rmd` (workspace
`32-test-for-trend/testtrend`). Reviewed at commit state with a clean
working tree.

## 1. Summary of the manuscript

The manuscript reports a small Monte Carlo study, framed by the
ADEMP structure of Morris, White and Crowther (2019). It compares
three tests for trend in a 4 x 2 table with 10 subjects per group:
the asymptotic Cochran-Armitage (CA) test, the unordered Fisher exact
test, and the exact conditional CA test from the `CATTexact` package.
For the exact test, the two-sided p-value is formed by doubling the
package's one-sided p-value. There are two data-generating
mechanisms: a null with p = 0.2 in every group, and one increasing
alternative, p = (0.05, 0.10, 0.20, 0.35). Each uses B = 2,000
replications. The authors report realized null rejection rates that
differ by about three-fold (0.043, 0.033, 0.015). They attribute the
low size of the exact trend test partly to the doubling rule. They
also sketch, without implementing it, a tree-based pruning scheme that
adapts companion Fisher-exact enumeration code to the trend
statistic. The stated contribution is the simulation comparison
itself.

## 2. Overall assessment and recommendation

**Recommendation: reject in current form (for a top-tier journal);
major revision for a methods-note or applied venue.**

The study is small (one null, one alternative, one design) and is
carefully hedged about what it does not do. That candor is a
strength. However, the central empirical finding rests on a
construction error. `CATTexact::catt_exact()` returns an upper-tail
p-value for an *increasing* trend only. Doubling an upper-tail
p-value does not yield a two-sided test. It yields a one-sided test
at level alpha/2 that can never reject a decreasing trend. The
'anomaly' that the Discussion attributes to the conservatism of
doubling under an asymmetric null distribution is, in large part, an
artifact of this directional error (Major comment 1). After that
correction the paper has no new method, no new theory and one design
point. That falls well short of the bar for JASA, JRSS-B,
Biometrika or *Statistics in Medicine*. A corrected and broadened
study could suit *Journal of Statistical Computation and Simulation*,
*Communications in Statistics: Simulation and Computation*, or
*Pharmaceutical Statistics*.

## 3. Significance and novelty

Exact conditional versus asymptotic versus exact unconditional trend
tests for small binomial dose-response designs have been compared
before, at larger scale and with more designs:

- Tang, Ng, Guo, Chan and Chan (2006), JSCS 76(10):847-859, compare
  exact conditional and unconditional CA tests across logit, probit,
  one-hit and extreme-value models. They find the conditional test
  conservative and less powerful in small samples. (This paper is
  already in the bibliography, but under the wrong authors; see
  Section 6.)
- Shan, Ma, Hutson and Wilding (2012), *Stat Med* 31:155-164, and
  Consiglio, Shan and Wilding (2014), *Int J Biostat* 10:221-230,
  compare conditional and unconditional exact trend tests.
- A paper in *J Mod Appl Stat Methods* 1(1), attributed in the
  search results to Corcoran and Mehta, reports exact level and power
  of permutation, bootstrap and asymptotic trend tests. It uses 40
  mice in four dose groups, which appears to be the present design.
  (This is from search snippets only. I did not confirm the title,
  year or pages.)

Against this literature the present study adds no new design point,
statistic or comparator. The tree-based pruning idea is the only
potentially novel element, and it is explicitly unimplemented. As a
result the manuscript has no claim to methodological novelty in its
current form. Its framing is also dated: there is no engagement with
mid-p values (Lancaster 1961; Hirji 2006, ch. 2), with unconditional
exact tests as comparators, or with robust trend tests that
combine several score vectors (MAX-type tests). My searches did not
locate a post-2015 exact-trend comparison study. The absence of
recent citations may therefore reflect the literature as well as the
manuscript.

## 4. Major comments

**1. The exact 'two-sided' p-value is an upper-tail test at alpha/2
(Methods, 'Operational two-sided p-value'; chunks `sim-null`,
`sim-alt`; Results 3.1; Discussion para 2). Verified by algebra
and rerun.**

- The `CATTexact` DESCRIPTION states that the package computes
  'one-sided p-values'. Inspection of `catt_exact()` (v0.1.1, the
  version pinned in `renv.lock`) and direct calls confirm the
  direction. For counts (0, 1, 2, 4) over scores 1:4 with 10 per
  group, `exact.pvalue` = 0.012. For the reversed table (4, 2, 1, 0)
  it is 0.996.
- Hence `pmin(2 * p, 1) < 0.05` is identical to `p < 0.025`. This is
  an upper-tail test at level 0.025. It never rejects for a
  decreasing trend. The 'two-sided' exact test therefore spends only
  half of alpha.
- I reran the manuscript's simulation exactly (same `RNGkind`, seed,
  and `sim_trend_test.R`, under R 4.6.0 with `CATTexact` 0.1.1). It
  reproduces the reported null rates (0.0425, 0.0330, 0.0150). The
  median doubled p-value is 1.000, and 61.8 percent of null doubled
  p-values equal 1. Both are signatures of a one-directional
  p-value, not of a symmetric two-sided construction.
- I replayed the manuscript's own random stream (same `RNGkind`,
  seed and draw order) and computed the lower tail for each
  replicate by reversing the group order. Reversal is valid here
  because the scores are equally spaced and the groups are of equal
  size. Under the null, the equal-tail construction
  2 min(p_upper, p_lower) rejects at 0.027, against 0.015 for the
  doubled-upper rule. That is, the directional error roughly halves
  the realized size. (An earlier independent run, with B = 3,000 and
  a different seed, gave 0.029 against 0.014.) The corrected rate,
  0.027, is still below Fisher's 0.033 and the nominal 0.05. The
  remaining gap is plausibly due to discreteness (inferred, not
  shown), so part of the anomaly is still unexplained after the
  correction.
- Under the increasing alternative, the equal-tail rejection rate on
  the same stream is 0.370, identical to the reported value. The
  lower tail is essentially never reached there. The reported power
  of the exact trend test is therefore unaffected. The damage is to
  the null size, to any decreasing alternative, and to the
  interpretation.
- Consequences: (a) the Discussion's explanation (doubling is
  conservative under asymmetric discrete nulls, with p-values 'just
  above alpha/2') is not the main mechanism and should be rewritten;
  (b) the size comparison in Table 1 sets a one-sided alpha/2 test
  against two two-sided tests; (c) the Methods statement that
  doubling is 'a standard heuristic' presupposes that the one-sided
  p-value is taken in the direction of the observed statistic, which
  the code does not do.
- Remedy: compute both tails (reverse the score order or the group
  order) and report 2 min(p_upper, p_lower). Alternatively, report a
  mass-based two-sided p-value from an exact linear-by-linear
  permutation test. `coin::independence_test()` with ordered scores
  and an exact distribution may provide this without new enumeration
  code. I could not confirm this, because `coin` failed to load in
  the review environment. Rerun the null scenario and rewrite
  Results 3.1 and the Discussion. Add a decreasing alternative to the
  design.

**2. The two-sided p-value definition uses an uncentered statistic
(Methods 2.1, display for T(t) and p). Inspected.**

The Introduction defines Armitage's T as centered
(sum d_i X_i - R_1 sum d_i n_i / N). Methods redefines
T(t) = sum d_i x_i(t), which is uncentered and always nonnegative for
positive scores. With that definition, |T(t)| >= |T_obs| is a
one-sided (upper-tail) criterion, not a two-sided one. The same
defect carries into the pruning rule in 2.3. Remedy: define
T(t) = sum d_i (x_i(t) - n_i R_1 / N), or state the rejection region
as |T(t) - E_0 T| >= |T_obs - E_0 T|. Note that this error in the
displayed definition parallels the directional error in Major
comment 1.

**3. The 'refinement of the reference set' argument is incorrect
(Results 3.1 para 1; Discussion para 2). Inspected.**

The Results say that conditioning 'on a strict refinement of the
reference set' should bring the size closer to nominal. The
Discussion calls the trend test's reference set 'a refinement of the
Fisher reference set'. Both tests condition on the same margins and
use the same reference set, as the manuscript itself states in
Section 2.2. They differ only in the ordering that defines the
rejection region. There is no a priori reason why the trend
ordering should give realized size closer to nominal. Size depends on
the granularity of the attainable p-values under each ordering.
Remedy: delete the refinement argument. If a size explanation is
wanted, tabulate the attainable p-values near alpha under each
ordering for these margins.

**4. Monte Carlo design and MCSE claims are inconsistent with the
results (Simulation size paragraph; Limitations item 5). Verified.**

The text justifies B = 2,000 by assuming power near 0.8 and states
that both MCSEs are 'well below 1 pp'. The realized powers are 0.27
to 0.43, and the reported MCSE for rejection under the alternative
is 0.0099 to 0.0111. Two of three exceed 1 percentage point. The
claim is false as stated. Remedy: either increase B (B >= 2,500 gives
MCSE <= 0.01 at any rate), or restate the target and acknowledge the
realized MCSE. With one alternative the comparison of 0.370 versus
0.426 is precise enough. However, a richer design will need a
pre-specified MCSE target per Morris et al.

**5. The design is too narrow to support any general conclusion
(ADEMP Data-generating mechanisms; Limitations item 1). Inspected.**

There is one k, one n, one score vector, one null p and one
alternative. Size of exact tests is a highly non-monotone function of
the nuisance probability. Reporting it at p = 0.2 alone can mislead
in either direction. Remedy: at minimum, (a) plot realized size over
a grid of the common null probability (e.g., 0.02 to 0.5); this needs
no Monte Carlo, since it can be computed exactly by summing over the
product-binomial sample space for k = 4, n = 10; (b) add unequal n_i,
k in {3, 5, 6}, and two or three alternative shapes (linear,
threshold, plateau); (c) include a decreasing alternative, which
would have exposed Major comment 1 immediately.

**6. Missing comparators that the manuscript itself names
(Introduction 1.5, 1.6; Limitations item 3). Inspected.**

The unconditional exact test of Shan et al. (2012) and the `coin`
exact linear-by-linear test are both discussed in the Introduction
but not run. The mid-p version of the conditional test is neither
discussed nor run. `coin` is an existing package. A mid-p value
needs the point probability of the observed statistic, which
`catt_exact()` does not return. It can, however, be obtained as
p_upper + p_lower - 1 from the two tail calls already needed for
Major comment 1. A study whose stated contribution is 'the
simulation comparison itself' needs these comparators. Without them it cannot
say anything new relative to Tang et al. (2006).

**7. The pruning sketch is incomplete as an algorithm (Section 2.3).
Inspected.**

The stated rule prunes a subtree only when no completion reaches the
rejection region. An efficient exact algorithm also needs the
complementary case: when every completion is in the rejection region,
the subtree's total probability should be added in closed form. This
is the analogue of find_max/find_min acceptance in the Fisher tree.
The greedy O(k - j) bound also needs a statement of the box
constraints 0 <= x_i <= n_i and the assumption that the scores are
sorted. Because the section is labeled as an unimplemented proposal,
this is a presentation issue rather than a correctness one. However,
the paper would be stronger without the section, or with the
section moved to an appendix with these gaps closed.

**8. Power interpretation omits the most salient result (Results 3.2;
Discussion para 3). Verified (numbers).**

The asymptotic CA test has the highest power (0.426) and a realized
size (0.043) closest to nominal. The text focuses on exact trend
versus Fisher. A reader would conclude from these data that, at
k = 4 and n = 10 with p near 0.2, the asymptotic test is adequate.
This is contrary to the motivating claim in Section 1.3. That claim
should be stated, and Section 1.3 should be qualified accordingly.

## 5. Minor comments

1. Section 1.2, line after the display for T: 'Here X_i is ...' should
   follow a display ending with a period, not a comma; or 'where'
   should follow the comma. The same issue occurs in Section 2.1.
2. Section 1.2: 'It is now standard in regulatory submissions' needs a
   citation (e.g., FDA 2001 carcinogenicity guidance, or
   Liu and Agresti 2005 for general use).
3. Section 1.3: 'by analogy, to general r x c tables' is out of place
   in a paper about r x 2 trend tests.
4. Section 1.4: 'He' is used for the three-author works
   @agresti1990exact and @mehta1992exact. Use 'They' or 'The
   authors'.
5. Section 1.5: '[@bartholomew1959test]' is attached to a sentence
   about Consiglio et al.'s finding. Move it to the first mention of
   Bartholomew's statistic.
6. Section 1.6: 'The `perm` package' is uncited. Cite Fay and Shaw
   (2010), *J Stat Softw* 36(2).
7. Section 1.6: the `CATTexact` DESCRIPTION says its algorithm is
   'following an idea by' Mehta et al. (1992). The manuscript says
   it 'implements ... using the Mehta-Patel network algorithm' and
   computes 'via network enumeration'. Soften to match the source.
8. Section 1.6: the claim about SAS PROC FREQ and SPSS is flagged as
   unverified. Either verify it (PROC FREQ `EXACT TREND` is
   documented) or delete it.
9. Section 1.6 para 3: 'we have shown that tree-based pruning can be
   competitive ...' needs a citation to the companion work.
10. Notation: 'group sizes' are called column totals in Section 1.5
    and row margins in Section 2.1. Choose one orientation. The code
    builds a 2 x k matrix, so groups are columns there.
11. Simulation-size paragraph: 'Monte Carlo SEs ... per Morris Table
    6' should read 'Morris et al. (2019), Table 6'.
12. ADEMP 'Estimands': a binary decision is not an estimand. The
    estimands are the rejection probabilities (size and power).
13. The 'Morris et al. (2019) ADEMP compliance' subsection reads as a
    change log (file line numbers, audit dates). Move it to a
    supplement or the repository README.
14. Table 2 caption: say 'the same two replications', since the code
    confirms the NA rows coincide (verified in my rerun).
15. Figure 3 (agreement scatter): because the exact trend p-value is
    one-directional, the scatter mixes the two tails. Re-draw after
    Major comment 1. Consider log scales.
16. Section headings: the Data and code availability section cites
    `renv.lock` and the companion repository by local path. Give a
    persistent identifier (Zenodo DOI or repository URL with commit).
17. The rendered `report.tex` and `report.pdf` in the repository are
    stale: their prose predates the current Rmd (e.g., 'should reject
    at approximately 5%'). Re-render before circulating.

## 6. Missing and questionable references

| Location | Issue | Suggested source |
|---|---|---|
| Bib key `lin2006exact` | Wrong authors and pages (verified via Crossref for DOI 10.1080/10629360600569519). | Tang M-L, Ng HKT, Guo J, Chan W, Chan BP-S. Exact Cochran-Armitage trend tests: comparisons under different models. *J Stat Comput Simul* 2006; 76(10):847-859. |
| Section 1.3, `lin2006exact` | Cited for 'asymptotic test can be anticonservative'. The abstract is about exact conditional versus unconditional tests and finds the *conditional* test conservative. Claim and source do not match (inspected from abstract only). | Corcoran C, Mehta C. Exact level and power of permutation, bootstrap, and asymptotic tests of trend. *J Mod Appl Stat Methods* 2002; 1(1), article 7 (pages to verify). |
| Bib key `williams1988tests` | DOI 10.2307/2347314 resolves to Rosenbaum (1988), a different paper (verified). | Correct DOI: 10.2307/2347316. |
| Section 1.3, `williams1988tests` | Cited as proposing an exact conditional trend alternative. The paper concerns tests for differences among several small proportions. Check that it supports the specific claim (suspected mismatch, not verified). | |
| Bib key `agresti1998order` | Volume and pages wrong (verified via RePEc and search). | Agresti A, Coull BA. *Comput Stat Data Anal* 1998; 28(2):139-155. |
| Bib key `edelmann2025catexact` | Year 2025; v0.1.1 was published on CRAN on 2020-06-25 (verified from installed DESCRIPTION). | Edelmann D. CATTexact, R package version 0.1.1, 2020. |
| Section 1.6, `perm` | Uncited software. | Fay MP, Shaw PA. Exact and asymptotic weighted logrank tests for interval censored data: the interval R package. *J Stat Softw* 2010; 36(2). (Cites `perm`; verify the preferred citation.) |
| Section 1.5 | Unconditional exact CA test predates Shan et al. | Tang M-L, Chan PS, Chan W. On exact unconditional test for linear trend in dose-response studies. *Biom J* 2000; 42(7):795-806 (details to verify). |
| Discussion | Mid-p as the standard remedy for discreteness. | Lancaster HO. Significance tests in discrete distributions. *JASA* 1961; 56:223-234. Hirji KF. *Exact Analysis of Discrete Data*, 2006, ch. 2. |
| Section 1.2 | 'Standard in regulatory submissions' unsupported. | US FDA (2001). Draft guidance: Statistical aspects of the design, analysis, and interpretation of chronic rodent carcinogenicity studies. |

## 7. Suggestions for strengthening the paper

1. Replace the Monte Carlo null with exact size computation. For
   k = 4 and n = 10 the sample space has 11^4 = 14,641 points, so size
   and power can be computed exactly for every method at every p.
   This removes MCSE from the comparison altogether.
2. Report size curves over the nuisance probability, in the style of
   Tang et al. (2006) and Consiglio et al. (2014).
3. Add the mid-p conditional test, the `coin` exact test, and the
   Shan et al. (2012) unconditional test as comparators.
4. If the tree-based enumeration is to remain the motivation,
   implement it and validate it against `coin` on a grid of tables
   before claiming any algorithmic contribution. Otherwise, remove it
   from the title-level framing.
5. Add a unit test in the workspace that checks the direction of
   `catt_exact()` on a decreasing table. This would have caught Major
   comment 1.

## 8. Scope of this review

**Verified (executed or re-derived):**

- Full rerun of both scenarios with the manuscript's seed and RNG
  (R 4.6.0, `CATTexact` 0.1.1 from `/opt/homebrew/lib/R/4.6/
  site-library`). It reproduces every number in Tables 1 and 2 to the
  reported precision. The NA rows under the alternative are identical
  for the asymptotic and exact trend methods.
- `CATTexact::catt_exact()` source and direction (upper tail only).
- Equal-tail versus doubled-upper rejection on a replay of the
  manuscript's own random stream: null 0.027 versus 0.015;
  alternative 0.370 versus 0.370. An earlier independent run
  (B = 3,000, different seed) gave 0.029 versus 0.014 under the null.
- Asymptotic CA variance formula in `sim_trend_test.R` (permutation
  variance with N - 1 denominator): matches the standard form.
- MCSE arithmetic in the Simulation size paragraph, and the B >= 2,500
  worst-case figure.
- Crossref metadata for the `lin2006exact`, `williams1988tests`,
  `consiglio2014comparison`, `shan2012efficient` and
  `mehta1992exact` DOIs.

**Inspected (read and judged):** all prose, displays, chunks and the
bibliography; `renv.lock` pins `CATTexact` 0.1.1;
`docs/morris-audit-2026-04-17.md` exists; `sim_trend_test.R` line
references cited in the text are accurate.

**Not checked:** the content of Williams (1988), Patefield (1982),
Mehta et al. (1984) and Liu and Agresti (2005) against their attached
claims; the full text of Tang et al. (2006) beyond its abstract;
the SAS/SPSS claim; figures (not rendered in this review). Whether
`coin` gives an exact two-sided linear-by-linear p-value for this
design (the package failed to load in the review environment). The
bibliographic details I suggest for Corcoran and Mehta, Lancaster
(1961), the FDA (2001) guidance, Fay and Shaw (2010) and Tang, Chan
and Chan (2000) come from memory or search snippets and need
checking. The attribution of the remaining size gap to discreteness
is an inference.

**Literature searches run (WebSearch / Crossref):** 'Exact
Cochran-Armitage trend tests: comparisons under different models'
(authorship); Agresti and Coull 1998 order-restricted (citation
details); comparisons of trend tests for binomial small samples;
exact CA small-sample simulation with two-sided and mid-p p-values,
2016 to 2024 (no direct post-2015 comparison found); Shan
unconditional exact trend tests, post-2012 (no newer Shan paper
located). Crossref bibliographic query for Williams (1988).
