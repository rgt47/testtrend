# Morris et al. (2019) ADEMP Audit: 32-test-for-trend
*2026-04-17 09:02 PDT*

## Scope

Files audited:

- `analysis/scripts/sim_trend_test.R`
- `analysis/report/report.Rmd`

## ADEMP scorecard

| Criterion | Status | Evidence |
|---|---|---|
| Aims explicit | Partial | narrative only |
| DGMs documented | Met | trend DGM parameterised |
| Factors varied factorially | Partial | null and alternative scenarios separate |
| Estimand defined with true value | Met | trend slope input to DGM |
| Methods justified | Met | trend tests listed |
| Performance measures justified | Partial | rejection rate only |
| n_sim stated | Met | `R = 2000` |
| n_sim justified via MCSE | Not met | no derivation |
| MCSE reported per metric | Not met | `summarize_sim()` at `sim_trend_test.R:64-87` omits MCSE |
| Seed set once | **Not met** | `set.seed(seed)` at `sim_trend_test.R:7` inside `sim_trend_test()` — reseeded on every call |
| RNG states stored | Not met | not stored |
| Paired comparisons | Partial | methods applied to same dataset within rep |
| Reproducibility | Partial | per-call seed present; RNGkind not pinned |

## Overall verdict

**Not compliant.**

## Gaps

- `set.seed(seed)` is inside the worker function at
  `sim_trend_test.R:7`, reset on every call. The Rmd invokes
  `sim_trend_test()` twice (`seed = 42` for null, `seed = 99` for
  alternative), so each scenario resets the RNG.
  Morris §4.1 requires one `set.seed()` at program start, not
  per-scenario reseeds. Using distinct scenario seeds also breaks
  cross-scenario paired comparisons.
- `summarize_sim()` at `sim_trend_test.R:64-87` returns rejection_rate
  but no MCSE.
- `R = 2000` not justified via MCSE.
- `RNGkind()` not pinned.
- No `.Random.seed` per-rep capture.

## Remediation plan

1. Remove `set.seed()` from inside `sim_trend_test()`
   (`sim_trend_test.R:7`). Set a single seed in the caller (the Rmd)
   before the first invocation.
2. Allow the caller to pass a pre-generated L'Ecuyer stream for the
   scenario, removing the integer `seed` argument semantics entirely.
3. Add MCSE cols in `summarize_sim()`:
   `mcse_reject = sqrt(r * (1 - r) / n)`.
4. Derive `R` from a target MCSE. For rejection-rate MCSE ≤ 1 pp at
   p = 0.5, need R ≥ 2500.
5. Pin `RNGkind("L'Ecuyer-CMRG")` at the top of the Rmd.
6. Store `.Random.seed` per rep.
7. Add ADEMP Methods section to `report.Rmd`.

## References

Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate
statistical methods. Stat Med 2019;38:2074-2102. doi:10.1002/sim.8086

---
*Source: ~/prj/res/32-test-for-trend/testtrend/docs/morris-audit-2026-04-17.md*
