# Voice Revision Log: report.Rmd
*2026-11-17 17:26 PST*

Pass: prose revision of `analysis/report/report.Rmd` into the
`/use_rgtvoice` register. Changes are in the working tree only. Nothing
was staged, committed or rendered.

## Rule applied

This pass applied the stricter 'Short, simple sentences' rule in the
current `use_rgtvoice.md`:

- one idea per sentence;
- a target of 15 to 20 words, with 25 as the ceiling;
- subordinate clauses split into separate sentences, each led by a
  plain connective;
- common words preferred over Latinate ones, except for statistical
  terms of art.

All content and hedges were kept. Long sentences were split rather
than shortened.

## Kinds of changes made

The edits cover about 33 prose passages, which appear as 39 diff
hunks against `HEAD`.

- **Long sentences split.** About 30 prose sentences over 25 words
  were split, measured by script before and after. The second
  sentence now opens with a plain connective, such as 'As a result',
  'That is,', 'There,', 'It also', 'Nevertheless,' or 'First, Second,
  Third'.
- **Trailing clauses moved out.** Trailing ', which ...' clauses and
  participial clauses became separate sentences. Examples:
  - the three 'Doing so ..., which ...' sentences in Methods 2.1,
    Results 3.2 and the Discussion;
  - the SAS/SPSS caveat;
  - the `coin` null-distribution options;
  - the degenerate-margin explanation;
  - the doubling-heuristic sentence;
  - the Morris worst-case MCSE sentence.
- **Lists given a short lead sentence.** Two list-heavy sentences were
  reordered so that a short head sentence comes first. These are the
  study objectives in Section 1.7 and Limitations item 3.
- **Inline lettering replaced.** In the Morris compliance paragraph,
  the inline enumeration '(a) ... (b) ... (c)' became 'First, Second,
  Third'.
- **Trade-off made explicit.** The trade-off for the unconditional
  approach is now an 'On the plus side ... On the minus side' pair.
- **Common words.** 'Demonstrated' became 'showed' in two places.
  'Problematic' is used once, in its calibrated sense.
- **Inline R lists moved.** The lists of inline R rates in Results 3.2
  and Discussion paragraph 2 now sit in their own sentences. The same
  inline expressions appear in the same order.

## Invariants checked by script against `HEAD`

The following are identical to `HEAD`:

- the YAML header and the abstract;
- every code chunk;
- all display math and all inline math;
- all inline R expressions;
- all `\@ref()` cross-references and all cite keys;
- all headings;
- the multiset of numerals in the body.

No em or en dashes were added.

After revision, the script counts four sentences over 25 words. None
of them is a real long sentence:

- Two are artifacts of the script. One runs into a display equation
  and one runs into a markdown list.
- One is the ADEMP 'Methods' enumeration, which was left as an
  enumeration.
- One is Future research item 6. The sentence itself is 24 words; the
  count includes the bold label.

Some sentences of 21 to 25 words remain. These are mostly definitional
sentences built around inline math, or list sentences in which
further splitting would leave fragments. Short prepositional uses of
'rather than' and one short 'because' clause in the Discussion were
also kept, because they do not stack clauses.

## Passages deliberately left alone

- **Abstract.** It is excluded by the skill and by instruction.
- **Chunk, table and figure captions, and figure content.** These are
  inside code chunks.
- **Headings and bold run-in labels.** The `\@ref()` anchors are
  derived from the headings.
- **Limitations '(1) to (5)' and the Future research numbering.**
  Other passages refer to these as 'item 1' and so on.
- **The ADEMP template (Aims, Estimands, Methods).** It is a
  structured reporting template, not narrative prose. Only the
  Simulation size sentences were split.
- **Substantive claims flagged in the referee report.** Per
  instruction, these were not corrected. Where a sentence was split,
  the claim's wording was preserved. The claims are:
  - the 'two-sided' doubled p-value framing (Major 1);
  - the uncentered T(t) definition (Major 2);
  - the 'refinement of the reference set' argument in Results 3.1
    and the Discussion (Major 3);
  - 'well below 1 pp' and the MCSE justification that assumes p near
    0.8 (Major 4);
  - 'He' for the multi-author works @agresti1990exact and
    @mehta1992exact (Minor 4);
  - the placement of the @bartholomew1959test citation (Minor 5).
- **'as a plausible explanation, not a demonstrated one'.** This was
  kept because 'demonstrated' here means 'proven'. Replacing it with
  'shown' would weaken the contrast.
- **The pre-existing suppressed-author citation
  '[-@morris2019using]'.** After the split it follows 'recommended
  by', as it did before. It was left as is because citation syntax
  was out of scope.
- **Rhetorical questions and worked examples.** None were added,
  because they would introduce new content rather than revise
  existing prose.
