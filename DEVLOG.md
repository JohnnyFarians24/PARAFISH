# DEVLOG

## 2026-07-24 — RULE: changing a test's FIXTURE is not changing its ASSERTION. And a floor with no boundary test has no proven value.

**The rule, because the question will come up again.** Task 3 added sequence validation
to the custom-DB upload, and three existing `DetectKindTests` broke: they used a 4 bp
`ACGT` body as a placeholder while `detect_fasta_kind` read only headers. The bodies were
grown to 80 bp and the tests passed again — which is the shape of the worst habit in
testing, so it needs to be told apart from it explicitly.

> **If the ASSERTION changed, suspect it. If only the INPUT changed and the assertion
> still demands the same thing, it was the fixture that was wrong.**

Here the assertions (`kind == "blast"`, `kind == "sintax"`) are untouched. What changed is
that the input became something a user could actually upload: a 4 bp subject can never
produce a blastn hit — the word size is 11, and megablast's is 28 — so the fixture was
exercising a scenario that cannot occur. It passed only because no validation existed;
once validation existed, 4 bp stopped being valid input at all. Weakening an assertion
lets the bug through; correcting an impossible fixture does not.

**And the reverse, which the review caught and I had not done.** The new floor
(`MIN_REFERENCE_LENGTH = 50`) had exactly one test, using an 8 bp sequence. The constant
could have been 10 or 100 and it would still have passed: **the floor had no proven
value, and the `<` versus `<=` choice was unproven.** It surfaced only because growing
those fixtures raised the question of where the boundary actually is.

Added: 49 rejected / 50 accepted, gaps excluded from the count (padding a 49-base record
with dashes must not buy it past a floor that exists to guarantee alignable sequence).
**And a second, deliberately literal test pinning the value 50 itself** — because the
boundary test derives its inputs *from* the constant, so it proves the comparison and not
the number. Verified by mutation: with the constant moved to 100, the first version still
passed; with the value pin added, 11 tests fail.

**The generalisation, which is the part worth carrying out of here:**

> **If you change the value and nothing goes red, the value was never pinned — it was
> only being used.** A scientific choice must force a test to change.

That is a sharper instrument than it looks, because it is *executable*. The first
boundary test looked rigorous, was written deliberately, and was checked into a suite of
540 green tests — and it was measured to be worthless by the simplest possible experiment:
change the number, re-run, see nothing fail. A test derived from the constant it means to
protect can only ever prove the code is self-consistent, which is the same failure mode as
`verify.py` check (f) (coherence is not correctness) and as `golden.py check` (a faithful
comparison of two stale files). Third time this project has met that shape, now in three
different instruments.

**On the list, not now: run the same mutation experiment on every other scientific
parameter** — `MIN_UNIQUE_SIZE`, PEAR's `min_overlap`, the clustering identity (0.97),
`sintax_cutoff` (0.8), and `trimmomatic.minlen`. Change each, re-run the suite, and record
which ones the suite fails to notice. Any that stay green are in exactly the position
`MIN_REFERENCE_LENGTH` was in, and every one of them is a documented, deliberate decision
that someone could currently alter without meeting any resistance at all. The 0.97 is the
one to check first: it is the parameter this pipeline's entire feature count rests on.

## 2026-07-24 — DECIDED (Task 4f, ahead of the work): the SINTAX columns are ALWAYS present, even when no SINTAX source ran.

After Task 1, a COI run legitimately has SILVA and PR2 marked `not_applicable`, so the
three `sintax_*` columns are empty in the normal case rather than the exceptional one.
The options were to omit them when no SINTAX source participated, or to carry a reason in
the file.

**Neither: keep them, always.** These tables are a deliverable other people will process
with their own scripts, and a schema that changes shape between runs breaks every one of
them — silently, since a missing column reads as a missing value. Schema stability is
worth more than self-explanation.

The explanation is not lost, and does not need to live in this file: `source_status.csv`
already records `skipped` / `not_applicable` per source **with the reason**, and the
run_config records which sources were requested. What is missing is the *link* being
documented, not a column disappearing.

## 2026-07-24 — KNOWN PROPERTY OF THE NET: `golden.py check` never re-runs the pipeline. "No regressions" means the frozen files have not changed — NOT that the code still produces them.

Found while deciding whether raising PEAR's `min_overlap` was safe. It is not, and the
way it is not-safe is the point: the fixture holds 83 merges with sub-50 bp overlap and
`_real` holds 1 649 (2.0% each), so the change would alter real output — **and
`golden.py check` would have stayed green through it.**

`check` reads a stored run directory (`runs/_fixture/current`) and diffs it against
`golden/fixture/`. Both sides are artifacts frozen in the past. Nothing in the loop
executes the pipeline, so **no change to code or parameters can make this check fail by
itself.** It fails only when someone re-runs the pipeline into that directory first —
which the discipline assumes and the tool does not enforce.

**The consequence, stated so the phrase cannot be misread again:** "golden: no
regressions" in any entry of this log means *the frozen files still match the run that
was last written beside them*. It does **not** mean the current code reproduces them. A
parameter changed without a re-run leaves the contract describing a pipeline that no
longer exists, with nothing to announce it.

Sixth incarnation of "looks alive, does nothing" in this project, and **the first in the
verification instrument rather than in the code being verified.** Sibling of what check
(f) of `verify.py` closed: there, internal coherence was not enough — a table can agree
with itself and still be wrong about the run. Here, faithful comparison is not enough —
two files can match perfectly and both be stale.

Nothing was silently wrong because of this: the min_overlap change was NOT made, and the
MINLEN change was verified safe by measuring the frozen inputs directly (zero pairs under
200 bp in either dataset) rather than by trusting the green.

**On the list, not now: a `--rerun` verification mode** that executes the pipeline over
the fixture and diffs the result against the frozen files. Cheap on the fixture (seconds,
26 MB) and expensive on `_real`. That would give the fixture a strong guarantee and leave
`_real` with the current weak one — **and the difference between the two would be
documented instead of implicit**, which is most of the value.

## 2026-07-24 — trim-minlen resolved DOWN at 50, against the first instinct, because the measurement said 200 was the wrong instrument. min_overlap is the right one, and it is OPEN.

The second divergence found beside `min_unique_size`: Setup offered 50, `RunConfig` said
200, so a run's read-length floor depended on whether it started from the UI or the CLI.
The instinct — mine and the reviewer's — was to align UP at 200. The reasoning was sound:
the Leray amplicon is ~313 bp and reads are 2×300, so a read trimmed to 50 leaves only
~11 bp of overlap at the merge, and a merge with 11 bp of support produces a sequence
that does not exist. On a detection platform that presents as a novel taxon — a false
positive wearing the clothes of a discovery.

**Measured on the COI run, 200 costs six times what it buys:**

| | |
|---|---|
| merges lost to MINLEN=200 | **14 536 (8.00% of all merges)** |
| …that are valid 300–330 bp amplicons | **92.82%** |
| …whose shorter mate was 175–199 bp (missed by a hair) | 50.41% |
| merges with genuinely frail overlap (<50 bp) | **2 306 (1.27% of all merges)** |

The reasoning was right about the *risk* and wrong about the *lever*. MINLEN acts on a
single read **before the overlap is known**, so it cannot distinguish a 190 bp mate that
merges with 120 bp of support from a 50 bp mate that merges with 11 bp. PEAR's
`min_length` (200) does not catch the frail case either — a 50 bp read merging with a
274 bp read yields 313 bp and sails through. The parameter that measures exactly the
quantity of concern is PEAR's **`min_overlap`, still at 10**.

Aligned at 50. Safe for the golden and verified before the change: **zero** pairs in
either frozen dataset have a read under 200 bp, so 50 and 200 produce byte-identical
output there. `KNOWN_DIVERGENCES` is now empty and a test asserts it stays empty.

**Open, and deliberately not done here: raise `min_overlap` from 10 to ~50.** It would
remove the 2 306 frail merges and keep the 12 230 good ones. It is not a free change:
the fixture has 83 merges under a 50 bp overlap and `_real` has 1 649 (2.0% of each), so
it **re-baselines the golden**. `golden.py check` would not notice — it compares a stored
run against frozen files and never re-runs — so the golden would stay green while
silently ceasing to be reproducible from the code. Green-and-stale is worse than red, and
it is the same shape as the five "looks alive, does nothing" findings. A re-freeze is a
decision of its own, and `_real` needs hardware this box does not have.

**Also on the list, not now:** a fixed MINLEN is wrong on a multi-marker platform. It
should derive from the expected amplicon length, exactly as the marker now derives from
the primer preset — same shape, same solution.

## 2026-07-24 — The same bug, fixed once already, in the CONSUMER instead of the SOURCE. Third recurrence of this pattern, and the first with documentary proof the earlier fix failed for that reason.

`min_unique_size` had four definitions. Three said 1 — `core/config.py`
(`PARAFISH_MIN_UNIQUE_SIZE`, with the measurements behind the choice), `core/models.py`
(delegating to the constant, explicitly commented as "not a literal here"), and
`.env.example`. The fourth, `ui/pages/page_setup.py`, hard-coded `value=2` into the
field, and the callback takes the field whenever it is non-empty — which it always is.
**Every UI run clustered at 2 while the project's documented, deliberate default was 1.**

**The part that makes this worth its own entry is the comment already sitting in
`ui/callbacks_setup.py`:**

> Empty control -> the configured default, not a literal typed here. This line's
> hardcoded 1 disagreed with the runner's hardcoded 2, so the policy depended on the
> entry point.

This bug was already found once, and already fixed once. The fix was applied **where the
value is read** — the callback learned to fall back to the constant — and the source that
was producing the wrong value, the page's `value=2`, was left alive. So the fallback
never fired (the field is never empty) and the same defect came back by another path.

Third recurrence of this pattern in this project (the phylo guard, the two
`latest_run_dir`, this), and the first where we can prove the earlier repair failed
*because of where it was applied*. **The transferable lesson: fixing at the point of
consumption leaves the source alive and guarantees the recurrence by a different route.
The two `latest_run_dir` were collapsed onto one primitive, not taught to agree — that is
the shape that holds.** `tests/test_setup_defaults_single_source.py` now parses the page
and fails if ANY Setup field disagrees with its `RunConfig` default, which closes the
class rather than this instance.

## 2026-07-24 — A defect class of its own: the interface teaching the OPPOSITE of the documented policy.

The "Min unique size" help text on the Setup page read:

> Singletons are overwhelmingly sequencing error and produce spurious OTUs. Set to 1 to
> keep them.

`core/config.py` says the opposite, and says it from measurement: the reads a floor of 2
uniquely discards are *not* sequencing errors — those sit within 3% of an abundant
sequence and get absorbed by clustering either way. What a floor of 2 discards is the
read that is globally unique AND >3% from everything else, **which is the exact profile
of a rare novel parasite** — the organism this platform exists to detect. The reference
run measured 1 and 2 producing identical output, with clustering time as the only real
cost.

So the field was not merely defaulted wrong. It was **arguing** for the wrong value, in
the user's own words, at the moment of choosing. A user who read that tooltip and acted
on it would raise the floor and delete the rare signal — and would believe they were
following the platform's advice while doing it. A wrong default is a bug; a wrong default
with documentation defending it is a bug that recruits the user against their own study.

Logged separately because the fix is separate: correcting `value=2` to `MIN_UNIQUE_SIZE`
does nothing about a tooltip that will still be read. Both were changed. Worth a sweep of
the other `_param_label` texts against what the code actually does — not done here.

## 2026-07-24 — Task 5 (measurement) + the bug it exposed: min_unique_size was 2 in every UI run, and the 90→112 rise is RECOVERY, not fragmentation.

**The prediction was wrong, and the reason is worth more than the prediction.** Cutting
the primers was expected to REDUCE the feature count: OTU1 and OTU2 differed only in the
primer regions and should have collapsed. Measured, OTU1 vs OTU2 have exactly **11
differences, all 11 inside the primers** (positions 2, 5, 8, 17, 20 and 344–359 — the
degenerate W/Y sites of mlCOIintF/jgHCO2198), bodies **100.00% identical**, global
identity 96.99% against a 97% threshold. They DID collapse — the new run's OTU1 *is*
that shared body, and across the 90 old OTUs this was the only pair to merge (90 → 89
distinct bodies). That is −1. The total went **+22**.

**What actually drives 90 → 112 is dereplication, and it is recovery.** Uniques fell
7 120 → 1 653 while the largest unique rose **284 → 139 938 reads**: before the cut, the
dominant sequence was shattered across thousands of "uniques" that differed only at the
primer's degenerate positions — synthesis noise, not biology. With `min_unique_size = 2`
in force, **45 365 of 179 061 filtered reads (25.3%) were discarded as singletons**;
after the cut, 8 087 (4.5%). **+37 195 reads entered clustering.** Real low-abundance
variants that had been scattered below the singleton floor now sum to ≥2 and seed their
own OTUs. The rise is not new fragmentation — it is variants the primer noise had been
destroying, coming back.

**Two hypotheses tested and REJECTED, both mine.** (1) A mixture of trimmed and
untrimmed reads causing new splitting: lengths got *more* homogeneous, not less (modal
band 90.3% → **96.40%**, and 365−313 = 52 = the two primers exactly). Untrimmed
survivors are **335 reads, 0.19%**. (2) A propagated-taxonomy estimate of host features
in the new run: I produced the number (4), then found it was an artifact of my own
fixed-offset 26/26 matching, which only reconstructs the body of OTUs in the modal band
— 21 of 112 matched, exactly the 21 in that band. **Withdrawn, not reported as a
result.** Q2 (host attribution, 61/90 before) stays OPEN: the new run was OOM-killed
during Databases and has no NCBI/BOLD taxonomy at all. It needs a run that finishes.

**The bug the measurement exposed, unlooked-for.** That `min_unique_size = 2` should not
have been there. The project decided **1** deliberately — a rare parasite is present as
few reads, and a floor of 2 deletes exactly that profile — and implemented it as
`PARAFISH_MIN_UNIQUE_SIZE` with default 1 (core/config.py), delegated properly by
core/models.py, documented as 1 in .env.example. But
**`ui/pages/page_setup.py` hard-coded `value=2` into the field**, and
`callbacks_setup.py` takes the field's value whenever it is non-empty — which it always
is. A fourth source, the only one anybody was reading. The phylo guard and the two
`latest_run_dir` again: two places defining one value, and the one nobody knew about
winning silently. The comment already sitting in callbacks_setup.py ("this line's
hardcoded 1 disagreed with the runner's hardcoded 2, so the policy depended on the entry
point") shows this exact bug was fixed once — in the callback, not in the page.

Collapsed onto `MIN_UNIQUE_SIZE`, and the field's help text — which still asserted
"singletons are overwhelmingly sequencing error" — replaced with what was actually
measured. `tests/test_setup_defaults_single_source.py` now parses the page and fails if
ANY Setup field disagrees with its RunConfig default, so the class is closed and not
just this instance.

**A second divergence, measured and NOT fixed:** `setup-trim-minlen` offers 50 where
`trimmomatic_params.minlen` is 200. Same shape, but which floor is right is a scientific
call, so it is recorded in `KNOWN_DIVERGENCES` with a test pinning both values — it
cannot drift further unnoticed, and it is on the list rather than decided by me.

**`--discard-untrimmed` enabled** (`PARAFISH_DISCARD_UNTRIMMED`, default on). The
argument is provenance, not signal: a read cutadapt could not trim cannot be shown to be
an amplicon of the targeted region, and a paper should be able to claim every feature
derives from a verified amplicon. Costs ~7.7% of INPUT reads (92.32% carried the forward
primer) but only 0.19% of what reaches clustering — expect read_tracking's input→trimmed
drop to widen by ~7.7 points, which is the flag working. Independent of
`min_unique_size = 1`: one keeps genuinely rare amplicons, the other removes material
that is not a verified amplicon at all.

**One regression caught in the making.** `--pair-filter=any` is rejected outright by
cutadapt in single-end mode, so the first version made the whole step fail and fall back
to untrimmed passthrough — primers silently KEPT, the exact opposite of the intent.
Caught by the existing `test_degenerate_primer_is_removed_from_every_read`, which failed
with `primers_removed 0.0 != 1.0`. The flag is now paired-end only.

522 tests pass (+8). Golden: no regressions on fixture, fixture_asv, real.

## 2026-07-24 — Task 1 (post-run COI series): the marker came from the FILE NAME, and "OTHER" was not inert — SILVA was cut with COI primers and the nothing it found was reported as an honest no_hits.

**The finding that matters is not the wrong marker. It is what the wrong marker made the
pipeline DO.** A COI/Leray run on samples named `aptl64b_01` resolved to marker `OTHER`,
because the marker was inferred solely from a regex over the FASTQ file stem
(`_infer_marker_from_name` looks for `coi`/`18s`/`its` tokens; that name carries none) and
the **primer preset — the one explicit statement the user makes about what was amplified —
was never consulted.** There was no preset→marker map anywhere; only the inverse
(`MARKER_DEFAULT_PRESET`), and only for 18S.

I expected `OTHER` to be inert — a missing value that would disable the marker filter.
Measured, it was the opposite. `marker_filter_applies("OTHER", Leray)` returned **True**:
`OTHER` is absent from `MARKER_LINEAGE`, so it fell past the lineage branch into the "has
primers" branch, which only asks whether two primer strings are set and cutadapt exists.
`build_marker_reference` skipped its own documented unknown-marker guard for the same
reason. **So the run took the COI Leray primers to SILVA and PR2 — small-subunit rRNA
references — and spent minutes cutting them to a "COI amplicon" that cannot exist there.**

**This is not dead code pretending to be alive.** That was the shape of the previous five
findings, and it is the wrong lens here. This is a *meaningless operation producing a
result-shaped output*: cutadapt ran, succeeded, wrote a filtered reference, and VSEARCH
loaded it and returned assignments-of-nothing. Every layer reported success. The empty
result then reached `source_status` wearing the word `no_hits` — which says "these taxa are
not in SILVA", a claim about biology, when the truth was "SILVA was asked a question it
cannot parse, using a reference we destroyed on the way in". **Had that shredded reference
returned even one spurious assignment, nothing in the run would have flagged it and I would
have believed it.** The danger was never the wasted minutes; it was that the output was
indistinguishable from a real one.

**Applicability was one-directional, and in the wrong phase.** `SOURCE_SUPPORTED_MARKERS`
held exactly `{"BOLD": {"COI"}}` — a COI source on an rRNA run was `not_applicable`, an rRNA
source on a COI run was not. And the decision lived in `_reconcile_source_status`, called at
line 1471, *after* the filter (1293/1321) and the SINTAX load. Adding SILVA/PR2 to the map
alone would have bought nothing: a post-mortem label with the minutes already spent. The
rule needed a second seat, before the work; reconciliation stays as the net for sources that
actually ran.

**Fixed:** `PRESET_MARKER` beside `PRIMER_PRESETS` (with a test that the key sets match, so a
new preset cannot be added marker-less); `resolve_run_marker()` — preset first, file names as
fallback, disagreement logged but not fatal; `marker_is_known()` treating `OTHER` and `""`
alike as "we could not tell", guarding both `marker_filter_applies` and
`build_marker_reference`; `SOURCE_SUPPORTED_MARKERS` gaining SILVA `{16S,18S}` and PR2 `{18S}`
(both shipped releases are SSU); `_marker_incompatible()` gating SILVA/PR2/BOLD before any
work.

**Two design decisions, both reruns of lessons this log already carries.**
`RunConfig.marker` is serialized but **never read back as input** — it is re-derived on every
load. A marker frozen before the preset changed would be a second source of truth
disagreeing with the first, which is the same failure the two `latest_run_dir` collapse (Task
7) was about. And `_not_applicable` emits **no `key` and no `raw`**: those are precisely the
fields `_reconcile_source_status` judges a source by, and a verdict about the *marker* must
not be re-litigated by machinery that judges *results*. Without a key the entry passes
through untouched, which is what makes it final.

**Addendum — the readout that agreed with nothing.** Checked whether Setup has a marker
control. It does not; it has a marker *display* — the `Marker` column of the sample table.
It rendered the literal **`NONE`**, because `sample.get("marker", "OTHER")` never fires its
default (the key exists holding `None`) and `str(None).upper()` is `"NONE"` — a value the
pipeline never assigns. A reader takes `NONE` for "not filled in yet", not for "resolved to
OTHER and about to act on it", which is why a COI run could be configured with nothing on
screen contradicting it. The DB-recommendation badge was worse: with samples on screen it
said *"Upload samples to see database recommendations"*, because every marker resolved to
`OTHER`, `OTHER` was discarded, and the empty set was reported as an empty upload. Both are
now fed by `resolve_run_marker` — the column shows the run's effective marker (COI, titled
with where it came from), and the badge names the marker and its recommended sources.

**Left open, deliberately.** `_infer_marker_from_name` still stands, and is still the only
source when no preset is chosen. There is no marker selector anywhere in the UI, so a run
configured without a preset falls back to a regex over file names **with no way for the user
to correct it** — the column now displays that guess honestly, but display is not control.
Logged, not fixed: it is a separate question from the one measured here.

502 tests pass (+14). Golden: no regressions on fixture, fixture_asv, real — baseline
confirmed green before the work started, so "no regressions" means something.

## 2026-07-23 — Tasks 1 & 5 (post-run UI series): dropdowns — the Task 1 diagnosis was REFUTED by the DOM; the real control is a native button.

**The refutation, stated plainly.** Task 1 "fixed" illegible dropdowns by styling react-select
classes (`.Select-*` v1, then `.Select__*` v5, then a `.Select div/span` catch-all). Measured now
against the real DOM and the installed bundle: **dcc.Dropdown in this Dash (4.1) is NOT react-select.**
It is a native `<button class="dash-dropdown">` (`dash/dcc/async-dropdown.js`, embedded
`src/components/css/dropdown.css`). Every react-select rule the theme carried matched NOTHING. The
cause named in Task 1 was CSS that never applied; the premise was inferred from a version number
("Dash 4 → react-select v5") instead of measured against the DOM; and so **two rounds of correction
could not have had any effect** — the dropdowns were only ever readable/unreadable for reasons the
theme never touched. The automated check went green because it read the FILE and confirmed the
selectors were present — never that they matched anything rendered. Fifth incarnation of
"looks alive, does nothing", this time inside the fix itself.

**The real mechanism, measured.** The button's stylesheet is keyed on `--Dash-*` design tokens whose
defaults are a LIGHT theme (`--Dash-Fill-Inverse-Strong: #fff`), and `.dash-dropdown { color: inherit }`
takes the value's colour from the ancestor — so on our dark theme the value read light-on-white, and
the contextual inheritance is what made it intermittent (the one thing Task 1 got right about the
*mechanism*, even while wrong about the culprit rule). Fixed the intended way — Via A, redefining the
tokens — not by overriding the component's own rules with `!important`:
- `:root:root { --Dash-* : <theme vars> }` — 11 tokens mapped to `--bg-*`/`--text-*`/`--blue`
  (the accent, not the bundle's purple). `:root:root` (0,2,0) out-specifies the bundle's injected
  `:root` (0,1,0) so it wins regardless of async-load order, without `!important`.
- `button.dash-dropdown { color: var(--text-0) }` — the value colour pinned explicitly (0,1,1 beats
  the component's 0,1,0), because leaving `color: inherit` is the exact fragility this bug taught.
- A small block ADDING option hover/`[aria-selected]` backgrounds — states the component ships bare;
  commented as additions, not licence to override.

**The scope, measured not guessed.** The `--Dash-*` tokens are shared by slider/datepicker, but the
theme has NO slider/datepicker rules (so no double-source), and the slider consumes only
`--Dash-Spacing` (left at its default) — so a global `:root` themes the dropdown with no surprise
elsewhere. Global was correct precisely because there was nothing to collide with.

**The guard changed, and this is the point.** The old check verified our own file. The new one
(`tests/test_dropdown_theme_tokens.py`) asserts the theme's defined `--Dash-*` set COVERS the set the
INSTALLED bundle consumes — it reads `async-dropdown.js`, so it fails the day a Dash upgrade renames
or adds a token, which is the failure mode that cost these two sessions and which no CSS-content check
could catch. But **closure is the user's visual confirmation, not this test** — green here has meant a
readable dropdown once before and been wrong.

**Task 5 folds in — no population bug.** The Results `results-stage-selector` / `results-sample-selector`
were logged as "empty / not populated". They have options (confirmed visually); they READ empty
because the selected value rendered invisible — the same contrast defect, not a missing-options bug.
That note is closed: there is no population work, only this contrast fix.

## 2026-07-23 — Task 7 (post-run UI series): one run selector — Results and the Eixo B pages can no longer pick different runs.

**Measured, not assumed.** Two functions named `latest_run_dir`: `ui.runtime_state` (Eixo B page
fallback) was hardened in Task 6 to order by run_id and skip incomplete runs; `core.run_registry`
(Results, via `run_dir_for`) still ordered by **mtime** with **no completeness gate**. Same
responsibility written independently, one hardened and one not — the phylo-guard class. And a
concrete consequence, not just debt: with an OOM orphan as the newest run, Results picked the orphan
while the Eixo B pages picked the last finished run — two different runs in ONE session, worse than
either being wrong alone. `run_registry.latest_run_dir` is called only by `run_dir_for`, whose
session-`output_dir` branch already serves the active/in-progress run, so its `latest_run_dir` branch
is purely the cross-session fallback — the same question Eixo B asks.

**Duplicate on one axis, legitimately different on another.** They shared the SELECTION question
(newest run to restore) but ran different implementations of it (the bug). They differ legitimately
on SCOPE: Results falls back across all projects (a user who never opened Setup still sees their last
run — pinned by `test_the_run_is_found_with_no_ui_state_at_all`), Eixo B is project-scoped. So the
collapse was of the selection, not the scope: a single primitive
`core.run_registry.newest_complete_run(run_dirs)` (newest complete by run_id, skips orphans) that
BOTH `run_registry.latest_run_dir` and `ui.runtime_state.latest_run_dir` call. Like
`MIN_TREE_SEQUENCES`: one definition, two call sites, not two definitions kept in sync. The third
selector, `restore_run_state` (newest NON-interrupted, includes in-progress), was left alone — it
answers a different question (resume a live session vs fall back to finished tables), a legitimate
difference, not duplication.

**Check seen to fail.** `tests/test_run_dir_single_source.py`: with an orphan newest by run_id AND
mtime, Results and Eixo B chose different runs (RED); now both pick the finished run. An older run_id
with a bumped mtime no longer beats the true newest run_id (the mtime hazard, gone). Full suite 490
passed; `golden.py check all` clean on fixture, real, fixture_asv. This closes the divergence Task 6
opened when it hardened only one of the two.

## 2026-07-23 — Task 6 (post-run UI series): Startup reconciliation — an OOM orphan must not masquerade as the active run and blank the restore.

**Numbering, so the record does not lie.** This is Task 6 of the post-run UI series
(1 dropdowns illegible, 2 stage-badge lag, 3 single Start trigger, 4 "server did not respond"
diagnosis, 5 the Results selectors reading empty, 6 this). Task 5 was gated on a human visual pass
when 6 landed; it later proved to be the SAME dropdown contrast defect as Task 1 (the selectors have
options, the value just rendered invisible), fixed in the Tasks 1 & 5 entry above — not a population
bug. Tasks were worked 1→2→3→4→6 because 5 was blocked, not skipped. (The
`Task 3` / `T0–T4` entries elsewhere in this DEVLOG belong to the EARLIER UI-unification series, a
different numbering context.)

**The problem, measured.** A pipeline runs in a thread inside the Flask process, so a run left
`stages: {"Analysis": "running"}` on disk with no `done` is one whose process the OOM killer took
mid-flight — `runs/` already held 7 such orphans. `latest_run_dir` (shared by Project and the four
Eixo B pages) picked the newest run by run_id regardless of completion, so if the newest was an
orphan the disk-restore rebuilt the Eixo B pages from its half-written records and they came back
BLANK — indistinguishable, to anyone looking, from the Eixo B restore bug having regressed.

**Two protections.** (1) Startup reconciliation (`core.run_registry.reconcile_interrupted_runs`,
wired in `app.main`, not `create_app`, so the test suite never mutates real runs): a run that claims
`running` with no live thread is marked `interrupted` (stages `running→interrupted`), deliberately
NOT `done`. (2) `latest_run_dir` returns only a COMPLETE run (`done`, not interrupted); when a
project has none it returns `None` — the pages then show their honest empty-state, because a
"marked orphan" would give the same blank with more steps (its CSVs never finished); only `None`
distinguishes "nothing here" from "broken", and that distinction is the whole task.

**The mid-task refinement — a test that resisted the wrong change.** The first cut applied a blind
`done: true` gate to `latest_run_dir`, which `restore_run_state` reused. An EXISTING test,
`test_a_run_still_in_progress_restores_config_but_no_stage_records`, went red: a genuinely live run
(no `done`, not interrupted, has a thread) must come back as config-only so a fresh tab reconnects
the poller. The fix was to REFINE the design, not weaken the test: two consumers hold opposite needs
over the same word "incomplete", so each got its own selector — `latest_run_dir` is complete-only
(finished tables for the page fallback), `restore_run_state` selects the newest NON-interrupted run
(resumes a live run, excludes a reconciled orphan). This is the exact inverse of the 401-tests
episode: there the tests shared the code's wrong belief; here a test held the line against the wrong
change. A test that stops the wrong thing is the whole return on the net.

**Proof.** `tests/test_startup_reconciliation.py`: orphan→interrupted, finished run untouched, live
run never reconciled, `latest_run_dir` prefers the complete run and skips a newer orphan (RED against
the old max-by-run_id — it returned the orphan), `None` when only orphans, failed-but-finished still
restores, in-progress resumes config-only, reconciled orphan declined. Eixo B restore against golden
intact (`test_qc/preprocess/analysis/databases_restore`, `test_project_flow` — 17 passed). Full suite
485 passed; `golden.py check all` clean on fixture, real, fixture_asv.

**Left standing, and named (Task 7).** A SECOND `latest_run_dir` in `core/run_registry.py` still
orders by mtime — the same class as the phylo guard: two functions with one responsibility written
independently, one hardened and the other not. It now has a concrete consequence, not just debt:
Results (via `run_dir_for`→mtime) can show one run while the Eixo B pages (via runtime_state→
complete-only) show another — internal divergence within a session, worse than either being wrong
alone. Task 7 collapses them to a single source, after measuring whether they truly duplicate.

## 2026-07-22 — verify.py check (f) now reads provenance: a Linnean absence that was CHOSEN is not-applicable, one that DIED still fails.

**The decision from the previous entry, implemented.** Making verify ASV-aware surfaced `(f/asv)`
failing on `fixture_asv` (SILVA+PR2-only → no Linnean backbone). The call: (f) should draw the same
line (e) already draws — an alarm that fires on an absence the user CHOSE (no Linnean source
requested) trains them to ignore it, and then it stops guarding the case it exists for, a requested
Linnean source that DIED. So exempt "never requested", keep the alarm for "requested and absent".

**The implementation trap, and where recording db_sources at init pays off.** To tell "not requested"
from "requested and died", (f) must read what was ASKED FOR, not what arrived — and that only exists
because the run's requested sources are written to `progress.json`'s `run_config` at init.
`_requested_sources(run)` reads `run_config.db_sources`; the decision is: a Linnean source
CONTRIBUTED to the consensus → PASS (unchanged); else no Linnean requested (and the request is
recorded) → SKIP (config choice, not a failure); else (requested-but-absent, OR the request is
unrecorded) → FAIL. The PASS-on-contribution comes first, so the OTU fixture — whose `run_config` is
empty in progress.json but whose consensus NCBI really fed — stays PASS without depending on the
provenance at all.

**Seen to fail on BOTH sides.** Config choice: the committed (pre-provenance) verify FALSE-ALARMS —
`(f/asv)` FAIL on `fixture_asv` (db_sources=silva|pr2); the new code SKIPs it (`config choice, not a
failure`). Death: a scratch run with `ncbi` added to db_sources but absent from the consensus →
`(f/asv)` FAIL (`requested NCBI but no Linnean source fed…`), exit 1 — the guard preserved. Both
locked in `tests/test_verify_asv.py` (config-choice→SKIP, requested-died→FAIL, contributed→PASS,
unrecorded→FAIL fallback).

`verify.py runs/_fixture_asv/current` now exits 0 (5 passed · 0 failed · 2 skipped — (f/asv)
not-applicable and (g) the OTU-only reader probe); no false alarm on the committed fixture. OTU
unchanged. Full suite 452 passed; `golden check all` clean. This truly closes the verify prerequisite
before the data phase, which will run datasets with and without Linnean sources and needs (f) telling
the two apart from the first run.

## 2026-07-22 — verify.py made ASV-aware: the blind instrument now sees an ASV run, and the blindness was seen over a real defect first.

**The prerequisite of the data phase, closed.** verify.py's seven checks all read `databases/otu/`
and `analysis/pooled/otu_table.tsv`, so on an ASV run every one SKIPped and the run exited 0 — a
blind instrument goes quiet, it does not scream, so an ASV run's "verify passed" validated nothing.
Measured on `runs/_fixture_asv/current`: 0 passed · 0 failed · **7 skipped · exit 0**. That is the
"green pretending to prove" one layer up, in the tool that judges the pipeline.

**Mapped the six before touching any (the checks fragment, they are not uniformly blind).** (a) read
totals, (c) class table, (d) banned ranks, (e) source_status, (f) Linnean backbone — all
OTU-hardcoded but fully ASV-able (`databases/asv/*` and the DADA2 abundance table exist). (b)
id_map — ASV-able too: DADA2 writes `analysis/dada2/asv_id_map.csv` with an `asv_id,size` shape. Each
now takes a `feature_set` and resolves the set's own paths through small helpers
(`_abundance_table`, `_id_map`, `_db_dir`); `main()` runs them PER SET actually present
(`_present_sets`), labelling each `(d/asv)` etc. Two seams stay OTU-only ON PURPOSE: the chimera leg
of (a) — DADA2 removes chimeras internally, there is no ASV chimera_stats — and check (g), the
Results-reader probe, which is chimera-based (`clustered_read_total`) and resolves its own paths.

**The fail-first trap, and it is subtle — honoured.** A blind verifier does not error on an ASV
defect, it stays silent, so "seen to fail" here is NOT "verify errors" — it is: plant a real defect
in an ASV run and prove the OTU-blind verifier does not catch it, then that the ASV-aware one does.
Planted an invented Linnean rank (Craniata as a Class) into a scratch `databases/asv/
taxonomy_consensus.csv`. The OLD verify (git HEAD): `(d)` SKIP, 0 failed, **exit 0 — the defect
uncaught**. The NEW verify: `(d/asv)` FAIL, **exit 1 — caught**. The instrument that verifies the
instrument, itself seen to fail before being trusted. Locked in `tests/test_verify_asv.py`.

**A real finding surfaced the moment the blindness lifted: `(f/asv)` FAILs on `fixture_asv`.** Not a
verify bug — a truth. `fixture_asv` is SILVA+PR2-only (NCBI was excluded for determinism, its remote
being non-reproducible), and SILVA/PR2 backbones are not Linnean, so Phylum..Family are genuinely
empty and (f) — "a run whose only Linnean source is absent is coherent and useless" — correctly
fires. In a PRODUCTION ASV run with a Linnean source (the COI+BOLD demo had BOLD) it would pass.
**Decision left to the domain owner, not silently coded:** whether (f) should distinguish a Linnean
source that was never SELECTED (a config choice, `skipped`/`not_applicable` in source_status) from
one that ran and DIED — the way (e) already distinguishes them. (f) is deliberately cause-agnostic
today; softening it is a real semantic change and yours to make, so it is flagged, not made.
**[RESOLVED 2026-07-22 — see entry above: (f) now reads provenance and draws exactly this line.]**

`verify.py runs/_fixture_asv/current` now exits 1 for that true reason (the tool is manual, in no
CI gate, so nothing automated regresses). OTU unchanged: `runs/_fixture/current` still 7 passed ·
exit 0. Full suite 448 passed; `golden check all` clean on all three (verify.py is not part of the
golden). This closes the last no-data item — the ruler is sharpened before the data phase begins.

## 2026-07-22 — Group 2 closed: the ASV export is contracted against the existing fixture_asv, no new data.

**The map's shortest close, done: `export/asv` now travels proven.** The databases stage was already
writing `databases/asv/{phyloseq,taxa_tables,table.biom}`, but `golden.py`'s export stage froze only
the `otu` paths, so the asv export shipped uncontracted — a change to it would have gone uncaught.
Added the six asv paths (mirroring otu) plus `databases/asv/table.biom` to the export spec, which is
now feature-set-adaptive by the same absent-guard as analysis/taxonomy: an OTU run freezes otu/*
(asv absent), an ASV run asv/* (otu absent). No aligner floor here — export is a deterministic
transform of the frozen taxonomy — so 3 features sufficed and it froze against the EXISTING
`fixture_asv`, no new dataset.

**Check seen to fail first, then the pin proven to bite.** Before: `check export --dataset
fixture_asv` SKIPped — no export layer, the asv export unprotected. After freeze: PASS. Perturbing
`golden/fixture_asv/export/.../taxa_table_genus.csv` makes `check` FAIL (`CONTENT CHANGED`, the asv
path named), restoring makes it pass — so the newly-added asv paths are really checked, where nothing
checked them before. Adding asv to the spec did NOT regress otu: `check all` PASS on fixture, real
and fixture_asv (the absent-guard skips asv paths on an otu run). Full suite 445 passed (golden data
only, no code change).

**Self-feeding order holds.** This contracts what `verify.py` ASV-aware will verify next: with the
asv export now frozen, the ASV-aware verifier (the data-phase PREREQUISITE) has a real export layer
to check. Group 2 first (contract), verify.py after (verify) — the second has its object because the
first froze it.

## 2026-07-22 — MAP (audit, no fixes): the ASV branches written but proven only for OTU, grouped by what proves each.

**What this is.** Task 3: an inventory, not a correction. Nothing is fixed here — fixing a branch
that cannot be proven is writing against a check that does not fail, which this project refuses.
The deliverable is a map so that when an ASV dataset lands the proof is mechanical: each site with
the check that fails, waiting for its artifact.

**Reconcile first — most of the old LIVING LIST (below, 2026-07-1x) is already CLOSED.** That list
predates the ASV harvest. `fixture_asv` (3 ASVs) has since exercised the analysis metrics, the
taxonomy, AND the Eixo-A/B restore + right-set of the harvest callbacks (`_render_hits`,
`_render_consensus`, `_render_occupancy`, the six Analysis residual callbacks) — those rows are
done and green against `fixture_asv`. So 3 features was enough for everything EXCEPT what needs
more. What remains is below, grouped BY WHAT PROVES IT — the requirement differs per group and they
must not be merged.

**PREREQUISITE OF THE DATA PHASE — ✅ CLOSED (2026-07-22, see entry above): `verify.py` is now
ASV-aware.** The six checks run per feature set present, an ASV defect is caught (proven over a
planted banned rank), and the fix surfaced a real finding — `(f/asv)` fails on the SILVA+PR2-only
`fixture_asv` for lack of a Linnean backbone, a decision-point flagged for the domain owner. Kept
below as the record of why it gated the phase. **The blindness it describes was real:** `verify.py`
checks `analysis/pooled/otu_table.tsv`, `otu_matrix.csv` and `databases/otu/*` — all OTU-hardcoded —
so an ASV run passes through it WITHOUT BEING VERIFIED, and reports green while proving nothing. That
is the exact class this project has chased the whole way (a checker that cannot see the thing it
claims to check: `conflict_report.csv`, `table.biom`'s date, the phylo guard) — now sitting in the
instrument itself. **Make `verify.py` ASV-aware BEFORE trusting any green from an ASV run.** Until
then, an ASV run's "verify passed" is validation against a blind instrument, worth nothing. This
gates the whole data phase; it is not scheduled among the groups below, it comes before them.

**GROUP 1 — needs a fixture ASV with ≥4 features.** `fixture_asv` has 3 ASVs; FastTree's floor is
4 (`MIN_TREE_SEQUENCES`), so the whole phylo/asv branch cannot run against it. Written, resolves
`asv`, proven for `otu` only.

| site | where | proven only for otu by | what proves the asv branch |
|---|---|---|---|
| phylo per-set writer | `ui/callbacks_progress.py:867,885` (`phylo/<key>/`) | golden phylo freezes `phylo/otu/*` (`scripts/golden.py:134–139`) | run phylo on a ≥4 ASV set → freeze `golden/fixture_asv/phylo` |
| phylo readers (per-set→flat) | `ui/feature_sets.py` `phylo_dir` / `phylo_tree_path` / `phylo_artifact` | `tests/test_phylo_golden.py` asserts `phylo_tree_path(..,"asv")` is **None** (no asv tree exists) | same ≥4 set → the None assertion flips to a real tree; extend the test with an asv reconstruction |
| phylo UI/results readers | `callbacks_phylo._phylo_kpis`/`_render_annotated_tree`/`_mapping_csv`/`_phylo_metric_paths`; `callbacks_ordination` (~145); `callbacks_results_summary._phylo_card` + `_EXPORTS` phylo rows | all hardcode/resolve `otu`; no asv phylo artifact to render | ≥4 asv phylo output on disk → the same resolvers exercised for `asv` |

**GROUP 2 — ✅ CLOSED (2026-07-22, see entry above).** The asv export is now contracted in
`golden.py` and frozen against `fixture_asv`; `check all` PASS on all three datasets, the pin proven
to bite on an asv path. Left here for the record of what it was.

| site | where | gap | what proved it |
|---|---|---|---|
| ~~export/asv~~ ✅ | `databases/asv/{phyloseq,taxa_tables,table.biom}` | was: otu-only export spec, no contract | added the six asv paths + `databases/asv/table.biom` checksum; froze `golden/fixture_asv/export` (3 features sufficed; biom `date` handled by the volatile-field guard) |

**GROUP 3 — related asv gaps with a DIFFERENT requirement — do NOT bundle under "≥4".** Listed so
the data phase does not mistake them for phylo-class items.

| site | where | needs |
|---|---|---|
| `feature_set_comparison` | `analysis/feature_set_comparison.py`; UI `databases-feature-comparison` | a run with BOTH methods in one output dir (otu AND asv) — neither `fixture` (otu-only) nor `fixture_asv` (asv-only) has both. Needs a dual-method fixture, not a bigger ASV one |
| `_render_dada2_tracking`, `_update_analysis_stats` | `ui/callbacks_analysis.py:237,527` | still store-only (blank on restore) — the deferred-whole rows from the LIVING LIST below; provable with existing `fixture_asv` once their Eixo B restore is wired (session-restore axis, not the ≥4 axis) |

(`scripts/verify.py` was here too, but it is not a feature to schedule — it is the PREREQUISITE
above: the verifier itself is OTU-blind, so it gates the data phase rather than sitting in it.)

**Consumption order when the data lands.** FIRST the prerequisite: make `verify.py` ASV-aware, or
every ASV-run green below is validated against a blind instrument. THEN Group 2 (no dataset needed —
contract the asv export against the existing `fixture_asv`, the shortest close on the map). A
≥4-feature ASV dataset closes Group 1 in one harvest (run phylo, freeze `golden/fixture_asv/phylo`,
flip the `None` assertion). Group 3's remaining items each need their own distinct artifact and
should be scheduled separately, not folded into the ≥4 event.

## 2026-07-22 — update_qc_stats restored: the QC KPI cards reconstruct the RIGHT values from the golden qc, a third restore form.

**The last store-only callback, left out by shape not by axis, now closes the session-restore
gap.** `update_qc_stats` computed its four cards — Total samples, Pass rate, Mean reads, Mean GC —
straight from `store-qc-outputs`. On a reloaded session that store is empty even though the run is
on disk, so the cards blanked to "—" while every sibling QC callback (`_render_qc_plots`,
`_render_qc_sample_table`) already reconstructed from disk through `_qc_metric_paths`. Measured,
the scope was exactly this one callback — no duplicate re-claims these outputs, so it alone blanked.

**A third restore form: semantic, not path.** The two axes restore the store's SHAPE (Axis A) or
resolve an output path (Axis B). Here neither fits: the store's per-sample dict is not on disk, but
its NUMBERS are — `qc_overview.csv` carries `total_reads` and `pct_gc` per sample, and a sample only
reaches that table if FastQC produced parseable output for it, so the run pass rate the store would
show is total/total. `_qc_stats` reads the populated store first and, when it is empty, reconstructs
the same four values from the overview via the existing `_qc_metric_paths` (with its `latest_run_dir`
fallback). One `_format_qc_stats` formats both paths, so store and restore cannot diverge.

**Strong-form check, because QC is the fifth golden layer.** `tests/test_qc_stats_restore.py`
reconstructs a run dir from the COMMITTED `golden/fixture/qc` and asserts the restored cards are the
RIGHT values — 2 samples, 100% pass, 5,000 mean reads, 54% mean GC — not merely "not blank"; and that
the restore path reproduces the store path for the same run. Seen to fail first for the documented
reason (reconstruction disabled → "restored session blanked instead of reconstructing from disk"),
green after. Full suite 445 passed; the callback smoke test confirms the new 3-arg signature fires;
`golden.py check all` clean on fixture/real/fixture_asv.

## 2026-07-22 — The phylo input floor, unified at 4: the honesty bug the phylo session flagged, now closed at its source.

**The UI guard said two sequences were enough for a tree; the aligner said four.** FastTree —
which the default HYBRID and FAST modes both route through — refuses below four tips (an unrooted
tree needs four to be non-trivial). The UI guard was `seq_count < 2`, so a 3-tip run PASSED the
guard, ran QC + preprocess + analysis, and only failed downstream inside phylo with "mínimo: 4":
a false "you may proceed" of the "server did not respond" family, revealed only at the end. The
floor of 4 was measured, not guessed — `fixture_asv` has 3 ASVs and hits it (which is also why the
phylo golden is OTU-only).

**Fixed at the source of the divergence: one floor, not two numbers written independently.**
`MIN_TREE_SEQUENCES = 4` now lives once in `phylo/__init__.py`; the tree runner
(`fasttree_runner`) and BOTH UI guards — the manual button (`callbacks_phylo`) and the pipeline
(`callbacks_progress`) — read it through a shared `phylo_input_error()` helper, so the number the
user is refused at cannot drift from the number the aligner enforces. The guard message names the
real floor (4).

**Check seen to fail first, for the documented reason.** `tests/test_phylo_input_guard.py` pins
that a 3-tip input is refused AT THE ENTRANCE naming 4, that a sufficient input passes, and that
the guard floor EQUALS the aligner floor — all red on the old code (guard 2 vs aligner 4, and no
shared source), green after. A wiring check fails if either guard site re-hardcodes its own
number. Full suite 442 passed; `golden.py check all` clean on fixture/real/fixture_asv (the frozen
layers untouched — the guard only refuses inputs the aligner would have refused anyway).

## 2026-07-22 — Phylo into the golden net: the trees ARE byte-reproducible, the freeze is OTU-only by nature, and "árvore não encontrada" was a per-set path the reader never learned.

**Phylo was the last territory without a net, and it entered on the full discipline — witnessed
failure, fix, freeze, contract test — never freeze first.** The order was kept: the reader
failure was seen with real artifacts before a line changed, the writer/reader pair was made
consistent, and only THEN was the golden frozen.

**The gating measurement: phylo trees are byte-reproducible.** Ran `run_phylo` twice on the
fixture and diffed. MAFFT and FastTree are deterministic and RAxML-NG is pinned by
`phylo_params.seed=42`, so `aligned.fasta`, `fasttree/tree.nwk`, `raxml.raxml.support`,
`bestTree`, `bestModel`, `bootstraps` and `startTree` came back byte-identical. The only churn
is `*.log` (wall-clock, thread scheduling) and `*.rba` (binary alignment cache with an embedded
timestamp) — non-content, the same trap as FastQC's zip banner and `table.biom`'s `date`, so
NOT frozen. This is what let the golden freeze the TREES themselves, not merely derived tables.

**The bug the freeze was blocked on: the tree lives under its RUNNER's subdir, never flat in
the set dir.** The pipeline writer already wrote per-set (`phylo/<set>/`), but the readers
looked for `phylo/otu/tree.nwk` and `phylo/tree.nwk` — neither of which the writer ever
produced; the real tree is at `phylo/<set>/raxml/raxml.raxml.support` (or `fasttree/tree.nwk`).
So a completed tree read as "árvore não encontrada": two mismatches at once, per-set-vs-flat AND
runner-subdir-vs-flat-in-set. Fix: `ui/feature_sets.py` gained `phylo_dir`/`phylo_tree_path`/
`phylo_artifact` (the per-set→flat ladder that knows the runner subdirs), and every reader —
KPIs, tree viewer, mapping, metrics (`callbacks_phylo`), ordination (`callbacks_ordination`),
Results `_EXPORTS` + `_phylo_card` (`callbacks_results_summary`) — now resolves through them.
The manual-run button in `callbacks_phylo`, the one writer still writing flat, now writes
per-set like the pipeline.

**OTU-only, and by structure not oversight.** The freeze is `phylo/otu/` alone because
`fixture_asv` has 3 ASVs and the aligner floor is 4 tips ("Alinhamento com apenas 3 sequências
(mínimo: 4)") — the ASV set cannot exercise phylo at all. `freeze()`'s absent-guard skips the
`asv/*` paths exactly as it skips a source a run did not produce. Frozen for `fixture` (travels
with the repo) and `real` (local-only, `golden/real` stays gitignored); the pin was proven to
bite (perturb the support tree → `CONTENT CHANGED`, restore → PASS). `tests/test_phylo_golden.py`
anchors on the committed golden subtree (clean-clone-safe): the reader resolves the per-set
support tree, and the old flat paths are asserted absent. Full suite 437 passed / 1 xfailed;
`golden.py check all` clean on fixture and real, phylo correctly SKIP on `fixture_asv`.

**Closed here. What remains is two small corrections, then product scope — not fixes:**
- **The phylo guard is a honesty bug, of the family we chased all along.** The UI tells the
  user two tips suffice (`seq_count < 2`); the aligner refuses below four. A user with three
  sequences passes the UI guard, starts phylo, and the stage blows up INSIDE with "mínimo 4" —
  after waiting through QC, preprocess and analysis. That is a false "you may proceed" revealed
  only downstream, the "server did not respond" pattern exactly. The real number is MEASURED
  this session (`fixture_asv` at 3 hits the floor of 4), not guessed. Next session: move the
  guard to 4 with a check that refuses the 3-tip case AT THE DOOR instead of at the aligner.
- **`phylo/asv/` is written but never exercised.** The per-set branch added this session was
  proven for `otu` only — no fixture ASV with ≥4 sequences exists to run it. This is the 717
  lesson pointed at this session's OWN work: a golden proves what it exercises, and this one
  does not exercise ASV. The `asv` resolver and writer path exist and look right, in the same
  state the deferred callbacks were in before their anchor — written, plausible, unmeasured.
  Recorded here with its true nature: **`phylo/asv/` written but unexercised — needs a fixture
  ASV with ≥4 sequences to prove.** It proves for free the day a larger ASV dataset lands.
- **`update_qc_stats`** — the last small, unblocked item, still open. With it and the guard
  closed, the data layer and the UI are essentially complete.
- Beyond those, product scope, not correction: the ASV living-list awaiting a real ASV dataset,
  validation on varied data, the push once the supervisors sign off, and the Applications Note
  benchmark.

## 2026-07-21 — The per-feature taxonomy summary: a NEW deliverable, single-sourced, and a cross-metric bug caught in the building.

**One row per feature, both evidence families side by side, no blended score.**
`taxonomy_summary.csv` is now written by `run_taxonomy_pipeline` next to hits/consensus: per
feature, the consensus call (LCA taxon + rank + reliability), then the best ALIGNMENT hit
(source, taxon, identity) and the best SINTAX hit (source, taxon, confidence) in SEPARATE
columns. There is deliberately NO cross-metric winner — identity and bootstrap are different
scales, and collapsing them into one "best" is the Phase-3 bug. When the two families disagree
(fixture OTU3: NCBI says *Aspergillus flavus* at 96%, SILVA says *Gadus morhua*), the reader
sees both and decides; the table never hides the disagreement behind a single number.

**Single-sourced, the Phase-2 principle: the pipeline writes the file, the UI renders THAT
file, the download hands over THAT file.** New card on the Databases page reads
`taxonomy_summary.csv` as-is (no recompute) and the download button sends the file itself
(`dcc.send_file`), so screen and export are the same bytes — and both resolve through
`_resolve_output_dir` with a `store-project` fallback, so a reconstructed session does not blank
them (the restore discipline the harvest just established, applied to a new pair at birth).

**A cross-metric bug, caught by building the fixture OTU set and reading it.** The first cut
split the two families by "has a confidence value" — and every OTU came back with
`sintax_source=NCBI, sintax_confidence=717`. That 717 is NCBI's BLAST *bitscore*: the record's
`confidence` field falls back to bitscore for aligners (`_normalize_source_row`), so keying on
its presence mislabels an alignment as a SINTAX call. The honest discriminator is the one the
`best_alignment_record` docstring already states — SINTAX never aligns, so it carries NO
identity. `best_sintax_record` now excludes any record with an identity, and NCBI/BOLD are
cleanly on the alignment side only. The bug was invisible in fixture_asv (SINTAX-only, no
aligner to mislabel) and only surfaced on an OTU set with NCBI — the measure-against-the-real-
artifact lesson again: the ASV anchor alone would have frozen the mistake.

**Contracted in the golden for all three sets, and the pin proven to bite.** Backfilled the
summary into the frozen run dirs by reconstructing records from `source_assignments.csv` (a
1:1 serialization of `records_by_asv`) and running the SAME builder — verified byte-identical
to the raw-collector path for fixture_asv, so it is the pipeline's own output, not a second
implementation. `taxonomy_summary.csv` added to `golden.py`, frozen for fixture/real (otu) and
fixture_asv (asv); perturbing a run-dir summary makes `check` FAIL (`CONTENT CHANGED`),
restoring makes it pass. Strong proof in `tests/test_taxonomy_summary.py`: records collected
from the raw SINTAX results by the pipeline's own collectors, run through the builder,
reproduce `golden/fixture_asv`'s summary byte-for-byte. Full suite 429 passed / 3 skipped;
`golden.py check all` clean on all three.

**Harvest and derivation done. Remaining: phylo (no-network writer, check-fails-first
discipline) and `update_qc_stats`.**

## 2026-07-21 — ASV harvest CLOSED. The conclusion is about the list, not the code: every site shrank when measured, and always for the same reason.

**Four sites of blocked ASV work were on the registry; measured against the fixture_asv
anchor, all four came in smaller than filed — and the cause was identical each time.** The
entries were written before the anchor existed, so each was a guess at scope with no artifact
to check it against. The guess ran wide every time:

- **Analysis residual** — filed "six residual callbacks", was ONE resolver swap. Behind the
  six sat two resolvers; only `_analysis_metric_paths` carried the `otu` assumption,
  `_find_abundance_matrix` was already set-adaptive.
- **Occupancy** — filed deferred-whole "restore + right-set proof", was restore only. The
  selector already passed its value straight to `analysis_metrics_paths`; the right-set half
  was never broken.
- **Hits** — filed as an Eixo A right-set selector problem, was the Eixo B restore only.
  `_load_hits_frame` already threaded the set through `taxonomy_file`.
- **Consensus** — same as hits: filed right-set, was restore only.

**The through-line: a list written without the artifact that could verify it will overstate,
because it is recording fear, not measurement.** Three of the four had a set-honoring path
that was already correct — the registry called it broken because, blind, "the selector might
be wrong" is the safe thing to write down. Only the anchor could turn "might" into "is" or
"isn't", and every time it said "isn't". This is the same lesson as the fixtures (author's
fixture inherits author's belief) pointed at a work list instead of a test: the estimate and
the check must come from the same measured artifact, or the estimate is just the belief.

**So the discipline for the registry, recorded as the harvest's conclusion:** a list of
blocked work exists to keep the SITES in view — so none is silently dropped — and never as an
estimate of the effort each will take. Effort is unknowable until there is an anchor to
measure against; until then the only honest field is "blocked, and why", not "how big". Read
the registry as a map of where to look, not as a forecast.

**What the harvest delivered.** `fixture_asv` is now a five-layer anchor (qc, preprocess,
analysis, dada2, taxonomy); the six Analysis callbacks, occupancy, hits and consensus all
reconstruct from disk and honor the selected set, each proven byte-identical to golden with
its wiring seen to fail first. `_download_hits` got the single-source fix alongside hits — the
export reads the same reconstructed dir as the screen, so the two cannot diverge on a restored
session (the Phase-2 two-ledgers bug in its UI form).

**Not in the harvest — the two that remain, with their disciplines:**
- **The derivation** (next session) — one row per feature, NEW build not a restore. Design and
  site fixed: a file WRITTEN BY THE PIPELINE, contracted in the golden, UI rendering that file
  rather than recomputing. Columns: consensus taxon (LCA), best alignment hit (source, taxon,
  identity), best SINTAX hit (source, taxon, confidence) — SEPARATE, no cross-metric winner.
- **Phylo** — a no-network writer, and it gets the full new discipline: golden `check` fails →
  change the WRITER → check passes → only THEN freeze. Never freeze first.
- **`update_qc_stats`** — still open.

## 2026-07-17 — ASV harvest 3/3, step 2: hits/consensus restored — and, measured, the scope was restore only.

**The last two deferred callbacks, and the same finding as occupancy: the right-set half was
never broken.** `_render_hits` and `_render_consensus` were the two SELECTOR callbacks
`test_databases_restore.py` left out — proving they reconstruct the RIGHT set was the Eixo A
"otu-hardcoded" problem, unprovable without an ASV golden. With `golden/fixture_asv/taxonomy`
now frozen (step 1), I could finally measure, and the set-honoring was already correct:
`_load_hits_frame` / `_find_consensus_csv` thread `feature_set` straight through
`taxonomy_file`, which resolves `databases/<set>/…` and returns nothing when the set is absent.
The only real gap was the Eixo B restore — both callbacks read `output_dir` off the per-session
`store-run-config`, so a fresh tab blanked them. So the fix is the plain restore wiring, nothing
about the selector: `_resolve_output_dir(run_config_data, project_data)` + a `store-project`
State on each, and the same on `_download_hits` so the export cannot diverge from the screen on
a reconstructed session.

**Proven, strong form.** `tests/test_databases_hits_consensus_asv.py`: the reconstructed run
dir feeds the selected set's hits and consensus tables — asv byte-identical to
`golden/fixture_asv` (hits and consensus), and `otu` on an asv-only run yields nothing (the set
is honored, not defaulted). Wiring seen to fail first (both render callbacks did not declare
`store-project`), passing after; the two logic proofs pass either way, which is exactly the
occupancy result — the resolver was sound, only the callback could not see the project. Full
suite 423 passed / 3 skipped; `golden.py check all` clean on fixture, real and fixture_asv.

**Harvest done; the derivation is next, and it is a NEW build, not a restore.** One row per
feature carrying the consensus taxon (LCA), the best alignment hit and the best SINTAX hit in
SEPARATE columns — no cross-metric winner (the Phase-3 bug). Open decision, filed for its own
session: it must live as a file WRITTEN BY THE PIPELINE (contractable in the golden like
`taxonomy_hits.csv`, and the file is what the supervisors take), with the UI RENDERING that same
file rather than recomputing — one source, two consumers, the Phase-2 single-source-of-abundance
principle, so the screen table cannot diverge from the download.

## 2026-07-17 — ASV harvest 3/3, step 1: fixture_asv taxonomy frozen (SILVA+PR2), and taxonomy_hits contracted.

**Ran fixture_asv through databases and it fits the box.** qc→preprocess→analysis(dada2)→
databases with SILVA+PR2, ~37 s, no OOM (551 MB free after) — the marker-filtered references
(25 MB / 76 MB, already cached on disk) keep SINTAX inside RAM even at 2.8 GB. The stage wrote
`databases/asv/` with the full taxonomy contract, including `taxonomy_hits.csv`.

**Determinism gate before freezing, and it passed.** Two runs, byte-identical across all six
taxonomy contract files (consensus, source_assignments, source_status, conflict_report,
taxonomy_hits, taxonomic_resolution). SINTAX is seeded, so ASV taxonomy is a contract, not a
sample — safe to freeze, unlike the unseeded case that started this whole discipline.

**Sources: SILVA+PR2 only (option A), and the asymmetry is recorded as debt, not precedent.**
The OTU golden carries NCBI, whose public remote is not reproducible; it was frozen before the
QC/DADA2 determinism gates sharpened the rule, and survives only because golden checks a fixed
run dir. Freezing the ASV set clean (deterministic sources) rather than matching OTU's NCBI is
the disciplined choice; the `run_config` in each run's provenance makes the source difference
legible between the two. Noted in `golden.py` at the taxonomy spec.

**Deliverable contracted.** `taxonomy_hits.csv` is now frozen for both feature sets (otu in
fixture/real, asv in fixture_asv) — it is the supervisors' evidence table, so it enters the
net. `golden.py check all` clean: fixture, real (six layers) and fixture_asv (analysis +
taxonomy).

**Next in the harvest:** restore hits/consensus (Eixo B, measured = restore only), then build
the one-row-per-feature derivation over both sets — with the metric design decided: NO
cross-metric winner (that is the Phase-3 bug again). One row per feature carrying the consensus
taxon (LCA), the best alignment hit (source, taxon, identity) and the best SINTAX hit (source,
taxon, confidence) in SEPARATE columns, so a reader sees the call and whether the two evidence
families agree.

## 2026-07-17 — BOLD's positive path is not untested — it is BROKEN against the real output format. Reported, not fixed.

**The cheap pre-meeting probe found more than "never ran".** BOLD has produced zero
attributions in this project (COI-vs-18S marker mismatch every run), so `_read_id_engine_hits`
and `select_bold_top_hits` — the parse of a positive BOLDigger3 result — never executed, in
production or in any test (the existing BOLD tests feed already-parsed dicts and stub files,
skipping both). The task was to test them against a synthetic BOLDigger3-format result.

**Validating the synthetic against the REAL schema is what caught it.** The first synthetic
DuckDB used capitalised rank columns (`Phylum`..`Species`) and passed — but that was the code
agreeing with itself. The real BOLDigger3 3.x `id_engine_results` table (checked in
`golden_baseline/.../sanitized_input.duckdb`) emits the ranks in **lowercase**
(`phylum`..`species`). `select_bold_top_hits` looks them up via
`_BOLD_RANK_ORDER = ["Phylum".."Species"]`, finds none, and **blanks every rank** — a 99.5%
hit comes back with its `pct_identity` (read correctly, lowercase) but no taxonomy at all.
Verified directly: all six ranks empty for a 99.5% synthetic hit.

**So "BOLD functional" is broken, not merely unproven.** A real positive BOLD identification
would be silently stripped of its taxonomy. `_read_id_engine_hits` (the store reader) is
sound; the rank-collapse is the break, and it looks like a column-case fix — but it must be
proven against the real format, not a hand-capitalised fixture.

**Not fixed, per instruction.** `tests/test_bold_positive_parsing.py` records both halves
against the real lowercase schema: the reader passes; the rank-collapse is
`expectedFailure` (xfail — suite stays green) so the bug is pinned in code and flips to an
unexpected success the day it is fixed. For the supervisor conversation: the honest line is
"the reader works; the taxonomy parse is broken by a column-case mismatch and needs a fix
plus a real COI demonstration" — not "the path never ran".

**Separately, on the deliverable table:** `_download_hits` exports `{feature_set}_best_hits.csv`,
but it is the ranked EVIDENCE table filtered by the UI controls ("export what you're looking
at"), NOT a one-row-per-feature best-hit derivation. The name oversells it. The evidence
table is real and correct (per-feature × per-source, with `identity` and `confidence` in
SEPARATE columns — the Phase-3 metric distinction already honoured); the curated
best-per-feature derivation for supervisors does NOT exist yet and is a build, not a rename.

## 2026-07-17 — ASV harvest 2/3: occupancy's scope, measured, was restore only — not restore + set-proof.

**The registry is hypothesis, not inventory — occupancy narrowed too.** It was filed
deferred-whole as "Eixo B restore + right-set proof". Measured against the anchor, the
right-set half was never broken: `_render_occupancy` already passes its
`analysis-occupancy-set` selector value straight to `analysis_metrics_paths(dir, set)`. The
only real gap was the Eixo B restore — it read `output_dir` off the session store with no
`store-project` fallback, blanking on a fresh tab. So the fix is the plain restore wiring the
other pages got, nothing about the selector.

**Occupancy had no contract at all, so it got one.** `feature_occupancy.csv` was frozen for
no dataset — occupancy was never in the net. Un-deferring it, I froze `feature_occupancy.csv`
(otu for fixture/real, asv for fixture_asv) via the same absent-guarded analysis spec, so the
proof can byte-match and the metric now regresses like the rest.

**Proven.** `tests/test_analysis_occupancy_asv.py`: the reconstructed run dir feeds the
selected set's occupancy — asv byte-identical to `golden/fixture_asv`, and `otu` on an
asv-only run yields nothing (the set is honored, not defaulted). Wiring seen to fail first
(callback did not declare `store-project`), passing after. `golden.py check` clean on
fixture, real (five layers) and fixture_asv.

**Harvest remaining:** `_render_hits` + `_render_consensus` — the last, and the only one that
still needs a NEW frozen layer: run fixture_asv through databases, gate taxonomy determinism
(SINTAX seeded; decide sources, likely the OTU set for comparability, made visible by the
run_config in provenance), freeze `fixture_asv/taxonomy`, then prove. Its own session.

## 2026-07-17 — ASV harvest 1/3: the Analysis residual was ONE resolver, not six callbacks.

**The registry said "six residual callbacks"; the anchor let me check, and it was narrower.**
Behind the six selector-free Analysis callbacks sit two resolvers, and only one carried the
`"otu"` assumption:

- **Resolver A — `_analysis_metric_paths`** — hardcoded `"otu"`: on the fixture_asv run it
  returned `{}` (metrics/otu absent), blanking alpha / plots / upset. Fixed: resolve the
  feature set from what the run produced (`resolve_feature_set(None, outputs, output_dir)` —
  recorded sets, or a disk scan when the store is empty), then read `metrics/<set>/`.
- **Resolver B — `_find_abundance_matrix`** — was ALREADY set-adaptive by its fallback order
  (otu_matrix → otu_table → asv_abundance), so otu-table / sample-table already rendered ASV.
  Locked with a test, not changed — it was never otu-hardcoded. Claiming otherwise and
  "fixing" it would have been a change with no failing check under it.

So the residual is a one-line resolution swap on resolver A, and `chimera_stats` stays
OTU-only by nature (DADA2 removes chimeras internally) — resolver A correctly yields no
`chimera_stats_csv` on an ASV run, so `_render_chimera_stats` shows empty, which is right.

**Proven both directions against the anchor.** `tests/test_analysis_asv_residual.py`:
resolver A reconstructs `metrics/asv/alpha_diversity.csv` and `rank_abundance.csv` byte-for-
byte against `golden/fixture_asv` (seen failing as `{}` before the fix); resolver B locked to
`asv_abundance.csv`. The OTU restore stays green — `resolve_feature_set` returns `otu` for an
OTU run, so `test_analysis_restore.py` is unchanged. `golden.py check` clean on fixture, real
(five layers) and fixture_asv.

**Harvest remaining:** `_render_occupancy` (deferred-whole, unblocked by this same anchor —
next), then `_render_hits` + `_render_consensus` (need a frozen fixture_asv *taxonomy* layer,
their own session with its determinism gate on the taxonomy sources).

## 2026-07-17 — ASV golden frozen (fixture_asv). The harvest's anchor exists; DADA2 is deterministic, so it was safe to freeze.

**Determinism gate first, before freezing anything — the SINTAX lesson.** Two DADA2 runs
over the fixture produced byte-identical contract files (`asv_table.tsv`, `asv_sequences.fasta`,
`metrics/asv/alpha_diversity.csv`, `rank_abundance.csv`). DADA2's denoising is deterministic
here — no unseeded distribution to accidentally enshrine — so a frozen ASV golden is a
contract, not one sample of many. (If it had differed, freezing would have been blocked
until seeded, exactly as SINTAX was.)

**One spec, either method — `golden.py` made feature-set-adaptive by the absent-guard it
already had.** `freeze()` skips files that do not exist, so the `analysis` stage now lists
BOTH the otu/* and the asv/* contract paths: a VSEARCH run freezes the former (asv absent),
a DADA2 run the latter (otu absent). No branching, no second stage definition. DADA2 removes
chimeras internally, so there is deliberately no asv `chimera_stats` to mirror.

**Frozen as its own dataset, `fixture_asv`, and committed** (`.gitignore` un-ignores
`golden/fixture_asv/` beside `golden/fixture/`) so the ASV contract travels with the repo —
the feature-set-resolved readers can be checked on a clean clone, same reason `golden/fixture/`
is in. `golden.py check analysis --dataset fixture_asv` PASSES; the otu datasets (fixture,
real, five layers) are untouched — `no regressions`.

**This is the anchor, not the harvest.** With the ASV contract frozen, the next sessions are
pure consumers, each provable against it: (1) the six Analysis *residual* callbacks — swap
the `"otu"` constant for the resolved set, check they reconstruct `metrics/asv/*` on the
fixture_asv run; (2) the *deferred-whole* selector callbacks (`_render_hits`,
`_render_consensus`, `_render_occupancy`) — Eixo B restore + right-set proof. One per session,
each with its own failing-first check, all before phylo.

## 2026-07-17 — Viability probe: both unblocks run at fixture scale. The phylo finding is wider than the old note said.

**The question was never "can the 2.8 GB box do a real DADA2/phylo run" — it is "can it do
one over the fixture", and the answer to both is yes.** An ASV registry harvest needs *any*
ASV artifact, not a deep one; a phylo contract needs a small tree as well as a large one.

- **DADA2 on the fixture: runs.** Full qc→preprocess→analysis(dada2), `overall=success`,
  ~36 s, no OOM (516 MB free after). It already writes the per-set layout: `analysis/dada2/*`
  and **`analysis/metrics/asv/`** — the exact directory `analysis_metrics_paths(dir, "asv")`
  expects. So an ASV golden fits the existing pipeline with NO writer change. (One env note:
  Trimmomatic needs `TRIMMOMATIC_JAR_PATH` pointed at the repo jar off-Docker; not a blocker.)
- **Phylo on the fixture: runs.** MAFFT+FastTree+RAxML over the 4 fixture OTUs in 0.4 s. Cost
  in phylo is scale (thousands of sequences), and scale is not what the contract needs.

**The phylo finding corrects the registry note, and it is worse than "flat vs per-set".**
The writer puts the tree at `phylo/fasttree/tree.nwk` (flat, by tool). The disk reader
(`callbacks_phylo.py:530`) looks at `phylo/otu/tree.nwk` then `phylo/tree.nwk` — **neither
matches**. It is not that the reader wants per-set and the writer is flat; it is that the
tree's location diverges outright, and the writer has NO feature-set namespace at all, so an
OTU tree and an ASV tree would collide in one directory.

**Layout decided (for the later phylo session): per-set `phylo/<set>/`, tree at
`phylo/<set>/tree.nwk`, reader per-set-then-flat** (the `feature_sets.py` idiom, proven three
times). It aligns writer with what `:530` already reaches for, matches
`analysis/metrics/<set>/` and `databases/<set>/`, and dissolves the OTU/ASV collision. The
flat fallback is required: runs already on disk have flat `phylo/` and would stop rendering
without it.

**Order and a new discipline.** ASV harvest first — mechanical: DADA2 already writes
`metrics/asv/`, so it is freeze-and-prove with no pipeline-code change. Phylo after: it needs
a WRITER change, and phylo has no golden — the first write-side change with no net under it.
Mandatory order there: check seen failing (flat writer + per-set reader → tree not found) →
change the writer → check passes → only THEN freeze the golden over the new layout. Freezing
before deciding would enshrine the wrong thing — the QC lesson, restated.

## 2026-07-17 — UI Eixo B, page 4/4: QC reconstructs from disk. Eixo B complete across every page that had the bug.

**The last page, and the cleanest.** QC is pre-clustering and per-sample — no feature-set
concept, so no Eixo A entanglement at all; the whole disk-reconstructing surface is provable
now that the QC golden exists. Same Preprocess shape: `_qc_metric_paths` gains `project_data`
and a `latest_run_dir` fallback; two render callbacks now declare `store-project` —
`_render_qc_plots` (`qc_length_distribution.csv` / `qc_per_base_quality.csv`) and
`_render_qc_sample_table` (`qc_overview.csv`).

**One callback left out, by shape not by axis.** `update_qc_stats` is store-only by design:
it computes from the per-sample outputs dict in the store, not from a metrics CSV, and its
numbers (pass/fail, mean reads, mean GC) would need a semantic reconstruction from
`qc_overview.csv` / `qc_module_status.csv` — a third restore shape, not the metrics-path
read. Deferred and noted here (NOT in the ASV registry — QC has no feature sets; this is a
plain follow-up, unblocked, whenever). The stat cards stay blank on a cold load; the plots
and the sample table — the page's substance — reconstruct.

**Strong form.** `tests/test_qc_restore.py`: empty session stores + `store-project`, the
mechanism reconstructs and the bytes match golden for `qc_length_distribution.csv` and
`qc_overview.csv`; plus the wiring assertion. Seen failing, passing after. `golden.py check
all` fixture+real: five layers, `no regressions`.

**Eixo B is now complete** on every page that carried the session-restore bug —
Preprocess, Databases (selector-free), Analysis (selector-free), QC. What remains is not
Eixo B: the ASV registry (Eixo A, one harvest under an ASV golden) and phylo (no contract
yet — the last hole in the net).

## 2026-07-17 — QC golden frozen. The "degenerate length histogram" was never a defect — and the net now covers a fifth stage.

**The remembered defect does not exist, so there was nothing to freeze around.** An early
audit filed `qc_length_distribution.csv` as degenerate — "all reads in one bin at 256.0,
`mean_length=256` a placeholder, not measured" — and it sat on the cosmetic list. Before
freezing QC I checked it against a real run, because freezing a known-bad output would
enshrine the defect as contract — exactly what Phase 3 refused when NCBI died. It turned
out there is no defect:

- The raw reads ARE uniformly 256 bp — `5000/5000` in the fixture, `101190/101190` in the
  real input, measured straight off the FASTQ. FastQC itself reports a single bin
  (`Sequence Length Distribution: 256 → count`) and `Basic Statistics: Sequence length 256`
  (a single value, not a range). The one-bin histogram is the *truthful* picture of
  fixed-length reads, not a broken parse.
- `mean_length` is measured, not hardcoded: `qc/fastqc_parser.py` computes
  `Σ(length·count)/Σcount` over the parsed module, with `_midpoint` expanding any `"a-b"`
  range. The refactor that introduced `qc_metrics.py` / `fastqc_parser.py` rewrote the old
  parser and fixed this in passing; the audit note predates it.

So the planned sequence (fix → freeze → Eixo B) collapsed: no fix, freeze directly.

**What was frozen, and what was deliberately not.** QC is the first per-sample stage to
enter the net, but its CONTRACT is still per-stage — the four `qc/metrics/*.csv` the reader
consumes — so it fits `golden.py` unchanged. The per-sample FastQC `.zip`/`.html` are NOT
frozen: they carry archive timestamps and a version banner (non-content that churns every
run, the table.biom trap) and are intermediate, not contract. No QC metrics CSV embeds a
date, so the four are byte-stable. `golden.py check all` now runs five layers, fixture and
real, all PASS.

**One hole left in the net, not two.** QC was the only pipeline stage with no contract at
all — the same condition that blocks phylo. Freezing it leaves phylo as the single
remaining gap. The QC Eixo B restore now has a golden to prove against, next.

## LIVING LIST — Eixo A sites pending an ASV golden (append as each page leaves selector callbacks behind)

**What this is.** Every Eixo B page fixes the session-restore of its *selector-free*
callbacks (provable against the OTU golden) and LEAVES its *selector-driven* callbacks —
the ones whose correctness depends on the resolved feature set — because proving they
reconstruct the RIGHT set is the "otu-hardcoded" Eixo A problem, unprovable without an ASV
golden. Rather than hunt those sites page by page on the day the ASV golden lands, they are
collected here as they are left behind, so that day is a harvest, not an excavation.

**The unblock is a single event:** generate a complete ASV/DADA2 golden run (itself a
viability question on the 2.8 GB box — if it will not run, that is a finding about the
machine, not the code). Once it exists, every row below gets a golden-anchored check.

**Two categories, and they must not be merged — the ASV-day work on each is different:**

**A. Deferred-whole** — NOT touched this session, because the proof path runs through the
feature-set selector or through ASV-only territory. These still need their Eixo B restore
(the `store-project` + `_resolve_output_dir` wiring), deferred WITH the set concern rather
than applied unproven. ASV-day work: BOTH session-restore AND set-resolution, at once.

| left by | callback | reads | why whole-deferred |
|---|---|---|---|
| Databases (2/4) | `_render_hits` | best-hits table, per selected set | `databases-feature-set` drives it — right-set proof needs ASV golden |
| Databases (2/4) | `_render_consensus` | `taxonomy_consensus.csv`, per selected set | consensus is per-set, selector-driven |
| Analysis (3/4) | `_render_occupancy` | `occupancy_csv`, per selected set | the ONE feature-set selector (`analysis-occupancy-set`) |
| Analysis (3/4) | `_render_dada2_tracking` | DADA2 read-tracking | ASV-only, store-only, no `output_dir` path to reconstruct |
| Analysis (3/4) | `_update_analysis_stats` | store-recorded FASTA/consensus paths | distinct restore shape (file-paths-in-store, not metrics dir) + mixed OTU/ASV content — reapproach whole with ASV golden |

**B. Residual** — Eixo B restore DONE and PROVEN this session, but the code carries an
`"otu"`-assuming resolver (`_analysis_metric_paths`), so on an ASV run these reconstruct OTU
rather than ASV. This is a latent Eixo A bug inside Eixo-B-complete code — the restore is
right and green. ASV-day work: only swap the hardcoded `"otu"` for the resolved set and
confirm the already-proven restore stays green. Less work, different in kind.

| left by | callbacks | residual |
|---|---|---|
| Analysis (3/4) | `_render_chimera_stats`, `_render_alpha_diversity`, `_render_otu_table`, `_render_analysis_sample_table`, `_render_analysis_plots`, `_render_upset` | reconstruct via `_analysis_metric_paths` / `_find_abundance_matrix`, both OTU-first; swap `"otu"` → resolved set when the ASV golden exists |

## 2026-07-17 — UI Eixo B, page 3/4: Analysis reconstructs from disk. Two categories opened in the registry; QC still last.

**The complex page, the finer line — drawn before touching it.** Analysis is the widest
render surface and the first where the provable/deferred line was not obvious. Mapped it
first, as instructed. Two resolvers feed the selector-free callbacks, and both got the disk
fallback: `_analysis_metric_paths` (otu-assuming metrics helper) and `_find_abundance_matrix`
(otu-first abundance table). Added a local `_resolve_output_dir` (session store → run_id-
ordered `latest_run_dir`) and wired `store-project` into six callbacks: chimera, alpha,
otu-table, sample-table, plots, upset.

**A presentation selector is not a feature-set selector.** `_render_otu_table` carries
`analysis-abundance-mode` (relative/absolute) — an orthogonal presentation axis, not
OTU-vs-ASV. Classifying it by the word "selector" would have deferred a callback that is
fully provable. The only feature-set selector on the page is `analysis-occupancy-set`
(`_render_occupancy`), and that one is deferred. Distinguishing the selector that picks a
feature set (blocked) from the one that picks a display mode (harmless) is what avoids both
the false-block and the false-advance.

**The registry now carries two categories, kept apart on purpose** (see LIVING LIST above):
- *Deferred-whole* — Eixo B restore NOT done: `_render_occupancy` (feature-set selector),
  `_render_dada2_tracking` (ASV-only, no disk path), `_update_analysis_stats` (distinct
  restore shape — store-recorded FASTA paths, not a metrics dir — plus mixed OTU/ASV; would
  cost a third restore pattern for half a callback, so deferred whole).
- *Residual* — the six above: Eixo B DONE and proven, carrying an `"otu"` constant. ASV-day
  work is a constant-swap, not a restore. Merging the two categories would, on ASV-golden
  day, hide which callbacks need a from-scratch restore and which need one line changed.

**Strong form, three frozen anchors.** `tests/test_analysis_restore.py`: with empty session
stores and only `store-project`, the mechanism reconstructs and the bytes match golden for
`chimera_stats.csv`, `alpha_diversity.csv` and `otu_matrix.csv`; plus a wiring assertion over
all six against the real dependency graph. Seen failing on the tree, passing after.
`golden.py check all` fixture+real: `no regressions`. Suite `415 passed` (412 + 3). The
callback-smoke test firing every rewired callback proves the six new signatures do not break.

## 2026-07-17 — UI Eixo B, page 2/4: Databases reconstructs from disk — the selector-free tabs, kept clear of Eixo A.

**Same class as Preprocess, wider surface.** The Databases render callbacks all read
`output_dir` off `store-run-config` (session) and blank on a fresh tab. Fixed through one
shared helper, `_resolve_output_dir(run_config_data, project_data)` — session store first,
then `latest_run_dir` — used by `_database_metric_paths` and by the inline `run_output_dir`
sites. Seven selector-free render callbacks now declare `store-project`: resolution,
confidence, conflicts, contaminants, prevalence, feature-comparison, db-visuals.

**Deliberately split along the axis boundary.** The two callbacks that carry the
`databases-feature-set` widget — `_render_hits`, `_render_consensus` — are LEFT OUT. Their
Eixo B restore is the same wiring, but proving they reconstruct the *right* set is the
"otu-hardcoded" Eixo A problem, unprovable without an ASV golden. Keeping them out is what
lets this session's check touch only session-restore and never the selector. The proof
path (`_database_metric_paths` / `_taxonomy_file` at the plain `otu` default) has the
widget nowhere in it.

**Strong form, three frozen anchors.** `tests/test_databases_restore.py`: with empty
session stores and only `store-project`, the shared mechanism reconstructs and the bytes
match golden for `taxonomic_resolution.csv` (resolution) and `conflict_report.csv`
(conflicts); `source_assignments.csv` (confidence) is the third frozen table the same path
serves. Plus a wiring assertion over all seven callbacks against the real dependency graph.
Seen failing on the tree (blank `{}`, `store-project` absent), passing after.

**The shared hardening paid out, as predicted.** This page inherited the run_id-ordered
`latest_run_dir` for free — no new selection logic, no re-proof of the mtime trap; the
Preprocess session already closed it under all four pages.

`golden.py check all` fixture+real: `no regressions`. Suite `412 passed` (409 + 3). The
callback-smoke test firing every rewired callback with empty inputs is what proves the
seven new signatures do not `TypeError` in a browser.

## 2026-07-17 — UI Eixo B, page 1/4: Preprocess reconstructs from disk on a fresh session. QC held for last.

**The bug, in one sentence a user would recognise:** close the app, reopen tomorrow, land
on Preprocess — blank, though the run is on disk. The stage stores are per-session; the
three render callbacks read `output_dir` straight off `store-run-config` (session) and
returned `{}` when it was empty. This is the second half of the Results-page bug — the
disk path vs the store path — and the OTU golden already pins the numbers the page must
show, so the check could be the STRONG form, not merely "not blank".

**Two gaps, each seen to fail first.** LOGIC: `_preprocess_metric_paths` had no way to
reach disk — no `project_data`, no `latest_run_dir` fallback. WIRING: the three render
callbacks did not even declare `store-project`, so the callback could not reach the run
on disk once the helper could. Confirmed the wiring gap by reading the deps *before*
touching logic — the fix is a wiring change (add `State("store-project")` to three
callbacks), not internal logic, exactly the distinction that decides whether a fix is
four lines or a callback-graph change.

**The fix.** `_preprocess_metric_paths(preprocess_data, run_config_data, project_data)`
falls back to `latest_run_dir(project_data)` when the session `output_dir` is gone — the
disk path the Results page already takes. Three render callbacks now pass `store-project`.

**Discipline.** `tests/test_preprocess_restore.py`: strong-form (reconstructs AND the
`read_tracking.csv` bytes equal the frozen golden), plus a wiring assertion against the
real dependency graph, plus the adversarial run-selection case (below). Seen failing on
the tree (blank `{}`, `store-project` absent from deps), passing after. `golden.py check
all` fixture+real: `no regressions` on all four layers. Suite `409 passed` (406 + 3).

**Order for the remaining three, decided.** Databases → Analysis → **QC last**. QC last
because it has no golden, so its check can only be the weak form ("not blank") unless a
QC golden is generated or QC is treated apart (it is per-sample, not per-feature-set) —
a three-way decision that is easier once the Eixo B pattern is proven on two golden-backed
pages. Databases before Analysis to keep the axes clean: Analysis is where Eixo B meets
Eixo A (the feature-set-selector widgets are "otu-hardcoded" territory, unprovable without
an ASV golden), and Databases is the more contained Eixo B page.

## 2026-07-17 — Hardening: latest_run_dir picks the newest run_id, not the newest mtime. The restore's silent-corruption source, removed.

**Filed on its own, not under the Preprocess page, because it is shared and outlives that
page.** `ui.runtime_state.latest_run_dir` is how a fresh session finds the project's run
on disk — it feeds Preprocess's restore, and equally the Project-page restore
(`restore_run_state` → `callbacks_project`), and it will feed the other three Eixo B pages.

**What was wrong.** It selected `max(runs, key=mtime)`. mtime is not identity: a re-freeze,
a copy, or an editor touching a file inside an OLD run bumps that run's mtime above the
genuine newest, and the restore then silently rebuilds the page from the stale/junk run.
This is the bug that does not appear until the day you copy an old run — and then it is a
page reconstructed from the wrong data, green, with no error.

**The fix removes the source rather than patching the symptom.** Select `max(runs,
key=name)`. The directory name is the run_id — the UTC creation timestamp (`timestamp_slug`)
— which is immutable and whose lexical order IS chronological. There is nothing left to
corrupt.

**Proven, not asserted.** A decoy run with an OLDER run_id but a deliberately newer mtime:
mtime-max picked the decoy's garbage before the fix, the real run after. Seen to fail, seen
to stop. The run-listing that surfaced this (ten dirs under `runs/`, several test/junk)
also confirmed the OTHER half is safe by construction: `latest_run_dir` is scoped to the
project's own slug, so cross-project orphans cannot be picked at all — the check pins that
too, with a cross-project decoy that is correctly ignored.

## 2026-07-17 — Gate CLOSED: a run now records its own config in progress.json. Two config-different runs are distinguishable.

**The mandatory pre-UI gate is closed.** The open item from the fixture entry below
("the run should write its own `run_config`") is done, and it was done first because
until it closed, a golden diff could not tell a config change from a regression.

**The contract, mapped before touching code.** `progress.json` is written at init in
`_run_pipeline_thread` (`ui/callbacks_progress.py`) and carried only `run_id`, per-stage
records and logs. `serialize_run_config` existed in `ui/runtime_state.py` and captured
the config, but **was never called from the thread nor written to disk.** So the run
recorded status but not recipe.

**The fix (`ui/callbacks_progress.py`, +36).** A new `_run_provenance()` writes a
`run_config` block into progress.json **at init** — before any stage — so it survives a
run whose tools never ran (same principle as the init write itself): `db_sources`,
`min_unique_size`, `analysis_method`, `branch`, `markers`, `primers` (fwd/rev), `ncbi_mode`.

**It records the REQUESTED db_sources, not the post-existence-filter set** the Databases
stage computes. Deliberate, and documented in the docstring: the requested set is
deterministic at init and survives a crash; the resolved set is only known mid-stage.
Do not "upgrade" this to the resolved set later — it would move the block out of the one
place that is crash-proof.

**Discipline held.**
- Check seen FAILING on the tree without the init line, for the documented reason:
  `KeyError: 'run_config'` — the real progress.json had no recipe.
- Check PASSING with the fix. `tests/test_run_provenance.py` — anchored on the real
  artifact: it drives `build_stores` (the helper `run_pipeline.py` uses) and the real
  `_run_pipeline_thread`, every stage disabled so no tool is needed, then reads the
  progress.json the pipeline wrote and asserts two runs differing only in `--db-sources`
  are distinguishable. Not a hand-written fixture.
- Full suite `406 passed` (was 403 + 2 new + 1). `golden.py check all` on **fixture** and
  **real**: `no regressions` on all four frozen layers — expected, since progress.json is
  not a frozen file and the change touches no runner or stage output.

**Known item, named not fixed (per instruction).** No `PROVENANCE.txt` is produced by a
real run — only `golden.py freeze` writes one (stage/dataset/source), and `golden_baseline/`'s
was hand-written. Any doc or overview that lists `PROVENANCE.txt` as a normal run output
(the repo overview did) is describing something the pipeline does not emit; it will mislead
a user who goes looking for it. Decide later whether the run should emit its own PROVENANCE
or whether the docs should stop promising one. Not touched now.

## 2026-07-17 — Pre-commit audit: a 2.3 GB near-miss, and the commit order is inverted by the import graph

**Nothing committed. The audit found one thing that had to be fixed before any `git add`,
and three entanglements that make the planned commit sequence impossible as specified.**

### Clean

No credentials (every hit is a reference to an env var, a test placeholder, or docs — `.env`
is ignored and stays ignored). No absolute paths, no username, no personal email anywhere in
the candidate tree. `runs/`, `data/`, `fixtures/`, `golden_baseline/` all correctly ignored.

### The near-miss: `git add .` would have staged 2.3 GB

    databases/bold/BOLD_Public.03-Jul-2026.parquet   2.3 GB   matched by NO pattern
    pr2_version_5.1.1_SSU_mothur.tax                  32 MB   `.tax.gz` was ignored; `.tax` was not
    databases/ncbi/lineage_cache.json                 12 KB   a rebuildable cache

GitHub refuses files over 100 MB, and **no later commit undoes this** — git history is
permanent. Now ignored (`*.tax`, `*.parquet`, `*.snappy`, `databases/bold/`,
`databases/ncbi/`); verified the four `databases/*.py` modules stay visible. **2.3 GB → 872 KB.**

**The root cause is structural and is on the list, not fixed here: `databases/` is both a
Python package and the default DB_ROOT**, so multi-GB downloads land inside the source tree.
Code and data sharing one directory is how a single missed pattern becomes a 2.3 GB mistake.

### The commit order is inverted, and it is the architecture saying so

The planned sequence was machinery → pipeline layers → UI fix → rest of UI. Tested rather
than assumed:

1. **`scripts/run_pipeline.py:121-122` imports `ui.callbacks_progress` and `ui.runtime_state`.**
   The runner that produces every golden **is** the UI. So "machinery first, UI last" cannot
   work: no commit before the UI can execute a run, and therefore **no pipeline-layer commit
   is provable by its own golden until the UI lands.** This dependency is deliberate and
   correct — *"a runner that took its own path could pass while the app stayed broken"* — so
   the ordering must change, not the code.

2. **`scripts/verify.py:455` imports `analysis.results_summary`** (check (g)). It does not
   degrade — proven by hiding the module:

       ModuleNotFoundError: No module named 'analysis.results_summary'

   verify.py **crashes**, so checks (a)-(f) die with it. The instrument is broken until the
   UI-fix commit.

3. **The pipeline layers share files and cannot be split by stage.**
   `ui/callbacks_progress.py` drives all of them in one file (cutadapt/trimmomatic →
   preprocess, vsearch → analysis, taxonomy_runner + build_database_metrics → taxonomy/export,
   phylo_runner → phylo). `databases/database_metrics.py` serves taxonomy
   (`taxonomic_resolution`, `taxon_prevalence`) **and** export (`write_biom`,
   `write_phyloseq_bundle`, `taxa_tables_all_ranks`). One commit per pipeline layer would
   mean splitting files, not selecting them.

**Corrected order, following the import graph rather than the narrative:**

    1. .gitignore hygiene            — tiny, protects every later commit
    2. scripts/golden.py + golden/fixture/  — no project imports; diffs a run dir standalone
    3. core/ + the stage runners     — by directory (they cannot split by stage)
    4. ui/ + app.py + run_pipeline.py — the commit where a golden becomes reproducible
    5. verify.py + results_summary + tests/ — (g) imports cleanly; everything now provable

### Open, and blocking a PUSH rather than a commit

`origin` is **https://github.com/JohnnyFarians24/PARAFISH.git** and its visibility could not
be determined here (`gh` unauthenticated). `golden/fixture/` carries the supervisors' data:
two sample ids (`aptl64b_02`, `aptl64f_01`), per-sample read counts, and the genera *Gadus*,
*Aspergillus*, *Melanogrammus*. That is small and not human data, but **it is theirs, and
publishing is not mine to decide.** Committing is local and reversible; pushing to a public
remote is neither. Resolve visibility + the supervisors' consent before any push.

## 2026-07-17 — Making the contract travel: the fixture IS deterministic, and the repo is the real hole

**The task was "confirm the fixture golden is byte-stable, then commit it so the checks run
on a clean clone". The first half is now proven. The second half turned out to rest on a
premise that is false, and that is the finding.**

### The fixture is byte-stable — proven, not assumed

An **independent re-run** of the fixture pipeline, diffed against the frozen golden:

    preprocess  identical
    analysis    identical
    taxonomy    identical      <- --randseed 42 holds across runs, as designed
    export      identical

The taxonomy being byte-identical on a second run is the SINTAX fix doing exactly what it
was for. Two things had to be separated first, and both were worth the run:

**1. A false alarm that was really a missing recipe.** The first re-run showed four taxonomy
files "changed". They had not: every biological number was identical (OTU1 *Gadus* 100.000,
OTU2 96.667, OTU3 *Aspergillus* 95.918, OTU4 *Melanogrammus* 99.638; `Consensus_Sources`
`NCBI|PR2|SILVA` throughout). The only delta was **BOLD's presence** — the golden was frozen
from a `silva,pr2,ncbi` run and the re-run used `run_pipeline.py`'s default
`silva,pr2,bold,ncbi,custom`. Re-running with the golden's own sources: byte-identical.

**So: `progress.json` does not record `run_config` at all, and `PROVENANCE.txt` records only
stage/dataset/source.** The run does not record how it was made, and the golden freezes an
output whose recipe is nowhere. A config difference is therefore indistinguishable from a
regression, and nobody — including us — can know which invocation reproduces a frozen run.
**Open item, and it is the next one:** the run should write its own `run_config.json`, and
`golden.py freeze` should carry it into `PROVENANCE.txt`.

**2. A real non-determinism: `table.biom` embeds a wall-clock timestamp.**

    frozen : "date": "2026-07-17T08:24:45Z"
    re-run : "date": "2026-07-17T09:51:33Z"

Its checksum could therefore **never** pass a re-run. This is the same defect
`conflict_report.csv` had — *a file that changes without its content changing cannot be
diffed, so a real change hides in the churn* — and it survived because **Phase 4's two-run
proof covered the nine COPIED contract files; `table.biom` is checksummed, not copied, and
was never in that diff.** The net had a hole exactly where it stopped copying.

**Fixed in the checker, not the writer, and the asymmetry is deliberate.** `conflict_report`'s
ordering was meaningless, so the writer was fixed. The BIOM date is *true*, and downstream
tools expect it; making an output lie about when it was written to keep a checker quiet is
the wrong trade for a platform whose argument is that its outputs are honest. So `golden.py`
hashes BIOM **content** (the `date` field dropped, keys canonicalised) and the field stays in
the file. Seen to ignore only what it claims to:

    date changed        -> IGNORED (same digest)
    a count 3960->3961  -> CAUGHT
    a taxon fabricated  -> CAUGHT

### The real hole: the code is not in the repo either

`golden/` being gitignored was the *visible* half. Checking what a clone would actually get:

    UNTRACKED  scripts/golden.py            <- the freeze/check machinery: the whole net
    UNTRACKED  scripts/verify.py            <- all seven checks
    UNTRACKED  analysis/results_summary.py  <- the module this session fixed
    UNTRACKED  tests/test_results_summary.py
    UNTRACKED  databases/database_metrics.py
    UNTRACKED  ui/feature_sets.py           <- the resolver that was right all along
    UNTRACKED  core/run_registry.py
    ... 55 untracked .py files, 49 modified

**The last commit is the initial import.** So committing `golden/fixture/` alone would ship a
contract for code that is not there: a clean clone gets neither the golden, nor `verify.py`,
nor `golden.py`, nor the module under test. "The checks run on a clean clone" was blocked by
something one layer under the golden, and the gitignore was a symptom.

This is the same shape as every other finding in this log — the visible defect was a symptom
of one nobody had looked for — and it is the one that matters most, because **reproducible on
clone was the project's goal.** Ten sessions of making the numbers trustworthy do not survive
a clone that contains none of the instruments that made them so.

`.gitignore` now tracks `golden/fixture/` (176 KB, 26 files, no absolute paths, no dates, no
tool versions, no hostname — verified) and still ignores `golden/real/` (size, and it is real
data). That is necessary and nowhere near sufficient. **What to commit, and in what order, is
a decision left open rather than taken unilaterally:** 55 untracked files spanning work this
log has not audited is not a commit to make while nobody is looking.

State: verify **7/7** real · **7/7** fixture · golden **no regressions** on both datasets *and*
on an independent re-run · **403 pass, 1 skipped** (`duckdb` absent, pre-existing).

## 2026-07-17 — UI session 2: check (g) seen failing, then the layout migration. Phylo held back.

**The gate, executed in order.** Check **(g) the Results page reads this run** was written
first and **seen failing on the frozen real golden run**, before a line of the fix:

    FAIL  (g) the Results page reads this run
          - clustered_read_total -> None, chimera_stats says 154223
          - species_classified_total -> None, taxonomic_resolution says 12
          - dashboard[aptl64b_02].richness -> 0, alpha_diversity says 9
          - dashboard[aptl64f_01].richness -> 0, alpha_diversity says 2
          - dashboard pct_species is 0 for every sample, yet 12 reads reached species

    6 passed · 1 failed

(g) is deliberately unlike (a)-(f): they compare two files, it compares **a file to the code
that reads it** — the gap where a correct pipeline and a blank screen coexist.

**The migration.** `database_metrics_paths(databases_dir, feature_set="otu")` now mirrors
`analysis_metrics_paths` exactly — per-set dir first, flat second — and a new
`taxonomy_table_path` gives non-UI readers the on-disk half of what
`ui.feature_sets.taxonomy_file` does for the UI (`analysis/` must not import `ui/`). The six
`results_summary` paths and the three `callbacks_databases` sites now route through them.
`callbacks_databases.py:265` needed **no change at the call site**: it always passed the
stage dir, and the resolver was the thing that did not know about feature sets.

**Nothing bled backwards** — the risk of "unifying" a resolver five callers already used
correctly:

    golden check all --dataset real     preprocess · analysis · taxonomy · export   no regressions
    golden check all --dataset fixture  preprocess · analysis · taxonomy · export   no regressions
    verify real 7/7 · verify fixture 7/7 · 402 pass, 1 skipped

(The 1 skip is `duckdb` missing in this env, pre-existing and unrelated.)

**The tests were the other half of the bug, and fixing the code did not fix them.** After
the migration `tests/test_results_summary.py` still passed — *because the resolvers fall
back to flat, so the lying fixture kept working*. The test would have passed with or without
the fix: **zero power to detect the defect it was meant to cover.** Fixed by moving the
fixture to the layout the pipeline actually writes, and by adding `GoldenLayoutTests`, which
takes its tree from `golden/` instead of from the test author's typing:

- `test_reader_sees_the_frozen_run` — the reader must return the frozen numbers.
- `test_the_fixture_layout_matches_the_pipeline_layout` — the hand-written fixture is pinned
  against golden, so it cannot drift back into agreeing only with its author.

**Both were seen failing before being accepted**, per the rule:

    revert the migration  -> FAIL: AssertionError: None != 154223
    flatten the fixture   -> FAIL: fixture writes taxonomy_consensus.csv at
                             'databases/taxonomy_consensus.csv', the pipeline writes it at
                             ['databases/otu/taxonomy_consensus.csv']

**A limitation, named: `golden/` is gitignored.** The contract is not in the repo — it is
real run output and lives only on this machine. So `GoldenLayoutTests` **skips** in a fresh
clone or in CI, and a skip there is *not* a pass: it means the contract went unchecked. The
suite that runs everywhere is still the one that can only prove the code agrees with itself.
How the golden contract travels (a committed miniature? checksums? a CI fixture derived from
golden?) is now an open item, and it is the difference between this discipline holding for
the project and holding only on this laptop.

**Phylo: held back, and this is an ordering rule, not a preference.**

> **No page is assessed before its stage has a golden.** The step before reading the phylo
> callbacks is to produce a run that exercises the stage and freeze it.

`callbacks_phylo` has the same flat-vs-per-set shape (`:39`, `:49`) and was **left untouched**.
It looks identical to what was just fixed, and "structurally identical, therefore the same"
is precisely the inference that produced the 401 passing tests — the fixture represented the
disk *because it ought to*. A reader without a contract is an unverified reader, however much
it resembles a verified one. If the pipeline cannot produce a real phylo run on this data,
that is a finding **about the stage**, and better learned now than halfway through assessing
its page.

## 2026-07-17 — PRINCIPLE: a hand-written fixture can only prove the code agrees with itself

**This is a rule, not the write-up of one incident.** It is filed here because it governs
every one of the ~93 callbacks still to be assessed, and because the incident that produced
it is the most dangerous class of defect this codebase has produced so far.

> **A test whose fixture is built by hand by the test's author proves only that the code
> agrees with itself. A contract test anchors on the real artifact the previous stage
> produces.**

The Results page read `databases/taxonomy_consensus.csv`. The pipeline writes
`databases/otu/taxonomy_consensus.csv`. `tests/test_results_summary.py::_run_tree` builds
its fixture at **the path the code expects** — so the test and the code shared one false
belief about the layout and agreed with each other perfectly. **401 tests passed while the
page was blank.**

**Why this is worse than having no tests at all.** No tests is honest: you know that you do
not know. A test whose fixture is invented by the author of the test is a lie with a seal of
approval — the instrument you consult to find out whether something is broken is actively
confirming that the broken thing is intact. The suite was not wrong by accident; it was
structurally incapable of being right, because it never touched anything the pipeline wrote.

**The asymmetry that makes golden files the contract:** a hand-written fixture *can* share
the code's hallucination — the same person, holding the same wrong belief, wrote both. The
pipeline's real output *cannot*, because it was not written by whoever wrote the reader.
That is the entire argument, and it is the tests' version of Rule 2 (*coherence is not
quality*): a fixture and a reader agreeing proves coherence between two things that came
from one head.

**Operationally, for the rest of the UI work:** a check that reads a path or a key must get
it from `golden/`, not from a tree the check itself mkdir'd. Where a fixture is genuinely
needed (speed, edge cases), it must be **derived from a golden artifact**, never typed from
the reader's expectations.

**Third incarnation (2026-07-17), and it proves the rule is not about CSVs — it is about
who built the fixture.** The BOLD positive-parsing probe: a synthetic DuckDB, hand-built to
test `select_bold_top_hits`, used capitalised rank columns because that is what the parser
reads — and it PASSED. The code agreed with itself, in a binary fixture this time. Only
checking the synthetic against the REAL BOLDigger3 DuckDB (`golden_baseline/…`), whose ranks
are lowercase, exposed that the parser blanks every real rank. Same failure as the 401 tests
and `test_results_summary` — author's fixture inherits author's belief — now shown to apply
to a DuckDB as squarely as to a CSV path. The fixture must be derived from, or validated
against, the real artifact; its file format does not exempt it.

## 2026-07-17 — UI Phase, session 1 addendum: the layout drift is systemic, and the store is what hides it

**The question: are the six broken paths specific to Results, or the same flat-vs-feature-set
drift elsewhere?** Swept every stage-directory path construction in `ui/`, `analysis/`,
`databases/`, `phylo/`, `core/`. **It is systemic — one incomplete migration, not a Results
bug.** The pipeline moved to one directory per feature set (`databases/<set>/`,
`phylo/<set>/`, `analysis/metrics/<set>/`); the readers did not all follow.

**Migrated, and correct** — `ui/feature_sets.taxonomy_file`, `analysis_metrics_paths`,
`callbacks_ordination` (:65-67, :91, :146-148), `callbacks_results_summary` (:353-355,
:535-537, :544, :556), `callbacks_analysis` (:62, :103). All resolve per-set first, flat
second.

**Left behind, flat only** — ~11 read sites across three modules:

| module | sites | store primary? | observed |
|---|---|---|---|
| `analysis/results_summary.py` | 6 (:88 :104 :129 :162 :213 :259) | **none** | **blank on golden, proven** |
| `ui/callbacks_databases.py` | 3 (:129 :265 :275) | yes, correct | latent |
| `ui/callbacks_phylo.py` | 2 (:39 :49) | yes, correct | latent, **unproven** |

**The shape of the answer is the useful part: the blast radius is decided by whether a store
masks the fallback.** `callbacks_databases` and `callbacks_phylo` are structurally identical
— read `metrics_files` off the stage record first (which the pipeline writes correctly),
fall back to a flat disk path (which is wrong). `results_summary.py` is **the only consumer
that takes `output_dir` and goes straight to disk with no store to mask it.** That is the
entire reason Results is the blank page and the others are not. Same bug, different exposure.

**So "latent" must not be read as "safe".** The masking store is a `session` store, and
`core/run_registry.py` exists precisely because those stores go empty on restart, in a new
tab, and on arriving at a page without having watched the run go by (DEVLOG 2026-07-14, T3).
**Every masked site becomes a blank site in exactly the scenario the store is known to fail
in.** The mask is made of the thing already documented as evaporating. These are not three
severities; they are one bug and an F5.

**A hole in the net, named rather than assumed.** The phylo claim above is **by construction,
not by evidence**: the writer is `output_dir / "phylo" / key` (`callbacks_progress.py:830`,
`:848`), the reader is `phylo/metrics/`. It could not be tested, because **the golden covers
preprocess/analysis/taxonomy/export and not phylo** — the real golden run has no phylo stage,
and the only runs on disk carrying `phylo/` predate feature sets and therefore have the flat
layout the reader expects. So phylo has **no contract at all**, and the check that would prove
this cannot be written against today's golden. That gap is now the thing to fix before phylo
is assessed, and it is recorded as unproven rather than counted as a finding.

**Consequence for the plan:** the fix is one migration applied to ~11 sites through resolvers
that already exist and are already correct in five places — not four page-by-page repairs. But
the sequence is unchanged and the gate holds: **a golden-anchored check, seen failing on
today's tree for the documented reason, before any path is touched.**

## 2026-07-17 — Why the golden net exists: the unseeded SINTAX, restated as method

The previous entry files "SINTAX was unseeded" under a list of fixes. That is the wrong
shelf, and the misfiling matters more than the fix. **This is not a bug we corrected; it
is the reason the method exists**, and it is recorded here as such so that nobody later
reads the golden discipline as ceremony.

What was actually true before `--randseed 42`: **the classifications were not results, they
were samples from a distribution.** The same reads, against the same reference, on the same
machine, produced different numbers every run. Those bootstraps set `Reliability` and are
compared against `--sintax_cutoff`, so a taxon sitting near the threshold could be **called
on one run and absent on the next, from byte-identical data**. Not a rounding artefact in a
report — a taxon present or absent. On a parasite-detection platform, that is the finding
itself changing its mind.

**Every way we had of looking would have missed it, and this is the whole point:**

- *Reading the code* could not see it. The code was correct in the only sense reading can
  check: it did what it said it did. `--randseed` was not passed wrong, it was not passed.
- *Unit tests* could not see it. There is no assertion to write. Each run's output was
  individually well-formed, internally coherent, and plausible.
- *A single run* could not see it — and this is the trap. Every run looked **perfectly
  good**. The defect had no appearance. It existed only in the relationship *between* two
  runs, which is a place no individual run can look.

It took **two runs over identical input and a byte-level diff** — 0.99 against 0.98, one
cell, one CSV. Nothing cheaper would have worked, because nothing cheaper compares a run to
another run.

**The counterfactual, stated plainly.** Had this project closed on the strength of code
review, passing tests, and runs that looked right — which is what closing normally means —
we would have handed the supervisors a platform whose numbers **nobody could reproduce, us
included**. A reviewer re-running it would get different taxa near the cutoff and would be
right to conclude the platform was unsound. The golden net is what stands between a
scientific instrument and a generator of plausible numbers, and the distinction is invisible
from inside any single run. That is the cost it saved; that is why the discipline is not
optional.

The corollary now governs the UI work: **a check you have never seen fail is not a check,
it is a hope.**

## 2026-07-17 — UI Phase, session 1: the data contract, and the Results page reads a layout that does not exist

**Not "score the page" — map the contract.** The recurring failure mode of this platform is
store/format mismatch: one module writes a shape, another reads a different one, `dict.get`
returns nothing, and the miss is silent. That is exactly what `class_` vs `class` was
(DEVLOG 2026-07-16), and the golden files are now a frozen contract, so the first question
is not "is the dropdown pretty" but **"do the paths and keys the UI reads exist in the
golden?"** They largely do not.

**The finding, by evidence rather than by reading.** `analysis/results_summary.py` run
against the frozen real golden run (`runs/_real/current`), which classified 9/9 OTUs with a
full Linnean backbone:

    clustered_read_total     -> None
    species_classified_total -> None
    sankey stages            -> ['Raw', 'Trimmed', 'Merged']        <- Clustered/Classified absent
    dashboard aptl64b_02     -> richness 0, shannon 0.0, pct_species 0.0
    dashboard aptl64f_01     -> richness 0, shannon 0.0, pct_species 0.0
    target 'Gadus morhua'    -> detected_in '—', prevalence None, reliability '—'

**Every number Phases 3 and 4 fixed is missing from the Results page**, and the cod the
pipeline correctly identified reads as *not detected*. The back end is right and the user
sees nothing — the exact way a working platform is judged broken.

**Root cause: the feature-set directory level.** The pipeline writes one directory per
feature set (`callbacks_progress.py:730` -> `databases/<set>/`,
`analysis_metrics.py:576` -> `analysis/metrics/<set>/`). `results_summary.py` reads the flat
pre-feature-set layout. Six of its eight paths cannot exist:

| pipeline writes (verified in golden) | `results_summary.py` reads | |
|---|---|---|
| `preprocess/metrics/read_tracking.csv` | same | ✓ |
| `analysis/otu_matrix.csv` | same | ✓ |
| `analysis/metrics/otu/chimera_stats.csv` | `analysis/metrics/chimera_stats.csv` | ✗ |
| `analysis/metrics/otu/alpha_diversity.csv` | `analysis/metrics/alpha_diversity.csv` | ✗ |
| `databases/otu/taxonomy_consensus.csv` | `databases/taxonomy_consensus.csv` | ✗ |
| `databases/otu/source_assignments.csv` | `databases/source_assignments.csv` | ✗ |
| `databases/otu/metrics/taxonomic_resolution.csv` | `databases/metrics/taxonomic_resolution.csv` | ✗ |
| `databases/otu/metrics/taxon_prevalence.csv` | `databases/metrics/taxon_prevalence.csv` | ✗ |

The two that work are the two written before feature sets existed. **The bug is not a typo;
it is a module that never learned a layout change**, and it degrades silently because
`_read_rows` returns `[]` for a missing file, which is indistinguishable from an honest
empty result — *the same sin as blastn exiting 0*, one layer up.

**A second site, and its severity is lower — recorded honestly rather than punchily.**
`ui/callbacks_databases.py:265` (`_database_metric_paths`) calls
`database_metrics_paths(Path(output_dir) / "databases")` without the feature set. Against
the golden run that returns `{}`; with `databases/otu` it returns both tables.
(`database_metrics_paths` is the one resolver that is not feature-set-aware at all — its
caller must supply the level, and this caller does not.) **But this is only the fallback.**
The primary path reads `metrics_files` off the stage record, and the pipeline writes that
correctly — verified in the golden `progress.json`:

    "taxon_prevalence_csv": ".../databases/otu/metrics/taxon_prevalence.csv"   <- correct

and `restore_run_state` hands those records back to the stores. So the prevalence card
works whenever the store is populated or restorable, and the broken fallback is reached only
when neither holds. It is a latent bug, not a blank page — unlike `results_summary.py`,
which takes `output_dir` and goes to disk **every time, with no store path to mask it**.
That difference is why Results is blank and Databases is not.

**The correct resolvers already exist and are used elsewhere**, which is what makes this a
drift rather than an oversight: `ui/feature_sets.taxonomy_file` (per-set first, flat
fallback) and `analysis_metrics_paths` (same pattern) both resolve correctly against the
golden, and `callbacks_results_summary.py:353` hand-rolls the same two-step fallback
correctly. Four conventions, three right, one module left behind.

**Why 401 tests pass while the page is blank — the important part.**
`tests/test_results_summary.py::_run_tree` builds its fixture at
`databases/taxonomy_consensus.csv` and `analysis/metrics/alpha_diversity.csv` — **the flat
layout the code expects and the pipeline does not write.** The test and the code share one
false belief, and agree with each other perfectly. The suite was not testing the pipeline's
output; it was testing a directory the test itself invented. Nothing ever compared either to
a real run. **A fixture written from the code under test can only confirm the code under
test** — this is the tests' version of "coherence is not quality" (Rule 2), and it is why the
golden files, not a fixture, are the contract.

**Not fixed in this session, deliberately.** The fix is small (route the six paths through
the existing resolvers, pass the feature set at `callbacks_databases.py:265`) and the
tempting move is to make it now. But the check that proves it must first be seen **failing**
on the current tree — otherwise it inherits the exact defect just documented: a fixture
agreeing with the code. Sequence: a golden-backed check that fails on today's `results_summary`
for the reason above, then the fix, then the check passes. Only then is the dropdown worth
looking at.

## 2026-07-16 — Phase 4 (export) CLOSED without code changes. Plan 0-4 complete.

**The question was whether `taxa_table_class`'s column offset was still a bug or a symptom
of the taxonomy Phase 3 fixed upstream. It was a symptom — but there were two distinct
problems wearing the same face, and separating them mattered.**

1. **A real, independent bug**: `database_metrics.TAXONOMY_RANKS` spells the rank `class_`,
   `community_analysis` keyed on `class`, the `dict.get` returned `{}` in silence and every
   read fell into `Unclassified`. **Already fixed before this plan started**
   (`canonical_rank` / `_RANK_ALIASES`), with tests pinning it.
2. **The symptom the plan describes**: `taxa_table_class` reading `Unclassified, 100%`. Not
   an offset — `tax_table.Class` was *genuinely empty*, because no Linnean source reached
   the consensus. Phase 3 fixed it upstream.

**The decisive evidence is check (c) passing on `golden_baseline`** — the pre-Phase-3 run
with the broken taxonomy. `taxa_table_class` already agreed with `tax_table` there; both
said empty. Had the offset still existed, the two would have disagreed. Diagnosis by
evidence, not by reading the code.

**Audited the instrument, not just the output.** Check (c) only ever tested *class* — so
all seven ranks were verified by hand against `tax_table`: kingdom, phylum, class, order,
family, genus, species all agree, all totalling 154 223. The rest of Phase 4's targets hold
too: `phyloseq/otu_table.tsv` matches the canonical abundance (9 OTUs, 154 223);
`table.biom` carries 154 223 and the harmonised lineage with `s__` empty where the LCA found
no agreement — NA, not invented.

**Nothing to fix, so nothing was touched.** Frozen: `golden/{real,fixture}/export`.

**Phases 0, 1, 2, 3 and 4 are closed.** The two that closed without code changes (2 and 4)
did so because the work had already been done upstream — which is what the phase ordering
was for.

## 2026-07-16 — Two permanent rules, the megablast trade-off, and a non-determinism

These are **rules, not fixes for one run**. Both were learned from blastn; neither is
about blastn.

### Rule 1 — judge an external tool by its OUTPUT, not by its exit code

`blastn` exits **0** having printed `(Severe Error) ... CPU usage limit was exceeded` for
every query and written an empty result file. `check=True` therefore saw success, and the
refusal travelled to `source_status` as `no_hits` — the one thing that column exists to
prevent. An exit code is a claim the tool makes about itself; the output is evidence.

The other half, which a naive version of this rule gets wrong in the opposite direction:
**an empty output is not automatically a failure.** A filter that legitimately passes
nothing and a crash produce the same empty file. Only the tool's own error output tells
them apart:

    empty output + error signature  -> FAILURE (refused / died)
    empty output + clean run        -> a real, honest, empty result

"Empty means broken" would turn every honest absence into a false alarm — the original sin
wearing the opposite mask. `core/process.py` encodes both halves (`run_tool`,
`find_failure_signature`, `ToolResult.ok`) so the rule lives in one place.

Applied beyond NCBI to `databases/custom_db.py`, which is the **exact twin**: same blastn,
same `check=True`, same blindness — and a 0-byte `custom_blast_raw.tsv` was already on the
fix plan's list. The remaining `check=True` call sites (vsearch, cutadapt, PEAR,
Trimmomatic, MAFFT, FastTree, RAxML, FastQC) are audited but **not** blanket-converted:
that would be fixing by assumption, which this plan forbids. Each needs its own evidence
of how it fails, and `run_tool` is what they move to as that evidence arrives. Two already
handle the signal case (`taxonomy_runner`, `dada2_runner` check for negative return codes
= OOM killer), which is the same rule discovered independently.

### Rule 2 — coherence is not quality

The five checks test whether the files agree with **each other**. They all passed on a run
whose only Linnean source had died: nothing contradicted anything, every rank was honestly
empty, `5/5`. **A coherent run can be a useless run.** Check (f) asks the question the
other five structurally cannot — *did anything actually deliver the goods?* — and any
future check suite needs at least one of its kind. Coherence checks can only ever say "no
lie was told"; they cannot say "something true was said".

### The megablast trade-off — the concrete number

Remote (default): **`task=megablast`**, which seeds on a **28 bp exact match**. Our
`perc_identity` floor is **90%**, where mismatches fall roughly **every 10 bp** — so a 28 bp
exact seed is unlikely to exist in a genuinely divergent hit, and megablast **can miss a
parasite at ~90% identity that blastn would find**. This is the sensitivity we are paying
for reliability, stated so nobody has to rediscover it.

It is paid because the alternatives do not run at all — measured, same 9 OTUs, same day:

    -task blastn -word_size 11   -> REFUSED (SIGXCPU on all 9)
    -task dc-megablast           -> REFUSED (SIGXCPU on all 9)
    -task megablast              -> 123 hits, 9/9 OTUs, ~3 min

**Where the sensitivity lives instead:** `PARAFISH_NCBI_LOCAL_DB=1` -> `task=blastn,
word_size=11`, on our own CPU, where no public quota can refuse it. That is the only place
the ~90% divergent hit is findable. **Do not "optimise" the remote default to blastn or
dc-megablast — it does not fail slowly, it fails completely.**

### SINTAX was unseeded — the taxonomy was not reproducible

The most serious of the three, and the golden net is the only reason it was ever seen.
Two fixture runs over byte-identical input, everything else matching:

    OTU4,SILVA,species,Gadus_morhua_(Atlantic_cod),,,0.99,high     <- run 1
    OTU4,SILVA,species,Gadus_morhua_(Atlantic_cod),,,0.98,high     <- run 2

SINTAX bootstraps by randomly subsampling k-mers, and vsearch's `--randseed` **defaults to
0, meaning "seed from a random data source"**. We never passed it. So the pipeline's
classifications were samples from a distribution, not results: **the same reads against the
same reference gave different numbers every run.**

This is not cosmetic. Those bootstraps set `Reliability` and are compared against
`--sintax_cutoff` — so a taxon sitting near the threshold can be **called on one run and
absent on the next, from identical data**. A reviewer re-running the pipeline could not
reproduce our numbers, and neither could we. Now `--randseed 42` (`SINTAX_RANDSEED`); any
fixed non-zero value would do, the point is that it is not 0.

Worth noting what this says about the method: a 0.99-vs-0.98 flicker in one cell of one CSV
is invisible to any amount of reading the code, and no test would have caught it — the code
was "correct" in the sense that it did what it said. It took two runs of identical data and
a byte-level diff.

### `conflict_report.csv` was non-deterministic

Caught by the golden discipline doing its job: two fixture runs over identical data, every
other taxonomy file byte-identical, and `conflict_report.csv` with the same rows shuffled.
`_build_conflict_report` walked `records_by_asv` in insertion order — i.e. the order the
sources happened to answer in — while the consensus beside it used
`sorted(records_by_asv.items())`. Now both sort. A file that changes without its content
changing cannot be diffed, which means nobody would notice it changing for a reason that
*did* matter.

### For the list (not now): the Linnean dependency is invisible

Removing SILVA/PR2 from `LINNEAN_SOURCES` was right — forcing their non-Linnean backbones
into Phylum/Class is what invented `Phylum=Opisthokonta`. The consequence to own: **the
Linnean backbone now rests on NCBI alone.** Without NCBI reachable there is no
Phylum..Family at all, and today a user discovers that from a failed run rather than being
told up front. Check (f) fails loudly, which is the floor, not the fix — the report/UI
should say it before the run, not after.

## 2026-07-16 — The NCBI refusal was our request, not NCBI's mood

**I blamed the server without auditing the request. The audit says the request was
wrong.** Same 9 OTUs, same query, same day, against remote `nt`:

    -task blastn -word_size 11   -> REFUSED, SIGXCPU on all 9      <- what we were sending
    -task dc-megablast           -> REFUSED, SIGXCPU on all 9
    -task megablast (default)    -> 123 hits, 9/9 OTUs, ~3 min

`word_size 11` against all of nt is a full blastn search — the most expensive thing you
can ask of a shared public quota, which is precisely what that quota exists to refuse.
Nobody chose this: `NCBIConfig.word_size` defaulted to 11 while `task` was never set, so
the search had blastn's cost under megablast's name. Remote now asks for `megablast` and
lets the task pick its own seed length.

**The rate-limit hypothesis was wrong too, and worth recording.** The requests were
already batched: `-query` is the whole multi-FASTA and blastn submits it as ONE job. The
nine identical errors are NCBI rejecting each query *inside* one job, not nine requests.
One run makes one request; retries are now the only way we could approach NCBI's
"once per 10 s", hence spaced growing backoff (30 s, 120 s, give up) rather than an
immediate retry — and giving up returns a **failure with the reason**, never a silent zero.

**blastn exits 0 when NCBI refuses it.** This is the load-bearing finding. It prints
`(Severe Error) ... SIGXCPU` for every query, writes an empty result file, and returns
**exit code 0** — so `check=True` never fires, `run_blast` returned `success=True,
hit_count=0`, and `taxonomy_runner` faithfully translated that into `no_hits`. The
refusal entered the system wearing the word for an honest absence. `run_blast` now reads
its own stderr for refusal signatures and fails loudly. Reconciliation's 0-byte rule was
catching this only by accident, downstream.

**Sensitivity, honestly.** megablast seeds on a 28 bp exact match; at ~90 % identity
(our `perc_identity` floor) mismatches fall every ~10 bp, so it *can* miss a divergent
parasite that blastn would find. dc-megablast is the textbook answer for that and the
public quota refuses it. So: remote (the default, per the supervisors' decision) uses
megablast because it is the only search that completes, and the **local db is where the
sensitive search is affordable** — `PARAFISH_NCBI_LOCAL_DB=1` opts in and gets
`task=blastn, word_size=11`, on our own CPU.

**Local BLAST is now opt-in and off by default.** It used to be selected automatically
whenever the files happened to be on disk, which silently swapped the source out from
under a run. The public remote nt is the source a reviewer can reproduce without our
filesystem; the local db is the emergency exit for when NCBI is down.

**Sixth check, mandatory: (f) a Linnean backbone reached the consensus.** The other five
test whether the files agree with each other, and they all passed on the run where NCBI
was refused — nothing contradicted anything, every rank was honestly empty, 5/5. That is
the one way this suite could mislead, so the question is now asked directly. On the
refused run:

    FAIL (f)  no Linnean source fed any feature's consensus (saw: PR2|SILVA) · NCBI=failed
              — Phylum..Family cannot be filled by SILVA/PR2, whose backbones are not Linnean

SILVA and PR2 are excluded from `LINNEAN_SOURCES` on purpose: their backbones are not
Linnean (PR2 has Supergroup/Division/Subdivision), and forcing them into those columns is
what invented `Phylum=Opisthokonta` in the first place.

## 2026-07-16 — Phase 3 (taxonomy) CLOSED and frozen

**NCBI cooperated with a well-behaved batched request.** `NCBI,ok` on the real dataset,
141 867 bytes of hits, `blast_no_hit.txt` **empty** — every OTU matched. The refusal was
never NCBI's mood; it was our `word_size 11` (see the entry above). Real dataset, 6/6:

    (a) read totals agree            all 11 sources = 154 223
    (b) otu_id_map.size              9/9 OTUs match
    (c) taxa_table_class             5 classes match
    (d) no invented middle ranks     no Opisthokonta/Amorphea/Obazoa/Craniata
    (e) source_status                NCBI=9 · PR2=9 · SILVA=9
    (f) Linnean backbone             NCBI · 9/9 features with Phylum..Family

The table this was all for:

    OTU1  Eukaryota  Chordata    Actinopteri      Gadiformes         Gadidae             Gadus
    OTU4  Eukaryota  Ascomycota  Eurotiomycetes   Eurotiales         Aspergillaceae      Penicillium coffeae
    OTU5  Eukaryota  Ascomycota  Dothideomycetes  Pleosporales       Massarinaceae       Pseudosplanchnonema phorcioides
    OTU6  Eukaryota  Chordata    Actinopteri      Gadiformes         Gadidae             Melanogrammus
    OTU7  Eukaryota  Arthropoda  Insecta          Lepidoptera        Saturniidae         Saturnia pavonia
    OTU8  Eukaryota  Chordata    Mammalia         Rodentia           Gliridae            Glis glis
    OTU9  Eukaryota  Ascomycota  Dothideomycetes  Mycosphaerellales   Mycosphaerellaceae  Nothopassalora personata

Every rank from Phylum to Family was empty before this; `taxa_table_class` said
`Unclassified, 100%`. It now resolves five real classes. `Consensus_Sources` is
`[NCBI|PR2|SILVA]` on all nine features — PR2 was harmonised, not replaced, and the table
says so rather than asking to be believed. Species stays empty where the LCA found no
agreement (OTU1/OTU2: *Gadus* genus, no species) — NA, not invented.

**Frozen:** `golden/{real,fixture}/taxonomy`. Golden rule held throughout — `preprocess`
and `analysis` PASS on both datasets; Phase 3 moved no number upstream. Tests: **399 pass**.

**Phases 0, 1, 2 and 3 are closed. Phase 4 (export) is next**, and its inputs are now a
frozen contract.

## 2026-07-16 — Phase 3 (taxonomy): the fixes, before NCBI would run on the real data

**The fabricated attribution — one line, and it explains the wrong species.**
`best_alignment_record` skipped records whose identity was `None`. BOLD writes its
non-identifications with `similarity = 0.0` — which is not `None`, so they scored; and
since SILVA/PR2 are SINTAX and carry no identity at all, BOLD's non-answer was the *only*
scoring record. Hence `Best_Source=BOLD, Best_Identity=0.0` on all 7 features while BOLD
had identified nothing and the species on the row had come from SILVA. **The consensus was
citing, as its evidence, the one source that had no evidence.** A record that named no
taxon, or aligned at 0, is not a hit.

**`Consensus_Sources`, a new column, because "PR2 was harmonised, not dropped" should be a
fact and not a promise.** `Best_Source` names the strongest *alignment*, so a SINTAX source
can never hold it once a BLAST source is present — SILVA and PR2 decided calls while being
invisible in the output. The column lists the sources that actually fed each feature's LCA.
On the fixture: `[NCBI|PR2|SILVA]` for every feature, `Best_Source=NCBI` at 100 %, 96.7 %,
95.9 %, 99.6 %.

This also killed a **false positive in my own check**: verify (e) had been asserting that a
source with assignments must reach the consensus, measured via `Best_Source` — so it
flagged "SILVA never reached the consensus", which is normal, not a bug. Participation is
now read from `Consensus_Sources`; attribution from `Best_Source`. Two questions, two
columns.

**Status honesty, second pass.** `_reconcile_source_status` only ever judged sources
claiming `ok`. NCBI labels itself `ok if hits > 0 else no_hits` up front, so a BLAST killed
server-side arrived already saying `no_hits` and was waved through — **a failure wearing
the words of an honest absence**, which is the exact sin this column was rewritten to stop.
It now judges every source that ran (only `skipped` is exempt), against assignments, marker
compatibility, and raw output. The real run now reads `NCBI,failed` off its 0-byte result
instead of `no_hits`. Also added: `not_applicable` for a marker the source cannot carry
(BOLD indexes COI; on 18S its 0/7 is not an absence, it is a question it could not answer).

**Linnean ranks, on the fixture, with NCBI alive:**

    OTU1  Eukaryota  Chordata    Actinopteri     Gadiformes  Gadidae         Gadus
    OTU3  Eukaryota  Ascomycota  Eurotiomycetes  Eurotiales  Aspergillaceae  Aspergillus
    OTU4  Eukaryota  Chordata    Actinopteri     Gadiformes  Gadidae         Melanogrammus

Phylum..Genus were empty for every feature before this. Species stays empty where the LCA
found no agreement — NA, not invented. No Opisthokonta, no Craniata.

**A limitation worth naming: verify.py passes 5/5 on the real run where NCBI died.** The
checks test *coherence*, not *quality*. A run with no Linnean source is perfectly coherent
and perfectly useless for Phylum..Family. The checks cannot be the only thing anyone looks
at, and `source_status` saying `failed` is what carries that signal.

**Open, and it is a platform question, not a bug:** the Linnean backbone currently depends
on a flaky remote service that rate-limits us. SILVA and PR2 were made to work locally by
cutting the reference to the marker + primers (788 MB -> 25 MB); the same trick is the
obvious candidate for a local BLAST subset, and would make Phase 3 reproducible instead of
dependent on NCBI's mood.

Tests: **396 pass**. Golden rule holds: `preprocess` and `analysis` PASS on both datasets —
Phase 3 moved no number upstream.

## 2026-07-16 — Singleton policy: retain, and it dissolves the second accounting

**Decision: `min_unique_size` defaults to 1 — singletons may seed OTUs.** Configurable via
`PARAFISH_MIN_UNIQUE_SIZE`; `core.config.MIN_UNIQUE_SIZE` is the one place it is decided.

**Why, precisely.** "Singletons are sequencing error" is folklore, and it is wrong about
the reads a floor of 2 actually discards. An error-derived read sits within 3% of the
sequence it is an error of, so 97% clustering absorbs it into that OTU whatever the floor
— it never seeds anything. What a floor of 2 *uniquely* discards is the read that is
globally unique AND >3% from everything else. That is not the error profile; it is the
rare-novel-taxon profile. On a parasite-detection platform that is the most expensive
false negative available, and it is unrecoverable after the fact. Time is the cheap
variable; a missed parasite is not repurchasable.

**This is not hypothetical — the fixture had one.** `aptl64b_02`, floor of 2:

    filtered 3969 -> mapped 3968, unmapped 1     (3 OTUs)
    filtered 3969 -> mapped 3969, unmapped 0     (4 OTUs)   <- floor of 1

That one read was globally unique and >3% divergent, so the floor barred it from seeding
and it then matched no OTU: it left the pipeline without appearing anywhere. With the
floor at 1 it founds OTU4 and gets counted and classified. The real dataset happened to
have zero such reads — the cost was luck, not design.

**The unexpected dividend: the two accountings collapse into one.**

    floor 2:  sum(size)=154 223   sum(cluster_size)=151 466   gap = 2 757
    floor 1:  sum(size)=154 223   sum(cluster_size)=154 223   gap = 0

With nothing excluded from clustering, the clustering count *is* the mapping count. The
147 010-vs-157 527 class of bug — the reason this whole plan exists — stops being
something we document and becomes structurally impossible. `cluster_size` stays as a
column: at floor 1 it is a tautology, at floor 2 it is the diagnostic.

**Cost, honestly:** clustering time, 57 s -> 247 s on the real run (4.3x), because 4 232
uniques are clustered instead of 1 475. A large host-dominated study may legitimately buy
that back with a floor of 2 — which is why this is a config parameter and not a code edit.
`read_tracking.csv` now carries `unmapped_reads`, so whatever the floor costs shows up as
a number instead of vanishing. That visibility is what makes the choice safe to expose.

**One default, not three.** It was hardcoded in `core/models.py` (2), re-hardcoded as a
fallback in `vsearch_runner` (2), and defaulted to 1 in the Setup page — so the "policy"
depended on which entry point you came in through. All three now read `MIN_UNIQUE_SIZE`.

**Correcting my own claim from the previous turn.** I reported the two floors as producing
"identical output". They do not, and I had compared aggregates (OTU count, read total)
rather than file contents. What is actually true: the *set of 9 OTU sequences* is
identical and the *multiset of abundances* is identical
(`[2,2,2,2,4,4,12,43,154152]`, total 154 223). What changed is `cluster_size` (converged
to `size`) and the OTU **ids**, permuted among equal-abundance ties because vsearch orders
by cluster size. Biology unchanged; labels and one column moved.

**And a method error.** I froze `golden/*/analysis` at the end of the previous session
*before* Phase 2's second deliverable — this policy — had been decided, then read the
resulting diff as a regression. It was not: it was the approved change landing on a
premature freeze. Analysis is re-frozen under the decided policy. The signal that
mattered was `preprocess`, which passed untouched. **A stage is not frozen until every one
of its deliverables is closed, not merely its checks.**

Tests: **392 pass**. `golden.py check all`: preprocess + analysis PASS on both datasets.

## 2026-07-16 — Phases 0 (safety net), 1 (preprocess) and 2 (analysis) closed

**The plan this follows was written against an older tree, and most of it was already
done.** Rather than take its word for it, Phase 0 built the net that answers the question
empirically. Result: on `golden_baseline` (the pre-fix reference run) four of the five
checks already passed. Phase 2's close criteria — `otu_id_map.size` == its otu_table row
sum, `chimera_stats.total_reads` == the otu_table total — hold at 157 527. Phase 4's
"column offset in the class collapse" was already fixed (`canonical_rank`). What the plan
called Phase 1's root cause ("integrate cutadapt before dereplication") was not the bug
either; see the previous entry.

**Phase 0 — four deliverables, all confirmed with real output.**

- `scripts/make_fixture.py` — 5 000 read-pairs/sample. seqtk is not installed; this
  samples with a fixed seed and *one* index set applied to both mates. Sampling the mates
  independently yields a fixture that silently fails to assemble; verified R1/R2 read ids
  are identical and the gzip output is byte-stable across rebuilds.
- `golden_baseline/` — the 120-file pre-fix run, kept out of `runs/`.
- `scripts/run_pipeline.py` — the pipeline had **no entry point outside the Dash UI**, so
  no run was reproducible. This drives the same `run_config_from_store` +
  `_run_pipeline_thread` the app drives, synchronously, to a fixed run id. A runner with
  its own stage order could pass while the app stayed broken. Fixture: **20 s**. Real: **57 s**.
- `scripts/verify.py` — the five checks, PASS/FAIL/SKIP (a SKIP is not a pass: a run
  without the Databases stage must not look like a clean bill of health).

**The net had a hole, and closing it found the biggest remaining bug.** Check (e)
originally walked `source_status`'s own rows, so a source that ran and was then dropped
entirely was *invisible* — it passed 5/5 on a run that names the wrong species. Judging
each source against its raw output on disk instead:

    FAIL (e)  NCBI: ran (blast_result.tsv, 9083 bytes) but is absent from source_status
              CUSTOM: ran (custom_blast_raw.tsv, 0 bytes) but is absent from source_status

`blast_result.tsv` holds 68 hits: OTU1 and OTU2 both *Micromesistius poutassou* at
99.296 % and 99.766 %. NCBI ran, produced exactly the right answer, and nothing consumed
it — while the consensus reports *Gadus morhua* from SILVA. The previous entry assumed
NCBI "is not being exercised here"; it was. **This is why phylum..family are empty** —
NCBI is the only Linnean source. That is Phase 3, and it is now the top of the list.

**Phase 1 — the primer fix, confirmed by a run rather than by reading the code.**

| | baseline | after |
|---|---|---|
| OTUs carrying TAReuk454FWD1 | **7/7** | **0/9** |
| OTU1 / OTU2 | 80 800/13 405 and 1 244/62 016 — one fish split in two | **collapsed**: OTU1 = 79 761 + 74 391 |
| `reads_with_primer` | empty | 99.98 % / 100 % |

The columns were never "unfilled": cutadapt passed through, so there was nothing to fill.
Same root cause as the split OTU.

**The ~22 k reads with no trace were the quality filter.** Not derep, not chimeras. It
runs in the *Analysis* stage, and `build_preprocess_metrics` only ever walked
`preprocess/`, so the waterfall stopped at "assembled" and the loss was never mentioned:

    assembled 179 605 − filtered 157 537 = 22 068   (+10 unmapped = 22 078)

`read_tracking.csv` now spans the whole pipeline. `build_preprocess_metrics` takes an
optional `analysis_dir` and the runner calls it a **second time** after Analysis — one
writer, two calls, rather than a second table that could drift. New columns:
`primer_trimmed_reads`, `filtered_reads`, `filter_discarded_reads`, `mapped_reads`,
`unmapped_reads`, `pct_filtered`, `pct_mapped`. `mapped_reads` is read from the otu_table
(the canonical basis) rather than recomputed. `failing_steps` now watches filtering — the
step that was eating the reads was the one nobody was watching.

**Deliberate deviation from the plan: no per-sample `dereplicated` / `non_chimeric`
columns.** The plan asks for the chain `... -> dereplicated -> non_chimeric ->
mapped_to_OTU`. Those are *pooled* operations, and `--minuniquesize 2` excludes reads
from clustering while `--usearch_global` counts them anyway. A per-sample "dereplicated"
column would therefore show a loss that does not exist — it would re-introduce the exact
two-accountings confusion the previous entry removed. That accounting is pooled and stays
in the vsearch pool log.

**Real dataset, every read accounted for:**

    input 191 410 -> trimmed 162 456 -> assembled 162 446 -> filtered 154 223 -> mapped 154 223 (+0 unmapped)
    per sample: filtered == mapped + unmapped   ✓ both

(`trimmed` drops further than the baseline's 188 355 because MINLEN=200 now applies to
reads that are ~20 bp shorter — the primer is gone. Legitimate, and the reason downstream
totals move.)

**Also deviating on `source_status`:** the plan wants BOLD's 0/7 called `failed`. It is
`no_hits`, and that is more honest — BOLD is COI-only and this is 18S data, so assigning
nothing is a true result. `verify.py` encodes the implemented rule (empty/absent raw
output ⇒ `failed`; real output ⇒ `no_hits`), not the plan's wording.

**Phase 2 closed without code changes**, which is the honest outcome: its criteria are
met on the regenerated real run — `otu_id_map.size` equals its otu_table row sum for
9/9 OTUs, and otu_table = otu_matrix = `chimera_stats.total_reads` = 154 223. The
`cluster_size` column stays: it is not the "second accounting" the plan wants deleted but
the diagnostic that makes the mapping/clustering gap visible, and the previous entry is
why. Nothing to fix, so nothing was touched.

**Frozen.** `scripts/golden.py freeze|check` — the contract files verbatim plus md5s for
the bulk FASTQ/FASTA, namespaced by dataset (fixture and real legitimately differ).
`preprocess` and `analysis` are closed for both datasets. Verified the check actually
fails on changed input rather than rubber-stamping it, which matters now: **Phase 3 must
not move a single number in either.**

Tests: **390 pass**. verify.py on the real run: (a) and (b) PASS at 154 223, (c)(d)(e)
SKIP (no Databases stage). **Phase 3 is next, and NCBI is the whole story.**

## 2026-07-16 — Singleton policy, decided and made visible (review #7, plan #5)

**The decision.**

> Singletons do not **found** OTUs, but they do **count** toward the abundance of OTUs
> founded by others.

`--minuniquesize 2` stays. This is not the contradiction it looked like. The review read
"derep has no size=1, so singletons were eliminated" and concluded the run had thrown
them away — but dereplication only decides what may *seed* a cluster. The mapping step
runs `--usearch_global` over the **full** pooled file, so those reads are still counted
against whichever OTU they match. On run 20260716T080902Z that is 10 527 reads: excluded
from OTU formation, present in `otu_table`.

**Why it serves the stated goal** (tell true absence from dropout): a taxon abundant in
one sample and seen once in another *is counted* in the second — 1, not 0. That is the
dropout case, and it works. What is deliberately given up is a taxon seen exactly once in
the entire run: it seeds no OTU and is lost. A globally-unique read is overwhelmingly
sequencing error, and letting it found an OTU is what produced the spurious singleton
OTUs this floor was added for.

**Coherence between the two halves is the thing to protect.** Raising `min_unique_size`
without knowing that the mapping still counts everything is precisely what let the
pipeline report 147 010 reads in some files and 157 527 in others. The policy is now an
output rather than folklore: `run_vsearch_pool_and_cluster` returns
`min_unique_size / pooled_reads / uniques_written / unique_seqs_discarded /
mapped_reads / unmapped_reads` and writes the policy into its log.

**Verified on the run — the accounting closes end to end:**

    pooled_reads           157 537
    unique_seqs_discarded   10 527   (kept out of clustering, still counted)
    mapped_reads           157 527   == otu_table total
    unmapped_reads              10   (matched no OTU at id=0.97)
    157 527 + 10 = 157 537 ✓

The 10 unmapped reads are the only ones that leave the pipeline uncounted, and they are
now named rather than missing. Tests: 390 pass.

**Note for the next run:** with primers removed (previous entry), the singleton picture
changes — part of the apparent uniqueness was the heterogeneous 5' spacer, not biology.
`unique_seqs_discarded` should fall.

## 2026-07-16 — Primers actually removed (review #3, #1)

**The diagnosis changed on inspection.** "cutadapt is outside the pipeline" is not what
was wrong. `run_cutadapt` *is* called (`callbacks_progress.py:221`), on the raw reads —
well before dereplication, so the ordering requirement was already met — and it uses
non-anchored `-g`/`-G`, which remove the primer **and everything 5' of it**, i.e. exactly
the heterogeneous spacer. The code was right.

**The actual root cause is one line of UI.** The preset dropdown shipped with
`value=""` — and `""` is the value of the option *"None — do not strip primers"*. So the
pipeline's default was *keep the primers*, `resolve_primers` returned `("", "")`,
`run_cutadapt` took its `_passthrough` branch, and the run recorded that as a success
with `primers_removed: 0.0`. An untouched control was indistinguishable from a deliberate
opt-out. This run: `primer_preset=None, primer_fwd=None, primer_rev=None`.

**Fix.** A new `"auto"` preset — now the dropdown's default — resolves the marker's
standard pair (`MARKER_DEFAULT_PRESET`, 18S → TAReuk). `""` stays a real opt-out and is
obeyed; an absent key (runs stored before "auto") is treated as auto. COI has no auto
default on purpose: Folmer and Leray are both standard and amplify different regions.
Guessing is safe *here* and nowhere else — cutadapt matches the primer as a sequence and
removes it only where it occurs, unlike the Analysis stage's blind length-strip.
`callbacks_setup` resolves "auto" the same way, so the preflight advises on the primers
the run will actually use.

**The prediction was falsifiable and it held.** Applying cutadapt to
`pooled_filtered.fasta` (the direct input to dereplication) and re-running
derep→cluster→chimera→map with the pipeline's own parameters:

    reads carrying the primer   157 537 / 157 537  (100%) — it was never removed
    OTU1 94 205 + OTU2 63 260 (split cod)  →  a single OTU with 157 439 reads
    run total                   157 527 → 157 527   (invariant (a) holds)

OTU1 and OTU2 were one *Micromesistius poutassou* split only by the 5' spacer
(GCATAGCA vs TCCTGAAC). They now collapse. The OTU *count* stays 7 because the other six
are the low-abundance contaminants, not because the merge failed.

**Not verified, and why:** `reads_with_primer`/`pct_primer` are filled from cutadapt's
`--json` report by `preprocess_metrics`, so they populate once cutadapt runs — but this
run's raw FASTQs are gone (`KEEP_INTERMEDIATES=0`; only the Trimmomatic-trimmed reads
remain), so preprocess cannot be re-run end-to-end here to show them filled. The
resolution logic is pinned by unit tests instead, including one that guards the actual
regression (the dropdown's default value). A fresh run is what will show those columns.

Tests: 390 pass.

## 2026-07-16 — Status = real result (review #2 rule 3, #5)

**Provenance.** Every source decided its own status from a *different* proxy, and none
of them was "did this source assign a taxon to anything" — the only question the column
claims to answer:

| source | old rule | proxy |
|---|---|---|
| BOLD | `"ok"` unconditionally | its parquet/xlsx **file existed** |
| NCBI | `ok if hits > 0` | BLAST's `hit_count` — rows, not usable assignments |
| custom BLAST | `ok if matched>0 else no_hits` | a match counter; 0-byte raw output still read as `no_hits` |
| SILVA/PR2 | `_classify_sintax_source` | counted real assignments — the only correct one |

That is how BOLD reported `ok, identificação BOLD concluída` on a run where all 7 of its
rows were `unassigned, identity 0.0`.

**Fix.** `records_by_asv` is the single source of truth — it is what
`source_assignments.csv` and the consensus are built from. New
`_source_assignment_counts` / `_reconcile_source_status` give it the last word: each
`source_status` entry is tagged with its record key and its raw output, and before the
CSV is written any source claiming `ok` with zero assignments is downgraded. Empty or
absent raw output ⇒ `failed` (a broken database, not an honest absence); real output ⇒
`no_hits`. `skipped` is left alone, so SILVA's out-of-memory message is not relabelled.

**Verified on run 20260716T080902Z** — BOLD went `ok` → **`no_hits`**, "correu mas não
atribuiu taxonomia a nenhuma das 7 features". That is correct rather than merely
different: BOLD is COI-only and this is 18S data. SILVA/PR2 now read `ok, 7/7 features
com taxonomia` and, thanks to the marker-filtered references, run locally (SILVA 788 MB
→ 25 MB).

**All five checks now hold at once on the regenerated run:**

- (a) otu_table = otu_matrix = chimera_stats = taxonomic_resolution = every taxa_table = **157 527**
- (b) every `otu_id_map.size` equals its otu_table row sum
- (c) `taxa_table_class` and `tax_table` agree
- (d) no invented middle ranks — Opisthokonta, Craniata, Amorphea, Obazoa all gone
- (e) `source_status` matches `source_assignments`: BOLD no_hits/0, SILVA ok/7, PR2 ok/7

Tests: 383 pass.

**Still open (untouched, later iterations):** cutadapt before dereplication (#3 — OTU1/
OTU2 are still one organism split by the 5' spacer), full read-tracking (#4), the
singleton policy (#5 — `--minuniquesize 2` still drops 10 527 reads from clustering
while the mapping counts them), and NCBI, which is the only Linnean source and is not
being exercised here, so phylum..family stay empty.

## 2026-07-16 — Harmonised taxonomic backbone (review #2 rule 2, #3, #4)

**Two bugs, one root: modules disagreeing about what a rank is.**

*(c) One rank, two spellings.* `taxa_tables_all_ranks` iterates
`database_metrics.TAXONOMY_RANKS`, which spells the rank `class_` (dodging the Python
keyword), and passes that string into `community_analysis._load_taxonomy`, whose
`_RANK_COLUMNS` is keyed `class`. A plain `dict.get` — so it never raised. It returned
`{}`, and `taxa_table_class.csv` reported 100% "Unclassified" for a run whose
`tax_table.tsv` said Craniata. Only `class` was affected: it is the one rank whose name
carries an underscore. Fixed with an explicit `_RANK_ALIASES`/`canonical_rank`, not by
duplicating the dict key, so `RANKS` stays the canonical seven.

*(d) Two non-Linnean backbones, both forced onto Linnean ranks.*

- PR2 v5 (`Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species`) was
  mapped `("d", None, "p", None, "c", "o", "f", "g", "s")` — Division→phylum,
  PR2's Class→class. That is where Phylum=Opisthokonta and Class=Craniata came from.
- SILVA's `tax_silva` **declares no ranks at all**, and `silva_rank_labels` claimed
  `p, c, o, f` positionally plus `g` on the last-but-one.

Fixing PR2 alone just handed the middle ranks to SILVA: the run came back with
Phylum=Amorphea (a supergroup), Class=Obazoa (a clade), Genus=Teleostei (an infraclass).
Same fabrication, other source. Both are now cut back to what they actually declare:

    PR2   → d (Domain), g (Genus), s (Species);  everything between = NA
    SILVA → d (first field), s (organism name);  everything between = NA

Plus PR2's unresolved-lineage placeholders (`Craniata_X`, `*_XX`, `*_sp.`) are **blanked,
not dropped** — the nine levels are positional and removing a field would slide Genus
onto Family. `SINTAX_FORMAT_VERSION` → v4 to rebuild the caches, which still held v2's
`p:Opisthokonta`.

**Why SILVA's `g` looked right for so long:** for `…;Sphaeropleales;Monoraphidium;
Monoraphidium sp.` the level before the organism genuinely IS the genus. It is the genus
often enough to look correct and not often enough to be correct; SILVA declares no ranks,
so the position cannot tell that case from the cod's. The existing test asserted the
happy case and is now updated to pin the contract instead.

**Verified on run 20260716T080902Z** (re-ran clustering→mapping, taxonomy and both metric
builders):

- (c) `taxa_table_class` and `tax_table` now agree (both "no class"). Because the run
  legitimately has no class any more, that check is vacuous on this data — so the alias
  fix is pinned by a positive unit test (`Actinopteri,150`) rather than by the run.
- (d) Not just Opisthokonta/Craniata gone: **Amorphea, Obazoa and Holozoa too**. Zero
  invented middle ranks from either source.
- (a)/(b) still hold after the rebuild: everything 157 527, every `otu_id_map.size`
  equal to its `otu_table` row.

**The cost, stated plainly.** The consensus is now Kingdom + Species and little else,
because the only two sources that ran do not declare the middle ranks. That is honest
rather than good: the Linnean ranks must come from a Linnean source, and NCBI — which
has one — is contributing nothing because of review item #5 (`source_status` says `ok`
on 0 assignments). Restoring phylum/class/order/family properly means either fixing NCBI
(next) or giving PR2/SILVA their own rank columns in the consensus schema (Supergroup/
Division/Subdivision, and SILVA's rank file `tax_slv_*.txt`), which is a schema change
across taxonomy_runner/database_metrics/community_analysis/phyloseq. Tests: 376 pass.

## 2026-07-16 — Single source of abundance truth (review #2), on run 20260716T080902Z

**Provenance first.** Two accountings existed, and the reference run reconciles them to
the read:

    pooled_filtered.fasta                          157 537 reads
      --derep_fulllength --minuniquesize 2   →  "10 527 clusters discarded"
      derep/cluster/otus.fasta ;size= sum     →     147 010   ← CLUSTERING basis
      --usearch_global --otutabout            →     157 527 mapped (+10 no-match)
                                                             ← MAPPING basis (canonical)

    157 537 − 10 527 = 147 010 ✓        157 527 + 10 = 157 537 ✓

So the "~10.5k reads with no trace" were never lost: `--minuniquesize 2` drops 10 527
singleton *reads* before clustering, and `--usearch_global` then maps the full pooled
file — singletons included — back onto the OTUs. Review items #2 and #7 are one bug.

**Only two files were on the wrong basis.** `otu_matrix`, `alpha_diversity`,
`rank_abundance`, `taxonomic_resolution` and every `taxa_table` already derived from the
mapping table (157 527, verified). The offenders:

- `analysis/pooled/otu_id_map.csv` — `_write_otu_id_map` read `;size=` out of the
  `otus.fasta` header, i.e. the `--cluster_size --sizeout` value.
- `analysis/metrics/otu/chimera_stats.csv` — `chimera_stats` used
  `sum(cluster_read_counts(clusters.uc))`. Its own docstring claimed that sum was "the
  reads actually clustered (157k)"; it is 147 010, and that stale assumption is what let
  the bug survive an earlier fix.

**Fixes at the source.** `_write_otu_id_map(otus, otu_table, out)` now takes `size` from
the mapping table and is called *after* the mapping step (it used to run before that
table existed — which is exactly why it could not have used it). The clustering value is
kept as a new `cluster_size` column rather than dropped: the gap between the columns is
where read redistribution shows up, and it is now visible instead of contradictory.
`chimera_stats(..., matrix_path)` quotes the canonical total.

**One subtlety, handled rather than papered over:** chimeric variants are removed
*before* mapping, so their reads are absent from `otu_table` by construction — the
cluster map is the only place they are counted. The denominator therefore depends on the
basis: on `otu_table` it is `mapped + chimeric`; on the cluster fallback the sum already
contains the chimeras and adding them back double-counts. An existing test caught me
doing exactly that.

**Verified by re-running pool/cluster/map + both metric builders on the run:**

- (a) `otu_table` = `otu_matrix` = `chimera_stats.total_reads` (basis now `otu_table`) =
  `taxonomic_resolution` = all 7 `taxa_table_*` = **157 527**, identical.
- (b) every `otu_id_map.size` equals its `otu_table` row sum; total 157 527. The new
  column shows why they used to disagree: OTU1 clustered 132 828 → maps 94 205;
  OTU2 clustered 14 144 → maps 63 260.

**Still failing, by design — later iterations:** (c) `taxa_table_class` = "Unclassified"
vs `tax_table` Class = Craniata; (d) Phylum=Opisthokonta / Class=Craniata; (e)
`source_status` reports BOLD/NCBI/PR2 `ok` with 0 real assignments and custom BLAST
`no_hits` on a 0-byte raw file. Root cause of (c) already located: `taxa_tables_all_ranks`
iterates `database_metrics.TAXONOMY_RANKS`, which spells the rank `class_`, while
`community_analysis._RANK_COLUMNS` is keyed `class` — the lookup misses, returns `{}`,
and every read falls into "Unclassified". It hits `class` only, because it is the one
rank whose name carries the underscore. Tests: 368 pass.

## 2026-07-16 — Missing deliverables in Downloads (plan #5)

**What was actually wrong.** Mostly *exposure*, not production. `taxa_tables_all_ranks`
(7 ranks × absolute+relative), `write_biom` and `alpha_diversity`/`rarefaction` were all
already written to disk by `build_database_metrics` / the Analysis metrics — they simply
were not in `_EXPORTS`, the declarative list that drives the Downloads table in
`ui/callbacks_results_summary.py`. So they existed in the run directory and in the run
.zip, and were invisible in the UI.

The one thing genuinely absent was phyloseq. Note also that recent runs had **no**
`taxa_tables/` or `table.biom` at all — not a bug in those writers, but
`build_database_metrics` returning early because there was no `taxonomy_consensus.csv`:
every reference had been skipped for lack of RAM. Fixed upstream by plan #3, same day.

**New: `write_phyloseq_bundle`** (`databases/database_metrics.py`, next to `write_biom`).
Writes `databases/<set>/phyloseq/` with `otu_table.tsv`, `tax_table.tsv`,
`sample_data.tsv` and `load_phyloseq.R`. It is *not* redundant with BIOM: PARAFISH's BIOM
writes `columns: [{id, metadata: None}]`, so the sample groups and the control flags —
the two things every downstream comparison needs — do not survive the BIOM route. They
survive this one. Counts are exported raw: rarefying/normalising is the analyst's call
and phyloseq has `rarefy_even_depth`/`transform_sample_counts` for it. Unassigned ranks
are blank rather than "unclassified", which `tax_glom()` would otherwise treat as one
real taxon and merge unrelated features under.

**Verified in real R 4.5.0** (phyloseq itself is not installed here, so the script's
final `phyloseq()` call is unexercised): the script parses, and the three tables load
with the alignment phyloseq silently depends on — `rownames(otu) == rownames(tax)`,
`colnames(otu) == rownames(sam)`, blanks arriving as `NA`, `sam$reads == colSums(otu)`.
A mismatch there does not raise; it yields an object with fewer taxa than the run had.

**Downloads (`_EXPORTS`) gained six rows:** alpha diversity, rarefaction, relative
composition at genus and species, BIOM, phyloseq. Missing rows keep the existing
greyed-out-with-a-reason treatment; taxonomy-derived rows get a reason that says "no
taxonomy" rather than "Databases did not run", which would send the user to fix the
wrong thing. The phyloseq row links the script alone and says so — its three tables are
siblings on disk and travel in the run .zip (confirmed: all four land together in
`<run>/databases/metrics/<set>/`, which is what the script needs to run in place).

**Known gap, pre-existing:** `_EXPORTS` hard-codes the `otu` feature set, so a DADA2-only
run shows these rows as missing. It predates this change (the consensus/hits rows already
did it) and is untouched here.

Tests: `PhyloseqBundleTests` in `tests/test_database_metrics.py` (10); 364 pass.

## 2026-07-16 — Marker-filtered references: SILVA/PR2 run locally again (plan #3, option 1)

**Problem.** This box has 2.8 GB of RAM and ~0.7 GB free. VSEARCH `--sintax` loads the
reference plus a k-mer index (~3x the FASTA), so full SILVA asked for ~2.5 GB and PR2 for
~1.0 GB; the memory gate in `core/preflight.py` skipped both on every run. That left only
the remote sources — BOLD (COI-only, so no hits on 18S) and NCBI (~1 min/sequence, and it
returned a species with no lineage above it). The local, fast, full-lineage path was dead.

**Change.** New `databases/reference_filter.py` cuts a reference down to the region the
run actually amplified, before the memory gate ever sees it:

    SILVA 824 MB --(keep Eukaryota)--> 110 MB --(cut to the V4 amplicon)--> 27 MB   (29x)
    PR2   334 MB --(already euk-only)--------> --(cut to the V4 amplicon)--> 85 MB  (4x)

The amplicon is cut with the run's **own** primers (`PRIMER_PRESETS` / `primer_fwd`+`rev`),
so reference and query describe the same stretch of gene — which is also why this makes
SINTAX *more* accurate, not just smaller: a 380 bp query is scored against 380 bp
references instead of 1 800 bp whole-SSU sequences whose bulk the query never spanned.

Wired through `run_taxonomy_pipeline(marker=, primers=)` from both call sites
(`callbacks_progress`, `callbacks_databases`); `preflight_notes(marker_filter=)` now says
"will be reduced to your marker" instead of the (no longer true) "will be skipped".

**Verified on the real run** `runs/pr/20260709T083137Z` (12 OTUs), with 300 MB free:
SILVA and PR2 both `ok`, 12/12 classified each, full kingdom→species lineage,
`Reliability: high`, 7.8 s end-to-end, 258 MB peak RSS. Identified *Gadus morhua*
locally — no network. Build costs ~4 s and 29 MB RSS, cached next to the source and
reused thereafter (cache keyed by source mtime + primers + format version).

**Degrades, never fails.** Unknown marker / missing cutadapt / primers that match almost
nothing (< 20% yield ⇒ wrong marker for this database) each fall back one rung:
lineage+amplicon → lineage-only → the original file (gate skips it, as before). A
mixed-marker run gets no filter — there is no one amplicon to cut to.

**Bug found and fixed on the way.** `_parse_sintax_assignment` stripped whatever sat in
the last parentheses as a bootstrap. 7.7% of SILVA's eukaryote names carry a common name
there (`Gadus_morhua_(Atlantic_cod)`), and the cutoff-filtered column it reads by
preference has no bootstraps at all — so it amputated the species to `Gadus_morhua_` and
set the bootstrap to the string `Atlantic_cod`. Now only a numeric parenthetical is
treated as a bootstrap (`_split_taxon_confidence`). Latent until now precisely because
SILVA was always skipped.

**Notes.** SILVA is filtered to Eukaryota, so bacterial contamination in an 18S run will
no longer be called by SILVA — that is the intended trade for a reference that fits, and
NCBI still covers it. Tests: `tests/test_reference_filter.py` (16) plus
`SintaxParentheticalNameTests`; 354 pass.

## 2026-07-15 — Language uniformisation (EN) + shared bar theme (plan #4)

**Language.** Every user-facing string in `ui/pages/*.py` is now English; the pages had
PT/EN mixed on the same screen ("Parâmetros usados" next to "Read tracking"). Converted
titles, card headers, subtitles, help notes, placeholders, radio/checklist labels and
buttons across Results, Databases, Analysis, QC, Preprocess and Phylo. Verified zero
accented/PT strings remain in the page layer.

**Safe by construction.** Only display text changed. Dropdown/radio **values** were
already English internal keys (`Consensus_Kingdom`, `otu`/`asv`, `shared`/`unique`,
stage keys) and were left untouched, so no callback logic depends on the translated
labels. All six page layouts render; 333 tests pass.

**Shared bar theme.** `ui/plots.py::_theme` — the single styling entry point every figure
already goes through — gained `bargap=0.32` / `bargroupgap=0.12`, so bar charts read as
bars instead of one solid wall. Ignored by non-bar figures, so it is safe to apply to
every figure and needed no per-page changes.

**Header typography.** `.page-sub` dropped `font-family: var(--font-mono)` in
`assets/parafish_theme.css`, so page descriptions render as a normal sans subtitle
rather than monospace.

**Still open for this item.** Long multi-line help paragraphs in the *callback* layer
(empty states, log/warning text) remain Portuguese — deliberately deferred: they are
lower-visibility and riskier to bulk-translate. A central string constant/dictionary
(the plan's anti-drift measure) is not yet introduced.

## 2026-07-15 — Primer removal with cutadapt (plan #2)

**Problem.** The reproducibility block showed `primer_fwd`/`primer_rev` empty and primers
were only ever handled by a fixed-length `vsearch --fastq_stripleft/stripright` in the
Analysis stage — blind (strips N bases whether or not the primer is there), no degenerate
matching, and a no-op unless the user filled the fields.

**Change.** New `preprocess/cutadapt_runner.py` removes the forward primer from R1 and the
reverse from R2 on the *raw* reads, before Trimmomatic, using cutadapt's native IUPAC
wildcard matching. Wired into the Preprocess loop in `ui/callbacks_progress.py`: cutadapt
runs, then the primer-trimmed reads are fed to Trimmomatic via `dataclasses.replace` on the
sample's paths. Presets (`analysis.vsearch_runner.PRIMER_PRESETS`: 18S V4 TAReuk, COI
Folmer/Leray) and the manual fwd/rev fields already in Setup are reused. `CUTADAPT_PATH`
added to `core.config`.

**No double removal.** When cutadapt ran, the pipeline sets
`vsearch_params["primers_removed_upstream"]`, and `_primer_strip_lengths` returns `(0,0)` —
so the Analysis-stage length strip never truncates real sequence on top of cutadapt.

**Degradation.** No primers configured → reads pass through untouched (`primers_removed
== 0.0`), behaviour identical to before. cutadapt missing or failing → same passthrough
with a clear warning, and the old length strip still applies as a fallback, so primers are
not silently kept.

**Read tracking.** cutadapt writes a `--json` report per sample; `preprocess_metrics` adds
`reads_with_primer` + `pct_primer` (how many reads actually carried the primer).

**Verified.** 18S V4 forward with degenerate S/Y removed exactly from concrete-base reads;
PE removes both primers; passthrough and cutadapt-absent paths confirmed. `cutadapt>=4.0`
added to requirements (installed 5.2). Tests: `tests/test_cutadapt_runner.py` (8). 333 pass.

## 2026-07-15 — Active-project visibility on Samples/Setup (project-tab UX)

**Problem.** A user created a project "new", but the samples and the run were attributed
to a different, still-active project ("projeto_1"). Audit found no attribution bug — the
code correctly saves to the active project and there is no code path that reverts
`active_id`. The real gap was UX: the active project showed only in the corner badge, so
it was easy to create a fresh project, upload into whichever one was still active, and run
the wrong data without noticing.

**Change (additive, no logic redefined).** A "Projeto ativo" strip at the top of Samples
and Setup (`samples-active-project`, `setup-active-project`) naming the active project and
its selected-sample count, and a line next to "Run Pipeline" (`setup-run-target`) stating
exactly which project and how many samples the run will use — or an amber warning when no
project is active / no samples are selected. Helpers `_active_project_banner` /
`_run_target_note` in `ui/callbacks_project.py`, driven by `store-project` + `store-samples`.
Reuses existing badge/pill classes (no CSS changes). Run logic untouched. 325 tests pass.

## 2026-07-15 — Multi-rank taxonomy: resolve NCBI lineage from taxids (plan #1)

**Symptom.** On a run whose hits were all NCBI BLAST, the consensus table had
`Consensus_Kingdom..Consensus_Family` blank and only Genus/Species filled. That
emptied the Databases sunburst (all inner rings "Unassigned"), the per-rank
breakdown (kingdom..genus = 0, species = 100%), and the Results community-composition
plot (a single solid Phylum bar). The "Nível taxonómico" selector had nothing above
species to offer.

**Cause.** `blastn -outfmt 6 ... sscinames` reports only the organism binomial. The
code recovered Genus (first word) + Species from it but never the lineage above genus.
BOLD would carry a lineage but found nothing for these sequences, and SILVA/PR2 (which
do) were skipped for RAM.

**Fix.** New `databases/ncbi_lineage.py` resolves each hit's `staxids` into a ranked
lineage via `Bio.Entrez.efetch(db="taxonomy")`, batched over unique taxids and cached
on disk (`<DB_ROOT>/ncbi/lineage_cache.json`) so re-runs never re-fetch. Wired into the
NCBI branch of `run_taxonomy_pipeline`, filling only blank ranks (a name BLAST already
asserted is never overwritten). NCBI's `superkingdom` and `kingdom` both map to the
`kingdom` column, with `kingdom` (Metazoa/Fungi/...) winning over `superkingdom`
(Eukaryota) as the more informative top level.

**Config.** `core.config.get_entrez_config()` reads `ENTREZ_EMAIL` (required by NCBI for
automated access) and an optional `NCBI_API_KEY`/`ENTREZ_API_KEY` (raises the rate limit
3→10 req/s). Both from the environment/.env.

**Degradation.** Best-effort by design: no email, no Biopython, or no network → returns
the cache plus a plain-language warning appended to the run log, and the upper ranks
stay blank instead of failing the stage. Verified on the real run's `blast_result.tsv`:
taxids `81636/100272/310288/2060155` captured, and with lineage supplied OTU1 resolves
Metazoa→Chordata→Actinopteri→Gadiformes→Gadidae→Micromesistius→M. poutassou.

**Note for existing runs.** The fix populates *future* runs; a past run's consensus CSV
is not rewritten. Set `ENTREZ_EMAIL` in `.env`, then re-run the Databases stage to see
the ranks fill. Tests: `tests/test_ncbi_lineage.py` (7).

## 2026-07-14 — UI unification, Results rework, ordination (T0–T4)

### T3 — the selector bug: what was actually wrong

The spec's first hypothesis ("no callback populates `options`") was **wrong** — the
callbacks existed (`ui/callbacks_results.py`). The real causes were the other two,
compounding:

1. **Options were gated on `status == "success"`.** Any stage that failed, was skipped,
   or ran partially never appeared in the list. On the run that actually needed
   debugging — analysis and databases both failed — the dropdown was literally empty and
   read as a dead control rather than a failed run.

2. **The source of truth was UI state, not disk.** Options derived from
   `store-{qc,preprocess,analysis,databases,phylo}-outputs`, all `storage_type="session"`.
   Those stores are only ever written by the progress poller, which *disables itself* when
   the run finishes. So: app restart → empty; new tab → empty; navigating to Results
   without having watched the run go by → empty. Meanwhile `<run>/progress.json` on disk
   had everything, and was never consulted.

**Fix:** `core/run_registry.py` reads the run from disk and is now the single source for
every selector. All five stages are listed with their real status (`Analysis ✕ falhou`),
and stages with no per-sample output say so (`"Este stage não tem output por amostra"`)
instead of offering an enabled, empty dropdown. Survives F5, new tabs and a server
restart, because the facts live where the pipeline put them.

### T4 — zero replacement in the CLR: the decision, and its cost

**Chosen: a fixed pseudocount of 0.5 on the counts**, applied after the prevalence
filter. Rejected `zCompositions::cmultRepl` (Bayesian multiplicative replacement) because
it pulls in an R dependency for a second-order refinement.

The cost is real and is documented in the code, in the UI, and in a test
(`test_pseudocount_breaks_exact_scale_invariance`): **`clr(x)` is exactly scale-invariant,
but `clr(x + 0.5)` is not**, because 0.5 means something different at a depth of 100 than
at 5,000. Measured distortion is <1% on the coordinates for a 50× depth swing — orders of
magnitude below the biological signal, and far smaller than the depth artefact it
replaces (on raw counts, PC1 captures sequencing depth almost entirely). But it is
non-zero, and the test exists so nobody later "fixes" it by tightening a tolerance.

Other guards worth knowing:
- `n < 3` → **no plot**, ever. An explicit empty state instead: with two samples there is
  one axis and the separation is an artefact of having two points.
- `3 ≤ n < 5` → plots, with a visible low-confidence warning.
- Degenerate distance matrices are refused (Jaccard on a table where every feature is in
  every sample gives all-zero distances; skbio then returns NaN axes).
- PERMANOVA reports its **p-value floor**. With two groups of three the smallest
  attainable p is 1/10 = 0.1, so a "p = 0.12" from that design is the floor of the test,
  not weak evidence — it now says so.

### New dependencies

`scikit-bio` (PCoA, Bray-Curtis/Jaccard/UniFrac, PERMANOVA) — already present in
`requirements.txt` and installed (0.7.3); no new pins were needed. `numpy`/`scipy` cover
the CLR-PCA via SVD, so no `sklearn` dependency was added.

### Other changes

- **T0.** One `data_table()` factory (`ui/components/table.py`); 20 direct
  `dash_table.DataTable` instantiations removed. QC's "Run overview" and Preprocess's
  "Read tracking" had passed only `{"padding": "8px"}`, inheriting Dash's near-black
  default on a dark card — ~2:1 contrast, effectively invisible. Now 13:1, enforced by a
  test that computes the WCAG ratio and by one that greps for stray `DataTable(` calls.
- **T1.** Phylo and Results moved onto the shared cards/KPIs/section headers. Error toasts
  consolidated into one banner with a count (`Servidor não respondeu · ×9 · 09:09:13`) and
  a `<details>` for the occurrences; the Dash dev-tools overlay that stacked ten identical
  cards is now off. Pages stratified into Layer 1 (status + KPIs + one hero plot, always
  visible) and Layer 2/3 (`<details>`, zero callbacks). Raw Illumina headers demoted and
  the ids cleaned. Distribution plots switch to labelled bars at `n ≤ 4`.
- **T2.** Results rebuilt as a run summary: run header with an aggregate badge, an
  end-to-end **read funnel** (input → trimmed → merged → clustered → classified, with any
  single step losing >20% flagged red), five clickable per-stage KPI cards, community
  composition with a rank selector and a host-exclusion toggle, compact ordination, then
  tables/exports/reproducibility. Stages that did not run render an explicit badge, never
  a blank.
- **Reproducibility block** reads tool paths from `core.config`, **not `$PATH`** — on this
  machine those differ (PATH has the broken blastn 2.9.0; the pipeline uses 2.17.0), and
  reporting the wrong one would put a version in a paper that never touched the data.

313 tests pass.
