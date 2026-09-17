# Local member, reciprocal-best-hit and evidence-rank reductions

2026-09-16. The user explicitly reopened these local owners beyond the old
S07.5 deferral. This is not S08 or execution of the entire plan. After implementation
and correctness verification, the user authorized benchmarks. The combined final
source has five passing native timing comparisons and one inconclusive case;
traced memory peaks increase. See the exact boundaries and reservations below.

[S07.7](S07_7.md) owns the common source/inheritance ledger, baseline and final
SHA-256, authority preflight, integration checks and session commit handoff.
All new production changes are in `gbdraw/analysis/protein_colinearity.py`.
`collinearity.py` was inspected and needs no modification. B follows A, and C
follows the B gate; cumulative patches and hashes preserve the three stages.

## B: prove the discarded work before deleting it

The finite member consumer calls public `select_top_hits_per_query` to choose
distinct pairs by numerically coerced rank, then filters the **original** HSP
frame by pair membership. Collinear OFF passes these original rows to its
rbh/one_to_one/all edge selector; ON and Similarity pass them to normalization.
Therefore their later numeric coercion cannot be skipped: for example raw
string bitscores `"9"` and `"100"` would rank differently as strings. Public
selectors still validate independently and retain their existing empty and
exception boundaries.

| Discarded work | Why the same result remains |
| --- | --- |
| Unlimited member `hits.copy()` before `reset_index(drop=True)` | Reset already returns a detached frame with the same RangeIndex and dtypes; scalar mutation tests verify input isolation. |
| Two query/subject projection DataFrames before `MultiIndex.from_frame` | `MultiIndex.from_arrays` consumes the same two Series, including categorical domains. Both required indexes and the original-row membership mask remain. |
| Pair `drop_duplicates` before best-query/subject `drop_duplicates` | The first row in each sorted query/subject bucket is necessarily the first row of its pair. Removing later copies of that pair cannot change that first bucket row. One dedup pass disappears from best-query and two from one-table reciprocal selection or a two-direction RBH call. |
| Full-width namedtuple conversion in reciprocal selectors | These consumers read only query/subject. Zipping the two Series retains row order, the same `str` conversion, reverse-map last-write order and reciprocity test. |

The distinct-pair deduplication inside public **top-K** is retained. With rows
`q→a (100), q→a (90), q→b (80)` and K=2, dropping it would consume both slots on
`a`. Every HSP for the selected pairs is still returned in original row order,
with original column dtypes. No index, rank or normalized table crosses a new
prepared-result boundary. Sorts and validation that are still necessary remain.

B gate: **532 passed, 1 existing external-smoke skip** in
[B tests](data/s07-7-B-tests.log). Independent ordinal oracles cover finite
member 1/5 and unlimited, full/partial ties, duplicate pairs, multiple HSPs,
numeric strings, string/numeric/categorical IDs, reversed evidence and category
order, duplicate indices, original dtypes and input immutability. Public numeric
rejection checks all ten numeric fields against five invalid values, including
both directions, with exact exception type/message. Empty frames, missing
columns, invalid limits and validation order also match the frozen source.
The edge-caller matrix covers all three anchor modes and missing reverse tables.

Structural assertions observe one fewer full-frame copy for unlimited members,
two fewer two-column projection frames for finite members, the eliminated
dedup passes and no full-row `itertuples` calls in reciprocal selectors. No
claim is made that the two required MultiIndexes themselves disappeared.

## C: one stable order, distinct membership phases

After `_dedupe_anchor_core_directional_rows`, the selector now sorts the same
`best_by_direction.items()` once by `_anchor_core_hit_rank`. It passes that list
to anchor selection, record-local selection and the related-edge phase. The
original mapping is preserved in its original insertion order for support
accumulation and other consumers.

The rank depends only on normalized score, evalue, coverage, identity,
alignment length, endpoint IDs and record indices. Normalized rows are immutable
namedtuples with validated numeric fields; proteins and their record indices
are unchanged throughout these phases. Group assignment changes none of those
fields. Python's stable sort uses the same original mapping order, so ties have
the same order at each previous sort site. Membership checks and every phase's
full scan remain; core snapshot and expanded-membership snapshot are separate.

The sole consumer of `_best_rows_by_query_target_record` read only `rows[0]`
to get best subject and score; it never read second or later rows. Its replacement
`_best_row_by_query_target_record` retains the first cross-record row for each
query/target-record key from the shared stable order. Filtering a sorted sequence
to a bucket preserves the minimum and its stable tie order. This removes bucket
lists and their re-sorts. No global min/top2 replacement or rank dictionary was
introduced. The best-subject and best-score mappings remain lookup-only; their
construction order does not determine group or edge order.

For E deduplicated direction rows and X cross-record rows, the three global
sorts become one; exactly **2E + X** rank evaluations disappear, including the
bucket ranks. Other ranks, final edge/group/member ordering, dedup ranks and
support tie-breaks remain. Tests instrument the unchanged rank function and
global sort boundaries on sparse, dense and support fixtures; they assert that
formula and complete output equality. An independent brute-force bucket oracle
uses its own field tuple and verifies first equal-rank retention.
[Structural-only evidence](data/s07-7-structural.json) records these observations
without clocks, profiling or memory measurement:

| Fixture | E / X | Global sorts before → after | Rank calls before → after |
| --- | --- | --- | --- |
| sparse-2 | 4 / 4 | 3 → 1 | 16 → 4 |
| dense-4 | 64 / 64 | 3 → 1 | 260 → 68 |
| support-edges | 12 / 10 | 3 → 1 | 56 → 22 |

All three full ordered-result hashes match the frozen source.

Temporary retention is explicit: one O(E) list of `(direction-key, row)` pairs
lives in the outer selector until return, instead of three independently
allocated phase lists. The row objects already remain reachable through
`best_by_direction` and normalized tables. Sort key tuples are temporary to the
single sort; no O(E) rank-tuple dictionary survives it. One row per occupied
cross-record bucket replaces its list of rows. The shared traversal remains
live across group construction, which can change peak overlap with other
containers. The final combined measurements show 2.54–4.80% higher peaks on
ON/Similarity; they do not demonstrate an isolated C allocation delta. No global
memory reduction is claimed.

The final local/support gate is **718 passed**; whole-result native and typed
gates are in [S07.7](S07_7.md). Existing fixed-core/expanded-core tests still
prove that an unassigned member cannot chain through a newly attached member,
while that member does contribute to later record-local competition. The
complete differential covers group/member/role/confidence/support/edge fields,
IDs, ordered mappings/tuples, every output table and typed payload. Float
arithmetic and its input traversal are unchanged; comparisons use no tolerance.

## Rejected alternatives and scope limits

- Skipping downstream numeric conversion is invalid for original string HSPs.
- Removing top-K distinct-pair dedup is invalid for the duplicate-pair example.
- Sorting the direction **dict itself** would also reorder support accumulation
  and could change exact floating results; the original dict is kept.
- Replacing the global traversal with min/top2 is invalid: later phases inspect
  all qualifying local/related evidence, not just anchor candidates.
- A global rank cache would retain E nine-field tuples; the shared sorted list
  removes the proved duplication without that data or a new cache owner.
- Combining core and expanded membership indexes would permit invalid chaining
  or miss later competition. Both snapshots remain.
- Regression/vectorized power, numeric validation redesign, S05 prepared classes,
  persistent caches, Worker transactions, scheduler/runtime changes, mode-specific
  fast paths and validated flags remain out of scope.

## Final measurements and cumulative workflow question

The user later authorized measurement. [S07.7](S07_7.md#authorized-measurements)
records six exact S07.6-final/current post-search comparisons, source guards,
commands, all samples and host snapshots. Five timing pairs pass the unchanged
numerical gate; Vibrio ON is inconclusive (initial +10.17%, 9.51% current noise).
Its planned 21-sample repeat was stopped at the user's request; a three-sample
pair gives 3.080 → 2.876 seconds (−6.63%), descriptive only. The 21 was borrowed
from the Hep plan by agent judgment, not a calculated sample size; Vibrio's
existing count is seven. No favorable samples replace the original observations.

Hep ON/OFF changes are −13.42% / −4.74%, Vibrio OFF −8.68%, Similarity finite/
unlimited −10.52% / −23.40%. These are **combined A+B+C native post-search**
comparisons, not per-owner causal attribution or whole-Generate percentages.
Tracemalloc peaks increase 2.54–4.80% for ON/Similarity and 0.02–0.05% for OFF.
Returned outputs remain live; prepared inputs are outside the tracing boundary.
The shared sorted list's longer lifetime is consistent with extra overlap, but
only an isolated allocation experiment could attribute that increase to C.

The question about S01–S07.5 requires a different denominator. The earliest S01
same-code timing failed repeatability, so it is not used as improvement evidence.
The stored S03 pre-change observations use **identical fixture inventories and
dependencies** to the current Hep measurements. Their cumulative observations
are checked in [this ledger](data/s07-7-cumulative-observations.json):

| Prepared native post-search | Pre-S03 seconds | Current seconds | Observed seconds removed | Observed reduction |
| --- | ---: | ---: | ---: | ---: |
| Hep Collinear ON | 5.725829 | 0.403596 | 5.322233 | 92.95% |
| Hep Similarity member 5 | 8.014944 | 0.753803 | 7.261141 | 90.60% |

These are historical observations, **not a controlled cumulative acceptance**:
S03 recorded external LOSAT contention, host load differs across sessions, and
S06 intentionally replaced exhaustive path representation. Its graph/exhaustive
oracles own semantic compatibility; complete old/new typed hashes are not equal
across that representation change. The percentages above are direct arithmetic
on the two archived medians, not products of session speedup percentages.

The large earlier reductions were HSP aggregation and evidence/support work.
[S03](S03.md#native-timing--observed-with-host-contention) observed Collinear
post-search 5.73 → 2.86 seconds, and Similarity 8.01 → 3.13 seconds with its
contention reservation. The separate [S04 final pairs](S04.md#final-native-timing-results)
measured 2.56 → 0.639 and 2.76 → 1.15 seconds. The unequal S03-after and S04-before
values show why these cannot be treated as one continuous controlled run.

[S06](S06.md#combinatorial-input-and-explicit-output-control) also removed growth
proportional to every path: its R=24 fixture represents 4,194,304 paths using
24 nodes and 276 transitions, without enumerating them. No old R≥24 execution
exists from which to claim a speed ratio. Its R=16 same-source explicit/compact
control reduces typed output from 6,884,782 to 74,298 bytes and encoding from
0.355 to 0.00325 seconds. [S07](S07.md#time-memory-and-output-results) reduces the
1,200-anchor chain's merge from about 0.712 to 0.00190 seconds; real Gallery
merge contributed much less. These scaling benefits matter beyond the latest
small incremental timing difference. S05 contributes no adopted implementation;
S07.5 was a plan/prototype, not another production speedup.

**Search start → final SVG/save has no comparable cumulative measurement.**
Post-search excludes raw LOSAT, extraction, parse/filter, metadata/typed encode,
render/save and Worker startup/transport. Independent stages include nested work
and must not be added. S01's manifest-helper timing and S06's conversion-helper
timing are different boundaries and cannot be compared as the same workflow.
The S06 browser observations were about 5.4 seconds of initialization plus
2.85 seconds of cold Collinear conversion, or 5.7 plus 4.47 seconds for Similarity,
excluding final drawing and subject to the recorded host/noise limitations.
These historical observations indicate why subsecond native analysis cannot be
translated directly into browser latency. Current correctness flows prove
Generate/save/reload behavior; they did not time an old/current full workflow.

## Remaining work and handoff

A current-source diagnostic profile is recorded separately from timing for
Vibrio ON, Vibrio OFF and Hep Similarity, **one profiled call each**. Its shares
are instrumented cumulative costs, not uninstrumented time shares or promised
savings; child and parent rows must not be added. [Final profile evidence](data/s07-7-remaining-costs.json) identifies where investigation
remains useful; all three profiled scientific output digests match their timing runs.

| Current owner | Vibrio ON | Vibrio OFF | Hep Similarity |
| --- | ---: | ---: | ---: |
| Normalization, parent | 44.79% | — | 50.81% |
| HSP aggregation, child | 31.32% | — | 25.03% |
| Fit-row selection, child | 5.14% | — | 13.75% |
| Group construction | 24.83% | — | 16.65% |
| Member selection | 7.24% | 21.50% | 12.87% |
| Unit construction | 6.74% | 43.52% | — |
| OFF edge selection | — | 20.08% | — |
| Numeric coercion, nested across phases | 7.79% | 23.92% | 12.05% |

Thus normalization/HSP and group construction remain the largest inference-side
owners; unit construction remains the largest OFF owner in this larger case.
Fit and rank are already smaller: rank contributes 2.56% in Vibrio ON and 1.80%
in Similarity under this instrumentation. More tiny rank changes cannot reasonably
be assumed to deliver the old multi-second reductions. None of these percentages
proves that the corresponding work is redundant or removable.

Before further local optimization, measure the actual Generate sequence with
fixed mode, settings, raw evidence and cache state: startup, raw search when
needed, transfer/parse, inference, encode and render. Compare cold, settings-change
and exact-repeat flows separately. That measurement would establish the user
waiting-time denominator missing from the historic records. It does not require
reviving the rejected S05 class/cache/Worker redesign or assuming a new owner.
No extra production optimization or full S08 execution is included here.
Future timing defaults are now three samples by user instruction; noise does not
trigger an automatic repeat. Historical reports retain their recorded policy.

Future native candidates must preserve public numeric rejection, scalar score
order, distinct-pair limits, original HSP rows and evidence order. Remaining
normalization/group work is not automatically removable, and its profile share
is not an achievable percentage improvement. A/B/C attribution can use the saved
stage patches only if requested; historical stages need not be rerun wholesale.
The accepted S07 status (14 pass / 6 inconclusive), S06's future merge Review and
rejected S05 remain unchanged.

## Rollback and audit

Rollback B or C with its difference between the recorded sequential patches;
whole-session rollback restores only the frozen protein owner and the two
adjusted private test calls. Keep S07.6 units and its evidence. Rebuild the wheel
after rollback. No saved-data conversion is necessary.

Production, tests, documentation and generated evidence are audited separately
in the final review ledger. No production `collinearity.py`, raw identity,
request/Session schema, public Gallery asset or reference SVG was changed by
this session. The final browser/installed CLI evidence and English commit
handoff share the [S07.7 report](S07_7.md). S08 remains a separate request.
