# Issue #598 — Final implementation handoff

S04 completes the existing implementation branch
`fix/issue-598-alignment-direction-reset-20260926` for
`satoshikawato <kawato@kaiyodai.ac.jp>`. The exact final commit SHA is reported
in the delivery response, avoiding self-reference in this file.

## Completed behavior and authority

The exact clicked reference and selected known-direction target anchors use
one exclusive Keep/right/left/Custom selection. A whole record reverses while
source bytes, identity and biological strands stay unchanged. Unknown, Skip,
missing and unusable anchors keep direction with visible reasons. A minority
reference can reverse independently. Its immediate pre-Align logical x and all
logical y positions are preserved; automatic composition/viewBox fitting and
screen pixels are separate.

Local choices start no Worker work. Explicit Apply validates once per attempt;
changed final facts require another acceptance. Failure retains editable intent
and the previous complete artifact. Reset defaults to positions only; combined
Reset restores absolute before-Align directions only for actual latest-Align
changes, including an affected reference and later manual edits on those
records. Both scopes consume plan/evidence; Undo restores them before another
scope. Historical absence differs from a valid empty change list. Fresh Load
retains admitted evidence and both scopes add zero LOSAT jobs. Full History,
pending form, readiness/finalization rollback and exact ribbon placement remain.

| Accepted concern | Actual dev authority | Preserved contribution |
| --- | --- | --- |
| PD-OI-027 | revision 5 | Exclusive scope, whole-record transform, source invariants and logical reference geometry. |
| PD-OI-029 | revision 3 | Both Reset scopes, absolute sparse receipt, later manual edits and complete History. |
| PD-OI-031 | revision 5 | Exact popup/drawer reference, automatic Keep, usable deterministic anchors and discoverable review. |
| PD-OI-034 | revision 5 | Local choices, explicit final validation, refreshed preview, rollback and editable retry. |
| PD-OI-039 | revision 2 | `EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH`: retire only the obsolete Match affordance/checkbox/flag; preserve its independent operation contributions. |

All nine receipt fields match the approved Packs against actual
`origin/dev@c922fc38aac78da9be83342c09ac0164ecef6ff6`. PD-OI-035 revision 3 remains
independently necessary, with canvas, keyboard, focus and Editor behavior
preserved. Contract SHA-256 is
`5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`.
Candidate authority is not used to authorize runtime. No new material Product
choice or BUG-17 behavior was needed.

The final dev advance adds independent PD-OI-046/047 authority and metadata/legend
architecture guards only. These retain the five complete Issue #598 receipts
and the independent PD-OI-035 section. The changed guards have a separate S04
run; no #600/#601 runtime implementation is introduced.

S04 normally integrates that inspected dev as the single commit's second
parent, retaining S00–S03 and prior writers' work. Its runtime is byte-identical to tested dev@9d967f1c; production equals that
runtime, apart from Gallery prose; the browser-test differential is
`af7aed11a44a93be878d5e3e8d9001cfba5b374ac90c106a10cbde3b9c8521a6`.
S00 `302dfa1136c95ab50ddef606ad835ec4087b95b2`, authority-only merge
`2edc00aebc74e01003da643dfc957b513d5dcfe5`, prior dev
`aab5ad4323d51fc8ce94577b900b1c35e5ad7ffe` and actual dev are required ancestors.

## Evidence and changes

[The S04 evidence](evidence/S04.md) records every A01–A16 condition separately,
owner/path differential, retired paths, fixtures, environment, commands,
observed geometry, limits and delivery verification. A01–A16 have the required
native, Node, real browser and documentation evidence within the stated scope.
The complete supported-version matrix and actual assistive-technology speech
are not claimed.

Existing Web/Session/Python references, release notes and the five-BGC
Tutorial/Gallery explain the finished behavior. No new public page is added.
The existing capture owner replays five original MIBiG GenBank inputs from an
empty GUI, produces six Tutorial images and checks optional minority-reference
directions, both Reset scopes and Undo. All five source files are byte-identical
to the authoritative version-5 downloads. Final figures retain labels, legend,
color rules, metadata, coordinate ruler and comparison context. All six images
were visually compared with the prior same-size images. Existing Gallery
popup/final media remains truthful and is retained; internal smoke/QA images
are not published. The owner-maintained social preview and tracked reference
SVGs are unchanged.

The only new acceptance-test contribution restores S03's reference-only Custom
center/pending-form test and strengthens retained canvas/error-clearing checks
alongside dev's current tests. No resolver, render/History/receipt/admission
owner, persistence policy/version, migration, privileged path or dependency is
added. Production, tests, docs and generated diffs were reviewed separately.

## Observed checks

| Check | S04 result |
| --- | --- |
| Required Python `pytest tests/ -v -m "not slow"` | 6268 passed, 17 skipped, 11 deselected; 1221.53s. Both real comparison browser shards pass. |
| Node | 712 passed on runtime@9d967f1c; latest-dev guard additions verified separately below. |
| Final architecture/Product guard boundary | 192 passed after latest dev integration, including both new guard cases. |
| Alignment + related History/palette browser | 14 passed. |
| Existing Rotate and region browser | 2 + 2 passed. Existing linear crop/reverse/fresh Load runs in comparison shards. |
| Gallery browser / media | 23 passed; 19 media / 19 operations checked. |
| Read-only reference SVG comparison | 16 passed before other tests; no regeneration. |
| Final documentation contracts | 41 passed after the final harness changes. |
| Final raw-input GUI recipe | PASS: 6 PNGs, 155 features; SVG 249,702 bytes; TSV 232 rows. |
| Ruff, local public links/images, JSON/JS syntax, manifest, whitespace | PASS. |
| Working-tree Web/Product/architecture policy | Gate PASS; Review CLEAR; ordinary, four hard architecture rules CONFORMING. |
| Mandatory post-commit policy against dev@c922fc38 | Gate PASS; Review CLEAR; ordinary, no blocker/review reason, four hard rules CONFORMING. Repeated on final amended HEAD. |

Regenerate the public example with:

```sh
python docs/capture/run_all.py --scenario T-GUI-04 --tier extended
```

[Capture README](../../capture/README.md) pins the existing environment and
source verification. S04 evidence lists all executable checks and isolated
checkout/venv/wheel/server details. Broad checks predate only final capture
sequencing/scroll adjustments; final contracts and real GUI replay cover that
changed boundary. Runtime, acceptance tests, fixtures and baselines are unchanged.
Slow/performance suites, remote CI, every Python/browser version, physical zoom
and assistive-technology speech are not newly executed. Local Chromium evidence
is not release, staging or comprehensive accessibility certification.

## Delivery and next boundary

The single commit uses the requested author/title. Mandatory post-commit check:

```sh
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
```

Gate and Review are separate: CLEAR does not grant publication/merge approval.
The prior S03 Review REQUIRED is not relabeled as approved. Immediately before
the normal push, fetch the same-name branch, confirm remote is an ancestor of
HEAD, then verify remote SHA=HEAD and clean tree using the commands in S04.
Unexpected advancement requires writer evidence and integration before retry;
force-push is forbidden.

Any later session must create a fresh dedicated checkout, inspect same-name
remote advancement/writer completion and actual latest dev authority, preserve
this work, and verify its own exact integrated candidate as needed. This task
authorizes only the named work-branch push. PR publication/merge, dev/main push,
deploy, tag and release remain separate unexecuted boundaries.

English commit title: **Document and verify alignment direction and reset behavior**.

Summary: Document exclusive alignment directions and selectable Reset scopes,
reproduce the five-BGC GUI workflow, and verify source, geometry, persistence
and History behavior on the integrated runtime.
