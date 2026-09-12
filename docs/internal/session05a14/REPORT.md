# SESSION 05A14 scoped remediation

All six defect families are **CLOSED in the local candidate**, with permanent
red/green tests and complete affected original programs. They are not yet
integrated into dev. This is not a clean full-acceptance declaration.

The unchanged failing witness, owner analysis, shared comparator and reproduction
command are in [RESULT_COMPLETION.md](RESULT_COMPLETION.md). Compact file hashes,
original-program hashes, operation counts, gate logs and package provenance are
in [PROOF.json](PROOF.json). Raw evidence remains at
`/tmp/gbdraw-session05a14-evidence`; historical A/B sources, fixtures and recorded evidence were not copied or
modified. The original seed manifest is verified by hash before each retained run.

| Defect | Status | Before correction | Permanent passing assertion | Complete original acceptance |
| --- | --- | --- | --- | --- |
| B-01 / A-01 | CLOSED | feature-color Redo left selected fill stale | Circular/Linear `color Redo publishes the completed SVG`; direct edit and every History checkpoint compare separately | J09, J10 PASS, including each of four Undo and four Redo steps and remaining operations |
| B-02 / A-02 | CLOSED | placement History changed unrelated default fills to `#cccccc` | `placement History preserves every unrelated feature`; defaults, custom default, type and specific/hash precedence controls | J13 PASS; four independent non-target History checks PASS |
| B-03 / A-03 | CLOSED | mobile Generate bar intercepted Apply Label pointer at 390x844 | Circular/Linear `mobile search Edit Apply Label accepts pointer input`; nonempty durable override, regenerated/restored text and control bounds | J43, J44 PASS through original pointer paths and downstream persistence/mode steps |
| B-04 | CLOSED | depth ancestor visibility changed without selected Result publication | Circular/Linear `depth OFF and ON publish parent visibility`, observed before subsequent edits | J27, J28 PASS through intervening edits and persistence |
| B-05 | CLOSED | delayed Web replay differed at legend by 1 px / 0.5 px and root width by 1 px | current-CLI single/composite replay after observable definition completion, plus fast control and controlled owner scheduler | supplemental C01 PASS, both single and two-source composite complete replays |
| B-06 | CLOSED | replacement scale/5 kbp parent translated by 44.5 px | J36 rejected FASTA recovery and label regeneration compares selected/mounted/downloaded scale | J36 PASS with original rejection/recovery prefix and intended label confirmed |

The common observer checks all scientific SVG nodes and relevant ancestors,
rather than only biological paths. It retains completed-action observations
before the actual export boundary, so later Save/Generate/export cannot erase an
earlier failure. An unrelated feature assertion rejects coherently wrong output.
This closes the old feature-only/group-blind, late-checkpoint and coherence-only
gaps independently. There is no claim that the six cases share a cause or were
introduced together.

## Mobile correction

The existing fixed feature/pairwise popups and their scope dialogs are teleported
to the document body so preview containment cannot trap their stacking context.
Their existing positioning constraints now supply the available height below
the actual popup top. The mobile Generate bar stays below the popup/dialog
stack, normal popup-body scrolling remains available, and long header/placement
content is constrained within the popup. No forced pointer input, injected CSS,
DOM click or Enter substitute is used for Apply Label.

Final source and packaged screenshots were inspected at 390x844. Apply Label,
Close and placement controls are visible and reachable in both modes. Existing
Editor-toggle and preview navigation checks pass. This uses the existing popup,
scroll and stacking owners; it adds no new state owner or interaction framework.

## Verification and review packages

The ten permanent browser cases first failed on unchanged production at their
intended target actions (`red-confirmed.log`), then passed unchanged semantic
assertions. The specifications were subsequently separated by review scope;
all ten pass again in `review-final.log`.

- Fast Node suite: 565 passed. After separating test files, the comparator and
  style/palette/layout owner checks passed again.
- Targeted and neighboring browser suite: 32 passed, protecting independent
  draft B / Result A / committed request A, multipart labels, source-bound
  intent, manual/category legends, scale, navigation and first-attempt recovery.
- Original programs: eight journeys plus C01 PASS, zero failed or blocked
  checkpoints. J13/J43/J44 were repeated after the final observer/layout refinements.
- Independent review commits: the first passes its seven browser cases and five
  owner/comparator checks without the palette/mobile changes; the second passes
  its palette browser case and five precedence controls without the mobile change.
- Existing PR smoke: 10 passed; Gallery publication parity: all 9 passed.
- Python browser gate: 23 passed initially; the remaining downloaded CLI replay
  passed with the candidate installed in an isolated environment. The first
  attempt had selected the older global CLI (session 40 versus current 41).
- Distribution wheel built successfully. Its six changed production files match
  source hashes. Four package browser cases passed from cold contexts with
  external requests blocked, covering desktop depth OFF/ON and mobile pointer
  editing, Generate, Save and fresh Load. The package has no analytics script.
- Architecture suite: 137 passed. Web policy Gate PASS; Review REQUIRED for the
  derived identity declaration. No maintainer review is represented as complete.

Review packages, in order:

1. Publish completed SVG edits and bind new Result geometry: B-01/B-04/B-05/B-06,
   carrying the shared comparator and retained-runner adapter.
2. Preserve palette defaults during placement History: B-02, using that comparator.
3. Keep mobile feature controls reachable: B-03, using the same comparator and
   carrying this combined scoped report.

CI routing, schemas, test-owned timeouts, the ten-case PR smoke budget, generated
Gallery artifacts and reference SVGs are unchanged. Production, tests,
documentation and ignored generated packages were reviewed separately. Earlier
setup failures and interrupted exploratory runs remain in the evidence directory;
they are not the red proof or final gate results. Two disposable Python bytecode
cache files created by the initial retained import were removed; the adapter
now disables bytecode writes before importing archive programs.

## J08 decision and remaining acceptance

The merged source/session contracts did not specify source-free Save behavior.
The user answered the precise Product question in this session:

> Save Session before any source file is loaded should preserve a settings-only session

That choice is recorded separately. It is not an existing base-branch decision
or a completed implementation. Source-free persistence and the corresponding
contract/fixture work remain a separate follow-up, as required by the task.
J08 is not waived and its historical failed checkpoint remains in the ledger.
No rationale, retirement intent or accepted residual risk has been inferred.

Publishing the prepared work branches/PRs and maintainer review remain pending.
After integration, targeted and neighboring regressions and required dev/Gallery
gates must run on the exact resulting dev. These local results do not substitute
for those gates. The full 48-journey run, retired dense Q01 workloads, S11 and S12
were not started; no technical/publication baseline SHA is assigned. The next
acceptance phase is one strengthened retained run after integration, with S11
only if its acceptance conditions are satisfied.
