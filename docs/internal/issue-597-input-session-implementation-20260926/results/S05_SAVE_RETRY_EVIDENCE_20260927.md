# Check Session Save busy outcomes and retry after diagram updates

Local delivery: READY FOR REVIEW, uncommitted evidence-only candidate.
S05 overall: INCOMPLETE. The branch
`fix/issue-597-s05-reflow-retry-evidence-20260927` starts directly from actual
`origin/dev` `968d08211ce2315c274be14d879ef9ceb332bb66`, with no upstream, in
`/tmp/issue597-S05-resume2-20260927-t65_9osr/evidence-view`.
No commit, push, PR, merge or deployment has been performed.

## Changed evidence and preserved acceptance

One protected file changes:
`tests/web/contracts/session-regenerate-intent.playwright.spec.js`.

- `saveCurrentSession` uses the existing public Save action. Every busy attempt
  must return the exact Updating diagram reason, preserve Result-array and
  committed-request identities and small editor-override maps, leave Save
  pending false and agree with the unified availability outcome. It waits for
  that existing availability to settle, checks reflow has no error and
  explicitly invokes Save again. All attempts share the original single
  180000 ms download deadline; it is never restarted.
- The divergent case polls its original pre-Save capture until existing
  availability is clear and the captured SVG still equals the current mounted
  SVG after asynchronous capture. Trusted dev does not expose the optional S04
  availability accessor: the wait makes no S04 settlement claim there.
  The S05 integration view uses the existing accessor. The test attaches raw
  pre-Save, saved-file and fresh-Load checkpoints before its strict comparison.

All 142 original `expect` lines remain. Original test deadlines (420000,
600000, 600000, 300000 ms), fixtures, reference outputs, Worker criteria and
performance budgets are unchanged. Draft, canonical request, preview,
overrides, catalog, History, Save/Load/Generate and strict SVG assertions
remain. No normalization or expected-output update is introduced. The Feature
fill body is byte for byte unchanged. The patch and assertion audit are in
`evidence/S05-save-retry-evidence-20260927/`.

## Existing authority and integration order

PD-OI-045 and the Web live-edit/availability invariants already select busy
plus explicit retry and draft/Result preservation. No Product choice,
retirement, risk acceptance or Python allowance is inferred. Current code and
tests are evidence, not authority.

`WEB_CHANGE_POLICY.md` requires prior evidence-only integration before
replacing an affected runtime's mapped contract. This delivery has no runtime,
guard, checker, workflow, privileged-permission, authority or generated-asset
change. Exact `path::test-name` references and checkpoint/execution/sensitivity
records stay identical, so no reference-only authority amendment is proposed.
Any later actual reference change still requires that separate step. This
candidate does not admit dependent S05 runtime in the same delivery.

Feature fill's original checkpoint now passes on trusted dev following the
already merged PR #626 repair (`03621c0ecc1498738535b0ffc79d71ee160e6c4f`). The
inherited caption color remains in the existing canonical-rules/legend/History
transaction. Explicit Legend color/stroke edits retain independent coverage
in the unchanged direct-edit journey. See the separate preserved-checkout
`FEATURE_FILL_REVIEW.md` under `evidence/S05-resume2-20260927-t65_9osr/`.

## Actual verification and remaining failures

Final trusted-dev candidate run: five contracts PASS, exit 0 (Feature fill,
no-draft continuity, direct edits, divergent round trips, bare legacy config).
Exact command, body/source hashes, environment, port and artifacts are in
`dev-candidate-final-invocation.json`; the runner log is retained.

Supplemental real-reflow probe: PASS, exit 0. A held existing Worker `run`
exercises this same final helper's busy path. A second transient busy is also
checked before successful Save. The disposable timing probe is not shipped
as another contract and is not mapped acceptance.

The provisional S05 view is a no-commit merge of `ee6b6207…` with trusted dev,
not an approved final source. Historical v4 had a passing batch followed by
 two strict divergent FAILs. A safe root-attribute observer found availability
became busy again during the separate wait/capture calls; geometry changed
before Save entry, while mounted and Result geometry agreed and Save entry
and return geometry stayed equal. That timing probe is not acceptance.

The final settled-capture candidate passes all five contracts on trusted dev
and provisional S05 (exit 0 for each), plus two sequential exact divergent
repeats on S05 (exit 0 each). All strict comparisons and original round-trip
assertions remain. Historical root-width/legend/composition FAILs are retained
rather than normalized. Final invocation/log files are retained here; these
bounded runs do not claim exhaustive race coverage.

Real transport responsiveness, full Load/Save responsiveness, Python-preview
Worker count 0, real Generate, strict CLI SVG agreement and native clone-byte
measurement also remain unresolved in independent current-source evidence.
The complete measurements, recipes and preserved FAILs are in the owned raw
root `/tmp/issue597-S05-resume2-20260927-t65_9osr/` and the S05 resume report.

## Conditions, review and authorization

Node/Python Playwright 1.61.0, Chromium 149.0.7827.55, Node v26.8.2, Python
3.13.3, Linux/WSL2. Generated wheel SHA
`2c6e00b82c6493ebb631bb68a42c804467448dd73f27a27cc63f28132a63b8e9` is reused;
all 228 corresponding source files match. Original Desktop Chrome config and
deadlines remain. An external config supplies owned cwd/testDir and
`reuseExistingServer=false`; an external reporter only retains attachments.
Host scheduling is not certified idle. Runner startup/selection failures
are separately retained and do not count as test verdicts.

Production diff: empty. Test diff: one helper and one scenario checkpoint,
without another Save implementation. Docs/evidence diff: report, patch,
provenance, logs and reviews. Generated diff: empty; wheel is ignored.
This is not architecture-bearing; production OE/PE/CB sets are unchanged.
Local policy/Node checks are separate from browser acceptance. Candidate
trusted CI and human review are pending external publication authorization.
`candidate-manifest.json` binds the reviewed local files for handoff.

Authorize commit, same-named branch push and an evidence-only PR to `dev` for
this concrete candidate. Required trusted CI/review and prior integration
remain necessary before dependent S05 runtime work. No S06+ work was started.
