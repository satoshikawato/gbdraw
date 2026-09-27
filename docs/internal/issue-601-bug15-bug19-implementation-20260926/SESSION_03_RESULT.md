# S03 result: operation errors and regex dialects

S03 connects structured causes to Generate, Align, rule/table/Label imports,
Session notifications and SVG/PNG/PDF export callers. Known causes retain their
code, stage, bounded context and cleanup facts. Public wording remains owned by
`services/error-normalization.js`. S04 has not started.

## Checkout and authority

- Independent clone: `/tmp/gbdraw-issue601-s03-4r0rIu/repo`.
- Starting local and remote S02 head: `6724475af540ed26c0c6b7395c443b04e2cdb6bf`.
- Branch: `fix/issue-601-bug15-bug19`; upstream:
  `origin/fix/issue-601-bug15-bug19`.
- Initial fetched dev: `88028fd242d263f0fe86aaf9da57b8dc9eb082f6`.
  Normal synchronization merge: `9f483bf7cb9dbef4505564cd944d9930da7db4fc`.
  Conflicts in documentation retained the newest Align description and both
  independent README additions.
- Final fetched dev: `d313b70b9f97c2c1d70f9ae885edbead80b62021` (PR #620).
  Its only additional path was `tools/web-change-policy.json`, authorizing the
  separate #597 Session Worker. Normal synchronization merge:
  `70c4e0eac843720c7494b98dd1d824ff71c94b91`.
- Final S03 head is the commit containing this record; its exact object ID,
  trusted final-head checks and matching published remote ID are recorded in
  the session evidence `publication.json` and reported at handoff. Recording
  its own future object ID in this file would change that object ID.
- The shared workspace and earlier S02 clone were not edited. No reset, rebase,
  force push, main/dev push, PR, tag, release, deployment or Issue closure.

Formal authority `744be7a5943d4a247d027369629058898aa3f33e` remains an ancestor.
[PD-OI-046](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md#pd-oi-046-guidance-with-bounded-diagnostics)
implements concern `web.errors.diagnostic-disclosure`, Choice A,
`GUIDANCE_WITH_BOUNDED_DIAGNOSTICS`.
[PD-OI-047](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md#pd-oi-047-rejected-color-rule-pattern-edit-recovery)
retains concern `web.rules.rejected-pattern-edit-recovery`, Choice A,
`KEEP_REJECTED_PATTERN_DRAFT`, for S04. Both complete nine-field receipts,
including source, rationale, preservation, retirement, risk, owner and date,
remain byte-identical. The formal contract SHA-256 remains
`5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`.
No concern reassignment, export receipt retirement, new authority or invented
BD number is introduced. Preflight outcome: implement existing authority.

The #602 accepted runtime commit `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7`
remains a work-branch ancestor and is not an ancestor of final fetched dev.
That runtime receipt is distinct from independent dev integration. #602 compact
Editor/canonical draft and #598 Align direction/Reset contracts are retained.
The export handoff was read at `f0c8128ac9d252bbc1b64074cf995f536f789d03`; PDF/font/asset/
export-snapshot implementation remains with that owner. No export implementation
or History snapshot implementation is edited.

## Caller migration and return audit

| Caller/boundary | Retained cause and UI/recovery owner |
| --- | --- |
| Ordinary Generate; reflow | `run-analysis` normalizes the supplied cause, returns `error`, and retains helper and engine code/stage/context. Reflow uses its safe summary. |
| Committed candidate Generate/Align | Internal and outer catches return the same safe cause. Candidate removal destructures only `generatedArtifactCandidate`, preserving other outcome fields. Align explicitly supplies root operation `align`. |
| Align review, direction revalidation, Reset | Uses `outcome.error`; never borrows global `errorLog`. Canceled/stale/superseded outcomes remain separate. Missing diagnostics are explicitly UNKNOWN. Review keeps choices and Apply retry. |
| Record discovery, FASTA/CDS, comparison analysis | Producer models are thrown directly instead of reconstructed with `new Error(String(model))`. Discovery retains models with separate summary projection for the selector. |
| Manual/new rules, presets, color preparation | Existing preparation/transaction owners remain. Resource/preparation failures keep their stage; syntax is classified only by the actual producer. Retry captures the attempted operation; Edit focuses its real current input. Preset selection changes invalidate retry. |
| Auxiliary color/priority/whitelist/blacklist files | Serialized live application retains its file and selection owner. Failed imports restore selection and return false to the existing History checkpoint; rejected upload/retry adds no entry. Corrected import adds one. |
| Label TSV | Keeps its actual input across Retry and Reselect, plus current revision, document snapshot and intent checks. Retry uses existing Label/History ownership. |
| Session import/save | Safe model notification replaces raw exception alert/console. Existing awaited rollback remains. Late read/save failures cannot replace a newer notification. |
| SVG/PNG/PDF export callers | Retain model code/stage/context and actual export root. Remove reconstructed `type/message/details` and duplicate export prefixes. Retry/Use SVG call existing export entry points, including a separate retained Interactive SVG retry target. |
| Warmup, preview, definitions, legend, palette, replay | Automatic error logs use only the safe model. No source-dependent public wording owner is added. |

All five literal `status: 'error'` return sites in `run-analysis.js` contain
`error`: shared `failOperation`, reflow discovery, reflow depth validation,
reflow engine failure and reflow catch. Ordinary discovery/depth/engine/catch,
committed engine/catch and both public-wrapper catches delegate to the same
failure owner. There is no remaining status-only error return or caller cause
inferred from a previous global alert. The removed wrapper reconstruction and
raw import/preset/export prefixes are not retained as alternate paths.

Generation recovery is a transient transaction fact, separate from the error
model: `no-result`, `preserved`, `restored`, `restore-failed`. Preservation checks
compare the captured artifact owner's Result identities. Restoration claims
require completed awaited restoration. A rollback failure keeps the initiating
cause and prevents a success claim; History reports failed restoration without
replacing that cause with the cleanup exception. The artifact restore wrapper
retains notification ownership when snapshot restoration clones presentation.
A held old Worker failure and a held Session read cannot overwrite a newer
export notification. History observer revision may advance on a failed attempt;
its Undo/Redo stacks and retained artifact do not change. Generation IDs and artifact-owner identity also guard
superseded work. No `resultPreserved` property is invented in the public model.

## Public diagnostics and concrete actions

The S02 model stays
`{code, operation, stage, context, secondary, summary, details:[{label,text}], actions:[id]}`.
Only enumerated codes/operations/stages/context fields/reasons and bounded
numbers are projected. Limits remain 1,000 summary characters, 4,000 characters
per detail text, 8 detail sections, 2 cleanup facts and integers at most 10^7.
Python regex positions remain zero-based Python character positions; `😀[` is
reported at position 1, with `python-character`, rather than a UTF-16 offset.

The shared OperationError component displays a finite operation title and safe
summary, initially closed native Details, a readonly diagnostic textbox,
explicit Copy and Select controls, and a polite status message. Copy reads only
the currently displayed model text after Details is opened. No automatic copy.
An unavailable/rejected Clipboard retains the cause, Details, manual selection
and recovery controls. Keyboard opening, Tab/focus, role/name, manual selection
and reachable controls were checked in Circular/Linear at desktop and 390px.

Private exception/message/cause/traceback/stdout/stderr/notes/pattern/file/record/
path/SVG fields are absent from normalized summary, Details and copied text.
The native render and Align tests use private pattern/file markers; the catch
fixture also supplies private cause/stdout/cleanup properties. Logged output
and displayed/copy models contain no private marker. Re-normalization preserves
finite facts and does not duplicate a public prefix.

Buttons bind to actual existing owners: Generate retry, review Apply, rule retry
and focus, table retry, Label retry/reselect, Save Session, reload and specific
SVG/PNG/PDF export or Use SVG. No generic action dispatcher or diagnostic-derived
row/document selection is created. Caller retry context is transient; it does
not introduce a persisted draft manager.

Color and Label guidance explicitly says case-insensitive Python regex, with
`(?i)NADH` and `(?P<name>...)`. Feature Search and downloaded standalone controls
say `Regex (JavaScript, i)`, with fixed invalid-JavaScript-regex guidance and
turning Regex off for the existing word search. Evaluators and matching grammar
are unchanged. The live and actual downloaded SVG choose the same seven NADH
product features for JavaScript named-group regex and ordinary word search.
The original BUG-19 audit pattern/build remains unavailable; this is not a
claim to reproduce that unpublished artifact.

## Verification and environment

Evidence root: `/tmp/gbdraw-issue601-s03-4r0rIu/evidence` (not committed).
Dedicated venv: Python 3.13.3, pytest 9.0.2; dedicated Node dependencies:
Node 26.8.2, @playwright/test 1.61.1, Chromium 153.0.8010.12/revision 1243.
Python Playwright 1.61.0 uses Chromium 149.0.7827.55/revision 1228. Both paths
were verified. A dedicated editable install targets this clone; shared
installations were not changed. Commands use this clone's PYTHONPATH, dedicated
venv/Node bins and the session's PLAYWRIGHT_BROWSERS_PATH.

The browser wheel was prepared before browser checks and was not rebuilt during
them. Ignored wheel SHA-256:
`d7bfea8d67c45802688409ed2533805ea412e34f59a240c9fd4b119162c84774`.
No cache-bust refresh, wheel commit, reference update or social-preview edit.

| Local required check/evidence | Result |
| --- | --- |
| Core + standard recipes + Gallery Python, `python-core-recipes-gallery.log` | 6,630 passed, 17 skipped initially; three reported failures were fixed and the relevant files rerun: 33 passed (`python-recipes-recheck.log`). Existing unaffected evidence is reused. |
| Node contracts, `node-contracts-complete.log`, `node-required-final.log` | 987 passed including architecture/CI contracts; the exact required fast set rerun on final production passed all 847 tests. Final changed owners separately: 11 + 57 + 55 passing tests. |
| Python browser, `python-browser.log` | 37 passed initially; the remaining required comparison shard passed after correcting the server environment (`python-browser-thread-recheck.log`, 8 native Node comparisons through one passing Python test). |
| Final native/focused browser, `browser-final-changed.log` | 6 passed: actual downloaded Search, real Label Retry/Reselect, Python Color/Label/invalid-TSV parity and native Align failure/Apply retry. |
| Auxiliary file imports, `browser-import-search-fifth.log` | Five successful import Undo/Redo cases and failed TSV Retry/History preservation passed. Twelve malformed Session cases passed in the subsequent focused run and final smoke. |
| Align/export affected run, `browser-affected.log`, `browser-errors-final.log` | All 18 received Align UI contracts passed; all seven export lazy-loading/capture/failure/retry tests passed after migration assertions were corrected; final caller checks also verify Interactive SVG retry produces a standalone file. |
| Gallery first-Generate, `browser-gallery-final.log` | 9 passed using real Gallery sessions and canonical SVG semantic comparison. |
| PR smoke, `browser-pr-smoke-complete.log` | 19 passed; cases cover #602 UI, #598/native Align, Circular/Linear errors at both widths, Session/depth, Search, export, History and editor paths. |
| Final export caller checks, `browser-export-complete.log`, `browser-caller-final.log` | All 8 export tests passed, including actual Interactive SVG Retry/download. Final Session caller regression also passed. |
| Lint, `lint-final.log` | `ruff check gbdraw/`: passed. |

Reproduction commands, run from the independent clone with the above environment:

```bash
python -m pytest tests/ -m 'not slow and not browser' --durations=30
node --test tests/web/*.test.mjs
python -m pytest tests/ -m 'browser and not slow' --durations=30
python -m pytest tests/test_linear_comparison_browser_contracts.py -k '2/2' -v
playwright test --config=playwright.pr-smoke.config.js
playwright test --config=playwright.gallery-publication.config.js
playwright test tests/web/export-lazy-loading.playwright.spec.js --config=playwright.functional.config.js --workers=1 --retries=0
playwright test tests/web/error-boundary.playwright.spec.js tests/web/similarity-alignment-ui.playwright.spec.js tests/web/python-rule-parity.playwright.spec.js --config=playwright.functional.config.js --workers=1 --retries=0 --grep 'live and downloaded|real reselect|native Python engine|Python-only color|Python label TSV|invalid color TSV'
ruff check gbdraw/
git diff --check
```

Failed checks were diagnosed without relaxing retention/History or test-owned
timeouts: actual native test inputs needed a valid resource cache token; copied
module fixtures needed the normalizer dependency; rejected uploads needed the
existing checkpoint's `shouldCommit` false result; Label retry needed its real
input retained; old tests still awaited raw alerts or status-only export errors;
interactive recipe payloads needed regeneration after embedded wording changed.
The standalone search counter is active-index/total, its controls start collapsed,
and its existing Search button applies the pending query. A full-page SVG
screenshot caused an adaptive viewport loop after all matching assertions passed;
the evidence capture now uses a viewport screenshot. No evaluator change.

The threaded WASI child Worker test failed on the stock Python HTTP server,
whose nested Worker responses lacked COOP/COEP. The same unchanged acceptance
criteria passed with the existing benchmark server's isolation headers supplied
on every response. `isolated-http-server.py` records that disposable environment.
The normalizer still reports the original unrecognized native failure as UNKNOWN;
it was not changed to make this test pass. No test execution was abandoned at
an early command timeout. The export-owned standalone dynamic import caches its
rejected promise; an
injected module-fetch failure still needs reload. It is outside the handed-off
S03 caller scope; no export loader fallback or cache policy is changed. Successful
Interactive SVG retry is verified after correcting its actual catalog input.
No remote CI or Python 3.11/Node 20 cross-version success is claimed by these local results.

Existing public artifacts H-CLI-13 and T-PY-08 were regenerated by their real
recipes, with only the three embedded Search wording/title lines changed:

```bash
python docs/recipes/run_cli_scenarios.py --scenario H-CLI-13
python docs/recipes/run_python_scenarios.py --scenario T-PY-08
```

SVG payload comparison passes. The two browser renders at 1400×1050 were visually
inspected: biological labels, legend, record metadata and existing quantitative
tracks are retained. Re-generated PDF/EPS/PS timestamp differences were restored;
no geometry/reference or export-owner changes were included. Final 390px error
crops were also inspected. All screenshots stay in the evidence root.

## Gate, review and S04 handoff

Trusted tools are archived from the fetched dev commit, rather than taken as
candidate authority. Initial base `88028fd242d263f0fe86aaf9da57b8dc9eb082f6`
and synchronization head `9f483bf7cb9dbef4505564cd944d9930da7db4fc` classified
`profile=pr / impact=full`, requiring `web-change-budget`, `core-pr`,
`recipes-standard`, `gallery`, `lint`, `web-contracts-pr`, `web-pr-smoke`.
Final base is `d313b70b9f97c2c1d70f9ae885edbead80b62021`; exact final S03 head
is checked after recording, with output retained as `policy-final-head.log` and
`ci-plan-final.log`. Working-tree Gate is PASS; Review is REQUIRED. These are
separate outcomes. A newly added privileged import was detected during development
and removed; no allowlist, checker, workflow or authority was changed to admit it.

Architecture review: one public wording/model owner before and after; existing
Worker producers, canonical renderer entry, artifact transaction and History
owners remain. Remaining reconstruction paths converge to the normalizer;
OperationError renders it, and concrete callbacks stay with current callers.
S03 adds no production module, evaluator, dispatcher, persistence format,
compatibility path or authority. The cumulative module addition is S02's existing
normalizer. Registered semantic-owner/canonical-entry results conform; dependency
cycles remain zero and privileged fan-out/permissions remain unchanged. Ordinary
non-increasing owner/path evidence applies; no architecture exception is sought.
Size/inventory/session-path review reasons remain visible in the trusted report.
Production, tests, documentation and generated diffs were reviewed separately.

S04 can consume the same finite field model, Python position units and diagnostic
component. It must use actual caller row/document ownership and preserve revision,
Result, History and Session semantics independently of the error model. Rejected
Color draft lifetime, Not applied, field Retry/Revert and their accepted PD-OI-047
outcome remain S04 work. This session retains the existing rejected-input behavior
and does not declare that lifecycle complete. Inherited metadata-free Session and
Legend override limitations remain. Export shared-file transfer still waits for
S05. S03 does not mean all BUG-15/BUG-19 implementation phases are complete.

Commit title: **Show actionable operation errors and clarify regex dialects**.
Summary: **Retain caller causes, safe diagnostics, recovery actions, and existing search semantics.**
