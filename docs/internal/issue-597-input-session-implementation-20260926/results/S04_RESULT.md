# S04 result — Session operation consistency

Status: **S04 complete**, with local validation below and required human Review / S05/S06 boundaries preserved. Scope: Issue #597 **BUG-02 / BUG-20, S04 only**; acceptance **S-02**, architecture **A-01**, delivery **W-01**. BUG-01/S02 were not reintroduced. This is not completion of Issue #597.

## Source and authority

Dedicated clone: `/tmp/gbdraw-issue597-S04.3Pw7Yo`; branch `fix/issue-597-input-session-20260926`, upstream the same-named remote branch. Shared checkout/environment/server were not changed. Starting remote/S03 HEAD: `e3684d496df7a3b2a7492443cad50168f6959b36`.

Fetched trusted `origin/dev`: `c922fc38aac78da9be83342c09ac0164ecef6ff6`. Clean dev integration, separate from the S04 implementation commit: `eed70f05d04f638a5a4ab6c6567982d45e774596`. S01 `15bcbca89392cbd63a88fe80dc44b71ad4868061`, S03 integration `ff02abbae948aaa870f18123805d7297a6c3495d`, S03 implementation and authority merge `af5d942af60353dda199aa487da9152a3576b3fe` are ancestors of the measurement HEAD. No merge conflict or authority/runtime self-authorization change.

The original `python tools/inspect_issue597_s01_contracts.py` was run and remains **FAIL** at its hard-coded revision-21 assertion. Latest trusted dev already has Contract revision 22. The script, checker and authority were not modified to make this pass. [S04_AUTHORITY.py](../evidence/S04_AUTHORITY.py) independently reproduces strict historical/current evidence: revision 21 at the actual authority merge; current revision 22 exact bytes equal fetched dev; PD-OI-044/045 nine fields equal their human and serialized receipts in both revisions; merge ancestry holds. Its result is **PASS**. This supplemental receipt proof is not a policy checker or new authority registry.

[Validation manifest](../evidence/S04-validation.json) records source SHA-256 values, runtime aggregate, actual measurement HEAD, fixtures/boundary files and environment. Measurement HEAD is the integration commit plus the recorded working source; it is not falsely labeled as the final implementation commit. Python 3.13.3, Node v26.8.2, Node/Python Playwright 1.61.0, Chromium 149.0.7827.55, WSL2 Linux. Node dependencies and browser wheel were prepared only in this clone; package/lock/runtime vendor artifacts were not changed.

## Implementation and ownership

[Mutation/checkpoint evidence](../evidence/S04_MUTATION_INVENTORY.md) and [action/native binding inventory](../evidence/S04-mutation-inventory.json) enumerate DOM, programmatic, asynchronous completion, owning state and availability checkpoints. Read-only viewport pan/zoom, scroll, search and popup/drawer navigation remain available. Persisted Result selection and diagram/legend/title/length-bar placement use semantic admission.

State derives one availability from existing Save/import flags. Generate/reflow and actual History/edit/import/definition resources supply preparation reasons. Native UI and existing action owners use the same predicate. Busy precedes validation, so invalid arguments cannot bypass it; the response contains a reason and retry continuation. No second lock ref, independent document snapshot owner or generic lock/queue framework.

`config.exportSession` owns the existing singleflight and exact duplicate promise join, before generic busy. The app's former singleflight was removed; app adapters handle title, paint and private record catalog preparation. Small mutable config/navigation state is copied immediately after title/pending publication. Adopted request/resources/catalog/cache/Result payloads remain under their existing owners; no new full graph hash/sign/base64 conversion. The private uncataloged Linear Save preserves one `cardinality: all` source request and exact two-record source bytes, without expanding live source rows during Save. Settings-only auxiliary rows retain their prior behavior.

Load reconstructs files into a private target through the existing reconstruction owner. Validation, normalization, resource/sequence/catalog/recovery preparation and SVG admission finish before reset/adoption. The old live legacy-recovery helper/path was removed. The existing transaction adopts main data and registers preview readiness before its first yield. Failure restores source/config/request/resources/Result/cache/evidence first, then reconciles transient UI. Immutable old references are reused. History baseline runs inside the rollback scope and checks lifetime after capture, before clearing stacks.

Generate's existing processing owner now includes source preflight and final UI reconciliation; automatic reflow's existing flag spans all iterations. Delayed discovery/file/rule/definition/style/legend completions recheck identity and availability. Discovery deferred during Session retries after settlement. Teardown cancels old operation handles and clears pending; stale Save cannot download or clear a newer Load.

Ordinary architecture evidence: one state availability owner and one config operation owner; existing canonical request/Result/resource/admission owners retained. App singleflight and live preparation paths removed in this change. No new compatibility path, parser/schema/renderer, Worker constructor/importer, active registry/checker or public default. Derived UI values and actual promise/timer handles do not create a second document lock. No independent ratchet exception was invoked; complete exception OE/PE/CB sets and waivers are not claimed.

## Verification and artifacts

All task artifacts are outside the repository under `/tmp/issue597-S04-evidence-3Pw7Yo`; tracked manifest holds checksums and regeneration commands. Diagnostic iterations are not acceptance PASS evidence. No sequence/full biological rows were intentionally logged by the added tests or evidence probe.

| Command / target | Result and log |
| --- | --- |
| `python tools/prepare_browser_wheel.py --no-build-isolation` | PASS; `wheel.log`; generated wheel ignored |
| `python tools/inspect_issue597_s01_contracts.py` | FAIL, obsolete revision assertion; `authority-historical-probe.log`; unchanged probe |
| `python docs/internal/issue-597-input-session-implementation-20260926/evidence/S04_AUTHORITY.py` | PASS; `authority-final.json` |
| `node --test tests/web/*.test.mjs` | 723 PASS, 0 failures, 78.923 s; `node-final-current.log` |
| Final combined functional browser command below | 39 PASS, 2.6 min; `browser-final-current.log`, `browser-final/` |
| Existing Vibrio performance command below | 1 PASS, 1.4 min; `browser-vibrio-final.log`, `browser-vibrio-final/` |
| `GBDRAW_WEB_TEST_PORT=43015 pytest tests/ -v -m 'not slow'` | 6262 PASS, 6 FAIL, 17 skipped, 11 deselected, 1500.33 s; `pytest.log`; all six failures corrected and rerun below; do not call this invocation PASS |
| `GBDRAW_WEB_TEST_PORT=43023 pytest tests/test_linear_comparison_browser_contracts.py -v` | 2 PASS, 213.91 s; `pytest-comparison-recheck.log` |
| Focused Python command below | 10 PASS, 0.94 s; `pytest-focused-final.log`; covers remaining four failing cases |
| Final changed History/preview/settings/mode/active-files/operation Node targets | 25 PASS; `node-final-changed.log` |
| Final style/palette/rule Node targets | 8 PASS; `node-style-final.log` |
| `ruff check gbdraw/` | PASS; `ruff-final.log` |
| `node tools/check-web-change-budget.mjs --base origin/dev` | Gate PASS / Review REQUIRED; `policy-final-before.log`; unchanged trusted-base checker |
| `node tools/check-web-change-budget.mjs --base origin/dev --head HEAD` | Required after implementation commit; actual result recorded at handoff in `policy-after-commit.log` |

Final single-observation Vibrio metrics: Save 5540 ms; heap delta 177459158 bytes; maximum 100 ms-heartbeat gap 962.6 ms; compressed 23738963 bytes. The existing test includes content comparison, fresh Load and CLI/Python replay. Passing its historical budgets does not meet the whole Issue #597 full-pipeline/100 ms responsiveness target.

Initial full-Python failures were two browser contract timeout/early-candidate failures, one fixture-builder executable-permission failure, and three source assertions. Both browser shards pass on corrected source with the existing 300 s test timeout unchanged. The fixture failed because its LOSAT executable mode was restored while the full run still needed it; the same builder passes with executable permission and the mode is restored only after child tests finish. Source assertions now cover dynamic uploader keyboard availability and suppression of the watcher that would invalidate already-imported mode profiles. Stable selectors remain intact. These are ordinary test adaptations, not active detector/authority or acceptance-budget changes. The additional History browser test's hard-coded 4173 origin was replaced with the existing configured port, preserving exact local-origin filtering and allowing this independent checkout's server.

Production, tests, plan evidence/docs and generated diffs are reviewed separately. No tracked reference SVG, Gallery/tutorial/public image, owner social preview, binary content, dependency manifest, dist/egg-info, cache-bust or vendor changes. Test-produced LOSAT executable mode is restored before staging. Stage uses explicit assigned paths only. The local Review requirement remains; it is not a waiver, trusted CI result or merge approval.

## Reproduction

From a dedicated clone at the delivered branch HEAD, prepare the existing wheel and local Playwright environment, and use separate ports. Preserve tests' timeout and performance budgets.

```bash
python tools/prepare_browser_wheel.py --no-build-isolation
python docs/internal/issue-597-input-session-implementation-20260926/evidence/S04_AUTHORITY.py
node --test tests/web/*.test.mjs
GBDRAW_WEB_TEST_PORT=43024 node node_modules/@playwright/test/cli.js test \
  --config=playwright.functional.config.js --workers=1 \
  tests/web/circular-record-presentation.playwright.spec.js \
  tests/web/record-display-discovery.playwright.spec.js \
  tests/web/settings-only-session.playwright.spec.js \
  tests/web/session-operation-consistency.playwright.spec.js \
  tests/web/session-save-lifecycle.playwright.spec.js \
  tests/web/session-loading-feedback.playwright.spec.js \
  tests/web/history-inputs.playwright.spec.js \
  tests/web/auxiliary-file-history.playwright.spec.js \
  tests/web/palette-history-preservation.playwright.spec.js
GBDRAW_WEB_TEST_PORT=43025 node node_modules/@playwright/test/cli.js test \
  --config=playwright.perf.config.js --workers=1 \
  tests/web/vibrio-session-save.performance.playwright.spec.js
GBDRAW_WEB_TEST_PORT=43023 pytest tests/test_linear_comparison_browser_contracts.py -v
pytest tests/test_tutorial_fixture_manifest.py::test_metazoan_mtdna_comparison_builder_is_byte_reproducible \
  tests/test_web_mode_profiles.py tests/test_web_ux_profile.py -v
GBDRAW_WEB_TEST_PORT=43015 pytest tests/ -v -m 'not slow'
ruff check gbdraw/
node tools/check-web-change-budget.mjs --base origin/dev
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
```

The authority reproduction pins the S04 fetched base and historical authority merge; future dev changes require their own current authority comparison rather than weakening this recorded assertion. The runtime aggregate in the manifest hashes sorted paths and SHA-256 values for tracked Python sources, `gbdraw/web/js/` and `gbdraw/web/index.html` using sorted compact JSON.

## Remaining boundary and next entry

S01's selected import transport remains **whole-object reply**. S04 adds no import Worker. Worker constructor/importer permission is a separate authority dependency awaiting dev integration; candidate permission does not authorize runtime. Full Issue #597 large Save/Load performance remains unresolved from S01 (**FAIL**, not replaced by a small Circular success or the narrower existing Vibrio test). Real large Linear preview Load's existing Python/helper Worker startup was not eliminated. S05 owns selected import transport/lifecycle; S06 owns measured full-pipeline stage/wall/100 ms heartbeat/long-task/main+Worker/process-heap/copy/content gates. Full new large-dataset S-03 matrix, repeated statistics and offline bundle audit are not claimed.

Next: **S05 only after required privileged authority has merged into dev**. Fetch/clone the latest same-named branch in a new dedicated checkout, read this result and manifest, [S01_RESULT.md](S01_RESULT.md), [S03_RESULT.md](S03_RESULT.md), [03_SESSION_OPERATIONS.md](../decisions/03_SESSION_OPERATIONS.md), [SESSION_WORKFLOW.md](../SESSION_WORKFLOW.md) and [S05_INSTRUCTION_PROMPT.md](../sessions/S05_INSTRUCTION_PROMPT.md). Verify actual trusted dev ancestry/receipts and exact Worker subjects; merge the required authority into a clean tree according to the workflow. If permission is absent, report that specific dependency and do not introduce the Worker. Keep one selected whole-object transport, remove main parse only within S05, preserve S04 availability/candidate/adoption/rollback checkpoints, and rerun S-02 plus S05 failure/stale/limits/transfer evidence. Do not start S06/S07/S08 in this session.

English commit title: **Keep Save and Load consistent across browser mutations**

English summary: Use shared Session availability across browser mutations, keep Save on one document, and prepare Load privately with source-first rollback while preserving browsing, History and supported Sessions.
