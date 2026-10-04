# S02 — Structured failure causes across Python and the diagram Worker

Date: 2026-09-27 (JST). Scope: Issue #601 BUG-15/BUG-19, **S02 only**.
This is an implementation and handoff record, not Product authority. S03–S05,
Details/Copy wiring, Generate/Align caller integration, regex dialect labels,
and rejected Color field drafts are not implemented here.

## Source, authority and history

| Boundary | Exact source |
| --- | --- |
| Dedicated clone / environment / evidence | `/tmp/gbdraw-issue601-s02-R33zEf/{repo,venv,nodeenv,evidence}`; independent clone, no shared objects, checkout, index or editable installation |
| Starting local and remote branch | `75bee53d2c98531be4a4d3184c11246b45a459dd` |
| Latest fetched `origin/dev` | `f5f86634459e0dcd46c1a452e9219fbba635d429` |
| Branch and upstream | `fix/issue-601-bug15-bug19` / `origin/fix/issue-601-bug15-bug19` |
| Separate ordinary dev synchronization merge | `6e9ef19706a9f0a44a41ccbdc181a368f7b47e78` |
| Formal authority integration | `744be7a5943d4a247d027369629058898aa3f33e` (PR #615), ancestor of both dev and this branch |
| Received #602 S07 | `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7`, ancestor of this branch; **not** an ancestor of dev |
| Published export owner transfer | `f0c8128ac9d252bbc1b64074cf995f536f789d03`; export `OWNER_HANDOFF_20260927.md` inspected from its remote ref |
| S02 commit | The commit containing this document: `git log -1 --format=%H -- docs/internal/issue-601-bug15-bug19-implementation-20260926/SESSION_02_RESULT.md` |

The latest start condition is [S01 handoff](SESSION_01_HANDOFF_RESULT.md):
`OWNERSHIP_TRANSFERRED / RUNTIME_RECEIVED / S02_READY`. Historical S01 blocked
records remain unchanged. #602 [S07](../issue-602-proposal-20260926/results/S07.md)
and #598 [dev integration](../issue-598-alignment-direction-reset-20260926/evidence/DEV_INTEGRATION.md),
[follow-up](../issue-598-alignment-direction-reset-20260926/evidence/DEV_INTEGRATION_FOLLOWUP.md)
and [S03 integration](../issue-598-alignment-direction-reset-20260926/evidence/DEV_S03_INTEGRATION.md)
were inspected. Receiving #602 in this branch does not publish its independent
dev integration. The sole merge conflict was the end of `session-request.test.mjs`;
independent #602 projection and dev pixel-gap cases were both retained and run.
No reset, rebase, force push or replacement with the dev tree was used.

Developer preflight: **IMPLEMENT_EXISTING_AUTHORITY**. OIPC revision 22,
PD-OI-046 (`GUIDANCE_WITH_BOUNDED_DIAGNOSTICS`) and PD-OI-047
(`KEEP_REJECTED_PATTERN_DRAFT`), scenario revision 1, already exist on dev.
Formal contract SHA-256 remains
`5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`.
The original concern keys and all receipt fields remain byte-identical. No
second active diagnostic authority, BD record, export receipt retirement,
checker/rule change or CI guard modification is introduced. The registered
Product map/store are unchanged; their pilot coverage does not by itself prove
these two independent outcomes complete.

## Finite diagnostic contract

Python `web_support/error_adapter.py` serializes native cause facts.
`services/error-normalization.js` exclusively owns public wording, permitted
context, bounded sections and action IDs. Worker and client delegate to that
same normalizer. Successful request/Session and outer Worker message schemas
remain unchanged; there is no parallel raw-error Worker protocol.

The producer envelope has `code`, `operation`, `stage`, `context`, and optional
`secondary`. The public model adds `summary`, `details: [{label,text}]` and
`actions: [id]`; it always includes a `secondary` list. Client errors expose the
same model with safe `message` and `stack`; raw exception names, traces, causes,
logs, notes and preformatted details are never copied. Re-normalization rebuilds
the model from its identities and context, without prefixes or cause loss.

- Codes: `UNKNOWN`, `VALIDATION_UNCLASSIFIED`, `INPUT_INVALID`, `INPUT_REQUIRED`,
  `FASTA_REQUIRED`, `INPUT_UNREADABLE`, `NO_RECORDS`, `RECORD_SELECTION`,
  `REGION_INVALID`, `DEPTH_INVALID`, `TABLE_INVALID`, `COMPARISON_INPUT`,
  `COMPARISON_IDENTITY`, `ANNOTATION_TARGET`, `TRACK_INVALID`, `TRACK_LAYOUT`,
  `REGEX_SYNTAX`, `RESOURCE_INVALID`, `HELPER_PROTOCOL`, `RUNTIME_INCOMPATIBLE`,
  `WORKER_INIT`, `FEATURE_METADATA`, `RESULT_INVALID`, `CLEANUP_FAILED`,
  `EXPORT_INPUT`, `EXPORT_DIMENSIONS`, `PNG_DPI`, `EXPORT_CONVERSION`,
  `PDF_LIBRARY`, `PDF_GLYPH`.
- Operations: `unknown`, `generate`, `align`, `feature-extraction`,
  `export-svg`, `export-png`, `export-pdf`, and the existing 17
  `DIAGRAM_HELPER_OPERATION_NAMES`. Native helper wrappers use the corresponding
  names, including `evaluateRules`, `resolveSimilarityAlignment` and `readPdfFont`.
  Domain parity is tested across Python, JS and the existing Worker protocol.
- Stages: `unknown`, `initialization`, `resource-staging`, `request-validation`,
  `helper`, `rule-validation`, `render`, `result-admission`, `cleanup`,
  `export-capture`, `export-conversion`, `font-validation`.
- Context strings: only registry `field` IDs and finite `reason` IDs;
  `positionUnit` is exclusively `python-character`. Arbitrary keys and values
  are dropped. Python identifiers are capped at 80 characters, and JS admits
  only members of its finite domains. The field domains are cross-tested.
- Context numbers: `position`, `row`, `inputOrdinal`, `recordIndex`,
  `seriesIndex`, `slotIndex`, `recordCount`, `columnCount`, `codepoint` must be
  nonnegative integers, at most 10,000,000; codepoint is at most U+10FFFF.
  Position is admitted only with the Python character unit. Position/indexes
  are zero-based; native table rows and input ordinals retain their native
  one-based units. Missing positions stay absent, never zero by assumption.
- Regex reason IDs: `UNTERMINATED_SET`, `UNTERMINATED_GROUP`, `UNBALANCED_GROUP`,
  `NOTHING_TO_REPEAT`, `MULTIPLE_REPEAT`, `FLAGS_POSITION`, `LOOKBEHIND_WIDTH`,
  `INVALID_ESCAPE`, `UNKNOWN_EXTENSION`, `CHARACTER_RANGE`, `GROUP_REFERENCE`,
  `GROUP_NAME`, `SYNTAX_ERROR`. Only an actual `re.error` or explicit `__cause__`
  supplies regex identity; raw JS/Python error text is not a regex classifier.
  Cause traversal is bounded to eight nodes and cycle-safe.
- Comparison reason IDs: `EMPTY_ENDPOINT`, `INDEX_ALIGNMENT`, `SOURCE_INDEX`,
  `SOURCE_VIEW_CONFLICT`. Other correction reasons are the finite `REASONS`
  registry in the wording owner: numeric/shape constraints, record/region
  selection, supported option domains, table fields and typed track issues.
- Cleanup: at most two `{code:'CLEANUP_FAILED', stage:'cleanup'}` secondary
  entries. Native Python notes stay available to native callers, while the Web
  envelope drops their free text. Render/helper cleanup never replaces a
  primary error. Cleanup-only failures have their own code and actual stage.
- Limits: summary ≤1,000 JS characters, each detail text ≤4,000, sections ≤8.
  Currently one fixed `Diagnostics` section is emitted. Caller limit overrides
  can reduce these bounds but cannot increase them.

Original patterns, sequences, file/record names, paths, SVG, free exception
text, stdout/stderr, notes and arbitrary detail labels/text are excluded.
`private_web_execution` discards native stdout/stderr without retaining them
and restores the prior native logging-disable level after each Web invocation.
Native validators, exceptions and explicit causes outside that scope remain
available to CLI/API callers. No regex compilation, input repair or second
validator is added to the diagnostic adapter.

## Known failure migration

The baseline [S00 inventory](SESSION_00_RESULT.md#known-failure移行表) is the source
of this table. Exact/anchored native validation clauses are projected into
finite correction facts; input values are not revalidated by the adapters.

| Native failure family | Finite destination and retained correction | Existing recovery action IDs |
| --- | --- | --- |
| GenBank/GFF3/FASTA required | `INPUT_REQUIRED` / `FASTA_REQUIRED`; input combination, actual input ordinal where supplied | `select-input`, `retry` |
| Input unreadable / record discovery / zero records | `INPUT_UNREADABLE` / `NO_RECORDS`; reselect readable input or input containing records | `select-input`, `retry` |
| Record selector range/absence/duplicates | `RECORD_SELECTION`; loaded count, `OUT_OF_RANGE`, `NO_MATCH`, `AMBIGUOUS`, `SELECTOR_FORMAT`; retain #index correction without selector text | `select-record`, `retry` |
| Region/crop/Align region facts | `REGION_INVALID` / `INPUT_INVALID`; both endpoints, integer/positive/order/bounds/crop-start constraints, fixed region/start/end fields | `edit-region` or `edit-input`, `retry` |
| Annotation target/track set/anchor | `ANNOTATION_TARGET` or native typed track issue → `TRACK_INVALID`; choose target/set, clear target or enable multi-record canvas, anchor/side/layer/order corrections; slot index without names | `edit-annotation` or `edit-track`, `retry` |
| Comparison BLAST/FASTA/no proteins | `COMPARISON_INPUT`; required outfmt 6/7/FASTA, empty FASTA ordinal, CDS protein requirement and identifier consistency | `edit-comparison`, `retry` |
| Pairwise Height and protein/collinear options | `INPUT_INVALID`; finite field, Auto/positive/integer/nonnegative/boolean or exact supported option choices; Pairwise max hits remains readable | `edit-input`, `retry` |
| Depth source/series/track mismatch | `DEPTH_INVALID` or typed `TRACK_INVALID`; source requirement, series/index where available, existing source-or-disable/remove guidance | `edit-depth`, `disable-track`, `retry`, or `edit-track` |
| Depth/GC min/max/window/step/ticks/font | `INPUT_INVALID`; fixed field, finite/positive/order correction | `edit-input`, `retry` |
| Circular track cannot fit inside | `TRACK_LAYOUT`; `CANNOT_FIT`: move track, reduce widths, disable conflicting labels or place outside; custom names removed | `edit-track`, `retry` |
| Python color/label/whitelist/visibility regex | `REGEX_SYNTAX`; reason, actual Python character position, native table row where present | `edit-pattern`, `retry` |
| Table columns/required values/color/action | `TABLE_INVALID`; row, column count, fixed field, supported colors or visibility actions; never classify these as regex syntax | `edit-table`, `retry` |
| Source/view comparison endpoint disagreement | Dedicated `ComparisonIdentityError(ValueError, GbdrawError)` → `COMPARISON_IDENTITY` with one of four reasons; same native message/capture behavior; no identity repair | `edit-comparison`, `retry`, `save-session` |
| Helper protocol/JSON/asset/runtime/resource | `HELPER_PROTOCOL`, `WORKER_INIT`, `RUNTIME_INCOMPATIBLE`, `RESOURCE_INVALID`; actual failure stage and safe correction; compatibility class capture retained | `retry`, `save-session`, `reload`, `select-input` as applicable |
| Result/feature metadata admission | `RESULT_INVALID` / `FEATURE_METADATA`; actual admission stage and reload/retry guidance | `retry`, `save-session` or `reload` |
| Existing SVG/PNG/PDF preparation, DPI, conversion, library, glyph | Finite export codes preserve direct normalizer guidance and bounded codepoint; PDF/font/snapshot producers are not edited; existing export caller still needs S03 integration | `generate`, `edit-dpi`, `retry`, `reload`, `use-svg` as applicable |
| Recognized validation class with unmapped clause, or unknown typed track issue | `VALIDATION_UNCLASSIFIED` rather than `UNKNOWN`; explicit classification omission, review inputs before retry | `review-input`, `retry` |
| Other unknown native cause | `UNKNOWN`; actual operation/stage when observed, otherwise `unknown`; no fabricated reason/position | `retry`, `save-session` |
| Cleanup secondary or cleanup-only | bounded `CLEANUP_FAILED` identities; primary cause wins | primary actions, or `save-session`, `reload`, `retry` |

Known correction fixtures exercise real native record/depth/table/visibility,
Collinear options, Align region and Web current-option validators. A missing
classification therefore fails the fixture instead of being hidden by a
universal UNKNOWN assertion. Malformed private context and newly unmapped
validation remain explicitly observable. This finite inventory does not claim
to explain every third-party runtime exception.

## Ownership, removed paths and architecture evidence

- Producer identity remains in `pairwise_match.py`. The five existing rejection
  sites use the dedicated ValueError-compatible type; geometry and successes
  are unchanged. Native validators keep their existing owner and causes.
- One Python serializer replaces render traceback dictionaries, helper
  `traceback.format_exc()` payloads and feature-metadata `str(exc)` payloads.
  `call_web_json_helper` wraps the existing helper function references before
  Pyodide can stringify an exception. It shares the same lazy runtime/cache.
- `request_render.py` annotates observed decode/render/admission/resource
  boundaries even when performance diagnostics are off. Native cleanup notes
  remain native; the Web envelope only receives finite secondary facts.
- Worker `serializeError` and client `deserializeWorkerError` delegate public
  semantics to one normalizer. Old exception-name/message/stack/notes transport,
  raw cleanup diagnostic concatenation and client fallback-message strings are
  removed. Helper rejection propagates the finite cause; success objects,
  prepared resource ownership, transfers, cancellation and stale guards remain.
- The wording owner removes free traceback/stdout/stderr/detail extraction,
  `safeErrorText`, and exception prefixes. Existing typed track issues are
  projected by their finite native codes, without copying private messages.
- The only `run-analysis.js` production change removes its superseded Circular
  raw traceback extractor. This is required removal of a migrated classifier,
  not S03 caller/UI implementation. Its successful/transaction paths are unchanged.

Changed capabilities are producer failure identity, native diagnostic projection,
public error wording, error transport, and migrated raw extraction. Each has
one determining owner and one canonical path after convergence. Existing Python
and JS validators own input semantics; adapters read failure facts, not inputs.
The old paths above are removed in this change. No privileged importer, reactive
state authority, compatibility reader, alternate Worker/runtime, schema, input
normalization or lifecycle is added. OE/PE do not increase; CB is unchanged.
No architecture-exception condition applies. Complete repository-wide totals
are not inferred from the passing detector. Ordinary rollback is a revert of
this S02 commit, leaving the preceding normal dev merge and received #602 history.

Registered Product requirements for current-result admission and canonical
render request remain independently satisfied. #602 Status, canonical projection,
Result/draft/History/Session and compact Editor/review owners are unchanged.
#598 Align/Reset/direction/retry owners remain intact. PD-OI-046 cause/guidance,
privacy, manual diagnostics/recovery and initial/previous-result conditions,
and PD-OI-047 regex semantics, atomic rule commit and field-draft lifecycle are
independent contributions: only their S02 boundary contributions are implemented.
Passing transport tests does not prove the remaining UI/draft contributions.

## Verification and trusted admission

Environment: Python 3.13.3, Node 26.8.2; Python Playwright 1.61.0 and dedicated
Node Playwright 1.61.1 confirmed. Venv uses existing system packages without
changing their editable install; tests use this clone's `PYTHONPATH`. Chromium
checks use the necessary local sandbox escalation because this host's sandbox
cannot mount `/mnt/wslg/distro`. No automatic approval rejection occurred.
All logs, ports, wheel/browser outputs and temporary recipe directories are
session-specific. Test-owned timeouts and acceptance conditions are unchanged.

| Check | Result / dedicated log |
| --- | --- |
| Complete non-slow Python: core, recipes, Gallery, browser and reference comparisons | 6,671 passed, 17 skipped, zero failed; 597.52s; `evidence/python-full-final.log` |
| Affected producer/native comparison set | 274 passed; `evidence/python-focused-final.log`; subsequent native Align-region tests included in full set |
| Required fast Web contracts | 846 passed, zero failed/skipped/canceled; `evidence/web-contracts-node-accepted.log` |
| Final finite contract additions and Worker/client consumers | 13 passed; `evidence/boundary-contract-accepted.log`; `normalizer-consumers-final.log` covers final public consumers |
| Architecture and CI contracts | 204 passed (139 architecture plus 65 CI), zero failed/skipped/canceled; `evidence/architecture-ci-final.log`; final architecture rerun 139 passed in `evidence/architecture-final.log` |
| Real Pyodide boundary and Python rule parity | 9 passed; `evidence/browser-boundary-accepted.log` |
| PR smoke on stable final wheel | 12 passed on stable wheel plus corrected Annotation case 1 passed; `evidence/pr-smoke-fixed.log`, `evidence/pr-smoke-annotation-final.log` (all 13 covered; no single-run whole-suite pass claimed) |
| Gallery first-Generate parity | 9 passed; `evidence/gallery-parity.log` |
| Earlier failed comparison shard, unchanged oracle rerun | 1 passed (eight browser contracts); `evidence/comparison-rerun.log`; both shards rerun in final Python set |
| Lint, whitespace, reference/social-preview immutability | pass; `evidence/ruff-final.log` and Git diffs |
| Normal generated browser wheel | pass; `evidence/wheel-final.log`; ignored asset, no cache-bust change |
| Trusted-base Web policy | Gate PASS / Review REQUIRED; `evidence/policy-final-working.log`, then exact committed base/head checks before push |

The real boundary tests run embedded Python helpers and a staged canonical
render, then compare code/operation/stage/context through Worker serialization,
client deserialization and re-normalization. Real browser tests use Python
syntax with empty/unrelated catalogs, Unicode position 1, failing then corrected
render/helper operations, private sentinel rejection from browser console/model,
one lazy Worker, Python rule targets/History/Session/stale isolation, responsiveness
and prepared reuse. Worker tests retain successful binary transfers/artifact
identity/cache ownership and cancel/stale/superseded recovery assertions.

Diagnostic failures found during validation were fixed without silent fallback:
Collinear search scope and Pairwise max hits needed finite corrections;
Annotation needed projection of its existing typed native issues. Old assertions
for arbitrary injected messages or private annotation names now assert codes,
correction facts, privacy and the same state-retention conditions. A fixture
failed initially because the clone's LOSAT binary was not executable; the
existing local test setup supplied its executable bit and reproduction passed.
The final handoff restores only that mode, with unchanged binary bytes.
The last Annotation UI assertion was corrected after a full smoke runner had
already loaded the spec; its separate corrected case passes with every original
state/History/preview assertion retained. A relative `PYTHON` path failed from a disposable recipe directory; final Node
execution uses an absolute path. One smoke attempt fetched the wheel while it
was being rebuilt and received a traced 404; the final run uses a fixed wheel.
Initial stale-wheel evidence and missing Node dependency environment are also
retained in the diagnostic logs. No earlier failed run is called a pass.

Trusted tools were extracted from exact latest dev into `../trusted/tools`,
not loaded from candidate checker/rules/authority. The cumulative diff includes
received #602 code, tests, docs and captures plus S02. Its full PR required set is
`web-change-budget`, `core-pr`, `recipes-standard`, `gallery`, `lint`,
`web-contracts-pr`, `web-pr-smoke`. The unchanged authority/CI contracts are
verified separately. This is local required-job evidence, not a claim of remote
GitHub PR success, complete integrated dev staging or a supported-version matrix.

Cumulative and S02-only trusted policy runs use different bases: latest dev vs
synchronization parent `6e9ef19706a9f0a44a41ccbdc181a368f7b47e78`.
The staged S02-only result is Gate PASS / Review REQUIRED for public exports
and resource-like declarations. The staged cumulative result is Gate PASS /
Review REQUIRED, additionally reporting inherited session/compatibility paths
and production gross churn 1,589 / net additions 637. Both have zero blocking
violations (`policy-staged-s02.log`, `policy-staged-cumulative.log`). Trusted
`required-jobs-staged.json` classifies the actual diff as valid and requires
the seven full PR jobs above. Gate PASS permits the deterministic boundary; Review REQUIRED remains human
review for architecture/scope (including inherited #602 session paths), not a
waiver or a new approval. The exact committed base/head and trusted impact plan
are evaluated before the authorized ordinary branch push. No candidate checker,
authority or CI changes are staged. No PR or dev merge publication is requested.

Reproduction commands from the dedicated repo (export absolute `PYTHON`,
`PYTHONPATH=$PWD`, dedicated Node `NODE_PATH` and venv/Node `.bin` on PATH):

```bash
python tools/prepare_browser_wheel.py
python -m pytest tests/ -m 'not slow' -n4 --dist loadfile --durations=30
mapfile -t web_tests < <(rg --files tests/web | rg '^tests/web/[^/]+\.test\.mjs$' | rg -v '(architecture-contracts|gallery-session-publication)\.test\.mjs$' | sort)
node --test "${web_tests[@]}"
node --test tests/web/architecture-contracts.test.mjs tests/ci/*.test.mjs
GBDRAW_WEB_TEST_PORT=47602 playwright test tests/web/error-boundary.playwright.spec.js tests/web/python-rule-parity.playwright.spec.js --workers=1
GBDRAW_WEB_TEST_PORT=47620 playwright test --config=playwright.pr-smoke.config.js
GBDRAW_WEB_TEST_PORT=47621 playwright test --config=playwright.gallery-publication.config.js
ruff check gbdraw/
WEB_ARCHITECTURE_CHANGE=true node ../trusted/tools/check-web-change-budget.mjs --base <exact-dev-SHA> --head HEAD
WEB_ARCHITECTURE_CHANGE=true node ../trusted/tools/check-web-change-budget.mjs --base 6e9ef19706a9f0a44a41ccbdc181a368f7b47e78 --head HEAD
git diff --check
git diff --exit-code -- tests/reference_outputs/ examples/gbdraw_social_preview.png
```

Production, tests, this non-normative document and generated diffs are reviewed
separately. Only enumerated task files are staged. The generated wheel and
LOSAT mode, dist/egg artifacts, reference outputs and owner-maintained social
preview are excluded. Publication proof (exact local/remote SHA, clean,
left/right counts and trusted postcommit logs) is kept in the dedicated evidence
root and reported at handoff; the enclosing commit supplies this document's SHA.

## S03 connection contract and limits

S03 may use the model directly for summary and optional, initially collapsed
keyboard Details. Manual Copy must use only the displayed bounded code,
operation, stage, context and secondary facts; free Error properties, source
inputs and Session content are not diagnostics. Clipboard unavailability keeps
manual selection and normal recovery. S03 must not create another normalizer,
raw prefix extractor or regex evaluator. A known native public cause should
retain this model instead of being reconstructed as `new Error(String(model))`.
The existing export caller reconstructs type/message/details and still loses
identity; it is explicitly pending S03, as are remaining raw caller/import/
preset/warmup console routes. They are not counted as migrated UI completion.

Action IDs describe existing operations and are not executable callbacks:
`retry`, `save-session`, `review-input`, `edit-input`, `select-input`,
`select-record`, `edit-region`, `edit-depth`, `disable-track`, `edit-table`,
`edit-comparison`, `edit-annotation`, `edit-track`, `edit-pattern`, `reload`,
`generate`, `edit-dpi`, `use-svg`. Bind each to the actual affected operation or
field using current caller/document/row ownership; never infer a target from a
private name in diagnostics. The caller may set the actual root action (e.g.
Align) while retaining the cause code/stage/context; helper identities alone
are not a user action label. Cancel/stale/superseded remain their own outcomes.
Unknown cause/stage/position stays unknown.

First failure vs preserved Result belongs to the actual transaction owner.
The model contains no invented `resultPreserved` flag. S03 must await rollback
and use its observed outcome; rollback/finalization failure cannot claim a
successful restoration. Preserve successful retry error clearing, canonical
request/draft/directions, Result, History, Save/Export, and #602 Status/compact
Editor/review behavior. Cover ordinary and committed-candidate Generate/Align
paths rather than reading global errorLog as a substitute for an operation result.

PD-OI-047 Color draft/Not applied/Retry/Revert remains S04. Its accepted rule,
Result/History/Session and stale/revision semantics remain independent of this
error model. Search semantics, TSV/new-rule/preset behavior and Session schema
are unchanged. PDF/font/snapshot implementation stays with export; shared
files transfer only after this plan's S05 handoff. No S03+ implementation,
release/tag/deployment, main/dev push, PR or Issue closure is performed here.
The original BUG-19 audit pattern/build is still unavailable; this session does
not claim to reproduce that unpublished artifact. Inherited metadata-free
Session and Legend-override limitations are not declared fixed.

S03 start condition: fetch the same branch, verify this result's enclosing S02
commit equals the published remote, clean checkout and current authority/base;
then follow [S03 instructions](sessions/SESSION_03_ERROR_UI_AND_DIALECT.md).
**S02 boundary implementation and verification complete; UI and overall BUG
completion remain pending.**

Commit title: **Preserve structured failure causes across the diagram worker**.
Summary: **Converge producer errors, transport, and bounded user-facing normalization.**
