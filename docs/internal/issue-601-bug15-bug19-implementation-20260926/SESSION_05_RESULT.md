# S05 — Integrated error recovery and regex acceptance

S05 completes every integrated acceptance ID for S02–S04 and updates the
existing Web and TSV references. No production correction was required.
The acceptance table below records independent requirements, inputs, assertions
and reproduction commands. Publication proof for this document's enclosing
commit is written after commit; its own SHA is not embedded through amend.
This record is evidence and handoff, not Product authority.

## Source and isolated environment

- Session root: `/tmp/gbdraw-issue601-s05-lPG7Gv`; independent HTTPS clone,
  objects/index, venv, Node dependencies, Chromium installations, wheel, server
  and evidence. Shared workspace and S02–S04 clones/environments/servers were
  read only. This session is the sole writer of its same-named remote branch.
- Start local/remote/tracking HEAD: `aa3b1207a31ec34e83a58044b3acdbc7f3ba224c`.
- Fetched dev / trusted base: `98c21116f6439d5721e7ea62ae1a21a7cf2d4319`.
  It is already an ancestor; no synchronization merge was needed at start.
- Branch: `fix/issue-601-bug15-bug19`; upstream:
  `origin/fix/issue-601-bug15-bug19`; initially clean, left/right 0/0.
- S04 final HEAD, its implementation `cf3ffe46f4422edbca0f4d3b0e11a51ff6bfb076`,
  authority integration `744be7a5943d4a247d027369629058898aa3f33e`, #602 S07
  `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7` and current dev are confirmed
  ancestors. Guidance and all S00–S05 plan documents match the read-only S04
  checkout byte for byte (`evidence/start-inventory.json`).
- Python 3.13.3 / pytest 9.1.1, Node 26.8.2 / @playwright/test 1.61.1;
  Python Playwright 1.61.0. The dedicated install uses Biopython 1.88,
  pandas 3.0.6, SVGwrite 1.4.3 and Ruff 0.15.12. Exact dependencies are in
  `evidence/python-freeze.txt` and npm's unchanged lockfile.
- Node CLI/version, Python Playwright import, Node module resolution and clone
  import path were checked separately. The dedicated server's dynamically
  selected port is **37739**, recorded in `evidence/port.txt`; every response has COOP
  same-origin, COEP require-corp and CORP same-origin. No external genome upload
  or external runtime dependency was introduced.
- The ignored browser wheel was prepared once before browser tests and never
  regenerated during verification (`wheel.log`, `wheel.sha256`). No cache-bust
  update. LOSAT's executable bit was prepared only in this clone and restored
  to its tracked mode before staging, with unchanged binary bytes.
- Bubblewrap fails before command execution on this host's `/mnt/wslg/distro`
  mount. The same local commands used required escalation; this is an
  environment limitation, not a runtime fallback or browser skip.

## Authorized behavior and owner review

Developer preflight: **IMPLEMENT_EXISTING_AUTHORITY**. OIPC revision 24 on
current dev already selects the complete PD-OI-046/047 outcomes. No material
outcome remains unresolved, no new judgment or evidence-dependent behavior is
selected, and no prohibited runtime/guard combination is proposed. The Product
map and BD store remain unchanged; no BD number is invented.

The complete sections for PD-OI-016/024/027/029/031/034/035/036/037/038/039/044/
045/046/047 are equal at start, trusted dev and candidate. PD-OI-046 and 047 each
retain all nine JSON fields exactly equal to their signed source, including
concern, revision, choice, rationale, preservation, retirement, risk, owner and
date (`receipt-verification.json`). The export handoff and original PDF/error
receipts were read at `f0c8128ac9d252bbc1b64074cf995f536f789d03`, without importing
that branch or changing/retiring/remapping its receipts.

S05 changes no production owner or path. Public wording/model remains in
`services/error-normalization.js`; Python producer/adapter and Worker/client
retain their finite transport; Generate/Align retain their orchestration and
artifact transaction owners. `feature-editor/rule-actions.js` owns rule mutation
and its focused `pattern-drafts.js` owns transient field display. Canonical
request, Result admission, History, Session and export remain separate accepted
owners. No classifier, translator/evaluator, generic draft manager, second
controller, schema, persisted rule ID, compatibility reader or runtime dependency
is added. OE/PE/CB do not increase; no exception condition applies.

S02 removed raw serializers/prefix extraction; S03 removed reconstructed caller
errors and operation-specific raw prefixes; S04 removed unconditional pattern
restoration and index-bound field edits. Those superseded paths remain absent.
S05 removes no additional runtime path because no runtime correction was needed.
Rollback is a revert of the S05 commit, retaining previous runtime and authority.

## Commands and evidence

All paths below are relative to the dedicated root's `evidence/` directory.
Tests ran from its `repo/` with the dedicated venv and Node CLI first on PATH,
absolute `PYTHON`, this clone's `PYTHONPATH`, dedicated `PLAYWRIGHT_BROWSERS_PATH`,
and the port from `port.txt` in `GBDRAW_WEB_TEST_PORT`. Browser retries were zero;
existing test timeouts, geometry/performance thresholds and failure assertions
were preserved. Long runs were monitored incrementally.

| Ref | Command | Result / log |
| --- | --- | --- |
| N | `node --test` on the top-level `tests/web/*.test.mjs`, excluding only `architecture-contracts.test.mjs` and `gallery-session-publication.test.mjs` (the required CI selection) | 871 passed; `node-required.log` |
| N2 | `node --test tests/web/error-normalization.test.mjs` after the new numeric/cleanup assertions | 1 file passed, all native and transport assertions; `node-boundaries-final.log` |
| A | `node --test tests/web/architecture-contracts.test.mjs tests/ci/*.test.mjs` | 204 passed (139 architecture + 65 CI); `node-architecture-ci.log` |
| C | `python -m pytest tests/ -m 'not slow and not browser' -n 4 --durations=15` with dedicated basetemp/cache | 6633 passed / 17 skipped / 10 warnings, 225.22 s; `python-core.log`; includes core, recipes, Gallery and read-only SVG comparison |
| B | Node Playwright: `error-boundary.playwright.spec.js python-rule-parity.playwright.spec.js specific-color-pattern-drafts.playwright.spec.js --config=playwright.functional.config.js --workers=1 --retries=0` | 31 passed, 5.3 min; `browser-integration.log` and `browser-integration/` |
| S | Node Playwright: `--config=playwright.pr-smoke.config.js --workers=1 --retries=0` | 19 passed, 3.6 min; `browser-pr-smoke.log` |
| PB | `python -m pytest tests/ -m 'browser and not slow' --durations=10` with dedicated basetemp/cache | 38 passed / 6661 deselected, 253.94 s; `python-browser.log` |
| RB | Node Playwright: `tests/web/similarity-alignment-ui.playwright.spec.js --config=playwright.functional.config.js --workers=1 --retries=0 --grep 'Gallery explicit directions\|exclusive directions include\|minority reference, source strands\|changed final facts\|compact review preserves'` | 7 passed, 3.0 min; `browser-received.log` and `browser-received/` |
| L | `ruff check gbdraw/` | PASS; `ruff.log` |
| F | `python ../verify-receipts.py`; `git diff --check`; read-only artifact diff | receipt equality PASS; final scope/whitespace reviewed before commit |

New blocking assertions run in existing required paths: numeric boundary and
secondary-cap assertions are in the top-level N file selected by
`web-contracts-pr`; real Align engine-error and catch assertions are in the
existing `@pr-smoke` case selected by `web-pr-smoke`. The inventory remains
19 smoke cases, all members of full functional acceptance, and 16 comparison
contracts in the separate Python browser wrapper. No tag, workflow, checker,
CI guard, threshold or upper-bound change was needed. PR #625's accepted ceiling
correction is reused as dev content rather than reimplemented.

S04's final logs and publication receipt were read and hashed
(`prior-evidence-review.json`). Their ownership, authority and inherited limits
are retained as historical facts. No earlier test count substitutes for S05
acceptance: new Python dependencies make broad environment reuse inappropriate,
so the relevant local gates, integrated browser scenarios and received
continuations are run in this session. Existing unchanged test assertions remain
the oracles; no oracle, reference SVG or timeout was regenerated to fit output.

Reproduction uses `node node_modules/@playwright/test/cli.js test` before the
browser arguments in B/S/RB above. N's literal selection was:

```sh
node --test $(rg --files tests/web -g '*.test.mjs' | awk 'index(substr($0,11),"/")==0 && $0 !~ /architecture-contracts.test.mjs|gallery-session-publication.test.mjs/')
```

## Acceptance by requirement

Each row is a separate contribution; matching an option ID or a suite count
does not substitute for its assertions. All commands refer to the table above.
No acceptance row relies only on a successful stubbed `runAnalysis`.

| ID | Actual input / operation | Required assertion and evidence |
| --- | --- | --- |
| E01 | B real render with a color-table `😀[` sentinel and EVALUATE_RULES Color/Label on empty/unrelated catalogs; Worker serialization/client deserialization; N native embedded Python oracle; S Align caller/UI/Copy | Exact REGEX_SYNTAX, evaluateRules/generate/align, rule-validation, Python position 1 and Label row survive; renormalized model equals original. B `real Python adapters…`, N2 transport and S real caller facts. |
| E02 | N2 real option/annotation/table/glyph/native validators plus second normalization | Existing known INPUT_INVALID/ANNOTATION_TARGET/TRACK_INVALID/PDF_GLYPH and finite correction/action retained; unclassified validation distinct from UNKNOWN; one public prefix, exact code/context/stage after transport. Native Python adapter tests are in C. |
| E03 | N2 repeated private text in name/message/cause/stack/traceback/stdout/stderr/notes/details/context and 20 cleanup entries; numeric values 0, maximum, maximum+1, fractions, strings, infinity/NaN; B/S real private pattern/file/clipboard/admission sentinels | No sentinel in summary/Details/Copy/automatic console; finite domains, summary ≤1000, details ≤4000, ≤8 sections, ≤2 cleanup facts, numbers ≤10^7 (codepoint ≤U+10FFFF). Missing/UTF-16 positions omitted; Python character position 1 retained. Original text stays only in its user-input field. N2, B four mode/width cases, S. |
| E04 | N real run-analysis/artifact owner with first failure, preactivation failure, postactivation adoption failure, awaited readiness rollback and rollback rejection; native helper cleanup and render cleanup | recovery no-result/preserved/restored/restore-failed are separately asserted; initiating cause survives cleanup; restore-failed never claims preservation. B first Circular generation vs existing Result, N2 proxy destroy/cleanup-only and N `run-analysis-simple-path` assertions. |
| E05 | S BGC Gallery → actual review → real native render failure → real successful render with one admission-owner exception → corrected Apply | REGEX_SYNTAX/rule-validation then COMPARISON_IDENTITY/result-admission plus cleanup fact reach Align UI/Copy. Every Result, committed Session/request, plan/receipt/direction and both History stacks match the prior artifact; local choices unchanged and error focused. Corrected retry clears error and adds exactly one History action. Only History's invalidation revision may advance, as documented by its owner. Producer inputs and admission throw are deliberate test injections; the real run-analysis, Worker, Python and transaction are not replaced by success stubs. |
| E06 | B held Worker replies across drawer, mode cycle, row removal/reorder, successful Session replacement, new edit/keystroke and cancellation; cold cancel/retry; N stale/superseded/cancel matrices | No invisible-mode, deleted/replaced row, newer edit/document commit; no new failure for cancel/stale/superseded; accepted state exact. Cold retry constructs a new actual Worker only after the canceled one is terminated. N/B races. |
| E07 | B/S keyboard Details, explicit Copy, unavailable/rejected Clipboard, Select diagnostics and real recovery buttons | Correct operation title; summary always visible, Details initially collapsed, Copy hidden until opened, only displayed bounded text copied, readonly manual selection and recovery remain. B helper `inspectSafeDetails`, field cases and S Align. |
| R01 | B regex_rules.gb corpus with `(?i)NADH`, `(?P<enzyme>NADH)`, `NADH\Z`, `\bβ`, Unicode `i`; native Color/Label in C/N | Expected exact locus targets R0/R1, R0, R2, R3/R4/R5; Color fills and Label text agree before and after Generate; native helper/renderer parity. B `Unicode regex corpus…`, C `test_web_rule_matching.py`, N rule preparation. |
| R02 | B invalid `[`/JavaScript named group on empty/unrelated catalogs; manual Color, specific TSV, Label TSV, Python preset, accepted History/Session/Generate | Syntax rejected independently of matching; TSV structure/resource errors distinct. Failed imports/Retry preserve accepted file/rules/Result/History; corrected imports add one entry. Manual/new/preset meanings unchanged. B adapter/parity/import/preset cases, N2 native validators, C. |
| R03 | B loaded mtDNA qualifier-value/product Search, JS `(?<enzyme>nadh)`, rejected Python `(?P<enzyme>NADH)`, word NADH, actual Interactive SVG download and file document | Exactly the same seven feature IDs in live and downloaded JS/word search. Actual bytes contain JS dialect and correction wording, runtime controls use same targets. No Python runtime added to standalone. B `live and downloaded standalone search…`. |
| R04 | B existing Color field rejection, failed Retry, Revert, successful Python correction; initialization and staging failures | Display text separately retained; accepted rule/Result/History exact; Not applied and syntax aria-invalid vs WORKER_INIT/RESOURCE_INVALID stage; four existing helper calls for two attempts, no added Worker; correction atomically adds one History entry, Revert no evaluation. B four mode/width and real runtime cases; N draft actions. |
| R05 | B drawer/mode retention, unrelated form History, target Undo/Redo, deletion, reset, successful Generate/Session, failed import preflight and post-reset rollback; held response races | Retain only same-document/surviving unchanged rule drafts; release replaced/deleted rows/documents; failed Session restores source before draft. Separate displayed text, canonical rules, whole Result and History assertions. B `real History…` and held edits; N actual-row/duplicate/current token tests. |
| R06 | B rejected NADH draft → Save → fresh context Load → SVG Export → Generate | Saved JSON excludes draft sentinel, retains accepted regex; fresh Load zero Worker, same accepted canonical and whole sanitized Result; actual Export XML equals current Result; Generate uses accepted rule and clears replaced-document draft. Export does not clear current draft before Generate. B `real History…`, downloaded boundary JSON/SVG. |
| U01 | B Circular/Linear at 1600 and 390px, keyboard edit/change/Retry/Revert, focus, field IDs, Details/Copy/manual selection; S narrow Align | aria-describedby points to the actual row cause/status; syntax-only aria-invalid, pending is Checking/Not applied; Retry/Revert/Copy reachable, successful edit/Revert focus retained. Real screenshots visually reviewed for Linear 390px field and generation alert; no public images were altered. |
| P01 | B 25000 features, alternating β-lactamase/other, one regex preparation then synchronous reuse; N draft 25k and keystroke/Revert | Exactly 12500 matches, one matching evaluation, one lazy Worker, reuse true, event loop responsive. S05 measured 477 ticks, 7662.87 ms preparation, 22.16 ms reuse. Field attempts retain existing captions+matching calls, no validate-only operation or keystroke/Revert call. |
| SCI01 | C native comparison/ValueError catch, native re.error/ParseError causes, complete read-only reference SVGs; B real targets/generated diagrams/Session/exports; N request and admission ownership | Native capture and scientific targets/geometry unchanged, current Session/request schemas and cache identities retained. No new runtime path or reference regeneration. Corpus targets and whole-artifact boundaries complement successful core comparisons. |
| G01 | N/N2/A/C/B/S/PB/L/F plus trusted exact dev checker and CI classification | Local required checks and evidence boundaries are recorded separately; no candidate authority/checker executes to admit runtime. Gate PASS and Review REQUIRED are distinct; remote CI, human review and dev staging are not inferred. |

## Independent receipt coverage and inherited limits

| Source / independent condition | Coverage |
| --- | --- |
| PD-OI-046 known corrections, actual stage, unknown stable code, private diagnostics/console, cleanup primary | E01–E04; N2/C/B/S. All migrated known failure families in S02's [migration table](SESSION_02_RESULT.md#known-failure-migration) remain in the same producer/adapter/normalizer; no known family is collapsed to UNKNOWN. |
| PD-OI-046 Result/request/local draft/orientation/History, retry, Save/Export, cancellation and first-result distinction | E04–E06/R06; actual transaction and artifact assertions, not diagnostic flags. |
| PD-OI-046 keyboard/manual displayed-only Copy, Clipboard fallback | E07/U01, independent of cause transmission. |
| PD-OI-047 Python semantics/priority/one preparation/Worker/valid targets/atomic commit/reuse | R01/R02/R04/P01/SCI01; helper targets agree with Generate and native owners. |
| PD-OI-047 display/accepted separation, Not applied, syntax/runtime, successful History only, field lifecycle and nonpersistence | R04–R06/E06/U01; each state boundary independently asserted. |
| Original export error receipt: optional initially collapsed Details, always-visible summary, keyboard/390px, native capture, local-only, actual stage, recovery/draft/cancel | E01–E07/U01/SCI01. These remain independent requirements despite its different option ID. Source receipt remains intact; it is not a second active authority. |
| Original PDF receipt: text/font/style/Unicode/geometry/snapshot/filename/DPI/local-only/lazy font/retry/glyph/packaging evidence | Preserved as export-owner obligations, not certified by S05. BUG-07/PDF producer/font/snapshot/distribution implementation is outside this acceptance. |
| #598 Apply/Reset/direction/receipt/retry and #602 compact Editor/Status/canonical projection/Session exclusion/independent History and Result | Existing independent authority sections unchanged; N/C/PB/S retain canonical and native contracts; RB freshly checks both Reset scopes, source strands, minority-reference direction, receipt/fresh Load/ribbon geometry, changed-fact separate Apply, and three compact entry scenarios. Their canvas, local no-Worker choices, focus, Editor availability and recovery assertions pass. No second controller, Match restoration, guessed direction or Session schema. |

The original BUG-19 audit pattern, field and build remain unidentified in the
current Issue and available source history. Current parity/recovery acceptance
is not reproduction of that original artifact or a deployed-site verification.
Standalone module-fetch failure still needs reload because its rejected import
promise is cached. Metadata-free Session and Legend-override limitations remain
unresolved. Firefox/WebKit, supported Python/Node matrix, physical OS keyboard,
physical zoom, remote S05 CI and exact-dev staging have not been certified.

## Separate diff reviews and delivery boundary

Production: no authored runtime change; existing owners and superseded-path
removal retained. Tests: numeric transport/privacy limits and real Align catch
assertions supplement existing accepted-state/History/retry checks; no timeout,
threshold, known classification or rollback assertion weakened. Docs: only
existing Web and input-schema owners gain descriptive specifications. No new public
manual/page, tutorial steps, screenshots or Gallery image were added. Generated:
ignored fixed wheel and private QA outputs only; references, Gallery, dist,
egg-info, owner social preview and tracked binary bytes/mode unchanged at stage.

Trusted tools are archived from the exact dev SHA above into `trusted-dev/`.
The start full-SHA diff is Gate PASS / Review REQUIRED (`gate-start.log`), and
`ci-start.json` is profile=pr / impact=full. Final cumulative working-tree check is also **Gate PASS / Review REQUIRED**
(`gate-final-working.log`); S05-only is **PASS / CLEAR**
(`gate-s05-working.log`). No guard files are touched. Exact committed checks
follow this record, using the same trusted tools and full base/head SHA. Required jobs are
web-change-budget, core-pr, recipes-standard, gallery, lint, web-contracts-pr,
web-pr-smoke; their local evidence is N/N2/A/C/B/S/PB/L. The exact committed
base/head report and final CI plan are in the publication handoff. Gate success
does not fulfill human Review REQUIRED.

See [export owner handoff](EXPORT_OWNER_HANDOFF.md) for shared-file acquisition
and preservation. Post-commit normal-push proof, exact local/remote/tracking SHA,
clean state and left/right are saved in the external
`evidence/SESSION_05_PUBLICATION_HANDOFF.json` and final session report. No amend
for a self-referential SHA. No runtime PR, dev/main integration, release, tag,
deploy, Issue close or external owner message is authorized or performed here.
The earlier #625 CI PR is completed base content, not authorization for a new PR.

Commit title: **Verify error recovery and regex editing across web workflows**

Summary: **Complete integrated acceptance and update existing behavior documentation.**
