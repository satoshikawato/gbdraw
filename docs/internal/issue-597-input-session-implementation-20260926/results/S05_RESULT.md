# S05 result — Resumed import Worker implementation

Status: **runtime implementation, owned measurements and result preservation complete; S05 acceptance / admission remains BLOCKED / FAIL**. The exact Worker permissions are integrated; they are no longer a missing prerequisite. Protected characterization delivery, existing mapped browser failures, real preview Python startup, real Generate and full-pipeline performance remain unresolved. Scope is BUG-02 / BUG-20, S05 only. BUG-01/S02 and S06+ were not started. The earlier prerequisite-only result is retained verbatim below as historical evidence.

## Checkout, integration and authority

Reused `/tmp/gbdraw-issue597-S05.ujOWyl`, initially clean on `fix/issue-597-input-session-20260926` with upstream `origin/fix/issue-597-input-session-20260926`. The repository, branch, ownership and saved prerequisite objects were checked before edits. Shared checkout and authority checkout were not used as implementation targets; unrelated changes and other sessions' caches/servers were preserved.

| Reference | Actual SHA / action |
| --- | --- |
| Start / initial actual remote work branch | `7209d7e2ad70ac71eaaab5409e75f2f7aab126b3` |
| Saved dev at intake | `88028fd242d263f0fe86aaf9da57b8dc9eb082f6` |
| Actual trusted remote dev / PR #620 merge | `d313b70b9f97c2c1d70f9ae885edbead80b62021` |
| Permission-only commit | `bb5b53f4b3d60e5952161d2e3a6ca05c4f6541bb` |
| This session's dev integration merge | `c11d4854a899349cd66863c859f6f871724ba3ab` |
| Measurement source | Integration HEAD above plus full runtime file hashes in [fingerprints](../evidence/S05-resume-fingerprints.json) |

Only dev was fetched because its actual object/ref was missing locally. The already-current work branch was not fetched/pulled, and no new clone, branch, main or tags fetch was needed. The clean-tree dev merge resolved six conflicts by retaining both S04 coordination and dev's annotation/color/History changes. Its commit is separate from the one implementation/tests/evidence/result commit. During validation the actual remote still matched these SHAs. Final inspection found dev advanced to `252986d096011fcf1a0f5564e940480d3b92844d` through PR #621: only signed Issue #619 decision text and revision 23, with PD-OI-044/045 unchanged. Only that new dev was fetched. To preserve the owned staged implementation and avoid a dirty-tree merge or stash/reset, its verified authority-only change is merged after the implementation commit into the clean same branch. This additional integration SHA, final implementation SHA, postcommit policy, nonforce push and actual remote equality are recorded at handoff in `/tmp/issue597-S05-resume-evidence`; no document self-SHA is required.

S01 `15bcbca89392cbd63a88fe80dc44b71ad4868061`, S03 `e3684d496df7a3b2a7492443cad50168f6959b36`, S04 integration `eed70f05d04f638a5a4ab6c6567982d45e774596`, S04 runtime `b092499d27705092d84037d5acd51032177c619a`, prerequisite evidence and PR #610 merge `af5d942af60353dda199aa487da9152a3576b3fe` are HEAD ancestors. PD-OI-044/045 human/serialized/static receipts match all nine fields each, with measured contract revision 22. The late separate trusted-dev amendment advances the final contract to revision 23 without changing either Issue #597 receipt. Product outcomes were not reselected or expanded.

[S05-resume-readiness.json](../evidence/S05-resume-readiness.json) records exact trusted-dev permissions: `allowedPrivilegedOwners["Diagram Worker"]` includes `services/session-import-client.js`; `allowedPrivilegedImporters["services/session-file.js"]` includes `workers/session-import-worker.js`. PR #620's prior CI / 192 focused PASS remains authority-delivery evidence only. No permission patch was reapplied. The readiness snapshot is explicitly **after integration, before runtime**, not a current no-runtime claim. At the measured checkpoint, active checker, policy, authority, guards and mapped contract bytes equal trusted dev d313. The late revision-23 authority arrives only through the separate clean dev merge; its before/after hashes and unchanged Issue #597 decision sections are recorded in the manifests.

## Runtime owners and preserved behavior

Selected transport: **File/Blob → JS import Worker → existing read/gzip/fatal UTF-8 codec → JSON.parse → whole-object reply → private untrusted main candidate**. Lifecycle is bounded; transport remains whole-object. No fallback, bounded sections, ACK, sampling pause, terminal flush or persisted protocol was added.

| Owner / path | Result |
| --- | --- |
| `services/session-import-client.js` | One owner for operation IDs, stale/duplicate/late replies, structured errors, cancel, crash, unreadable reply, single settlement, listener cleanup and termination. Every success/error terminates before Promise continuation reaches preflight/adoption. Unsupported Worker or startup/post failure is explicit; no main parser fallback. |
| `workers/session-import-worker.js` | Reuses `session-file.js` for read/decompression/fatal decode and performs JSON.parse. Returns the entire graph plus character count/observational timing. Imports no config, DOM/Vue, request, History, SVG or Python owner. Parse errors do not quote uploaded JSON. |
| `services/config.js` | Replaces/removes the main full-text read/JSON.parse path in the same change. Validates unsafe keys before candidate use. Adds AbortController to the existing operation lifetime, not another lock. Existing preflight, migrations, adoption, rollback and History owner remain. Save retains its existing codec/streaming path. |
| Settings-only export | Both-mode browser tests exposed dev's new empty `annotationWarnings` in runMetadata rejecting settings-only Save. Settings-only export now emits its required empty runMetadata; full-document metadata remains unchanged. No validator was weakened. |

The S04 availability owner still governs Save/Load/Generate/reflow, same-promise Save join, cancellation and post-await current checks. Safe adoption and source-first rollback preserve the old request/resources/Result/History on rejected candidates. Ordinary architecture evidence: one lifecycle owner replaces the main parse execution path; one codec and one semantic operation owner remain. The two new modules match exact preauthorized paths. No alternative parser, compatibility branch, global framework or ratchet exception was introduced. The inert guard patch is an external delivery candidate, not active admission authority.

## Reproduction and fingerprints

[Reproduction](../evidence/S05_RESUME_REPRODUCTION.md), [measurement recipe](../evidence/S05_MEASURE.py), [validation](../evidence/S05-resume-validation.json), [fingerprints](../evidence/S05-resume-fingerprints.json) and [artifact hashes](../evidence/S05-resume-artifacts.json) bind commands, environment, source/input/recipe hashes, logs and generated outputs. Raw owned outputs are `/tmp/issue597-S05-resume-evidence`; large biological payloads, logs, wheel and caches are not committed. Some test shell strings were not separately retained: the validation manifest explicitly distinguishes recovered reproduction selections from exact measurement/semantic commands. It does not invent shell flags or reinterpret failures.

Environment: Python 3.13.3, Node v26.8.2, Ruff 0.15.12, Playwright 1.61.0, Chromium 149.0.7827.55, WSL Linux / i9-14900HX, 32 logical CPUs, 33518579712 bytes RAM. CLI and Python Playwright paths were both checked; Node Playwright was installed only in this checkout with an owned npm cache. The ignored browser wheel was prepared without cache-bust refresh. Default bubblewrap failed mounting `/mnt/wslg/distro`; the same commands ran with escalation. A read-only approval timed out once; a simpler retry succeeded.

The unchanged S01 fixture recipe regenerated the same six pinned GBFF sources, first two distinct chromosomes each: **12 distinct DNA records, 62089 biological features, 132 nonself + 12 self = 144 protein cache entries**. Combined source SHA-256 `93fc27fc1ba27f21e513079b506244b072d5f836889acc0dfb5feca5fe7df67c`. Owned LOSAT/cache generation completed normally. Regenerated gzip SHA-256 `d3cafef9664bff958aec4932eea6e264ae872cce89762a7cbe617447af975388`, **134471323 compressed / 466767669 expanded bytes**. Historical S01 sizes **134471293 / 466767621** are retained separately; this is a new run, not an overwritten fixture. Large positive import is gzip only; expanded plain exceeds 200 MiB. Tracked LOSAT mode remains **0644**, unchanged, and its fixture child processes completed.

## Focused checks and required gates

| Invocation | Actual result / qualification |
| --- | --- |
| New client + real production-Worker codec Node tests | **8 PASS**: plain/gzip magic/Unicode, malformed JSON/gzip, fatal UTF-8, unsafe own keys, exact caps, success/error/crash/messageerror, stale/late/duplicate, cancellation/20 repeat cycles, unsupported/start/post failures; no fallback. |
| Full top-level Web Node invocation after adapter fixes | **1000 tests: 998 PASS / 2 FAIL**, 176.48 s. Both failures are unchanged protected Worker inventories omitting the authorized new constructor. Earlier invocation **980 PASS / 20 FAIL** remains historical; 18 adapter failures were fixed with a production Worker Node adapter, not a fake main parser. Corrected import-focused rerun **45 PASS**. |
| Python focused | **459 PASS**, 10 warnings, 55.26 s; API/session/settings-only/codec/current-released compatibility and record/profile scope. |
| Read-only `TestOutputComparison` | **16 PASS**, 33.90 s; references unchanged. |
| Ruff production / measurement recipe | PASS. |
| Browser initial / corrected | Initial **36 PASS / 3 settings-only FAIL**; corrected **24 PASS**, 11.2 m, including composite resources, current/released CLI compatibility, settings-only and Worker/coordination cases. Unaffected initial S03 discovery/disclosure/record, auxiliary, palette, loading and Save lifecycle PASS is reused with its scope identified. |
| New actual Chromium Worker assertions | **3 PASS**: five alternating JSON/gzip loads; termination precedes preflight; zero retained import Worker/listeners/object URLs; rejected unsafe/malformed/fatal/crash preserves document/History; pending teardown cancels before adoption and retry loads. |
| PR smoke / CI Node / sequence source coverage | **13 / 65 / 3 PASS**, respectively. Tobacco and Hepatoplasmataceae Session projection commands also exit 0. |
| Packaging non-slow/non-browser | **34 PASS / 12 deselected**, 13.64 s; not a full build/offline/browser suite claim. |
| Unchanged mapped browser combined | **1 PASS / 2 FAIL / 2 not run**; completed runner report, outer tool exit **143**, no shell exit marker. Feature fill commits no History; direct edit reaches its unchanged 600000 ms timeout. |
| Mapped Feature fill on trusted dev / observational diagnosis | Same FAIL and same alert on dev/candidate: `Cannot apply feature style: Legend entry "tRNA" already exists with a different color.` Existing `syncFileLegendEntries` rejects this conflict; outside the import transport repair. |
| Separate unchanged divergent / legacy | Divergent **FAIL**: busy `Updating diagram. Retry after the update finishes.` Legacy **PASS**. Three inert draft experiments also FAIL; last no-draft PASS, divergent FAIL, legacy skipped. Rejected patch stays outside git. S04's passing coordination contract still requires reflow busy behavior. |
| Inert protected inventory candidate | Disposable candidate-source **4 PASS**, trusted-dev-source **4 PASS**, including two unchanged mapped Node bodies. Not applied to active tests; not external review or CI. |
| Unchanged Vibrio single-Save budget | **1 PASS**, narrow budget: wall 4930 ms, heap delta 95291768 bytes, max heartbeat **916.9 ms**. Overall 500 ms target still FAIL. |
| Unchanged Hepatoplasmataceae regeneration performance contract | **1 PASS**, 36.8 s; exact five-record fixture generates twice with bounded work and reuses one Python Worker. Structural contract PASS is not Issue #597 whole-pipeline performance PASS. |
| Local base checker before/integrated/runtime unstaged/staged | All **exit 0 / Gate PASS / Review REQUIRED**, against actual dev `d313b70...`, unchanged checker/authority bytes. Final precommit and postcommit `--head HEAD` run at handoff; neither is trusted CI or waiver. |

No protected guard, mapped body, budget, timeout or Worker-free criterion was weakened. [S05_GUARD_DELIVERY.md](../evidence/S05_GUARD_DELIVERY.md) records the ready exact inert inventory patch and the separate unresolved browser boundaries. The guard-only patch needs separate review/checks/merge into trusted dev before runtime admission. Any mapped test/reference update follows its required evidence-only → authority-ref-only sequence. This request does not authorize those external PRs/branches or a new Product outcome. Pushing the authorized work branch does not trigger the dev/PR-only required workflow, and no remote CI PASS is claimed.

## Runtime and whole-pipeline measurements

Runs were **sequential within this session**, using owned server/browser/CDP/cache. The shared host is not certified idle. S01's 100 ms heartbeat, long tasks, main heap, worker used/backing memory and independent 100 ms RSS sampling are reused. A test-served observational codec/Worker trace records stage spans, URLs and constructor/post stacks without contents. Its initial smoke route failed to expose the codec probe; corrected smoke passed before accepted measurement. Completion excludes S01's 200 ms observer settlement; graph equality is checked afterwards. Overlapping read/decompression spans are not additive.

| Operation, three runs each | Completion / wall ms | Heartbeat p95 / maximum ms | Conclusion |
| --- | --- | --- | --- |
| Small plain JS transport | 15.6–18.5 | 99.9–100 / 99.9–100 | PASS transport responsiveness |
| Vibrio gzip JS transport | 555.8–684.0 | 150.6–196.4 / 150.6–196.4 | PASS transport responsiveness |
| Real gzip JS transport | 2966.6–3065.4 | about 100.1 / **512.2–526.0** | **3/3 FAIL max ≤500**; whole graph equal |
| Vibrio full Load | 17197.5–18196.9 | 16577.7–17534.1 / same | **FAIL** |
| Real full Load | 118212.1–146263.0 | 100.1–100.2 / **99702.6–124875.5** | **FAIL**; p95 alone hides the long stop |
| Vibrio full Save | 5038.1–5542.4 | 195.5–213.6 / 195.6–223.8 | Responsiveness PASS in these three observations |
| Real full Save | 21413.7–25258.2 | 189.0–217.8 / **1779.4–2146.3** | **FAIL** |

[Transport](../evidence/S05-resume-transport.json) records 9/9 whole-payload equivalence, one import Worker termination each, zero retained listeners. [Limits](../evidence/S05-resume-limits.json) records actual 466767669-byte plain-file rejection and 513 MiB gzip-expanded rejection (compressed 522893 bytes), both with termination/cleanup. A small plain positive is separate from these negatives.

Real transport sampled main heap is about **468.8 MB**, import Worker used heap **906.75 MB**, aggregate Chromium RSS **2.275–2.333 GB**. Real full Load main sampled peak reaches **1.628 GB**, import Worker used **906.75 MB**, import backing **554–586 MB**, Python used about **22–24 MB** / backing **106–108 MB**, aggregate RSS **3.264–3.507 GB**. Real Save main sampled peak reaches **2.066 GB**. These are sampled lower bounds across different notions of memory, not summed exclusive RAM or certified retained size. Closed-Worker CDP sampling errors stay recorded; unobserved short-lived isolates remain unknown. Native structured-clone **wire bytes, copy bytes and receiver task duration are null**, not inferred from file size. Explicit transferred-buffer bytes are zero because production uses no transfer list, not because native copying is zero.

[Pipeline](../evidence/S05-resume-pipeline.json) records the exact real preview Python trigger in every run: `config.js::importSessionDocument` → `validateUnmanagedConfigOverrides` → the app-setup registered validator → `runDiagramHelperOperation("validateConfigOverrides")` → `runAuxiliaryWorkerRequest` → `ensureWorkerInitialized` → `diagram-generation.js::getWorker` → `workers/diagram-generation-worker.js`. Constructor/post stacks and helper operation corroborate it. This is separate from the terminated JS import Worker. **Real preview Python count 1 FAILS required 0**; Vibrio preview count 0 passes. The prior sequence-recovery hypothesis is not substituted for this causal observation. Critical config validation was not skipped or retired to obtain zero.

[Semantic evidence](../evidence/S05-resume-semantic.json): all **nine checks pass for each of six source-versus-Save documents**, including request/resources/catalog/cache/manifest/overrides and saved SVG. Python readers reconstruct 4/12 records and both dataset CLI replays exit 0. [Journey](../evidence/S05-resume-journey.json): fresh Vibrio Load → Generate → Save passes and generated output separately replays. The real fixture's Generate fails with **`Sequence #1: Missing GenBank file.` on both trusted dev and candidate**. It contains full committed resources/request, but schema-2 active `webFiles.bindings.linearSeqs` is empty. Empty draft bindings were preserved, not silently replaced with committed inputs; the retained preview saved after the error is not a generated positive. Source-versus-itself checks on the generated Vibrio replay are not counted as Save fidelity.

Additional CLI SVG comparison reports **three strict FAILs**, each solely CLI root `baseProfile="full"` absent from Web SVG. An additional disposable comparison excluding only this root attribute passes all remaining scientific geometry/text/styles/metadata. Both results and the reproduction are preserved; the existing asymmetric verifier was not edited. CLI execution PASS is not silently upgraded to strict SVG PASS.

## Acceptance and remaining boundaries

| ID | Current disposition |
| --- | --- |
| S-01 | Worker JSON/gzip, caps/fatal decoding/unsafe keys and current/released/settings-only import checks PASS. Full real fresh-Load/Generate/replay acceptance remains incomplete because real Generate fails; mapped divergent/Feature-edit cases remain unresolved. |
| S-02 | Transport lifetime/crash/cancel/stale/late/duplicate/termination and actual-browser repeated cleanup PASS. S04 busy/adoption/rollback checks pass; no claim that all mapped edit/Save continuations pass. |
| S-03 | Required real fixture and Vibrio measurements complete; native clone bytes/task duration remain null and CDP coverage is explicitly partial. |
| S-04 | Six real source/Save semantic comparisons and reader/CLI execution PASS; Vibrio generated journey PASS. Real Generate and strict CLI SVG comparison remain FAIL as recorded. |
| S-05 | **FAIL**: real transport max heartbeat, full Load/Save responsiveness and real Python-free preview are unmet. Historical S01 whole-pipeline FAIL remains open. |
| D-01–D-04 | S03 browser discovery/disclosure/record presentation passes in this resumed runtime's initial browser invocation; unchanged paths retained. |
| A-01 | Exact permissions/receipts present and local Gate PASS / Review REQUIRED. Required protected characterization and mapped browser evidence remain blocking; no self-authority or waiver. |
| W-01 | Existing correct checkout reused, targeted dev fetch/merge, explicit staging and one runtime commit; nonforce same-target push/remote equality checked at handoff. |

Historical S04 full Python **6262 PASS / 6 FAIL / 17 skipped / 11 deselected** and the obsolete revision-21 probe FAIL are retained. The earlier full Node and settings-only browser failures, first codec smoke FAIL, failed mapped runs and rejected draft experiments are distinct invocations, never relabeled. Full Python, full functional/offline suites and external required CI were not claimed rerun/passing.

Read [guard delivery](../evidence/S05_GUARD_DELIVERY.md) first for the admission dependency, then the readiness/fingerprint/validation manifests and measurement JSONs for current evidence. Complete separate protected guard delivery and scoped mapped browser resolution without weakening accepted reflow busy behavior or choosing new legend outcomes in S05. Real Python-free preview and fixture Generate require a separately scoped resolution retaining config validation and draft/committed meaning. These boundaries are not missing PR #620 permissions.

**S06 was not started.** Its concrete later entrance is: real current preflight **8409.7–9875.3 ms** (catalog **6193.3–7047.6**, artifact validation **1923.3–2575.5**); restore-before-SVG **6666.9–8152.5**; SVG admission **1062.4–1172.9**; DOM preview mount **92914.5–117150.2**; Save projection **386.8–482.5**, compression **14796.3–18033.0 ms**. Vibrio DOM mount is **10954.7–11420.8 ms**. Stage spans overlap and do not explain wall time by simple summation. S05 transport termination already precedes these stages; they must not be hidden by reclassifying transport performance.

Production, tests, internal docs/evidence and generated diffs were reviewed separately. No public Gallery/tutorial/reference SVG/social preview, package/dependency manifest, dist/egg-info, tracked generated wheel/cache-bust or LOSAT mode changed. The ready inventory patch is inert; rejected draft patch and full biological output stay outside git.

English commit title: **Move Session import parsing to a bounded browser Worker lifecycle**

English summary: Move whole-object Session read and parsing into a dedicated JS Worker, centralize its lifetime and cancellation, preserve existing Session adoption and settings-only behavior, and record reproducible correctness checks and unresolved admission/performance failures.

---

The following prerequisite-only result is historical. Its old missing-permission conclusion and checkout workflow do not describe the resumed runtime above.

# S05 result — Import Worker prerequisite check

Status: **前提待ち / BLOCKED**. Independent readiness evidence and focused checks are complete; **S05 runtime is not started or complete**. Scope is Issue #597 BUG-02 / BUG-20, S05 only. BUG-01/S02 remain excluded. S06 and later were not started.

## Checkout and authority

Dedicated clone: `/tmp/gbdraw-issue597-S05.ujOWyl`; branch `fix/issue-597-input-session-20260926`, upstream `origin/fix/issue-597-input-session-20260926`. The requested clone/fetch/ff-only pull and explicit dev/main/tag fetch completed with a clean tree. Shared checkout, environments, servers and other sessions' artifacts were not changed.

| Reference | SHA / disposition |
| --- | --- |
| Starting HEAD / fetched target | `b092499d27705092d84037d5acd51032177c619a`; no additional remote commits at intake |
| Latest fetched trusted dev | `88028fd242d263f0fe86aaf9da57b8dc9eb082f6` |
| Fetched main | `4556e04e929a4a85ad28d1833ce7304bd764881c` |
| Product authority dev merge, PR #610 | `af5d942af60353dda199aa487da9152a3576b3fe`; ancestor of HEAD and latest dev |
| S01 | `15bcbca89392cbd63a88fe80dc44b71ad4868061`; HEAD ancestor |
| S03 | `e3684d496df7a3b2a7492443cad50168f6959b36`; HEAD ancestor |
| S04 dev integration | `eed70f05d04f638a5a4ab6c6567982d45e774596`; HEAD ancestor |
| S04 implementation | `b092499d27705092d84037d5acd51032177c619a`; HEAD ancestor |
| Import Worker authority merge | **Absent** in the fetched trusted dev |
| S05 dev integration merge | **None**; the required authority is absent, so no runtime integration was performed |

[S05-readiness.json](../evidence/S05-readiness.json) and its read-only [reproduction script](../evidence/S05_READINESS.py) compare the actual latest dev, historical Product merge and starting HEAD. Human/serialized/authority receipts for PD-OI-044 and PD-OI-045 match in all nine fields each. Historical Contract revision is 21; HEAD/latest dev are revision 22 and byte-identical. Product selection needs no new approval.

The unchanged `tools/inspect_issue597_s01_contracts.py` still **FAILS** its revision-21 assertion. The unchanged historical/pinned `S04_AUTHORITY.py` **PASSES** against its S04 base. Separate current-dev evidence also confirms the receipt/ancestry/byte comparisons; it does not reinterpret either older script or grant permission.

The decisive missing subjects in `origin/dev:tools/web-change-policy.json` are:

| Kind | Exact required subject | Actual trusted dev |
| --- | --- | --- |
| Worker constructor operator | `allowedPrivilegedOwners["Diagram Worker"]` → `services/session-import-client.js` | Absent |
| Codec import edge | `allowedPrivilegedImporters["services/session-file.js"]` → `workers/session-import-worker.js` | Absent |

The active policy is byte-identical between starting HEAD and latest dev. `git apply --check` and the existing `--inert` candidate verifier pass, but neither is runtime authorization. The candidate was applied only to a disposable policy copy. No active authority, checker, mapped contract or permission patch was edited/applied in this checkout.

The trusted Product map contains 14 contract references, recorded with dev file hashes in the readiness JSON. Twelve references have identical HEAD/dev file bytes; the two references to `session-request.test.mjs` share a file changed by the intervening dev track-gap fixes. The diff was inspected and preserved as dev-only work, not copied into this prerequisite-only session. These observations do not certify every mapped checkpoint or replace the prerequisite gate. Any needed mapped-contract update must retain evidence-only → authority-ref-only → runtime delivery.

## Preserved owners and preparation

Selected transport remains **File/Blob → JS Worker read/decompress/fatal UTF-8/JSON.parse → whole-object reply → private untrusted main candidate**. No size-based fallback, bounded sections, ACK, sampling pause or terminal flush is introduced.

| Responsibility | Existing / planned owner and disposition |
| --- | --- |
| File/gzip/UTF-8/limits and streaming Save | Existing `services/session-file.js`; unchanged. One codec, 200 MiB file cap / 512 MiB gzip-expanded cap, magic detection and fatal decoding remain. |
| Import transport/lifecycle | Planned `services/session-import-client.js` and `workers/session-import-worker.js`; neither exists in production. Operation IDs, errors, stale/cancel/crash/teardown, settlement and termination remain unimplemented. |
| Import coordination / validation | Existing `services/config.js::importSessionDocument`; still `readSessionText(file)` → main `JSON.parse(text)` → `assertSafeObjectKeys` → existing preflight/migrations/authority validation. **No main parse path was removed** because its replacement is blocked. |
| Atomic adoption / rollback | Existing config private candidate preparation, Result admission, source-first rollback and transient reconciliation; unchanged. History baseline lifetime checks remain in the existing app/History owners. |
| Availability / Save coordination | Existing state availability and config operation owner; duplicate exact-promise Save join, cross-operation/Generate/reflow reasons, cancellation and browsing remain. No second lock or snapshot owner. |

When authorized, dependency direction remains config → client → Worker → codec; config's Save codec import stays. Worker must not import config, History, request, SVG admission or Python helper owners. Terminate on reply/error settlement before preflight/restore/DOM work; validate unsafe keys before candidate use. Keep immutable adopted payload ownership and rollback checkpoints from [S04 mutation evidence](../evidence/S04_MUTATION_INVENTORY.md).

Ordinary architecture evidence for this delivery: runtime owners, paths and compatibility branches are unchanged; no new production module, dependency, protocol, parser, lifecycle or public default. No superseded runtime owner/path is removed in a prerequisite-only change. The replacement/removal is deferred together until permission is integrated. The read-only evidence script is not connected to application routing or the policy checker. No ratchet exception or review waiver is claimed.

## Verification and fingerprints

Raw outputs are in `/tmp/issue597-S05-evidence-ujOWyl`. [S05-validation.json](../evidence/S05-validation.json) stores environment, source/input/recipe fingerprints, log hashes, commands and scope. Runtime aggregate is `20f51b6c90b46575c7ef20bad43f57e0eca75f644c5f32cb1bd88a05d29cea73` over 372 tracked runtime files, identical to S04. Every individual S04 production/test/input/boundary hash also matches. Prior browser/performance evidence is historical context, not a new S05 PASS.

Python 3.13.3, Node v26.8.2, Ruff 0.15.12, Python Playwright 1.61.0. CLI and Python Playwright are available; checkout-local Node `@playwright/test` is absent. No dependency installation, wheel build, browser server or performance run was needed for this blocked runtime scope. The ordinary sandbox failed at bubblewrap's `/mnt/wslg/distro` mount; identical dedicated-checkout commands ran with sandbox escalation.

| Check | Actual result |
| --- | --- |
| S05 pinned readiness script | **exit 2 / BLOCKED_MISSING_PRIVILEGED_PERMISSION**; receipt and ancestry assertions hold |
| Original revision-21 probe | **exit 1 / FAIL**, unchanged obsolete revision assertion |
| S04 pinned authority evidence | exit 0 / PASS; historical and S04 pinned-base checks only |
| Candidate apply check / `verify_candidate.py privileged "$PWD" --inert` | exit 0 / PASS; shape/detector evidence only |
| Focused Node command in validation manifest | **104 PASS**, 0 fail/skip; existing codec, request, active files, record/discovery metadata, operation consistency, Save lifecycle, settings-only, backing, cache and recovery |
| `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$PWD" pytest tests/test_output_comparison.py::TestOutputComparison -v --basetemp=<task-dir>/pytest-output` | **16 PASS**, 19.95 s; references read-only |
| `ruff check --no-cache gbdraw/` and S05 evidence script | PASS |
| Local trusted-base policy before commit | **exit 0 / Gate PASS / Review REQUIRED**; unchanged checker/authority bytes equal latest dev; not trusted CI or a review waiver |
| Same policy command after commit with `--head HEAD` | Run at handoff; external `policy-after-commit.log`; not trusted CI or merge approval |

The full Python/browser suites and large fixture regeneration were not rerun for unchanged runtime and a blocked transport. S04's full-Python invocation remains **6262 PASS / 6 FAIL / 17 skipped / 11 deselected**, with separately passing corrected reruns; it is not relabeled PASS. No timeout, performance budget or Worker-free assertion was weakened. No LOSAT executable mode was changed in S05.

## Acceptance and remaining measurements

| ID | S05 status |
| --- | --- |
| S-01 | Existing codec/settings/historical admission assertions pass in the focused tests. Worker-backed JSON/gzip/released/current/draft/committed browser acceptance is **not run**. |
| S-02 | Existing operation/lifecycle assertions pass. New import Worker success/error/crash/stale/cancel/teardown/termination/repeat-resource-leak acceptance is **not implemented or measured**. |
| S-03 / S-04 transport | S05 stage/wall, heartbeat, long tasks, main/Worker/process memory, transfer/copy, content equivalence, fresh Load → Generate → CLI/Python replay are **not measured**. |
| S-05 | S01 full-pipeline performance **FAIL remains open**; no new S05 performance PASS. |
| D-01–D-04 | S03 discovery/disclosure/record presentation preserved byte-for-byte; no new S05 browser certification. |
| A-01 | Required authority subjects verified absent; runtime admission blocked. Ordinary unchanged-owner/path evidence recorded. |
| W-01 | Dedicated clone, explicit staging and one evidence/result commit; nonforce same-target push/local-remote equality verified at handoff. |

S01 whole-object transport feasibility stays separate from runtime correctness and full-pipeline performance. The existing Vibrio S04 single Save (5540 ms, heap delta 177459158 bytes, max heartbeat 962.6 ms) is its narrower historical budget result, not the Issue #597 overall target. S05 repeat-import Worker/object URL leaks and native structured-clone wire/copy bytes are **unknown/null**, not zero.

Real large Linear preview Worker-free acceptance remains **FAIL from S01**, not replaced by small Circular success. Observed trigger: fresh preview Load of the real 12-chromosome/full-pairwise gzip fixture. Observed Worker path: `workers/diagram-generation-worker.js`; constructor owner: `services/diagram-generation.js::getWorker` (`new Worker`), initialized by `ensureWorkerInitialized`. Existing helper/feature-extraction requests use `runAuxiliaryWorkerRequest`. The historical observer captures constructor URLs/counts, not the exact originating helper operation/call stack; that causal trigger must still be traced on the real fixture. S05 adds no new measurement or fix. Future probes must distinguish this Python Worker from the proposed JS import Worker while retaining Python count **0** acceptance.

Real fixture expanded JSON is 466767621 bytes, above the 200 MiB plain file cap; its positive case is gzip (134471293 bytes), not plain JSON. Source/fixture recipe and existing measurements remain in [S01 reproduction](../evidence/S01_REPRODUCTION.md) and [S01 transport evidence](../evidence/S01-final-transport.json). Regenerate in an owned environment after authorization; do not access another session's fixture/cache. Keep S01's metric definitions and sequential performance runs.

## Required external step and S05 resumption

A maintainer must deliver the two exact permissions through a **separate permission-only PR from latest dev**, review/merge it into dev and record its merge SHA, following the existing [candidate instructions](../authority-candidates/README.md). Any required mapped evidence/reference updates must precede runtime in their separate deliveries. This request does not authorize creating/pushing/merging that authority PR, editing active policy/checkers or granting blanket Worker permissions. Product outcomes are already approved.

Resume **S05**, not S06: fetch latest dev and the same-named work branch into a new dedicated checkout; read this result; recheck exact permissions, required mapped evidence, ancestry and current receipts; then merge the authorized dev into a clean tree. Record that integration merge separately from the runtime implementation commit. Implement the selected whole-object transport, remove main parse in the same runtime change and complete the requested lifecycle/browser/performance matrix. A candidate/`--inert` PASS or local Gate PASS alone cannot unlock runtime.

S06 remains dependent on completed S05. Its later entrance is the measured Save projection/preflight/restore/SVG admission/DOM mount bottlenecks recorded in S01; no S06 implementation was attempted here. S07 public docs and S08 integration remain later work.

Production, tests, plan docs/evidence and generated diffs were reviewed separately. Only the new S05 readiness script/JSON, validation manifest and this result are staged. Public Gallery, reference SVGs, social preview, vendor/dependency manifests, dist/egg-info, generated wheel and cache-bust are unchanged.

English commit title (requested): **Move Session import parsing to a bounded browser Worker lifecycle**

English summary: Record S05's missing trusted-base Worker permissions and reproducible readiness evidence, preserve approved Product receipts and unchanged runtime, and leave implementation and measurements pending separate authority integration.
