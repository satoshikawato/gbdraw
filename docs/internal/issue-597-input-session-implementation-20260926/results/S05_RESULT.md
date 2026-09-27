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
