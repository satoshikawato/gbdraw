# Issue #561 — S04 Web Worker and controller instruction prompt

```text
Implement Session S04 of gbdraw Issue #561: add one Web orchestration owner that
uses the existing diagram Worker to call the shared Python resolver and manages
an immutable alignment draft. Do not finish popup/dialog styling or lifecycle
invalidation in this session.

Repository and branch:
- Use issue-561-similarity-alignment.
- Verify accepted authority and complete S01-S03 ledger entries.
- Preserve unrelated changes; do not create another runtime branch.

Read before editing:
1. AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md
2. master plan and accepted six Product Decisions
3. S01-S03 contracts and tests
4. gbdraw/web/js/app/app-setup.js, orthogroups.js, run-analysis.js
5. gbdraw/web/js/services/diagram-generation.js,
   diagram-worker-protocol.js, session-request.js, history.js,
   history-snapshot.js, feature-identity.js
6. gbdraw/web/js/workers/diagram-generation-worker.js
7. gbdraw/web_support/orthogroup_metadata.py
8. existing helper-operation, Worker lifecycle, orthogroup identity, Generate,
   cancellation, and stale-result tests

Implement:
1. Add one focused createSimilarityAlignmentActions-style module under
   gbdraw/web/js/app/. It owns alignment orchestration, not group browsing,
   candidate ranking, rendering, Session encoding, or generic history.
2. Add one typed helper operation to the existing lazy diagram Worker. It calls
   the S01 Python resolver. Do not create a Worker, a JavaScript resolver, or an
   argv-shaped Python bridge.
3. Build the exact reference from current canonical feature identity. Popup
   actions must reject a click that cannot resolve uniquely. Group-drawer calls
   require a selected exact reference; a group ID is never submitted as the
   reference.
4. Send only validated current group members, direct edges, record/crop/display
   facts, and explicit choices to the helper. Treat helper responses as typed
   resolved/ambiguous/skipped data and validate unknown values before use.
5. Keep unresolved choices in an ephemeral draft. Candidate hover/selection may
   request existing feature highlighting but must not mutate canonical state,
   Result, Session, or history before Apply.
6. Cancel and Worker error discard the draft and retain the current Result and
   committed request. Stale helper completion cannot replace a newer draft or
   action.
7. Apply builds the canonical plan and invokes generation inside one existing
   undoable artifact replacement. Do not set plan state first and call Generate
   afterward. Commit plan/request/Result/history together only on successful
   admission.
8. Reuse compatible committed protein evidence and prove the helper/action does
   not schedule LOSATP or group inference.
9. Replace current popup and drawer functions that write
   selectedOrthogroupAlignmentFeature. Do not leave two active action paths.
10. Keep UI rendering minimal in S04: expose controller state/actions required
   by S05 without adding a second state owner in index.html.

Testing:
- exact popup reference survives into helper input;
- drawer refuses absent exact reference;
- auto-resolved path commits one artifact transaction;
- ambiguous path creates a draft and commits nothing;
- explicit choice/Skip is sent back to the same resolver;
- Cancel, helper error, render failure, cancellation, supersession, and stale
  completion preserve the current Result/request/history;
- one existing Worker is lazily created and reused;
- no LOSATP dispatch;
- no group-ID-only normal action remains.

Design constraints:
- SOLID: controller coordinates; Python resolves; request codec persists;
  renderer renders; history commits.
- KISS: one draft, one plan, one helper operation.
- DRY: no JS candidate-ranking copy.
- YAGNI: no global workflow engine, additional Worker, smart alignment, or
  future Collinear abstraction.

Finish:
- Run focused Web unit tests, Worker protocol tests, run-analysis simple-path
  tests, orthogroup identity tests, and relevant Python resolver tests.
- Run architecture contracts and change-budget checks without weakening them.
- Update master-plan section 16 and state whether S05 may start.
- Provide an English proposed commit title and summary; do not push or create a
  PR unless separately authorized.
```
