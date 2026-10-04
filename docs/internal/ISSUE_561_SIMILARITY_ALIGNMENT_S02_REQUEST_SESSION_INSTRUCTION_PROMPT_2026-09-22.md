# Issue #561 — S02 request and Session instruction prompt

```text
Implement Session S02 of gbdraw Issue #561: make the resolved alignment plan and
stable per-record X/Y translations part of the one canonical typed request and
Session path, with a bounded legacy reader. Do not implement final renderer
geometry, CLI behavior, or Web interaction in this session.

Repository and branch:
- Use issue-561-similarity-alignment; do not create another runtime branch.
- Verify origin/dev contains the six accepted similarity-alignment Product
  Decisions and that S01 is complete in the master-plan ledger.
- Inspect git status and preserve unrelated changes. Do not commit to dev/main.

Read before editing:
1. AGENTS.md, CLAUDE.md, and gbdraw/web/CLAUDE.md
2. docs/internal/ISSUE_561_SIMILARITY_ALIGNMENT_MASTER_PLAN_2026-09-22.md
3. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md
4. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
5. docs/internal/PRODUCT_IMPACT_RATCHET.md
6. docs/SESSION_COMPATIBILITY.md
7. S01 implementation and tests
8. gbdraw/api/requests.py, gbdraw/api/options.py,
   gbdraw/session_request_codec.py, gbdraw/session_io.py
9. gbdraw/web/js/services/session-request.js,
   gbdraw/web/js/services/config.js,
   gbdraw/web/js/services/history-snapshot.js, and state.js
10. relevant request/session/documentation contract tests and released fixtures

Required architecture:
- One current typed request owns the SimilarityAlignmentPlan.
- One Linear layout value owns finite X/Y base translations keyed by recordKey.
- RecordPresentation remains the base orientation owner; do not add a second
  source-orientation boolean.
- The Session projects the same canonical plan and translations. Do not add a
  top-level Session-only copy or a Web-only duplicate.
- Current writers never emit align_orthogroup_feature or
  alignOrthogroupFeature.
- Supported old schemas retain a named reader-only adapter for their string.

Implement:
1. Re-audit current Session and request versions. Allocate the next current
   versions only if required and available. Update all Python, Web, docs, and
   test version authorities atomically. Do not retain a branch-only intermediate
   reader or migration.
2. Add strict codec support for the typed plan and record translations using
   current naming conventions. Reject unknown fields, duplicate record keys,
   non-finite translations, invalid identity combinations, and inconsistent
   record coverage.
3. Keep translation scope deliberately small: x and y only; no scale, rotation,
   generic matrix, constraints, animation, or Circular extension.
4. Replace the current writer projection from protein comparison settings. The
   alignment is display/layout state, not LOSATP pipeline identity.
5. Preserve supported legacy request schemas. Decode their string into an
   explicit private legacy-alignment representation and route it only to the
   later compatibility materializer. Current request construction and typed API
   must not accept that private type.
6. Add a representative positive fixture proving that the legacy contract was
   released or present in first-parent main history. Add malformed and
   ambiguous legacy negatives.
7. Ensure load of a saved preview does not construct the Worker. Validation or
   materialization can remain deferred to the first operation that needs Python.
8. Audit save-before-materialization. If a valid released legacy Session lacks
   enough metadata to write a resolved current plan, stop this edge case and
   prepare a narrow Product Decision Pack. Do not silently drop intent, keep a
   current legacy writer, or invent a save block.
9. Remove superseded current writer/state string owners in the same change only
   when all current projections and tests use the new plan. Keep the legacy
   decoder visibly isolated by schema/version.

Testing:
- Python and Web current request round trip;
- current Session save/fresh Load/regenerate round trip;
- finite keyed X/Y translations, reordered records, duplicate keys, unknown
  fields, and invalid plan combinations;
- supported legacy request/session reads and reproduces a private legacy value;
- current writer contains no legacy alignment field;
- future/unsupported schemas fail;
- load-only saved preview remains Worker-lazy;
- request equality and canonical current-state effects remain intact.

Design constraints:
- SOLID: codecs validate and translate; they do not resolve anchors or render.
- KISS: one current writer, one plan schema, one bounded old-schema reader.
- DRY: Session and Generate project the same canonical request.
- YAGNI: no downgrade writer, no parallel JSON store, no generalized transform
  format, and no support for branch-only version numbers.

Finish:
- Run focused session_request_codec, session_io, Web session-request, lazy-load,
  documentation-contract, and architecture-contract tests.
- Run the Web change-budget check against the appropriate base and record its
  output; do not weaken its acceptance criteria.
- Review production, migration fixtures, tests, and docs separately.
- Update master-plan section 16 and state whether S03 may start.
- Provide an English proposed commit title and summary; do not push or create a
  PR unless separately authorized.
```
