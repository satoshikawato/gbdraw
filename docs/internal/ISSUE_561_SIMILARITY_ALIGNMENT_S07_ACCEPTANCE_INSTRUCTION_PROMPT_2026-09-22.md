# Issue #561 — S07 final acceptance instruction prompt

```text
Complete Session S07 of gbdraw Issue #561: independently audit the finished
Similarity Group alignment implementation against all six accepted Product
Decisions, close integration gaps, run the full required verification, update
public documentation, and prepare the implementation handoff. Do not broaden
scope or weaken tests to obtain a pass.

Repository and branch:
- Use issue-561-similarity-alignment.
- Verify origin/dev contains the accepted authority and S01-S06 are recorded as
  complete in the master-plan ledger.
- Inspect all branch commits and the full diff from origin/dev. Preserve
  unrelated changes and do not create another runtime branch.

Read before editing:
1. AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md
2. the complete master plan and all S01-S06 ledger entries
3. all six accepted authority records
4. PRODUCT_IMPACT_RATCHET.md, ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md,
   WEB_CHANGE_POLICY.md, and SESSION_COMPATIBILITY.md
5. all changed production, test, fixture, and documentation files
6. docs/DOCS.md, docs/CLI_Reference.md, and the existing Web technical owner
   page before choosing documentation locations

Acceptance audit:
1. Trace every Issue #561 acceptance criterion to a Product Decision, production
   owner, and executable test or explicit manual evidence.
2. Confirm exact non-representative reference, missing member, other valid
   member without representative, inparalog selection, direct RBH, multiple
   RBH ambiguity, Skip, and unselected-copy preservation.
3. Confirm Align preserves orientation/Y and Align & orient follows current
   displayed strand, whole-record reverse, unknown-strand position-only, text
   readability, and persistent derived rev indication.
4. Confirm same-row records, negative/positive translations, canvas bounds,
   ruler/definitions/tracks/labels/annotations/comparison edges, and repeated
   Generate/idempotence.
5. Confirm active plan persistence, stable reorder, every invalidation trigger,
   stale repair, Align A/B, Reset, Undo/Redo, Cancel, and failed/canceled/stale
   Result isolation.
6. Confirm current request/Session writer has no legacy string and old supported
   schema access is reader-only and isolated. Search for accidental current-path
   uses, duplicated owners, and stale documentation.
7. Confirm Web, typed Python, and strict CLI share the resolver. Confirm no JS
   ranking implementation and no new LOSATP/group-inference work.
8. Confirm deferred scope remains absent: anchor TSV, Collinear UI, smart
   alignment, synteny propagation, and multi-hop inference.

Documentation:
- Update existing technical documentation and CLI reference rather than adding
  one page per surface. Explain exact reference, candidate priority, both
  actions, Skip, active persistence, invalidation, Reset, Session compatibility,
  strict CLI ambiguity error, and deferred Collinear/TSV scope.
- Add a FAQ entry only if it answers a distinct troubleshooting question.
- Add or change Gallery material only if a finished existing Gallery workflow
  needs it. Do not create a minimal smoke figure as a public showcase.
- Keep documentation examples executable and verify literal commands.

Required verification:
- Run all focused tests listed in master-plan section 13.
- Run:
    ruff check gbdraw/
    python -m pytest tests/ -v -m "not slow"
    node tests/web/architecture-contracts.test.mjs
    node tools/check-web-change-budget.mjs
    python -m build
- Run representative desktop and narrow real-browser journeys. Use Node
  Playwright when available or Python Playwright otherwise. Rerun with required
  escalation if Chromium is blocked by sandboxing.
- Run output-comparison tests read-only. Regenerate tracked references only for
  an intentional reviewed geometry change, inspect the SVG diff, then rerun.
- Inspect production, tests, docs, fixtures, generated artifacts, and full
  architecture diff separately.

Architecture/Product handoff:
- Record concise owner/path evidence and show that the normal implicit resolver
  and string writer were removed while only the bounded legacy reader remains.
- Verify canonical request and current Result admission subjects still satisfy
  every mapped requirement; do not use matching option names as sole evidence.
- If any Architecture or Product Impact exception condition is reached, stop
  and obtain the required independent decision instead of declaring completion.

Finish:
- Update master-plan section 16 with commands, exact results, browser evidence,
  residual risks, and completion status.
- Do not claim complete while a required check, authority, surface, or accepted
  behavior is missing.
- Provide an English proposed commit title and concise summary for the complete
  implementation session. Do not push, create/merge a PR, tag, release, or
  deploy unless separately authorized.
```
