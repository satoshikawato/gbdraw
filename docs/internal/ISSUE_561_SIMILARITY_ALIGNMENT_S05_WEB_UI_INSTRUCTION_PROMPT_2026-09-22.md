# Issue #561 — S05 Web UI and accessibility instruction prompt

```text
Implement Session S05 of gbdraw Issue #561: expose the complete Web Similarity
Group alignment journey using the S04 controller, including exact reference
selection, two explicit actions, ambiguity resolution, summary, plan inspection,
candidate preview, and accessibility. Do not duplicate resolver rules in UI.

Repository and branch:
- Use issue-561-similarity-alignment.
- Confirm authority and completed S01-S04 ledger entries.
- Preserve unrelated changes; do not create another runtime branch.

Read before editing:
1. AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md
2. master plan, especially sections 4, 7, 8, and 13
3. accepted similarity-alignment authority
4. S04 controller and tests
5. gbdraw/web/index.html, app-setup.js, orthogroups.js, components.js
6. existing modal, focus-management, toast, feature-popup, drawer,
   highlighting, and accessibility tests

Implement:
1. Feature popup offers distinct `Align` and `Align & orient` actions for the
   exact clicked feature. Labels, descriptions, and titles make the orientation
   effect explicit.
2. Group drawer supplies an exact reference record/feature selector before
   enabling either action. Do not retain the current group-ID-only Align button.
3. Open the selection dialog only for records reported ambiguous by the shared
   resolver. Show feature identifier, coordinates, displayed strand,
   representative/role status, and direct evidence. Do not show score as an
   automatic ranking or imply that representative status is preferred.
4. Each ambiguous record requires Select or Skip. Apply remains disabled until
   all are resolved. Cancel is a complete no-op.
5. Candidate hover/focus uses existing highlight capabilities and restores the
   prior highlight on close. Preview highlighting is not canonical selection.
6. After successful Apply, show a concise live summary with aligned, unchanged,
   explicitly skipped, no-candidate, and reversed counts.
7. Add an active-plan inspector showing exact reference and each record's anchor
   or Skip plus the stored rationale. It is read-only except for the existing
   Reset/clear action.
8. Show a persistent `rev` indicator by deriving effective orientation relative
   to source. Do not persist a separate indicator flag.
9. Use existing local icons/styles and no network dependency or build step.

Accessibility acceptance:
- semantic button and form controls with accessible names;
- modal label/description association;
- focus moves into the dialog and returns to the invoking control;
- Tab/Shift+Tab containment while modal;
- Escape performs Cancel and commits nothing;
- arrow/radio/select behavior follows the chosen native pattern;
- candidate highlighting is not the only indication of selection;
- summary is announced without stealing focus;
- disabled reasons are visible and programmatically associated;
- narrow viewport keeps actions and dialog controls reachable.

Testing:
- popup and drawer exact-reference journeys;
- automatic path has no ambiguity dialog;
- multiple candidate Select and Skip;
- Cancel/Escape no-op;
- keyboard order, focus trap/return, accessible names, live summary;
- candidate highlight cleanup;
- inspector rationale values and rev derivation;
- desktop and narrow real-browser checks with a representative inparalog fixture.

Design constraints:
- UI renders controller state; it does not resolve or mutate the plan directly.
- Keep markup in index.html and focused logic in the existing module layout.
- Do not add a component framework, build step, general wizard, anchor TSV, or
  Collinear UI.

Finish:
- Run focused JS tests and real-browser verification. If Node Playwright is not
  installed, use Python Playwright; if Chromium sandboxing fails, rerun with
  required escalation.
- Capture screenshots only as disposable review evidence unless public docs
  require them in S07.
- Review markup, behavior, accessibility tests, and visual evidence separately.
- Update master-plan section 16 and state whether S06 may start.
- Provide an English proposed commit title and summary; do not push or create a
  PR unless separately authorized.
```
