# Session 03 instruction prompt: single Align action and review palette

## Mission

Deliver the visible Issue #586 Web workflow on top of the completed controller:
one `Align…` entry point, immediate loading feedback, and one accessible
floating review palette that shows all records, recommendations, candidates,
per-record orientation intent, and effective outcomes. Preserve existing canvas
pan/zoom and candidate picking.

## Branch and prerequisites

Use only `issue-586-similarity-alignment-ux-20260924`. Fetch, inspect status,
verify Sessions 00–02 are ancestors, and preserve unrelated changes.

When this prompt is supplied as the session request, it authorizes one focused
Session 03 commit and push to the same-named branch. It does not authorize a PR,
merge, release, or direct push to `dev`/`main`.

## Read before editing

Read `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`, the master plan, active
Product decisions, architecture/product ratchets, Sessions 01–02 diffs and
ledger evidence, the current HTML/CSS/templates, feature-popup action owner,
right drawer, similarity controller, and browser/unit tests.

Keep the SPA build-free. Keep `create*` entry points in top-level `app/*.js` and
split a focused private view helper only if the existing module would otherwise
mix independent responsibilities. Do not add a framework or dependency.

## Entry-point UI

- Replace the separate Web `Align` and `Align & orient` actions with one
  `Align…` action in both the feature popup and Similarity Groups drawer.
- Both surfaces pass the exact selected record/feature identity to the same
  controller operation.
- On activation, disable every conflicting alignment start control and display
  `Resolving…` (or equally clear text/spinner) immediately. Expose `aria-busy`
  and an appropriate live status without creating repeated announcements.
- Restore focus deliberately on error/Cancel; when review opens, move focus to
  its heading or first meaningful control.

## Review palette UI

The existing nonmodal/floating palette remains the review surface. Do not bring
back a backdrop or focus-trapped modal.

Render:

- exact reference biological name/ID, record, coordinates, and strand;
- one stable displayed-record row for every non-reference record;
- unchanged status for missing/unusable records;
- selected anchor plus a `Recommended` indication and plain-language reason;
- every eligible candidate with biological name or feature ID as the primary
  label, coordinates, strand, representative status, direct-evidence detail,
  and internal ID only as secondary detail;
- Select controls and Skip for each applicable row;
- `Match reference direction` per applicable record, default off;
- the effective outcome before Apply: orientation preserved, whole record will
  reverse, or orientation preserved because strand is unknown;
- an immediately enabled Apply button when the Python-provided default draft is
  valid, plus Cancel and retryable error/status text.

Do not claim scientific confidence or superiority. `unique representative` and
`deterministic candidate 1` are convenience recommendation bases. A changed
user selection should remain visibly selected without falsely retaining a
`Recommended` badge on the replacement.

Unknown-strand handling must be truthful: a match request cannot display a
future reversal when either required displayed strand is unknown. Follow the
domain projection; do not recalculate strand semantics in the template.

## Canvas and responsive behavior

- Preserve manual canvas pan/zoom while the palette is open.
- Candidate number/guide clicks and row controls dispatch the same controller
  draft mutation.
- A candidate without a unique drawable location remains selectable from the
  list and does not receive a misleading marker.
- Preview-only guides/labels never enter Result SVG, downloads, Session, or
  History.
- Palette drag/reposition remains available where currently supported.
- At 390 CSS pixels, every record, Select/Skip control, orientation control,
  reason, Apply, Cancel, and error is reachable without horizontal page loss.
- Keyboard users can traverse and operate all choices; visible focus remains
  clear; radio grouping and accessible names include record/candidate context.

## Tests and visual review

Update focused DOM/unit tests and Playwright coverage. Required browser journeys:

1. popup exact non-representative reference -> immediate busy state -> review;
2. drawer exact reference -> same review path;
3. uniquely resolved targets still open an immediately applicable review;
4. ambiguous targets display recommendations/reasons, replacement, and Skip;
5. independent per-record orientation outcomes, including unknown strand;
6. canvas selection parity with row selection while pan/zoom remains usable;
7. Apply failure keeps choices and supports correction/retry;
8. Cancel and stale completion preserve the prior Result/History;
9. keyboard-only operation and 390-pixel viewport.

Run at minimum:

```bash
node --test \
  tests/web/similarity-alignment-actions.test.mjs \
  tests/web/right-drawer.test.mjs

npx playwright test \
  tests/web/similarity-alignment-ui.playwright.spec.js \
  --project=chromium --workers=1

node tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
git diff --check
```

Prepare the browser wheel first if required by the harness. Visually inspect at
desktop and 390-pixel widths at readable scale. Store only evidence expected by
the repository; do not replace `examples/gbdraw_social_preview.png`.

## Architecture evidence, commit, and handoff

Record the exact removed actions, common controller path, accessibility/browser
results, and any justified private view decomposition in the ledger. Review
production, test, and visual/generated diffs separately.

Use an English commit title such as:

```text
Unify the similarity alignment review UI
```

Push the branch and report commit, focused/browser results, visual findings,
and lifecycle/artifact work remaining for Session 04.
