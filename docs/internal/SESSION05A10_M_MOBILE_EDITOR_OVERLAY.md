# 05A4-06: Mobile Editor overlay accessibility

Starting `origin/dev`: `c9ba546a1a83d51f231853291c4944ca5228f860`.
Work was isolated from the user's dirty checkout. Only 05A4-06 is in scope.

## Original reproduction

Before editing CSS, load the retained J43 `HmmtDNA_ATskew.gbdraw-session.json`
or J44 `lambda_basic_linear.gbdraw-session.json` at exactly 390 × 844.
Use the retained journey's `scrollIntoView({block: 'center'})` on the toggle.
Both sessions reproduce the same geometry. Keyboard Enter was used only to
measure the original open state because pointer activation was blocked.

Coordinates below are `(left, top, width, height)` in CSS pixels, relative to
the Result Preview rectangle; page scroll does not affect this comparison.

| State | Preview size | Toggle | Preview controls | Drawer | Toggle center |
| --- | --- | --- | --- | --- | --- |
| Before, closed | 374 × 416 | (337, 274, 35, 60) | (10, 306, 354, 100) | (372, 2, 360, 412), hidden | Feature-search status intercepts |
| Before, open | 374 × 416 | (-23, 274, 35, 60) | (10, 306, 354, 100) | (12, 2, 360, 412) | Outside Preview |
| After, closed | 374 × 477.5 | (336, 10, 36, 60) | (10, 367.5, 354, 100) | (372, 2, 334, 473.5), hidden | Toggle descendant |
| After, open | 374 × 477.5 | (2, 10, 36, 60) | (38, 367.5, 354, 100), behind drawer | (38, 2, 334, 473.5) | Toggle descendant |

The original closed toggle/toolbar intersection is **27 × 28 = 756 px²**.
The original open intersection is **2 × 28 = 56 px²**. Both become zero-area
intersections. The new original-case test also fails against the unchanged
starting HTML because the toggle center is blocked.

## Cause and correction

Several independent absolute overlays shared the same small Preview area:

- The two-row toolbar grew into the toggle's fixed `bottom: 5rem` reservation.
- The search panel, at z-index 30, intercepted the z-index 20 toggle center.
- Translating the toggle by a 360 px drawer width exceeded the available
  Preview width. The toolbar at z-index 20 also painted over the drawer at 10.
- Hidden overflow could scroll horizontally during focus/scroll-into-view,
  shifting the open drawer and toggle out of their intended bounds.

`index.html` continues to own all layout. A container query uses the actual
Preview pane width (up to 40rem), where the full toolbar needs to wrap, so it
also covers the 844 × 390 landscape layout. Search occupies a scrollable grid
row above the naturally sized toolbar. The toggle has a 2.25rem side lane at
the top; drawer width is capped by both 360 px and the remaining Preview width.
A local stacking context puts the narrow Preview overlays behind the open
drawer. Their complete hit targets move behind its edge, without exposed
button fragments. `overflow: clip` prevents programmatic horizontal scrolling
of that region. The canvas-padding panel shares the search row without
covering the persistent toolbar.

Closed-drawer controls are all retained. The nearly full-width open drawer
obscures its underlying toolbar; closing it restores access. This is the
open-drawer exception specified in the task, not a new collapsed toolbar.
The wider existing mobile rules still apply outside the compact container;
the desktop overlay oracle is unchanged. Existing safe-area handling stays
with the fixed Generate bar; no new platform-specific inset is needed.

## Verification

- J43 and J44 at 390 × 844, 375 × 667, 430 × 932, and 844 × 390: closed/open
  geometry, positive-area overlap, hit testing, real controls, Editor list
  scrolling, Enter/Space/Escape, title and ARIA preservation, mounted SVG
  feature geometry/text, selected Result identity, and zero browser errors or
  external requests. JSON geometry and screenshots are written to test output.
- No Result, Circular → Linear → Circular, real GenBank source replacement,
  desktop → mobile → desktop, exactly one drawer/toggle, and no inline toggle
  geometry: one additional browser case.
- Unchanged Issue #461 matrix: 1024 × 768, 1366 × 768, 1920 × 1080, and
  2541 × 1409, with controls operable in both drawer states.
- Neighboring browser suites: `mode-transition-editor-state`,
  `mode-transition-result`, `multipart-label-regeneration`,
  `composite-placement-regeneration`, `source-visibility-reconciliation`, and
  `source-legend-reconciliation`, plus the rest of `right-drawer`.
- Drawer and visibility unit tests; focused offline packaging tests; built
  wheel browser verification at 390 × 844 and 1366 × 768, external networking
  blocked, including cold and repeated generation.
- PR smoke remains exactly **10 cases in 7 files**. New cases are selected by
  full functional CI. No workflow, routing, or timeout policy changes.

The test scrolls Preview into the usable space between the existing sticky
header and fixed Generate bar before hit testing. This is test navigation,
not a runtime geometry watcher. Reset layout's intentional layout/viewport
changes are allowed; biological feature geometry and SVG text are preserved.

Local detailed evidence is retained under
`docs/internal/gbdraw_v014_session05a10_m_evidence_2026-09-12/` in the user's
original checkout (outside this commit): original closed/open measurements,
screenshots, execution logs, final geometry, and packaged verification.
The tests regenerate the after screenshots and measurements in their output
directory. No Gallery or reference screenshot is updated.

## Product and architecture

Product preflight: **IMPLEMENT_EXISTING_AUTHORITY**. The task's explicit
operable-visible-toggle contract and `gbdraw/web/CLAUDE.md`'s reactive
availability invariants select the outcome. The available actions, close and
Escape transitions, Result authority, and saved intent remain intact. There
is no unresolved alternative outcome, missing authority, or prohibited
change; EVIDENCE_REQUIRED, PRODUCT_DECISION_REQUIRED, and NOT_ALLOWED do not
apply. No new BD record or candidate authority is introduced.

Architecture: ordinary non-increasing layout correction; this is not
architecture-bearing. Layout remains in `index.html`; `right-drawer.js` still
owns every visibility/tab transition. No semantic owner, canonical path,
compatibility path, state, watcher, dependency, or resource signal is added.
The geometry oracle remains in `right-drawer.playwright.spec.js`. Schemas,
History, persistence, and scientific rendering are unchanged. Rollback is a
revert of this single implementation commit.

Other 05A4 findings are unchanged. The complete adversarial rerun, S11
reacceptance, S12, and baseline/publication-candidate assignment are outside
this change. Closure requires the exact merged dev checks and staging gates;
local passing evidence alone does not establish that acceptance.
