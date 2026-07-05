---
name: web-gallery-screenshot-maintenance
description: Maintain gbdraw Web Gallery tutorial screenshots, manual text, and structured tutorial content. Use when auditing, planning, capturing, replacing, or reviewing screenshots, captions, alt text, tutorial JSON, operation registers, or documentation for `gbdraw/web/gallery`, especially operation crops for dropdowns, segmented controls, file uploads, text/number inputs, track slots, color rules, tables, toolbars, popups, generated previews, and post-generation editor instructions.
---

# Web Gallery Screenshot Maintenance

## Overview

Use this skill to keep `gbdraw` Web Gallery tutorials visually and textually truthful: every operation screenshot should show the actual UI the reader operates, every caption should name the action or selected value, tutorial text should match the real workflow, and repeated-field content should use structured display.

## Read First

Before editing or reviewing Gallery tutorial screenshots, writing, or structured content, read:

- `CLAUDE.md`
- `gbdraw/web/CLAUDE.md`
- `docs/WEB_GALLERY_OPERATION_SCREENSHOTS_PLAN.md`, if present
- `docs/WEB_GALLERY_REAL_UI_SCREENSHOT_ROLLOUT_PLAN.md`, if present
- `docs/WEB_GALLERY_DROPDOWN_SCREENSHOT_PLAN.md`, if present
- `docs/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`, if present
- `docs/WEB_GALLERY_SCREENSHOT_RECAPTURE_PLAN.md`, if present

Then inspect the target tutorial JSON and media directory:

- `gbdraw/web/gallery/tutorials/<example-id>.json`
- `gbdraw/web/gallery/media/<example-id>/`

## Audit Workflow

1. Audit screenshots, tutorial writing, and structured content together. Do not treat a screenshot pass as complete until the tutorial also satisfies the Writing Rules and Structured Content Rules below.
2. Enumerate only media referenced by tutorial JSON first. Treat unreferenced media as cleanup candidates after visible screenshots are fixed.
3. Compare each referenced screenshot with its operation text and caption.
4. Check adjacent operation screenshots for parent/child duplication. If a broader crop already shows the child controls clearly and the child crop adds no new action context, remove the child operation/media instead of recropping the broader image.
5. Hash or otherwise compare referenced screenshots after the visual pass. Pixel-identical screenshots inside the same tutorial are acceptable only when they are intentionally reused across non-adjacent, clearly different contexts; otherwise consolidate the operation, recrop one image, or delete the duplicate media.
6. Classify each screenshot as:
   - `keep`: real UI/result crop, readable, caption matches.
   - `recrop`: correct state but poor extent, blurry text, missing context, or too much page.
   - `replace`: final preview or stale/empty crop used for an input/edit/upload/setting operation.
   - `add`: operation needs an extra crop to avoid one overloaded screenshot.
7. Review tutorial writing and structure: remove generic Requirements, avoid unnecessary Files tab steps, keep pre-generation setup, generated-result inspection, and post-generation edits separate, and omit redundant operation `title`/`body` when a clear step plus media caption already conveys the action.
8. Review structured content: convert repeated-field setup lists, file mappings, record metadata, track recipes, and color rules to tables when rows share the same fields; keep simple one-column checklists as bullets.
9. For sweeping audits, make a temporary contact sheet of referenced media for each tutorial and inspect the updated contact sheet after fixes.
10. Record screenshot, writing, and structured-content decisions in the operation register or the active screenshot plan before replacing files. If no active register exists, create or update `docs/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`.
11. After tutorial references are fixed, scan `gbdraw/web/gallery/media/` for WebP files unreferenced by tutorial JSON. If tests or docs still refer to an otherwise stale media fixture, update them to use current referenced media before deleting the stale file.
12. Keep batches small. Do not migrate every tutorial in one risky change.

## Screenshot Rules

Use the real UI the user operates.

- Tutorial media must render as appropriately scaled thumbnails in the Gallery, while click-to-open behavior must show the original image itself at natural bitmap size. Do not use a framed page-like modal, captions, headers, or extra white shell around the opened image.
- If a screenshot is unreadable in the Gallery, fix the source capture density, crop extent, or thumbnail/lightbox behavior. Do not rely on browser zoom, CSS upscaling, or a larger wrapper around the same low-resolution bitmap.
- Mode selector screenshots must show the real Circular/Linear control with the active choice readable.
- File uploads must show the real uploader/input card the user operates, not a reconstructed list of file chips. The crop must include the uploaded file and the neighboring per-record controls that are visible in the actual row.
- When the session or command uses per-record metadata that affects the figure, such as `--record_label`, `--record_subtitle`, `definition`, `record_subtitle`, organism/strain, subtitle/title, region, reverse complement, or LOSAT filename, include those fields in both the tutorial table/text and the upload-row crop. Do not split file upload and record label/title entry into disconnected screenshots when the values live in the same input card.
- For tall multi-file upload workflows, keep the full file order and metadata as a tutorial table and use one representative real row crop with the file chip and relevant metadata fields, rather than one unreadable crop of every row.
- Dropdown/select operations must show the select in context plus a capture-only opened-option overlay with the chosen option highlighted.
- Text and number inputs must show labels, values, units, and nearby related controls. Highlight changed fields when multiple inputs are visible.
- When one operation sentence lists several pre-generation settings, include an overview crop that shows the actual text inputs, selects, toggles, or selected chips, even if separate dropdown crops also show individual menu choices.
- Plot title operations must include the full Plot Title input value and the Plot Title Position control in the same readable context. If the default sidebar width clips the title text, widen the sidebar or crop wider for capture instead of accepting a cut-off title input.
- Checkboxes, toggles, and radio groups must show selected state in context.
- Track slots, conservation rings, record labels, record positions, and color rules must show the actual row/panel where values are typed or edited.
- Feature-specific color rule operations must show the real `SPECIFIC RULES (-t)` controls, not a generated plot preview.
- When several color rules are entered before one final generation, use one focused crop showing all entered rule rows if it is readable; otherwise use focused per-row input crops. Do not place a full-plot screenshot after each rule.
- If a control screenshot looks too wide because rows have large unused horizontal gaps, treat the capture source as wrong. Compare sibling operation images in the same example, then recapture the real UI with a compact sidebar/panel width that matches their density. Do not fix this by changing Gallery display CSS or by cropping only the already-wide bitmap.
- Show the generated plot only once after all required pre-generation inputs and rules have been entered, unless a later step is explicitly about inspecting a popup or visual result.
- Toolbar actions must show the actual toolbar, not an isolated or reconstructed button.
- Crops with highlighted toolbar buttons, track slots, drawer rows, or other position-sensitive controls must include enough surrounding controls above, below, left, or right for the highlight to convey where the target sits. Do not crop exactly to only the highlighted box when neighboring items define the action context.
- Popups must come from real clicks in the restored session whenever possible.
- Popup crops must be tight enough that the highlighted clicked feature, match ribbon, orthogroup ribbon, or collinear block and the popup text are readable in the rendered Gallery. Do not use a full generated-preview crop when the inspected target is only a small part of the figure.
- If the popup covers the highlighted feature, orthogroup ribbon, or collinear block, move the popup during capture and then crop the real UI state. Do not accept a crop where the popup is readable but the highlighted diagram target is hidden.
- Before recapturing a popup from a restored session, verify the saved SVG and session metadata agree with the tutorial's intended label rules. Stale session artifacts can make popups fall back to `product` even when the rendered figure or tutorial expects `gene`; fix the session/result metadata first, then recapture.
- Feature editor popup crops must keep the relevant affected item(s) or multi-target controls readable, including apply-to-all/scope icons or dialogs when they are the point of the operation. Crop or move the popup so these controls are not cut off at the image edge.
- Judge readability at the rendered Gallery size, not at the source bitmap's full size. A source image that looks readable only when opened standalone is a `recrop` if Gallery CSS downscales it enough to make popup text, highlighted targets, or affected item(s) hard to read.
- Generated preview crops are appropriate only for final result checks, visual inspection, legends, popups, or rendered-output comparisons.
- Generated preview screenshots with floating controls must leave enough capture-time breathing room around the SVG content so zoom/original/search controls do not cover top records, labels, or legends, and bottom titles do not visually collide with legends. If the visible crop is cramped, recapture with a larger viewport or SVG canvas instead of accepting overlap.
- Do not show the same generated preview or popup twice in immediate succession as both step-level media and operation media. If the operation already carries the result or popup crop, omit step-level media or make the two crops visibly different and purposeful.
- More generally, omit step-level media when the operation media already shows the same UI state, Files tab, generated preview, or popup without adding new context.
- When one operation crop already includes another operation's entire control area, keep the broader crop if it remains readable and remove the narrower duplicate. Prefer a single truthful panel crop over consecutive screenshots that repeat the same controls.
- Post-generation editor screenshots must come from the exact restored session for that example. Before capturing drawers such as Legend, Features, or Orthogroups, verify the restored editor state matches the example-specific generated result. A generic or stale drawer state, such as a BGC legend editor showing only `CDS`, is a `replace`, even if it is a real drawer crop.
- Drawer screenshots must show the named active drawer tab and must not be overlapped by a popup or another drawer's controls. Legend drawer crops should be scrolled so the entries named by the caption are visible.
- If a restored session's saved editor state and rendered SVG disagree, fix the app restore behavior or refresh the session artifact before capturing. Do not document the broken intermediate state as the tutorial screenshot.
- If a restored session's UI controls disagree with the tutorial command or rendered output, fix the app restore behavior or refresh the session artifact before capturing. This includes settings that may exist only in `cliInvocation.args` in older sessions, such as multi-record `--multi_record_position` tokens.
- When capturing a drawer, keep the actual drawer controls visible, but exclude unrelated floating preview controls if they visually overlap the drawer; use temporary capture-only CSS rather than permanent app changes.

Reject preview screenshots for non-final input/edit operations even if the generated result is visually correct.

## Capture Rules

Use Playwright against the current app or Gallery state. If `@playwright/test` is missing, use Python Playwright for focused browser checks. If Chromium fails with a sandbox permission error, rerun the same check with the required sandbox escalation.

Capture standards:

- For restored-session captures, assert critical control values in the browser before saving the crop when the operation depends on exact row/order/selector state. Multi-record Record Order crops, for example, should verify `#1@1, #2@1, #3@2, #4@2, #5@2, #6@2` rather than trusting the restored page by sight.
- For restored-session label or popup captures, assert representative label text in both `results[0].content` and the relevant feature/orthogroup metadata before saving. For qualifier-priority examples, check `qualifierPriorityRules`, feature `label`/`display_label`, and orthogroup overrides or candidates instead of trusting a visual crop alone.
- For representative upload-row crops, assert both the file order and the row's metadata values in app state before capture; do not trust visible labels alone because input values may be clipped in the DOM text.
- Use device scale factor 2 or higher; use 3 for dense forms and popup/detail screenshots. Treat requests such as "150 dpi" as a demand for a higher-density source bitmap: increase device scale factor or crop tighter, never upscale an existing low-resolution image.
- When the request is only to improve resolution, preserve the same semantic crop. Increase device scale factor, capture density, or source viewport as needed, but do not include extra sibling slots, following sections, or unrelated context. For stack/list panels, compute the clip bottom from the same last visible row rather than from the whole container.
- Save temporary PNGs under `/tmp` and commit WebP only.
- Use WebP quality 92-95 for UI controls and dense forms.
- Do not upscale a crop. Recapture at higher scale instead.
- Keep labels and selected values readable at the rendered Gallery size.
- For multi-row form captures such as `SPECIFIC RULES (-t)`, preserve the intended vertical scope when recapturing. If the previous screenshot included the surrounding uploader rows and whitespace, keep that top/bottom extent, but reduce excess horizontal whitespace by changing the source panel width before capture.
- Use temporary capture-only CSS or DOM changes to hide unrelated fixed bars, footers, or following cards only when they intrude into the preserved crop area. Do not commit those UI changes unless the app itself is wrong.
- For highlighted track-slot rows or toolbar buttons, prefer a selector crop with asymmetric padding or a broader panel crop so adjacent relevant items stay readable.
- When overlapping SVG paths intercept a real popup target, use a forced click only for the capture action and then visually verify the opened popup corresponds to the intended feature, orthogroup, ribbon, or block.
- When a replacement screenshot must be reviewed immediately in a running Gallery tab or hosted static cache, change the media filename or add a deliberate cache-bust path change; do not rely on same-URL image replacement being visible after a normal refresh.
- Prefer multiple focused crops over one tall unreadable crop.
- Do not add permanent tutorial-only UI to the app; injected capture overlays must be temporary.

For broad "all screenshots" refreshes, distinguish "all capture-defined operations" from all referenced media. If a tutorial still has manual-only screenshots, either add capture metadata for high-risk controls/popups before recapturing or explicitly report which referenced media remain manual.

## Writing Rules

Keep tutorial text aligned with the actual workflow.

- Use user-facing labels in tutorial text. Avoid implementation/UI-pattern jargon such as `segmented control`; write the action the reader should take, such as `Select Circular`.
- Omit generic Requirements sections when they only restate that the reader needs a browser or the listed input file.
- Do not add a Files tab step unless checking Gallery artifacts is part of the workflow. Input filenames can usually live in `downloads` or the upload step.
- Do not add Color rules or Post-generation edits sections when the example does not require color-rule setup or real post-generation editing.
- Prefer screenshots over dense setup lists when the UI can show the values clearly. Keep any remaining text to the values the reader must type or choose.
- Do not describe automatically generated output as a manual editor task. If `Generate Diagram` creates legend entries, tracks, labels, or previews, say they are generated and reserve drawers/editors for review or optional tweaks.
- Use action text for operations and state text for generated results. Avoid vague instructions such as `keep visible` when the UI already produced the state.
- Omit operation `title` and `body` when a single media item already sits under a clear step title/body and the extra operation text would only create a redundant bold subheading.
- Keep pre-generation setup, generated-result inspection, and post-generation edits as separate concepts.
- When a step already has a `table`, do not repeat the same fields as slash-delimited operation text or captions. Refer to the table row or category instead.

Write captions and alt text as action/state descriptions:

- Good: `Select Circular.`
- Good: `Set Dinucleotide to GC, Window to 500, and Step to 50.`
- Good: `Enter CDS / product / tyrosine recombinase in the SPECIFIC RULES (-t) row.`
- Bad: `The rendered web app crop shows ...`
- Bad: `Confirm ...`
- Bad: `Web app session view.`

When replacing a screenshot, update the caption and alt text in the same change unless they are already concrete and accurate.

## Structured Content Rules

Use structured display when tutorial content has repeated fields.

- Prefer a table over slash-delimited bullets for color rules, track-slot recipes, file mappings, record metadata, or any list where each row shares the same fields.
- When one slash-delimited or comma-packed structured list is found, scan sibling tutorial JSON for the same pattern before stopping; convert true repeated-field content and leave simple one-column checklists as bullets.
- Use clear column headers that match the UI or data model, such as `Feature`, `Qualifier`, `Value pattern`, `Color`, and `Legend caption`.
- Keep caveats such as regex specificity in a separate `note`; do not hide them inside a dense row.
- If the Gallery renderer only supports bullet lists, add a small reusable renderer field such as `table: { columns, rows }` rather than encoding table markup in JSON strings.
- Add compact CSS with horizontal overflow for narrow viewports and update focused browser tests to assert table headers and representative rows.
- Keep existing `items` bullets for simple one-column checklists where a table would add noise.

## Verification

Run focused checks after each tutorial or batch:

```bash
python -m json.tool gbdraw/web/gallery/tutorials/<example-id>.json >/tmp/<example-id>.json.check
python tools/capture_gallery_tutorial_screenshots.py --example <example-id> --check
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
```

If Node Playwright is unavailable, use Python Playwright for equivalent media and viewport checks.

When adding or removing tutorial media references, update focused browser tests that assert image counts or exact `src` values.
When removing a duplicate operation media reference, delete the now-unreferenced WebP and assert the removed `src` is absent if a focused test already covers that tutorial.

Manual review must check desktop and mobile Gallery views, image readability, caption specificity, writing/section structure, table suitability, preview toolbar/title/legend overlap, and that no input/edit operation is represented only by a generated preview.
