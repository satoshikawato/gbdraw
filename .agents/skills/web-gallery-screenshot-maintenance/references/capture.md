# Gallery capture requirements

Read for screenshot audits, recaptures, or capture-metadata changes. These
requirements preserve example identity and the documented operation.

## Screenshot Rules

Use the real UI the user operates.

- Treat a good existing screenshot as the semantic-crop baseline. Prefer the smallest truthful crop that keeps the operated control, selected value, and only the neighboring context needed to locate it. A larger crop is not an improvement merely because it contains more real UI.
- Tutorial media must render as appropriately scaled thumbnails in the Gallery, while click-to-open behavior must show the original image itself at natural bitmap size. Do not use a framed page-like modal, captions, headers, or extra white shell around the opened image.
- If a screenshot is unreadable in the Gallery, fix the source capture density, crop extent, or thumbnail/lightbox behavior. Do not rely on browser zoom, CSS upscaling, or a larger wrapper around the same low-resolution bitmap.
- Mode selector screenshots must show the real Circular/Linear control with the active choice readable.
- Do not reuse another example's data-dependent crop. A shared mode selector
  or equally data-independent control may be reused only when the operation
  declares `genericMedia: true` and the crop contains no example-specific
  file, metadata, setting, or result.
- File uploads must show the real uploader/input card the user operates, not a reconstructed list of file chips. The crop must include the uploaded file and the neighboring per-record controls that are visible in the actual row.
- When the session or command uses per-record metadata that affects the figure, such as `--record_label`, `--record_subtitle`, `definition`, `record_subtitle`, organism/strain, subtitle/title, region, reverse complement, or LOSAT filename, include those fields in both the tutorial table/text and the upload-row crop. Do not split file upload and record label/title entry into disconnected screenshots when the values live in the same input card.
- For tall multi-file upload workflows, keep the full file order and metadata as a tutorial table and use one representative real row crop with the file chip and relevant metadata fields, rather than one unreadable crop of every row.
- Dropdown/select operations must show the select in context plus a capture-only opened-option overlay with the chosen option highlighted.
- Text and number inputs must show labels, values, units, and nearby related controls. Highlight changed fields when multiple inputs are visible.
- When one operation sentence lists several pre-generation settings, include an overview crop that shows the actual text inputs, selects, toggles, or selected chips, even if separate dropdown crops also show individual menu choices.
- Plot title operations must keep the Plot Title input and Plot Title Position control in the same readable context. If a long title cannot fit at the normal compact sidebar width, put the exact copyable value in tutorial text or a table and let the screenshot identify the real input without implying that its visible substring is complete. Widen only enough to make the control recognizable; do not turn a compact operation crop into a page-wide image solely to expose every character of a long title.
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
- When one operation crop already includes another operation's entire control area, keep the broader crop only if it remains compact at Gallery size and every extra area helps explain the action. Otherwise keep the narrower crop or recrop the shared panel. Prefer a single truthful panel crop over consecutive screenshots that repeat the same controls, but do not trade away compactness merely to reduce the image count.
- Post-generation editor screenshots must come from the exact restored session for that example. Before capturing drawers such as Legend, Features, or Orthogroups, verify the restored editor state matches the example-specific generated result. A generic or stale drawer state, such as a BGC legend editor showing only `CDS`, is a `replace`, even if it is a real drawer crop.
- Drawer screenshots must show the named active drawer tab and must not be overlapped by a popup or another drawer's controls. Legend drawer crops should be scrolled so the entries named by the caption are visible.
- If a restored session's saved editor state and rendered SVG disagree, fix the app restore behavior or refresh the session artifact before capturing. Do not document the broken intermediate state as the tutorial screenshot.
- If a restored session's UI controls disagree with the tutorial command or rendered output, fix the app restore behavior or refresh the session artifact before capturing. This includes settings that may exist only in `cliInvocation.args` in older sessions, such as multi-record `--multi_record_position` tokens.
- Treat capture metadata as executable documentation even when the current bitmap is retained. Its active tab, selected values, target selector, and crop composition must reproduce the visible screenshot state; do not leave a recipe that would silently replace a Details image with Qualifiers or otherwise change the documented action on the next refresh.
- For ordered or movable controls such as track slots, comparison series, and record rows, replay the tutorial steps against the reset/default state and compare the resulting order with the saved Gallery session. Document every required move explicitly instead of relying on a pre-arranged restored session.
- When capturing a drawer, keep the actual drawer controls visible, but exclude unrelated floating preview controls if they visually overlap the drawer; use temporary capture-only CSS rather than permanent app changes.

Reject preview screenshots for non-final input/edit operations even if the generated result is visually correct.

## Capture Rules

Use Playwright against the current app or Gallery state. If `@playwright/test` is missing, use Python Playwright for focused browser checks. If Chromium fails with a sandbox permission error, rerun the same check with the required sandbox escalation.

Capture standards:

- Treat capture metadata as a declarative test contract. A
  `dataDependent: true` operation must use the example's explicit
  `capture.session`, declare exact `capture.assertAppState` paths, and list the
  expected `capture.visibleControls` and/or `capture.visibleText`. The capture
  tool must verify that those controls, selected values, checked states, and
  identity text are inside the final crop before writing the WebP.
- For restored-session captures, assert critical control values in the browser before saving the crop when the operation depends on exact row/order/selector state. Multi-record Record Order crops, for example, should verify `#1@1, #2@1, #3@2, #4@2, #5@2, #6@2` rather than trusting the restored page by sight.
- For restored-session label or popup captures, assert representative label text in both `results[0].content` and the relevant feature/orthogroup metadata before saving. For qualifier-priority examples, check `qualifierPriorityRules`, feature `label`/`display_label`, and orthogroup overrides or candidates instead of trusting a visual crop alone.
- For representative upload-row crops, assert both the file order and the row's metadata values in app state before capture; do not trust visible labels alone because input values may be clipped in the DOM text.
- Use device scale factor 2 or higher; use 3 for dense forms and popup/detail screenshots. Treat requests such as "150 dpi" as a demand for a higher-density source bitmap: increase device scale factor or crop tighter, never upscale an existing low-resolution image.
- When the request is only to improve resolution, preserve the same semantic crop. Increase device scale factor, capture density, or source viewport as needed, but do not include extra sibling slots, following sections, or unrelated context. For stack/list panels, compute the clip bottom from the same last visible row rather than from the whole container.
- Before overwriting an existing tutorial image, preserve or render the previous image for a same-size side-by-side review. Do not accept the new capture until the action remains at least as easy to identify, the crop is no less compact without a documented reason, and any clicked feature or ribbon remains visibly associated with its popup.
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
