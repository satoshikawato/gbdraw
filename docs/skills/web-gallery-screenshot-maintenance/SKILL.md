---
name: web-gallery-screenshot-maintenance
description: Maintain gbdraw Web Gallery tutorial screenshots and manual text. Use when auditing, planning, capturing, replacing, or reviewing screenshots, captions, alt text, tutorial JSON, operation registers, or documentation for `gbdraw/web/gallery`, especially operation crops for dropdowns, segmented controls, file uploads, text/number inputs, track slots, color rules, toolbars, popups, and generated previews.
---

# Web Gallery Screenshot Maintenance

## Overview

Use this skill to keep `gbdraw` Web Gallery tutorials visually truthful: every operation screenshot should show the actual UI the reader operates, and every caption should name the action or selected value.

## Read First

Before editing or reviewing Gallery tutorial screenshots, read:

- `CLAUDE.md`
- `gbdraw/web/CLAUDE.md`
- `docs/WEB_GALLERY_OPERATION_SCREENSHOTS_PLAN.md`
- `docs/WEB_GALLERY_REAL_UI_SCREENSHOT_ROLLOUT_PLAN.md`
- `docs/WEB_GALLERY_DROPDOWN_SCREENSHOT_PLAN.md`
- `docs/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`
- `docs/WEB_GALLERY_SCREENSHOT_RECAPTURE_PLAN.md`, if present

Then inspect the target tutorial JSON and media directory:

- `gbdraw/web/gallery/tutorials/<example-id>.json`
- `gbdraw/web/gallery/media/<example-id>/`

## Audit Workflow

1. Enumerate only media referenced by tutorial JSON first. Treat unreferenced media as cleanup candidates after visible screenshots are fixed.
2. Compare each referenced screenshot with its operation text and caption.
3. Classify each screenshot as:
   - `keep`: real UI/result crop, readable, caption matches.
   - `recrop`: correct state but poor extent, blurry text, missing context, or too much page.
   - `replace`: final preview or stale/empty crop used for an input/edit/upload/setting operation.
   - `add`: operation needs an extra crop to avoid one overloaded screenshot.
4. Record the result in the operation register or the active screenshot plan before replacing files.
5. Keep batches small. Do not migrate every tutorial in one risky change.

## Screenshot Rules

Use the real UI the user operates.

- Mode segmented controls must show the real Circular/Linear control with the active segment readable.
- File uploads must show the uploader and uploaded file rows/chips in the correct order.
- Dropdown/select operations must show the select in context plus a capture-only opened-option overlay with the chosen option highlighted.
- Text and number inputs must show labels, values, units, and nearby related controls. Highlight changed fields when multiple inputs are visible.
- Checkboxes, toggles, and radio groups must show selected state in context.
- Track slots, conservation rings, record labels, record positions, and color rules must show the actual row/panel where values are typed or edited.
- Toolbar actions must show the actual toolbar, not an isolated or reconstructed button.
- Popups must come from real clicks in the restored session whenever possible.
- Generated preview crops are appropriate only for final result checks, visual inspection, legends, popups, or rendered-output comparisons.

Reject preview screenshots for non-final input/edit operations even if the generated result is visually correct.

## Capture Rules

Use Playwright against the current app or Gallery state. If `@playwright/test` is missing, use Python Playwright for focused browser checks. If Chromium fails with a sandbox permission error, rerun the same check with the required sandbox escalation.

Capture standards:

- Use device scale factor 2 or higher; use 3 for dense forms.
- Save temporary PNGs under `/tmp` and commit WebP only.
- Use WebP quality 92-95 for UI controls and dense forms.
- Do not upscale a crop. Recapture at higher scale instead.
- Keep labels and selected values readable at the rendered Gallery size.
- Prefer multiple focused crops over one tall unreadable crop.
- Do not add permanent tutorial-only UI to the app; injected capture overlays must be temporary.

## Writing Rules

Write captions and alt text as action/state descriptions:

- Good: `Select Circular in the mode segmented control.`
- Good: `Set Dinucleotide to GC, Window to 500, and Step to 50.`
- Good: `Enter CDS / product / tyrosine recombinase in the SPECIFIC RULES (-t) row.`
- Bad: `The rendered web app crop shows ...`
- Bad: `Confirm ...`
- Bad: `Web app session view.`

When replacing a screenshot, update the caption and alt text in the same change unless they are already concrete and accurate.

## Verification

Run focused checks after each tutorial or batch:

```bash
python -m json.tool gbdraw/web/gallery/tutorials/<example-id>.json >/tmp/<example-id>.json.check
python tools/capture_gallery_tutorial_screenshots.py --example <example-id> --check
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
```

If Node Playwright is unavailable, use Python Playwright for equivalent media and viewport checks.

Manual review must check desktop and mobile Gallery views, image readability, caption specificity, and that no input/edit operation is represented only by a generated preview.
