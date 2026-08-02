# Web Gallery coordinate-scale capture plan

Date: 2026-07-31
Status: completed

## Scope

Recapture only these four Linear Scale Style operations after the Web
coordinate-scale control and session projection are stable:

- `lambda_basic_linear/manual-03-03-scale-style-ruler.webp`
- `hepatoplasmataceae_orthogroup/manual-05-01-layout-gc-skew-ruler.webp`
- `hepatoplasmataceae_collinear/manual-05-01-layout-gc-skew-ruler.webp`
- `BGC0000708-BGC0000713/manual-05-02-scale-style-ruler.webp`

The previous crops showed the old Axis & Scale card. The new checkbox moves
**Scale Style** below the axis controls, so the old 54 px top padding no longer
reproduces the same card context.

## Capture contract

Each operation must:

1. Load its own Gallery session.
2. Set and assert `mode=linear`, `form.show_scale=true`, and
   `form.scale_style=ruler`.
3. Open **Axis & Scale** and the **Scale Style** menu.
4. Keep the checked **Show Coordinate Scale** control and selected
   **Ruler (Ticks)** value inside the final crop.
5. Use device scale factor 3 and WebP quality 94.

The final crop uses 145 px top and 130 px right padding. This keeps the card
heading, checkbox, intervening axis controls, and open Scale Style menu
locatable while excluding the unrelated collapsed Labels card.

## Review

Before replacing a WebP, copy the current image to `/tmp` and compare the old
and new images at the same displayed size. Accept a replacement only when:

- the checked checkbox and selected ruler are readable;
- the open menu remains associated with **Scale Style**;
- no unrelated card or floating preview control enters the crop;
- the crop is no larger than the new control layout requires.

## Validation

Run:

```bash
python -m json.tool gbdraw/web/gallery/tutorials/lambda_basic_linear.json
python -m json.tool gbdraw/web/gallery/tutorials/hepatoplasmataceae_orthogroup.json
python -m json.tool gbdraw/web/gallery/tutorials/hepatoplasmataceae_collinear.json
python -m json.tool gbdraw/web/gallery/tutorials/BGC0000708-BGC0000713.json
python -m json.tool gbdraw/web/gallery/tutorials/HmmtDNA_ATskew.json
python tools/capture_gallery_tutorial_screenshots.py --example <example-id> --check --strict
node --check gbdraw/web/gallery/gallery.js
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium
```

Review the affected tutorials at desktop and mobile Gallery widths after the
media replacements.

All four replacements were captured at 990×768. Their full asserted controls
remain inside the crop, and same-width comparison against the preserved old
images confirmed that the new checkbox and changed control order are shown
without unrelated card content.

## Audited without recapture

- `HmmtDNA_ATskew` documents explicit authority in the custom `Ticks` row and
  does not show Axis & Scale. Its existing WebP is retained; its capture
  metadata now names the exact session and asserts the final tick slot.
- `vibrio-harveyi-group-collinear` and `WSSV_genome_comparison` mention scale
  settings in tables, but their operation crops do not include Axis & Scale.
  Their default-visible scales and existing media remain correct.
