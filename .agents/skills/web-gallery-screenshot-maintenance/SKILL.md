---
name: web-gallery-screenshot-maintenance
description: Maintain or audit tutorial JSON, text, and screenshots in gbdraw/web/gallery. Does not cover general Web runtime changes or all public documentation.
---

# Maintain Gallery tutorials

Keep each tutorial faithful to the actual workflow: show the operated control,
use the example's own data and settings, and make captions name the action or
result. A final preview cannot stand in for an input or editing operation.

## Select the relevant guidance

Inspect the target `gbdraw/web/gallery/tutorials/<example-id>.json` and its
referenced media before editing. Read only the references needed for the request:

- Screenshots or capture metadata: [capture.md](references/capture.md).
- Text, captions, or structured tables: [content.md](references/content.md).
- Runtime, restored state, or Gallery generation ownership: the relevant sections
  of `gbdraw/web/CLAUDE.md`.
- Broad screenshot audits: the active operation register at
  `docs/internal/WEB_GALLERY_OPERATION_SCREENSHOT_REGISTER.md`, if present.
  Consult older plans only to resolve an otherwise unclear requirement.

A caption correction does not require a new capture if the existing image and
capture contract remain accurate. A resolution-only refresh preserves the same
semantic crop. A full tutorial audit covers images, prose, and structured content
in the requested inventory; batching does not reduce that coverage.

## Capture and audit decisions

Start with media referenced by tutorial JSON. For an audit, classify each as
`keep`, `recrop`, `replace`, or `add`, comparing the image with its instruction
and caption. Record broad-audit decisions in the active register (create it at
the path above if absent); keep a local edit's rationale with the change.

Data-dependent operations require `dataDependent: true`, example-local media,
an explicit `capture.session`, exact `capture.assertAppState` paths, and expected
`capture.visibleControls` and/or `capture.visibleText`. `genericMedia: true` is
only for crops with no example-specific file, metadata, setting, or result.
Prove these assertions inside the final crop before writing it.

Use the existing `tools/capture_gallery_tutorial_screenshots.py` capture owner.
For ordered controls, also replay the documented moves from reset/default state;
a correctly arranged restored session does not prove the instructions work.
A Gallery inspection crop may restore its example session. That does not replace
fresh-input GUI evidence when a tutorial teaches the reader to build the figure.

Compare every replacement with the old image at the same displayed size.
For broad audits, inspect contact sheets and check adjacent/identical images for
redundancy without sacrificing compactness or readability. After changing
references, remove stale media only after checking other callers. Report manual-only
media separately when an "all screenshots" refresh cannot reproduce them.

Generator-owned Gallery outputs follow `AGENTS.md`'s overwrite rule. Generate
from declared inputs; keep unrelated files and the owner-maintained social
preview out of scope.

## Verification and completion

For changed tutorial JSON and capture contracts, run:

```bash
python -m json.tool gbdraw/web/gallery/tutorials/<example-id>.json
python tools/capture_gallery_tutorial_screenshots.py --example <example-id> --check
```

For Gallery renderer changes, check `gallery.js` and the relevant cases in
`tests/web/gallery-tutorial.playwright.spec.js`; use the Web guide's Python
Playwright fallback if Node Playwright is unavailable. Update existing assertions
for changed media references or table behavior.

Inspect changed images and layout at their published size on relevant desktop
and mobile views. Check readable controls, caption accuracy, intact short table
tokens, and unobscured targets. Reuse unchanged evidence. Finish with the affected
examples, captures/checks performed, and any specific unresolved capture path.
