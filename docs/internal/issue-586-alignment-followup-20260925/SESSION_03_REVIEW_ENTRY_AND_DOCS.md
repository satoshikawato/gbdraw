# Session 03 instruction prompt — explicit review entry and public wording

Paste this entire file into a new Codex session. It is self-contained.

## Mission and branch

The selected Choice A for gbdraw Linear Similarity Group alignment gives the normal Align… action a one-action resolved path. Users who want to inspect or change an anchor, Skip, or Match reference direction before commitment must have a visible Review alignment options… action at both the exact-feature popup and Similarity Groups drawer. Ambiguous resolutions automatically open the existing full Select alignment anchors palette.

Use exactly fix/similarity-alignment-auto-apply-definition-column-20260925. Confirm the Choice A authority for PD-OI-026/027/031/034 is merged into origin/dev and the implementation branch, and that Session 02's controller branch/review entry mode is present. If either is missing, do not invent a parallel entry controller. Check branch/upstream and preserve unrelated changes.

## Read and implement

Read AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md, 01_MASTER_PLAN.md, the merged Product Contract, and the current Web reference. Inspect gbdraw/web/index.html and gbdraw/web/js/app/app-setup.js for popup/drawer buttons, busy state, focus restoration, Escape, palette visibility, canvas overlay, and narrow viewport behavior. Inspect tests/web/similarity-alignment-ui.playwright.spec.js and focused action tests.

Expose a named Review alignment options… action next to the existing Align… action in both entry points. Each passes the exact same selected reference and group into the one controller with explicit review intent. Give both actions clear enabled, busy, and accessible names/states. Avoid putting the only review route behind hover, modifier keys, or a pointer-only menu. Ensure the palette opens promptly for explicit review even when Python returns resolved; it retains every target row, candidate details, Select/Skip, orientation effects, canvas guides/badges, pan/zoom, keyboard, Apply/Cancel, and retry. Focus moves into the palette and returns to the initiating control on close. Normal resolved Align does not flash the palette. The Result summary and errors are understandable without seeing the palette.

Update docs/REFERENCE/web-app.md to distinguish normal resolved automatic application, ambiguous review, and explicit review. Audit existing help text and any Gallery tutorial that actually documents this control. If Gallery tutorials/screenshots need edits, read the web-gallery-screenshot-maintenance skill and regenerate through its supported capture process; do not hand-edit generated assets. Do not add new public pages when the existing Web reference answers the user question.

## Verification and review

Add focused desktop and 390-pixel browser assertions for both entry points, exact reference identity, immediate busy feedback, keyboard discovery and focus, orientation/Skip on explicit resolved review, ambiguous review, automatic success, and failure/retry. Confirm Session/result/download do not include transient guide or draft data. Run focused Node/Python tests and browser checks, architecture contracts, and git diff --check. Check both Node and Python Playwright installations; rerun a Chromium sandbox-only failure with escalation. Review production, tests, docs, and generated diffs separately. Record concise owner/path evidence; commit on the fixed branch only after verification, and push only if the session authorizes it.

## Handoff

Report changed files, browser and accessibility observations, test commands/results, documentation updates, and commit. At the end, print the entire contents of SESSION_04_ACCEPTANCE_AND_HANDOFF.md in a copyable code block as the exact prompt for the user's next Codex session.
