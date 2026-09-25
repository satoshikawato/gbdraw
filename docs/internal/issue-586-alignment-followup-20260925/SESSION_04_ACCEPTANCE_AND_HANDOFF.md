# Session 04 instruction prompt — integrated acceptance and handoff

Paste this entire file into a new Codex session. It is self-contained.

## Mission and prerequisites

Finish and verify the selected Choice A follow-up on fix/similarity-alignment-auto-apply-definition-column-20260925. The work includes: automatic application of Python-resolved Linear Similarity Group alignment; a full explicit Review alignment options… action at popup and drawer; full automatic review for ambiguity; and one locked Definition column after unequal record translations. Confirm Sessions 01–03 are present and the complete Choice A authority for PD-OI-026/027/031/034 is merged into origin/dev and the fixed branch. Do not treat a planning file or unmerged authority candidate as sufficient.

Read AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md, docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md, docs/internal/PRODUCT_IMPACT_RATCHET.md, and 01_MASTER_PLAN.md. Verify branch and upstream, preserve unrelated files, and inspect all in-scope production, test, documentation, and generated diffs separately once before repeating changed checks.

## End-to-end acceptance

Use deterministic fixtures to verify:
1. Normal Align with only one usable candidate, with a unique direct RBH among multiple candidates, and with missing/unusable targets returns a validated resolved plan, commits once without palette flash, and preserves default direction.
2. Genuine ambiguity opens the full preselected review; user can change a candidate, Skip, and set Match reference direction. Explicit Review alignment options… opens the same full review for a resolved plan at both popup and drawer.
3. An automatic render failure opens a retained draft, reports an actionable error, and supports corrected retry. A failed Apply retains the draft. Cancel, stale completion, repeated starts, and supersession preserve the last Result/History; success creates one action. Reset and Undo/Redo work from the committed result.
4. Save/load and regeneration preserve the typed active plan, source-relative orientation, and rendered result semantics; exports contain no transient guide, badges, or draft.
5. With Lock Definition Column on, use Chromium to measure the transformed SVG left edge of every row/side Definition after unequal positive and negative record translations. All share a left edge and retain the configured gap from the leftmost sequence. With Lock off, definitions follow rows; record-local labels follow their sequences. Cover ordinary and mixed multi-record rows, alignment, manual translation, save/load, and regeneration.

Render and visually inspect the user's five-record Streptomyces-style result or an equivalent realistic Gallery-quality recipe at readable scale, plus a 390-pixel viewport. Confirm text remains legible, the palette/actions and candidate list remain reachable, and closing the palette restores the preview. The approved PD-OI-035 revision-2 outcome permits the measured 390 px palette coverage and retires simultaneous canvas selection and manual pan/zoom while it is open at that width; desktop canvas interaction remains required. This revised criterion is authorized by PD-OI-035 revision 2, merged into `origin/dev@cfe28a86` through PR #593 and included in the fixed branch at `f26eb1f2`. Keep smoke fixtures in tests rather than public figures. Do not modify examples/gbdraw_social_preview.png.

## Required gates

Run focused tests first and fix concrete failures. Then run the applicable repository gates:
- python -m pytest tests/test_linear_definition_alignment.py tests/test_linear_multi_record_layout.py tests/test_similarity_alignment.py tests/test_similarity_alignment_rendering.py -v
- node --test tests/web/similarity-alignment-actions.test.mjs and any changed Web unit specs
- npx playwright test tests/web/similarity-alignment-ui.playwright.spec.js tests/web/linear-multi-record.playwright.spec.js --project=chromium --workers=1, when Node Playwright is installed; otherwise run equivalent targeted checks with Python Playwright
- ruff check gbdraw/
- node tests/web/architecture-contracts.test.mjs
- node tools/check-web-change-budget.mjs
- python tools/update_cli_reference_help.py --check
- python -m pytest tests/ -v -m "not slow"
- python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
- python -m build
- git diff --check

Prepare the gitignored browser wheel if packaging or browser tests need it. Allow at least 30 minutes for test commands and monitor incrementally. Do not update tracked reference outputs unless geometry intentionally changed, then review SVG diffs and rerun comparison. Run an offline browser audit only if runtime dependencies, privacy, bundle composition, or Worker lifecycle changed. Check both Playwright installation paths and rerun a Chromium sandbox-only failure with escalation.

## Architecture and final record

Confirm one Python resolver/validator, one Web controller, one shared Definition placement, one current request/Session writer, and one Result/History admission route. Remove superseded unconditional-review behavior and duplicate definition offsets. Check every jointly required Product effect rather than matching decision IDs alone. Provide concise before/after owner/path evidence, or full OE/PE/CB sets only if the architecture ratchet's exception criteria apply. Do not weaken a failing gate.

Record exact commits, merged authority SHA, tests and results, measured browser geometry, visual review, documentation and reference-output disposition, remaining limitations, and an English proposed commit title and short summary. Commit fixes on the fixed implementation branch after verification. Push only if authorized for this session, and do not create a PR or merge without the relevant authorization.

This is the final implementation session. In the final response, state that no next-session prompt is needed only if every acceptance criterion and required gate is complete. If work remains, create a self-contained continuation prompt naming the exact remaining task, branch, evidence, and next gate, then display its full text in a copyable code block for the user to paste into a new Codex session.
