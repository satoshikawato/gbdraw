# Session 02 instruction prompt — automatic application of resolved alignment

Paste this entire file into a new Codex session. It is self-contained.

## Mission and gate

gbdraw Linear Similarity Group alignment has a Python helper response with status resolved or ambiguous. The current Web controller always opens Select alignment anchors. The selected Choice A outcome requires normal Align to apply a validated resolved plan directly, while ambiguity still opens the full review. A named explicit review action will be added in Session 03. Preserve the ability to request full review by exposing a controller entry mode now; do not add a second alignment controller.

Use exactly fix/similarity-alignment-auto-apply-definition-column-20260925. Before changing dependent runtime, verify that the complete Choice A authority for PD-OI-026/027/031/034 has merged into origin/dev. Fetch origin/dev, inspect the merged records and SHA, and merge that base into the fixed branch. If authority is not merged, stop this dependent runtime session and report the exact gate; the independent Definition correction from Session 01 may remain on this branch. Never treat the Decision Pack or unmerged authority candidate as runtime authorization.

## Read and map

Read AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md, 01_MASTER_PLAN.md, 00_DECISION_PACK.md, the merged Product Contract, docs/internal/PRODUCT_IMPACT_RATCHET.md, and docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md. Inspect gbdraw/web/js/app/similarity-alignment.js, its app-setup.js call sites, the Python response adapter, and the current tests/web/similarity-alignment-actions.test.mjs and tests/web/similarity-alignment-ui.playwright.spec.js. Confirm response.status is derived by Python; do not use unresolvedCount, selected candidate count, or visual badges as a substitute.

## Implementation

Change the existing controller entry to distinguish normal Align from an explicit review request. Both must resolve the same exact feature reference and use the same Python helper. After validating the helper response and current artifact:
- Normal Align plus resolved: use the response's validated plan, enter applying, and call the existing applyPlan / runAnalysis path once. Preserve direction by default. Publish Result, History, summary, and focus/status feedback only after successful generation. Do not mount the palette for a frame or make a redundant second resolver call.
- Normal Align plus ambiguous: build the existing complete draft and enter reviewing. Keep recommendation reasons, Select, Skip, canvas selection, orientation, and batch Apply.
- Explicit review plus either status: build the same full draft and enter reviewing without commit. Session 03 will supply the visible action.
- Automatic generation failure: retain a usable resolved draft and present the existing review with an actionable error and Apply retry. The current Result, active plan, and History remain intact. Resolver failure without a validated response reports an error and leaves state intact.
- Stale, canceled, or superseded work cannot install a plan, Result, or History entry. Preserve immediate busy feedback, no-duplicate-start behavior, and one action per successful alignment.

Do not change the Python recommendation algorithm, typed plan schema, Session writer, Worker lifecycle, or the CLI. Keep one state machine and remove obsolete unconditional-review assumptions.

## Verification and review

Add controller tests for resolved only-usable and unique-direct-RBH cases, genuinely ambiguous cases with preselected recommendation, missing/unusable targets, explicit review of a resolved response, error/retry, stale completion, and one Result/History change. Verify no second helper operation on the resolved default and no palette mount before success or failure. Extend a focused browser test with real Worker/rendering behavior; use Python Playwright if Node's test runner is unavailable. Run relevant Node/Python tests, ruff if Python changes, architecture contracts, and git diff --check. Review production and test diffs separately and record owner/path evidence. Commit verified implementation on the fixed branch; push only if authorized for this session.

## Handoff

Report merged authority SHA, branch/upstream, changed files, tests, observed resolved/ambiguous/retry behavior, and commit. At the end, print the entire contents of SESSION_03_REVIEW_ENTRY_AND_DOCS.md in a copyable code block as the exact prompt for the user's next Codex session.
