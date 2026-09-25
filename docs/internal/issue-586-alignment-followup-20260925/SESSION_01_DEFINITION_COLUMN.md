# Session 01 instruction prompt — repair the locked Definition column

Paste this entire file into a new Codex session. It is self-contained.

## Mission and branch

In gbdraw Linear diagrams, Lock Definition Column must align row/side Definition blocks on one common left edge. A Similarity Group Align can assign different horizontal translations to records; the current drawing and collision code add those translations even to the locked Definition origin, so the blocks drift. This is a regression against accepted PD-OI-024 and needs no new Product choice.

Use exactly fix/similarity-alignment-auto-apply-definition-column-20260925, created from origin/dev. Do not create a replacement implementation branch. Before editing, verify branch, upstream, git status, and the state of origin/dev. Preserve unrelated changes, especially the pre-existing untracked docs/internal/gbdraw_meet_gbdraw_implementation_v2/ directory. Do not edit Web auto-apply runtime in this session; that change waits for separate authority.

## Read and trace

Read AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md, docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md, PD-OI-024 in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md, docs/internal/issue-586-alignment-followup-20260925/01_MASTER_PLAN.md, and relevant code/tests.

Trace place_linear_definition in gbdraw/layout/linear.py, final_translations and _record_collision_bands in gbdraw/diagrams/linear/assemble.py, and add_record_definition_group in gbdraw/diagrams/linear/builders.py. The existing keep_left path computes a common x, then both paint and collision add record-specific translation_x. Inspect ordinary one-record rows and multi-record row headings separately. A record-local label above its sequence must continue to follow the record.

## Implementation and acceptance

Use one final shared column origin after final record translations. Ensure the locked Definition blocks have one left edge and enough configured gap from the leftmost displayed sequence, including negative translations. Feed the same placement into paint and collision planning; avoid copied formulas and new layout owners. Keep Lock-off definitions following their record rows, and keep local labels over the corresponding sequence. Preserve explicit text-anchor behavior within its existing scope.

Add focused deterministic Python coverage and an actual browser SVG geometry assertion with unequal nonzero positive and negative record translations. Measure transformed Definition left edges and the nearest sequence edge; also cover one-record rows and mixed multi-record rows, Lock off, and save/load plus regeneration where the existing browser harness makes that practical. A simple assertion of the text-anchor attribute is insufficient. Adapt the existing tests/test_linear_definition_alignment.py and tests/web/linear-multi-record.playwright.spec.js rather than duplicating fixtures without need. Treat tracked reference_outputs as read-only until an intentional output change is reviewed.

Run focused pytest and browser checks, ruff on changed Python, and the relevant output-comparison test. Check Node and Python Playwright before choosing the browser harness; if Chromium fails only because of the sandbox, rerun with the required escalation. Diagnose a failing geometry boundary rather than lowering tolerance to hide it. Review production and test diffs separately, record concise owner/path evidence under the architecture ratchet, and commit a verified correction only on the fixed implementation branch. Push only if authorized for that action in the session.

## Handoff

Report exact changed files, test commands/results, measured invariant, any effect on SVG references, commit, and remaining work. Session 02 cannot begin until the Choice A authority-only change is merged into origin/dev and merged into the fixed implementation branch. At the end, print the entire contents of SESSION_02_AUTO_APPLY_CONTROLLER.md in a copyable code block as the exact prompt for the user's next Codex session, and state the authority prerequisite immediately before that block.
