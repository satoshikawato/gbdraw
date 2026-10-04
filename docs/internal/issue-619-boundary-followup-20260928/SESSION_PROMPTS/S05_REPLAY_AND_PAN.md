# S05 instruction prompt — semantic replay and stable pan measurement

You have no prior chat context. Read [MASTER_PLAN.md](../MASTER_PLAN.md), S01–S04 results under `SESSION_RESULTS/`, repository/Web guidance, and the public [Session compatibility contract](../../../REFERENCE/session-and-request-compatibility.md). OS Japanese IME is excluded.

Use `<repo>/.worktrees/issue-619-boundary-followup-20260928` and **reuse S01's branch** `fix/issue-619-boundary-followup-20260928`. Inspect status/HEAD/upstream, fetch and fast-forward its same-named remote branch. Do not create another implementation branch/worktree or alter another session's checkout.

## Native replay

Use the existing Web Save and `python -m gbdraw.cli ... --session ... --session_output ...` path on a generated Session and a Session whose draft differs from its committed Result. Pin inputs and record hashes. Verify source bytes, biological record IDs, scalar values and units, request grouping, resolved geometry, and drawing element/text/transform structure. Use exact equality for scientific identities and values; use the already documented narrow absolute numeric tolerance for Web wasm/native floating-point SVG attributes where needed. Record the compared attribute set, maximum difference, tool/font versions, and any omitted/default native sidecar fields. Preserve the existing public statement that cross-version/font SVG byte identity and full sidecar JSON equality are not promised. Do not add a writer normalization or new Session schema merely to force byte equality.

If a semantic mismatch appears, diagnose the first authoritative boundary (typed Session decode, request projection, planner, or renderer) and fix it in that existing owner; do not hide it by widening the numeric tolerance. Add or refine a focused regression in the existing Session/replay test owner. If only serialization differs, document the proven semantic scope and leave runtime untouched.

## Pan retry

The hosted first attempt in `preview-navigation.playwright.spec.js` had pan state `(60,80)` and zero scroll delta, but x screen displacement `60.5802`; the before/after SVG 100-unit screen distance changed from about `30.0439` to `29.9996`. The existing `settle()` only matches CSS transform to pan/zoom state. Reproduce with retries disabled. If a pre-sample layout race is confirmed, wait for stable SVG CTM scale and wrapper geometry before capturing `before` in this **test's existing helper**; retain its `toBeCloseTo(..., 1)` oracle, timeout, and mouse gesture. Avoid a global arbitrary sleep or a new generic polling framework. If the geometry is stable and displacement remains wrong, fix the existing preview-navigation runtime owner instead and prove the user effect.

Run several focused no-retry repetitions after the final change, retaining every attempt in the result. Record exact commands, attempts, pan/scroll/CTM measurements, semantic replay evidence, owner/path effects, and rollback in `SESSION_RESULTS/S05.md`. Review production/test/docs/generated diffs separately; run applicable policy/architecture checks. Commit scoped changes and the result with an English title, verify branch/upstream, and push only the same-named remote work branch without force. A no-runtime parity result still receives a commit/push. Do not open a PR or merge.
