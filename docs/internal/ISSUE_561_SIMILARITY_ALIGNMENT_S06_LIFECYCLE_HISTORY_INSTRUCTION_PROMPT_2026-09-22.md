# Issue #561 — S06 lifecycle and history instruction prompt

```text
Implement Session S06 of gbdraw Issue #561: complete active-plan persistence,
materialization, invalidation, Reset, sequential Align, Undo/Redo, and the bridge
from existing composition deltas to stable record-keyed base translations.

Repository and branch:
- Use issue-561-similarity-alignment.
- Confirm merged authority and complete S01-S05 ledger entries.
- Preserve unrelated changes; do not create another runtime branch.

Read before editing:
1. AGENTS.md, CLAUDE.md, gbdraw/web/CLAUDE.md
2. master plan sections 4.3-4.5, 7.2, 7.5, 8.3-8.5, and 13
3. accepted lifecycle, reset/history, and compatibility decisions
4. S01-S05 code and tests
5. gbdraw/web/js/services/history.js and history-snapshot.js
6. gbdraw/web/js/app/legend-layout/composition-actions.js and diagram-drag.js
7. source/crop/selector/orientation/reorder mutation owners in linear-sources,
   watchers, app-setup, reset, and config
8. run-analysis Result admission and stale/cancel handling

Implement one coherent state machine:
1. The base record translation map and base RecordPresentation values are the
   immediate pre-align baseline. The active plan is an overlay.
2. Ordinary Generate after style, label, or canvas-size edits retains and
   revalidates the active plan.
3. Stable record reorder remaps by recordKey and retains the plan.
4. Starting a new Align materializes the old plan's effective translations and
   orientations as base values, removes the old plan, and installs the new plan
   in one successful artifact transaction.
5. Reset Align removes the active plan and renders the stored base values. It
   does not revive an older active plan. Normal Undo may restore that artifact.
6. Manual record drag or manual orientation change first materializes the active
   effective state and clears the plan with a visible reason, then applies the
   edit in the same user operation.
7. Source replacement, crop change, and selector change clear the plan with a
   visible reason through their existing semantic mutation owners. Avoid a
   broad deep watcher that guesses why state changed.
8. Convert current per-record composition deltas to/from stable keyed base
   translations at one adapter boundary. Prevent double application after
   Generate. Leave legend, title, length-bar, and whole-diagram deltas under
   their existing owners.
9. Validate before Generate. Stale reference blocks and offers Reselect/Clear.
   Stale target offers Select/Skip. Keep the last successful Result and saved
   preview during repair.
10. Apply, Reset, and plan-clearing manual edits each create one history entry.
    Failed, canceled, superseded, stale, or no-op operations create none.
11. Save/fresh Load preserves base translations, active plan, rationale, and
    Reset behavior. Load-only preview remains Worker-lazy.

Testing:
- style/label/canvas Generate retains alignment;
- stable reorder retains it by recordKey;
- source/crop/selector/orientation/drag clears it with the correct reason;
- Align A then Align B, Reset B, Undo/Redo sequence;
- Reset restores immediate geometry and clears the plan;
- record drag after alignment starts from effective rather than old base
  geometry and does not double-translate;
- same-row records retain separate keyed deltas;
- legacy composition array order is not used as persisted identity;
- stale reference/target recovery and current Result retention;
- failed/canceled/superseded/stale/no-op history counts;
- Session round trip and lazy preview;
- no LOSATP invocation during lifecycle-only operations.

Design constraints:
- Reuse existing artifact transactions; do not add an alignment history stack.
- Use explicit callbacks at semantic mutation owners; do not scatter resolver
  logic or duplicate the plan.
- Keep translations to X/Y and Linear mode only.
- If an operation cannot preserve a Product effect without a new behavior
  choice, stop and prepare a Decision Pack rather than inventing a fallback.

Finish:
- Run focused history, composition, session, source-mutation, Generate, and
  browser lifecycle tests plus architecture/change-budget checks.
- Review production/tests/session fixtures/generated output separately.
- Update master-plan section 16 and state whether S07 may start.
- Provide an English proposed commit title and summary; do not push or create a
  PR unless separately authorized.
```
