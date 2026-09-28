# S03 instruction prompt — numeric History guard and rule preparation

You have no prior chat context. Read [MASTER_PLAN.md](../MASTER_PLAN.md), S01 and S02 results at `../SESSION_RESULTS/S01.md` and `../SESSION_RESULTS/S02.md`, repository guidance, Web guidance, the accepted `PD-OI-048`–`050` scalar outcomes, and the Product/architecture ratchets. OS Japanese IME is excluded.

Use `<repo>/.worktrees/issue-619-boundary-followup-20260928` and **reuse S01's branch** `fix/issue-619-boundary-followup-20260928`. Inspect status/HEAD/upstream, fetch and fast-forward its same-named remote branch. Do not create another branch/worktree or alter another session's checkout. Preserve unrelated changes.

## Numeric draft

S01 must establish whether nonfinite JavaScript numbers can enter through `app/circular-track-slots.js::updateCircularTrackSlotMeasure()`. Existing `services/json-clone.js` uses JSON, which turns such numbers into `null` in a History snapshot. Keep invalid/incomplete **text** draft support, including `"Infinity"`, `"1e"`, and unsupported suffixes; these must remain visible and be rejected on Generate without changing the prior Result/request. Add the narrowest guard at the existing owning action before slot mutation/History capture, using scalar semantics from `app/track-slot-validation.js`. Reject nonfinite numeric leaves only; never coerce them to null, zero, or text. Do not replace the generic JSON clone with a broader clone or add a second scalar parser. If S01 proves this action is unreachable from supported user/API flows, document the supported boundary and add a focused regression at the actual reachable edge instead of a speculative global fix.

Test bare and typed nonfinite numbers, finite numbers, Auto, valid typed lexemes, invalid text, no-op edits, Undo/Redo, Save rejection, and unchanged request/Result/History after rejection. Verify the existing current writer and native typed request still reject invalid canonical scalars.

## History restore and rules

Use S01's call-count evidence for `app/app-setup.js::historySnapshots.setAfterApplyHistoryIntent`. A config-domain restore may need one `rulePreparation.prepare()` to reproduce effective rules. Keep it when required. If a particular call is redundant, narrow the existing callback's trigger or reuse the already prepared result; remove the superseded invocation in the same edit. Do not skip preparation for a domain whose rule-derived SVG changes. Test Undo/Redo of rule-changing and unrelated config edits, exact visible rule outcome, call count, no extra Worker construction, and a subsequent independent edit.

Review production/test/docs/generated diffs separately and run focused Node/browser and policy/architecture checks. Write `SESSION_RESULTS/S03.md` with the measured requirement for each call, precise owner/path effects, commands, Product classification, and rollback. Commit only justified changes and that result with an English title; verify branch/upstream and push only the same-named remote work branch without force. A no-runtime disposition still receives an evidence commit and push. Do not open a PR or merge.
