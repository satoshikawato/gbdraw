## Summary

<!-- What changes, and why? -->

## Verification

<!-- Exact commands and outcomes. -->

## Behavior-preserving Web refactor

Complete this section when existing Web behavior is moved, centralized,
split, replaced, or removed.

- [ ] I classified the task as behavior-preserving, intentional behavior change,
      bug fix, performance change, or mixed.
- [ ] Base characterization tests were established before production edits.
- [ ] The new implementation is not used as its own behavioral oracle.
- [ ] Renderer-owned and editor-owned state provenance was audited.
- [ ] Large borrowed owners were audited for clone, serialization, traversal,
      retention, and cleanup.
- [ ] Mounted SVG or selected Result ownership does not cross `await`, or an
      explicit lease is revalidated and race-tested.
- [ ] Superseded production paths were removed.
- [ ] Base/head differential evidence was recorded.
- [ ] Structural metrics observe production code rather than local dummy values.
- [ ] Circular, Linear, and multi-Result applicability was assessed.
- [ ] Failure, Cancel, stale settlement, Undo, Redo, Save, Load, and replay were
      covered where applicable.
- [ ] A separate adversarial review inspected the completed diff.
- [ ] `npm run test:web:refactor-guards` passes.
- [ ] Every required CI failure is green or reproduced and classified against
      the exact base.

## Artifacts and compatibility

- [ ] No unintended session, request, or feature-catalog schema change.
- [ ] No unintended Gallery session or reference SVG change.
- [ ] No generated browser wheel was tracked.
- [ ] No runtime network dependency was added.
