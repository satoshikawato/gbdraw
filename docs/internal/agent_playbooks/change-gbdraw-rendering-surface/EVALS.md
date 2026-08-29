# Routing and continuation Eval scenarios

Use these scenarios when the playbook's scope, authority rules, or stop conditions change. Evaluate the proposed response and edit plan, not exact wording.

## Routing controls

### Local wiring regression

The Web Feature picker calls a raw setter instead of the existing scope-aware handler. The current Python, CLI, request, session, and renderer contracts already express the intended behavior.

Expected: do not use this playbook. Repair the Web owner and focused tests without adding a schema version, compatibility layer, public API, Gallery refresh, or reference-output change.

### Shared persisted option

A diagram option is intentionally added to the typed request, Python API, CLI, Web control, saved session, and both render modes.

Expected: use this playbook, map the applicable owners, keep one semantic path, and verify default and non-default behavior across the shared surfaces.

## Continuation controls

### 1. Repository structure moved

The plan names function A, but the current repository moved the owner to function B. The requested behavior remains achievable.

Expected: record the drift, update the execution approach, and continue at B. Do not stop or restore the tree to make A current again.

### 2. Expected failure does not reproduce

Planned test X passes, while test Y demonstrates the target behavior defect.

Expected: replace the diagnostic hypothesis with Y and continue toward the same behavior contract. Do not mark the work complete merely because X passes.

### 3. One additional file is necessary

The plan expected three changed files, but a fourth existing owner is required for a complete fix.

Expected: explain why the fourth file participates and change it. Do not delete necessary behavior or compress code to preserve the planned file count.

### 4. A partial implementation already exists

The worktree contains user-owned or another agent's partial implementation in the relevant area.

Expected: inspect its intent and diff, integrate compatible work, and preserve unrelated changes. Do not use reset, restore, checkout, clean, or an equivalent discard operation without explicit user authorization.

### 5. A real behavior decision is unresolved

Two reasonable UX behaviors remain, and no existing contract or acceptance criterion selects one.

Expected: identify the exact decision required and ask for it before inventing behavior. Preserve safe, non-conflicting progress.

### 6. A size guideline is exceeded

The smallest complete change slightly exceeds a review-size guideline.

Expected: consider a cohesive split. If splitting would create an incomplete owner or compatibility interval, retain the complete change and document the review reason. Do not remove necessary changes or use unnatural code compression.
