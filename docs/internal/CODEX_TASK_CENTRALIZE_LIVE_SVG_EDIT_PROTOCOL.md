# Codex task: centralize live SVG edit commit and replay

Repository: `satoshikawato/gbdraw`

## Prerequisite

Start from the latest `main` containing PR #336:

```text
Preserve loaded session intent on regeneration
```

Do not continue on the PR #336 branch and do not append this refactor to PR #336.

Before editing, record:

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git diff --stat
```

Preserve unrelated worktree changes. Do not reset, discard, broadly reformat, commit, push, or open a PR unless explicitly instructed later.

## Mission

Replace the distributed live-edit commit protocol with two explicit owners:

1. one owner for mutating the mounted SVG and synchronizing the selected Result; and
2. one owner for replaying canonical editor overrides onto a newly admitted SVG.

The user-visible behavior must remain unchanged:

```text
Feature fill, stroke, visibility, label text, label visibility,
legend fill, and legend stroke update immediately without Generate.

The selected Result is synchronized before the action returns.

Save without Generate preserves the edit.

Load restores the edit without constructing the diagram Worker.

Generate or label reflow reapplies the canonical edit state.

Undo and Redo restore both the canonical overrides and the visible Result.
```

Do not begin with another broad audit or performance investigation.

Do not stop after proposing an architecture. Implement the focused refactor, remove the superseded paths, and prove the ownership boundary with executable contracts.

## Why this is the next maturity task

PR #336 exposed four failures in one conceptual operation:

```text
Feature visibility changed the DOM but did not flush the selected Result.
Direct label visibility unnecessarily initialized the diagram Worker.
Label overrides could be cleared during generated-SVG admission.
Legend color overrides were omitted from generated-SVG admission.
```

These were not four unrelated defects. They resulted from one live-edit operation being implemented across multiple modules with separate decisions about:

```text
canonical override mutation
mounted DOM mutation
preview dirty state
selected-Result serialization
History intent
optional label reflow
new-SVG override replay
failure behavior
```

The repository guidance already requires a live edit to commit canonical state and the mounted Result before returning. Convert that prose invariant into a single executable ownership path.

## Read first

Read only the nearest relevant material:

```text
AGENTS.md
CLAUDE.md
gbdraw/web/CLAUDE.md

gbdraw/web/js/app/preview-runtime.js
gbdraw/web/js/app/candidate-render.js
gbdraw/web/js/app/run-analysis.js

gbdraw/web/js/app/feature-editor.js
gbdraw/web/js/app/feature-editor/color-actions.js
gbdraw/web/js/app/feature-editor/visibility-actions.js
gbdraw/web/js/app/feature-editor/label-actions.js
gbdraw/web/js/app/feature-editor/svg-actions.js

gbdraw/web/js/app/legend.js
gbdraw/web/js/app/legend/entry-actions.js
gbdraw/web/js/app/legend/stroke-actions.js

gbdraw/web/js/services/svg-serialization.js
gbdraw/web/js/services/svg-result-ingestion.js
gbdraw/web/js/services/history.js
gbdraw/web/js/services/history-snapshot.js

tests/web/preview-runtime.test.mjs
tests/web/feature-color-actions.test.mjs
tests/web/feature-visibility-actions.test.mjs
tests/web/history.test.mjs
tests/web/architecture-contracts.test.mjs
tests/web/contracts/session-regenerate-intent.playwright.spec.js
tests/web/contracts/current-session-lazy-materialization.playwright.spec.js
```

Inspect additional direct dependencies only when required.

Do not create another dated audit, roadmap, or incident document.

## Current ownership to verify

The current implementation has a useful low-level runtime in:

```text
gbdraw/web/js/app/preview-runtime.js
```

but direct-edit callers still independently decide when to:

```text
mark the mounted SVG dirty
serialize the SVG
replace the selected Result
set skipCaptureBaseConfig
invalidate indexes
queue label reflow
```

Override replay for a fresh generated or reflowed SVG is also distributed across:

```text
candidate-render.js
label-actions.js
legend helpers
run-analysis.js
session restoration
```

Verify the actual current call graph before editing.

Do not retain multiple equivalent commit paths merely for compatibility with current internal call sites.

## Required architecture

### A. Mounted SVG commit owner

Make `preview-runtime.js`, or one narrowly named adjacent module, the only production owner for committing a mounted-DOM edit into the selected Result.

Introduce an explicit operation such as:

```js
previewRuntime.commitDomEdit({
  reason,
  mutate,
  invalidateIndexes
})
```

The exact name and shape are your decision.

Required semantics:

1. Resolve the active mounted SVG and selected Result exactly once.
2. Invoke the supplied synchronous DOM mutation.
3. Treat a returned `false` or zero change count as a no-op.
4. Invalidate only the indexes named by the operation.
5. Mark the active Result dirty.
6. Serialize the SVG exactly once.
7. Replace the selected Result exactly once.
8. Clear the dirty state after a successful flush.
9. Return a structured result such as:

```js
{
  changed,
  flushed,
  resultIndex,
  reason
}
```

10. Never construct or invoke the diagram Worker.
11. Never call Python.
12. Never create a Generate artifact-replacement History entry.

A compound action affecting several DOM elements must still perform one Result serialization and replacement.

A no-op action must perform zero Result serialization and replacement.

### B. Canonical editor-state projection owner

Introduce one pure or narrowly stateful owner for applying canonical editor state to a newly admitted SVG.

A suitable module name is:

```text
gbdraw/web/js/app/editor-svg-projection.js
```

The exact name is your decision.

It must own the deterministic replay of:

```text
Feature fill overrides
Feature stroke overrides
Feature visibility overrides
Label text overrides
Label visibility overrides
Legend color overrides
Legend stroke overrides
Specific-rule-derived Feature fills
Any current suppression flag that changes admitted SVG content
```

The projection owner must operate on an SVG plus explicit input data. It must not read arbitrary global state internally.

Use shared low-level DOM mutators where practical so direct edits and replay agree on:

```text
Feature identity
fill targets
stroke targets
visibility semantics
label identity
legend identity
paint normalization
```

Do not make a fresh SVG pass through UI action functions.

### C. Action-module responsibility

After the refactor, Feature, Label, and Legend action modules own only:

```text
scope resolution
validation of user input
canonical override-state mutation
History intent boundaries already owned by the application
optional geometry-reflow policy
calling the mounted-SVG commit owner
```

They must not independently:

```text
serialize the full mounted SVG
replace state.results[selectedResultIndex]
implement another dirty/flush protocol
invoke Generate to apply an ordinary direct edit
apply persisted editor state to a newly generated SVG
```

### D. Optional reflow semantics

Preserve the established live-edit rule:

```text
directly applicable edit
→ canonical override commits
→ mounted DOM and selected Result commit immediately
→ optional reflow may follow
```

If optional reflow fails:

```text
the already committed direct edit remains visible and saved
the error is reported
no rollback to the pre-edit SVG occurs
```

When automatic reflow is disabled and the mounted target was updated successfully:

```text
no Worker construction
no Worker initialization
no Worker run
```

Only edits that genuinely require creating, removing, or repositioning geometry may queue the owning rerender.

Do not use a hidden Generate as the common direct-edit mechanism.

## Required migrations

Migrate at least these production operations to the single mounted-SVG commit owner:

```text
single Feature fill
bulk Feature fill
single Feature stroke
bulk Feature stroke
single Feature visibility
bulk Feature visibility
label text
label visibility
legend fill
legend stroke
legend rename when it mutates the mounted SVG
```

Where an operation currently has a fallback path without `previewRuntime`, either:

1. route the fallback through the same owner; or
2. prove it is test-only or unreachable in production and remove it.

Do not leave two production serialization paths.

Migrate these fresh-SVG paths to the canonical projection owner:

```text
normal Generate candidate admission
label-reflow admission
current-session preview restoration when editor overrides must be replayed
trusted History artifact restoration only if it currently replays rather than restores an owned artifact
```

Do not route trusted immutable artifact restoration through untrusted SVG admission.

## Structural architecture contracts

Extend `tests/web/architecture-contracts.test.mjs` with durable assertions.

Required properties:

```text
Direct editor modules do not assign Result SVG content directly.
Direct editor modules do not call serializeCleanSvg directly.
Only the mounted-SVG commit owner calls flushActiveResult for live edits.
Only the canonical projection owner applies the complete persisted editor override set to a fresh SVG.
Feature/Label/Legend action modules do not import diagram-generation or run-analysis owners.
A direct edit has no Worker or Python call path.
```

Use import-graph and source-owner assertions consistent with the existing architecture test style.

Do not assert fragile exact occurrence counts for incidental helper names when an ownership/import assertion is sufficient.

## Executable direct-edit matrix

Create a reusable test matrix rather than adding another long one-off sequence.

The matrix must cover:

```text
Feature fill
Feature stroke
Feature visibility
Label text
Label visibility
Legend fill
Legend stroke
```

For each case prove:

```text
canonical override changed
mounted DOM changed immediately
selected Result SVG changed before the operation returned
Result replacement count = 1
SVG serialization count = 1
Worker construction delta = 0
Worker initialization delta = 0
Worker run delta = 0
Python helper/render delta = 0
Feature catalog owner unchanged
extracted Feature owner unchanged
biological Feature owner unchanged
orthogroup owner unchanged
```

For a no-op repeat of the same value prove:

```text
Result replacement count = 0
SVG serialization count = 0
History count unchanged
```

For one bulk edit affecting multiple Features prove:

```text
multiple DOM targets changed
one Result serialization
one Result replacement
```

## Save, load, Generate, and History contracts

Retain the PR #336 browser contract, but refactor reusable setup and evidence helpers out of the monolithic spec where doing so reduces duplication.

Do not weaken the existing semantic or visual assertions.

Required workflow:

```text
Load current session
→ perform every representative direct edit without Generate
→ Save without Generate
→ fresh Load
→ confirm every edit immediately
→ first Generate
→ confirm every edit is replayed
→ Undo Generate
→ confirm loaded edited artifact restored
→ Redo Generate
→ confirm regenerated edited artifact restored
```

Before the first Generate, prove:

```text
Worker construction = 0
Worker initialization = 0
Worker runs = 0
resource transfer = 0
```

During direct-edit Undo/Redo, prove:

```text
canonical overrides restore
mounted DOM restores
selected Result restores
Worker/Python remain idle
no artifact-replacement checkpoint is created
```

## Label-context preservation

Preserve and explicitly test the PR #336 rule:

```text
A changed SVG context key does not clear Feature-scoped label overrides
when their stable target Features still exist in the admitted SVG.
```

Required cases:

```text
same targets, different generated geometry → retain overrides
target removed by active biological/render filtering → remove or mark unresolved according to existing semantics
ambiguous target identity → fail closed
reflow failure → retain the already committed direct edit and canonical override
```

Do not preserve stale label overrides by rendered positional ID alone.

## New-SVG replay equivalence

For each direct-edit domain, prove that these two paths produce equivalent committed SVG state:

```text
Path A:
Generate
→ direct edit mounted SVG

Path B:
same active canonical override state
→ Generate or reflow
→ replay override onto fresh SVG
```

Compare the relevant DOM attributes and canonical override state. Use full semantic/raster comparison only for the focused combined fixture, not for every unit case.

## Test maintainability

The current session-regeneration Playwright contract is over 2,000 lines.

While adding the new matrix:

- extract stable, reusable helpers into a clearly named test helper module;
- keep scenario-specific assertions in the spec;
- do not create a generic test framework;
- do not duplicate the large evidence-capture function in another file;
- reduce or hold the total line count of the main spec where practical.

Do not mix unrelated test cleanup into this PR.

## Static API documentation

Add concise JSDoc typedefs for the new commit and projection boundaries.

Document:

```text
owner
inputs
mutation behavior
return value
no-op behavior
failure behavior
Worker/Python prohibition
```

Do not convert the application to TypeScript in this PR.

## Verification order

Run focused Node tests first:

```bash
node --test \
  tests/web/preview-runtime.test.mjs \
  tests/web/feature-color-actions.test.mjs \
  tests/web/feature-visibility-actions.test.mjs \
  tests/web/history.test.mjs \
  tests/web/architecture-contracts.test.mjs
```

Add focused tests for the projection owner and direct-edit matrix.

Then run:

```bash
npx playwright test \
  tests/web/contracts/session-regenerate-intent.playwright.spec.js \
  --workers=1 \
  --retries=0

npx playwright test \
  tests/web/contracts/current-session-lazy-materialization.playwright.spec.js \
  --workers=1 \
  --retries=0

node --test tests/web/*.test.mjs

npm run test:web:functional-smoke

ruff check gbdraw/

git diff --check
```

Run the real Vibrio generation contract once at the end:

```bash
npm run test:web:vibrio-generate
```

Do not repeatedly rerun it to select favorable timing.

## Correctness requirements

Preserve all established behavior:

```text
Current session loads without resource materialization.
Saved preview loads without constructing the diagram Worker.
Direct SVG edits work before first Generate.
Direct edits survive Save and fresh Load.
Active session intent remains authoritative for the next Generate.
Generate uses immutable artifact-replacement History.
Failed and canceled Generate restore exact pre-Generate owners.
Prepared biological input caches remain valid and bounded.
Unchanged warm Generate creates no History entry.
Current-session and frozen-v39 compatibility remain intact.
Generated SVG remains sanitized.
```

## Explicit non-goals

Do not implement:

```text
pre-Worker unchanged-Generate short-circuiting
generation-key design
final SVG caching
prepared biological input cache changes
session active-config schema extraction
session schema changes
canonical request schema changes
feature-catalog schema changes
History redesign
render queue redesign
Canvas or WebGL
main-thread Pyodide
new UI features
scientific behavior changes
```

Do not change a direct edit into a Generate operation.

## Scope and diff discipline

The refactor must remove superseded paths in the same change.

Do not merely wrap the existing duplicated logic while leaving all callers intact.

Prefer:

```text
one mounted-SVG commit owner
one fresh-SVG editor projection owner
thin action modules
structural ownership tests
```

Avoid adding another general-purpose utility module.

Production line growth should be modest. A net reduction is preferred, but correctness and clear ownership are mandatory.

## Durable guidance

Update `gbdraw/web/CLAUDE.md` only with a concise invariant:

```text
Live editor actions own semantic scope and canonical overrides.
The preview runtime exclusively commits mounted DOM edits into Results.
The editor SVG projection exclusively reapplies canonical editor state to a fresh admitted SVG.
```

Do not add a long implementation narrative.

## Compatibility and artifact checks

Verify explicitly:

```text
SESSION_VERSION unchanged
canonical request schema unchanged
feature-catalog schema unchanged
no Gallery session changed
no reference SVG changed
no social-preview image changed
no dependency-file change
no generated browser wheel tracked
no runtime network dependency added
no public Python or CLI signature changed
```

## Completion report

Do not return only an audit or proposal.

Implement the ownership refactor and report:

1. starting SHA and worktree state;
2. previous direct-edit commit paths removed;
3. mounted-SVG commit owner and API;
4. fresh-SVG projection owner and API;
5. action-module responsibilities after the change;
6. optional reflow and failure semantics;
7. files changed;
8. production lines added, removed, and net change;
9. architecture assertions added;
10. per-domain direct-edit matrix outcomes;
11. no-op and bulk-edit serialization counts;
12. Worker/Python deltas for direct edits;
13. Save/load/Generate replay result;
14. direct-edit and Generate Undo/Redo result;
15. label-context retention/removal results;
16. exact commands and test outcomes;
17. Vibrio compatibility result;
18. schema and generated-artifact checks;
19. remaining risks;
20. proposed focused PR title and concise PR body.

Suggested PR title:

```text
Centralize live SVG edit commit and replay
```
