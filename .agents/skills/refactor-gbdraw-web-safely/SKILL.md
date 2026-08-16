---
name: refactor-gbdraw-web-safely
description: >
  Use for behavior-preserving gbdraw Web refactors, ownership centralization,
  or moving existing logic across mounted SVG, Results, editor replay, sessions,
  History, Workers, caches, or asynchronous boundaries. Establishes independent
  base characterization before production edits, preserves state provenance,
  forbids accidental copying of large owners, requires differential evidence,
  and ends with a separate adversarial review.
---

# Refactor gbdraw Web safely

Read `AGENTS.md` and the repository-level `CLAUDE.md` first. Read
`gbdraw/web/CLAUDE.md` for every task that touches the Web application.

Use this Skill when the task intends to preserve user-visible behavior while
moving, centralizing, splitting, replacing, or deleting existing Web code.
Examples include:

- moving live SVG edits behind a single owner;
- consolidating fresh-SVG replay;
- extracting a session or History boundary;
- changing Result, Worker, cache, or rollback ownership;
- replacing several call paths with one abstraction;
- moving an operation across a synchronous or asynchronous boundary.

Use `$change-gbdraw-rendering-surface` as well when the task intentionally
changes a setting, renderer behavior, topology, SVG semantics, Python/CLI
surface, or persisted format.

## Classify the task before editing

State explicitly which category applies:

```text
behavior-preserving refactor
intentional behavior change
bug fix
performance change
mixed task
```

For a mixed task, separate the behavior-preserving portion from the intentional
change. Do not describe an intentional semantic change as a refactor.

Record:

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git diff --stat
```

Identify the exact base commit that defines prior behavior.

## Establish the oracle before production changes

A behavior-preserving refactor begins with characterization, not implementation.

Before modifying production code:

1. Run the nearest existing tests on the base.
2. Add focused characterization tests for behavior not already fixed by an
   independent oracle.
3. Prefer a test-only commit or clearly separable test-only diff before the
   production refactor.
4. Do not rewrite characterization expectations merely because the new design
   produces different output.
5. If the old behavior is intentionally wrong, reclassify that part as a bug
   fix and state the new semantic contract explicitly.

A valid oracle is one or more of:

```text
base implementation output
an existing released fixture
a real renderer-produced SVG
a persisted session produced by the current writer
an independently specified state-transition table
a base-versus-head differential harness
```

The new implementation is not an independent oracle for itself.

Do not prove equivalence by comparing two paths that both call the same new
helper. Such a test proves shared implementation, not preserved behavior.

## Map state provenance

Before consolidating code, list every affected state domain and its owner.

For each value, record:

```text
producer
validator
mutation owner
persistence owner
rollback owner
fresh-artifact replay owner
default/reset semantics
```

Distinguish provenance even when two states use the same DOM attribute.

Examples:

```text
renderer-owned display="none"
editor-owned visibility override
History-owned previous artifact
session-owned active intent
mounted Result content
temporary candidate state
```

Clearing an editor override restores the renderer baseline. It must not infer or
overwrite that baseline from the current DOM after the editor has changed it.

For tri-state controls, write a transition table before implementation. At
minimum cover:

```text
base/default
explicit on
explicit off
reset/default
Undo
Redo
Save
Load
fresh Generate or reflow
```

## Audit data ownership and cost

Classify every input crossing the new boundary.

Use these categories:

```text
large immutable borrowed owner
small mutable copied value
ownership transfer
derived compact index
temporary candidate
```

Large values include, where applicable:

```text
Feature catalog
extracted Features
biological Features
orthogroups
sequence payloads
SVG strings or DOM trees
Results
LOSAT raw or derived payloads
session resource tables
```

Rules:

- Pure does not mean cloned. A pure function may borrow validated immutable data
  by reference.
- Do not use `JSON.parse(JSON.stringify(...))` as defensive copying for large
  owners.
- Do not serialize a large owner merely to compare, sign, estimate, freeze, or
  prove immutability.
- Copy only the small mutable domain that the boundary may modify.
- Prefer a compact index containing only fields consumed by the new boundary.
- Document who releases retained references and what invalidates them.

Add structural evidence for costly operations. Prefer:

```text
fullFeatureGraphCloneCount = 0
fullFeatureGraphSerializationCount = 0
unusedProjectionBuildCount = 0
fullSvgCloneCount = 0
```

over a tight wall-clock threshold.

## Keep mounted ownership synchronous by default

Mounted SVG and selected Result ownership is synchronous unless explicitly
proved otherwise.

Preferred lifecycle:

```text
resolve current target
validate every target and normalized value
mutate synchronously
serialize
commit to the same Result owner
release target
start optional asynchronous work afterward
```

Do not retain a mounted SVG, Result object, or Result index across `await`
without an explicit lease.

If asynchronous ownership is unavoidable, the lease must include enough
identity to reject stale settlement, for example:

```text
Result index
Result object identity
mounted SVG identity
document or runtime revision
generation key
```

Revalidate the lease after every `await` and immediately before commit. Fail
closed if any owner changed.

Test deterministic races:

```text
Result selection changes while suspended
same-index Result owner is replaced
Generate commits while suspended
session Load replaces the document
History restore replaces the artifact
mounted SVG is remounted
```

Do not use sleeps as synchronization.

## Make mutation failure semantics real

Do not claim that an operation “commits nothing” unless both the DOM and Result
are transactional.

Before mutation:

- resolve and validate all targets;
- normalize and validate all values;
- detect missing or ambiguous bindings.

Where a later failure is still possible, use a compact undo log for only the
attributes, text nodes, or child positions touched by the operation.

Do not clone the entire SVG for every edit.

After any failure, exactly one of these states is allowed:

```text
DOM unchanged and Result unchanged
DOM explicitly dirty with a tested recovery/retry path
```

Never leave:

```text
DOM changed
Result stale
dirty false
```

## Preserve one real owner

A centralization refactor is incomplete if callers still perform the operation
outside the new owner and invoke the owner only afterward.

Remove superseded paths in the same change.

Reject patterns such as:

```js
mutateDomElsewhere();
commitDomEdit({ mutate: () => true });
```

The owner must receive the target and perform the mutation inside its transaction
boundary, or the API must require and validate an explicit pre-existing mutation
lease.

Architecture tests should enforce capabilities and dependency direction, not
only the presence of function names.

## Build an operation matrix

Do not treat a broad domain name such as “Legend” or “visibility” as one case.

Inventory user-visible operations separately. For live SVG work, typical rows
include:

```text
Feature fill: single and bulk
Feature stroke: single and bulk
Feature visibility: Feature, orthogroup, product, protein, rule edit, reset
Label text: direct and textPath
Label visibility: on, off, default, renderer-hidden baseline
Legend: add, remove, rename, fill, stroke, reorder, duplicate caption
layout drag and reposition
```

For each applicable row, cover:

```text
direct mounted edit
selected Result synchronization
no-op repeat
Undo
Redo
Save
fresh Load
fresh Generate
reflow
failure
Cancel or stale settlement where applicable
Circular
Linear
multi-Result isolation where applicable
```

## Use production-connected evidence

Counters and owner assertions must observe production code.

Invalid evidence includes:

```js
const execution = {
  workerRuns: 0,
  pythonCalls: 0
};
assert.equal(execution.workerRuns, 0);
```

Use actual runtime hooks, spies on production boundaries, browser lifecycle
events, Result owner identity, or independently derived output.

Test-only structural metrics must be:

- emitted by the production boundary;
- scoped to the operation under test;
- absent from persisted sessions, SVG, or public APIs.

## Run differential verification

For a behavior-preserving refactor, compare base and head on the same inputs.

Capture the smallest sufficient evidence:

```text
canonical state
normalized SVG semantics
relevant DOM attributes and order
Result content
History counts
Worker/Python calls
resource transfer
structural copy/traversal/serialization metrics
```

Byte equality is preferred where the old contract is byte-stable. Otherwise use
a documented canonical or semantic comparison.

Run real Circular and Linear fixtures for shared rendering behavior. Add a
multi-Result fixture when Result ownership is involved.

## Separate implementation from adversarial review

After implementation and focused tests pass, use a separate agent or session to
review the completed diff.

The reviewer is not asked to confirm compliance. It is asked to find a reason
not to merge.

The adversarial review must search for:

```text
new deep clones or full serialization
new full traversals
renderer/editor provenance loss
ownership crossing await
stale Result or session settlement
new parallel pipelines
placeholder commit callbacks
partial mutation on error
tests that use the implementation as their oracle
dummy counters or dummy owner objects
unclassified CI failures
```

Do not let the implementation session silently waive findings.

## Split large refactors

Default split threshold:

```text
more than 15 production files
or
more than 1,000 changed production lines
or
more than one independent state machine
```

Above the threshold, split unless the completion report explains why atomicity
requires one PR.

A safe sequence is often:

```text
characterization tests
mounted synchronous owner
fresh replay for one domain
fresh replay for the next domain
session/History/reflow migration
superseded path removal and architecture enforcement
```

## Required verification

Run focused tests first, then the affected broad gates.

For Web ownership refactors, always include:

```bash
npm run test:web:refactor-guards
node --test tests/web/*.test.mjs
npm run test:web:functional-smoke
ruff check gbdraw/
git diff --check
```

Run the applicable real-session, LOSAT, performance, Gallery, and reference
gates when their owners are touched.

Do not declare completion while any required CI job is red or unclassified.
Reproduce failures against the exact base before calling them pre-existing.

## Completion report

Return:

1. task classification;
2. base SHA and worktree state;
3. characterization tests established before production edits;
4. state provenance table;
5. large-data ownership and cost table;
6. synchronous or leased ownership model;
7. superseded paths removed;
8. base/head differential evidence;
9. production-connected structural metrics;
10. Circular, Linear, and multi-Result coverage;
11. failure and race coverage;
12. exact commands and outcomes;
13. all CI failures and their classification;
14. separate adversarial-review findings and dispositions;
15. production/test line counts;
16. remaining risks;
17. proposed commit title and summary.

A behavior-preserving refactor is incomplete if its independent oracle,
ownership audit, adversarial review, or required CI classification is missing.
