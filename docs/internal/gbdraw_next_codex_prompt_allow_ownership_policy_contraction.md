# Codex task: make the Web ownership gate ratchetable downward

Repository: `satoshikawato/gbdraw`

## Base revision

Start from the fetched `origin/main` whose reference revision for this task is:

```text
43c14ddec49835458f6bac88aa134a7d755c89b1
Merge pull request #338: Add Web change-budget and ownership gate
```

Before editing, run and record:

```bash
git fetch origin main
git status --short --untracked-files=all
git rev-parse HEAD
git rev-parse origin/main
git diff --stat
```

Work from a clean branch based on that `origin/main`. Preserve unrelated worktree changes. Do not reset, discard, broadly reformat, commit, push, or open a pull request unless explicitly instructed later.

If `origin/main` has moved beyond the reference revision, inspect only the three in-scope files listed below. Continue only if the current checker still rejects every `tools/web-change-policy.json` contraction with `privileged capability allowlists may only expand`. Otherwise stop and report that this task has already been superseded.

## Mission

Make the Web ownership policy capable of **ratcheting downward** after runtime ownership has been reduced.

The current gate correctly prevents a pull request from introducing a new privileged owner or importer by editing runtime code and its allowlist together. Preserve that protection.

The current gate also rejects removal of any allowlist entry, including an unused stale entry. Change only contraction handling and the policy-key validation needed to distinguish a path-entry contraction from policy shape damage:

```text
privileged expansion
→ still requires guard-only preauthorization on the base revision

privileged contraction
→ allowed in a later guard-only pull request
→ only when the proposed policy still covers every privileged owner and importer
   present in the pull-request head runtime
→ removes path entries only; every capability and import-target key from the
   base policy remains present
```

This is a bootstrap change to the guard itself. It must not contract the real repository policy yet.

## Why this is the next task

PR #338 established a trusted-base checker, runtime/guard separation, change budgets, and privileged ownership allowlists. Those controls currently prevent ownership expansion but cannot record ownership reduction because `policyContractions` are always added to `integrityViolations`.

Without safe contractions, `tools/web-change-policy.json` remains a historical inventory rather than a ratcheting architectural boundary. Future narrow runtime migrations could remove a distributed owner, but the gate could not subsequently remove that path from the allowed set.

Do not begin another broad architecture audit. Do not modify Web runtime behavior. Implement this focused guard correction and prove it using the existing temporary-repository tests.

## Read first

Read only the nearest relevant files:

```text
AGENTS.md
CLAUDE.md
tools/check-web-change-budget.mjs
tools/web-change-policy.json
tools/web-change-source.mjs
docs/internal/WEB_CHANGE_POLICY.md
tests/web/architecture-contracts.test.mjs
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
```

The workflows are read-only context for this task.

## Strict scope

Production/runtime files under these paths must not change:

```text
gbdraw/web/index.html
gbdraw/web/js/
gbdraw/web/vendor/
```

Expected changed files are limited to:

```text
tools/check-web-change-budget.mjs
tests/web/architecture-contracts.test.mjs
docs/internal/WEB_CHANGE_POLICY.md
```

Do not change:

```text
tools/web-change-policy.json
tools/web-change-source.mjs
.github/workflows/test.yml
.github/workflows/web-base-policy.yml
package.json
package-lock.json
pyproject.toml
```

Do not add a new module, dependency, workflow, policy schema, command-line option, label, or general policy framework.

## Required semantics

### 1. Preserve base-policy authority for expansion checks

The base revision's allowlist must remain authoritative when deciding whether the pull-request head introduces a privileged owner or importer.

A pull request must still fail when it does either of the following in one change:

```text
add a new privileged runtime owner and add it to the proposed policy
add a new privileged importer and add it to the proposed policy
```

The existing two-PR preauthorization flow must remain intact:

```text
guard-only policy expansion merged first
→ later runtime implementation against that updated base
```

An `architecture-change` label must not waive this integrity rule.

### 2. Permit safe policy contractions

When `tools/web-change-policy.json` changes, parse the proposed head policy as data.

Removing an entry from either section must no longer be an integrity violation by itself:

```text
allowedPrivilegedOwners
allowedPrivilegedImporters
```

Keep reporting removed entries in the step summary. A successful contraction should therefore be visible but should not fail merely because it is a contraction.

Contraction in this task means removing path entries from an existing allowlist
array. The proposed policy must retain every capability and import-target key
present in the base policy. Deleting or renaming one of those keys is a policy
shape change, not a safe contraction, even when its current array is empty or
the head runtime has no matching use. Treat a missing base-policy key as a
distinct, non-waivable integrity violation. Do not rely on the pull-request
checkout's architecture test to preserve these keys.

### 3. Validate the proposed contracted policy against head runtime

Do not rely only on the pull-request checkout's architecture test to establish contraction safety. The trusted-base checker must be able to validate the proposed policy while treating pull-request files strictly as data.

Reuse the current source inventory, import parsing, capability specifications, and masked-source matching. Add a second coverage check when a proposed policy exists:

```text
base policy versus head runtime
→ detects unpreauthorized expansion

proposed policy versus head runtime
→ detects an over-contraction that excludes an active owner or importer
```

A proposed contraction must fail if the head runtime still contains the removed privileged owner or importer.

Use a distinct diagnostic, for example:

```text
proposed privileged capability allowlist excludes active owners or importers
```

The exact wording may differ, but it must clearly distinguish an unsafe contraction from an unpreauthorized expansion.

This failure is an integrity violation and must not be waivable with `architecture-change`.

Policy-key validation and proposed-policy runtime coverage are separate checks.
A proposed policy that retains every active owner and importer must still fail
if it deletes or renames a capability or import-target key from the base policy.

### 4. Keep trusted-base execution properties unchanged

The checker may read the pull-request head's policy and runtime files through Git, but it must not check out or execute pull-request code.

Do not alter the trusted-base workflow. Do not add package installation or third-party parsing dependencies.

### 5. Preserve runtime/guard separation

The existing rule must remain unchanged:

```text
Web runtime files and Web guard/CI files cannot change in the same pull request
```

A future ownership reduction therefore remains a sequence of narrow pull requests:

```text
runtime implementation PR removes the distributed owner
→ guard-only contraction PR removes the now-unused allowlist entry
```

## Implementation guidance

Keep the implementation local to the existing policy comparison and capability coverage logic.

A suitable structure is:

1. Parse the base policy exactly as today.
2. Parse the proposed head policy only when the policy file changes.
3. Confirm that every section key in the base policy remains present in the
   proposed policy; report a missing key as an integrity violation.
4. Compute and report removed path entries.
5. Evaluate head runtime against the base policy to preserve expansion protection.
6. Evaluate head runtime against the proposed policy to reject unsafe contractions.
7. Do not add path-entry contraction itself to `integrityViolations`.

Avoid duplicating the complete owner/importer scan. A small helper that evaluates one policy against the already-built `allProductionSources` is appropriate. Do not introduce a class, registry, manager, policy AST, or generalized rule engine.

Keep the current JSON policy format unchanged.

## Required tests

Use the existing temporary Git repository harness in:

```text
tests/web/architecture-contracts.test.mjs
```

The new contraction and policy-shape cases must exercise the same revision-based
path as the trusted-base workflow. For each case, commit the proposed change,
retain the baseline SHA, and invoke the checker with explicit revisions through:

```js
execute({ base: baseSha, head: headSha })
```

Do not use only `runChangeBudgetCase()` or another working-tree-only invocation
as acceptance evidence for these cases. A working-tree smoke case may remain,
but it does not replace explicit `--base`/`--head` coverage. The assertions must
therefore exercise `git show` and `git ls-tree` reads of pull-request head data.

Revise the current test named approximately:

```text
privileged allowlists cannot contract
```

and add focused cases proving all of the following.

### A. Unused owner contraction passes

The fixture policy contains `app/future-owner.js` as an allowed `Diagram Worker`
owner. Add that file to the baseline fixture with harmless executable code that
does not import a privileged target and does not match the `Diagram Worker`
owner expression. The positive case must therefore prove contraction of an
existing but inactive source path, not only removal of a path whose file does
not exist.

Remove only that unused owner from the proposed policy.

Prove:

```text
checker exit status = 0
Result: PASS
removed entry is reported
no integrity violation is reported
```

### B. Active owner contraction fails

Remove `app/editor.js` from the proposed `Diagram Worker` owner allowlist while the head runtime still calls `runDiagramGeneration()` there.

Prove:

```text
checker exit status = 1
diagnostic identifies the active owner excluded by the proposed policy
failure is an integrity violation
```

Also prove that `WEB_ARCHITECTURE_CHANGE=true` does not waive it.

### C. Active importer contraction fails

Remove `app/editor.js` from the proposed importer allowlist for `services/diagram-generation.js` while the import remains in head runtime.

Prove:

```text
checker exit status = 1
diagnostic identifies the active importer excluded by the proposed policy
failure is an integrity violation
```

Also prove that `WEB_ARCHITECTURE_CHANGE=true` does not waive it.

### D. Same-PR expansion remains blocked

The current suite has a runtime-only unapproved-expansion test; it does not yet
prove that a proposed policy can never self-authorize the runtime change. Retain
that runtime-only test and add an explicit revision-based case that, in one head
commit, both:

```text
adds a privileged runtime owner and importer
adds the same path to the proposed policy
```

Prove that the case reports both the runtime/guard co-change violation and the
base-allowlist owner/importer diagnostic. An exit status of 1 caused only by
runtime/guard separation is insufficient evidence for base-policy authority.
Also prove that `WEB_ARCHITECTURE_CHANGE=true` does not waive either integrity
failure.

Do not weaken or delete the existing two-step preauthorization test.

### E. Runtime/guard co-change remains blocked

Retain the existing test proving that runtime and guard files cannot change together.

### F. Reporting remains useful

For a passing contraction, retain a report section for removed allowlist entries. Rename headings or internal variables only if that improves accuracy; do not remove the evidence.

### G. Policy keys cannot disappear

Give the fixture policy at least one import-target key whose allowlist is empty,
for example:

```json
"workers/diagram-generation-worker.js": []
```

Delete that key from the proposed policy while leaving head runtime unchanged.
Prove:

```text
checker exit status = 1
diagnostic identifies the missing base-policy import-target key
failure is an integrity violation
WEB_ARCHITECTURE_CHANGE=true does not waive it
```

This case must fail because the key disappeared, not because an active importer
was excluded. It closes the gap that runtime-coverage checks alone cannot see.

## Documentation change

Update `docs/internal/WEB_CHANGE_POLICY.md` so it distinguishes expansion from contraction.

Document this exact operating sequence:

```text
Expansion:
1. guard-only preauthorization PR adds an owner/importer
2. later runtime PR uses it

Contraction:
1. runtime PR removes the owner/importer use without editing the policy
2. later guard-only PR removes the stale allowlist entry
```

State that the proposed contracted policy is checked against head runtime and
that runtime plus guard changes remain prohibited. Also state that contraction
removes path entries only: capability and import-target keys present in the base
policy cannot disappear, even when their arrays are empty.

Do not add a dated incident report, execution plan, roadmap, or new policy document.

## Explicit non-goals

Do not:

- modify any Web application behavior;
- change the current real allowlist contents;
- address Feature visibility, SVG commits, History, session ownership, or Worker lifecycle in this PR;
- alter ordinary file/line/module/dependency budgets;
- alter the `architecture-change` label behavior except to preserve non-waivable integrity checks;
- refactor unrelated checker code;
- add stronger lexical heuristics for cache, token, handle, journal, protocol, manager, or session fields;
- create another ownership abstraction.

## Verification

Run in this order:

```bash
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
git diff --check
git diff --stat
git diff -- tools/web-change-policy.json
git diff -- gbdraw/web
git status --short --untracked-files=all
```

Required results:

```text
architecture-contract tests pass
new contraction cases pass through explicit base/head revision reads
same-PR policy plus runtime expansion fails against the base allowlist
active owner and importer contractions remain non-waivable
deleting an empty base-policy capability/import-target key fails
working-tree Web change-budget check passes
tools/web-change-policy.json has no diff
gbdraw/web has no diff
no dependency or workflow file changed
```

Do not run browser, Playwright, Pyodide, rendering, Gallery, or full Python suites for this guard-only change unless a focused failure demonstrates that they are unexpectedly relevant.

## Diff audit before handoff

Audit the final diff separately for:

1. checker behavior;
2. temporary-repository tests;
3. policy documentation.

Confirm that:

```text
base-policy expansion protection remains intact
proposed-policy contraction coverage is enforced
safe contractions pass
unsafe contractions fail
base-policy capability and import-target keys remain present
runtime/guard separation remains intact
real policy contents remain unchanged
```

Prefer deletion or replacement of the unconditional contraction failure over layering an exception around it.

## Completion handoff

Report:

1. the exact base SHA used;
2. changed files;
3. the previous behavior and new behavior;
4. the positive and negative contraction cases and policy-key integrity case added;
5. verification commands and results;
6. confirmation that no Web runtime, workflow, dependency, or real policy content changed;
7. final diff statistics;
8. any residual risk.

Provide this proposed commit title and a concise English commit summary:

```text
Allow Web ownership policy to ratchet downward
```

Do not continue into the Feature visibility runtime migration in this session.
This pull request is complete when the guard can safely accept a later
path-entry contraction while preserving the base policy's keys and all existing
expansion protections.
