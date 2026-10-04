# Session 02 instruction prompt: Web controller and local review draft

## Mission

Refactor the Similarity Group alignment controller into one explicit operation
state machine. It must give immediate busy feedback, create an immediately
applicable preselected review draft from Python-owned facts, edit that draft
locally without Worker calls, and submit one atomic validation/regeneration job
on Apply. Do not finish visual layout in this session.

## Branch and prerequisites

Use only `issue-586-similarity-alignment-ux-20260924`. Fetch the remote, verify
the branch/upstream, and confirm Session 00 authority and Session 01 schema-2
domain work are ancestors. Preserve unrelated changes.

When this prompt is supplied as the session request, it authorizes one focused
Session 02 commit and push to the same-named branch. It does not authorize a PR,
merge, release, or direct push to `dev`/`main`.

## Read before editing

Read `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`, the master plan, the
active Product decisions, both architecture/product ratchets, Session 01's diff
and ledger evidence, and the current controller/tests. Trace the existing
Worker helper, Result admission, active-plan repair, History, Session draft,
and stale-operation paths before changing state.

Likely owners include:

- `gbdraw/web/js/app/similarity-alignment.js`;
- `gbdraw/web/js/app/python-helpers.js`;
- existing popup and drawer action adapters;
- existing Worker/generation controller and current Result admission owner;
- `tests/web/similarity-alignment-actions.test.mjs` and related state tests.

## Controller contract

Keep one controller and one draft. Use a small explicit state machine equivalent
to:

```text
idle -> resolving -> reviewing -> applying -> idle
                     ^             |
                     |--- failure--|
```

Required behavior:

1. Capture the exact initiating reference identity.
2. Enter `resolving` synchronously before awaiting the Python helper. Expose one
   derived `busy` state to every initiating control and reject duplicate starts.
3. Admit only the current operation token and compatible committed-input
   snapshot. A stale or superseded completion is a no-op.
4. Convert the Python projection into one local draft row for every displayed
   non-reference record:
   - only-usable/direct-RBH outcomes are selected with their Python reason;
   - ambiguous rows select the Python recommendation with its reason;
   - missing/unusable rows are explicitly unchanged;
   - orientation policy defaults to `preserve`.
5. Enter `reviewing` for every successful initial resolution, including when no
   ambiguity exists. Never auto-apply.
6. Candidate Select, Skip, orientation changes, and canvas candidate clicks
   mutate the same local draft function and start no Worker.
7. `Apply` sends all target rows as one explicit batch. It invokes the shared
   Python resolver once for validation and then the canonical generation path.
8. Admit the new diagram/plan/history only after the complete operation
   succeeds. Do not partially apply record transforms.
9. An Apply validation or generation failure returns to `reviewing` with the
   same editable draft and actionable error.
10. Cancel closes the draft without changing active plan, current diagram, last
    successful Result, or History.
11. Existing active-plan regeneration and repair use the same plan/admission
    owners. Do not add a hidden shortcut, alternate Worker, or second history
    stack.

Keep UI rendering passive: templates/buttons consume state and dispatch
controller actions. They must not determine recommendations, eligibility, or
effective orientation. Popup and drawer must call the same start operation.

## Tests

Use deterministic fakes around the Worker boundary and assert job counts.
Cover at least:

- busy state is observable before the first awaited helper result;
- repeated invocation while resolving/applying starts one job;
- exact reference survives popup and drawer initiation;
- no-ambiguity and ambiguity responses both enter review;
- all initial selections and reasons come from Python data;
- local candidate, Skip, orientation, and canvas edits create zero jobs;
- initial draft is immediately applicable;
- Apply produces one validation/generation request with complete explicit rows;
- Apply failure keeps draft and prior Result/History;
- Cancel changes no committed state;
- stale/superseded resolve or Apply completion changes no committed state;
- successful Apply creates one history action and preserves existing active-plan
  repair/regeneration semantics.

Run the focused Web unit tests, including at minimum:

```bash
node --test \
  tests/web/similarity-alignment-actions.test.mjs \
  tests/web/right-drawer.test.mjs \
  tests/web/diagram-generation-worker.test.mjs \
  tests/web/history-config-restore.test.mjs \
  tests/web/session-draft-authority.test.mjs

node tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
git diff --check
```

Adjust the exact test list to changed owners, never by weakening existing
assertions.

## Architecture evidence, commit, and handoff

Record controller states, Worker job-count evidence, stale-result owner, and
removed parallel paths in the master ledger. Confirm popup/drawer are adapters,
not semantic owners. Review production/test diffs separately.

Use an English commit title such as:

```text
Build the similarity alignment review draft
```

Push the branch and report commit, tests, architecture/change-budget results,
and the presentational work intentionally left to Session 03.
