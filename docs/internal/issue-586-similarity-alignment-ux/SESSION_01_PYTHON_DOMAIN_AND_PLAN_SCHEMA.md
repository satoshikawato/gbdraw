# Session 01 instruction prompt: Python domain owner and plan schema 2

## Mission

Implement the shared Python semantics required by Issue #586: deterministic,
disclosed recommendations for ambiguous Similarity Group members; explicit
per-record orientation intent; and a fully resolved schema-2 alignment plan.
Keep Python as the single biological and transform-semantic owner. Do not build
the Web palette in this session.

## Required branch and authority gate

Use only:

```text
issue-586-similarity-alignment-ux-20260924
```

Fetch `origin` and verify that the Session 00 authority replacements are merged
into `origin/dev`. Merge that `origin/dev` into the implementation branch if it
is not already an ancestor. Do not continue if `PD-OI-026`, `PD-OI-027`,
`PD-OI-031`, and `PD-OI-034` revision 2 are absent from the base authority.

When this prompt is supplied as the session request, it authorizes one focused
Session 01 commit and a push to the same-named implementation branch. It does
not authorize a PR, merge, release, or direct push to `dev`/`main`.

## Read before editing

Read completely:

1. `AGENTS.md`, `CLAUDE.md`, and `gbdraw/web/CLAUDE.md`;
2. the master implementation plan in this directory;
3. the four active revision-2 Product decisions and the still-active plan
   lifecycle, reset/history, Session compatibility, and canvas decisions;
4. `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`;
5. `docs/internal/PRODUCT_IMPACT_RATCHET.md`;
6. current owners in `gbdraw/layout/similarity_alignment.py`,
   `gbdraw/api/record_planning.py`, `gbdraw/diagrams/linear/assemble.py`,
   `gbdraw/web_support/similarity_alignment.py`, and
   `gbdraw/session_request_codec.py`;
7. their focused tests.

Inspect the working tree before editing and preserve unrelated changes.

## Pre-edit persisted-format audit

The plan was written when similarity plan schema 1 existed on `dev` but not in
first-parent `main` history or any release tag. Repeat that audit using Git
history and release tags. Inventory schema-1 examples and fixtures, including
Gallery Sessions and `tests/test_inputs`.

- If schema 1 remains unreleased, replace it in place with schema 2. Remove the
  superseded current reader, tests, fixtures, and documentation rather than
  adding a compatibility branch.
- If schema 1 is now released or covered by an accepted compatibility promise,
  stop before schema edits. Report the evidence and prepare a bounded
  reader-only migration proposal. Do not silently create dual current schemas.

Record the result in the master plan ledger.

## Domain model to implement

Keep the shared resolver as the only owner of all rules below.

### Selection and recommendation

Retain current final selection precedence:

1. validated explicit Select or Skip;
2. only usable candidate;
3. one distinct direct RBH candidate to the exact reference;
4. unresolved ambiguity.

For an unresolved ambiguity, return one transient recommendation and its finite
reason:

1. exactly one representative candidate -> `unique representative`;
2. otherwise candidate 1 in the existing canonical identity-based order ->
   `deterministic candidate 1`.

Also expose the existing `only usable candidate` and `unique direct RBH`
rationales as review reasons. A recommendation is not a plan decision. It must
remain replaceable, skippable, and subject to final validation.

Prove with tests that recommendation is independent of viewport, scroll,
visibility, ribbon geometry, confidence score, supporting-edge count, and
multi-hop evidence. Do not add a score or generalized ranking abstraction.

### Per-record orientation

Replace global orientation mode in the current typed representation with one
small per-record enum, conceptually `preserve` and `match_reference`.

- preservation is the default;
- reverse a whole target only for explicit `match_reference` and known opposite
  displayed anchor strands;
- same known strands or either unknown strand preserve current orientation;
- reference, missing, unusable, and skipped records preserve orientation;
- target Y positions never change;
- anchor-center X alignment is computed after effective orientation and remains
  exact and idempotent;
- retain both requested policy and effective reverse-complement outcome in the
  final plan.

Extend the existing record choice rather than adding a second request model if
that keeps one clear contract. Enforce complete record coverage and exact stable
record/feature identity. Reject invalid policy/choice combinations explicitly.

### Plan schema 2

Advance `SIMILARITY_ALIGNMENT_PLAN_SCHEMA` to 2. The current writer and typed
plan use only per-record orientation policy plus effective outcome; remove the
global mode from the current representation and its normal consumers. Update
exports, renderer/materializer inputs, Python codec, and Web JSON projection as
needed, but do not implement Web state or markup here.

The initial Web-resolution projection must contain enough Python-owned facts for
one row per target: candidate identities/details, automatic outcome or
ambiguity, recommended choice, recommendation reason, and orientation facts.
The adapter only serializes; it must not rank or decide orientation.

Preserve strict CLI behavior: its exact reference remains position-preserving,
and unresolved ambiguity remains an actionable error rather than accepting a
recommendation automatically. Typed Python accepts only fully resolved plans.

## Tests

Add or update tests for at least:

- exact reference preservation and stable identity validation;
- only-usable and unique-direct-RBH automatic outcomes;
- unique-representative recommendation;
- canonical candidate-1 recommendation with record/candidate reorder
  invariance;
- explicit replacement of a recommendation and explicit Skip;
- zero/missing/unusable candidate preservation;
- every orientation-policy/strand combination;
- effective source-relative `rev`, text readability, crop coordinates, exact
  post-reversal centers, and repeated-apply idempotence;
- schema-2 exact round trip and rejection of malformed/incomplete plans;
- strict CLI ambiguity rejection;
- adapter output containing reasons/facts without independent policy.

Run at minimum:

```bash
python -m pytest \
  tests/test_similarity_alignment.py \
  tests/test_similarity_alignment_rendering.py \
  tests/test_similarity_alignment_web_adapter.py -v

python -m pytest \
  tests/test_session_request_codec.py \
  tests/test_session_compat.py \
  tests/test_api_session.py -v

ruff check gbdraw/ tests/test_similarity_alignment.py \
  tests/test_similarity_alignment_rendering.py \
  tests/test_similarity_alignment_web_adapter.py

git diff --check
```

Do not regenerate public Gallery artifacts yet; Session 04 owns coordinated
artifact regeneration.

## Architecture evidence, commit, and handoff

Document the semantic owner, canonical Python-to-adapter path, schema audit, and
removed superseded paths in the master ledger. Review production and test diffs
separately. There must be no JavaScript recommendation algorithm and no second
plan model.

Use an English commit title such as:

```text
Add per-record similarity alignment plans
```

Push the implementation branch and report commit, changed files, tests, schema
audit result, and remaining work for Session 02.
