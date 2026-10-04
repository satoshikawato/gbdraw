# Session 04 instruction prompt: persistence, lifecycle, artifacts, and docs

## Mission

Complete the cross-workflow integration for the schema-2 per-record alignment
plan. Preserve Session, regeneration, active-plan repair, Reset Align,
Undo/Redo, record reorder, failure/retry, CLI, and typed Python meanings.
Regenerate branch-owned artifacts through their owners and update user
documentation truthfully.

## Branch and prerequisites

Use only `issue-586-similarity-alignment-ux-20260924`. Fetch, verify branch and
upstream, confirm Sessions 00–03 are ancestors, inspect status, and preserve
unrelated changes.

When this prompt is supplied as the session request, it authorizes one focused
Session 04 commit and push to the same-named branch. It does not authorize a PR,
merge, release, or direct push to `dev`/`main`.

## Read before editing

Read the repository guidance, master plan, active Product decisions,
architecture/product ratchets, all prior session diffs/ledger evidence, and the
current owners for:

- Python/Web Session request codecs and authority;
- released Session compatibility;
- active alignment plan regeneration and repair;
- Reset Align and immediate pre-align baseline;
- Undo/Redo and current Result admission;
- record reorder/stable key handling;
- Gallery Session generation/publication;
- Web, CLI, typed Python, Session-compatibility, and release-note documentation.

If editing Gallery tutorial screenshots, captions, alt text, tutorial JSON, or
operation registers, load and follow the repository's
`web-gallery-screenshot-maintenance` skill before acting. Do not use that skill
merely because a generated Session file changes.

## Persistence and lifecycle requirements

- Current writers emit schema 2 only, with exact per-record requested
  orientation and effective outcome.
- Current readers reject malformed, partial, mismatched, or unsupported current
  plans explicitly.
- Retain only compatibility paths justified by released first-parent/tag
  evidence. Do not resurrect the unreleased schema-1 format.
- A loaded current Session regenerates through the canonical plan/render path
  without asking the resolver to guess an anchor.
- A successful new Apply replaces the active plan and immediate pre-align
  baseline as one history operation.
- Reset Align restores the immediate pre-align baseline already defined by the
  active Product Contract.
- Undo/Redo restores the plan, translations, requested orientation, effective
  orientation, diagram, and inspector state together.
- Record reorder resolves by stable record key and biological feature identity;
  array position is not authority.
- Active-plan repair after a supported edit retains one canonical plan owner.
  An invalid/stale repair cannot overwrite the last successful Result.
- Apply validation/generation failure keeps the local review draft. Cancel and
  stale/superseded work commit nothing.
- Session save never includes floating palette state, candidate markers,
  transient recommendation badges, or other preview overlays.

Do not bump a broader Session/request schema merely because the nested alignment
plan schema changed unless the exact codec contract requires it. If a bump is
required, document why and update both Python and Web canonical writers/readers
together; do not introduce dual normal writers.

## Generated artifacts

Inventory every checked-in occurrence of similarity plan schema 1. Known
candidates include:

- `gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json`;
- Session fixtures under `tests/test_inputs/`;
- executable documentation fixtures.

Regenerate Gallery-owned Session JSON with the existing owner, such as
`tools/refresh_gallery_sessions.py` or the narrower current supported command.
Never hand-edit generator-owned JSON. Review the session, manifest, and any
derived media diffs independently; retain only outputs required by the owner.
Do not regenerate or edit `examples/gbdraw_social_preview.png`.

## Documentation

Update existing owners rather than adding a parallel guide:

- Web reference: one `Align…`, always-open review, recommendation reasons,
  per-record Match reference direction, local edits, Apply/retry/Cancel;
- typed Python: schema-2 fully resolved plan and per-record policy;
- CLI: exact reference, position preservation, and strict ambiguity rejection;
- Session/request compatibility: current schema 2 and only proven released
  legacy readers;
- release notes and navigation where existing contracts require them.

State explicitly that recommendations are convenience heuristics and that
Collinear alignment controls, anchor TSV, scored inference, support-count
ranking, and multi-hop automatic selection are unsupported.

Keep literal code examples executable. Update documentation contract tests
rather than weakening them.

## Tests

Cover at least:

- Python/Web schema-2 exact round trips;
- current malformed/partial plan rejection;
- proven released legacy Session load and current-format save;
- regeneration and active-plan repair;
- immediate-baseline Reset Align;
- Undo/Redo across selection and per-record orientation;
- record reorder with stable identity;
- failure/cancel/stale isolation and review-draft retry;
- Gallery Session regeneration/publication parity;
- literal typed-Python example and CLI reference text.

Run the focused suites appropriate to changed owners, including:

```bash
python -m pytest \
  tests/test_session_request_codec.py \
  tests/test_session_compat.py \
  tests/test_api_session.py \
  tests/test_session_io.py \
  tests/test_documentation_contracts.py \
  tests/test_documentation_reference_contracts.py -v

node --test \
  tests/web/session-request.test.mjs \
  tests/web/session-authority.test.mjs \
  tests/web/session-draft-authority.test.mjs \
  tests/web/history-config-restore.test.mjs \
  tests/web/gallery-session-migration.test.mjs \
  tests/web/gallery-session-publication.test.mjs

python tools/update_cli_reference_help.py --check
node tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
git diff --check
```

Run the repository's supported Gallery reproduction/parity command for every
regenerated artifact.

## Architecture evidence, commit, and handoff

Record current writer/reader paths, exact released compatibility evidence,
schema-1 removals, generated commands/diffs, and lifecycle test results in the
ledger. Review production, tests, docs, and generated artifacts separately.

Use an English commit title such as:

```text
Persist per-record similarity alignment choices
```

Push the branch and report commit, tests, generated artifacts, documentation,
and any residual concern for Session 05.
