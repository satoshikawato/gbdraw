> Reference playbook. Use only when a user or task-specific instruction explicitly requests this playbook. Do not auto-route to it from file paths, changed modules, or task category.

# Change a gbdraw rendering surface

Read `AGENTS.md` and `CLAUDE.md` first. Read `gbdraw/web/CLAUDE.md` when the Web app, session projection, browser worker, or Gallery is in scope. Use the [maintain Python API playbook](../maintain-python-api/PLAYBOOK.md) for public API contract details and the [Web Gallery screenshot maintenance playbook](../web-gallery-screenshot-maintenance/PLAYBOOK.md) for Gallery media work.

## Decide whether this playbook applies

Use this playbook when the requested outcome changes a shared or persisted rendering contract. Typical triggers are:

- a typed schema, option, default, or validation rule shared by more than one public surface;
- the same option semantics crossing Python, CLI, Web request construction, or the renderer;
- session serialization or compatibility behavior;
- output topology, geometry, SVG semantics, or export behavior shared by multiple modes.

A local bug can touch several implementation and test files without changing a shared contract. If one surface is merely wired to the wrong existing handler, projects an existing value incorrectly, or fails to replay an already-defined behavior, inspect the adjacent boundary and keep the fix local. Do not add schema versions, migrations, public APIs, Gallery assets, or reference-output changes to complete a surface matrix.

If evidence found during a local fix shows that the shared contract itself is inconsistent, state that evidence and apply the remaining workflow only to the owners that actually participate. A surface audit is not authorization to expand the requested behavior.

## Map the current owner

After confirming that the playbook applies, trace the concept through the surfaces it actually uses:

| Surface | Inspect |
|---|---|
| Typed core | config model, request model, validator, planner, prepared value |
| Python | package-root facade, `gbdraw.api`, option dataclass, export contract |
| CLI | parser, normalization, shared request construction, help text |
| Web | reactive state, visible control, canonical request projection, worker payload |
| Persistence | current writer, supported readers, migration evidence, fixtures |
| Rendering | Circular and Linear planners, drawers/groups, SVG semantics, export |
| Documentation | CLI reference, Python recipes, Gallery tutorial, release notes |

Record the default, non-default value, validation owner, final consumer, and persistence field. Mark non-applicable surfaces explicitly; do not create an owner merely to fill the matrix.

## Keep one behavior path

- Put shared behavior in the typed core or planner that owns it.
- Keep adapters limited to translating surface-specific input.
- Project Web state once in `gbdraw/web/js/services/session-request.js`.
- Keep fresh generation, session replay, Python, and CLI paths convergent.
- Do not add an argv-shaped Web contract, parallel config model, or duplicate renderer.
- Add an abstraction only when it unifies at least two active paths and remove the superseded paths in the same change.

When a setting controls topology, preserve explicit empty positions and record alignment. When it controls output naming, preflight every target from the same plan used for rendering.

## Implement the lifecycle

Apply only the steps owned by the changed contract:

1. Add or change the typed model and default when the shared model changes.
2. Normalize and validate once at its existing owner.
3. Carry the resolved value through the prepared request when rendering consumes it.
4. Consume it in the real rendering owner for each applicable mode.
5. Update public API and CLI adapters only when they expose the contract.
6. Bind or update the Web control and request projection only when the Web surface exposes it.
7. Update the current session writer and evidence-backed readers only when persisted semantics change.
8. Update help, recipes, Gallery material, or reference outputs only when users would otherwise receive stale evidence.
9. Remove paths superseded by the complete implementation; do not remove a partially correct path merely to match a planned file set or change-size target.

Defaults must preserve previous behavior unless the request explicitly changes the default.

## Test the boundary

Start with a matrix of applicable modes and surfaces. Cover:

- unchanged default behavior;
- one non-default value reaching the final consumer;
- validation and error semantics;
- Python and CLI agreement;
- Web request and session round trips;
- Circular and Linear behavior when the setting is shared;
- output naming and grouping when topology changes;
- SVG semantic IDs or metadata when interaction depends on them.

Run focused tests first, then the relevant broader gate. Use Python Playwright when Node Playwright is unavailable. Rerun Chromium with sandbox approval when the local sandbox blocks it.

## Handle visual and generated evidence

Treat `tests/reference_outputs/` as read-only until a reviewed geometry change requires regeneration. Compare output first. If regeneration is intentional, use the repository command, inspect the SVG diff, and rerun comparison tests.

Gallery sessions, source SVGs, examples, thumbnails, and `examples.json` are generator-owned. Regenerate them from declared inputs. Preserve the labels, legend, metadata, tracks, and comparison context that make a public example useful. Never replace a Gallery-quality example with a minimal test diagram.

Render and inspect the final artifact at a readable scale. Check clipping, label collisions, legend placement, quantitative tracks, comparison ribbons, and interactive metadata.

## Finish

Audit production code, tests, documentation, and generated output as separate diffs. Report the surface matrix covered, commands run, browser cases checked, reference-output status, and any unsupported surface. Provide the requested commit title and summary.

When changing this playbook's routing or execution rules, verify the lightweight scenarios in [EVALS.md](EVALS.md).
