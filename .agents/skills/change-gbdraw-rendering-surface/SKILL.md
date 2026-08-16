---
name: change-gbdraw-rendering-surface
description: Implement or review a gbdraw rendering option across every surface that owns or persists it. Use when adding, changing, renaming, or removing diagram settings, track behavior, labels, legends, geometry, comparison data, output topology, Python APIs, CLI flags, Web controls, typed render requests, saved sessions, Gallery examples, or reference SVGs.
---

# Change a gbdraw rendering surface

Read `AGENTS.md` and `CLAUDE.md` first. Read `gbdraw/web/CLAUDE.md` when the Web app, session projection, browser worker, or Gallery is in scope. Use `$maintain-python-api` for public API contract details and `$web-gallery-screenshot-maintenance` for Gallery media work.

## Map the current owner

Before editing, trace the concept through the surfaces it actually uses:

| Surface | Inspect |
|---|---|
| Typed core | config model, request model, validator, planner, prepared value |
| Python | package-root facade, `gbdraw.api`, option dataclass, export contract |
| CLI | parser, normalization, shared request construction, help text |
| Web | reactive state, visible control, canonical request projection, worker payload |
| Persistence | current writer, supported readers, migration evidence, fixtures |
| Rendering | Circular and Linear planners, drawers/groups, SVG semantics, export |
| Documentation | CLI reference, Python recipes, Gallery tutorial, release notes |

Record the default, non-default value, validation owner, final consumer, and persistence field. Omit a surface only after confirming that the concept does not apply there.

## Keep one behavior path

- Put shared behavior in the typed core or planner that owns it.
- Keep adapters limited to translating surface-specific input.
- Project Web state once in `gbdraw/web/js/services/session-request.js`.
- Keep fresh generation, session replay, Python, and CLI paths convergent.
- Do not add an argv-shaped Web contract, parallel config model, or duplicate renderer.
- Add an abstraction only when it unifies at least two active paths and remove the superseded paths in the same change.

When a setting controls topology, preserve explicit empty positions and record alignment. When it controls output naming, preflight every target from the same plan used for rendering.

## Implement the lifecycle

1. Add or change the typed model and default.
2. Normalize and validate once at its owner.
3. Carry the resolved value through the prepared request.
4. Consume it in the real rendering owner for each applicable mode.
5. Update public API and CLI adapters without duplicating policy.
6. Bind the Web control and project it into the current render request.
7. Persist user-owned state in the current session writer and update only evidence-backed readers.
8. Update help, recipes, and Gallery material that users rely on.
9. Remove retired fields, branches, documentation, and branch-only migrations.

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
