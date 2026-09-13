# Chloroplast feature placement regression

Base: `991f8eef` (`origin/dev`, fetched on 2026-09-09).

The Session 00 delivery reuses this fix on the same current `origin/dev` base.
See [delivery verification](DELIVERY.md) for provenance and the fresh-browser
restoration check.

The published `tobacco-chloroplast.gbdraw-session.json` reproduced both reported
validation errors and the overlapping legend in a local Chromium instance.
The screenshots below were captured during that reproduction; none were copied
from the user's attachments.

## Reproduction and causes

1. Load the Gallery tobacco chloroplast session and Generate. Separate Strands
   is enabled. The old feature placement capability list permits only Auto/Main.
2. Clear Separate Strands and Generate. Set multipart CDS features `clpP`,
   `rpl16`, and `rpoC1` to Inward lane 1. Generate with Label Mode Both.
   The renderer reports `required_additional_px=12.943`.
3. Change Label Mode to Out and Generate. The legend moves onto the circle.
4. Enable Separate Strands and Generate. The retained directional placements
   are rejected as unsupported in the split slot.

The failures have five related causes:

- Circular split slots unnecessarily reject directional placement whenever
  strands are separate. A saved direction therefore also makes strand toggling
  fail. Lane 1 now means one lane beyond the corresponding nominal strand lane;
  the override remains bound to the complete biological feature.
- Adding an inward lane shifts the combined-strand Main lane by half a lane
  step. Main now stays at the slot anchor. Separated lanes also retain their
  levels when another lane is empty, and resolved slot direction takes priority
  over a preset hint.
- A clockwise-only label packing sweep can cross leaders whose features lie
  at different radii. This was repeatedly reported as a radius deficit. When
  leaders collide, the packer balances the clockwise and counterclockwise
  solutions before asking for more space. It checks the resulting text and
  leader geometry through the existing collision checks.
- Corner-legend placement treats external labels as obstacles but omits the
  record body. With outer labels only, it can place a legend inside the circle.
  A rendered legend now also avoids the record/track/definition bounds. The
  same obstacle is serialized for the Web editor's later layout reflow.
- Session import/export drops the resolved track geometry, including the
  available placement targets. The session service now saves and restores this
  existing artifact field and clears it when resetting the session baseline.
  The refreshed Gallery session offers placement immediately after import;
  saved sessions retain it after reloading. Older sessions that do not contain
  this metadata acquire current targets when Generate is run.

Canvas dimensions still follow visible content and the space needed for labels
and the legend. The fix preserves nominal feature lanes and prevents overlap;
it does not freeze canvas dimensions across different label modes.

| Before | After |
| --- | --- |
| [Inner-label validation failure](before-inner-label-error.png) | [Combined strands with inner labels](after-inner-labels.png) |
| [Legend over the circle](before-legend-overlap.png) | [Outer labels with a clear legend](after-outer-labels.png) |
| [Strand-toggle validation failure](before-strand-toggle-error.png) | [Separated placement after session restore](after-restored-session.png) |

## Verification

Final Python gate: `pytest tests/ -q -m "not slow" -n 4` completed with
**3969 passed, 17 skipped**. The 16 tracked SVG comparisons pass without changing
their references. The dedicated documentation/placement suite passes all 393
cases, and the session/history JavaScript checks pass all 15 test files. Ruff
and `git diff --check` also pass.
All 14 browser scenarios pass: 12 passed together, then the two existing
Circular/Linear saved-draft scenarios passed after replacing their old
"Main is disabled until Generate" expectation with the restored behavior.
The scenarios still check Undo, saved drafts, and regeneration afterwards.

`tests/test_chloroplast_placement.py` uses the published session, including its
color rules, labels, region annotations, GC track, and multipart source features.
Its eight cases cover combined/separate strands, inner/outer label modes, and
placements both with and against the biological strand. They assert physical
lane centers, unchanged nominal lanes, label anchors, no text/leader collisions,
and legend clearance. All eight fail against the four untouched base production
modules and pass with the fix.

Restoring each production module independently to its base version also makes
the new regression test fail:

| Restored module | Failing cases out of eight |
| --- | --- |
| Feature placement | 4 |
| Circular lane geometry | 8 |
| Radial label packing | 2 |
| Circular composition obstacles | 8 |

`tests/web/chloroplast-placement.playwright.spec.js` loads the real session in a
cold browser context with external requests blocked. It edits placement controls
immediately after import and again immediately after reloading a saved session,
toggles strand/label settings, regenerates repeatedly, checks mounted SVG bounds,
saves and reloads a session, and checks the restored result at desktop and narrow
mobile widths. Its screenshots and SVGs are emitted into the Playwright output
directory. The before screenshots use a 1600 × 1000 viewport at 100% zoom; the
after screenshots use the same viewport at 70% so more of the figure is visible.
Before the session-service fix, this browser test failed at its first Inward
selection: the imported session contained resolved targets but the UI's geometry
state was null.

Run from an environment with this checkout installed, so child CLI processes
also use the candidate code:

```bash
python -m pytest tests/test_chloroplast_placement.py tests/test_feature_placement.py tests/test_circular_radial_layout.py tests/test_circular_radial_labels.py tests/test_circular_composition.py -q
python -m pytest tests/test_output_comparison.py::TestOutputComparison -q
python tools/prepare_browser_wheel.py --no-build-isolation
node node_modules/@playwright/test/cli.js test tests/web/chloroplast-placement.playwright.spec.js
```

The existing Circular and Linear rotation/placement browser smoke tests also
exercise downloaded Source recipe and Exact replay commands. Use a dedicated
virtual environment: the workspace's global editable installation points at a
different, older checkout and does not implement the current placement CLI.

The Gallery session was regenerated with
`python tools/refresh_gallery_sessions.py --session tobacco-chloroplast`.
Its render request, source resources, and editor state are unchanged; only the
result and resolved run metadata change. The Gallery source, interactive SVG,
and artifact manifest were regenerated through that same command.

The affected documentation outputs were regenerated through the existing
CLI recipes `T-CLI-05`, `T-CLI-06` and Python recipes `T-PY-02`, `T-PY-11`,
`H-PY-03`. For example:

```bash
python docs/recipes/run_cli_scenarios.py --scenario T-CLI-06
python docs/recipes/run_python_scenarios.py --scenario T-PY-02
```

The Chloroplast and quantitative-map outputs retain identical rendered pixels
at 1200 px width; their differences are composition metadata and floating-point
serialization. `H-PY-03` also changes visible lane geometry: its explicit custom
slot direction now takes precedence over the separate-strand preset. Both
versions were rendered and inspected. All records, labels, colors, quantitative
tracks, and annotations remain present. No reference output or owner-maintained
social preview was regenerated.

Documentation progress: framing, existing-page ownership, pinned inputs,
smoke proof, existing execution harnesses, regeneration, retained public pages,
and artifact review are complete. No public page or capture harness was added.

## Scope and ownership review

The user's request explicitly includes moving whole multipart features with
Separate Strands enabled and retaining those choices through strand toggles.
The circular split restriction documented on the base branch is expanded for
that requested combination. One-sided slots and Linear separate-strand slots
retain their existing restrictions. Main, automatic allocation, fixed-placement
conflicts, source identity, biological strand, and session wire formats retain
their meanings. No saved override is silently removed to make a render succeed.

The same owners remain before and after: `features/placement.py` validates and
assigns feature placements; `diagrams/circular/radial_layout.py` measures their
lanes; `labels/circular_radial.py` packs radial labels;
`diagrams/circular/assemble.py` supplies composition obstacles. The Web UI still
consumes Python's supported-target metadata and uses the existing typed render
request/Worker path. `web/js/services/config.js` remains the session import/export
owner and now retains the already-defined run metadata. There is no new owner,
alternate render path, compatibility
reader, schema, or dependency. Rollback is the single implementation diff; no
persisted data migration is required.
