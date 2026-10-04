# Issue #561 — S03 Python rendering, API, and CLI instruction prompt

```text
Implement Session S03 of gbdraw Issue #561: apply the typed plan through the
existing Python record-transform and Linear renderer paths, expose it through
the typed Python API, and make CLI resolution strict and non-interactive.

Repository and branch:
- Use issue-561-similarity-alignment.
- Confirm merged Product authority and completed S01/S02 ledger entries.
- Preserve unrelated changes and do not create another runtime branch.

Read before editing:
1. AGENTS.md and CLAUDE.md
2. docs/internal/ISSUE_561_SIMILARITY_ALIGNMENT_MASTER_PLAN_2026-09-22.md
3. accepted similarity-alignment decisions in OPTION_INTEGRITY_PRODUCT_CONTRACT.md
4. S01/S02 code and tests
5. gbdraw/api/record_planning.py, request_render.py, requests.py, options.py,
   diagram.py, and package exports
6. gbdraw/layout/record_coordinates.py
7. gbdraw/diagrams/linear/orthogroup_alignment.py, assemble.py, builders.py,
   positioning.py, and precalc.py
8. gbdraw/linear.py and CLI tests/documentation
9. PD-OI-016, PD-OI-018, and canonical request contracts

Implement rendering:
1. Materialize an effective per-record orientation at request-planning time.
   Base orientation stays in RecordPresentation. The plan may override the
   effective boolean for aligned targets; do not add another source orientation
   flag or reverse the same record twice.
2. Project every anchor through the existing RecordDisplayTransform after the
   effective orientation is known.
3. Compute absolute target translations from the reference world X and the
   target post-transform center. Preserve every Y translation. Keep the
   reference base translation unchanged. Keep skipped records' base X/Y and
   orientation unchanged.
4. Feed the same final translations to records, definitions, tracks, rulers,
   labels, annotations, comparison geometry, bounds, and composition metadata.
   Do not apply an SVG-only feature-group patch after layout.
5. Remove the current rejection of several records in one row and support each
   displayed record independently.
6. Recalculate canvas left/right extents from final record bounds and avoid
   clipping under normalize-length, center alignment, negative translations,
   and same-row placement.
7. Prove idempotence: repeated planning/rendering of the same request produces
   the same translations and geometry.
8. Legacy string input may reach only the old-schema adapter. Remove normal
   current renderer selection of representatives or scores after the adapter is
   covered.

Implement public surfaces:
1. Export the typed plan through the existing gbdraw.api boundary without
   exposing renderer implementation classes.
2. Replace the current typed public str field with the plan. Public invalid
   values fail before rendering.
3. Keep --align_orthogroup_feature as the initial CLI spelling unless the
   current CLI contract requires an alias. It now identifies an exact feature
   or protein reference, not permission to select by group ID.
4. CLI calls the shared resolver. It succeeds without prompting only when each
   record is resolved uniquely or has no usable member. Ambiguity reports the
   affected record and exact candidate identifiers and exits with an actionable
   error.
5. CLI/Python alignment reuses available group evidence and does not invoke a
   new LOSATP job merely to resolve an active/supplied plan.
6. Do not add an anchor TSV, Collinear UI, smart alignment, or another plan
   builder.

Testing:
- formula-level center alignment, negative/positive offsets, normalize-length,
  center-alignment, crop, and reverse display;
- strand matrix: same, opposite, unknown reference, unknown target;
- position-only orientation preservation and position-and-orientation whole
  record reverse;
- reference fixed, skipped/missing fixed, Y preserved;
- same-row multi-record, definitions/tracks/ruler/comparison agreement, canvas
  bounds, and idempotence;
- typed API success and validation;
- CLI exact unique success, group/ambiguous error, missing-member preservation,
  and no interactive prompt;
- instrumentation or mocks proving no new LOSATP invocation;
- legacy old-schema output parity through the isolated adapter.

Architecture:
- Keep the existing RecordDisplayTransform and renderer as coordinate owners.
- The shared resolver selects anchors; renderer never ranks candidates.
- Remove superseded normal implicit selection paths in the same change.
- Record concise owner/path evidence in the master-plan ledger.

Finish:
- Run focused Python, CLI, request render, multi-record layout, orthogroup, and
  output-comparison tests. Treat reference outputs as read-only unless a
  reviewed intentional geometry change requires regeneration.
- Run ruff on gbdraw and review production/tests/generated diffs separately.
- Update ledger section 16 and state whether S04 may start.
- Provide an English proposed commit title and summary; do not push or create a
  PR unless separately authorized.
```
