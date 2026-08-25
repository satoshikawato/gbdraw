# Web GUI historical evidence routing

- **Evidence source:** `GBDRAW_CONCERN_DRIVEN_BEHAVIOR_HARDENING_ARCHITECTURE_REV4_2026-08-22.md`, Section 14
- **origin/dev observed:** `d96247a7f5f6c6918be30398f1410bcd4aa0fe54`
- **Curator:** Codex implementation agent
- **Date:** 2026-08-26
- **Status:** NON-NORMATIVE

This overlay routes historical evidence. It does not define correct behavior,
assert current failure, assert current conformance, certify current test
coverage, or close any concern. Code, tests, issues, commits, and historical
repairs remain evidence until the applicable authority adopts their meaning.

The source extraction contained 86 table rows and 86 unique HREG aliases before
routing. Each alias appears once below. Priority is qualitative (`high`,
`medium`, or `low`); no weighted score is used. `S0 required` means that current
reachable characterization is needed before product or contract work can rely
on the historical evidence.

## Semantic families

- `feature-legend-edit-continuity`: Feature and Legend style intent across direct editing and lifecycle boundaries.
- `active-result-edit-transactions`: edit scope, active artifact updates, and transaction ownership.
- `biological-identity-selection`: stable biological identity and record or Feature targeting.
- `interactive-metadata-sequence-sources`: popup metadata, coordinates, sequence sources, and ordering.
- `track-layout-state`: track and slot intent, occupancy, capability, and resolved placement.
- `render-cardinality-cross-surface`: record cardinality, output order, mode semantics, and surface parity.
- `transient-ui-continuity`: input drafts, focus, dialog timing, selection, scroll, pan, and zoom.
- `visible-bounds-composition`: labels, Legend, ticks, spacing, and canvas composition.
- `scientific-output-semantics`: domain values, strand meaning, biological connectors, and scientific output.
- `session-normalization-recovery`: input normalization, Session migration, restore, and rollback.
- `history-async-artifact-lifecycle`: History boundaries, Generate replacement, async epochs, and settlement.
- `performance-large-artifact-guards`: resource lifetime, Worker use, parsing, copying, and large-data budgets.
- `artifact-provenance-export`: source and runtime artifact provenance, download, and export settlement.
- `insufficiently-specified-evidence`: the source does not establish a defensible reachable context.

## Routing overlay

| HREG | Historical status | Semantic family | Evidence route | Priority | S0 required | Canonical source |
|---|---|---|---|---|---|---|
| `HREG-LEG-001` | insufficient-evidence | `insufficiently-specified-evidence` | `INSUFFICIENT_EVIDENCE` | low | no | Section 14.1, row 1 |
| `HREG-LEG-002` | mapped | `feature-legend-edit-continuity` | `ENGINEERING_GUARD` | medium | no | Section 14.1, row 2 |
| `HREG-LEG-003` | mapped | `feature-legend-edit-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.1, row 3 |
| `HREG-LEG-004` | mapped | `feature-legend-edit-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.1, row 4 |
| `HREG-LEG-005` | mapped | `feature-legend-edit-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.1, row 5 |
| `HREG-LEG-006` | mapped | `feature-legend-edit-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.1, row 6 |
| `HREG-LEG-007` | mapped | `feature-legend-edit-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.1, row 7 |
| `HREG-LEG-008` | mapped | `feature-legend-edit-continuity` | `DERIVABLE_INVARIANT` | high | no | Section 14.1, row 8 |
| `HREG-EDIT-001` | mapped | `active-result-edit-transactions` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.1, row 9 |
| `HREG-EDIT-002` | mapped | `active-result-edit-transactions` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.1, row 10 |
| `HREG-EDIT-003` | mapped | `active-result-edit-transactions` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.1, row 11 |
| `HREG-EDIT-004` | mapped | `active-result-edit-transactions` | `ENGINEERING_GUARD` | medium | no | Section 14.1, row 12 |
| `HREG-ID-001` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 1 |
| `HREG-ID-002` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 2 |
| `HREG-ID-003` | mapped | `biological-identity-selection` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.2, row 3 |
| `HREG-ID-004` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 4 |
| `HREG-ID-005` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 5 |
| `HREG-ID-006` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 6 |
| `HREG-ID-007` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 7 |
| `HREG-ID-008` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 8 |
| `HREG-ID-009` | mapped | `biological-identity-selection` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 9 |
| `HREG-META-001` | mapped | `interactive-metadata-sequence-sources` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 10 |
| `HREG-META-002` | mapped | `interactive-metadata-sequence-sources` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 11 |
| `HREG-META-003` | mapped | `interactive-metadata-sequence-sources` | `DERIVABLE_INVARIANT` | medium | no | Section 14.2, row 12 |
| `HREG-META-004` | mapped | `interactive-metadata-sequence-sources` | `DERIVABLE_INVARIANT` | high | no | Section 14.2, row 13 |
| `HREG-TRACK-001` | mapped | `track-layout-state` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.3, row 1 |
| `HREG-TRACK-002` | mapped | `track-layout-state` | `DERIVABLE_INVARIANT` | high | no | Section 14.3, row 2 |
| `HREG-TRACK-003` | mapped | `track-layout-state` | `DERIVABLE_INVARIANT` | medium | no | Section 14.3, row 3 |
| `HREG-TRACK-004` | mapped | `track-layout-state` | `DERIVABLE_INVARIANT` | high | no | Section 14.3, row 4 |
| `HREG-TRACK-005` | mapped | `track-layout-state` | `ENGINEERING_GUARD` | medium | no | Section 14.3, row 5 |
| `HREG-TRACK-006` | mapped | `track-layout-state` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.3, row 6 |
| `HREG-TRACK-007` | mapped | `track-layout-state` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.3, row 7 |
| `HREG-TRACK-008` | mapped | `track-layout-state` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.3, row 8 |
| `HREG-TRACK-009` | mapped | `track-layout-state` | `DERIVABLE_INVARIANT` | medium | no | Section 14.3, row 9 |
| `HREG-TRACK-010` | mapped | `track-layout-state` | `DERIVABLE_INVARIANT` | medium | no | Section 14.3, row 10 |
| `HREG-TRACK-011` | insufficient-evidence | `insufficiently-specified-evidence` | `INSUFFICIENT_EVIDENCE` | low | no | Section 14.3, row 11 |
| `HREG-TRACK-012` | mapped | `track-layout-state` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.3, row 12 |
| `HREG-RENDER-001` | mapped | `render-cardinality-cross-surface` | `DERIVABLE_INVARIANT` | high | no | Section 14.3, row 13 |
| `HREG-RENDER-002` | mapped | `render-cardinality-cross-surface` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.3, row 14 |
| `HREG-TRACK-013` | umbrella | `track-layout-state` (canonical family) | `DUPLICATE_OR_UMBRELLA` | low | no | Section 14.3, row 15; delegates to the track-state rows |
| `HREG-UI-001` | mapped | `transient-ui-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.4, row 1 |
| `HREG-UI-002` | mapped | `transient-ui-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.4, row 2 |
| `HREG-UI-003` | mapped | `transient-ui-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.4, row 3 |
| `HREG-UI-004` | mapped | `transient-ui-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.4, row 4 |
| `HREG-LAYOUT-001` | mapped | `visible-bounds-composition` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.5, row 1 |
| `HREG-LAYOUT-002` | mapped | `scientific-output-semantics` | `DERIVABLE_INVARIANT` | high | no | Section 14.5, row 2 |
| `HREG-LAYOUT-003` | mapped | `visible-bounds-composition` | `DERIVABLE_INVARIANT` | high | no | Section 14.5, row 3 |
| `HREG-LAYOUT-004` | mapped | `visible-bounds-composition` | `DERIVABLE_INVARIANT` | medium | no | Section 14.5, row 4 |
| `HREG-LAYOUT-005` | mapped | `visible-bounds-composition` | `DERIVABLE_INVARIANT` | low | no | Section 14.5, row 5 |
| `HREG-LAYOUT-006` | mapped | `visible-bounds-composition` | `DERIVABLE_INVARIANT` | medium | no | Section 14.5, row 6 |
| `HREG-LAYOUT-007` | mapped | `visible-bounds-composition` | `DERIVABLE_INVARIANT` | medium | no | Section 14.5, row 7 |
| `HREG-LAYOUT-008` | mapped | `visible-bounds-composition` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.5, row 8 |
| `HREG-LAYOUT-009` | mapped | `feature-legend-edit-continuity` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.5, row 9 |
| `HREG-LAYOUT-010` | mapped | `visible-bounds-composition` | `ENGINEERING_GUARD` | medium | no | Section 14.5, row 10 |
| `HREG-LAYOUT-011` | mapped | `scientific-output-semantics` | `DERIVABLE_INVARIANT` | high | no | Section 14.5, row 11 |
| `HREG-LAYOUT-012` | mapped | `session-normalization-recovery` | `DERIVABLE_INVARIANT` | high | no | Section 14.5, row 12 |
| `HREG-LAYOUT-013` | mapped | `visible-bounds-composition` | `PRODUCT_SEMANTICS_CANDIDATE` | medium | yes | Section 14.5, row 13 |
| `HREG-LAYOUT-014` | mapped | `visible-bounds-composition` | `ENGINEERING_GUARD` | medium | no | Section 14.5, row 14 |
| `HREG-LAYOUT-015` | mapped | `visible-bounds-composition` | `PRODUCT_SEMANTICS_CANDIDATE` | low | yes | Section 14.5, row 15 |
| `HREG-LAYOUT-016` | insufficient-evidence | `insufficiently-specified-evidence` | `INSUFFICIENT_EVIDENCE` | low | no | Section 14.5, row 16 |
| `HREG-LAYOUT-017` | mapped | `scientific-output-semantics` | `DERIVABLE_INVARIANT` | high | no | Section 14.5, row 17 |
| `HREG-LAYOUT-018` | insufficient-evidence | `insufficiently-specified-evidence` | `INSUFFICIENT_EVIDENCE` | low | no | Section 14.5, row 18 |
| `HREG-LAYOUT-019` | mapped | `scientific-output-semantics` | `DERIVABLE_INVARIANT` | high | no | Section 14.5, row 19 |
| `HREG-LAYOUT-020` | mapped | `scientific-output-semantics` | `DERIVABLE_INVARIANT` | high | no | Section 14.5, row 20 |
| `HREG-SESSION-001` | mapped | `session-normalization-recovery` | `DERIVABLE_INVARIANT` | high | no | Section 14.6, row 1 |
| `HREG-SESSION-002` | mapped | `session-normalization-recovery` | `DERIVABLE_INVARIANT` | high | no | Section 14.6, row 2 |
| `HREG-SESSION-003` | mapped | `session-normalization-recovery` | `ENGINEERING_GUARD` | medium | no | Section 14.6, row 3 |
| `HREG-HISTORY-001` | mapped | `history-async-artifact-lifecycle` | `ENGINEERING_GUARD` | high | no | Section 14.6, row 4 |
| `HREG-HISTORY-002` | mapped | `history-async-artifact-lifecycle` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.6, row 5 |
| `HREG-ASYNC-001` | mapped | `history-async-artifact-lifecycle` | `ENGINEERING_GUARD` | high | no | Section 14.6, row 6 |
| `HREG-ASYNC-002` | mapped | `performance-large-artifact-guards` | `ENGINEERING_GUARD` | high | no | Section 14.6, row 7 |
| `HREG-PERF-001` | mapped | `performance-large-artifact-guards` | `ENGINEERING_GUARD` | high | no | Section 14.6, row 8 |
| `HREG-PERF-002` | mapped | `performance-large-artifact-guards` | `ENGINEERING_GUARD` | medium | no | Section 14.6, row 9 |
| `HREG-PERF-003` | mapped | `performance-large-artifact-guards` | `ENGINEERING_GUARD` | high | no | Section 14.6, row 10 |
| `HREG-PERF-004` | mapped | `performance-large-artifact-guards` | `ENGINEERING_GUARD` | medium | no | Section 14.6, row 11 |
| `HREG-PERF-005` | mapped | `performance-large-artifact-guards` | `ENGINEERING_GUARD` | medium | no | Section 14.6, row 12 |
| `HREG-SOURCE-001` | mapped | `interactive-metadata-sequence-sources` | `DERIVABLE_INVARIANT` | high | no | Section 14.6, row 13 |
| `HREG-EXPORT-001` | mapped | `artifact-provenance-export` | `DERIVABLE_INVARIANT` | medium | no | Section 14.7, row 1 |
| `HREG-EXPORT-002` | mapped | `artifact-provenance-export` | `DERIVABLE_INVARIANT` | high | no | Section 14.7, row 2 |
| `HREG-EXPORT-003` | mapped | `artifact-provenance-export` | `DERIVABLE_INVARIANT` | high | no | Section 14.7, row 3 |
| `HREG-EXPORT-004` | mapped | `artifact-provenance-export` | `PRODUCT_SEMANTICS_CANDIDATE` | high | yes | Section 14.7, row 4 |
| `HREG-SCI-001` | mapped | `scientific-output-semantics` | `DERIVABLE_INVARIANT` | high | no | Section 14.8, row 1 |
| `HREG-SCI-002` | mapped | `scientific-output-semantics` | `DERIVABLE_INVARIANT` | high | no | Section 14.8, row 2 |
| `HREG-SCI-003` | mapped | `scientific-output-semantics` | `ENGINEERING_GUARD` | medium | no | Section 14.8, row 3 |
| `HREG-SCI-004` | duplicate | `interactive-metadata-sequence-sources` (canonical family) | `DUPLICATE_OR_UMBRELLA` | low | no | Section 14.8, row 4; primary evidence is Section 14.2, row 10 |
| `HREG-SCI-005` | duplicate | `scientific-output-semantics` (canonical family) | `DUPLICATE_OR_UMBRELLA` | low | no | Section 14.8, row 5; primary evidence is Section 14.5, row 19 |

## Allowed claim

> Historical regression evidence routing is complete. All 86 aliases are
> represented once as non-normative evidence.

This statement is limited to historical routing. It does not say that any alias
is current, fixed, conforming, protected, normatively correct, or closed.
