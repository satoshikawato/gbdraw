# Circular Track Slots Root Fix Plan

## Status

Retired on 2026-05-16.

This older root-fix plan described the transitional
`resolve_circular_track_slots()` helper and annulus-style circular slot
geometry. Those APIs were removed by the unified circular Track Slots resolver
migration.

The detailed historical design notes in
`docs/CIRCULAR_TRACK_SLOTS_UNIFIED_RESOLVER_PLAN.md` predate the retired
`strict`/`compress`/`reserve` flags. The current implementation treats circular
Track Slots as the single radial layout input model, uses slot-level `radius`,
`width`, `spacing`, `side`, and `z` fields, auto-compresses implicit inside
numeric/depth slots, and rejects obsolete `ri`/`ro` and `gap`/`gap_after` slot
grammar.
