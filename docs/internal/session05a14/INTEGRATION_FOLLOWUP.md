# Preserve restored previews after the 05A14 series

The final integrated dev `884b0182df79ecbe1e4d32cdcdcd54084d4489c9`
passed the six-family scoped acceptance, but its full Tests run 34691796434
failed four additional cases. This corrective change addresses those failures.

- Result replacement notified the preview watcher for both the new Result and
  its template ref, requesting the same binder twice. Run that existing watcher
  after Vue's render batch so each completed replacement requests one binding.
- Mode profile changes triggered the live stroke watcher against a retained
  Result. Batch this existing watcher after profile updates and exclude mode
  transitions and suppressed restore operations. Ordinary stroke edits still
  update the mounted SVG and saved Result.
- Palette reconciliation treated split-feature outline paths as fill targets,
  changing `fill="none"` to a palette color during file-replacement Undo.
  Exclude the renderer's `__outline` paths in the existing fill predicate while
  preserving their block identity for stroke edits. Transparent filled bodies
  remain editable; connectors remain excluded.

The no-op stroke setter experiment was discarded: it did not address the
confirmed causes and adds no bytes to this patch. No assertions, test timeouts,
CI routing, structural counters, smoke budget or generated wheel bytes change.
No J08 implementation is included.

## Product and architecture preflight

Classification: `IMPLEMENT_EXISTING_AUTHORITY`. The base Web guidance requires
immutable generated artifact replacement, automatic live-edit completion and
preservation of saved Result versus active draft. The base Option Integrity
contract's OIPC-C05/C06 protects valid intent and explicit replacement; Session
compatibility keeps the last committed render distinct from per-mode drafts.
No new default, compatibility retirement or alternative Product outcome is
selected. The observed color/stroke corruption violates these existing outcomes.
There is no unresolved evidence-dependent Product choice and no non-waivable
constraint is relaxed. J08's separate settings-only behavior is unchanged.

Owners and entry paths remain `setupWatchers` -> the existing preview runtime,
`createSvgStyles` -> the existing Result persistence action, and the existing
feature-DOM fill predicate. No reactive ref, watcher, export, module, fallback
owner or canonical path is added. The superseded pre-render callback scheduling
and overly broad fill predicate are replaced in place. Architecture Gate and
Review results are retained with the exact candidate commit in the handoff proof;
maintainer review remains required regardless of machine Review classification.

## Verification

The unchanged three functional failures pass on the corrective source. The
unchanged saved-session performance case passes with one mount, one binder and
one accepted readiness receipt per Generate, and no rejected duplicate request.
New Circular and Linear controls cover repeated Generate, mode round trips and
live stroke edits with selected/mounted/export coherence. The fill-target unit
regression fails on the integrated-base helper and passes with the correction.
The final 32-case neighbor run and exact source/commit hashes are recorded in the
durable handoff. Existing failing tests are unchanged.

These local results do not mark `884b0182...` or its failed dev CI as passing.
The candidate is prepared for separate review. After authorized integration,
refresh the required checks and acceptance against actual resulting dev.
