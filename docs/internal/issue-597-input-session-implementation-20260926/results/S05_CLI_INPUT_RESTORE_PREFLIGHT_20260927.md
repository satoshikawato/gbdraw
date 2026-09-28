# S05 CLI Session input restoration — developer preflight

Status: IMPLEMENT_EXISTING_AUTHORITY, local bug fix; no new Product option or
authority receipt. Base: `c2e374cdc2454fb05b41e7dd243ebabf1fb03e39`.

## Trigger and user clarification

The unchanged real records-table CLI Session contains 12 GenBank resources and
no saved Web `config`. The CLI writer's initial `linearSeqs: []` replaced the
existing typed projection during Load. The saved SVG still displayed, but
Generate read one empty placeholder and returned `Missing GenBank file`.
This is not evidence of a researcher clearing input controls.

The user clarified: 「GenBank入力(GFF3+FASTAも含む)とか、TSVとか、描画に必要な
ファイルとかデータがすべてSession ファイルに入っているべき」 and identified the
failure as a bug. This selects correction of incomplete Session restoration;
it does not select the withdrawn A/B recovery actions or approve publication.

## Existing authority and ownership

- `docs/SESSION_COMPATIBILITY.md`, Current request ownership: resource bytes
  stored once; Web bindings point to those resources; committed render and
  separately persisted Web draft remain distinct.
- `gbdraw/session_request_codec.py::_encode_record`: existing typed writer
  embeds GenBank or GFF3 + FASTA sources. Typed drawing options/comparisons own
  their required data and resource references. Record-table selection/presentation
  data is retained in the typed request; no argv/table replay is introduced.
- PD-OI-044/045: preserve exact source identity, draft/artifact separation,
  preview-only Load and existing Generate/error/retry continuations.
- Original S05 resume instruction, acceptance 4: resolve at the existing typed
  request/resource/inventory boundary; never automatically refill an explicit
  saved empty Web draft from the committed artifact.
- Code/tests describe the fault and explicit-empty semantics; they do not
  independently authorize a different user outcome. No candidate BD is cited.

## Complete local outcome

The config import owner identifies CLI-only initial input setup only when the
admitted document has no saved `config` and has the existing CLI writer marker.
It requests initialization at the existing canonical projection owner. This is
not a trust bypass: original bindings, resources, request and catalogs are still
validated. No resource is invented, no scientific source is reconstructed from
SVG, and no new source selector/button/confirmation is added.

For that initial CLI projection, empty writer slots do not erase typed inputs.
Real original CLI input bindings, particularly original GFF + FASTA, remain
usable. Auxiliary typed file/data references are retained. A stored Web config
or an ordinary projection keeps explicit null/empty bindings authoritative,
even when the committed Result still has sources.

All requests use the existing builder, resource views, Python Worker renderer,
SVG admission and History. Load does not automatically Generate. The imported
Session document is not modified by the initial projection. On later Web Save,
the ordinary writer records the now initialized active file bindings.

## Effect review

Affected continuation: CLI Save -> Web Load -> ordinary Generate. Expected
change: existing embedded sources reach Generate, replacing the missing-input
error caused by initial slots. A saved Web draft retains its original failure
or explicit clearing. Source bytes/identity, record selectors/order, scientific
comparison/depth data, committed preview and normal editor/export/History
semantics remain required. No new compatibility default for GUI drafts,
source synthesis, argv builder, renderer, authority/map/checker, or mapped
contract body is introduced.

The old real FAIL is retained, and the original fixture is used after the fix.
Focused checks cover GB, GFF + FASTA, depth/comparison TSV, immutable imported
bindings and explicit empty Web inputs. End-to-end real Generate and final
source checks are recorded in the new S05 result; passing a focused check does
not waive any remaining S05 condition.
