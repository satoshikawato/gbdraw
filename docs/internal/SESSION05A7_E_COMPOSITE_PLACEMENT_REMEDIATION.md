# SESSION 05A7-E: composite placement remediation

PR B starts at exact `origin/dev` `bed447b813d07581764b92373da62bcf464e8517`,
the squash merge of PR #508 after final-head human review. The user's original
worktree remains untouched. This PR addresses only 05A4-08.

## Original failure and independence

Both defects reproduced before any production changes on
`ac1bb944cd26441cf7552823769a54ced7c53881`. The retained two-source
mitochondrial/lambda composite and full 11-record Vibrio Circular construction
lost every placement option after Generate with unchanged source/configuration.
PR A repaired the multipart label and left this placement failure intact.
The PR B base also fails the new composite transition contract while its single
record browser control passes. The roots are independent.

Generate changes the committed source representation: separate component
resources become one concatenated GenBank resource. The UI retains the same
composite file view. `record-display-options.js::matchesSavedSource` compared
each component payload with the concatenated payload. That could never prove
current-source identity, so `isCurrentFeature` rejected every placement including
Auto. The catalog, record identities, feature type/strand, track preset, and
Separate Strands state were not lost or reinterpreted. Save retains the ordered
component bindings, so the same comparison also failed after a fresh Load.

## Existing owner correction

`session-resource-backing.js` already owns immutable ordinary and composite file
views, ordered component bytes, and their LF concatenation rule. Its new
`matchesSessionResourceDescriptor` operation verifies the committed descriptor
against either an original component or the complete ordered composite bytes.
`record-display-options.js` delegates its former payload comparison to that owner.
The exact source object and paired-source guards stay in the record display
owner; a real new File/view cannot inherit an old feature's capability.

The comparison checks declared sizes, all bytes, order, repeated/empty components,
CRLF and appended LF boundaries. On a successful whole-composite comparison,
the existing backing's descriptor slot references that immutable serialized
representation. Repeated checks reuse the verified descriptor and retain no
decoded content. This adds no Result/editor cache, second catalog, placement
owner, identity namespace, or persisted metadata. Component file bindings and
source metadata remain unchanged. The old editor-local payload comparison is
removed. Native File behavior and ordinary resource matching are retained.

`placement-actions.js` remains the capability/action owner. Capability follows
the draft feature slot, using the same track-slot resolver as request projection.
`valueFor` reads the current placement override independently. The popup renders
`choices`, and `setPlacement` rejects unavailable choices through that same
predicate. Middle permits Auto/Main/Outward/Inward; tuckin and spreadout permit
Auto/Main. An existing outward value does not make an unsupported lane available.

An initial two-line candidate used the existing payload-owner association and
fixed Generate, but its privileged import was disallowed and its ephemeral
association did not survive fresh Load. That candidate was removed. No allowlist,
CI policy, compatibility field, or schema was changed to admit it.

## Architecture and Product evidence

Ordinary non-increasing correction: the source-backing owner gains one operation;
the editor delegates its source payload comparison there. There is one canonical
source comparison path and one placement capability owner. No parallel owner,
compatibility reader, production module, reactive state, or watcher is retained
or introduced. Privileged permissions and import fan-out do not increase; no
first-party cycles are added. The new public export requires ordinary human
review, not an architecture exception. Rollback is a revert of this PR.

Product classification: `IMPLEMENT_EXISTING_AUTHORITY`. The base Web live-edit
and reactive availability contracts in `gbdraw/web/CLAUDE.md` require accepted
editing to remain generatable, enabled rendering and action admission to share
a predicate, and editing to remain bound to current source data. The fix preserves
that affordance across Generate, Save/Load, and mode changes; it selects no new
placement semantics or persisted compatibility outcome. The original failing
reproductions and the existing slot owner establish the implementation defect.
No unresolved Product choice requires EVIDENCE_REQUIRED or PRODUCT_DECISION_REQUIRED,
and the final correction needs no forbidden architecture or persistence exception.

Session 41, request 7, feature catalog 3, and file-binding schema 2 remain unchanged.
PR smoke remains ten cases. The new browser test is discovered only by the existing
full-functional selection; CI policy and generated Gallery artifacts are unchanged.

## Regression coverage

- `record-display-options.test.mjs`: exact component-to-combined transition through
  the real file view and placement actions, same-name/same-byte replacement rejection,
  and unsupported draft lanes independently of the selected value.
- `session-resource-backing.test.mjs`: exact composition, order, empty/repeated parts,
  LF/CRLF boundaries, mismatched same-length bytes, replacement, and no repeated decoding
  after successful verification or decoded-content release.
- `composite-placement-regeneration.playwright.spec.js`: P1 single, P2 composite,
  P3 repeat Generate, P4 label plus placement, P5 Save/fresh Load/Generate, P6 mode
  round trip, P7 actual same-byte File replacement, P8 preset/Separate Strands matrix.
  It compares semantic enabled sets, current value, DOM options, invalid-action
  rejection and renderer placement targets, preserving both records and their order.

The real CLI writer constructs the two-source test session from the existing
Gallery source sessions. Cold browser contexts block external requests and require
no page errors. The retained original full 11-record Vibrio construction remains a
separate acceptance case; it is not reduced to the two-source regression fixture.

Evidence lives under `/tmp/gbdraw-session05a7-e-evidence/`. Initial policy and P5
failures are retained as diagnosis and are excluded from passing evidence. The
large-source allocation measurement motivated comparing decoded binary strings
without unnecessary Uint8Array construction. Chromium verifies all 71,480,904 bytes
in 154–170 ms across three trials; 1,000 subsequent checks take 0.2–0.4 ms. These are
local source-comparison measurements, not a new performance baseline or an output
cache. The existing performance and exact-dev staging gates remain authoritative.

Deterministic verification: 37 focused source/placement contracts, 398 fast Web
contracts, and 137 architecture contracts pass. The local final Web checker reports
Gate PASS / Review REQUIRED solely for the new source-backing public operation.
