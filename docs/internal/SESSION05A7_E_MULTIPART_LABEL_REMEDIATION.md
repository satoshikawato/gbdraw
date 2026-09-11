# SESSION 05A7-E: multipart label remediation

Starting `origin/dev`: `ac1bb944cd26441cf7552823769a54ced7c53881`.
The remote was fetched and matched before creating the isolated worktree.
The user's existing worktree and uncommitted changes were preserved.

## Original reproductions

Both findings reproduce on the exact starting revision, before production edits.

- **05A4-01:** load the retained tobacco Gallery session, Generate, edit the
  multipart `rps12` label, Generate. The second Generate returns an error:
  `Sanitized SVG content is missing or ambiguously binds an editable Label.`
  Label-only is sufficient. The previous Result remains usable and exportable.
- **05A4-08:** load the retained CLI-resumed mitochondrial/lambda composite,
  inspect placement, Generate, inspect the same feature. All four options change
  from enabled to disabled. The original 11-record Vibrio construction also
  reproduces the failure, retaining every record and its order.

The tobacco target is record `NC_001879.2`, CDS `NitaCp049`, protein
`NP_054549.2`. Its biological and rendered IDs are `f6a2826aa`; its editor ID is
`record-1\u0000f6a2826aa`. Its three biological parts, in source order, are
`72224..72337`, `100624..100855`, and `100062..100087`, all on the negative
strand. Three SVG block elements belong to that logical feature.

## Boundary diagnosis

The roots are independent. The label repair leaves the composite reproducer
failing with every placement option disabled after a successful Generate.

**05A4-01 classification: wrong SVG occurrence identity.** The catalog retains
the correct biological identity and all three parts. Native Circular rendering
emits one label for this biological feature, anchored to its selected segment.
It did not emit the label's feature binding. The editor grouped every feature's
geometry into a centroid and assigned labels by distance with a one-use rule.
The two native `rps12` labels became bound to `f8ec4b35c` and `f5ba16ffb`,
leaving `f6a2826aa` unbound. A logical label override could be accepted, but
the subsequent mounted-label validation found zero bindings. It did not find
multiple legitimate labels or unstable renderer IDs.

The original native SVG still contained `rps12` after the edit's request reached
the renderer. A second fault in the same required regeneration chain reused
compiled label rules after an explicit label table replaced their input. The
canonical TSV correctly selected `NitaCp049`; the renderer used the old rules.
Supplying native identity alone made Generate succeed but exposed a Result /
mounted-text disagreement. The completed repair also invalidates the compiled
label-override rules when their table is replaced. It does not invalidate
unrelated state or introduce a cache.

**05A4-08 classification: composite source identity projection.** Generate
replaces the committed component resource descriptors with one combined
resource. The unchanged UI source remains a composite file view.
`record-display-options.js::matchesSavedSource` compares each component's
payload with the combined payload and rejects the current source.
`placement-actions.js` consequently disables even Auto. Neither the feature
catalog nor the track configuration loses its placement semantics.

Generate therefore changes the representation used to prove current-source
identity, while the placement owner still derives valid targets from the draft
track configuration. It does not derive capability from the selected value.
The label case collapses multipart geometry to a centroid for inference; it
does not need an artificial single-part biological identity or a relaxation of
the missing/ambiguous-label gate.

## Owners and change

- `gbdraw/labels/circular.py`: project the existing feature hash and existing
  duplicate-source instance ID into the prepared label.
- `gbdraw/render/drawers/circular/labels.py`: emit the already-supported
  `data-label-feature-id` on the label text, for embedded, horizontal, and
  radial labels.
- `gbdraw/api/diagram.py`: keep the label reference aligned when existing
  copied-record identity binding namespaces the corresponding feature;
  invalidate compiled label rules when an explicit replacement table is attached.

No placement production code is changed in this PR. The catalog, biological
ID namespace, editor override map, label reconciliation gate, Generate admission,
Result owner, and session writer remain the existing owners. The renderer
continues to decide label count and geometry. No label is copied to arbitrary
parts, and no feature name or fixture is special-cased in production.

Architecture: ordinary non-increasing correction at existing projection owners.
No new catalog, owner, production module, public export, compatibility reader,
reactive state, or Result/editor cache. Current rendered Circular labels carry
the identity the editor already accepts; existing handling of saved SVG and
Linear labels is retained. Rollback is a revert of this PR.

Product: `IMPLEMENT_EXISTING_AUTHORITY`. The base Web live-edit contract requires
accepted editor intent to survive Generate and keep Result, preview, and export
consistent. `docs/SVG_SEMANTIC_HOOKS.md` explicitly requires semantic identity
for compound features and forbids inferring it from geometry. The existing
regeneration contract preserves post-generation editing. No alternative product
outcome or persisted compatibility disposition is selected.

Schema disposition: unchanged session version 41, request schema 7, feature
catalog schema 3, and file-binding schema 2. The SVG attribute already exists in
the current editor and sanitization profile. No required persisted field or
version is added. No reference SVG or Gallery artifact is regenerated.

## Verification

`tests/test_circular_label_identity.py` adds deterministic coverage for all three
label rendering policies, horizontal/radial labels, multipart membership,
colliding source hashes, copied records using existing record-instance binding,
and repeated rendering with replacement label tables. The latter test failed
with old label text before the rule invalidation was added.

`tests/web/multipart-label-regeneration.playwright.spec.js` covers M1–M7 with the
retained tobacco session: single-part label, multipart label, placement-only,
placement plus label, fresh Save/Load, mode round trip, and Undo/Redo. It checks
the actual outgoing request against the selected Result's committed request,
logical part membership, durable override maps, semantic Result/mounted/export
agreement, and the unchanged second `rps12` feature. The new browser cases enter
full functional discovery; the ten-case PR smoke budget is unchanged.

Final local results:

- Core PR selection: 3,738 passed, 17 skipped (Python 3.13.3, 182 seconds).
- Focused Python labels/request suite: 142 passed, including the 13 new cases.
- Existing read-only SVG comparisons: 16 passed; tracked references unchanged.
- Final browser matrix and recent-fix controls: 31 passed (Chromium, 4.4 minutes).
- Circular warning and Auto controls: 5 passed; the final 31-case run also
  includes both nullable-region History tests.
- Four focused JavaScript contract files passed. Missing/ambiguous label
  rejection and invalid placement controls remain unchanged.
- Ruff and whitespace checks passed. Local Web policy: Gate PASS, Review CLEAR,
  with no registered architecture or Product delta.

Browser contexts block all external requests. Tests use the packaged local
runtime assets and freshly generated browser wheel, cold session loads and
repeat generation, 1600×1000 desktop and 390×844 mobile. Final desktop and mobile
screenshots and actual exported SVGs are retained for visual review. Hosted
integration gates remain separate from these local results.

Evidence is retained under `/tmp/gbdraw-session05a7-e-evidence/`: original
reproducers `M01`, `M18`, `J05`; raw Worker output and complete target catalogs in
`D01-R2`; intermediate renderer/request diagnosis in `D02`; repaired original
cases `M01-A4`, `M02-A4`; and independent placement failure in `M18-A4`.
Original and candidate SVGs, source order, component sizes, requests, catalogs,
override maps, downloads, screenshots, and runtime logs are retained there.

Intermediate attempts remain recorded. `M01-A`/`M02-A` used the old wheel after
a build-dependency fetch failed; `M18-A` overlapped wheel replacement and saw a
404. They are excluded from candidate evidence. `A2` exposed svgwrite's strict
SVG-1.1 attribute validation; `A3` exposed stale compiled label rules. The first
combined control run lost a Playwright-owned shared server and was stopped;
the final run uses an independently served loopback origin. None of these
attempts is counted as a passing regression.
