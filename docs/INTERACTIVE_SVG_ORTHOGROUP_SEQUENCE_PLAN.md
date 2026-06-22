# Interactive SVG Orthogroup Sequence Access Plan

## Goal

Interactive SVG output should let users copy and download nucleotide and amino-acid FASTA for orthogroup members from Feature, Orthogroup, and pairwise/collinearity match popups, matching the GUI experience.

The sequence payload must remain single-sourced. Feature metadata already owns the sequence data, so orthogroup and match popups should resolve sequences through member `featureSvgId` references instead of duplicating sequence strings in orthogroup or match metadata.

## Current State

- GUI orthogroup drawer enriches orthogroup members from extracted feature sequence metadata and then builds FASTA on demand.
- Interactive SVG export already embeds feature metadata in `<metadata id="gbdraw-interactive-feature-metadata">`.
- Rich interactive feature payloads include `nucleotide_sequence`, `amino_acid_sequence`, `nucleotide_fasta`, and `amino_acid_fasta`.
- Interactive SVG match payloads are precomputed during export, but their member rows currently carry table details only, not reusable sequence actions.
- Interactive SVG feature popups already show per-feature sequence blocks with Copy support.

## Design Principles

### SOLID

- Single Responsibility: keep sequence resolution, FASTA assembly, filename generation, clipboard, and download behavior as separate small helpers.
- Open/Closed: extend popup rendering through reusable action helpers rather than special-casing every popup type.
- Interface Segregation: popup renderers should consume a narrow member-row shape containing identifiers and display fields, not full feature objects.
- Dependency Inversion: orthogroup and match renderers depend on `featuresById` lookup and stable member identifiers, not on where the feature metadata came from.

### KISS

- Do not introduce new dependencies, build steps, or external files for standalone SVG behavior.
- Use the existing embedded metadata and embedded script model.
- Prefer existing FASTA strings on feature payloads. Generate FASTA from raw sequence only as a fallback.

### DRY

- One helper should build group/member sequence FASTA for both feature orthogroup sections and match popups.
- One helper should render Copy/DL action controls for nucleotide and amino-acid FASTA.
- One download helper should handle all text file downloads.

## Data Contract

Feature payload is the only sequence source:

- `feature.svg_id`
- `feature.nucleotide_fasta`
- `feature.amino_acid_fasta`
- `feature.nucleotide_sequence`
- `feature.amino_acid_sequence`

Orthogroup and match payloads may carry identifiers and display fields:

- `member.featureSvgId` / `member.feature_svg_id`
- member display fields such as record, coordinates, protein ID, role, confidence, product/note

Orthogroup and match payloads must not carry duplicated `ntFasta`, `aaFasta`, `memberNtFasta`, or `memberAaFasta` strings in standalone interactive SVG metadata.

## Implementation Plan

### 1. Add Standalone Runtime Helpers

Target: `gbdraw/web/js/services/standalone-interactivity.js`

Add small pure helpers inside the embedded standalone runtime:

- `memberFeatureSvgId(memberOrRow)`
- `featureForMember(memberOrRow)`
- `featureFasta(feature, sequenceKind)`
- `memberFasta(memberOrRow, sequenceKind)`
- `memberSequenceFilename(memberOrRow, sequenceKind, orthogroupId)`
- `groupFasta(memberRows, sequenceKind)`
- `groupSequenceFilename(orthogroupId, displayName, sequenceKind)`
- `downloadText(filename, text, mimeType)`

Rules:

- `featureFasta()` first returns existing `feature.nucleotide_fasta` or `feature.amino_acid_fasta`.
- If FASTA is absent, it may fall back to raw `feature.nucleotide_sequence` or `feature.amino_acid_sequence` using a minimal FASTA formatter.
- Empty sequence data should produce no button rather than an empty download.

### 2. Preserve Member Identifiers in Runtime Rows

Target functions:

- `orthogroupMemberTableRows()`
- `normalizeMatchMemberRows()`
- `buildStandaloneMatchMemberRows()`
- `normalizeMatchBlockOrthogroups()`

Extend member-row objects with:

- `featureSvgId`
- `orthogroupId` when available
- `displayName` or parent orthogroup display label where useful for filenames

Do not attach sequence strings to these rows.

### 3. Add Reusable Sequence Action Renderer

Add a renderer that takes member rows plus optional orthogroup metadata and returns action buttons:

- Copy nt
- DL nt
- Copy aa
- DL aa

The renderer should include counts when useful, for example `Copy nt (4)`, where the count is the number of members with available FASTA.

Use data attributes that point to internal action indexes, similar to existing `data-copy-index`, so large FASTA strings are kept in runtime arrays rather than duplicated into HTML attributes.

Suggested runtime arrays:

- `copyValues`
- `downloadValues`, each item `{ filename, text, type }`

### 4. Extend Feature Popup Orthogroup Members Section

Target:

- `renderOrthogroupMembers()`

Add sequence actions to the Orthogroup members block title or block action row.

Behavior:

- The table remains focused on member metadata.
- The new actions aggregate sequences for all orthogroup members that resolve to a feature with sequence data.
- Individual member sequence actions are optional for the first pass; group-level actions are the required minimum.

### 5. Extend Match and Collinearity Popups

Targets:

- `renderMatchMemberTable()`
- selected block orthogroup rendering in `renderMatchSections()`

Add the same sequence actions to member tables inside:

- Orthogroup match popup summary
- Pairwise match popup orthogroup section
- Collinearity block selected orthogroup section

Behavior:

- For a direct orthogroup match, actions operate on the orthogroup member table.
- For a collinearity block, actions operate on the selected orthogroup once the user selects a row from the covered orthogroups table.
- If no selected orthogroup exists, no group sequence action is shown for that selected section.

### 6. Add Download Click Handling

Target:

- popup click handler near existing `data-copy-index` handling

Add handling for `data-download-index`:

- Resolve indexed download payload.
- Create a Blob.
- Trigger an `<a download>`.
- Revoke the object URL.

Clipboard fallback should remain unchanged.

### 7. Keep Export Metadata Compact

Targets:

- `buildStandaloneOrthogroupPayloads()`
- `buildStandaloneMatchPayloads()`
- `buildStandaloneMatchMemberRows()`

Audit these functions after changes to ensure they do not serialize aggregated FASTA strings into `metadata.textContent`.

Acceptable:

- identifiers
- table display fields
- filenames if small and useful

Avoid:

- per-member FASTA strings outside feature payloads
- group-level aggregated FASTA strings

## UI Details

Use concise action labels in popup blocks:

- `Copy nt`
- `DL nt`
- `Copy aa`
- `DL aa`

Disable or omit unavailable actions rather than showing buttons that do nothing. Omission is simpler and avoids explaining unavailable states in the SVG UI.

Suggested filenames:

- Orthogroup group nt: `<orthogroup_id>_<display_name>_nt.fna`
- Orthogroup group aa: `<orthogroup_id>_<display_name>_aa.faa`
- Single member nt: `<orthogroup_id>_<member_id>_nt.fna`
- Single member aa: `<orthogroup_id>_<member_id>_aa.faa`

Use the existing safe filename normalization style.

## Testing Plan

### Unit-Style JS Tests

Extend `tests/test_pairwise_match_popup.py` or add a focused Node-based test for standalone helper behavior if the current test harness can import or evaluate the relevant module section.

Coverage:

- Orthogroup member rows retain `featureSvgId`.
- Group FASTA is assembled from `featuresById` feature payloads.
- Missing member sequences are skipped.
- No group-level FASTA is serialized into match payload metadata.

### Packaging Tests

Extend `tests/test_web_packaging.py`:

- Interactive SVG metadata still contains feature sequence data.
- Orthogroup metadata does not duplicate sequence data.
- Match metadata does not duplicate aggregated sequence data.
- Standalone script contains download handling for sequence action buttons.

### Manual Verification

Generate an interactive SVG from a linear orthogroup comparison and verify:

- Feature popup Sequence tab can copy/download nt and aa FASTA.
- Feature popup Orthogroup members section can copy/download group nt and aa FASTA.
- Orthogroup match popup can copy/download member group nt and aa FASTA.
- Collinearity block popup can select a covered orthogroup and copy/download that selected group's nt and aa FASTA.
- Exported SVG file size does not grow by a second copy of the sequence payload.

## Risks and Mitigations

- Risk: some orthogroup members may not resolve to a feature payload.
  Mitigation: skip unresolved members and show actions only when at least one sequence is available.

- Risk: GFF or partial-region inputs may lack amino-acid sequences.
  Mitigation: nucleotide and amino-acid actions are independently available.

- Risk: browser clipboard permissions may fail in standalone SVG contexts.
  Mitigation: keep the existing prompt fallback.

- Risk: popup code grows further inside a large standalone script.
  Mitigation: add small helpers with clear names and reuse them from all popup renderers.

## Acceptance Criteria

- No sequence strings are duplicated into orthogroup or match metadata.
- Orthogroup-related popup sequence actions resolve data from feature payloads by `featureSvgId`.
- Copy and download work for nucleotide and amino-acid FASTA where data exists.
- Existing feature popup sequence behavior remains intact.
- Existing interactive SVG search, hover, panning, and popup interactions continue to work.
