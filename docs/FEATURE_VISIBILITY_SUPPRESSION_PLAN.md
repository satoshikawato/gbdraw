# Feature Visibility And Suppression Plan

## Goal

Support two related behaviors without turning feature filtering into a second
pipeline:

- A feature matched by the feature-specific color table (`-t/--table`) is drawn
  even when its feature type is not listed in `-k/--features`.
- A user can suppress an individual feature through
  `--feature_visibility_table`, including a feature whose type is listed in
  `-k`; suppressed features are excluded from LOSAT/BLASTp pairwise links,
  orthogroup inference, and collinearity calculation.

Keep the implementation small. `-k` remains a drawing default, not an analysis
filter.

## Principles

- **KISS:** keep one draw predicate and one analysis predicate. Do not add a new
  filtering subsystem.
- **DRY:** reuse the existing feature rule matching helpers for hash, location,
  record location, and qualifiers.
- **SOLID:** separate responsibilities:
  - input loading chooses candidate feature types,
  - drawing decides visible features,
  - analysis decides suppressible features.
- **YAGNI:** do not add broad per-feature policy objects, new table formats, or
  global analysis filters until a concrete use case needs them.
- **Debt control:** every added helper must remove duplicated logic or prevent a
  real bug. Avoid cosmetic refactors while implementing this change.

## Terminology

- `-k/--features`: default feature types to draw.
- `-t/--table`: feature-specific color table. Matching rows can force drawing of
  individual features outside `-k`, but this table does not control
  suppression.
- `--feature_visibility_table`: per-feature visibility table. This controls
  explicit `show`, `hide`, and `suppress` overrides. `--feature_table` is the
  currently implemented legacy name and should be replaced in user-facing docs,
  generated Web arguments, and examples.
- `show`: draw the feature even if `-k` would omit it.
- `hide`: do not draw the feature, but keep it available for analysis.
- `suppress`: do not draw the feature and do not include it in protein
  comparison, orthogroup, or collinearity inputs. Only
  `--feature_visibility_table` can request suppression.

## Behavior

### Drawing

Draw a feature when the first matching rule says `show`, or when no explicit
visibility rule matches and either:

- the feature type is in `-k`, or
- the feature matches a `-t/--table` color rule.

Do not draw a feature when the first matching visibility rule says `hide` or
`suppress`.

### Analysis

Include CDS features in protein analysis unless a matching visibility rule says
`suppress`.

Do not use `-k` to decide protein analysis membership. For example, removing
`CDS` from `-k` hides CDS by default, but pairwise/orthogroup/collinear analysis
still uses CDS unless the feature is explicitly suppressed.

### GFF3 Loading

When parsing GFF3/FASTA, load a candidate superset instead of only `-k`:

```text
candidate_types =
  -k feature types
  + concrete feature_type values from -t/--table
  + concrete feature_type values from --feature_visibility_table
```

If either table uses `feature_type="*"`, load all feature types. Final drawing
and suppression are still decided after parsing.

This keeps GFF loading light while avoiding false negatives for table-driven
features.

## Implementation Steps

1. Rename the CLI/API surface from `--feature_table` to
   `--feature_visibility_table`.
   - Keep `-k/--features` unchanged as the type-level drawing default.
   - Prefer `feature_visibility_table` for variable names, docs, examples,
     session metadata, and Web-generated arguments.
   - Keep `--feature_table` as a deprecated or hidden alias because it is already
     implemented in CLI/API/Web/session paths. It should write to the same
     destination and should not be documented as the primary option.

2. Extend visibility action parsing in `gbdraw/features/visibility.py`.
   - Normalize `show/on/display/include/true/1` to `show`.
   - Normalize `hide/off/false/0` to `hide`.
   - Normalize `suppress/exclude` to `suppress`.
   - Keep first-match-wins behavior.
   - Reject or ignore suppression actions from any source other than
     `--feature_visibility_table`; `-t/--table` remains color-only plus drawing
     promotion.

3. Add a small analysis predicate in `gbdraw/features/visibility.py`.
   - Suggested name: `should_include_feature_in_analysis(feature, rules, record_id=None)`.
   - Return `False` only when the first matching visibility rule action is
     `suppress`.
   - Export it from `gbdraw/features/__init__.py`.

4. Add specific-color matching support without duplicating regex logic.
   - Reuse or lightly extract the existing color-rule matching logic from
     `gbdraw/features/colors.py`.
   - Add an optional `specific_color_rules` or `color_map` argument to
     `should_render_feature()`, or wrap it with a drawing-only helper.
   - The helper should answer only: "does this feature match any specific color
     rule?"

5. Update drawing call sites to pass the color-rule match context.
   - `create_feature_dict()`
   - `precompute_used_color_rules()`
   - `check_feature_presence()`
   - interactive metadata extraction, so clickable feature payloads match the
     rendered feature set.

6. Add a candidate feature-type resolver for GFF3 loading.
   - Suggested location: near CLI/API input setup, or a small helper in
     `gbdraw/features/visibility.py` if both CLI modes need it.
   - Use it in `circular.py`, `linear.py`, and API paths that load GFF3/FASTA.
   - Pass candidate types to `load_gff_fasta()` while preserving the original
     `selected_features_set` for drawing.

7. Apply suppression to protein extraction.
   - Add optional `feature_visibility_rules` to `extract_cds_proteins()` and
     `extract_web_stable_cds_proteins()`.
   - Skip CDS features when `should_include_feature_in_analysis()` returns
     `False`.
   - Pass the same filtered extraction into pairwise, orthogroup, and collinear
     code. Collinearity units will then be built from the filtered extraction.

8. Update Web wiring only where needed.
   - Add `suppress` as a feature visibility override mode.
   - Serialize suppress rows to `/web_feature_visibility_table.tsv`.
   - Include suppress state in protein extraction cache keys.
   - Pass the generated feature visibility table to
     `extract_cds_protein_fasta()` so browser LOSAT cache generation and Python
     diagram generation use the same CDS set.

## Tests

Add focused tests rather than broad snapshot churn:

- `-t` row for a feature type outside `-k` draws only matching features.
- GFF3 candidate loading includes `-k` plus concrete
  `-t`/`--feature_visibility_table` feature types.
- `feature_type="*"` in either table makes GFF3 loading keep all features.
- `hide` removes a feature from drawing but leaves it in `extract_cds_proteins()`.
- `suppress` removes a feature from drawing and from `extract_cds_proteins()`.
- `suppress` is accepted only from `--feature_visibility_table`; `-t/--table`
  remains a color table and has no suppression action.
- pairwise/orthogroup/collinear paths receive filtered protein extraction.
- Web-generated `/web_feature_visibility_table.tsv` preserves `show`, `hide`, and
  `suppress` distinctly.

Reference SVG updates should only be needed if existing fixtures intentionally
exercise the changed visibility semantics.

## Compatibility

Existing `show` and `hide` behavior remains drawing-compatible. The only
intentional semantic change is that `suppress` and `exclude` now mean analysis
exclusion as well as hidden drawing. This matches the new requirement and should
be documented in CLI help and web UI labels for `--feature_visibility_table`.

Because `--feature_table` is already implemented in CLI/API/Web/session paths,
keep it as a compatibility alias for one release cycle or as a hidden alias.
New docs, sessions, and Web-generated invocations should use
`--feature_visibility_table`.

## Non-Goals

- Do not make `-k` an analysis filter.
- Do not remove `-k/--features` or require a TSV file for ordinary type-level
  feature selection.
- Do not add a new table format.
- Do not add `suppress` semantics to `-t/--table`.
- Do not implement non-CDS analysis suppression until another analysis path
  needs it.
- Do not refactor LOSAT, orthogroup, or collinearity internals beyond passing a
  filtered `ProteinExtractionResult`.
