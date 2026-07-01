# Feature Visibility Matching Exclusion Plan

## Goal

Change feature visibility controls so users can remove problematic CDS entries
from orthogroup and pairwise-match inputs without changing the diagram
annotation visibility policy.

The motivating case is SARS-CoV-like polyprotein annotations such as
`YP_009725295.1 ORF1a polyprotein`: when the feature would otherwise be drawn,
it should stay drawn, but it should not seed LOSATP/BLASTP pairwise matches when
it duplicates a larger overlapping polyprotein feature.

## Principles

- **SOLID:** keep drawing policy and analysis-membership policy as separate
  predicates. UI code should build rules; Python code should interpret them.
- **KISS:** keep the existing five-column feature visibility table. Do not add a
  second table or a nested policy object.
- **YAGNI:** implement only the requested scopes: this feature, current
  orthogroup members, exact product, and exact protein ID.
- **DRY:** use the existing rule-matching machinery for hash, qualifiers,
  record IDs, and feature types. Add helpers only when they prevent duplicate
  JS/Python interpretation.
- **Debt control:** every added line must pay for itself. Prefer renaming and
  tightening current behavior over adding compatibility layers, except for the
  narrow deprecated `suppress` migration needed by existing saved/manual inputs.

## Target Behavior

### User-Facing Modes

The Rich Feature Popup should expose these feature visibility modes:

| UI mode | Serialized action | Drawn in diagram | Used for protein matching |
| --- | --- | --- | --- |
| Default | no rule | Follows feature filters | Yes |
| On | `show` | Yes, bypasses feature filters | Yes |
| Off | `off` | No | No |
| Exclude from matching | `exclude_matching` | Follows feature filters unchanged | No |

`Off` keeps its name but changes meaning to "hide and exclude". The old
`suppress` action is treated as a deprecated input token that normalizes to
`exclude_matching`; new UI/session/TSV output must not serialize it.

`Exclude from matching` does not force the feature visible. It only leaves the
normal drawing decision untouched while removing the feature from protein
matching inputs.

### Drawing

`should_render_feature()` should use first-match-wins:

- `show`: return `True`.
- `off`: return `False`.
- `exclude_matching`: ignore the rule for drawing and return the normal base
  drawing decision.
- no matching rule: return the normal base drawing decision.

The normal base drawing decision remains unchanged: draw selected feature types,
plus features promoted by specific color rules.

### Analysis

`should_include_feature_in_analysis()` should return `False` when the first
matching rule action is either:

- `off`
- `exclude_matching`

It should return `True` for `show` and for features with no matching rule.

`-k/--features` remains a drawing default. It must not become a protein-analysis
filter.

## Rule Format

Keep the current TSV shape:

```text
record_id<TAB>feature_type<TAB>qualifier<TAB>value<TAB>action
```

Canonical action values:

- `show`
- `off`
- `exclude_matching`

Accepted aliases should stay minimal:

- UI mode `on` serializes to `show`.
- UI mode `off`, plus optional manual TSV aliases `hide`, `false`, and `0`,
  normalize to `off`.
- Deprecated manual TSV/session token `suppress` normalizes to
  `exclude_matching` with a Python warning. Web import/session migration may
  accept it only to migrate old data; web serialization must emit
  `exclude_matching`.
- `exclude_matching` is the only analysis-exclusion-with-unchanged-drawing
  action.

Invalid after this change:

- bare `exclude`

Rejecting bare `exclude` keeps the vocabulary narrow while preserving the
already-used `suppress` token long enough for older saved/manual inputs to load.

## Scope Application

The popup should support applying a selected mode with these scopes.

### This Feature

Use the existing exact hash rule:

```text
*    *    hash    ^<escaped svg_id>$    <action>
```

This remains the default scope and keeps current behavior simple.

### Current Orthogroup Members

Use the current orthogroup member list and create one exact hash rule per member.

Do not add an `orthogroup_id` pseudo-qualifier to the TSV. Orthogroup IDs are
derived analysis output, not source annotations, and will disappear after the
matching-excluded features are removed from analysis.

### Exact Product

Create one qualifier rule:

```text
*    <feature_type>    product    ^<escaped product>$    <action>
```

Use exact regex matching so `ORF1a polyprotein` does not accidentally match
`ORF1ab polyprotein`. Exact means anchored regex matching under the existing
case-insensitive Python rule compiler.

### Exact Protein ID

Create one qualifier rule:

```text
*    <feature_type>    protein_id    ^<escaped protein_id>$    <action>
```

This is the preferred scope for examples like `YP_009725295.1`.

Exact protein ID matching also uses anchored regex matching under the existing
case-insensitive Python rule compiler.

### Scoped Rule Priority

Feature visibility remains first-match-wins. Popup-created editor rules should
use this priority order:

1. exact hash rules for this feature and orthogroup members;
2. exact qualifier rules for product/protein ID scopes;
3. existing manual/file rules in their current table order.

This keeps the most specific editor hash rule above broader product/protein ID
rules. A user can still change the final order in the advanced Feature
Visibility table.

## Implementation Steps

1. **Normalize the action vocabulary.**
   - Update `gbdraw/features/visibility.py`.
   - Replace canonical `hide`/`suppress` output with canonical `off` and
     `exclude_matching`.
   - Keep `suppress` as a deprecated input alias for `exclude_matching` and log
     a warning when Python sees it.
   - Reject bare `exclude`.
   - Update validation messages to mention the new canonical actions and the
     deprecated `suppress` compatibility token.
   - Keep first-match-wins behavior.

2. **Split drawing and analysis semantics.**
   - Update `should_render_feature()` so `exclude_matching` preserves the normal
     drawing decision.
   - Update `should_include_feature_in_analysis()` so `off` and
     `exclude_matching` both exclude CDS entries from extraction.
   - Keep all call sites using these two predicates; do not add a parallel
     analysis-filter path.

3. **Update Web action mapping.**
   - Update `gbdraw/web/js/app/feature-visibility.js`.
   - Use UI modes `default`, `on`, `off`, and `exclude_matching`.
   - Serialize canonical actions only.
   - Migrate legacy `suppress` parser/session values to `exclude_matching`, but
     remove `suppress` from new override cache values and mode conversion.

4. **Update instant preview.**
   - Update `gbdraw/web/js/app/feature-editor/svg-actions.js`.
   - Hide SVG elements only for `off`.
   - Leave `exclude_matching` visually unchanged.
   - Keep `on` and `default` visually visible unless the generated SVG lacks
     that feature.
   - Do not try to remove pairwise/orthogroup match edges instantly; regenerated
     analysis updates those tracks.

5. **Add scoped rule construction.**
   - Add small helpers in `feature-visibility.js` for:
     - exact hash rule,
     - exact qualifier rule,
     - mode-to-action conversion,
     - escaped exact regex values.
   - Reuse these helpers from popup actions and tests.
   - Do not introduce a new table schema.

6. **Add a visibility scope dialog.**
   - Add minimal state in `gbdraw/web/js/state.js`, modeled after existing color
     and label scope dialogs.
   - Trigger it from the Feature visibility select when multiple useful scopes
     exist.
   - Offer buttons for this feature, orthogroup members, exact product, and
     exact protein ID when each scope has a usable value.
   - Apply immediately to this feature when no broader scope is available.

7. **Apply scoped rules.**
   - For this feature and orthogroup members, upsert or remove exact hash rules.
   - For exact product/protein ID, upsert or remove the exact qualifier rule.
   - Insert editor exact hash rules above editor exact qualifier rules; insert
     editor exact qualifier rules above manual/file rules.
   - After applying rules, rebuild the compatibility override cache and update
     current SVG preview for affected currently extracted features.
   - Keep the advanced Feature Visibility table as the escape hatch for manual
     cleanup.

8. **Update UI labels.**
   - Popup select:
     - `Default (follow feature filters)`
     - `On (force show)`
     - `Off (hide and exclude from matching)`
     - `Exclude from matching (leave visibility unchanged)`
   - Advanced rule table actions:
     - `Show`
     - `Off`
     - `Exclude from matching`
   - Help text should state that `Off` hides and removes the feature from
     protein comparison inputs, while `Exclude from matching` removes it from
     protein comparison inputs without changing drawing.

9. **Update session and history compatibility only for current schema.**
   - Session JSON should preserve `featureVisibilityRules` with new canonical
     actions.
   - Migrate old `suppress` override/rule values to `exclude_matching`.
   - Keep old `featureVisibilityOverrides` fallback only if it costs little, but
     map `on`, `off`, and deprecated `suppress`; do not preserve `suppress` in
     newly saved sessions.

10. **Pass visibility rules to Rich Feature Popup extraction.**
    - Update `gbdraw/web/js/app/run-analysis.js` so generated feature
      visibility TSV is passed to feature extraction and included in the
      feature-extraction cache key.
    - Update `gbdraw/web/js/services/diagram-generation.js` and
      `gbdraw/web/js/workers/diagram-generation-worker.js` so feature
      extraction requests can stage the generated TSV and pass its path.
    - Update `gbdraw/web/js/app/python-helpers.js` and
      `gbdraw/web_support/feature_metadata.py` so
      `extract_features_from_genbank_json()` compiles the same visibility rules
      and passes them to `extract_features_from_records_payload()`.
    - This keeps popup metadata, editable feature lists, and generated SVG
      feature visibility in sync. `exclude_matching` still leaves extraction
      visibility unchanged because drawing uses `should_render_feature()`.

11. **Update tests.**
    - Python tests:
      - `off` hides and excludes.
      - `exclude_matching` stays drawable but excludes from
        `extract_cds_proteins()`.
      - `suppress` warns and normalizes to `exclude_matching`.
      - bare `exclude` is rejected.
      - first-match-wins still controls both predicates.
    - Web helper tests:
      - action parsing and serialization use `show`, `off`, and
        `exclude_matching`.
      - legacy `suppress` imports/session values migrate to `exclude_matching`
        and are not serialized.
      - popup scope helpers build exact hash, product, and protein ID rules.
      - instant preview hides only `off`.
    - Packaging tests:
      - index options and helper strings no longer offer `suppress` as a new
        action.
      - `/web_feature_visibility_table.tsv` is still passed to Python protein
        extraction.
      - Rich Feature Popup extraction receives the generated feature visibility
        table.

## Files Expected To Change

- `gbdraw/features/visibility.py`
- `gbdraw/analysis/protein_colinearity.py` tests around existing extraction
  behavior, if needed
- `gbdraw/web/js/app/feature-visibility.js`
- `gbdraw/web/js/app/feature-editor/visibility-actions.js`
- `gbdraw/web/js/app/feature-editor/svg-actions.js`
- `gbdraw/web/js/app/python-helpers.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/diagram-generation.js`
- `gbdraw/web/js/state.js`
- `gbdraw/web/index.html`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/workers/diagram-generation-worker.js`
- `gbdraw/web_support/feature_metadata.py`
- `tests/test_feature_visibility.py`
- `tests/web/feature-visibility.test.mjs`
- relevant assertions in `tests/test_web_packaging.py`

Avoid touching reference SVGs unless an existing fixture intentionally exercises
the changed feature visibility actions.

## Verification Commands

Run focused checks first:

```bash
python -m pytest tests/test_feature_visibility.py -q
node tests/web/feature-visibility.test.mjs
python -m pytest tests/test_web_packaging.py -k "feature_visibility or web_run_analysis" -q
```

Then run broader checks only if the focused changes pass:

```bash
python -m pytest tests/ -m "not slow" -q
```

## Non-Goals

- Do not add a new TSV schema.
- Do not make `-k/--features` affect protein extraction.
- Do not serialize the old `suppress` token in new web/session/TSV output.
- Do not add broad bulk-edit workflows beyond the four requested scopes.
- Do not refactor orthogroup or pairwise-match internals; filtered protein
  extraction should be enough.
