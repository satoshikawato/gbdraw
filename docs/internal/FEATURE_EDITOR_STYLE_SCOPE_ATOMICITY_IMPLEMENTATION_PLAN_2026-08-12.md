# Feature editor style scope and atomic legend consistency implementation plan

- Date: 2026-08-12
- Status: implementation paused at the 2026-08-13 session handoff; Phases 1–5 are largely implemented, but completion is blocked by a real Similarity-group projection regression and the final broad verification gates
- Planning baseline: `bugfix-20260812` / `c71d998b5362e556ce003ab5c615127dc7fca577`
- Audit revision: 2026-08-12; Circular batch scope, exact record-instance identity, History/watchers, Generate/save concurrency, v40 recovery, and phase slicing revised before execution
- Regression origin: `7aab08f8f919cb941a015143e836acccb20a3dc2`
- Primary surfaces: Web Feature editor, legend editor, PreviewRuntime, History, canonical request/session projection

Related contracts and plans:

- [Web application maintenance contract](../../gbdraw/web/CLAUDE.md)
- [Session and request compatibility history](../SESSION_COMPATIBILITY.md)
- [Performance regression remediation execution plan](PERFORMANCE_REGRESSION_REMEDIATION_EXECUTION_PLAN_2026-08-11.md)
- [Responsive Preview Phase 1 implementation plan](RESPONSIVE_PREVIEW_PHASE_1_IMPLEMENTATION_PLAN_2026-08-11.md)
- [gbdraw v0.14.0 roadmap](gbdraw_v0.14.0_codex_roadmap.md)

This document is an implementation plan. At creation time it changes documentation only. Production code, tests, Gallery assets, browser wheels, and reference SVGs have not been changed by this plan.

Execution correction (2026-08-12): `HEAD` is `68bbc8d047314670cffc9f66177a9cb12e492ab6`; its tracked tree matches the planning baseline. The named Gallery session has SHA-256 `5e566c305535bb4caccf50723158462471ecf82ce5f71ecc2fda57f91e80f937`. Its canonical specific-color table already matches the four stored `config.rules`; the stale resource is the canonical default-color table, which differs from `config.colors` and `ui.appliedPaletteColors`. Recovery therefore validates both color resources and repairs only a mismatch explained by the saved catalogue, overrides, and Result. The named file cannot serve as a positive specific-table mismatch fixture.

The earlier overlapping Web changes are committed in the planning baseline. The plan file itself is untracked at this revision. Phase 0 must still inspect the current worktree and re-resolve every owner before implementation; symbol names and line locations below describe the planning baseline, not a license to overwrite later work.

## 0. 日本語要約

今回の回帰は、`Auto / None / Color` の汎用部品を Feature editor に導入した際、既存の scope-aware な入口ではなく個別変更用 setter へ UI を接続したことが直接原因である。修正は旧 handler へ戻すだけで終わらせず、迂回経路そのものを作れない構造へ変更する。

- UI が呼べる fill mutation は `requestFeatureFillChange()` だけにする。raw setter、rule 書き換え、Feature/凡例 DOM 書き換えは module-private にする。
- fill の変更要求から pure planner が matching rule、同じ feature type、同じ凡例、rendered label、source annotation label、Similarity group、selected features、single の候補と対象を計算する。group 候補がある場合、選択確定前は何も変更しない。
- semantic scope と Result extent を別の値として扱う。Circular batch の `matching-rule` など広域の意味操作は、該当する全 Result を一つの transaction で更新する。active Result だけへの広域変更は暗黙の縮退ではなく、将来追加できる明示的な `current-result` snapshot exception とする。本計画で追加する record-instance selector は `single` と明示 selection の再生成に使い、広域 snapshot exception は rule 爆発と別の product semantics を伴うため本計画では提供しない。
- `Auto` は Feature editor から外し、実効色と継承元を表示しながら picker を常に操作可能にする。`Use inherited` と `No fill` は同じ planner/command を通る明示操作にする。
- `manualSpecificRules` と applied palette/default-color state を renderer-visible fill の意味的 owner とする。`featureColorOverrides`、Feature由来の凡例、popup state、live SVG、保存済み Result は derived projection/cache とする。Specific Rules editor、preset、TSV import、palette apply が live semantics を変える場合も、同じ all-affected-Result command を通し、単なる revision increment で済ませない。
- scope 確定後、全 async prerequisite と Result ごとの legend layout を active Result の外で preflight する。rules、derived cache、active DOM、全 affected `results[].content`、凡例 model、popup、semantic revision を一つの同期 commit で反映し、一つでも準備または commit に失敗すれば全状態を元へ戻す。
- editor の凡例追加・rename は main-thread Pyodide を起動せず、同じ font/style を使う offscreen SVG で測定する。既存 swatch の色変更には測定自体を行わない。
- scope dialog と cancel は History に入れず、確定した編集だけを一つの undoable command とする。History の intent capture と stack transition まで transaction 境界に含め、post-apply capture が失敗した場合は補償操作を行う。Undo/Redo も同じ atomic projection と補償契約を使う。
- `manualSpecificRules` の deep watcher は command 後の第二の mutation owner にしない。command 中は抑止または撤去し、同じ明示的 reconciler を一度だけ呼ぶ。History の補償 rollback は通常の Undo と区別し、revision を進めず exact before-state へ戻す。
- cached canonical request が編集後の rules または applied palette/default colors より古い状態で session 保存されないよう、両 color resource を copy-on-write で overlay する semantic revision/fingerprint を導入する。Generate と session export は開始時の immutable style snapshot と revision を保持し、adopt/download 直前の不一致を reject、再起動、または coverage を証明した rebase にする。
- first-parent `main` 上の v40 Gallery session は default-color resource の stale overlay を示しているため、v40 recovery は必須とする。specific-color 側は保存済み `config.rules` と一致している。両 color resource を検査し、catalogue、derived override、saved Result で差異を説明できる場合だけ canonical resource を修復する。
- reproducible single/selection editing のため、描画前に共有 `FeatureInstanceIdentityPlan` を構築し、reserved exact TSV selector `__gbdraw_instance_hash__` を導入する。これは public persisted vocabulary であり、session writer を v41、canonical request を schema 6 に上げる。既存 schema-3 catalogue は additive に読み、v40 から安全に identity を補完できない場合は public hash へ縮退せず exact scope を無効化して再Generateを求める。
- fill を最初の vertical slice として完成させる。label/stroke が二つ目の実利用者になった段階で、共通の scope/command protocol だけを抽出する。先に汎用 style framework は作らない。

テストは、architecture contract、pure planner、atomic command、legend projection、failure injection、actual popup/list DOM、Circular single/grid/batch と Linear の save-load-regenerate、session compatibility、PreviewRuntime/performance の各層で実装する。テストを設計ルールの代用にはせず、raw setter 非公開、scope-before-mutation、all-affected-Result atomic command、derived ownership という構造を検証する防波堤として扱う。

## 1. Outcome

Restore the Feature editor workflow in which a user can edit one feature or explicitly choose a semantic group such as the matching specific-color rule, the same Feature type, the same legend item, the same rendered label, the same source annotation label, or a Similarity group. Remove the `Auto`-as-lock interaction from Feature fill and stroke controls. Make every accepted edit update the feature, rule intent, legend, mounted SVG, saved Result, History, and session regeneration basis as one observable operation.

The completed implementation must provide all of the following:

1. An inherited Feature color is displayed as an editable effective color, not as a disabled `Auto` picker.
2. Changing fill from the popup or Feature list always enters the same scope-aware request path.
3. If more than one meaningful scope exists, no mutation occurs until the user chooses one.
4. A group edit updates every target and its existing legend entry in the same commit.
5. A single-feature edit preserves the old group color and creates or reuses a distinct legend caption in the same commit.
6. `None` and `Use inherited` are explicit semantic actions routed through the same planner and command as a color value.
7. A failed or stale edit leaves rules, overrides, Feature DOM, legend DOM/model, Result content, and History unchanged.
8. Save, load, Undo/Redo, fresh regeneration, and candidate replay preserve the accepted edit in Circular and Linear diagrams.
9. Raw style mutation primitives are not reachable from templates or the public app surface.
10. Tests enforce the ownership boundary, atomicity, actual UI wiring, failure rollback, session parity, and performance invariants.
11. In Circular batch, a semantic group or bulk rule/palette edit updates every affected Result in one History command; switching Result, saving, loading, and regenerating cannot reveal a stale sibling output.
12. A successful Generate or session export cannot adopt or download a mixed rule/palette revision.
13. Single and explicit-selection edits remain exact across duplicate public record IDs and same-record duplicate Features, or fail closed with a Generate-to-enable diagnostic for a legacy catalogue that cannot be upgraded safely.

## 2. Confirmed regression and current constraints

### 2.1 Scope behavior still exists but is bypassed

The current scope-aware entry point is `requestFeatureColorChange()` in `gbdraw/web/js/app/feature-editor/color-actions.js`. It still computes a matching rule, same-caption siblings, same rendered-label features, and same source-label features. `handleColorScopeChoice()` still applies all dialog choices, and the modal remains in `gbdraw/web/index.html`.

The active controls do not call that entry point:

| Current control | Current call | Required call |
| --- | --- | --- |
| Popup header fill | `setFeatureColorValue()` | scope-aware fill request |
| Popup Edit-tab fill | `setFeatureColorValue()` | the same scope-aware fill request |
| Feature list fill | `setFeatureColorValue()` | the same scope-aware fill request |
| Popup label text | `updateClickedFeatureLabelText()` direct override | scope-aware label request |

Before `7aab08f8`, popup fill called `updateClickedFeatureColor()`, which delegates to `requestFeatureColorChange()`, and the Feature list called `requestFeatureColorChange()` directly. The tri-state control replacement changed the template wiring without deleting the scope subsystem. This is a wiring regression, not an intentional product removal.

### 2.2 `Auto` represents inheritance but behaves like a lock

`ColorValueControl` maps `null` to `Auto`, `none` to `None`, and other values to `Color`. It disables the native color picker unless the mode is `Color`. `getFeatureColorValue()` returns only an explicit override and therefore returns `null` for an ordinary palette- or rule-derived feature, even though `getFeatureColor()` can resolve its visible color.

The generic tri-state contract remains valid for nullable settings. It is not the correct Feature-editing contract. Feature fill and stroke need an effective-style view model that separates:

- the value currently rendered;
- the origin of that value;
- whether a local explicit value exists;
- whether the user can reset to inheritance.

### 2.3 Feature and legend updates are not atomic

`setFeatureColor()` currently writes `featureColorOverrides` and changes the Feature DOM before it awaits legend work. Adding a legend entry may initialize main-thread Pyodide, calculate text dimensions, generate an SVG fragment, reflow the legend, and then re-extract `legendEntries` from the DOM. A missing runtime or `legendReflow` metadata can leave the feature and rule changed while the legend remains old.

This ordering creates an observable partial commit:

```text
feature override/rule -> Feature DOM -> await legend/Pyodide -> legend DOM/model -> Result flush
```

The target ordering is:

```text
intent -> pure scope plan -> user choice -> async preflight/staging
       -> one synchronous semantic/DOM commit -> one Result flush
```

No semantic or Result mutation may occur before preflight succeeds.

### 2.4 The same intent has too many writable representations

The current implementation can independently mutate:

- `manualSpecificRules`;
- `featureColorOverrides`;
- live Feature SVG attributes;
- `legendEntries` and `legendColorOverrides`;
- live legend SVG nodes;
- `results[].content`;
- `clickedFeature` display state;
- the cached canonical render request and its specific-color resource.

These are not equal owners. The target authority model is fixed in Section 3.4.

### 2.5 Public mutation APIs permit bypass

`feature-editor.js` and `app-setup.js` currently expose both orchestration methods and primitives such as `setFeatureColorValue()` and `setFeatureColor()`. Each is separately wrapped as an undoable action. A template change can therefore choose the wrong level without a build, type, architecture, or runtime failure.

### 2.6 Existing tests validate disconnected halves

`feature-color-actions.test.mjs` directly calls `handleColorScopeChoice()`. The existing Linear Playwright regression directly invokes `app.requestFeatureColorChange()` and `app.handleColorScopeChoice()` from `page.evaluate()`. Global control tests correctly verify the generic `Auto / None / Color` contract but do not exercise Feature controls.

All those tests can pass while the real popup is wired to a raw single-feature setter. The new test plan must go through the rendered controls and must also enforce the module boundary statically.

### 2.7 Circular batch has shared semantics but only one mounted Result

Circular batch owns one logical Result per input record, while the canonical request owns one shared specific-color resource. `featureStateFromCatalog()` currently flattens the Features from every catalogue item into one `extractedFeatures` list, but `PreviewRuntime` indexes and flushes only the mounted `selectedResultIndex`.

Changing a shared regex rule only in the active Result is therefore not a valid optimization: inactive Result strings remain old, while fresh regeneration applies the new rule to them. Converting a broad rule edit into active-Result exact hash rules is also not an equivalent operation. It freezes the current match set, expands the rule list, and cannot distinguish duplicate record instances because renderer hash and `record_location` selectors use public record IDs rather than the catalogue's unique `recordKey`.

The fixed contract is semantic rather than mount-driven: every Result affected by the selected durable intent participates in the transaction. D4 adds a renderer-recognized record-instance selector so `single` and explicit selections remain reproducible, but that does not turn a broad group into an active-only snapshot. Such a snapshot is a separate product capability and is not part of this work package; do not simulate it with ambiguous hashes, mass-generated exact rules, or persisted `resultIndex` values. The same rule applies to the Specific Rules editor, preset application, and TSV import: if they mutate live rules, they must stage every affected Result or be classified as Generate-time drafts that do not claim to update current Results.

### 2.8 Reactive and History work currently extends past the apparent commit

The deep `manualSpecificRules` watcher reapplies rules, refreshes caches, and reconciles legend state after the array changes. It can therefore create a second projection and dirty mark outside a nominal synchronous command. Watcher execution is not an invariant boundary and cannot remain an independent mutation owner.

`history.runUndoableCommand()` currently applies the command, awaits a new intent capture, and only then pushes the undo entry. Undo and Redo likewise mutate before their post-operation capture and stack transition. A capture failure can leave a visible mutation without the corresponding History movement. The transaction design must cover these post-apply seams, not only the Feature/legend DOM steps.

### 2.9 Generate and save can cross style revisions

Generation builds a canonical request before awaiting the render worker, later applies then-current live rules to the candidate Results, and finally adopts the earlier canonical basis. A style edit or Undo during a successful generation can therefore pair newer artifacts with older canonical semantics. Session export also crosses asynchronous resource assembly after reading Result and editor state.

`generationKey` prevents an old generation response from replacing a newer generation; it does not prove that style semantics stayed unchanged. Generate adoption and session download need an immutable start snapshot plus a semantic revision/fingerprint comparison at their final boundary.

## 3. Fixed product and architecture decisions

### D1. The Feature editor uses an effective-style control

Introduce a Feature-specific control/view model; do not change the meaning of the generic `ColorValueControl` used by global nullable settings.

The fill/stroke view model is:

```text
{
  effectiveColor: "#54bcf8" | "none",
  explicitValue: null | "none" | "#rrggbb",
  origin: "palette" | "specific-rule" | "local" | "svg-default",
  originLabel: "Inherited from palette" | "Inherited from rule: CDS" | ...,
  canReset: boolean,
  allowNone: boolean
}
```

Required interaction:

- The swatch/picker always displays `effectiveColor` and remains enabled while inherited.
- Picking a color promotes the effective value to an explicit color intent in one interaction.
- `No fill` emits `{ kind: "none" }`.
- `Use inherited` emits `{ kind: "inherit" }` and is disabled when already inherited.
- The word `Auto` is not used as a disabled Feature-editing mode.
- The popup header may keep a compact quick editor and the Edit tab may keep a full editor with origin/reset actions. Both are variants of the same Feature fill component, consume the same derived view model, and emit the same request action; neither owns state or a handler. The Feature list uses the same compact variant.
- Native picker `input` events update only the component-local draft swatch. The committed request is emitted once on `change`/picker acceptance. Opening or dragging a picker must not create multiple History entries or partial single-feature previews before scope selection.

Keep UI and persistence color domains separate. The native picker emits normalized six-digit hex. `RuleColor` accepts and round-trips the already supported 3/4/6/8-digit hex forms plus the explicit literal `none`; normalization must not reject a legacy value merely because the picker cannot author it. A no-fill Feature has no Feature-derived legend membership. A manual-only entry with the same caption remains because it has a different owner.

### D2. Semantic scope and Result extent are separate domain values

Define a stable scope vocabulary for fill first:

```text
single
matching-rule
feature-type
legend-caption
rendered-label
source-annotation-label
similarity-group
selected-features
```

Define the Result extent independently:

```text
current-result
all-affected-results
```

Each candidate contains a stable ID, semantic scope, Result extent, user-facing label, total target count, affected Result count, `targetsByResult`, stable target feature keys, and any rule identity needed for persistence. The planner returns candidates; the dialog renders the returned list. The template must not contain a separately maintained hard-coded inventory of scopes or infer extent from the mounted Result.

The transaction domain equals the durable semantic target domain:

- `matching-rule` changes the shared rule and therefore uses `all-affected-results`, including every Circular batch Result containing a match.
- `feature-type`, `legend-caption`, `rendered-label`, `source-annotation-label`, and `similarity-group` use every Result reached by that selected semantic group. They must not be silently narrowed because only one Result is mounted.
- `single` uses the origin Result and one stable Feature target. Its durable exact rule uses the renderer-recognized record-instance Feature selector defined in D4, not a public-record hash that can hit a duplicate input.
- `selected-features` uses the Results represented by the explicit selection. The current selection UI is active-Result-local, but the contract does not infer that from `selectedResultIndex`.
- In non-batch topologies, `all-affected-results` may contain one Result; keep the same contract rather than adding mode-specific planner branches.
- `current-result` plus a broad semantic scope means “materialize a snapshot exception for this output.” That is a separate, optional product action. This work package does not expose it: the D4 exact selector is reserved for `single` and explicit selections, and mass materialization would change future-match semantics and expand rules with target count. Do not approximate the omitted action with public record hashes, `record_location`, output position, or persisted `resultIndex`.

Rules:

- `single` is available only when the target has a validated `instanceHash`. A v40 catalogue that cannot be upgraded safely keeps viewing and group scopes but disables `single` and `selected-features` with `Generate to enable exact feature edits`; it never falls back to public `hash`.
- Group candidates are included only if they contain at least one other target.
- Two scopes may remain visible when they target the same current features but encode different durable intent, for example editing one regex rule versus materializing exact feature rules.
- Candidate ordering is deterministic: matching rule, Feature type, legend caption, rendered label, source annotation label, Similarity group, selected features, single. Within a future scope that supports more than one Result extent, batch-wide semantic intent precedes the explicitly named current-output snapshot exception.
- If any group candidate exists, the request returns `needs-scope`; it must not silently default to `single`.
- Cancel closes the dialog and restores the local control draft without changing semantic or artifact state.
- Existing-caption color reuse is a conflict-resolution choice, not a hidden scope. It is shown only when the requested caption already belongs to a different color.
- A candidate that spans multiple Results names both totals, for example `Update rule across batch — 84 features in 12 outputs`. A Result-local candidate names the Result and states that it creates exact exceptions. The UI never says only `Apply to all (N)` when more than one Result exists.

### D3. Request, planning, confirmation, and commit are separate phases

Use serializable, DOM-free data contracts. The first implementation is fill-specific; proposed JSDoc shapes are:

```text
FeatureFillValue =
  { kind: "inherit" }
  | { kind: "none" }
  | { kind: "color", color: "#rrggbb" }

FeatureFillIntent = {
  targetFeatureKey,
  value,
  requestedCaption,
  source: "popup" | "feature-list" | "selection-toolbar" | "legend-editor",
  originResultIndex,
  resultGenerationKey,
  documentEpoch,
  styleFingerprint
}

FeatureFillScopeCandidate = {
  id,
  semanticScope,
  resultExtent: "current-result" | "all-affected-results",
  targetCount,
  affectedResultCount,
  targetsByResult: [{ resultKey, resultName, featureKeys }],
  durableRuleIntent
}

FeatureFillPlan = {
  token,
  status: "noop" | "ready" | "needs-scope" | "conflict" | "invalid",
  intent,
  candidates,
  diagnostics
}

ResolvedFeatureFillPlan = {
  ...FeatureFillPlan,
  status: "ready",
  selectedCandidateId,
  semanticScope,
  resultExtent,
  affectedResultKeys,
  targetsByResult,
  semanticBefore,
  semanticAfter,
  legendProjectionByResult,
}
```

`planFeatureFillChange()` and `resolveFeatureFillPlan()` are pure: neither writes reactive state, DOM, Result content, History, or dialog state. A `needs-scope` plan has no single `semanticAfter`; resolving one candidate produces the immutable `ResolvedFeatureFillPlan` consumed by the command builder. Inputs include normalized catalogue items with durable Result keys plus current indexes, Result identity, normalized rule/palette state, document epoch, and style fingerprint. Tests deep-freeze them and verify no mutation. Persisted rule intent and History commands key Results by durable identity; array index is presentation/navigation data only.

The presentation facade exposes only request/resolve/cancel operations. Only the internal command builder can mutate. A reasonable final surface is:

```text
requestFeatureFillChange(intent)
resolveFeatureFillScope(planToken, scopeId)
cancelFeatureFillScope(planToken)
requestFeatureStrokeChange(intent)
requestFeatureLabelTextChange(intent)
```

The number of presentation methods is not the invariant. The invariant is that every accepted fill mutation reaches one internal `buildFeatureFillCommand(plan)`/commit boundary and no primitive setter is public. When label or stroke becomes the second proven consumer in Phase 5, extract only the common target/scope protocol and keep property-specific planning and persistence in their current owners.

### D4. Semantic ownership is explicit

Use the following authority table:

| Data | Authority | Treatment |
| --- | --- | --- |
| Base feature fill | Applied palette/default-color state | User-owned semantic input |
| Explicit feature fill/caption | Normalized `manualSpecificRules` | Sole semantic owner for renderer-visible explicit color intent |
| Effective feature fill/caption | Palette plus ordered specific-rule resolver | Derived value |
| `featureColorOverrides` | Stable-identity replay cache | Derived from semantic rules; never an independent UI write target |
| Feature-derived legend entries | Actually used effective `(caption, color, feature IDs)` projection | Derived; not independently authored |
| Manual legend-only entries | Document-global legend intent | User-owned editor state, projected into every Result in the document |
| Legend order and caption-keyed stroke/style | Document-global legend intent | User-owned intent; Result-specific transforms are derived layout |
| Live Feature/legend SVG | Mounted Result projection | Artifact |
| `results[].content` | Serialized per-Result projection | Artifact/cache; every affected entry must share the committed style revision |
| `clickedFeature` fields | Effective-style view model | Derived presentation state |
| Canonical specific/default-color resources | Projection of normalized rules and applied palette/default colors at the request/session boundary | Persisted renderer contract |

Extend the existing `specific-color-rules.js` normalization and legend-intent helpers rather than introducing a parallel rule model. The resolver order is a nested contract: for the exact Feature-type rule set and then the `*` rule set, try reserved instance literal, legacy/public-record `hash`, `record_location`, source qualifiers in source order, then location. Thus an exact-type legacy hash precedes a wildcard instance rule. Add a used-feature projection that applies the same resolver to committed catalogue items and returns groups per Result.

Before visibility or color resolution, build one immutable `FeatureInstanceIdentityPlan` over every materialized record and all source Features, including hidden/nested Features. A shared core owner, not Web-support post-processing, consumes normalized `record_key`, public record ID, source-coordinate stable Feature ID, and authoritative source-feature index. It applies the catalogue's existing collision rule to produce `biological_feature_id` (including the source-index discriminator for same-record duplicates), then freezes:

```text
instanceHash = "fi1_" + base32(
  sha256(domain("gbdraw-feature-instance-v1")
         + lengthFrame(record_key)
         + lengthFrame(biological_feature_id))[0:16]
)
```

Visibility, Feature drawing, used-rule/legend calculation, presence checks, Feature objects, interactive metadata, and catalogue construction consume that same plan. Ambiguous source identity fails closed. Ordinary rendering without an instance rule remains unchanged when no identity context is available; an instance rule requires the context and produces a capability diagnostic rather than being treated as a source qualifier.

The persisted TSV pseudo-qualifier is `__gbdraw_instance_hash__`, not `instance_hash`, so existing biological `/instance_hash=...` qualifier rules keep their old meaning. It is a newly reserved public extension: case-sensitive literal equality, fixed `fi1_...` grammar, machine-generated/opaque, and read-only in the generic Specific Rules UI. Python and Web parsers/normalizers/serializers distinguish this literal rule from regex rules; canonical requests containing it require schema 6. The token is stable within the originating canonical request/session record-key domain. Typed Python replay is portable when `RecordInput.record_key` is preserved; an ad-hoc CLI invocation is not promised to reproduce Web UUID record keys. Exact-row TSV download names the minimum compatible gbdraw version and warns about the record-key domain; older standalone CLI behavior cannot be made safe by the session schema gate. No CLI flag, public Python option, request field, or per-Result color-table is added, but the public TSV reference and compatibility notes are updated.

Keep feature-catalog schema 3 as an additive container: old catalogues may omit `instanceHash`, while the v41 writer requires it for every biological Feature in a non-empty Result set. Normal generation persists it once and Web never recomputes it. The v40 compatibility upgrader may derive it only from the already canonical `(recordKey, biologicalFeatureId)` pair using the frozen helper; if either value is missing/ambiguous, exact scopes stay disabled. Session v41 and canonical request schema 6 gate the reserved vocabulary so older v40/schema-5 consumers reject rather than silently interpret it as a source qualifier. Existing schema-3 Gallery/standalone assets remain readable; newly generated interactive metadata bytes intentionally change.

Manual-only legend intent is document-global, matching the existing single editor-state owner. Every Result projection merges that global intent with its own used Feature groups and computes its own layout. `legendEntries` remains only the selected Result's derived view. v40 import treats the saved editor legend intent as global; conflicting inactive artifact-only legend state is not promoted to authority and requires regeneration before batch editing/saving rather than being guessed.

The legend editor must classify the selected entry before mutation. Changing the color or caption of a Feature-derived entry emits a fill intent with `legend-caption` scope and therefore changes the semantic rules through the same command. Only an entry with no Feature members and no used rule may use the manual legend-only editor owner. `legendColorOverrides` must not remain an alternative color owner for Feature-derived entries.

Inventory every writer before cutover. Popup, Feature list, selection toolbar, Feature-derived legend editing, the Specific Rules editor, preset application, TSV import, applied palette/default-color actions, session hydration, Reset, and History restore must each be classified as either:

- a live semantic edit routed through the owning command and style revision;
- a document-replacement/hydration owner that installs an already validated snapshot; or
- a Generate-time draft that is deliberately outside live Result semantics.

Specific Rules field acceptance/add/remove/reorder, preset apply, TSV import, and applied palette/default-color acceptance are live semantic edits. They use one bulk semantic-snapshot command that computes before/after effective style and visibility for all catalogue Features, stages every affected Result, updates both canonical color resources as applicable, and creates one History entry. Input/picker typing stays component-local. Session hydration and Reset are document replacements. No watcher, import callback, palette helper, preset helper, or legacy setter may remain an unclassified live writer.

### D5. One atomic command owns semantic and visual convergence

`buildFeatureFillCommand(plan)` returns the existing History command shape: `apply()`, `revert()`, compensation restore, metadata, and an exact upper-bound byte estimate. The bulk rule/palette writer builds the same projection command from complete before/after semantic snapshots. Apply and revert re-run preflight against durable Result identities and the Result mounted at that moment; they do not retain DOM nodes or assume that the apply-time Result remains selected.

Preflight must complete before the first active mutation:

1. Verify document epoch, Result generation key, expected before/after style fingerprint, stable target identities, affected durable Result keys, and current mounted-root ownership. A mounted root is required only when its Result is affected.
2. Re-resolve `targetsByResult` against the current catalogue items and reject a stale, missing, or cross-Result target.
3. Build the complete next normalized rule and applied palette/default-color snapshot.
4. Derive the next override replay cache and before/after visibility/effective-style diff.
5. Derive used Feature legend groups per affected Result and merge them with the document-global manual legend intent.
6. Stage every Result-specific legend diff and layout measurement outside the active Result.
7. If the currently mounted Result is affected, build forward/inverse patches and a non-mutating mounted serialization candidate. For each other affected Result, parse its already admitted SVG into a detached root, resolve Feature/legend targets, apply the staged projection, and serialize one candidate. Do not mount it, install handlers, or re-sanitize it.
8. Use the Result-ingestion owner to create fresh admitted metadata for changed Result objects; do not share the private mutable admission/mount record through object spread. Stage one complete next `results` array plus a separate `validatedStyleFingerprintByResultKey` ledger, retaining unaffected Result objects/private metadata and the exact previous objects/runtime ledger for compensation and History accounting.
9. Validate rule precedence, before/after visibility coverage, target coverage, caption/color uniqueness, Result count/order/name, per-Result legend membership, admission ownership, and post-commit fingerprint invariants. Results whose bytes do not change are still proven compatible with the new semantic snapshot and may retain object identity.

The planner may await detached-Result preparation and browser font readiness. Routine fill changes do not rerender. If a bulk rule/palette change adds or removes rendered Feature geometry or admitted SVG metadata cannot support the required patch, the owning worker may render complete candidates during preflight; no semantic change is committed until every candidate succeeds. This required-geometry transaction is the explicit exception to the existing Web live-edit fallback, so update `gbdraw/web/CLAUDE.md`: its “keep a direct edit after optional rerender failure” rule continues for non-atomic optional reflow, but not for this all-or-nothing style command. The UI reports progress and remains cancellable until all candidates are ready. A newer style request, Generate, Result replacement, session load, or Reset invalidates the whole staged transaction. Performance pressure must never silently narrow `all-affected-results` to the mounted Result.

Commit is synchronous and bounded:

1. Replace/update the normalized semantic rules and/or applied palette/default-color state. Remove fill/legend projection from the deep rule/palette watchers; synchronous suppression flags are not an invariant because Vue may drain them in a later microtask.
2. Replace the derived override cache.
3. If the current mount is affected, apply its prepared Feature and legend patches and serialize it once into a candidate without writing `state.results`.
4. Replace the selected `legendEntries` view only when the selected Result is affected; otherwise keep its projection unchanged.
5. Assign the fully staged `results.value` array exactly once, then finalize PreviewRuntime and Result-admission metadata against the committed object. No call to the existing state-writing `flushActiveResult()` occurs inside this transaction.
6. Run the explicit rule/cache/legend reconciler once; no deep watcher replays the transition later.
7. Refresh the popup/list effective-style view.
8. Advance the semantic revision/fingerprint and atomically replace the transient Result-validation ledger: changed Results contain new bytes, while unchanged Results were proven unaffected. Do not mutate private admission metadata on unchanged Result objects merely to stamp a revision.

No `await` occurs between steps 1 and 8. If an unexpected exception occurs, apply inverse mounted-DOM operations and restore semantic snapshots, the prior Result array/object and admission-metadata references, Result-validation ledger, active legend state, PreviewRuntime dirty/index/mount state, and exact prior revision before returning failure. Do not publish a partial Result array. Serialization may clone the active SVG internally; the fast-path invariant is zero active-root parse, sanitize, remount, Feature-index rebuild, or handler setup. With `r` affected Results and `a` equal to 1 only when the current mount is affected, apply/revert performs at most `a` mounted serializations and `r-a` detached parse/index/serialize passes.

### D6. Editor legend updates do not initialize main-thread Pyodide

Existing generated legend geometry remains authoritative. Ordinary editor changes must use browser-owned measurement in the legend module:

- Updating an existing swatch color requires no text measurement.
- Adding or renaming an entry uses a reusable offscreen SVG measurement host with the same computed font and `getComputedTextLength()`/`getBBox()`.
- The measurement host is outside the active Result and is removed or reused without entering saved SVG content.
- Font readiness is awaited before measurement; the host must be connected in the same document/style context long enough for browser SVG metrics to be authoritative.
- Existing `legendReflow` metadata is validated independently for every affected Result during preflight.
- A missing font/measurement/metadata result is an explicit preflight error; it does not fall back to a partial feature edit.
- Python rendering and Worker-owned generation remain unchanged.

The editor path must initialize main-thread Pyodide zero times for fill, `None`, inherit/reset, caption-group edits, and single-feature legend creation.

### D7. History starts only after scope resolution

Request, dialog display, cancel, and local picker draft changes are not undoable semantic operations. Remove individual `undoableAction()` wrappers from raw request/setter/choice functions.

After the user confirms a scope, call:

```text
history.runUndoableCommand("Change feature fill", () => buildFeatureFillCommand(resolvedPlan))
```

Requirements:

- One successfully committed interaction creates exactly one History item. Before apply, reject a command whose own retained before/after data exceeds History's byte cap; show a specific “edit is too large to undo” error and mutate nothing. Older entries may be evicted only after the new command itself is known to fit.
- Cancel and no-op create none.
- `history.js` is a required owner in the atomic-command phase. The History transaction does not end when `command.apply()` returns; it ends only after post-apply intent capture succeeds and the stack transition is committed.
- If post-apply intent capture or any stack/file-retention/bookkeeping step fails, `runUndoableCommand()` invokes compensation mode, not normal Undo. Compensation restores the exact prior semantic revision, fingerprints, Result/admission references, current intent, retained-file IDs, byte totals, and both stacks. If compensation itself fails, block further History actions and report an explicit integrity error rather than claiming that the edit failed cleanly.
- Undo and Redo update rules, derived caches, Feature DOM, per-Result legend projections, all affected Result content, and semantic revision through the same atomic path. If their post-operation capture/bookkeeping fails, use compensation mode to restore the exact pre-attempt state before leaving the stacks unchanged.
- Normal Undo validates the command's expected after-fingerprint and normal Redo validates its before-fingerprint; each advances revision as a new accepted transition. Session load, Reset, Generate adoption, and Result replacement advance document epoch and clear/prune commands from an older epoch so a stale top entry cannot block all later Undo.
- Selection-toolbar bulk edits use an explicit `selected-features` scope and the same command builder.
- `estimateBytes()` includes retained before/after content and admission metadata for every affected Result plus semantic deltas. History's action/byte limits remain authoritative; broad batch edits may evict older entries with the existing visible diagnostic, never by under-reporting their size.

### D8. Canonical request and session output cannot remain stale

`services/session-request.js` remains the only Web-to-canonical projection boundary. Add a narrow copy-on-write projection for editor-owned color semantics. It overlays current normalized specific-color rules and applied palette/default colors onto the last committed canonical basis without adopting unrelated drafts or rerunning analysis.

Track a transient semantic revision/fingerprint:

- successful style command, Undo, and Redo advance it;
- a separate transient ledger validates every durable Result key against the style fingerprint; only Results with changed bytes receive new objects/admission metadata;
- successful Generate/adoption aligns the committed canonical revision only to the exact style snapshot projected into both the adopted canonical resource and every candidate Result;
- save/export-to-session and regeneration obtain an effective canonical basis with both current color-resource overlays;
- an unchanged revision reuses the existing canonical request/resources;
- only the changed specific/default-color resource and its reference are replaced;
- Linear comparison resources, sequence resources, LOSAT evidence, topology, and unrelated draft settings retain identity.

Generate and session export use the revision as a concurrency token, not only as a dirty bit:

1. Capture an immutable style snapshot at operation start: normalized rules, applied palette/default colors, fingerprint/revision, document/Result generation keys, and canonical basis identity.
2. Use that snapshot consistently for request construction and initial candidate projection; do not reread live rule arrays midway through an asynchronous operation.
3. At Generate adoption, compare with the latest committed style revision. If it advanced, first prove candidate coverage under the newer rules/palette, including Feature visibility. When topology is unchanged, stage the same copy-on-write overlay against candidate Results and both canonical color resources, then recheck immediately before synchronous adoption. If required Feature nodes differ, discard/restart generation from the newer snapshot; never pair a newer resource with an unpatchable older SVG. A second advance invalidates and restarts staging.
4. At session export, build Results, editor state, canonical overlay, and resources from one captured revision. Recheck before download. On mismatch, discard the envelope and restart from the newer snapshot or fail explicitly; never download a mixed-revision session.
5. Failed, cancelled, or stale Generate/export work does not advance or align revisions.

Two explicit compatibility changes are in scope:

1. **Confirmed v40 stale-overlay recovery.** The first-parent `main` Gallery session `gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json` is v40 and supplies concrete evidence: four broad `config.rules`, 66 derived `featureColorOverrides`, saved SVG agreement, and a stale canonical default-color resource. Its canonical specific-color table matches `config.rules`. Phase 0 minimizes the default-color mismatch into a frozen fixture; this branch is mandatory. A separate specific-table mismatch is tested only as a derived compatibility case, without attributing it to the Gallery producer.
2. **Exact-selector capability gate.** The new writer emits session v41 and canonical request schema 6. Request schema 6 reserves `__gbdraw_instance_hash__`; readers of schema 5 continue treating unknown keys as ordinary source qualifiers. Catalogue schema 3 remains additive as defined in D4. Existing v40/schema-5 sessions are promoted in memory only after the rules/identity checks below.

The v40 stale-overlay recovery rule is:

1. Detect a v40 envelope whose canonical specific-color table differs from stored editor `config.rules`, or whose canonical default-color table differs from the stored applied palette/default-color state.
2. Normalize the complete editor rules, including broad and exact selectors, and the applied default colors.
3. Re-resolve the rules and defaults against the persisted catalogue and compare every derived override, used caption/color, and saved Result Feature fill. Accept an overlay only when all observed state is explained, no target is ambiguous, and no caption/color conflict exists; do not require one exact rule per derived cache entry.
4. Materialize the recovered resource and canonical request identically in Web and the full-envelope typed Python boundary (`session_to_request()`/session compatibility), then promote the in-memory document to the current contract.
5. Separately upgrade missing v40 `instanceHash` values only from validated canonical `(recordKey, biologicalFeatureId)` pairs. If this fails, preserve viewing and group-rule editing, disable exact scopes, and require Generate before a v41 save containing Results.
6. If stale-rule validation fails, keep the prior canonical request and saved Result for viewing and emit a non-mutating compatibility diagnostic. Do not guess.

On load, persisted numeric revision stamps do not exist. Mark all Results untrusted, validate/reproject them against the recovered rules/palette and admission metadata, and only then stamp the current runtime fingerprint and enable live edits. Session export captures one logical fingerprint; it does not persist transient revision counters as a second authority.

### D9. Fill is the first vertical slice; label and stroke close the same bypass class

Complete fill end to end before generalizing. The fill slice includes popup, Feature list, selection toolbar, Feature-derived legend color/caption actions, Specific Rules editor field/add/remove/reorder, preset/TSV import, and applied palette/default-color acceptance. The latter bulk writers use the semantic-snapshot command from D4/D5. These consumers move before raw fill/rule/palette setters are deleted; none is deferred to the later label/stroke phase.

Then use the proven scope planner/command protocol for the two adjacent regressions:

- Popup label text must call the existing scope-aware label request rather than writing a per-feature override directly. Rendered-label/source-label target resolution uses the shared scope-plan representation; label-specific persistence remains in the label owner.
- Feature stroke uses the effective-style control so inherited stroke is editable without `Auto`. The existing same-legend and selected-feature stroke operations become explicit scopes handled by a stroke command. Stroke remains editor-artifact state unless and until the renderer contract owns it; do not invent a Python/CLI stroke option.
- Feature visibility already follows a scope-aware command pattern. Treat it as the reference and keep its public request/resolve boundary; do not refactor it merely for naming symmetry.

Do not build a generic command framework before at least fill and one of label/stroke use the same abstraction. Any new shared helper must replace the superseded paths in the same phase.

### D10. Async intent is generation-bound

Every pending plan and preflight carries the origin Result key plus navigation index, ordered affected durable Result keys, document epoch, `resultGenerationKey`, starting style fingerprint, and a monotonic request token.

- Result switch, session load, Reset, new generation, popup target change, and disposal invalidate pending plans.
- A newer style request supersedes an older preflight.
- A stale or out-of-order preflight result is discarded before commit and cannot apply rollback data over a newer committed command.
- The dialog stores stable data, not mutable rule or feature object references.
- The UI reports a concise failure and keeps the prior effective value when a request becomes stale or preflight fails.
- Result switch may close a pending dialog, but switching after a committed batch-wide edit does not undo or narrow that edit. Apply/revert re-resolves the current mount: Undo/Redo after switching to an affected or unaffected Result restores all command Results without changing the selection.

### D11. Performance complexity is part of correctness

For `R` total Results, `r` affected Results, `n` catalogue Features, `m` ordered rules, `k` changed Features, total serialized bytes `S_all` across affected Results, and `a ∈ {0,1}` indicating whether the current mount is affected:

- without a proven rule index, planning/effective-style comparison has the honest upper bound `O(n*m)`; exact literal selectors should use a map and regex rules should be partitioned by Feature type, with structural counters recording actual evaluations;
- a `matching-rule` edit changes the durable shared rule; it must not materialize `O(k)` exact rules merely to preserve active-only runtime performance;
- mounted DOM updates are `O(k_active + changed active legend entries)` using the existing Feature DOM index;
- the Result transaction performs one shallow `O(R)` array copy/assignment, retaining unaffected object identity;
- it performs at most `a` mounted SVG clone/serializations and `r-a` detached parse/index/legend/serialize passes, for total work bounded by `O(R + n*m + k + S_all)` before any explicitly required geometry rerender;
- ordinary patchable edits perform zero admission/sanitization/remount/handler-setup passes, zero mounted Feature-index rebuilds, and zero main-thread Pyodide initializations;
- staged and retained History bytes for broad batch edits are measured honestly. Unchanged Result objects retain identity, and no second full candidate copy is kept after commit beyond what rollback/History requires;
- picker-local `input` events perform no catalogue scan or Result serialization;
- one accepted `change`/scope confirmation performs one command.

The implementation must not regress the direct-edit ownership and Result-ingestion guarantees established by the performance remediation plan. Batch-wide preflight exposes progress and remains cancellable before commit. Structural counters are mandatory; wall-clock measurements diagnose scaling and set evidence-backed thresholds after the representative batch fixture exists.

## 4. Non-goals

- No new CLI flag, public Python option, request field, or per-Result color-table model. The reserved public TSV pseudo-qualifier, request schema 6, session v41, and additive catalogue metadata defined in D4/D8 are the only renderer/persistence extensions.
- No second render path and no return of heavy generation to main-thread Pyodide.
- No general-purpose arbitrary SVG patch engine.
- No routine regeneration for patchable fill/stroke edits. Required Feature-geometry changes may use the existing Worker during atomic preflight and must leave all state unchanged on failure.
- No full rewrite of the legend editor or History service.
- No change to generic `ColorValueControl` semantics for nullable global settings.
- No visual/geometry reference regeneration merely for the editor fix. New interactive catalogue metadata bytes intentionally change and require focused metadata fixture/reference review; existing Gallery SVGs remain valid and need not be regenerated solely to add optional metadata.
- No Gallery session/media regeneration solely because the editor interaction changed. Audit documented popup screenshots and recapture only those that visibly become incorrect.

## 5. Target module boundaries

### 5.1 Proposed new focused modules

Names may be adjusted to current repository conventions, but responsibilities must remain separate.

| Module | Responsibility | Forbidden responsibility |
| --- | --- | --- |
| `app/feature-editor/fill-scope-plan.js` | Pure stable-target fill-scope discovery and plan normalization | Vue state, DOM, History, serialization |
| `app/feature-editor/fill-view-model.js` | Resolve fill effective value, origin, reset capability | Mutation or dialog policy |
| `app/feature-editor/fill-command.js` | Build/preflight/apply/revert one fill command | Template state and generic settings |
| `app/feature-editor/style-snapshot-command.js` | Build bulk rule/palette before/after commands over all affected Results | UI field drafts or renderer identity |
| `app/feature-editor/fill-result-projection.js` | Stage mounted patches and detached candidates for the affected Result set | Rule semantics, Result admission, or History |
| `app/feature-editor/fill-color-control.js` | Feature-specific Vue fill control and local draft | Scope discovery or semantic mutation |
| `app/legend/feature-fill-projection.js` | Derive used Feature legend groups per Result; stage legend projection | Owning manual color rules or committing Results |
| `features/instance_identity.py` | Build the pre-render `FeatureInstanceIdentityPlan` shared by renderer/catalogue consumers | Web UI or session policy |

Prefer extending `specific-color-rules.js`, `preview-runtime.js`, and existing legend helpers over copying their normalization, identity, or flush logic into these modules.

### 5.2 Public surface after completion

`createFeatureEditor()` and `app-setup.js` may expose presentation queries and orchestration actions, but not mutation primitives.

Allowed examples:

- `getFeatureFillViewModel`
- `requestFeatureFillChange`
- `resolveFeatureFillScope`
- `cancelFeatureFillScope`
- `requestFeatureStrokeChange`
- `requestFeatureLabelTextChange`
- existing visibility request/choice actions

Forbidden examples:

- `setFeatureColor`
- `setFeatureColorValue`
- `updateClickedFeatureColor`
- direct rule upsert/remove helpers
- direct feature/legend DOM attribute setters
- direct per-feature label override setter

Low-level functions may remain as module-private implementation details only when the command builder is their sole caller.

### 5.3 Successful fill data flow

```text
FeatureColorControl / Feature list
  -> requestFeatureFillChange(intent)
  -> planFeatureFillChange(snapshot, intent)
      -> noop/invalid
      -> needs-scope -> dynamic dialog -> resolve scope
      -> ready
  -> history.runUndoableCommand(...)
  -> buildFeatureFillCommand(resolved plan)
  -> preflight rules/palette + targetsByResult + derived cache
     + mounted DOM patch + detached inactive Result candidates
  -> synchronous semantic/DOM/Result commit
  -> zero/one non-mutating mounted serialization + one array-level Result swap
  -> History intent capture + stack transition, or compensating revert
  -> effective canonical semantic overlay marked current
```

## 6. Surface and authority matrix

| Surface | Planned change | Required evidence |
| --- | --- | --- |
| Typed Python renderer | Build shared pre-render instance identities; consume schema-6 reserved literal TSV selector | Visibility/fill/legend/catalogue agree; legacy selectors and schema-5 behavior remain unchanged |
| Python/CLI session compatibility | Recover the confirmed v40 broad-rule stale overlay; read v40 and write v41/schema 6 | Main-backed frozen fixture produces the same canonical rules as Web; unsafe recovery fails closed |
| Feature catalogue/interactive SVG | Keep schema 3 additive; persist optional/required-by-v41 `instanceHash` and read old/new assets | v40 upgrade, v41 required-field, standalone schema-3, and intentional metadata-byte diffs tested |
| Web state | Replace raw dialog payload with stable plan/token/status; add document epoch, style revision/fingerprint, and Result-keyed validation ledger | Reset/session/Generate invalidate pending plans; transient fields are not persisted |
| Feature UI | Dedicated effective-style control; one popup editor; inline list routes to same request | Real DOM interaction shows inherited editable value and opens expected scopes |
| Scope planning | Pure stable-target plan with semantic scope, Result extent, and `targetsByResult`, shared by fill and later label/stroke | Table-driven, mutation-free unit matrix including Circular batch |
| Semantic color rules/palette | Normalize once; Feature and bulk snapshot commands mutate; no watcher owns convergence | Cross-runtime precedence, exact identity, bulk rule/preset/import/palette batch gates, complete writer allowlist |
| Override replay cache | Rebuild/validate from rules after command/load | Cache cannot diverge from normalized rules |
| Legend | Used Feature groups are derived per Result; document-global manual intent is merged; staged browser measurement | All affected SVGs agree; zero Pyodide initialization; legacy conflict fails closed; rollback on any Result failure |
| PreviewRuntime/Result projection | Non-mutating mounted candidate plus detached candidates and fresh admission metadata, committed by one array assignment | `a` mounted serializations, `r-a` inactive passes, no partial swap or shared private metadata |
| History | Confirmed command plus post-operation capture/compensation | One edit = one entry; capture failure restores artifacts and stacks; Undo/Redo parity across Results |
| Canonical request | Schema-6 copy-on-write specific/default-color overlay with style snapshot/revision | Changed resources only; topology coverage checked before Generate rebase |
| Session | v41 one-fingerprint writer envelope; confirmed v40 recovery and exact-capability upgrade | Edit-during-export cannot download mixed state; v40/v41 parity and diagnostics |
| Gallery/docs | Audit only visible documented workflow | No generated asset change unless a named capture is intentionally refreshed |

## 7. Implementation phases

### Phase 0: Freeze evidence and re-audit the shared baseline

Production changes: none.

Tasks:

1. Record `HEAD`, branch, `git status --short`, and production/test/documentation/generated diffs separately.
2. Confirm that the right-drawer and lazy-export changes now present at the planning baseline remain authoritative; inspect any newer overlapping edits in `index.html`, `app-setup.js`, `right-drawer.js`, `run-analysis.js`, `config.js`, `architecture-contracts.test.mjs`, and Playwright specs before editing.
3. Re-run the lambda Linear browser reproduction through actual controls and capture machine-readable before-state:
   - inherited CDS control presents `Auto` and disabled picker;
   - fill change does not open scope dialog;
   - Feature/rule changes precede legend convergence;
   - main-thread Pyodide initialization count;
   - History entries, PreviewRuntime flushes, root mounts, and Result updates.
4. Add a two-output Circular batch diagnostic reproduction. Record catalogue target ownership, scope counts, the mounted Result, inactive Result bytes, Result switching, saved session Results, and regeneration. Prove the current shared-rule/active-runtime seam without modifying production code; a focused test harness may supplement, but not replace, the visible-control reproduction.
5. Use controlled promises/failure hooks to freeze the current History post-apply capture seam, the deep-rule-watcher second projection, successful Generate plus intervening style edit/Undo, and session export plus intervening style edit.
6. Freeze the confirmed v40 evidence from first-parent `main`: minimize `BGC0000708-BGC0000713.gbdraw-session.json` while retaining the default-color mismatch, broad rules, derived overrides, and saved SVG agreement; record producer commit/build, hash, and expected recovered request. Add one deliberately ambiguous negative mutation. Keep the already-consistent canonical specific-color table in the provenance record.
7. Run current focused tests and record that they pass despite the browser regression.
8. Verify Node and Python Playwright availability and use the repository functional/performance configs. Escalate Chromium sandbox permission when required.

Phase gate:

- The reproduction uses real DOM events, not `window.__GBDRAW_APP__` mutation methods.
- The Circular batch evidence names every affected Result and demonstrates why mounted-only projection is insufficient.
- The main-backed v40 default-color positive and derived ambiguous-negative fixtures are reproducible. A specific-color mismatch case, when tested, is labeled as derived rather than main-backed.
- History, watcher, Generate, and export interleavings have reproducible failure boundaries or explicit evidence that the suspected seam is already fixed at the execution baseline.
- Unrelated worktree changes are inventoried.
- No production file is modified.

### Phase 1: Add pure effective-style and scope planning

Primary files:

- new `gbdraw/web/js/app/feature-editor/fill-view-model.js`
- new `gbdraw/web/js/app/feature-editor/fill-scope-plan.js`
- `gbdraw/web/js/app/feature-editor/rule-actions.js`
- `gbdraw/web/js/app/specific-color-rules.js`
- `gbdraw/web/js/app/file-imports.js`
- `gbdraw/web/js/services/feature-catalog.js`
- `gbdraw/web/js/services/standalone-interactivity.js` and generated source owner `standalone-interactivity-assets.js`
- `gbdraw/web/js/state.js`
- new `gbdraw/features/instance_identity.py`
- `gbdraw/core/record_metadata.py`, `gbdraw/features/ids.py`
- `gbdraw/features/selector_values.py`
- `gbdraw/features/colors.py`, `gbdraw/features/visibility.py`, `gbdraw/features/factory.py`, `gbdraw/core/sequence.py`
- Circular/Linear precalc and assembly call boundaries
- `gbdraw/render/interactive_context.py`, `gbdraw/render/interactive_svg.py`
- `gbdraw/web_support/feature_metadata.py`, `gbdraw/web_support/feature_catalog.py`, and request-render metadata owner
- `gbdraw/session_io.py`, request/session codecs, and authority validators for v41/schema 6
- `docs/REFERENCE/input-formats-and-tsv-schemas.md` and `docs/REFERENCE/interactive-svg-and-semantic-hooks.md`

Tasks:

1. Factor catalogue collision identity into the pre-render `FeatureInstanceIdentityPlan` and freeze `fi1_...` token construction. Build it once for all source Features and thread it through visibility, fill, used legend, presence, Feature-object, interactive, and catalogue consumers.
2. Add reserved literal `__gbdraw_instance_hash__` parsing/matching/serialization under request schema 6; preserve biological `instance_hash` qualifier and every legacy schema-5 selector result. Extract the pure Web resolver with exact-type-then-wildcard nesting and cross-runtime oracle.
3. Return origin metadata separately from explicit value.
4. Define stable intent/value/semantic-scope/Result-extent/plan normalization, separate picker/RuleColor validation, and stable feature/Result keys.
5. Preserve durable Result key plus current navigation index on every flattened catalogue Feature and compute `targetsByResult` from immutable feature/rule snapshots; never persist array index as selector identity or retain mutable objects in dialog state.
6. Add deterministic target-set and rule-identity signatures.
7. Add conflict planning for an existing caption with a different color.
8. Apply the fixed extent policy: semantic group scopes cover every affected Result; `single` and the current selection use their explicit target set. Do not emit a broad `current-result` exception candidate in this work package.
9. Persist additive catalogue `instanceHash`; implement v40 canonical-pair upgrade and fail-closed exact-scope capability. Update v41/schema-6 codecs, session authority, standalone readers, TSV warnings, and references without regenerating old Gallery SVGs.
10. Replace the field-heavy `colorScopeDialog` runtime representation with `pendingFeatureFillPlan`, token, candidates, progress, and status. Keep it transient and reset it on lifecycle boundaries.
11. Do not route production UI through the planner until the fill command in Phase 2 can commit every outcome atomically.

Phase gate:

- Pure planner unit matrix passes.
- Web/Python rule-precedence oracle cases agree, including exact-type versus wildcard, reserved literal, duplicate records, and same-record duplicate Features; visibility, fill, used legend, metadata, and catalogue resolve the same identity.
- v40 missing-token upgrade succeeds only from canonical pairs; otherwise exact scopes are disabled. v41/schema-6 round trip and schema-3 old/new standalone reads pass.
- Circular batch plans report exact Feature and Result counts and never flatten away target ownership.
- Deep-frozen inputs remain unchanged.
- No new mutable color owner is introduced.
- No visual geometry/fill changes when the reserved selector is absent; intentional catalogue/interactive metadata-byte changes are isolated and reviewed. No editor UI behavior changes yet.

### Phase 2: Implement staged legend projection and atomic fill command

Primary files:

- new `gbdraw/web/js/app/feature-editor/fill-command.js`
- new `gbdraw/web/js/app/feature-editor/style-snapshot-command.js`
- new `gbdraw/web/js/app/feature-editor/fill-result-projection.js`
- new `gbdraw/web/js/app/legend/feature-fill-projection.js`
- `gbdraw/web/js/app/feature-editor/color-actions.js`
- `gbdraw/web/js/app/legend/entry-actions.js`
- `gbdraw/web/js/app/preview-runtime.js`
- `gbdraw/web/js/app/watchers.js`
- `gbdraw/web/js/services/history.js`
- `gbdraw/web/js/services/svg-result-ingestion.js`
- `gbdraw/web/CLAUDE.md`

Tasks:

1. Build next rule/palette/default-color snapshots without mutation; preserve rule order and provenance.
2. Implement `color`, `none`, and `inherit` as plan values. Inherit removes only the rules owned by the selected plan/scope and reveals the next effective palette/rule value; it must not remove an unrelated broader rule.
3. Derive the stable override replay cache from the next rules and catalogue without losing Result ownership.
4. Derive used Feature legend groups with Feature IDs per Result; merge document-global manual-only intent and compute Result-local layout.
5. Split current legend reconciliation into prepare and commit stages and await real browser font readiness before measurement.
6. Replace editor legend text sizing/generation through main-thread Pyodide with an offscreen browser SVG measurement owner.
7. Build forward/inverse patches for whichever affected Result is mounted at apply/revert time plus detached candidates for the other affected Results, keyed by durable Result identity.
8. Add a non-mutating mounted serializer and one Result-array transaction API. Use the ingestion owner to replace admitted content/metadata without sharing its private mutable runtime record; preserve unaffected identity and restore mount/dirty/index/admission state on failure.
9. Remove fill/legend projection from deep rule/palette watchers. Hydration/replacement calls the same explicit reconciler; no microtask suppression flag owns correctness.
10. Implement one synchronous semantic/DOM/Result commit and exact rollback on injected failure in any Result.
11. Make History intent capture, push/clear/eviction, byte totals, and retained-file bookkeeping compensatable for apply/Undo/Redo. Pre-reject a single oversized command and epoch-prune stale commands.
12. Route single/group/selection-toolbar/Feature-derived-legend internals through the command behind private temporary adapters. Do not delete adapters needed by the old UI until Phase 4.
13. Route Specific Rules field/add/remove/reorder acceptance, preset, TSV import, and applied palette/default-color acceptance through the bulk semantic-snapshot command; stage geometry rerender only when before/after visibility requires it.
14. Clarify the repository Web contract so required atomic geometry failure rolls back, while unrelated optional-reflow behavior retains its existing fallback.

Phase gate:

- Command/legend/failure unit tests pass.
- The new fill implementation has one internal mutation owner and no new caller uses `setFeatureColor()`. Any legacy public adapter still needed by the pre-cutover UI is an explicitly unlanded transition and is removed before the non-separable Phase 2–4 slice can pass.
- An existing legend update and a newly suffixed single-feature legend both commit without Pyodide.
- Circular batch Feature scopes and bulk rule/preset/import/palette edits change every affected Result and no unaffected Result in one command; Result switching exposes no stale output.
- Undo/Redo after switching to an affected or unaffected sibling re-resolves the mount and preserves selection.
- A forced failure in any mounted/inactive Result, admission metadata, watcher drain, array assignment, or post-operation History bookkeeping preserves exact before-state, revision, retained files, bytes, and stack depths.
- A command larger than the History cap is rejected before mutation.
- Single and batch structural/performance counters remain within D11, and History byte estimates include batch candidates.

### Phase 3: Establish canonical/session revision safety before UI exposure

Primary files:

- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-authority.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/session.py`, `gbdraw/api/session_compat.py`, `gbdraw/session_io.py`, and request/session codec owners
- `docs/SESSION_COMPATIBILITY.md`

Tasks:

1. Separate document epoch, monotonic semantic revision, content fingerprint, and a transient durable-Result-key validation ledger; normal commands/Undo/Redo advance revision, compensation restores it exactly, and load starts untrusted.
2. Add copy-on-write specific- and default-color overlays in `session-request.js`.
3. Capture immutable rule/palette snapshots for Generate and session export; construct request, candidate projection, Result envelope, and both canonical color resources from one snapshot.
4. Revalidate at Generate adoption. Rebase only when Feature topology/coverage is unchanged; otherwise restart/reject from the latest snapshot.
5. Revalidate before session download. Discard/restart or fail a stale assembly; a cached committed request must not bypass a newer semantic revision.
6. On session load, mark all Results untrusted, restore normalized rules/palette, recover/upgrade the catalogue where safe, validate or reproject every Result, then stamp runtime fingerprints and enable edits.
7. Implement the confirmed broad/exact v40 stale-overlay recovery at both Web and full-envelope typed Python boundaries. Promote safe documents to v41/schema 6; diagnose rather than guess on mismatch.
8. Update session/request authority, v41/schema-6 codecs, compatibility docs, and exact-TSV portability docs without adding duplicate semantic fields.
9. Verify inactive mode drafts, Linear comparison plans, LOSAT resources, sequence resources, and unchanged Result object identity remain untouched by the style overlay.

Phase gate:

- Internal command harness: edit -> save -> load -> regenerate preserves rules, Feature fills, captions, and per-Result legends in Linear, Circular single/grid, and Circular batch.
- Edit/Undo during successful Generate and edit during export cannot produce mixed revisions.
- Web and typed Python projection agree on the frozen main-backed v40 broad-rule fixture and both reject its ambiguous mutation.
- The v41/schema-6 writer emits one-fingerprint rule/palette/Result envelope or refuses save with a specific error; it never silently saves stale semantics or missing exact capability.

### Phase 4: Replace UI wiring and remove `Auto` locking

Primary files:

- new `gbdraw/web/js/app/feature-editor/fill-color-control.js`
- `gbdraw/web/js/app.js`
- `gbdraw/web/index.html`
- `gbdraw/web/js/app/feature-editor.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/js/app/legend.js` and `gbdraw/web/js/app/legend/entry-actions.js`
- `gbdraw/web/js/app/right-drawer.js` or the current lifecycle owner

Tasks:

1. Register the Feature-specific color control without changing the generic `ColorValueControl`.
2. Replace popup header, Edit-tab fill, Feature-list fill, and selection-toolbar fill with compact/full variants or actions that call the same fill request owner.
3. Bind every Feature presentation and Feature-derived legend color/caption action to `requestFeatureFillChange()`; derive displayed values from the same effective fill/legend projection after commit/cancel.
4. Render scope buttons from the planned candidate array. Show Feature count, affected Result count, active Result name for local targets, and batch-wide progress; never hide Result extent behind `Apply to all (N)`.
5. Bind confirmation/cancel by plan token and candidate ID.
6. Show inherited origin and effective hex color. Keep the picker enabled. Add explicit `No fill` and `Use inherited` actions.
7. Close/reset pending plans on target change, Result switch, generation, Reset, session load, and component disposal.
8. Bind Specific Rules editor, preset, TSV import, and applied palette/default-color controls to the bulk semantic-snapshot request owner; commit on acceptance, not each local input event.
9. Remove template/app exports and undo wrappers for raw color/rule/palette setters, selection bulk fill, Feature-derived legend fill/caption, and old reset branches only after all consumers use the commands.
10. Keep global nullable color controls and their `Auto / None / Color` behavior unchanged.

Phase gate:

- Static architecture contracts reject raw setters in templates/public app exports.
- Real popup/list DOM tests open the scope dialog.
- Cancel is mutation-free.
- Group and single choices reach immediate converged states; a Circular batch group choice updates every affected Result before command completion.
- Result switching and Undo/Redo after a batch edit preserve the selected Result and reveal no stale sibling; live rule/preset/import/palette paths have the same guarantee.
- Generic global `Auto / None / Color` tests still pass unchanged.

### Phase 5: Close adjacent label/stroke bypasses

Primary files:

- `gbdraw/web/js/app/feature-editor/label-actions.js`
- `gbdraw/web/js/app/feature-editor/color-actions.js` or a focused stroke command module
- `gbdraw/web/js/app/feature-editor.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/index.html`

Tasks:

1. Route popup label editing through the existing label scope request and confirmed command; delete the direct popup-only override path.
2. Reuse stable scope candidates for rendered/source label groups while retaining label-specific persistence and reflow behavior.
3. Present inherited stroke through the Feature-specific effective-style control.
4. Convert same-legend stroke and selected-feature stroke edits into explicit-scope History commands with preflight/apply/revert.
5. Keep visibility on its existing command boundary and add it to the structural allowlist rather than rewriting it.
6. Remove every superseded label/stroke public raw mutation export and template call. Fill, selection fill, and Feature-derived legend fill were already closed in Phases 2–4.
7. Add the second-consumer capability matrix here, not in Phase 1: stroke exposes no matching-rule scope without a renderer contract; rendered/source label scopes use stable identities and every affected Result.

Phase gate:

- Fill, stroke, and label public surfaces contain request/resolve operations only.
- Popup label group choice and stroke group choice are real-DOM tested.
- Each accepted operation is one History entry and one affected-Result transaction where applicable.
- Label/stroke capability cases pass only after their property-specific commands exist; Phase 1 remains fill-specific.
- No generic abstraction remains with only one consumer.

### Phase 6: Documentation, visual audit, and final cleanup

Tasks:

1. Add concise user-facing help: inherited color is editable; reset returns to palette/rule; group options state Feature count and Result extent; Circular batch rule edits update every affected output before completion.
2. Update release notes with the restored scope options, unlocked inherited color editing, affected-Result atomicity, and legend consistency fix.
3. Document session v41/request schema 6, confirmed v40 recovery limits, additive catalogue `instanceHash`, the reserved exact TSV pseudo-qualifier, typed-request stability, standalone CLI/download warning, and Generate-to-enable legacy diagnostic.
4. Audit Gallery tutorials/manual screenshots for the old `Auto` control or missing scope dialog. If a visible workflow is affected, use the Gallery screenshot maintenance workflow and regenerate only named captures from declared sessions.
5. Delete orphaned dialog fields, old reset branches, raw setters, duplicate legend reconciliation helpers, unused imports, and tests that directly exercise retired public APIs.
6. Audit production, tests, docs, compatibility fixtures, visual references, and metadata-bearing interactive artifacts as separate diffs.
7. Confirm visual/geometry reference outputs are unchanged. If geometry changes unexpectedly, stop and investigate rather than updating references; review intentional catalogue/interactive metadata diffs separately.

Phase gate:

- No orphaned scope implementation or dual mutation API remains.
- User help, TSV/interactive references, compatibility docs, and release notes match the implemented labels and version gates.
- Any screenshot update has reproducible capture evidence and visual review.
- No unrelated generated artifact changes.

## 8. Test implementation plan

Tests are enforcement for the architecture decisions above. They are not a substitute for private primitives, a single command boundary, derived state, or atomic commit.

### 8.1 Test layers and ownership

| Layer | Target file | Purpose |
| --- | --- | --- |
| Static architecture | `tests/web/architecture-contracts.test.mjs` | Make bypass structurally illegal |
| Pure planner | new `tests/web/feature-fill-scope-plan.test.mjs` | Scope/value/identity matrix without DOM |
| Effective view model | new or planner test file | Inherited effective color and reset semantics |
| Fill command | refactor `tests/web/feature-color-actions.test.mjs` | Semantic + DOM + History atomicity |
| Legend projection | extend `tests/web/legend-sync.test.mjs` | Used groups, manual entries, measurement, rollback |
| Result transaction | extend `preview-runtime.test.mjs` and `svg-result-ingestion.test.mjs` | `a` mounted candidates, `r-a` detached candidates, fresh admission ownership, one array swap |
| Session/request | existing request/authority/export tests | One-fingerprint rule/palette overlay, schema 6, Generate/export concurrency |
| Compatibility | Frozen Web/Python v40 fixtures plus v41/catalogue tests | Confirmed broad-rule recovery, ambiguous rejection, exact-capability promotion |
| Real UI | new `tests/web/feature-style-editor.playwright.spec.js` | Popup/list wiring, scope UI, immediate convergence |
| Mode/regeneration | new spec plus existing Linear spec refactor | Circular single/grid/batch and Linear save-load-regenerate parity |
| Performance | existing Web performance spec | Single-Result fast path and affected-Result-proportional batch path; no Pyodide/remount/input storm |

### 8.2 Static architecture contracts

Add contracts that inspect reachable production sources and `index.html`:

| ID | Contract |
| --- | --- |
| ARCH-STYLE-01 | `index.html` contains no call to `setFeatureColor`, `setFeatureColorValue`, direct rule mutation, or direct popup label/stroke mutation. |
| ARCH-STYLE-02 | `feature-editor.js` and `app-setup.js` do not return/export raw Feature style setters. |
| ARCH-STYLE-03 | Feature fill/stroke controls use the Feature-specific component; generic `ColorValueControl` remains allowed only on non-Feature nullable settings. |
| ARCH-STYLE-04 | Feature-editor writes to `featureColorOverrides` occur only in the derived projection/restore owners, never in templates or presentation actions. |
| ARCH-STYLE-05 | Scope dialog markup iterates planned candidates; it does not enumerate a second hard-coded scope inventory. |
| ARCH-STYLE-06 | `app/file-imports.js` remains the one shared specific-rule parser/serializer owner; session, Generate, and download callers do not add a parallel serializer. |
| ARCH-STYLE-07 | Main-thread Pyodide is not imported/called by the Feature fill command or editor legend projection. |
| ARCH-STYLE-08 | The deep `manualSpecificRules` watcher cannot apply a second fill/legend projection after the command; all live writers appear in a named allowlist. |
| ARCH-STYLE-09 | Selection fill and Feature-derived legend color/caption have no raw mutation path after the fill UI cutover. |
| ARCH-STYLE-10 | Circular batch group candidates carry `all-affected-results` and `targetsByResult`; no active-only fallback or persisted `resultIndex` selector is present. |
| ARCH-STYLE-11 | Specific Rules, preset/import, and applied palette/default-color live paths call the bulk snapshot command; no template/public raw writer remains. |
| ARCH-STYLE-12 | `__gbdraw_instance_hash__` is literal/reserved only under schema 6; biological `instance_hash` and schema-5 qualifier behavior remain intact. |

Use owner allowlists rather than brittle whole-repository occurrence counts where config hydration, History restoration, or import actions are legitimate writers. The allowlist must name owners, not arbitrary paths added to make a test pass.

### 8.3 Pure scope planner and effective-style matrix

Add table-driven Node tests:

| ID | Scenario | Expected |
| --- | --- | --- |
| PLAN-01 | Palette-derived CDS with no local override | effective hex, `origin=palette`, editable picker, `canReset=false` |
| PLAN-02 | Regex-rule-derived CDS | effective rule color/caption, `origin=specific-rule` |
| PLAN-03 | Exact local rule | effective local color, reset reveals next rule/palette |
| PLAN-04 | Same caption has other features | `legend-caption` plus `single` candidates with exact counts |
| PLAN-05 | Matching regex, rendered label, and source label all exist | deterministic complete candidate order and stable target keys |
| PLAN-06 | No group candidate with exact capability | plan is `ready` with `single`; no dialog |
| PLAN-07 | Requested caption conflicts with an existing color | explicit reuse/separate-caption conflict choices |
| PLAN-08 | `none` and `inherit` | values remain distinct; inherit removes owned intent, none persists explicit no-fill |
| PLAN-09 | Duplicate public record IDs and same-record identical Features | record key plus biological collision discriminator selects only the intended Feature |
| PLAN-10 | Visibility before/after changes | planner includes hidden source Features for semantic comparison; only after-state rendered members enter derived legend, and topology change requests preflight rerender |
| PLAN-11 | Stale generation/result key | plan is invalid; no targets |
| PLAN-12 | Deep-frozen inputs | no input mutation; identical input gives identical normalized plan |
| PLAN-13 | Circular batch each group scope | matching-rule, legend-caption, rendered-label, and source-label candidates each carry exact per-Result counts/keys |
| PLAN-14 | Circular batch single/selection | only explicit targets participate; no unrelated Result is inferred from the shared catalogue |
| PLAN-15 | Broad current-Result snapshot request | invalid/unsupported; never mass-materialized through exact or ambiguous hashes |
| PLAN-16 | Web/Python precedence oracle | exact type then wildcard; within each, reserved literal, hash, record-location, source-order qualifier, location |
| PLAN-17 | Legacy exact capability | safe v40 canonical-pair upgrade enables exact scope; failed upgrade disables it with Generate diagnostic |
| PLAN-18 | RuleColor round trip | 3/4/6/8-digit hex and `none` round-trip; picker emits six-digit hex; no-fill derived legend is omitted |
| PLAN-19 | Reserved-key versioning | schema-6 canonical exact row uses literal matching; schema-5 biological `instance_hash` and legacy qualifiers keep regex/source semantics |
| PLAN-20 | Identity stability domain | crop/reverse preserve token, record-key change changes token, typed request with same key replays, legacy CLI hash behavior is unchanged |

Include a cross-product over `color / none / inherit`, `single / group`, and `palette / regex / exact override` for invariants, not a hand-picked happy path only.

### 8.4 Atomic command and History tests

Refactor the current color-action harness so legend actions and PreviewRuntime are realistic staged dependencies rather than always-success mocks.

| ID | Scenario | Required assertions |
| --- | --- | --- |
| CMD-01 | Apply same-caption group color | all target rules/cache/Feature fills update; existing legend swatch updates; one affected-Result transaction |
| CMD-02 | Apply single different color | one exact rule; original group unchanged; deterministic suffixed legend created immediately |
| CMD-03 | Reuse existing caption color | no duplicate caption; feature joins existing group |
| CMD-04 | Explicit `none` | rule/cache/DOM/legend membership follow renderer semantics consistently |
| CMD-05 | Inherit single | only owned exact override removed; next broader rule/palette becomes effective |
| CMD-06 | Inherit group/rule | selected rule intent removed or reset without affecting unrelated precedence |
| CMD-07 | No-op | no rule/cache/DOM/legend/Result/History change |
| CMD-08 | Cancel before command | no command built; zero History entries |
| CMD-09 | Apply then Undo then Redo | every semantic/artifact domain, affected Result, and caption membership round-trips exactly |
| CMD-10 | Stale target before apply/revert | command returns false; History stacks and every Result remain unchanged |
| CMD-11 | Selected-feature scope | one command/History entry/transaction for all selected targets |
| CMD-12 | Picker event storm | many local `input` events plus one accepted `change` produce one plan and one command |
| CMD-13 | Circular batch shared rule | every affected Result and only those Results receive staged Feature/legend content; selected Result stays fixed |
| CMD-14 | Circular batch Undo/Redo | all Result contents, semantic revision, and mounted view restore as one History transition |
| CMD-15 | Out-of-order preflights | older completion cannot commit or apply inverse data over a newer command |
| CMD-16 | History post-operation capture failure | apply/Undo/Redo compensate, artifacts and both stacks equal before-state |
| CMD-17 | Switch Result then Undo/Redo | affected and unaffected current mounts both preserve selection and receive `a`/`r-a` processing |
| CMD-18 | Bulk writer matrix | rule field/add/remove/reorder, preset, TSV import, and palette/default apply each update all affected Results once |
| CMD-19 | Oversized batch | estimate exceeds max bytes before apply; no semantic/DOM/Result/History change |
| CMD-20 | Stack bookkeeping failure | failures after capture, push, redo clear, limit eviction, or file retention compensate intent/stacks/bytes/files/revision exactly |

Failure-injection matrix:

- missing active SVG;
- missing target Feature node;
- missing or malformed target in one affected inactive Result;
- detached parse/serialization failure in the first, middle, and final affected Result;
- missing legend reflow metadata;
- browser text measurement failure;
- font readiness/measurement failure for one inactive Result;
- staged legend conflict;
- stale plan token after preflight;
- injected exception after semantic replacement;
- injected exception after Feature patch;
- injected exception after legend patch;
- non-mutating mounted serialization, Result-array assignment, or admission-metadata finalization failure;
- queued legacy deep-rule/palette watcher (architecture test proves no projection callback remains);
- post-apply/revert/redo intent-capture and each stack/file-bookkeeping failure;
- required-geometry worker failure and History-cap preflight rejection.

For each seam, snapshot normalized rules/palette, override cache, active `legendEntries`, global manual legend intent, relevant DOM, every Result object/content/admission record, Result-validation ledger, selected Result, PreviewRuntime mount/dirty/index state, document epoch, semantic revision/fingerprint, current History intent, retained files/bytes, and both stacks. Assert exact equality after failure. The transaction does not call the state-writing flush API; any exception after mounted patching must use compensation before reporting failure.

### 8.5 Legend projection tests

Extend legend tests with:

| ID | Scenario | Expected |
| --- | --- | --- |
| LEGEND-01 | Existing Feature group color update | update only swatch/rule projection; no text measurement |
| LEGEND-02 | Single-feature split | new caption/color/featureIds inserted and layout reflowed in same prepared patch |
| LEGEND-03 | Feature rejoins existing group | old unused derived entry removed; existing entry gains Feature ID |
| LEGEND-04 | Unused specific rule | no Feature-derived legend entry unless actually used by committed features |
| LEGEND-05 | Document-global manual legend-only entry | projected into each Result and preserved while each Result derives only its own Feature groups/layout |
| LEGEND-06 | Caption/color conflict | rejected or resolved by the explicit plan; never two colors under one caption |
| LEGEND-07 | Undo/Redo | order, transform, stroke, feature IDs, and colors restore exactly |
| LEGEND-08 | Measurement path | same-font browser measurement used; main-thread Pyodide call count remains zero |
| LEGEND-09 | Missing metadata/measurement | preflight fails before active mutation |
| LEGEND-10 | Circular batch group edit | each affected Result derives only its used members; inactive legend SVGs are staged before commit |
| LEGEND-11 | Browser font readiness | connected measurement host waits for the intended font and produces deterministic per-Result layout |
| LEGEND-12 | Legacy batch manual-state conflict | v40 inactive artifact conflict is diagnosed; no artifact is promoted or guessed before regeneration |

Do not assert only `legendEntries`; compare the reactive model, live `g[data-legend-key]` nodes, swatch fills, Feature IDs, order/transforms, and serialized Result.

### 8.6 Canonical request and session tests

Add or extend the following:

| ID | Target | Scenario |
| --- | --- | --- |
| SESSION-01 | `session-request.test.mjs` | normalized rules and applied palette serialize once into their canonical color resources under schema 6 |
| SESSION-02 | `session-draft-authority.test.mjs` | style revision overlays only changed color resources; unrelated drafts remain unchanged |
| SESSION-03 | `session-export-validation.test.mjs` | cached committed request cannot bypass a newer style revision |
| SESSION-04 | `session-authority.test.mjs` | rule authority/cache/artifact inventory and mismatch rejection |
| SESSION-05 | v41 config round trip | direct edit -> writer -> reader preserves rules/palette, global manual legend intent, exact identity, and derived cache |
| SESSION-06 | frozen main-backed v40 fixture | broad+exact rule re-resolution repairs the missing canonical table and yields expected effective request |
| SESSION-07 | ambiguous v40 mutation | recovery refuses to guess and produces a non-mutating diagnostic |
| SESSION-08 | resource identity | sequences, comparisons, LOSAT evidence, and unchanged resources retain identity/bytes |
| SESSION-09 | failed/cancelled generation | semantic revision and last valid canonical basis remain correct |
| SESSION-10 | Python/Web parity | full-envelope readers produce the same recovered v40 request and v41/schema-6 projection |
| SESSION-11 | edit/Undo during successful Generate | candidate Results and both color resources rebase only with proven topology; otherwise generation restarts |
| SESSION-12 | edit during asynchronous export | stale envelope is discarded/restarted or fails before download; no mixed revision is emitted |
| SESSION-13 | Circular batch writer/reader | every saved Result and shared rule/palette resources represent one fingerprint; switch/load/regenerate preserves all outputs |
| SESSION-14 | v40 exact capability | canonical-pair upgrade enables exact scopes; missing/ambiguous pair disables them and v41 save requests Generate |
| SESSION-15 | version gates | v41/schema 6 rejects silently incompatible old consumers; schema-3 catalogue and old standalone assets remain readable |

SESSION-06/07/10 use the minimized main-backed producer with commit/build/hash and expected output; the negative fixture is derived by one named ambiguity mutation. Keep the original Gallery session as provenance, not as the only test input.

### 8.7 Real browser tests through visible controls

Create `tests/web/feature-style-editor.playwright.spec.js`. Use role/label/test-ID selectors and actual picker/select/button DOM events. Do not call scope or mutation methods through `page.evaluate()` to perform the edit.

| ID | Mode | Workflow | Assertions after the command promise/status settles |
| --- | --- | --- | --- |
| UI-01 | Linear lambda | Open untouched CDS | shows effective hex and inherited origin; picker enabled; no disabled `Auto` mode |
| UI-02 | Linear lambda | Change fill, then Cancel | complete scope list visible; no Feature/rule/legend/Result/History change |
| UI-03 | Linear lambda | Choose same legend group | all group Feature fills and existing legend update; no suffixed entry; popup/model/Result agree |
| UI-04 | Linear lambda | Choose single feature | only target changes; original legend remains; new suffixed entry appears before completion |
| UI-05 | Linear lambda | Use inherited and `No fill` | explicit actions use scope planner and converge atomically |
| UI-06 | Linear lambda | Undo/Redo UI-03 and UI-04 | visual and session-owned state round-trip; one History item per edit |
| UI-07 | Feature list | Edit the same kind of feature | same scope dialog and command behavior as popup |
| UI-08 | Popup label | Edit shared rendered/source label | scope dialog opens through actual input/button; selected group changes once |
| UI-09 | Popup stroke | Edit inherited stroke and apply same legend | picker enabled; explicit group command; one Result transaction/History entry |
| UI-10 | Circular fixture | Group and single fill | same semantic/legend behavior on Circular SVG structure |
| UI-11 | Lifecycle | Open scope dialog then switch Result/load session/Reset | pending plan invalidated; confirmation cannot mutate new Result |
| UI-12 | Runtime | Perform first editor color change before main Pyodide ready | command succeeds; Pyodide initialization counter stays zero |
| UI-13 | Legend editor | Change a Feature-derived entry, then a manual-only entry | Feature entry uses semantic command; manual intent projects globally with Result-local layout |
| UI-14 | Rich/simple popup | Repeat inherited display and one group edit with Rich Feature Popup on and off | both modes use the same editor control and command path |
| UI-15 | Duplicate presentation | Edit through header, then Edit tab, then Feature list | every instance opens the same scopes and all instances resynchronize after commit/cancel |
| UI-16 | Circular batch | Choose matching-rule and one non-rule group spanning outputs | buttons name Feature/Result totals; every affected output updates; one History item each |
| UI-17 | Circular batch lifecycle | Cancel progress or inject one inactive-Result failure | no Result/rule/History mutation; later confirmation from the stale token cannot commit |
| UI-18 | Bulk writers | Accept one rule edit, preset/import, and palette change | each shows progress, updates all affected outputs atomically, and is one History item |
| UI-19 | Legacy exact capability | load safe and unsafe v40 catalogues | safe upgrade permits single; unsafe case disables exact controls with Generate diagnostic |

"Immediate" is a state-boundary assertion, not a generous sleep:

1. Wait for the exposed command status/token to become `committed` or for the dialog to close after the actual click.
2. In the next synchronous evaluation, assert Feature DOM, legend DOM/view, normalized rules/palette, popup view, mounted-serialization/array-swap counters, admission ownership, and every affected Result show the same fingerprint. Assert unaffected Result object identity is unchanged.
3. Do not wait for Pyodide readiness, a later regeneration, or a fixed multi-second timeout to make the assertion pass.

### 8.8 Save-load-regenerate mode matrix

At least these end-to-end cases are required:

| ID | Mode/topology | Edit | Round trip |
| --- | --- | --- | --- |
| RT-01 | Linear multi-record | rendered-label or caption group | edit -> save -> new context load -> regenerate |
| RT-02 | Linear multi-record | single split caption | edit -> save -> load -> regenerate |
| RT-03 | Circular single | matching rule group | edit -> save -> load -> regenerate |
| RT-04 | Circular grid when shared catalogue/rules apply | one semantic group | edit -> save -> load -> regenerate |
| RT-05 | Circular batch, at least three outputs | matching-rule and rendered/source/legend group cases | edit -> switch every Result -> save -> new context load -> regenerate every output |
| RT-06 | Circular batch with duplicate public IDs and same-record duplicate Features | single exact edit | only intended instance changes; save/load/regenerate never leaks |
| RT-07 | Either mode | inherit after explicit edit | edit -> save/load/regenerate -> effective fallback parity |
| RT-08 | Either mode | `none` | edit -> save/load/regenerate -> no-fill and legend parity |
| RT-09 | Circular batch | rule import/preset and applied palette | each bulk edit -> switch -> save/load/regenerate; all affected outputs and resources agree |

For every case compare normalized specific rules, target Feature fills, non-target fills, caption/color/feature membership, legend order, and Result count/topology. Do not accept “a magenta legend exists” without proving it represents the intended targets and that no stale old group remains.

### 8.9 Performance and ownership guards

Extend the existing Playwright performance lane with a representative large session and test hooks/counters already used by PreviewRuntime/ingestion tests.

| ID | Guard |
| --- | --- |
| PERF-STYLE-01 | Ten picker-local `input` events cause zero semantic scans, History items, serializations, and root mutations before acceptance. |
| PERF-STYLE-02 | One accepted single edit causes one planner pass, one command, one mounted serialization and array swap, zero active parse/sanitize/remount/index rebuild/handler work, and zero Pyodide init. |
| PERF-STYLE-03 | One caption-group edit records rule evaluations, performs one reconciliation per affected Result, and does not rebuild the mounted Feature index. |
| PERF-STYLE-04 | Repeating no-op color assignments performs no Result update or History work. |
| PERF-STYLE-05 | Undo/Redo performs one command projection each and does not enter admission/sanitization/remount; only transaction-owned admitted-content metadata is replaced. |
| PERF-STYLE-06 | A batch apply/revert performs `a` mounted serializations, `r-a` detached passes, and one shallow Result-array swap; unaffected Results have zero parse/serialize and retain identity. |
| PERF-STYLE-07 | Updating a shared matching rule does not generate one exact TSV rule per matched Feature; rule growth follows durable intent, not target count. |
| PERF-STYLE-08 | `estimateBytes()` accounts for Result/admission/semantic retention; a command above max bytes is rejected before apply, while a fitting command may evict only older entries. |
| PERF-STYLE-09 | Representative `n`, `m`, `R`, `r`, and `S_all` counters support the declared `O(R + n*m + k + S_all)` upper bound; exact lookups use the literal map. |

Use structural counters as the primary gate. Record elapsed time for diagnosis, but do not replace ownership assertions with a flaky wall-clock-only threshold.

### 8.10 Focused and broad verification commands

Run focused tests after each phase, then the justified broad gate. Exact filenames may be added by the implementation, but the intended command set is:

```bash
node --test tests/web/feature-fill-scope-plan.test.mjs
node --test tests/web/feature-color-actions.test.mjs
node --test tests/web/feature-catalog.test.mjs
node --test tests/web/file-imports.test.mjs
node --test tests/web/legend-sync.test.mjs
node --test tests/web/preview-runtime.test.mjs
node --test tests/web/svg-result-ingestion.test.mjs
node --test tests/web/history.test.mjs
node --test tests/web/session-request.test.mjs
node --test tests/web/session-draft-authority.test.mjs
node --test tests/web/session-export-validation.test.mjs
node --test tests/web/session-authority.test.mjs
node --test tests/web/architecture-contracts.test.mjs

npx playwright test tests/web/feature-style-editor.playwright.spec.js \
  --project=chromium --workers=1
npx playwright test tests/web/linear-multi-record.playwright.spec.js \
  --project=chromium --workers=1

pytest tests/test_feature_selector_values.py tests/test_color_table_parsing.py \
  tests/test_feature_catalog_source_identity.py tests/test_web_feature_metadata.py \
  tests/test_web_feature_catalog.py -v
pytest tests/test_api_request_render.py tests/test_session_request_codec.py \
  tests/test_interactive_svg_cli_format.py -v
pytest tests/test_session_compat.py tests/test_api_session.py tests/test_session_io.py -v

node --test tests/web/*.test.mjs
npm run test:web:functional-smoke
npm run test:web:perf-smoke
pytest tests/ -v -m "not slow"
pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
git diff --check
```

Prepare the ignored browser wheel when browser/session tests require it. Do not refresh the cache-bust token unless preparing a deployable bundle. Allow long-running pytest commands at least 30 minutes and report progress at intervals shorter than 60 seconds.

## 9. Implementation order and change slicing

Execute the numbered phases in order, but use four internally consistent landing slices. Phases 2–4 form one non-separable fill vertical slice: they may be reviewed and tested as separate diff sections, but must not land or be marked complete independently. Exposing the new command before canonical revision safety, or exposing the UI before every fill writer is cut over, would recreate a stale-session or dual-owner interval.

1. **Exact identity, pure planning, and effective-style model (Phase 1)**
   - pre-render shared identity plan, schema-6 reserved selector, v41/catalogue compatibility, Web/Python precedence oracle, planner/view-model modules, and unit tests;
   - no editor wiring or visual geometry/fill change when the selector is absent; catalogue/interactive metadata diffs are intentional and reviewed separately.
2. **Complete fill vertical slice (Phases 2–4; one landing slice)**
   - all-affected-Result staging, global-manual/per-Result legend projection, watcher removal, admission-safe Result swap, History compensation, and failure/performance evidence;
   - rule/palette-revision-safe Generate/session writer and mandatory v40 recovery;
   - popup/list/selection/Feature-derived-legend plus rule/preset/import/palette UI cutover, then raw writer deletion;
   - the slice is acceptable only when command, canonical/session, architecture, and real-browser gates all pass together.
3. **Label/stroke convergence (Phase 5)**
   - close the adjacent bypasses using the proven protocol; remove only the superseded label/stroke paths.
4. **Documentation and cleanup (Phase 6)**
   - user wording, v40/v41/schema-6 compatibility documentation, visual audit, orphan removal, and separate artifact review.

If implemented as one pull request, preserve these as clearly separated diff sections or commits only where the boundary above permits it. Do not keep an old and new mutation owner active as a fallback, and do not publish an intermediate build between Phases 2 and 4.

## 10. Risk register and mitigations

| Risk | Failure mode | Mitigation/evidence |
| --- | --- | --- |
| Rule precedence changes | A reset exposes the wrong color or removes a broader rule | Pure precedence/inherit matrix; normalized before/after rule snapshots |
| Stable identity collision | Wrong duplicate record or same-record Feature instance changes | Shared pre-render plan, reserved literal token, source-index collision discriminator, and duplicate save/regenerate tests; never persist `resultIndex` |
| Active-only scope narrowing | Shared rule updates mounted output but leaves sibling Results stale | `resultExtent`, `targetsByResult`, Circular batch transaction and switch/save/regenerate gates |
| Exact-rule explosion | Broad semantic edit becomes thousands of snapshot rules | Shared-rule intent remains shared; planner/performance guard rejects target-count-proportional rule growth |
| Legend owner overreach | Manual legend-only entries disappear or inactive artifact becomes authority | Document-global manual intent, Result-local derived layout, v40 conflict diagnostic, and preservation tests |
| Async stale commit | A dialog edits a different Result after load/switch | Result/generation/token binding and lifecycle invalidation tests |
| Partial multi-Result rollback | One Result, mounted DOM, admission record, or semantic state commits alone | Stage all Results first; one array swap; failure injection and exact object/runtime/DOM compensation |
| Pyodide removal changes text width | Legend layout clips or drifts | Same-font SVG measurement tests and readable-scale visual inspection |
| Large batch regression | Too many Results/rules are evaluated or retained | `R/n/m/r/S_all` counters, `a`/`r-a` projection bounds, size pre-rejection, progress/cancel |
| History post-capture/bookkeeping failure | Mutation succeeds without a coherent stack/intent/file state | Compensation-mode tests at every post-operation seam and epoch pruning |
| Generate/export revision race | Newer rule/palette artifacts pair with older canonical resources | Immutable snapshots, topology coverage, adopt/download recheck, controlled interleavings |
| v40 over-migration | Reader invents broad-rule intent or exact identity | Main-backed fixture, full resolver agreement, canonical-pair-only identity upgrade, ambiguous rejection |
| Version/TSV divergence | Old consumer treats exact pseudo-key as a biological qualifier | v41/schema-6 gate, reserved literal parser, schema-5 qualifier regression, CLI portability diagnostic |
| Cross-surface divergence | Visibility/fill/legend/catalogue or Web/Python differ | One pre-render identity plan and complete precedence/fixture parity |
| Generic color regression | Other nullable settings lose Auto/None | Preserve component; rerun existing global color tests unchanged |
| Shared dirty files | Unrelated user changes are overwritten | Phase 0 diff audit and separate final diff review |

## 11. Completion criteria

Implementation is complete only when all statements are true:

### Product behavior

- Untouched Feature fill and stroke show their effective value and are directly editable.
- Popup and Feature list expose all applicable semantic scopes with target counts.
- Circular batch candidates expose affected Result counts; shared semantic scopes update every affected output and never silently become current-output snapshots.
- Cancel changes nothing.
- Group edits update the intended features and one existing legend group atomically.
- Single edits preserve the old group and add/reuse a distinct legend group atomically.
- Single edits remain exact across duplicate public record IDs and same-record duplicate Features by using validated stable record-instance identity; unsafe legacy catalogues fail closed with a Generate diagnostic.
- Specific Rules, preset/import, and applied palette/default-color live edits update all affected Results atomically.
- `No fill` and `Use inherited` use the same plan/command boundary.
- Popup label text no longer bypasses its scope flow.

### Architecture

- Templates and public app exports cannot call raw style setters.
- Scope candidates come from one pure planner.
- Semantic scope and Result extent are separate values; flattened catalogues retain per-Result ownership.
- `manualSpecificRules` plus applied palette/default state own renderer-visible fill intent; global manual legend intent is separate from Result-derived Feature legends.
- Feature override cache, Feature-derived legend entries, popup state, DOM, and Result content are projections.
- No active mutation occurs before async preflight succeeds.
- One confirmed edit is one atomic History command and one affected-Result transaction with zero/one mounted serialization and one Result-array assignment.
- History capture/bookkeeping failure restores exact revision, Results/admission, intent, stacks, bytes, and retained files; oversized commands fail before apply.
- No deep rule/palette watcher can replay a second projection.
- Editor color/legend work does not initialize main-thread Pyodide.

### Persistence and convergence

- Effective canonical rules and default colors cannot be older than the saved editor style fingerprint.
- Generate adoption and session download cannot cross rule/palette revisions or rebase across a Feature-topology change.
- The v41/schema-6 writer and reader are deterministic; the confirmed v40 broad-rule recovery agrees in Web/Python, and unsafe exact-identity upgrade fails closed.
- Linear and Circular single/grid/batch save-load-regenerate parity passes.
- Undo/Redo and candidate replay preserve the same state.

### Verification and artifacts

- Static, selector/precedence, planner, command, legend, session, applicable compatibility, real-UI, mode, concurrency, failure, and performance tests pass.
- Existing global `Auto / None / Color` behavior passes unchanged.
- Fast Python tests, Web functional/performance suites, output comparison, Ruff, and `git diff --check` pass.
- Visual/geometry references remain unchanged unless separately reviewed; intentional catalogue/interactive metadata-byte diffs are reviewed separately.
- Production, tests, documentation, fixtures, and generated-artifact diffs are reviewed separately.
- Any Gallery screenshot change has regeneration and visual-inspection evidence.

## 12. Execution evidence ledger template

Keep this section updated during implementation rather than marking phases complete from code inspection alone.

| Work item | Status | Production evidence | Test/command evidence | Browser/visual evidence | Artifact impact | Remaining issue |
| --- | --- | --- | --- | --- | --- | --- |
| Phase 0 baseline/reproduction | complete | HEAD/tree/worktree and all style mutation owners were audited; the disabled inherited picker, raw setter bypass, active-only projection, History compensation gap, and Generate/export revision seams were reproduced | Initial focused Python and Node gates established the pre-change behavior | Real DOM probes reproduced the inherited-picker lock and raw fill path | none | none |
| Phase 0 confirmed v40 fixtures | complete | The named Gallery SHA-256 and both canonical color resources were inspected; the real stale authority is the default-color resource, not the specific table | Minimized positive, consistent, and ambiguous fixtures exercise Web/Python recovery without mutating their inputs | All six tracked v40 Gallery sessions were probed; tobacco empty-rule authority remains compatible | compact compatibility fixtures added | Final broad rerun still required |
| Phase 1 identity/schema/scope/effective-style model | implemented; focused verified | Pre-render `FeatureInstanceIdentityPlan`, `fi1_` literal selector, schema 6/v41 gates, semantic `fs1:` selectors, effective fill/stroke view models, pure scope planners, and Result ownership are wired in Python and Web | Identity, precedence, catalogue, parser, planner, request, and session focused suites passed; shared Python/Web selector vectors pass | A real majanivirus click exposes Feature type, source-label, Similarity-group, legend, and single choices with target counts | intentional catalogue/interactive metadata additions | The new end-to-end Similarity acceptance currently fails after scope commit; see Section 13 |
| Phase 2 affected-Result fill/bulk-writer/legend/History command | implemented; focused verified | Atomic fill, bulk-style, Result/legend projection, manual legend, legend order, and compensated History commands are present; deep fill/rule replay owners were removed | Command, rollback, inactive-Result, History compensation, ledger, legend projection, and architecture suites passed | Circular batch and immediate Feature-derived legend paths passed before the latest Similarity acceptance was added | no geometry-reference change expected | Similarity semantic intent is committed but its current Results/legend are not recolored in the real Gallery case |
| Phase 3 v40/v41/schema-6 canonical/session safety | implemented; focused verified | v41/schema-6 writers/readers, copy-on-write color authority, deterministic v40 recovery, style revision/fingerprint/Result ledger, and Generate/export rechecks are present | Focused Web authority/recovery/export and Python session/request/compatibility suites passed | Source-label save-load-regenerate passed in real Linear Chromium | compact fixtures and current-writer metadata changes | Final full Python/functional rerun and current-writer artifact publication remain |
| Phase 4 Feature/rule/palette/derived-legend UI cutover | blocked in final acceptance | Popup, Feature list, selection, Specific Rules, TSV, preset, palette/default, non-Feature legend, and derived legend paths are routed through request/command facades; managed exact/semantic rows are opaque | All 86 Web Node suites passed at an intermediate revision; focused style suites pass at the handoff revision | The seven-case Feature style suite passed before the new full Similarity apply/save/load/regenerate assertion; that assertion now fails with zero recolored Result/legend targets | docs capture audit found no required screenshot replacement | Fix the real Similarity Result projection and the manual-legend stale failure, then rerun the full suite |
| Phase 5 label/stroke convergence | implemented; focused verified | Fill-equivalent label and stroke planners/actions/atomic Result commands are wired; Label TSV rollback was added | Focused label, Label TSV, stroke, History, and adjacent fill suites passed | Visible label/stroke editor regressions passed in the Feature style browser suite | none expected | Broad functional rerun remains; Label TSV remains current-Result scoped by its existing product contract |
| Phase 6 documentation/visual cleanup | in progress | Reference, Web app, compatibility, interactive metadata, release-note, capture, and recipe contracts were updated | Documentation/recipe focused tests and semantic artifact comparators passed | Existing captures were audited; no visual recapture was required for the unchanged UI geometry | five current-writer recipe groups remain stale; see Section 13 | Publish and review current-writer session/interactive metadata artifacts after production gates are green |
| Focused verification | mostly complete | n/a | Latest recorded evidence includes 86/86 Web Node suites at an intermediate revision, 351 focused Python semantic/session tests, 16/16 output comparisons, and focused browser gates | Feature style 7/7 before the extended Similarity assertion; Linear source-label save-load-regenerate 1/1; right drawer/non-Feature palette 4/4 | reference SVG geometry remained unchanged | Rerun all focused gates after the Similarity/manual-legend fixes |
| Broad verification | blocked | n/a | Latest functional smoke was 86 passed / 7 failed; two v40 expected-value assertions were patched to v41 but not rerun, two manual-legend stale failures and three parallel BGC timeouts remain unclassified. Full fast Python, performance smoke, Ruff, and final diff-check have not passed on the handoff tree | New Similarity acceptance fails with a committed rule but zero magenta Features/legend entries | reference outputs remain read-only | Execute the ordered gate list in Section 13 |
| Final diff/visual review | pending | Production/test/docs/fixture/generated-asset separation has not received a final root review | No final all-green evidence bundle yet | No final visual sign-off after the extended Similarity test | unrelated untracked report must be preserved | Do not commit until all blockers and broad gates are closed |

## 13. Session handoff — 2026-08-13

This execution session is intentionally closed before the plan is complete. The working tree contains the implementation, tests, documentation, and fixtures described above, but it has not been committed, pushed, tagged, or published. Do not reset or replace the dirty tree when starting the next session. Preserve the unrelated untracked file `reports/codex-session-knowledge-audit-2026-08-12.md`. The browser wheel was rebuilt locally with `python tools/prepare_browser_wheel.py --no-build-isolation`; it is an ignored generated test artifact and must not be committed. The temporary HTTP server was stopped at handoff.

### 13.1 What is working now

- Clicking a real Feature opens one planner-owned scope inventory. Applicable Feature type, legend caption, rendered label, source annotation label, Similarity group, selected-feature, matching-rule, and single-feature choices are derived from the catalogue rather than hard-coded in the template.
- A single edit uses validated `instanceHash`; broad scopes persist as bounded schema-6 semantic rules rather than expanding into one exact rule per current target.
- Inherited fill/stroke controls show their effective value and remain editable. `No fill` and `Use inherited` use the same request/command boundary.
- Fill, rule/preset/TSV, palette/default-color, label, stroke, manual legend, and legend-order implementations have atomic command and rollback coverage, including inactive Results and History post-capture failures.
- Feature-derived legend projection is Result-local, updates without main-thread Pyodide, and preserves per-Result order/style in the covered command tests.
- v41/schema-6 persistence, v40 color recovery, exact-identity fail-closed behavior, Generate/export style fingerprints, and source-label save-load-regenerate have focused passing evidence.
- Existing reference SVG geometry still matches in the recorded output-comparison run.

### 13.2 Confirmed unfinished work and current failures

1. **Similarity group is the release-blocking product regression.** The new real-browser acceptance at `tests/web/feature-style-editor.playwright.spec.js:104` selects `similarity-group`, commits exactly one rule `* / __gbdraw_semantic_scope__ / fs1:similarity-group:og_1 / #ff00ff`, and creates no exact rules, but the current saved Results contain zero magenta Feature targets and zero magenta legend entries. The captured failure is `test-results/feature-style-editor.playw-677d2-milarity-and-single-choices-chromium/error-context.md`. This must be fixed in the semantic target-to-Result projection boundary; weakening the assertion or restoring post-render fill replay is not acceptable.
2. **Manual legend composition has a stale-intent failure.** In the broad functional smoke, the Circular and Linear `composition-layout-real.playwright.spec.js` cases fail in `addAndRenameLegendEntry()` with `manual-intent-command.js` reporting `intent became stale before command preparation`. Reproduce serially, then fix the legitimate revision/fingerprint contract or caller snapshot timing. Do not bypass the stale guard globally.
3. **The broad browser run is not green.** The last functional smoke reported 86/93 passing. Two assertions expecting session v40 were changed to v41 at `tests/web/depth-track-session.playwright.spec.js:1425` and `:1880` immediately before handoff, but were not rerun. Three BGC session tests timed out under the 16-worker run and must be rerun with one worker before deciding whether they are functional failures or parallel-startup pressure.
4. **Broad non-browser verification is incomplete.** The latest full Python run predates the final semantic-selector fixes and reported 2,909 passed / 63 failed. Many failures were subsequently classified and fixed, but the full fast suite has not been rerun. The performance smoke, repository-wide Ruff, final output comparison, and final `git diff --check` also remain.
5. **Current-writer documentation artifacts are intentionally not yet published.** The stale groups are H-CLI-12, H-PY-05, H-CLI-13, T-PY-08, and T-CLI-11. Their visual SVG geometry is unchanged; the intended changes are v41/schema-6 session data and interactive catalogue metadata. Regenerate only after the production gates are green, then review the generated diff separately.
6. **The final root audit was interrupted.** Raw-writer ownership, final session convergence, production/test/docs/fixture separation, and plan completion criteria require one final review. The implementation must not be described as complete until that review and all required gates pass.

### 13.3 Ordered next-session work

1. Start with `git status --short`; preserve all shared edits and the unrelated report. Read this handoff and the captured Similarity failure before editing.
2. Fix the Similarity target projection so the confirmed command updates every affected Result and its Feature-derived legend immediately. Keep Python-rendered Generate output authoritative; do not reintroduce the obsolete candidate fill overwrite. Run the single real majanivirus test through apply, save, fresh-page load, and Generate until it passes.
3. Reproduce and fix the manual-legend stale-intent composition cases with one worker. Add a focused regression at the boundary that was wrong.
4. Rerun the two v41 expectation tests and the three BGC tests serially. Retain a timeout increase only if measured startup evidence justifies it; do not mask a load failure.
5. Rebuild the ignored browser wheel after any Python renderer change, without refreshing the cache-bust token, and rerun the focused Feature style, Linear, Circular batch, right-drawer, non-Feature palette, and composition browser tests.
6. Run all Web Node tests, then functional and performance smoke. Run the focused Python identity/selector/session suites, the fast Python suite, output comparison, Ruff, architecture guards, and `git diff --check`.
7. When production is green, publish the five current-writer recipe groups with their official scenario commands, inspect metadata/session diffs separately from visual outputs, and rerun recipe/documentation contracts.
8. Perform the final owner/dead-API/diff audit. Update this ledger with exact commands and counts. Only then change the plan status to complete and prepare the single-session commit handoff.

### 13.4 Commands to resume with

```bash
# Release-blocking Similarity acceptance
TMPDIR=/tmp TEMP=/tmp TMP=/tmp npx playwright test \
  tests/web/feature-style-editor.playwright.spec.js \
  --project=chromium --workers=1 \
  --grep "real majanivirus Feature click exposes applicable type"

# Manual legend composition regression
TMPDIR=/tmp TEMP=/tmp TMP=/tmp npx playwright test \
  tests/web/composition-layout-real.playwright.spec.js \
  --project=chromium --workers=1 \
  --grep "real render preserves composition edits"

# Session version and BGC serial classification
TMPDIR=/tmp TEMP=/tmp TMP=/tmp npx playwright test \
  tests/web/depth-track-session.playwright.spec.js \
  --project=chromium --workers=1 \
  --grep "v40 layout preferences|P3 Custom Track drafts|BGC session"

# Focused non-browser gates after the fixes
node --test tests/web/*.test.mjs
python -m pytest tests/test_feature_instance_identity.py \
  tests/test_feature_semantic_selectors.py \
  tests/test_session_request_codec.py tests/test_session_compat.py \
  tests/test_api_session.py -q
python -m pytest tests/test_output_comparison.py::TestOutputComparison -q
ruff check gbdraw/
git diff --check
```

After those pass, run the repository functional/performance and fast Python gates from `AGENTS.md`; do not treat the intermediate counts above as final evidence.

### 13.5 Deferred artifact publication commands

```bash
python docs/recipes/run_cli_scenarios.py --scenario H-CLI-12 --output-root docs/images
python docs/recipes/run_python_scenarios.py --scenario H-PY-05 --output-root docs/images
python docs/recipes/run_cli_scenarios.py --scenario H-CLI-13 --output-root docs/images
python docs/recipes/run_python_scenarios.py --scenario T-PY-08 --output-root docs/images
python docs/recipes/run_cli_scenarios.py --scenario T-CLI-11 --output-root docs/images
```
