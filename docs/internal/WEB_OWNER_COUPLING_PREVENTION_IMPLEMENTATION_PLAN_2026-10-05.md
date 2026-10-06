# Web owner-coupling prevention implementation plan

Status: proposed, 2026-10-05. Baseline: `dev` at `8f450194` (merge of #812).
Scope: `gbdraw/web/js/**`, `gbdraw/web/CLAUDE.md`, `tools/check-web-change-budget.mjs`
and its registries, `tests/web/**`, `.github/workflows/**`.

This plan follows the
[architecture fitness-function ratchet](ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
protocol (detector-first onboarding, trusted-base evaluation, shrink-only
frozen baselines, authority-only PRs). It adds no new protocol; it adds rules,
detectors, and registries the existing protocol can carry.

## 1. Objective

Stop cross-owner coupling in the Web app from growing through channels the
current gate does not measure, and make every existing instance of it a
shrink-only accepted violation that a runtime PR can only remove.

"Cross-owner coupling" here means one owner module reaching another owner
through a path that is not a static import of a service: a whole owner object
injected at composition, a closure that resolves to an owner created later, a
function assigned into `state`, a second copy of a projection, or an API
parameter that carries another owner's concern.

## 2. Evidence

Audit of the 159 first-parent merges since 2026-09-01 that touched
`gbdraw/web/js` (measured per merge; the tool is reproduced in PR B1):

| Metric | 2026-09-01 | `8f450194` | Direction |
| --- | ---: | ---: | --- |
| Static import cycles (`gbdraw/web/js`) | 0 | 0 | flat; the only hard structural check, never broken in 315 PRs |
| State fields written by ≥2 non-bulk owner modules (pairs) | 284 | 233 | down (#665, #670, #728, #741, #757) |
| Repairing watchers | 48 | 37 | down |
| Generate-path draft assignments (R7 baseline) | 62 | 39 | down (#666) |
| Owner-object parameters taken by `create*` factories | 32 | 58 | **up** |
| Forward-reference closures to a later-created owner | 0 | 4 | **up** (#617, #622, #737, #764) |
| Functions assigned into `state.<name>` outside `state.js` | 0 | 3 | **up** (#474, #641, #807) |
| `rulePreparation` trigger sites outside its owner | 0 | 15 (9 modules, 6 call shapes) | **up** (#538 → #795) |
| Projection call shapes, feature visibility + label | 2 | 10 sites in 3 orchestration shapes | **up** (#692, #789, #807) |
| Projection call shapes, legend order | 1 | 4 | **up** (#670, #673, #739, #742) |
| Projection call shapes, palette + specific rules | 1 | 3 | **up** (#538, #617, #670) |

Coupling moved from the measured channels (state co-writes, watchers, import
graph) into the unmeasured ones (runtime injection, lazy closures, state
backdoors, duplicated projections, concern-carrying parameters). Pairs found
beyond feature/label: History ⇄ Legend/LegendLayout (#622, #737, #742),
Legend ⇄ Feature rules (#617, #622, #692), `rulePreparation` fan-out
(#538–#795), Generate → editor state (#670, #761, #784), TrackSlots ⇄ Placement
(#768, #805), run-analysis ⇄ record-display (#510, #574), `state.js` ⇄
app-setup availability cycle (#641).

## 3. Root causes

1. **Unmeasured channel.** `tools/web-architecture-rules.json` carries two rule
   kinds (`single-canonical-entry-edge`, `single-semantic-owner`), both over
   static imports and definitions. No detector observes composition
   (`app/app-setup.js:1190-1191,1298`, `app/feature-editor.js:31`), state
   function assignment (`app/app-setup.js:1226,2506,2507`), or call shapes.
2. **Self-authorization.** `gbdraw/web/CLAUDE.md` and every guard test the
   R1–R12 sections name (34 paths, for example
   `tests/web/track-layout-transition.test.mjs`,
   `tests/web/option-input-integrity.test.mjs`) are absent from `guardPaths`,
   `authorityPaths`, and `governancePaths` in `tools/check-web-change-budget.mjs`
   (lines 380-452). PR #805 rewrote R10 (`gbdraw/web/CLAUDE.md:355-363`),
   added its own allowlist entries
   (`tests/web/track-layout-transition.test.mjs:187-195`), and changed the
   runtime in one STANDARD PR with Gate PASS, Review CLEAR.
3. **Undefined layering.** R10 says a transition "calls the owner's reconcile
   explicitly" but not through what. The cheapest reading is to inject the whole
   owner object and call it, which produced the bidirectional pairs above.
   Creation order in the composition roots then forced lazy closures.
4. **Decisions assign behavior, not the reaction owner.** F-3 ("a hidden
   feature hides its label"), Owner Q2 ("Show feature and label"), Q3 (slot
   inputs through placement) each require owner A to react to owner B. The
   records do not say who owns the reaction or through which channel, so each
   implementation chose locally and differently.

## 4. Definitions

- **Owner module**: a module listed in the Module ownership table of
  `gbdraw/web/CLAUDE.md`, or a `create*` factory module under `app/*/`.
- **Composition root**: a module whose job is to create owners and connect
  them. Initial registry: `app/app-setup.js`, `app/feature-editor.js`,
  `app/legend.js`, `app/legend-layout.js`. PR B1 characterizes
  `app/feature-editor/svg-actions.js`, `app/annotations.js`, and
  `app/linear-comparisons.js` (each creates ≥3 factories) and the registry
  records the outcome.
- **Owner object**: the value returned by an owner's `create*` factory, or a
  service instance with mutable state (`history`, `historySnapshots`,
  `previewRuntime`, `rulePreparation`, `legendLayout`, `*Actions`).
- **Port**: a single named function injected into an owner by a composition
  root. A port has one direction: the owner that receives it calls it; it never
  returns an owner object.
- **Injection edge**: `root: consumer <- provider` when a composition root
  passes an owner object (not a port) into `consumer`'s factory.
- **Forward-reference closure**: a closure passed at composition that
  references an owner instance assigned later in the same file, including
  `let x = null; … x = createX(...)` (`app/app-setup.js:461,1313`). A closure
  at any depth counts, and so does a late-bound function variable that the
  root reassigns after the owner exists (detector
  `owner-graph.forward-closure.v2`, added after E11).
- **State backdoor**: `state.<name> = <function or owner member>` outside
  `state.js`.
- **Projection domain**: a derived display the Result shows from editor intent
  (feature visibility, label intent, legend order, palette + specific rules,
  strokes, composition offsets). R3 names one projection function per domain.
- **Call shape**: the normalized text of a call site of a projection function
  (arguments and surrounding orchestration on the same statement). R3 holds
  when a domain has one call shape outside its owner.
- **Trigger site**: a call of a heavy derived value's producer
  (`rulePreparation.prepare|prepareDrawn|prepareVisibility|prepareCandidate|evaluate|run|runDrawn`),
  including a call of a producer port (a producer method handed to an owner as
  a function, such as `rulePreparation.run`) in the module that receives it
  (detector `heavy-derived.trigger-site.v2`, added after E11).

## 5. Target authority layout

| Path | Role | Added by |
| --- | --- | --- |
| `tools/web-owner-graph.json` | Registry: composition roots, owner object names, allowed injection edges, ports, projection domains with their projection functions, heavy derived producers | B2 |
| `tools/web-owner-graph-detectors.mjs` | Versioned detectors `owner-graph.*.v1`, `projection.call-shape.v1`, `heavy-derived.trigger-site.v1` | B1 |
| `tools/web-design-rule-guards.json` | Registry: each R-rule id → rule text anchor in `gbdraw/web/CLAUDE.md`, guard test paths, allowlist/baseline locations inside those tests | A1 |
| `tools/report-web-owner-graph.mjs` | Maintainer CLI: per-commit trend table over a revision range (the audit tool) | B1 |
| `tests/web/owner-graph-detectors.test.mjs` | Fixture corpus + untouched-base characterization at `8f450194` | B1 |
| `tests/web/owner-graph-replay.test.mjs` | Historical replay: the merges that introduced each coupling produce NEW observations against their first parent | B1 |
| `tests/web/architecture-contracts.test.mjs` | Extended: every `Guards:` path in CLAUDE.md R-sections is registered; layering section present | A1, C1 |
| `.github/workflows/web-structure-audit.yml` | Daily report on `dev` | F1 |

## 6. Required implementation sequence

Phases A and B are independent and may run in parallel. C and D depend on
nothing but should land before E. E depends on B3 (frozen store exists so each
fix contracts it). F depends on B1.

### Phase A — close the self-authorization hole (GOVERNANCE)

#### PR A1 — register design-rule authority and guards (checker-only)

Purpose: make `gbdraw/web/CLAUDE.md` R1–R12 text and their guard tests
registered authority, and refuse a runtime PR that expands them.

Modify:

- `tools/web-design-rule-guards.json` (new): one entry per rule id R1–R12 with
  `ruleAnchor` (heading text), `guardPaths` (the 34 paths the sections name
  today), and `allowlistLocations` (`{ path, symbol }` for each in-test
  allowlist or baseline constant, for example
  `tests/web/track-layout-transition.test.mjs#TRACK_LAYOUT_WRITERS`,
  `tests/web/option-input-integrity.test.mjs#BASELINE`).
- `tools/check-web-change-budget.mjs`:
  - load the registry from the trusted base; add `gbdraw/web/CLAUDE.md`, the
    registry itself, and every registered guard path to `guardPaths` and
    `authorityPaths`;
  - new Gate rule `design-rule.co-change`: a diff that changes a registered
    rule anchor's section text, a registered guard, or a registered allowlist
    **and** any production path under `gbdraw/web/js/` or `gbdraw/web/index.html`
    fails the Gate, unless every registered change is a contraction (allowlist
    entry removed, baseline constant decreased, guard assertion strengthened by
    deletion only). Expansion requires a separate authority-only PR, as the
    ratchet's "Authority directions" table already requires for
    `tools/web-architecture-rules.json`.
  - report the registered changes under "Governance and authority".
- `tests/web/architecture-contracts.test.mjs`: assert every `tests/...` path
  that appears in a `Guards:` paragraph of an R-section is in the registry, and
  that every registry path exists.
- `tests/web/architecture-ratchet-fixtures.test.mjs`: fixture built from the
  #805 merge diff (`1808d20d^1..1808d20d`, paths
  `gbdraw/web/CLAUDE.md`, `tests/web/track-layout-transition.test.mjs`,
  `gbdraw/web/js/app/feature-editor/placement-actions.js`) expecting Gate FAIL
  with reason `design-rule.co-change`; a contraction-only fixture expecting
  PASS.

Verification:

```bash
node --test tests/web/architecture-contracts.test.mjs tests/web/architecture-ratchet-fixtures.test.mjs
node tools/check-web-change-budget.mjs --base 1808d20d^1 --head 1808d20d   # expected: Gate FAIL design-rule.co-change
node tools/check-web-change-budget.mjs --base 8f450194 --head HEAD         # this PR: Review REQUIRED (checker changed), Gate PASS
```

Acceptance: the #805 replay fails; the current `dev` tip passes; no runtime
path changed; `docs/internal/WEB_CHANGE_POLICY.md` gains the rule in one
paragraph.

Proposed commit title: `governance: register Web design rules and guards as authority`

#### PR A2 — PR template and change-class consistency (authority-only)

Modify `.github/pull_request_template.md`: under GOVERNANCE evidence add
"Design-rule or guard change: rule id, direction (contraction | expansion),
and the separate runtime PR it precedes". Add to
`docs/internal/WEB_CHANGE_POLICY.md` the two-PR sequence for a rule change.
No checker change.

### Phase B — make owner coupling measurable (detector-first)

Revision after B1 (2026-10-05): schema version 1 of
`tools/web-architecture-rules.json` caps the registry at four rules, admits
two kinds, forbids a generic rule kind, and `evaluateArchitectureRuleResult`
in `tools/web-architecture-evaluation.mjs` rejects `FROZEN` (the frozen-store
mechanics are declared, not implemented). Rather than a schema plan, the
baseline uses the design-rule guard mechanism of Phase A: R13 (Phase C)
records the B1 characterization as two registered allowlist literals in
`tests/web/owner-graph-baseline.test.mjs` (`OWNER_GRAPH_BASELINE`, a count
map of `detector|subject`; `PROJECTION_SHAPE_BASELINE`, a count map per
domain). The test fails on a new subject and on a fixed subject that is still
recorded, so the baseline is exact and shrink-only; `design-rule.co-change`
blocks an expansion in a runtime diff. B2 and B3 below are superseded by
PR C1 plus the R13 entry in `tools/web-design-rule-guards.json`;
`tools/web-owner-graph.json` and the rules-registry route stay available for
a later schema plan.

#### PR B1 — owner-graph detectors, fixtures, replay, report CLI (checker-only)

Purpose: observe the unmeasured channels with versioned detectors, report-only,
without touching authority.

Modify / add:

- `tools/web-owner-graph-detectors.mjs` exporting, through
  `WEB_ARCHITECTURE_DETECTORS`-compatible entries:
  - `owner-graph.injection-edge.v1`: in each composition root, for every
    `<ident> = create<X>(…)` call (both `const` and deferred `let`), each
    argument value that is an owner object name → subject
    `root|consumer<-provider`.
  - `owner-graph.forward-closure.v1`: closure arguments whose body references
    an owner instance assigned later in the same file → subject
    `root|consumer.port->provider.method`.
  - `owner-graph.state-backdoor.v1`: `state.<name> = ` with a function, arrow,
    or owner-member right-hand side outside `state.js` → subject `path|name`.
  - `owner-graph.whole-object-port.v1`: factory parameters outside composition
    roots whose name is an owner object name → subject `path|factory|param`.
  - `projection.call-shape.v1`: for each registered domain, the set of
    normalized call sites of its projection functions outside the owner →
    subject `domain|shape`.
  - `heavy-derived.trigger-site.v1`: trigger sites per module → subject
    `producer|path`.
  Owner object names and projection functions are read from
  `tools/web-owner-graph.json`; until B2 lands, the detector module carries the
  same list as a default so B1 can characterize the base.
- `tools/report-web-owner-graph.mjs`: `--range <from>..<to> --first-parent`
  prints the per-merge table of Section 2 (the audit tool), and `--at <sha>`
  prints the subject lists.
- `tests/web/owner-graph-detectors.test.mjs`: synthetic fixtures under
  `tests/web/fixtures/owner-graph/` for each detector (positive, negative,
  deferred-`let` case, port-vs-object case), and untouched-base
  characterization at `8f450194` with these exact expectations:
  - injection edges: 18 distinct (`app/app-setup.js` 13, `app/feature-editor.js` 4,
    `app/legend.js` 1; the audit heuristic reported 20 because it counted the
    two `historyFileStore` edges twice — the detector's exact count is the
    characterization);
  - forward closures: 4 (`app/app-setup.js:1190,1191,1298`,
    `app/feature-editor.js:31`);
  - state backdoors: 3 (`app/app-setup.js:1226,2506,2507`);
  - trigger sites: 15 (`app/app-setup.js` 3, `color-actions.js` 2,
    `rule-actions.js` 1, `svg-actions.js` 1, `app/legend.js` 1,
    `app/svg-styles.js` 1, `app/run-analysis.js` 2, `visibility-actions.js` 2,
    `label-actions.js` 2);
  - call shapes: feature visibility+label 3 orchestration shapes
    (`visibility-actions.js`, `feature-editor.js:64-68`,
    `app-setup.js:2685-2689`); legend order 4; palette+rules 3; strokes 4.
  Characterize the three candidate roots and record whether each is a root.
- `tests/web/owner-graph-replay.test.mjs`: for each of `#617 f5f86634`,
  `#622 27939fae`, `#641 b1744142`, `#692 942fa954`, `#737 cafed980`,
  `#764 0bf9034a`, `#805 1808d20d`, `#807 4cf78f5d`: evaluate detectors at
  `<merge>^1` and `<merge>` and assert at least one NEW subject in the detector
  the Section 2 row attributes to it. This is the acceptance test for the whole
  plan: each historical coupling PR becomes observable.
- `tools/check-web-change-budget.mjs`: no behavior change; print the six
  detector outputs under a new report-only heading "Owner graph (report)".

Verification: `node --test tests/web/owner-graph-*.test.mjs`;
`node tools/report-web-owner-graph.mjs --range e11e1d02..8f450194 --first-parent`
reproduces Section 2 within the stated numbers.

Acceptance: all characterizations exact; replay passes for all eight merges; no
registry, policy, or runtime change.

Proposed commit title: `tools: add owner-graph detectors, replay test, and trend report`

#### PR B2 — register the owner graph (authority-only, report-only)

Modify:

- `tools/web-owner-graph.json` (new): composition roots (from B1's
  characterization), owner object names, the current injection edges (B1 characterization) as
  `allowedEdges`, projection domains with projection functions
  (`feature-visibility`: `projectFeatureVisibility`/`reconcileFeatureVisibility`,
  `applyFeatureVisibilityToLabels`; `label-intent`: `reconcileLabelOverrides`,
  `projectLabelIntent`; `legend-order`: `orderLegendEntries`,
  `reconcileLegendEntries`; `palette-rules`: `applyPaletteToSvg`,
  `applySpecificRulesToSvg`; `strokes`: `reconcileStrokeOverrides`,
  `applyStrokeOverridesToSvg`; `composition`: `reconcileCompositionUserDeltas`),
  heavy producers (`rulePreparation`).
- `tools/web-architecture-rules.json`: six rules, `enforcement: report-only`,
  `baselineEligible: true`, referencing the B1 detector ids.
- `tools/web-architecture-evaluation.mjs`: accept the new rule kinds
  (`allowed-edge-set`, `forbidden-subject`, `shrink-only-subject-set`).
- `tools/check-web-change-budget.mjs`: add the three new paths to `guardPaths`
  and `authorityPaths`.

Acceptance: Gate PASS on `dev`; report lists exactly the B1 characterization.

Proposed commit title: `governance: register owner-graph rules as report-only`

#### PR B3 — freeze the owner graph (authority-only)

After at least five merged runtime PRs have produced reports with no detector
error, tighten:

- `owner-graph.forward-closure`, `owner-graph.state-backdoor`: `frozen`; the
  4 + 3 current subjects go into `tools/web-architecture-violations.json`.
- `owner-graph.injection-edge`: `frozen` against `allowedEdges`; a new edge is
  a NEW violation (Gate FAIL) until an authority-only PR adds it.
- `owner-graph.whole-object-port`: `frozen` per `path|factory|param` (58
  subjects).
- `projection.call-shape`: `frozen` per domain as a shrink-only count.
- `heavy-derived.trigger-site`: `frozen` per `producer|path` (15 subjects).

Acceptance: Gate PASS on `dev` with zero NEW; a fixture PR that adds one
injection edge fails; a fixture PR that removes one accepted subject without
contracting the store fails (`FIXED` without authority contraction).

Proposed commit title: `governance: freeze owner-graph baselines`

### Phase C — define layering and ports in the authority text (GOVERNANCE)

#### PR C1 — "Owner layering and ports" section in `gbdraw/web/CLAUDE.md`

Add, after "Module ownership", a section stating:

- layers: `state.js` → `services/` → owner modules → composition roots →
  template; a module depends only on lower layers;
- an owner never imports, holds, or calls another owner object; it receives
  ports; a port is one function, one direction, named for the reaction it
  performs (`applyFeatureVisibilityToLabels`), never for the provider
  (`labelActions`);
- a reaction of owner A to a change in owner B is wired by the composition
  root: either as a port injected into B, or as the root calling A's single
  projection after B's transition returns; never by B calling A's object and
  never by A calling back into B;
- a forward-reference closure is a composition error: reorder creation, or
  register the port after both owners exist (`historySnapshots.registerCapture(domain, fn)`);
- `state.<name> = fn` is forbidden; a value an owner needs from a higher layer
  arrives as a port or on an existing data path (Result metadata, catalog);
- each projection domain has one projection function and one call shape
  outside its owner (R3 restated as a shape, not an existence);
- a heavy derived producer is triggered from one place per flow: the
  composition root's projection, the Generate compiler, or the owner that
  owns the input; an owner that only reads the result awaits the producer's
  `pending` rather than triggering it.

Amend R10's first paragraph: "calls the owner's reconcile explicitly, through
the port the composition root injected or by returning to the root, which
calls the reconcile (see Owner layering and ports)". Leave the #805 paragraph
and mark it with the accepted-violation id from B3; it is re-decided in E7.

Guards: `tests/web/architecture-contracts.test.mjs` asserts the section and the
registry file names appear; `tools/web-design-rule-guards.json` gains the
section as rule id R13 with `tests/web/owner-graph-detectors.test.mjs` as its
guard.

Proposed commit title: `docs(web): define owner layering and ports; restate R3 and R10 as shapes`

### Phase D — decision records name the reaction owner (GOVERNANCE)

#### PR D1 — reaction fields in decision records

Modify `docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md` ("Choice" sections)
and the record format in `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`:
for a decision in which one editor domain follows another, the record must
carry `Reaction owner: <module>`, `Channel: port | root-projection | event`,
and `Projection: <function>`. `tests/web/product-impact-ratchet-fixtures.test.mjs`
checks that every decision dated after this PR whose text matches
`(follows|hides|shows|with it|same action)` carries the three fields.
Backfill F-3 (label follows feature), Q2 (Show feature and label), Q3 (slot
inputs through placement) with the owners chosen in Phase E.

Proposed commit title: `governance: decisions that couple owners name the reaction owner and channel`

### Phase E — contract the frozen baselines (STANDARD, one PR each)

Each PR removes the accepted subjects it fixes from
`tools/web-architecture-violations.json` in the same diff (the protocol's
"remove each accepted entry in the same runtime PR that fixes its subject"),
and may not add any. Order by blast radius; E1 unblocks the open #810.

| PR | Subject | Change | Removes |
| --- | --- | --- | --- |
| E1 | feature visibility ⇄ label | One `projectFeatureVisibility` owned by `visibility-actions.js` that projects features and then calls the injected port `applyFeatureVisibilityToLabels`; `feature-editor.js:64-68` and `app-setup.js:2685-2689` call only it. The rerender request leaves the label API: a `requestAutomaticRerender` port owned by the composition root replaces `applyFeatureVisibilityToLabels({ rerender })` (`label-actions.js:759-765`, `visibility-actions.js:496`); `watchers.js:270-287` consume one `automaticRerenderRequestSeq`. `confirmLabelOn` (`label-actions.js:1035,1046`) stops calling `setFeatureVisibility`; the composition root provides `showFeatureForLabel` as a port that runs the visibility owner's transition and returns a restore function. `visibility-actions.js:266-271` requires `previewRuntime`. Rebase #810 onto E1 or fold it in. | forward closure `feature-editor.js:31`; call shapes feature-visibility 3 → 1; trigger sites `label-actions.js` 2 → 0 |
| E2 | History ⇄ Legend/LegendLayout | `historySnapshots.registerCapture('legend', fn)` and `('composition', fn)` called from `app-setup.js` after `legendActions`/`legendLayout` exist; `reconcileLegendEntries` takes `{ entries, entryOwners, from }` built by the root from the History step, and the B19 decision (`app-setup.js:2711-2718`) moves into the legend owner as one function; `entry-actions.js:796-831` routes through `orderLegendEntries`. | forward closures `app-setup.js:1190,1191`; legend-order shapes 4 → 1 |
| E3 | Legend ⇄ Feature rules | `app-setup.js:1298` shim replaced by a port `commitLegendRowRules` registered after both owners exist; `retiredLegendIntents` leaves `rule-matching.js:213`, computed in `rule-actions.js` from the candidate result; `entry-actions.js:947` no longer executes the rule mutation inside the legend transaction (the root runs the rule transition, then the legend sync). | forward closure `app-setup.js:1298`; trigger-site param leakage |
| E4 | `rulePreparation` fan-out | `run-analysis.js:2136-2140/4475/4482/4774` and `:5322-5324/5345/5348/5390` share one `prepareAndAdmitCandidate`; palette+rules projection has one caller (`projectMountedEditorIntent`), `svg-styles.js:490-492` and `rule-actions.js:126-129` call it; single-caller parameters `onProgress`, `{drawn}`, return fields `previousIntents`, `featureColorOverrides`, `changes` removed or generalized. | trigger sites 15 → ≤9; palette-rules shapes 3 → 1 |
| E5 | Generate → editor state | `run-analysis.js:1670-1677` → label port `closeLabelTextScopeDialog`; `:5385-5387` → `previewRuntime.selectResult`; `:4874-4876` → placement owner `restorePlacements` transition in `services/feature-placement.js`; `:2132,5316` warning clears → label owner port. | no detector subject; R7/R2 guard baselines decrease |
| E6 | state backdoors | `committedDiagramOptions` travels on the admitted catalog (`services/feature-catalog.js`); `recordDisplayRows` becomes an explicit argument of `buildCanonicalRenderRequest` (both callers already construct the argument); `sessionPreparationBusyReason` is composed in `app-setup.js` and injected into `services/config.js` and `history` as `mutationAvailability`. | state backdoors 3 → 0 |
| E7 | TrackSlots ⇄ Placement | Owner decision required (see Section 8): either placement receives slot-owner setters as ports (slot owners stay the writers of `form.track_type`, `linear_track_layout`, `separate_strands`), or the #805 direction stands and `trackLayoutActions` becomes a registered port rather than a wrapper of the slot editors. | whole-object port `circular-track-slots.js:1409`, `linear-track-slots.js:822` |
| E8 | run-analysis ⇄ record-display | `isCurrentFeature` becomes a pure function in `services/feature-identity.js` taking the record rows; `feature-record-rotation.js:101` requests the run through a root-provided port. | injection cycle |
| E12 | `rulePreparation` trigger sites outside the allowed places | Triggers sit in the root projection, the Generate compiler, and the owner of the input (R13). E12a: the popup (`svg-actions.js`) and Label On (`label-actions.js`) ask the visibility owner's `prepareDrawnFeatureMatches` port (the owner's `prepareDrawn`, `strict` as a run); Label TSV import takes the root's `evaluateLabelRules`; `runDrawn` and the owner's `evaluate` leave `rulePreparation`. E12b: colour actions take `runWithRuleMatches` from the rule owner (`rule-actions.js`), which prepares through its one `prepareCandidate` funnel; `run` leaves `rulePreparation`. | trigger sites `svg-actions.js` 1, `label-actions.js` 2 (E12a), `color-actions.js` 2 (E12b); the remaining five sites are the allowed places |

Each E PR: STANDARD; "Architecture impact: architecture-bearing (ordinary,
non-increasing)"; evidence is the contracted store plus the detector report
before/after; no `gbdraw/web/CLAUDE.md` or guard change (A1 would fail it).

### Phase F — continuous structural audit

#### PR F1 — daily report on `dev`

`.github/workflows/web-structure-audit.yml`: `schedule: cron '15 20 * * *'`
(05:15 JST) and `workflow_dispatch`; checks out `dev`; runs
`node tools/report-web-owner-graph.mjs --at HEAD --json` and
`--range <previous-run-sha>..HEAD --first-parent`; uploads the JSON as an
artifact; writes the table to the job summary; fails the job when any frozen
metric is above the previous artifact's value (trend check independent of the
per-PR gate). `permissions: contents: read`. Add the job name to
`docs/internal/WEB_PERIODIC_AUDIT.md`.

Proposed commit title: `ci: daily Web owner-graph structural audit`

## 7. Acceptance for the plan as a whole

- `tests/web/owner-graph-replay.test.mjs` passes: all eight historical coupling
  merges are observable as NEW subjects (B1).
- The #805 replay produces Gate FAIL `design-rule.co-change` (A1).
- After B3, a PR that adds an injection edge, a forward closure, a state
  backdoor, a trigger site in a new module, or a second call shape for a
  registered domain fails the trusted-base Gate from `dev`'s checker, with the
  PR unable to amend the registry in the same diff.
- After E1–E8, the frozen store holds zero owner-graph subjects and the daily
  audit's frozen metrics equal the Section 2 "2026-09-01" column or better.

## 8. Open owner decisions

1. **E7 direction.** R10's #805 paragraph makes `placement-actions.js` the
   writer of three settings inputs and wraps the slot editors. Keep that
   direction (then register `trackLayoutActions` as a port and move the three
   setters into placement) or restore slot-owner writers with placement as a
   port consumer. Both satisfy the layering; the first keeps #805's product
   behavior with less code churn.
2. **Rerender channel name and owner (E1).** `runLabelReflow` is today the
   automatic rerender (R1(c)) for every domain. Proposal: rename to
   `runAutomaticRerender`, owned by `app/run-analysis.js`, requested through
   one port; label reflow becomes one caller.
3. **Freeze timing (B3).** Five merged runtime PRs with clean reports, or a
   fixed date; the plan assumes five PRs.

## 9. Non-goals

- No TypeScript migration; the registries and detectors give the same
  visibility for this codebase without a build step.
- No big-bang refactor; every E PR is independently revertible and contracts
  the store.
- No change to product decisions F-3, Q2, Q3; E1–E7 change where the reaction
  lives, not what the user sees. E7's direction is the one exception and is
  listed in Section 8.
- No watcher-removal campaign; repairing watchers are already decreasing and
  are covered by R10 as written.

## 10. Rollback

Each PR is revertible alone. Reverting B3 returns the rules to report-only;
reverting B2 removes them; reverting A1 returns the design rules to
unregistered authority (the state this plan corrects). An E PR revert must also
restore the accepted subjects it removed, which the frozen-store protocol
requires the revert PR to carry.
