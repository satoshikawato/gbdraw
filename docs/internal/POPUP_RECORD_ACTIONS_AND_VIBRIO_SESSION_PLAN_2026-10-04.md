# Feature popup record actions and Vibrio Session plan (2026-10-04)

Status: approved by the Product Decision Owner `satoshikawato` on 2026-10-04.
Base: `dev` at `0bf9034a`. Line numbers below are from `f2a7f8d4`/`0bf9034a`;
re-locate them by symbol before editing.

## 1. Owner request and decisions

The Owner loaded the Gallery Session "Vibrio harveyi group" on gbdraw.app
(`v0.14.0+fe6861f`) and reported, verbatim:

> ExamplesのVibrio Harveyi groupの.jsonファイルを読み込むと、
> - "File order is unavailable because Record Layout is custom. Use Advanced comparison and layout → Record Layout to restore one row per File with no shared rows." なぜこのような警告が出る？"Generate Diagram"しないと"RECORD ACTIONS"ができない。
> - レコードごとに.gbkファイルが分割されている！これはおかしい。もともと1ファイルだったものは1ファイルに複数レコードを保持しなければならない。
> - フィーチャーのポップアップ画面がごちゃごちゃしてわかりづらい。
> - "Apply and regenerate"以外に、ただ単にアクションを予約というかセーブするというか、複数のフィーチャーについて"RECORD ACTIONS"を適用して、最後に一気に描画するオプションが欲しい。ApplyとOKの違いというか。あるでしょ画像レタッチソフトとかで。

The Owner asked for a fix plan that respects SOLID, KISS, DRY, and YAGNI, then
answered the three open questions and approved the plan:

| ID | Question | Owner answer |
| --- | --- | --- |
| OD-1 | Move Record actions from the top of the Edit tab into a Layout group? | 「移しましょう。」 |
| OD-2 | Name the staging button **Apply on Generate**? | 「推奨通りでお願いします。」 |
| OD-3 | Keep **Apply and regenerate** target-only (`PD-OI-032` item 4), so staged changes on other records stay pending? | 「推奨通りでお願いします。」 |
| OD-4 | Proceed | 「計画書をまとめた後、実装に移ってください。最初に私が貼り付けたスクリーンショットをファイルとして保存したり、コミットすることは禁止します。」 |

OD-4 constraint: the Owner's screenshots are never saved as files or committed.
Evidence images, if needed, are captured fresh from local builds.

After R1 merged, the Owner added a request about the rotation controls
(OD-5), verbatim:

> Record actionsで、Ｐｌａｃｅ　ｔｈｉｓ　ｆｅａｔｕｒｅ　ａｔ　ｔｈｅ　ｅｎｄがでかすぎませんか？数あるオプションの一つでいいというか、むしろどちらかというとPlace this feature at the start とかのほうが優先されるべきでは？Anchorとの関係性がわかりにくい。ここらへんはもう抜本的にデザインしなおして下さいわかりやすいように。

The Owner delegated the redesign. Section R4 "Rotation controls" records the
design. It re-presents the operations that `PD-OI-032` item 2 and the
`PD-OI-033` revision 2 receipt preserve (anchors, signed offset, orientation,
feature-end placement) and adds a "start of the record" shortcut over the
existing anchors; it changes no persisted format.

## 2. Findings

### F1. The Vibrio Gallery Session stores each chromosome as its own File

- `gbdraw/web/gallery/examples.json` (entry `vibrio-harveyi-group-collinear`)
  declares `featureSources` = the two GBFF files and `inputSummary` = "2
  multi-record GenBank files; 4 chromosomes". Its command uses
  `examples/vibrio-harveyi-group-linear-records.tsv`, which points both rows of
  each species at the same unsplit GBFF under `tests/test_inputs/`.
- The stored Session (`gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz`,
  version 44) has four resources, `record-1-genbank` .. `record-4-genbank`,
  each holding one LOCUS, with original names `NC_004603.1__GCF_...gbff` etc.
- Origin: the first Session (`229b2576`, 2026-07-18) already had per-record
  resources; the CLI codec then wrote one resource per record. `c9426c77`
  (2026-10-03) made the codec keep one source file as one resource, but it
  landed after the last Gallery refresh. `tools/refresh_gallery_sessions.py`
  replays the stored Session, and grouping keys on source paths, so replay
  cannot merge the split resources. The split persists until the Session is
  rebuilt from its declared inputs.

### F2. The custom-layout warning is a consequence of F1

`planLinearSourceRowMove` (`app/linear-record-layout.js:52-80`) returns
`custom-layout` when two Files share a row. Four Files on rows 1, 1, 2, 2 hit
that rule, so `linearSourceMoveBlockedReason` (`app/app-setup.js`) shows "File
order is unavailable because Record Layout is custom." Two Files with two
records each fall under the default "one row per File" and show no warning.

### F3. Record actions do not work after Session Load until Generate

- `776a2f93` made a current-schema Session Load read no record bytes: "the first
  Generate or an explicit record action reads them". `importSession` sets
  `sessionResourceDiscoveryDeferred`, discovery watchers are suppressed, and
  nothing restarts discovery afterwards.
- The popup rotation target is built only from discovered record rows.
  `targetForFeature` (`app/record-display-options.js:402-408`) filters
  `committedRows` by `record_key` and throws "The popup feature target is stale
  or ambiguous." when there is not exactly one row. After Load there are zero
  rows. Generate works because `prepareGenerate` runs discovery first.
- The popup never starts that "explicit record action" read, and the reason
  text is wrong: the Result matches the inputs; the records are only unread.
- `committedRows` is a plain variable and the draft is computed once on popup
  open, so nothing recomputes when discovery finishes later.

### F4. One failure is shown four times

The `recompute` catch in `app/record-display/feature-record-rotation.js` copies
`error.message` into `anchorChoices[].message`, `orientCapability.message`,
`featureEndCapability.message`, and `disabledReason`; the template renders all
four, plus "Record: Unavailable" and "New display start: Unavailable".

### F5. The popup repeats itself and mixes apply timings

- Two identical fill color inputs (`aria-label="Feature fill color"`): one in
  the header, one under Fill Color.
- The similarity group appears in the header (`og_10 - 2 members`) and in the
  Similarity group section.
- The disclosure button "Record actions · Rotate record using this feature" is
  followed by an inner heading "RECORD ACTIONS / Rotate record using this
  feature…", a state badge (`placementLabel`), and Record/Feature rows that
  repeat the header.
- Feature placement (applies on Generate) sits above the tabs; the Edit tab
  starts with a "Live edit: …" paragraph; controls with three timings (live,
  Apply Label, on Generate, regenerate now) are interleaved.

### F6. The Generate-time queue already exists

- `recordDisplayDrafts` in `app/record-display-options.js` holds one draft per
  `[scope, sourceUid, selector]`; the sidebar's record display controls edit it
  and Generate applies every draft.
- Apply and regenerate already ends in `writeResolvedTransform`, which writes
  the resolved start/orientation into that draft (single-record Linear cards
  write `seq.region_reverse`).
- `hasPendingChanges` drives the existing "Record rotation, feature placement,
  or tolerance has changes pending Generate." notice in the Input Genomes card.
- `captureTargetDraft` / `restoreTargetDraft` already give a History checkpoint
  for one target's draft.

So staging needs no new queue, store, or History engine.

### F8. The rotation controls hide what they do (OD-5)

- Every operation sets one value: which source base becomes the first base of
  the displayed record. The form does not say so. It shows an Anchor select
  (5′ end, midpoint, 3′ end), an Offset field, and a separate full-width
  "Place this feature at the end" button that switches to a different mode
  (`placement: 'feature-end'`) and disables the Anchor select.
- The common goal, "put this feature at the start of the record" (for
  example, dnaA at position 1), has no control. The 5′-end anchor does it
  only when the feature reads forward in the resulting display. For a
  minus-strand feature without orientation change, the 5′ end is the
  feature's rightmost base, so the feature is split between both ends of the
  display. The display-first base in that case is the 3′ end.
- `resolveFeatureAnchor` (`app/record-display/feature-anchor.js`) already
  computes the displayed strand after the operation and the feature-end
  boundary. The persisted provenance accepts only `placement` `anchor` or
  `feature-end` (`app/record-display-options.js`, anchor intent validation).

### F7. Closing the popup during Apply and regenerate leaks state

`closeFeaturePopup` resets the draft without aborting the run; when the
in-flight `apply` resolves, it writes its status into the reset or the next
popup's draft (`feature-record-rotation.js`, `workflow.apply`).

## 3. Rules for every PR

- Reuse the owners named in F3 and F6. Add no new store, queue, History engine,
  Worker path, request owner, or record rotation engine (`PD-OI-032` item 10).
- One reason source: a failure is computed once and rendered once.
- Group popup controls by when they apply; the group heading states the timing
  once (`02_DECISION_PACK.md`: every setting is truthfully labelled Live edit
  or Applies on Generate).
- Load still reads no record bytes (`776a2f93`); only an explicit record action
  may start discovery.
- Delete superseded paths, constants, and tests in the same PR.
- Follow `gbdraw/web/CLAUDE.md`, `docs/internal/WEB_CHANGE_POLICY.md`, and the
  PR template. Branch from the latest `origin/dev`.

## 4. PRs

| ID | Branch | Class | Scope | Depends on |
| --- | --- | --- | --- | --- |
| R0 | `docs/popup-record-actions-plan` | STANDARD (docs) | This plan | none |
| R1 | `docs/product-contract-popup-record-actions` | GOVERNANCE | Product Contract revision 30: `PD-OI-033` scenario revision 2 and `PD-OI-085` from Appendix A | R0 |
| R2 | `fix/web-record-actions-after-session-load` | STANDARD | F3, F4 | none |
| R3 | `fix/gallery-vibrio-multi-record-files` | STANDARD | F1, F2, recurrence guard | none |
| R4 | `feat/web-feature-popup-layout-groups` | STANDARD | F5, F8, OD-1 (`PD-OI-033` rev 2), OD-5 | R1, R2 |
| R5 | `feat/web-record-rotation-apply-on-generate` | STANDARD | OD-2, OD-3 (`PD-OI-085`), F7 | R1, R4 |

R2 and R3 run in parallel. R4 and R5 both edit the popup template and
`feature-record-rotation.js`, so they run in sequence. Open PR #762 (Feature
placement) also touches `record-display-options.js` and `run-analysis.js`;
rebase on it once it merges.

### R1. Product Contract revision 30

- Add Appendix A's two receipts to
  `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` in the form of
  `PD-OI-084`: receipt text block, nine-field JSON, UTF-8 SHA-256 of the
  receipt (final newline excluded), and a Revision 30 entry that quotes the
  Owner answers in section 1 verbatim.
- `PD-OI-033` scenario revision 2 supersedes revision 1
  (`A / EDIT_DISCLOSURE`). `PD-OI-085` is new; add it to the approved IDs.
- Authority only; no runtime.

### R2. Record actions after Session Load (F3, F4)

1. Give the target failure a kind. "Records not discovered yet" (no discovered
   rows for the feature's source) is distinct from a real identity mismatch.
   Only the latter says "stale or ambiguous".
2. When the Record actions section is open and the failure kind is "not
   discovered", start the existing discovery for the mode
   (`linearRecordSelector.refresh()` in Linear, `refreshCircularRecordOrder()`
   in Circular), show one "Reading records…" status with the form disabled,
   and recompute the draft when discovery settles. A discovery error shows its
   own reason. This is the explicit record action that `776a2f93` allows; Load
   itself still reads nothing.
3. F4: `recompute` puts a whole-target failure only in `disabledReason`.
   Per-control messages are set only for control-specific limits (for example,
   orientation unavailable for a strandless feature). When there is no target,
   the template shows the reason line and hides the form.
4. If the File card's Record options are empty after Load for the same reason
   (check in a browser), opening them uses the same entry point. Otherwise
   leave them.
5. Tests: a unit test that the draft reports "not discovered" (not stale)
   with zero rows; a Playwright test that loads a current-schema Session,
   opens a feature popup and Record actions without Generate, and sees the
   record ID and new display start; the existing "Load reads no record bytes"
   tests still pass.

### R3. Vibrio Session with two multi-record files (F1, F2)

1. Rebuild the Vibrio Gallery Session from its declared inputs: the records
   table and the two GBFF files. Result: two GenBank resources, each with both
   LOCUS records, records bound by record-ID selectors, File 1 on row 1 and
   File 2 on row 2, no custom-layout warning, File moves enabled.
2. Change the generator, not only the bytes: after this PR,
   `tools/refresh_gallery_sessions.py` must reproduce the two-file shape from
   current code and declared inputs. Remove superseded per-record constants
   (`_VIBRIO_RECORD_KEYS`, per-record pair expectations, `VIBRIO_RAW_ENTRY_COUNT`
   if it changes meaning) and size limits that no longer apply.
3. Guard: one test asserts that every Gallery Session's resource original
   names match its entry's declared inputs (`featureSources` plus any other
   declared input files). This prevents per-record splits in any entry.
4. LOSAT uses the file-level database (`PD-OI-018` revision 4), so comparison
   output can change. Regenerate the entry's SVG, thumbnail, source figure,
   and tutorial captures with the owner tools
   (`.agents/skills/web-gallery-screenshot-maintenance/SKILL.md`), inspect them
   at readable scale, and mark Review REQUIRED for scientific output.
5. Update tests that pin four files (`tests/web/contracts/vibrio-sequence-source-coverage.test.mjs`,
   `tests/web/session-request.test.mjs` near 5072,
   `tests/web/vibrio-session-save.performance.playwright.spec.js`
   `resourceCount`, `tests/test_refresh_gallery_sessions.py`,
   `tests/web/contracts/vibrio-full-generation.serial.spec.js`). Grep main-only
   and dispatch suites as well.

### R4. Popup layout groups (F5, OD-1)

Target layout (rich popup; the simple popup shows the Edit content without
tabs):

```text
■ <feature label>                                     [×]
  <record ID>: <location> · <gene>
 [Edit] [Details] [Qualifiers] [Sequence]
 Appearance · updates the current Result
   Fill Color, Stroke, Label text + Label visibility (Apply Label),
   Feature visibility, Legend name
 Layout · applies on Generate
   Feature placement
   ▸ Rotate record using this feature   (disclosure, collapsed)
 Similarity group <id> (<n> members)     (existing actions)
```

- Remove the header fill color input (keep the swatch), the header similarity
  group lines, the "Live edit: …" paragraph, and the per-control "Applied on
  the next successful Generate." caption.
- Move Feature placement into the Layout group.
- In the rotation section, remove the inner heading, the `placementLabel`
  badge, and the Record/Feature rows; the header carries the record ID.
- Keep keyboard operation and a 390 px viewport working; update tests that
  depend on the old role names and texts.
- Docs: `docs/REFERENCE/web-app.md` popup description, `CHANGELOG.md`.

#### Rotation controls (F8, OD-5)

One question, one answer: where does this feature go in the record. The
presets name outcomes; Custom exposes the raw rule "the record starts at a
reference point of this feature, shifted by an offset".

```text
▾ Rotate record using this feature
  Put this feature at
    (•) Start of the record
    ( ) End of the record
    ( ) Custom position
          Record starts at [this feature's 5′ end ▾] shifted by [ 0 ] bp
            reference choices: 5′ end · midpoint · 3′ end · just after the feature
  [ ] Show this feature on the forward strand
      (reverse-complements the record when the feature is on the − strand)
  NC_004603.1 will start at 7,680 · orientation unchanged       (i)
  [Cancel]                       [Apply and regenerate]   (+ R5 button)
```

1. The position radio group replaces the Anchor select, the Offset field at
   the top level, the full-width "Place this feature at the end" button, and
   the `placementLabel` badge. Default: Start of the record.
2. Start of the record: the feature's first base in display order becomes
   base 1. In `resolveFeatureAnchor`, a `feature-start` intent resolves to
   the five-prime anchor when the feature reads forward in the resulting
   display (displayed strand after the operation is `+`) and to the
   three-prime anchor otherwise, with offset 0. Its provenance is recorded as
   `placement: 'anchor'` with that resolved anchor, so no persisted format
   changes. Unavailable, with its own reason, for features without a known
   strand.
3. End of the record: the existing `feature-end` placement with offset 0.
4. Custom position: reference select (5′ end, midpoint, 3′ end, just after
   the feature) plus the signed offset, counted along the feature's strand
   (help tip). "Just after the feature" is `feature-end` with the offset. This
   keeps every operation of `PD-OI-032` item 2, including feature-end with an
   offset.
5. Orientation checkbox: unchanged meaning, clearer label. Its effect on
   "Start of the record" is resolved by item 2, so the feature lands at the
   start with or without it.
6. Preview: one sentence with the record ID, the new start, and the
   orientation; the displayed-strand change appears only when it changes.
   "Coordinates refer to the original record." moves into a help tip.
7. An unavailable choice is disabled with its own short reason next to it.
   A whole-target failure shows one reason line (R2) and hides the form.
8. The domain math stays in `feature-anchor.js`; the workflow in
   `feature-record-rotation.js` only maps the radio choice to an intent.
   Remove `placeAtFeatureEnd`, `exactFeatureEnd`, and `placementLabel` when
   nothing else uses them.
9. Tests: `feature-start` resolution for +/− strands, with and without
   orientation change and with an already reversed record (feature occupies
   display bases 1..n); unstranded is unavailable; Custom covers all four
   references with offsets; Playwright for keyboard operation and 390 px.

### R5. Apply on Generate (OD-2, OD-3, F7)

1. Buttons: `Cancel`, `Apply on Generate`, `Apply and regenerate`. Both apply
   buttons use the same validation (`canApply`).
2. `Apply on Generate` writes the resolved transform through
   `writeResolvedTransform` into the target's draft and records one History
   step with the existing target-draft checkpoint. It runs no candidate and
   does not touch the Result. The status line says the change applies on
   Generate.
3. Re-staging the same record replaces its draft. When the popup opens on a
   record with a pending draft, the preview shows the pending start and
   orientation.
4. `Generate Diagram` applies all staged drafts (existing behavior). The
   existing pending notice signals staged changes; no badge is added to the
   Result or the Generate button (`PD-OI-037` revision 2).
5. `Apply and regenerate` stays target-only (`PD-OI-032` item 4, OD-3).
6. F7: an in-flight Apply and regenerate whose popup was closed or retargeted
   writes nothing into the current draft.
7. Tests: stage two records, then Generate once and see both rotated; Undo and
   Redo of a staged step; staged record B survives Apply and regenerate on
   record A; F7.
8. Docs: `docs/REFERENCE/web-app.md`, `CHANGELOG.md`.

## Appendix A. Product receipts for R1

Both receipts record the Owner answers in section 1 (OD-1 to OD-4) on
2026-10-04.

```text
PRODUCT_DECISION
Concern: web.feature-popup.record-actions-presentation
Scenario revision: 2
Choice: B / LAYOUT_GROUP_DISCLOSURE
Rationale: popup の Edit を、現在の Result に反映する Appearance と、Generate で適用する Layout の 2 つに分ける。record 回転は Feature placement と同じ Layout グループの開閉セクションに置く。重複した表示と、同じ理由文の繰り返しをなくし、popup を読みやすくする。
Must preserve: 開いた feature だけを対象とする回転、既存の anchor・offset・orientation・feature-end 操作、適用前 preview、操作できない理由の表示（1 か所）、rich と simple の両 popup での同じ操作、keyboard と 390 px での到達性、既存 sidebar 操作、成功時の一体的な Result と Undo/Redo、Cancel・失敗時の直前 Result と record transform、Feature placement の次回 Generate 適用、fill color・stroke・label・feature visibility・legend name・similarity group の操作。
May retire: Edit タブ先頭の Record actions 配置、セクション内の重複見出し・Record と Feature の行・状態バッジ、同じ理由文の複数表示、ヘッダの fill color 入力（Fill Color と重複）、ヘッダの similarity group 行（Similarity group セクションと重複）、タブより上の Feature placement 配置とその個別注記、Edit 上部の「Live edit: …」説明文（グループ見出しで置き換える）。
Accepted residual risk: record 回転は Edit の上部から下へ移るため、見つけにくくなる。Layout 見出しと開閉ボタンの名前で補い、keyboard と 390 px で到達できることを確認する。
Owner: satoshikawato
Decision date: 2026-10-04
```

```text
PRODUCT_DECISION
Concern: web.feature-popup.record-rotation-apply-on-generate
Scenario revision: 1
Choice: A / APPLY_ON_GENERATE
Rationale: Record actions に Apply on Generate を加える。開いた feature の record の表示開始位置と、指定したときの向きを、sidebar と同じ未適用の record display 設定に書き込むだけで、Result は再生成しない。Generate Diagram を 1 回押すと、予約したすべての record がまとめて描画される。複数 record の回転を予約してから一度に描画したいという要望（画像編集ソフトの「適用」と「OK」の区別）に応える。
Must preserve: PD-OI-032 の対象特定・anchor・offset・orientation・feature-end・適用前 preview・理由表示。Apply and regenerate は最後の committed request から対象 record だけの candidate を作り、他の record の予約は予約のまま残す（PD-OI-032 item 4）。sidebar の record display 操作とその意味。Generate Diagram による未適用設定の一括適用。同じ record への再予約は後の値が有効。予約は 1 回の Undo/Redo で戻せる。popup を開き直すと予約済みの値が分かる。既存の pending Generate 通知のほかに、Result と Generate への常時 Pending/Applied 表示を加えない（PD-OI-037 revision 2）。
May retire: none
Accepted residual risk: 予約したまま Generate しないと、表示中の Result と設定がずれたままになる。既存の pending Generate 通知と popup での予約値表示で補う。
Owner: satoshikawato
Decision date: 2026-10-04
```
