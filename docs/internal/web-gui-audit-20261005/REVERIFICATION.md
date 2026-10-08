# Web GUI 再監査（2026-10-05）の確かめ直しと、修正の計画

`FINDINGS.md` の 59 件を、2026-10-07 の `dev` `9b0c56bd` で確かめ直した。mode、Result、Load、凡例に関わるものは、E1 の head（`fix/web-per-mode-results` `76b45de9`）でも確かめた。計画の元は `PLAN.md`。順序と範囲は、下の「PR の計画」で置き換える。

## 結果の要約

- 別の PR ですでに直っていた: TK-05（#805）、UJ-03/FL-01（#857、#870）、FL-03 と FL-06（#892）、FL-04（#894）、FL-13（#812）。
- E1 で直る: UJ-02（E1 の受け入れの case に足した）。
- 別のセッションに渡した:
  - 凡例の FL-02、FL-05、FL-07 は Phase E の LEGEND の PR に渡した（OV-152、OV-154、OV-127）。
  - UI-10 は v0.15.0 の multi-drawing に渡した（Owner の決定、2026-10-07）。
- そのほかは、再現するか、症状が変わって残っている。

## PR の計画（Owner の指示、2026-10-07: PR の本数を減らす。0.14.0 に要らない作業は止める）

| PR | 中身 | 時期 |
|---|---|---|
| A-out（#939、マージ済み 3fbfac97） | TK-01（P1）、TK-02、外側の行の順序（Owner の決定）、入力と比較のエラー（CI-04、CI-06 ほか）、Gallery の refresh tool（UJ-04、GX-03）、この記録 | Phase E の PER-MODE の前 |
| A-ui（#941、マージ済み 3f09501c） | UI-01 の刷新（提案）、UI-05、UI-06、Fit ボタン（UI-11、GX-02） | Phase E の PER-MODE の前 |
| B = Lane B（`fix/web-gui-audit-b`） | TK-04、TK-06〜10、TK-12、TK-13、TK-15、UI-07、CI-02、CI-03、UJ-01、UJ-06、UJ-08〜10、FL-14、観察 2 件（mode ごとの検索、Linear の Auto の通知）、GX-01（Duplicate）、GX-04、GX-05、GX-10、GX-17〜23。TK-03、UJ-02、FL-02、FL-05、FL-07、FL-10 の Web 側、GX-16 は別の PR で直り、`29c631e2` で確かめた | PER-MODE（#945）の後 |

P3 のうち上にないもの（アクセシビリティ G04、補足の文字の統一、Interactive SVG、注釈、その他）は、0.14.0 の後に回した（Owner の指示、2026-10-07）。

## 確かめ直しの表

| ID | P | Lane/PR | dev 9b0c56bd | E1 76b45de9 | Note |
|---|---|---|---|---|---|
| TK-01 | P1 | A/G05 | 再現する (CLI: features cannot fit inside between 274.9px and 193.0px) | n/a | probes/tk01 |
| TK-02 | P2 | A/G05 | 再現する (Web RENDER_FAILED; raises at _resolve_circular_radial_layout L1943, _validate_same_side_order L1600/1606, _place_outside_auto L833) | n/a | test_web_error_producer_coverage.py L82 is in PY-B (v0.15.0 branch) |
| TK-03 | P2 | B/G06 | 症状が変わった (phantom series still created; now DEPTH_INVALID "Depth series 3. Supply the required value." since OV-89 #922) | 症状が変わった (same) | app-setup.js ensureDepthTrackEditableConfigCount: PR-1 plan §5 OV-109 fixes the getter write → hand case to Phase E |
| TK-04 | P2 | B/G07 | 再現する | - | circular/linear-track-slots.js (PR-0c) |
| TK-05 | P2 | - | 直っている (#805, 654cbff4) | - | close |
| TK-06 | P3 | B/G06 | 再現する | - | session-request.js + error-normalization.js (reserved) |
| TK-07 | P3 | B/G06 | 再現する | - | circular-track-slots.js (PR-0c) |
| TK-08 | P3 | B/G06 | 再現する | - | app-setup.js + index.html help text |
| TK-09 | P3 | B/G07 | 再現する | - | circular-track-slots.js (PR-0c) |
| TK-10 | P3 | B/G06 | 再現する | - | circular-track-slots.js (PR-0c) |
| TK-11 | P3 | A/G03 | 再現する | - | index.html CSS |
| TK-12 | P3 | A? (services/track-slot-validation.js free) | 再現する | - | parseOptionalCircularScalar; measure-editor field name |
| TK-13 | P3 | A? (validation free; feature-rendering PR-0c) | 再現する | - | Owner: match CLI |
| TK-14 | P3 | A/G03 | 再現する | - | index.html details open or docs |
| TK-15 | P3 | B/G07 | 再現する | - | display free; estimate in circular-track-slots.js (PR-0c) |
| TK-16 | P3 | A/G04 | 再現する (+ row buttons repeat names) | - | index.html |
| UI-01 | P2 | A/G01 | 再現する (plain <style>, 18 @apply) | n/a | orchestrator grep |
| UI-02 | P3 | A/G04 | 再現する (7 dialogs; 6 lack role/label; no focus-in; Escape closes popup behind) | n/a | rv/ui.md; openers in PR-0b files → shared mount/unmount component (new dialog-focus.js) + ui.js |
| UI-03 | P3 | A/G04 | 再現する | n/a | ui.js handleEscapeKey, app-setup closeFeaturePopup (reserved) |
| UI-04 | P3 | A/G04 | 症状が変わった (#809 named icon buttons; legend caption inputs, similarity search/sort, Help×66, per-row dup names, drawer tabs remain) | n/a | |
| UI-05 | P3 | A/G01 | 再現する (lang="ja") | n/a | |
| UI-06 | P3 | A/G01+G02 | 再現する (35 failing text nodes) | n/a | |
| UI-07 | P3 | B (run-analysis.js reserved) | 再現する (Pyodide starts before assertActiveModeInputs; 5.7 s warm) | n/a | cause in run-analysis.js (reserved) |
| UI-08 | P3 | A-out (aout-input) | 再現する | n/a | components.js FileUploader (free) |
| UI-09 | P3 | A/G04 | 再現する (19/28 disabled w/o reason) | n/a | index.html + app-setup.js predicate |
| UI-10 | P3 | 0.15.0 | handed to gbdraw-ff | n/a | handover/ui-10 |
| UI-11 | P3 | A/G20 | 再現する (997x817 in 904x406) | n/a | |
| UI-12 | P3 | A/G04 | 再現する (Generate = 35/384) | n/a | index.html only |
| UI-13 | P3 | A/G02/G03 | 症状が変わった (h-scroll only with sections open: 343/320; 70 px dead space; title 67/113 px) | n/a | |
| CI-01 | P3 | A-out (G11) | 再現する (Web RENDER_FAILED + false CLI text) | same (code) | crop_genbank.py check_start_end_coords (free); diagnostic= fixes Web too |
| CI-02 | P2 | B (G10) | 再現する | same (code) | app/run-info.js appendInputArgs records-table order (PER-MODE plan) |
| CI-03 | P2 | B (G10) | 再現する (Circular replace/Remove + Linear single→single) | same (code) | app-setup.js setCircularRecordPresentationSelector / setLinearSeqPrimaryFile, run-analysis runCircularRecordRefresh, watchers.js (reserved) |
| CI-04 | P2 | A-out (G12) | 再現する (message only) | same (code) | error_adapter.py _classify_native (free); error-normalization.js (E1) wording |
| CI-05 | P3 | A-out (G12) | 症状が変わった (5d3ade20 fixed Interval/FASTA; Query/Subject span rows still search frame) | same (code) | pairwise-match-popup.js buildMatchSpans (free); match-sequences.js (E1) |
| CI-06 | P2 | A-out (G12) | 再現する (uncoded throw → UNKNOWN) | same (code) | services/session-request.js buildComparisons L2127 (E1, PER-MODE) → exception hunk or hand over |
| CI-07 | P3 | A-out (G12) | 再現する a/b/d; c matches docs (Owner-delegated: explain) | same (code) | warnings not located yet; d also CLI |
| CI-08 | P3 | A-out (G11) | 再現する | same (code) | error_adapter.py: UnicodeError test before validation flag (free) |
| obs Reset alignment dialog | - | A-ui (G17) | 再現する | same | index.html ~L7491 fixed right-3 top-3 |
| obs drawer @1600 | - | A-ui (G17) check | 再現する (may be intended) | same | |
| obs popup vs footer | - | A-ui (G17) | 再現する (3 px) | same | CSS |
| obs region+rc order docs | - | A-out (docs) | 再現する (docs gap) | same | CLI_Reference.md / input-formats docs (free) |
| UJ-01 | P2 | B | 再現する | - | label-actions.js applyLabelVisibilityPreview (PER-MODE) |
| UJ-02 | P3 | E1 | 再現する | E1 で直る | hand case to gbdraw-ff (E1 acceptance) |
| UJ-03/FL-01 | P2 | - | 直っている (#857 8354d9f6; #870) | - | parity spec lacks popup Legend-name rename + Reset fill cases → Phase E LEGEND |
| UJ-04 | P3 | A-out (G18a) | 再現する (7/10 `out`) | - | tools/refresh_gallery_sessions.py _refresh_one_session; publication finalize keeps config (reserved) |
| UJ-05 | P3 | A-out (aout-ci) | 再現する | - | comparison-ui.js intentKeyForPlan (free) |
| UJ-06 | P3 | B (G21) | 再現する | - | Owner decision |
| UJ-07 | P3 | A-out (aout-input) | 再現する | - | record-discovery.js normalizeSequenceRecords (free); wording error-normalization.js (E1) |
| UJ-08 | P3 | B | 再現する | - | run-analysis.js runCircularRecordRefresh / watchers.js (reserved) |
| UJ-09 | P3 | B (G14) | 再現する | 再現する | app-setup importSession / config.js |
| UJ-10 | P3 | B (G14) | 再現する | 再現する | run-analysis.js cancel (reserved) + index.html text |
| obs search query across modes | - | B (G14) | 再現する | 再現する | feature-search/preview-actions.js (PER-MODE) |
| obs Linear Auto notice | - | B | 再現する | 再現する | app-setup.js linearLabelAutoDisclosure |
| FL-02 | P2 | Phase E LEGEND (OV-152) | 再現する | 再現する | retireSupersededLegendColors misses removed rules |
| FL-03 | P2 | - | 直っている (likely #892 OV-63) | 直っている | told Phase E to drop OV-153 |
| FL-04 | P3 | - | 直っている (#894 OV-62: Merge only same type) | 直っている | told Phase E to drop OV-151 |
| FL-05 | P3 | Phase E LEGEND (OV-154) | 再現する | 再現する | |
| FL-06 | P2 | - | 直っている (#892 OV-63) | 直っている | told Phase E to drop OV-155 |
| FL-07 | P2 | Phase E LEGEND (OV-127) | 再現する | 再現する | compactLegendEntries y=rectSize/2; moveLegendEntryToAnchor |
| FL-08 | P3 | paused (G15) | 再現する | n/a | standalone-interactivity-assets.js keydown (free) |
| FL-09 | P3 | paused (G15) | 再現する | n/a | .gfs-qualifier CSS |
| FL-10 | P3 | paused | 再現する | n/a | requests.py message unclassified; session-request safePrefix (reserved) |
| FL-11 | P3 | paused | 再現する (code) | n/a | rule-actions.js (PER-MODE), index.html wording |
| FL-12 | P3 | paused (explain) | 再現する | n/a | selector_values.py; Owner-delegated: explain |
| FL-13 | P3 | - | 直っている (#812 OV-27) | n/a | new bug GX-07 (tab-only line) |
| FL-14 | P3 | paused | 再現する | n/a | annotations.js (PER-MODE) + annotation-state.js |
| obs Interactive SVG popup title gene vs product | - | paused (G15) | 再現する | n/a | renderPopup display_label first |
| obs renamed Legend row moves to end | - | Phase E (told) | 再現する | 再現する | |
| obs Feature Edits TSV no fill/stroke | - | - | 仕様（help text documents it） | | close |
| obs scope choice closes popup | - | paused | 再現する | | |

## 新しく見つけた不具合（GX）

| ID | Found by | Summary | Where | 扱い |
|---|---|---|---|---|
| GX-01 | rv-tracks | Track-row controls without `:disabled="!semanticMutationAvailable"` (Circular Depth legend label, skew colors, annotation Lane gap/Padding, Duplicate; Depth legend title; Linear Height, Lane gap/Padding, Depth legend label, skew colors, Duplicate) accept input during a Session operation, the setter refuses, value reverts silently; canDuplicate*TrackSlot ignore running ops | index.html (bindings) + circular/linear-track-slots.js canDuplicate* (PR-0c) | 修正済み。`index.html` の bindings は #941、Duplicate の predicate は Lane B（`fix/web-gui-audit-b`）（b-tracks） |
| GX-02 | G20 agent | Zoom in passes 500% (template `zoom+=0.1` unclamped); in/out drift off 0.1 grid | index.html zoom group | 修正済み（#941） |
| GX-03 | rv-journeys | Gallery Sessions HmmtDNA_basic_circular and lambda_basic_linear carry maintainer-local paths in cliInvocation.args (`--gbk /mnt/c/...`, `-o /tmp/gbdraw_beginner/<id>`) | tools/refresh_gallery_sessions.py _canonicalize_recorded_cli_invocation | 修正済み（#939） |
| GX-04 | G05 agent | Explicit outer_gap_px on the row after a pinned inside row → "order cannot be honored" (now TRACK_LAYOUT row 2); same rows without the pin work. Cap after a pinned/group-boundary row uses only that row's inner gap, _validate_same_side_order requires the larger facing gap | radial_layout.py | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405）。隣り合う二行は向き合う gap の大きい方を保つ。基準 SVG は変わらない |
| GX-05 | G05 agent | Rows above a pinned row don't compress: pin 0.6013 (0.007 px above Auto) fails in Tuckin; "(auto)" note rounds, so typing the shown value fails where the note rounds up | radial_layout.py / track-slot-display | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405）。resolver は元から Auto と同じく縮める。GX-17 で例のコマンドが描ける |
| GX-06 | G05 agent | diagrams/circular/assemble.py:2281 "radial inner labels cannot fit the fixed circular geometry" raises without diagnostic= | assemble.py (R1/R2 branch file) | 修正済み（#939） |
| GX-07 | rv-features | Specific colors TSV with a tab-only line: CLI exits 1 "Missing values … at line 2" (bad line is 4; Web skips it). read_color_table reports table row+1, ignoring skipped lines | gbdraw/io/colors.py read_color_table | 修正済み（#939） |
| GX-08 | rv-features | Interactive SVG Details tab shows "Protein ID TRNF" for tRNA-Phe (displayProteinId falls back to gene) | standalone-interactivity-assets.js ~2541 | 修正済み（#939） |
| GX-09 | aout-ci agent | LOSAT FASTA path `applyLosatSequenceTransforms` (run-analysis.js) still throws the false "Start position … must be less than end position" for a region beyond the record (CI-01 fixed crop_genbank.py only) | app/run-analysis.js (PER-MODE) | 修正済み（#939） |
| GX-10 | aout-misc agent | getFeatureCaption fallback prints a 0-based start ("tRNA at 576..647" for 577..647); the runtime copy too; same string is a Legend caption key | services feature caption + runtime | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx10、b-labels）。caption は 1 始まり、strand は付けない（D-B01） |
| GX-11 | aout-misc agent | Interactive SVG search bar (390 px + 12 px margin) cut off on the right at 390 px viewport | standalone-interactivity-assets.js CSS | 修正済み（#939） |
| GX-12 | g02g03 agent | Settings controls outside track rows without busy `:disabled` (plot title, default font size, definition-line style, Linear scale font, Linear default definition/subtitle, record row number, annotation Duplicate set + stroke/fill colours, Preset scheme select, collinear evidence scope select) | index.html | 修正済み（#941） |
| GX-13 | g02g03 agent | Circular single-record group dimmed with opacity-60 dims its help tips (2.30:1) | index.html (G04 UI-09 lines) | 修正済み（#941） |
| GX-14 | g02g03 agent | Gallery tutorial caption "Choose Reset to Middle" (gallery/tutorials/HmmtDNA_ATskew.json) now names a button "Middle" under "Reset to preset" | gallery/** (PER-MODE) | 修正済み（#941） |
| GX-15 | g02g03 agent | Pairwise popup CSS still uses #94a3b8 text | index.html CSS (aui-misc region) | 修正済み（#941） |
| GX-16 | G04 agent | Color Change Scope Cancel (button/Escape) can take >5 s to close after a Session load: handler runs only after history.runUndoable returns from an intent capture | color-actions.js / history (PER-MODE) | 別の PR で修正済み（#945、OV-161）。`29c631e2` で case が通る |
| GX-17 | gx0405 agent | A row pinned at its own Auto radius fails when Auto had compressed it (MG1655 Tuckin default: Auto gc_content 66.09 px at 0.706555 R; r=0.706555 keeps 74.1 px → ticks cannot fit 316.5-386.1 px). A radius disables compression (`normalize_circular_track_slots` sets compress only when r is blank). Options: (A) keep + document; (B, recommended) a pinned numeric row with an Auto width compresses like Auto, centred on its radius, down to Auto's minimum (also makes the GX-05 command render) | radial_layout.py / tracks/circular.py normalize_circular_track_slots | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405）。Owner の決定（2026-10-08）の option B |
| GX-18 | gx0405 agent | Ticks row "(auto)" Radius note shows the tick band centre (`radiusFactor` = center_radius_px in assemble.py) but `r` pins the tick anchor: typing the shown value shifts ticks by half the tick length; fails in Tuckin and Middle on both records | assemble.py radiusFactor / track-slot-display.js | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-depth）。Web の注が tick の anchor を示す |
| GX-19 | gx0405 agent | MG1655: the CLI default stack (Use custom stack off) and the Web preset custom stack with nothing typed give different figures (default keeps numeric rows at preset anchors, Middle gc_content 0.75 R; custom stack packs them under ticks, 0.759 R). HmmtDNA identical. Pre-existing | presets.py / radial_layout | 修正済み。Lane B（`fix/web-gui-audit-b`）（gx0405、D-B02）。CLI の既定の stack を基準にした |
| GX-20 | b-depth agent (2026-10-08) | After a Session load, Legend Name Scope Cancel/Escape took > 5 s to close (OV-161 fixed only Color Change Scope); found when the 60 s wait in accessibility spec was removed | app-setup.js / feature-editor.js / color-actions.js | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-depth）。Cancel はすぐ閉じ、History の step を作らない |
| GX-21 | b-labels agent (2026-10-08) | Auto Reflow on: Undo of a label Off that the rerender left out does not bring the label back live (Generate draws it); same cause as UJ-01 | app/feature-editor/label-actions.js applyStoredVisibilityOverridesToSvg | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-labels） |
| GX-22 | b-session agent (2026-10-08) | Releasing a popup drag past its clamp (e.g. on the footer) closed the popup: the following click landed outside it | app/ui.js | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-session） |
| GX-23 | b-session agent (2026-10-08) | test_build_py_copies_offline_gui_assets (slow, main only) required web/js/app/record-discovery.js, moved to services/ by layering D (c8f00000) | tests/test_web_packaging.py | 修正済み。Lane B（`fix/web-gui-audit-b`）（b-session） |
