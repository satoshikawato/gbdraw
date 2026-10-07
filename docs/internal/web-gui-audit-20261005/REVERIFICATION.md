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
| A-out | TK-01（P1）、TK-02、外側の行の順序（Owner の決定）、入力と比較のエラー（CI-04、CI-06 ほか）、Gallery の refresh tool（UJ-04、GX-03）、この記録 | Phase E の PER-MODE の前 |
| A-ui | UI-01 の刷新（提案）、UI-05、UI-06、Fit ボタン（UI-11、GX-02） | Phase E の PER-MODE の前 |
| B | Lane B の P2（TK-03、TK-04、UJ-01、CI-02、CI-03）と記録 | PER-MODE の後 |

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

| ID | Found by | Summary | Where | Plan |
|---|---|---|---|---|
| GX-02 | G20 agent | Zoom in passes 500% (template `zoom+=0.1` unclamped); in/out drift off 0.1 grid | index.html zoom group | fixed in G20 (A-ui) |
| GX-03 | rv-journeys | Gallery Sessions HmmtDNA_basic_circular and lambda_basic_linear carry maintainer-local paths in cliInvocation.args (`--gbk /mnt/c/...`, `-o /tmp/gbdraw_beginner/<id>`) | tools/refresh_gallery_sessions.py _canonicalize_recorded_cli_invocation | A-out with G18a (tool fix; PR-1 regenerates) |
| GX-04 | G05 agent | Explicit outer_gap_px on the row after a pinned inside row → "order cannot be honored" (now TRACK_LAYOUT row 2); same rows without the pin work. Cap after a pinned/group-boundary row uses only that row's inner gap, _validate_same_side_order requires the larger facing gap | radial_layout.py | not fixed (changing gap rule moves rendering stacks); report |
| GX-05 | G05 agent | Rows above a pinned row don't compress: pin 0.6013 (0.007 px above Auto) fails in Tuckin; "(auto)" note rounds, so typing the shown value fails where the note rounds up | radial_layout.py / track-slot-display | report (residual) |
| GX-06 | G05 agent | diagrams/circular/assemble.py:2281 "radial inner labels cannot fit the fixed circular geometry" raises without diagnostic= | assemble.py (R1/R2 branch file) | report; small, candidate for A-out if free |
| GX-07 | rv-features | Specific colors TSV with a tab-only line: CLI exits 1 "Missing values … at line 2" (bad line is 4; Web skips it). read_color_table reports table row+1, ignoring skipped lines | gbdraw/io/colors.py read_color_table | paused (P3) |
| GX-08 | rv-features | Interactive SVG Details tab shows "Protein ID TRNF" for tRNA-Phe (displayProteinId falls back to gene) | standalone-interactivity-assets.js ~2541 | paused (G15) |
| GX-01 | rv-tracks | Track-row controls without `:disabled="!semanticMutationAvailable"` (Circular Depth legend label, skew colors, annotation Lane gap/Padding, Duplicate; Depth legend title; Linear Height, Lane gap/Padding, Depth legend label, skew colors, Duplicate) accept input during a Session operation, the setter refuses, value reverts silently; canDuplicate*TrackSlot ignore running ops | index.html (bindings) + circular/linear-track-slots.js canDuplicate* (PR-0c) | index.html part in G03; canDuplicate part Lane B G07 |
