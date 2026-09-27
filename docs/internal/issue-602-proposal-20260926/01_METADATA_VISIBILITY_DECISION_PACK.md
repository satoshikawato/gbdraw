# Product Decision Pack — Linearメタデータ表示のfresh defaultとAutoの説明

状態: **SIGNED — B / AUTO-FRESH-RESET-WITH-DISCLOSURE**。署名者: **satoshikawato**、署名日: **2026-09-26**。

本Packは選択肢と判断根拠を保存する。採択はB、Aは非採択。正式authorityへの反映はS00の責務で、runtime baseへのmergeは未完了。[署名済み回答と機械表現](06_SIGNED_PRODUCT_DECISIONS.md)に採択全文を保存している。

## Identity

- Concern key: `linear.record-label-auto-visibility`
- Scenario revision: `2`
- Discovery lane: developer preflight
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Prepared for proposed branch/head: `fix/issue-602-linear-live-edit-20260926`、計画書をcommitし、実装を継続するwork branch。
- Prepared by: Codex、2026-09-26
- Related issue: [#602](https://github.com/satoshikawato/gbdraw/issues/602)
- 実装・証拠・全体依存: [00_MASTER_PLAN.md](00_MASTER_PLAN.md)

## Trigger / Authority search

BUG-11。共有行へ配置しただけでAccession/Lengthが図全体から消える。fresh defaultを変更するか、Autoを維持して理由を説明するかを決める。

| Source inspected | Result | Conflict or gap |
| --- | --- | --- |
| Product Impact map / BD | mapはrequest/Resultの境界を保護。tools/web-product-decisions.jsonのdecisionsは空 | 基準SHAでこの表示結果を選択するBDは存在しない |
| Static Product contract / Web reference | docs/REFERENCE/web-app.mdはfresh Auto・any shared row時のdiagram-wide Hiddenを明示。#562計画はINDEPENDENT-AUTO-SHOW-HIDEと当時の承認を記録 | fresh Showへの変更は既存記述を更新する必要がある。独立durable recordがなかったという記録も保持 |
| Domain / integrity | 生物学的identityと座標、ユーザー編集の整合性を維持 | A/BのUI選好を科学規則で選べない |
| Released compatibility | docs/SESSION_COMPATIBILITY.md、現在のreader/migratorとpositive fixtures | supported formatsだけを維持。新しい互換性約束はしない |
| Eligible exact-head decision | なし。PR-localはmappedなAFFORDANCE_PRESERVEDに限定 | 新しいdurable結果は正式base authorityへ反映する。候補だけでruntimeを自己承認しない |
| Current code/tests | 00のbase evidenceおよび本Packのobserved behavior | 現在のcode/testsだけではProduct authorityにならない |

Decision result: **SIGNED — B / AUTO-FRESH-RESET-WITH-DISCLOSURE**。
Procedural next action: **DURABLE_AUTHORITY_REQUIRED**。選択は完了している。S00で署名全文を既存authorityへ正確にserializeし、baseへのmerge後に依存runtimeを実装する。A/Bの比較資料や署名だけでGate/test failureを免除しない。

## User journey / Current observed behavior

Actor: Linear図を作る研究者。Journey ID: J602-METADATA。Entry: fresh起動/Reset、Arrange in rows、Record Labels。Preconditions: 同一/異なるソースから複数recordを選択。Goal: 配置を変えても識別情報の表示意図を理解する。Steps: 生成 → 共有行に変更 → 実効表示を確認 → 必要ならShow → Generate → Save/Load。Checkpoints: fresh/reset、draft共有行、Generate成功/失敗、Session読込、再生成、export。Failure/recovery: 失敗なら旧Resultと選択を保持して修正・再試行。Next actions: Generate、表示選択変更、共有行の解除、Undo。

現在のfresh Accession/Lengthは各Auto。実際にrenderするshared rowが一つでもあると各Autoは図全体でHiddenになる。独立Show/Hideは維持され、休眠layoutは無視される。Record Labelsのselect内にAuto · Hiddenはあるが、レイアウト操作場所に理由と修正先がない。既存pure testはこの挙動を確認している。

## Non-waivable constraints

- Architecture: SPA/no build step、単一canonical request / Worker / Result admission / History ownerを維持。新たなrender path、active authority store、互換readerを追加しない。
- Security/privacy: genome/Resultをbrowser内に保持し、同一originの資産と既存SVG sanitizationを維持。新たな外部依存を追加しない。
- Scientific correctness: sequence、annotation、source座標、feature identity、比較根拠を変更しない。表示文言は自動決定を科学的な確実性として説明しない。
- Persisted compatibility: 現在サポートするSession/requestを維持し、読み込みだけで保存Resultを生成し直さない。新schema/branch-only migrationを作らない。
- Performance/resource safety: Status/画面配置を理由にPython/LOSATを呼び出さない。追加のgenome複製・全量hash・SVG cloneを行わない。
- Required evidence: 最終runtimeでfocused testsとbrowser acceptanceを実施。Product選択は既存Gate/test failureの免除ではない。
- 明示Autoの意味はA/Bとも従来どおり。手入力SubtitleとRepliconは別責務。Showを有効にした代表入力でテキストの欠落/切断/誤座標を認めない。

## Choice A / Choice B / Comparison matrix

ファイルや実装ownerではなく製品結果を比較する。表の各行は各Choiceの本文であり、共通制約も両方に適用する。

| Dimension | Choice A（非採択） | Choice B（採択） |
| --- | --- | --- |
| Stable outcome ID | A / SHOW-FRESH-RESET-WITH-AUTO-DISCLOSURE | B / AUTO-FRESH-RESET-WITH-DISCLOSURE |
| Complete normative outcome | fresh/resetのAccessionとLengthは独立Show。Autoは保持。Auto+shared rowのfieldについてLinear Layoutに理由、diagram-wide範囲、次回Generateの効果、Record Labelsへの変更先を表示 | fresh/resetのAutoを保持。同じ場所・内容のAuto非表示説明と変更先を表示 |
| Preserved effects | Auto/Show/Hide独立、明示値、保存Result、共有行解除でAuto再表示 | Aと同じ |
| Added effects | fresh共有行でメタデータを保持。Auto使用時の説明と修正先 | Auto使用時の説明と修正先 |
| Lost / retired effects | fresh/resetがAutoを選ぶ既定値だけ退役 | 既定値/表示意味の退役なし |
| Entry / discoverability / accessibility | LayoutとLabelsから理解・keyboard移動。polite statusはfieldと理由を特定 | Aと同じ |
| Immediate feedback | 共有行のdraft変更で次回Generateの表示を更新。現在Resultと混同しない | Aと同じ。freshでAuto非表示を経験しやすい |
| Canonical state update | fresh factoryのみShow。説明はresolverから派生。変更先ボタンはfocus/scrollだけ | fresh factoryはAuto。説明とボタンはAと同じ |
| Undo / Redo | 説明にHistoryを作らない。表示変更と共有行変更は既存履歴どおり | Aと同じ |
| Session round trip / compatibility | 明示値を維持。supported old omission/boolean migrationは従来の意味。読込Result不変 | Aと同じ |
| Regeneration | 選択値と実際のshared rowからbooleanを一度投影。SaveのAutoをShowへ変えない | Aと同じ |
| Export / artifact | 現在Resultを出力。fresh/resetから新生成した図はShow | 現在Resultを出力。fresh shared rowから新生成した図はAutoでHidden |
| Validation / error | 未知のmodeは既存validation error。失敗は旧Resultを維持 | Aと同じ |
| Failure / recovery | 選択を保持し、Show変更またはlayout修正からGenerate再試行 | Aと同じ |
| Scientific-output consequence | 科学データは不変。識別/長さテキストのfresh表示だけ増える | 科学データ不変。fresh Autoの表示結果も不変 |
| Cache / provenance | 検索・render cache規則は不変。表示選択は既存configへ保存 | Aと同じ |
| Performance | 説明は純粋なresolver。Showによる既存text/geometry増加を代表入力で確認 | 説明のみ追加。通常render量は現行 |
| Architecture | 既存default factoryと同じvisibility resolverを使用。第二のlayout policyを作らない | Aと同じ |
| Evidence available / missing | pure resolverとSession evidenceあり。fresh Showの長文/shared-row browser geometryと新説明が未検証 | Auto baselineあり。新説明のbrowser evidenceが未検証 |
| Residual risk | fresh Showで図の文字量が増え、非常に密な配置は手動調整が必要になりうる。黙って隠さない | 説明があっても共有行でメタデータが消えること自体は残る |
| Route | DURABLE_AUTHORITY_REQUIRED | DURABLE_AUTHORITY_REQUIRED（説明を新しい保証として記録） |
| Next action | 非採択。default Showへの変更は実装しない | 既定Autoと説明のauthorityをmerge後、V01–V03を実装/確認 |

## Evidence-first option / Engineering recommendation

必要な比較入力は単一/共有/混在行、layout disabled、長いaccession、短いrecord、大きなフォント、Auto/Show/Hideの混合。出力は選択値、実効boolean、実SVGテキスト/幾何、save/load/regeneration。測定は表示品質を検証するが、識別情報と図の簡潔さの優先順位を自動決定しない。evidence-only作業でdefault/authorityを変更しない。全状況の無衝突保証はない。

採択結果: **B / AUTO-FRESH-RESET-WITH-DISCLOSURE**。実装はこの結果の全寄与を満たす。未選択の結果は実装対象に含めない。

## Product Decision Owner response

採択Bの全文は署名済み。Aは比較用に保存した非採択回答。採択されたRationale、維持/退役範囲、riskは改変しない。正式authorityのserializationとbaseへのmergeはS00で行う。署名は外部操作の許可を含まない。

### A — 非採択回答

```text
PRODUCT_DECISION
Concern: linear.record-label-auto-visibility
Scenario revision: 2
Choice: A / SHOW-FRESH-RESET-WITH-AUTO-DISCLOSURE
Rationale: 共有行への配置は識別情報を隠す意思表示ではないため、初めて作る図ではAccessionとLengthを表示し、自動非表示は利用者が選べる形で残す。
Must preserve: 独立したAuto/Show/Hide、明示Autoのdiagram-wide非表示条件と共有行解除時の再表示、disabled layoutの休眠行除外、既存Sessionの明示値と対応済み省略意味、読込時の保存Result、手入力Subtitle、Repliconの独立性、Undo/Redo、GenerateとExportの区別。Auto非表示時はLayoutにも理由・対象field・図全体の範囲・次回Generateの効果・Record Labelsへの変更先を表示する。
May retire: Web fresh/resetがAccessionとLengthの初期値にAutoを選ぶことだけ。明示Autoの意味は退役しない。
Accepted residual risk: fresh Showは文字量を増やし、非常に密な図は利用者の間隔やフォント調整が必要になりうる。代表的な長文/短record/共有行で欠落・切断・誤座標を認めず、自動的な非表示や省略で回避しない。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名なし。非採択。

### B — 署名済み回答

```text
PRODUCT_DECISION
Concern: linear.record-label-auto-visibility
Scenario revision: 2
Choice: B / AUTO-FRESH-RESET-WITH-DISCLOSURE
Rationale: 共有行の図では簡潔な既定表示を維持し、情報が非表示になる理由とShowへの変更先を配置操作の場所で明示する。
Must preserve: fresh/resetの独立Auto、Show/Hideの明示値、diagram-wide Auto解決、休眠行除外、既存Sessionと保存Result、Undo/Redo、GenerateとExportの区別。Auto非表示時はLayoutにも理由・対象field・図全体の範囲・次回Generateの効果・Record Labelsへの変更先を表示する。
May retire: なし。
Accepted residual risk: 共有行でAccession/Lengthが非表示になる結果自体は残る。説明を見落とす可能性があるため、layout操作場所とLabelsの両方で実効値と変更先を示す。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。
