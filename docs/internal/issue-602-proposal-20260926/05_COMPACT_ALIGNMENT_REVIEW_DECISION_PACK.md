# Product Decision Pack — 狭いPreviewでAlignment reviewとcanvasを同時に使う

状態: **SIGNED — A / DOCKED-COMPACT-ALIGNMENT-REVIEW**。署名者: **satoshikawato**、署名日: **2026-09-26**。

本Packは選択肢と判断根拠を保存する。採択はA、Bは非採択。正式authorityへの反映はS00の責務で、runtime baseへのmergeは未完了。[署名済み回答と機械表現](06_SIGNED_PRODUCT_DECISIONS.md)に採択全文を保存している。

## Identity

- Concern key: `web.similarity-alignment.review-presentation`
- Scenario revision: `1`
- Discovery lane: developer preflight
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Prepared for proposed branch/head: `fix/issue-602-linear-live-edit-20260926`、計画書をcommitし、実装を継続するwork branch。
- Prepared by: Codex、2026-09-26
- Related issue: [#602](https://github.com/satoshikawato/gbdraw/issues/602)
- 実装・証拠・全体依存: [00_MASTER_PLAN.md](00_MASTER_PLAN.md)

## Trigger / Authority search

BUG-16のAlignment review。現在の非モーダルpaletteが390pxの図を覆いうる。既存のcanvas操作・Select/Skip・Apply/Cancel・retryは既に必須。その意味を変えず、狭い画面の固定dockと浮動panelのどちらを選ぶかを決める。

| Source inspected | Result | Conflict or gap |
| --- | --- | --- |
| Product Impact map / active BD | mapはcanonical request/Resultの継続を保護。BD storeは空 | 基準SHAでこの新しい表示結果を決めるBDはない |
| Static Product contract | PD-OI-031/PD-OI-034はresolved自動Apply、明示review、local draft、1回のbatch validation、canvas interaction、failure/retryを要求。PD-OI-027/28/29/30/35がtransform/plan/Historyを要求 | これらを退役しない。狭い画面のdrag退役とEditorを閉じる優先順位は新判断 |
| Static Web runtime contract / UI | reviewはrole=dialogだがaria-modalなし。paletteのdrag/clamp/focusはapp-setup、選択/適用はsimilarity-alignment owner | PD-OI-035 scenario 2は390px遮蔽と同時canvas操作退役を明示的に許容。採択Aの導入ではこの旧mobile例外を限定supersedeする。非モーダル、identity/keyboard/Skip/非描画候補/focus/overlay非保存は維持 |
| Domain / integrity | 正確な適用状態、source/feature identity、atomic commitmentを維持 | A/Bの表示方式を一意に選ばない |
| Released compatibility | docs/SESSION_COMPATIBILITY.mdと既存restore fixtures | supported writer/readerと保存Result/draft分離を維持 |
| Eligible exact-head decision | なし。PR-local routeはmapped AFFORDANCE_PRESERVEDに限定 | candidate authorityでruntimeを自己承認しない |
| Current code/tests | 各Packのobserved evidenceと00のbase検証 | 現在の挙動がそのまま正解になるわけではない |

Decision result: **SIGNED — A / DOCKED-COMPACT-ALIGNMENT-REVIEW**。
Procedural next action: **DURABLE_AUTHORITY_REQUIRED**。選択は完了している。S00で署名全文を既存authorityへ正確にserializeし、baseへのmerge後に依存runtimeを実装する。A/Bの比較資料や署名だけでGate/test failureを免除しない。

## User journey / Current observed behavior

Actor: スマートフォンでSimilarity alignmentを確認する研究者。Journey ID: J602-ALIGNMENT-COMPACT。Entry: 明示Review、ambiguous応答、automatic render失敗のretry。Preconditions: 比較済みResultと対象候補。Goal: canvas上で候補を見比べ、Select/SkipとMatch reference directionを選び安全に適用する。Steps: review開始 → 対象行/候補を選択 → canvasをzoom/pan → draft編集 → ApplyまたはCancel → failure correction/retry。Checkpoints: open、候補marker、local編集、resize、Apply開始/成功/失敗、Cancel、drawer reopen。Next actions: 候補変更、Skip、Apply、Cancel、error修正/retry、終了後Editor再open。

現行paletteはTeleport body、最大26rem、狭い場合100vw−1.5rem、最大100dvh−1.5rem。ドラッグしてwindow境界へclampする。canvas overlayと選択controllerは別で、局所draftの編集はWorkerを呼ばない。既存browser testsは390pxでpanel幅/Apply到達等を検証しているが、計画基準SHAでのreviewの実遮蔽率は追加計測していない。Editorの実測をreviewの測定として転用しない。

## Non-waivable constraints

- Architecture: SPA/no build step、単一canonical request / Worker / Result admission / History ownerを維持。新たなrender path、active authority store、互換readerを追加しない。
- Security/privacy: genome/Resultをbrowser内に保持し、同一originの資産と既存SVG sanitizationを維持。新たな外部依存を追加しない。
- Scientific correctness: sequence、annotation、source座標、feature identity、比較根拠を変更しない。表示文言は自動決定を科学的な確実性として説明しない。
- Persisted compatibility: 現在サポートするSession/requestを維持し、読み込みだけで保存Resultを生成し直さない。新schema/branch-only migrationを作らない。
- Performance/resource safety: Status/画面配置を理由にPython/LOSATを呼び出さない。追加のgenome複製・全量hash・SVG cloneを行わない。
- Required evidence: 最終runtimeでfocused testsとbrowser acceptanceを実施。Product選択は既存Gate/test failureの免除ではない。
- PD-OI-031/034と現行transform/plan/reset/historyの全寄与を維持。1つのMatch reference directionと各targetの結果方向を維持。resolved自動Applyは通常どおり、ambiguous/明示reviewはApply前draft。
- canvasは非モーダルで操作できる。Select/Skipの変更はlocal、Applyの共有Python検証とatomic Result/Historyを維持。失敗/Cancel/stale/supersededは以前のResult/orientation/Historyを維持し、retry draftを失わない。

## Choice A / Choice B / Comparison matrix

表は各Choiceの完全な本文であり、共通制約も両方に適用する。実装ファイルの選択ではない。

| Dimension | Choice A（採択） | Choice B（非採択） |
| --- | --- | --- |
| Stable outcome ID | A / DOCKED-COMPACT-ALIGNMENT-REVIEW | B / FLOATING-COMPACT-ALIGNMENT-REVIEW |
| Complete normative outcome | 狭いPreviewでは図を上、同じreviewを下段に固定。local draftのまま、候補listをscroll、Apply/Cancelを到達可能にする。wideはdrag可能なpalette。狭い場合review開始時にEditorをownerで閉じ、終了まで理由付きでopenをdisable。tabは保持、終了後は明示reopen | 狭いPreviewでは同じreviewを高さ制限のfloating下側panelとし、上側canvas/toolbarを操作可能にする。ドラッグを維持し、canvasの必要部分はpanで確認。Editorの優先処理はAと同じ。wideは現行palette |
| Preserved effects | resolved自動適用、ambiguous/明示review、全候補/方向/Skip、canvas、validation/retry | Aと同じ |
| Added effects | reviewと可視canvasの非重複領域。panel位置を調整せず候補比較できる | 遮蔽を制限したcanvas領域。panel位置の手動調整を維持 |
| Lost / retired effects | 狭いPreviewでreviewを自由dragする操作だけ。review中のEditor同時openを理由付きで制限 | 狭いreviewの高さ占有を制限。review中のEditor同時openを理由付きで制限。dragは退役しない |
| Entry / accessibility | 既存entryとrole=dialog。非モーダル、focus移動/復帰、Escape/Cancel。canvasへkeyboardで移動可能 | Aと同じ。dragに依存せず初期位置とApply/Cancelが到達可能 |
| Immediate feedback | local choiceが候補markerと結果方向へ反映。canonical ResultはApply成功まで変えない | Aと同じ |
| Canonical state update | draft/controller不変。画面配置はcanonical planを変更しない。drawerはclose ownerのみ | Aと同じ |
| Undo / Redo | 成功Applyは既存1操作。panel resize/drag/closeはartifact履歴を増やさない | Aと同じ |
| Session round trip | plan/Resultの既存意味。transient draft/位置を新しく保存しない | Aと同じ |
| Regeneration | 既存typed plan/path。viewportやpanel位置で候補resolve/renderを追加しない | Aと同じ |
| Export / artifact | accepted Result/planのみ。候補markerとpanelをexportしない | Aと同じ |
| Validation / error | local編集はWorkerなし。Applyで既存batch validation。理由とerrorを狭い画面でも表示 | Aと同じ |
| Failure / recovery | Apply失敗はdraftとerror保持して修正/retry。Cancel/staleは旧artifact保持。focus復帰、Editor明示reopen | Aと同じ |
| Scientific output | 正確なreference/biological identity、方向policy、Python eligibilityを維持。viewportで候補を選ばない | Aと同じ |
| Cache / provenance | 検索cache、plan provenance、source座標不変 | Aと同じ |
| Performance | 配置切替でWorker呼出ゼロ。候補再計算/第二SVGゼロ | Aと同じ |
| Compatibility | 新schemaなし。現在のplan/review/Reset意味を維持 | Aと同じ |
| Architecture | review presentationと選択/適用を分離。既存controller一つ、旧narrow clampを除去 | 既存controller一つ、height制限と可視viewportへのclampを使用 |
| Evidence available / missing | 既存controller/retryと390px幅テストあり。review遮蔽baseline、dock/resize/focus/draft継続の証拠が必要 | 同じbaselineが必要。floating下側可視canvas/drag/Apply hitの証拠が必要 |
| Residual risk | review listのscrollが増え、狭い幅では自由位置調整ができない。review中Editorは一旦閉じる | 下側図の部分遮蔽と位置調整/panが必要。review中Editorは一旦閉じる |
| Route | DURABLE_AUTHORITY_REQUIRED | DURABLE_AUTHORITY_REQUIRED |
| Next action | authority merge後、既存controllerを維持してreview表示とM01–M05、failure/retryを確認 | authority merge後、bounded floating表示と同じ受入を確認 |

## Evidence-first option / Engineering recommendation

実装前にresolved明示review、ambiguous複数target、automatic failure retryを390×740/844、500px高、landscapeで実測する。候補Select/Skipと方向、canvasのelementFromPoint/pan/zoom、Apply/Cancel/フォーカス復帰、開いたままresize、Editor優先close/tab保持を確認。390×740/844のcanvas可視幅/高さ200px条件はEditor Packと同じ測定式を使うが、draft/Apply証拠は別である。evidence-onlyでresolver/default/authorityを変更しない。A/Bはdrag継続と図の空間配分の製品判断で、計測だけで選ばない。

採択結果: **A / DOCKED-COMPACT-ALIGNMENT-REVIEW**。実装はこの結果の全寄与を満たす。未選択の結果は実装対象に含めない。

## Product Decision Owner response

採択Aの全文は署名済み。Bは比較用に保存した非採択回答。採択されたRationale、維持/退役範囲、riskは改変しない。正式authorityのserializationとbaseへのmergeはS00で行う。署名は外部操作の許可を含まない。

### A — 署名済み回答

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.review-presentation
Scenario revision: 1
Choice: A / DOCKED-COMPACT-ALIGNMENT-REVIEW
Rationale: 狭いPreviewでもalignment候補をcanvasで確認できるよう、reviewを図の下段に固定し、候補比較へ操作を集中させる。
Must preserve: PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、1つのMatch reference directionと各targetの結果方向、canvas操作、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。狭いreview開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisable、終了後は明示reopen可能。wideのdragと非モーダルcanvasを維持する。
May retire: 狭いPreviewでreviewを自由にdragする操作、およびreview中にEditorを同時openする継続だけ。候補や方向の選択、Apply前draft、failure/retryは退役しない。
Accepted residual risk: 狭いreviewではlist scrollが増え、自由に位置を動かせなくなる。開始時Editorは閉じるがtabは保持し、終了後再openできる。位置変更で候補draftやResultを変えないことをbrowserで確認する。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### B — 非採択回答

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.review-presentation
Scenario revision: 1
Choice: B / FLOATING-COMPACT-ALIGNMENT-REVIEW
Rationale: 狭い画面でも自由にreview位置を調整できる継続を残し、panelの高さと初期位置を制限してcanvasを使えるようにする。
Must preserve: resolvedの通常自動Apply、ambiguous/明示reviewのlocal draft、独立Select/Skip、正確なreference/候補根拠、1つのMatch reference directionと結果方向、非モーダルcanvas、local no-Worker編集、共有Python batch validation、atomic Result/History、failure draft/error/retry、Cancel/stale/supersededの旧artifact、transform/plan/reset/history、Session/regeneration/Export、focus復帰。狭いfloating panelの上に利用可能幅全体と390×844/740で高さ200px以上のcanvas/toolbarを残し、dragなしでもApply/Cancelへ到達可能にする。Editorはreview開始でowner経由close、tab保持、review中理由付きopen disable、終了後明示reopenとする。
May retire: 狭いreviewの全面に近い高さ占有、およびreview中のEditor同時open継続だけ。drag、候補/方向選択、Apply前draft、failure/retryは退役しない。
Accepted residual risk: panelが図の下側を覆い、候補を見るためpanや位置調整が必要になりうる。Editorは一旦閉じる。可視canvas/toolbar、drag-independent Apply/Cancel、resize時のdraft保持をbrowserで確認する。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名なし。非採択。
