# Product Decision Pack — 狭いPreviewで即時Editorと図を同時に使う

状態: **SIGNED — A / DOCKED-COMPACT-EDITOR**。署名者: **satoshikawato**、署名日: **2026-09-26**。

本Packは選択肢と判断根拠を保存する。採択はA、Bは非採択。正式authorityへの反映はS00の責務で、runtime baseへのmergeは未完了。[署名済み回答と機械表現](06_SIGNED_PRODUCT_DECISIONS.md)に採択全文を保存している。

## Identity

- Concern key: `web.editor.compact-presentation`
- Scenario revision: `1`
- Discovery lane: developer preflight
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Prepared for proposed branch/head: `fix/issue-602-linear-live-edit-20260926`、計画書をcommitし、実装を継続するwork branch。
- Prepared by: Codex、2026-09-26
- Related issue: [#602](https://github.com/satoshikawato/gbdraw/issues/602)
- 実装・証拠・全体依存: [00_MASTER_PLAN.md](00_MASTER_PLAN.md)

## Trigger / Authority search

BUG-16のEditor。幅390pxでdrawerが図とtoolbarを実質的に遮蔽する。即時編集の意味は維持し、図の空間を再配分するか、部分遮蔽のsheetにするかを選ぶ。Alignment reviewのApply/Cancelや優先順位はPack 05の責務。

| Source inspected | Result | Conflict or gap |
| --- | --- | --- |
| Product Impact map / active BD | mapはcanonical request/Resultの継続を保護。BD storeは空 | 基準SHAでこの新しい表示結果を決めるBDはない |
| Static Web runtime contract | right-drawer.jsが唯一のvisibility/tab owner。単一編集はlive commit。閉じることでcommitしない | 狭い画面の同時Preview空間の保証は未規定 |
| Existing responsive evidence | 既存40rem container CSSとright-drawer mobile tests | 幅390pxのテストはtoolbar遮蔽を許容。テスト成功はProduct選択ではない |
| Domain / integrity | 正確な適用状態、source/feature identity、atomic commitmentを維持 | A/Bの表示方式を一意に選ばない |
| Released compatibility | docs/SESSION_COMPATIBILITY.mdと既存restore fixtures | supported writer/readerと保存Result/draft分離を維持 |
| Eligible exact-head decision | なし。PR-local routeはmapped AFFORDANCE_PRESERVEDに限定 | candidate authorityでruntimeを自己承認しない |
| Current code/tests | 各Packのobserved evidenceと00のbase検証 | 現在の挙動がそのまま正解になるわけではない |

Decision result: **SIGNED — A / DOCKED-COMPACT-EDITOR**。
Procedural next action: **DURABLE_AUTHORITY_REQUIRED**。選択は完了している。S00で署名全文を既存authorityへ正確にserializeし、baseへのmerge後に依存runtimeを実装する。A/Bの比較資料や署名だけでGate/test failureを免除しない。

## User journey / Current observed behavior

Actor: スマートフォン/狭いPreviewで図を調整する利用者。Journey ID: J602-EDITOR-COMPACT。Entry: Editor toggle。Preconditions: CircularまたはLinearのResult。Goal: 編集の変化を図で確認する。Steps: Editorを開く → 色/ラベル/visibilityを編集 → 図をzoom/pan → tab変更 → Close/Escape → 幅変更。Checkpoints: open/close、編集、camera操作、keyboard、resize、Result置換。Failure/recovery: live rerender失敗は既存契約どおり、直に適用された編集を保持してerrorを示す。Next actions: 修正、Undo、panel scroll、camera、Close。

390×844の実測ではPreview幅374px、drawer幅334px、高さ473.5px。図の大部分が見えず、toolbar6操作の中心hitが遮られる。既存テスト1件は成功し、この遮蔽を許容している。画像を目視済み。wide drawerの可用性/close/tab invariantは既存unit testsでも確認済み。

## Non-waivable constraints

- Architecture: SPA/no build step、単一canonical request / Worker / Result admission / History ownerを維持。新たなrender path、active authority store、互換readerを追加しない。
- Security/privacy: genome/Resultをbrowser内に保持し、同一originの資産と既存SVG sanitizationを維持。新たな外部依存を追加しない。
- Scientific correctness: sequence、annotation、source座標、feature identity、比較根拠を変更しない。表示文言は自動決定を科学的な確実性として説明しない。
- Persisted compatibility: 現在サポートするSession/requestを維持し、読み込みだけで保存Resultを生成し直さない。新schema/branch-only migrationを作らない。
- Performance/resource safety: Status/画面配置を理由にPython/LOSATを呼び出さない。追加のgenome複製・全量hash・SVG cloneを行わない。
- Required evidence: 最終runtimeでfocused testsとbrowser acceptanceを実施。Product選択は既存Gate/test failureの免除ではない。
- 同じSVG・同じEditor content・同じselection/open stateを用いる。mini-previewやmobile版の第二editorは作らない。
- Close/Escapeはvisibilityだけを変更。tab可用性の同期reconcile、live History/Export/Session、失敗時の復旧を維持する。

## Choice A / Choice B / Comparison matrix

表は各Choiceの完全な本文であり、共通制約も両方に適用する。実装ファイルの選択ではない。

| Dimension | Choice A（採択） | Choice B（非採択） |
| --- | --- | --- |
| Stable outcome ID | A / DOCKED-COMPACT-EDITOR | B / PARTIAL-OVERLAY-COMPACT-EDITOR |
| Complete normative outcome | 狭いPreviewで図を上、Editorを下に配置して図の領域を確保。wideは既存side drawer。Editorのlistは独立scroll、Close/headerとcanvas toolbarを到達可能にする | 狭いPreviewでEditorを下側の高さ制限付きsheetとして重ねる。上の図領域とtoolbarは常に操作可能。図はpanで必要部分を可視領域へ移す。wideは既存side drawer |
| Preserved effects | 全tab、単一live edit、close/tab selection、camera、History | Aと同じ |
| Added effects | 編集中にも図の利用可能幅全体と可視高さを確保。重ならない操作領域 | 編集中の上側可視図と操作領域。図の下側は遮蔽しうる |
| Lost / retired effects | 狭いPreviewで横からslideして全面高さを覆う配置だけ | 同じ旧配置だけ。下側の部分遮蔽は明示して残る |
| Entry / accessibility | 同じtoggle/tab/Close/Escape。mobileで折返しとscroll。screen readerの順序は図→Editor | 同じ操作。sheetの上側可視図とtoolbarへkeyboardで移動可能 |
| Immediate feedback | live変更が同じ図へ反映。自動rerender中/失敗は既存意味 | Aと同じ。ただし対象が下側ならpanが必要 |
| Canonical state update | open/tabは既存ownerのみ。responsive layoutはstateの複製を持たない | Aと同じ |
| Undo / Redo | 編集の既存History。open/resize/scrollはartifact履歴を増やさない | Aと同じ |
| Session round trip | 新しいpanel位置/heightを保存しない。既存Result/editor stateを維持 | Aと同じ |
| Regeneration | 既存Generate意味。表示切替はrender/LOSATを起動しない | Aと同じ |
| Export / artifact | 同じResultのみ。panel配置やcameraを科学的図のgeometryへ混ぜない | Aと同じ |
| Validation / error | 各editorの既存validation/error。図を見えるようにしても未適用をAppliedとしない | Aと同じ |
| Failure / recovery | 直に適用された編集を保持、既存rerender errorとretry/Undo。panelはClose可能 | Aと同じ |
| Scientific output | sequence/coordinates/Result SVGをresponsive配置で変更しない | Aと同じ |
| Cache / provenance | 不変。viewportはPython/LOSAT cacheの入力にしない | Aと同じ |
| Performance | CSS/container配置中心。geometry変更に限る既存camera調整。SVG cloneゼロ | 高さ制限中心。SVG cloneゼロ |
| Compatibility | 新schemaなし。desktop操作を維持 | Aと同じ |
| Architecture | 同一markup/contentと既存right-drawer owner。新visibility flagを追加しない | Aと同じ |
| Evidence available / missing | 390pxの失敗品質を実測済み。新配置のM01–M05とactual live updateが必要 | 現行overlayの測定あり。sheetと可視canvas hit/scrollの証拠が必要 |
| Residual risk | 図の縦の領域とEditor listの両方が短くなり、list scrollが増える | 図の下部をsheetが隠すため、対象をpanする追加操作が必要 |
| Route | DURABLE_AUTHORITY_REQUIRED | DURABLE_AUTHORITY_REQUIRED |
| Next action | authority merge後、同一drawerのresponsive配置とM01–M05を確認 | authority merge後、制限sheet/toolbarの配置と可視図操作を確認 |

## Evidence-first option / Engineering recommendation

390×844/740では可視canvas幅をPreview利用可能幅、高さ200px以上とし、header/Generate barを差し引く。390×500、320px、landscape、200% zoom、soft keyboardでは全操作へ到達/scroll回復を確認。矩形の収まりだけでなくelementFromPointと実際のzoom/pan/live編集を確認する。baseline画像は00のevidence。追加の比較prototypeは製品選好の前提ではない。evidence-onlyでruntime/default/authorityを変更しない。mini-previewを使って達成したとは扱わない。

採択結果: **A / DOCKED-COMPACT-EDITOR**。実装はこの結果の全寄与を満たす。未選択の結果は実装対象に含めない。

## Product Decision Owner response

採択Aの全文は署名済み。Bは比較用に保存した非採択回答。採択されたRationale、維持/退役範囲、riskは改変しない。正式authorityのserializationとbaseへのmergeはS00で行う。署名は外部操作の許可を含まない。

### A — 署名済み回答

```text
PRODUCT_DECISION
Concern: web.editor.compact-presentation
Scenario revision: 1
Choice: A / DOCKED-COMPACT-EDITOR
Rationale: 狭いPreviewでも即時編集の変化を図で確認できるよう、図とEditorを上下の領域へ配置する。
Must preserve: 同じSVGとEditor、全tabと同期可用性、canonical live commitと必要時rerender、既存History/Session/Export、camera操作、keyboard、Close/Escapeのvisibility-only意味、選択tab、Result置換/失敗復旧。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、Editor内容を独立scrollさせ、Close/headerとtoolbarを操作可能にする。短いviewport/soft keyboardでは全操作へscrollで到達できる。wideのside drawerを維持する。
May retire: 狭いPreviewでEditorが横から全面高さを覆う表示配置だけ。編集機能や保存意味は退役しない。
Accepted residual risk: 上下分割で図とEditor listの縦領域が短くなり、list scrollが増える。実操作のpointer/keyboard/browser検証を必須とし、複製Preview・SVG clone・第二editorによる回避は受け入れない。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### B — 非採択回答

```text
PRODUCT_DECISION
Concern: web.editor.compact-presentation
Scenario revision: 1
Choice: B / PARTIAL-OVERLAY-COMPACT-EDITOR
Rationale: 図のviewport寸法を保ちながら狭い画面の遮蔽を制限し、編集パネルの上側で図を確認できるようにする。
Must preserve: 同じSVG/Editor、全tab、live commit/rerender、History/Session/Export、Close/Escape/tab意味と失敗復旧。狭い画面では高さ制限sheetの上に利用可能幅全体と390×844/740で高さ200px以上の操作可能canvasを残し、toolbarを覆わない。内容scroll、keyboard、pan/zoom、短いviewport/soft keyboardで全操作への到達、wide side drawerを維持する。
May retire: 狭いPreviewでEditorが横から全面高さを覆う表示配置だけ。
Accepted residual risk: sheetは図の下側を覆い、編集対象を見るためpanが必要になりうる。可視領域のpointer操作、toolbar非遮蔽、Closeへの到達を必須とする。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名なし。非採択。
