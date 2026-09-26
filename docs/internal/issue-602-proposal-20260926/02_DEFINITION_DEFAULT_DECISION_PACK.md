# Product Decision Pack — Web fresh/resetのDefinition列をLock ONにするか

状態: **SIGNED — A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A**。署名者: **satoshikawato**、署名日: **2026-09-26**。

本Packは選択肢と判断根拠を保存する。採択はA、Bは非採択。正式authorityへの反映はS00の責務で、runtime baseへのmergeは未完了。[署名済み回答と機械表現](06_SIGNED_PRODUCT_DECISIONS.md)に採択全文を保存している。

## Identity

- Concern key: `linear.definition-display`
- Scenario revision: `2`
- Discovery lane: developer preflight
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Prepared for proposed branch/head: `fix/issue-602-linear-live-edit-20260926`、計画書をcommitし、実装を継続するwork branch。
- Prepared by: Codex、2026-09-26
- Related issue: [#602](https://github.com/satoshikawato/gbdraw/issues/602)
- 実装・証拠・全体依存: [00_MASTER_PLAN.md](00_MASTER_PLAN.md)

## Trigger / Authority search

BUG-12。行をalign/offsetしたときfresh Lock OFFではDefinitionが段状に追従する。Webの初期値変更と、既に修正されたLock ON座標バグを区別して選択する。

| Source inspected | Result | Conflict or gap |
| --- | --- | --- |
| Product Impact map / BD | mapはrequest/Resultの境界を保護。tools/web-product-decisions.jsonのdecisionsは空 | 基準SHAでこの表示結果を選択するBDは存在しない |
| Active static Product contract | PD-OI-024 scenario 1: D1-A / D2-P / D3-A | D1-AはLock=false既定を明示。AはD1のWeb fresh/resetだけを改訂。D2/D3の独立寄与を保持する必要がある |
| Domain / integrity | 生物学的identityと座標、ユーザー編集の整合性を維持 | A/BのUI選好を科学規則で選べない |
| Released compatibility | docs/SESSION_COMPATIBILITY.md、現在のreader/migratorとpositive fixtures | supported formatsだけを維持。新しい互換性約束はしない |
| Eligible exact-head decision | なし。PR-localはmappedなAFFORDANCE_PRESERVEDに限定 | この新規/改訂のdurable結果を候補だけでruntimeを自己承認しない |
| Current code/tests | 00のbase evidenceおよび本Packのobserved behavior | 現在のcode/testsだけではProduct authorityにならない |

Decision result: **SIGNED — A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A**。
Procedural next action: **DURABLE_AUTHORITY_REQUIRED**。選択は完了している。S00で署名全文を既存authorityへ正確にserializeし、baseへのmerge後に依存runtimeを実装する。A/Bの比較資料や署名だけでGate/test failureを免除しない。

## User journey / Current observed behavior

Actor: publication用のLinear比較図を作る研究者。Journey ID: J602-DEFINITION。Entry: fresh起動/ResetとLinear Layout。Preconditions: unequal translationを含む単一/共有/混在行。Goal: Definitionの読み比べに適した列を得る。Steps: recordを配置 → Align/offset → Generate → Lock切替 → Save/Load → Generate。Checkpoints: default、入力、最終座標、読込、再生成、download。Failure/recovery: 失敗で以前の図を保持。次はLock切替、位置変更、Generate、Undoが可能。

fresh Lock=false。D1-Aどおり各rowの移動にDefinitionが追従する。Lock=trueの共通左列は最新baseでfinal translations後の共通originを用い、非browserのDefinition tests 45件が成功している。初期値変更はこの修正の代替ではない。

## Non-waivable constraints

- Architecture: SPA/no build step、単一canonical request / Worker / Result admission / History ownerを維持。新たなrender path、active authority store、互換readerを追加しない。
- Security/privacy: genome/Resultをbrowser内に保持し、同一originの資産と既存SVG sanitizationを維持。新たな外部依存を追加しない。
- Scientific correctness: sequence、annotation、source座標、feature identity、比較根拠を変更しない。表示文言は自動決定を科学的な確実性として説明しない。
- Persisted compatibility: 現在サポートするSession/requestを維持し、読み込みだけで保存Resultを生成し直さない。新schema/branch-only migrationを作らない。
- Performance/resource safety: Status/画面配置を理由にPython/LOSATを呼び出さない。追加のgenome複製・全量hash・SVG cloneを行わない。
- Required evidence: 最終runtimeでfocused testsとbrowser acceptanceを実施。Product選択は既存Gate/test failureの免除ではない。
- PD-OI-024の手入力/保存Subtitle、継承規則、row共通とrecord固有の区別、Replicon/Organelle制御、既存text_anchor受入範囲は両案で維持する。
- CLI/Pythonの省略defaultは変更しない。Lock=falseの明示選択は共通幅中央＋row追従を維持する。

## Choice A / Choice B / Comparison matrix

ファイルや実装ownerではなく製品結果を比較する。表の各行は各Choiceの本文であり、共通制約も両方に適用する。

| Dimension | Choice A（採択） | Choice B（非採択） |
| --- | --- | --- |
| Stable outcome ID | A / D1-B-WEB-LOCKED-FRESH | B / D1-A-WITH-PROMINENT-GUIDANCE |
| Complete normative outcome | Web fresh/resetのLockはtrue。最終表示で共通左端とconfigured gapを使う。OFFの選択は従来のrow追従。CLI/Python省略は従来どおり。D2-P/D3-Aを維持 | Web fresh/resetはfalse。Linear Layoutでrow追従とLock ONの共通列を常時説明。描画とCLI/Pythonは現行。D2-P/D3-Aを維持 |
| Preserved effects | Lock ON列、OFF配置、全Subtitle/Replicon寄与、既存Session | Aと同じ |
| Added effects | freshでalign/offset後もDefinition共通左列。Layoutで選択の効果を常時説明 | Layoutで選択の効果を常時説明 |
| Lost / retired effects | Web fresh/resetがLock OFFを選ぶ既定だけ | 既定/描画意味の退役なし |
| Entry / discoverability / accessibility | Linear Layoutのcheckboxと常時短い説明。tooltipだけに依存しない | Aと同じ |
| Immediate feedback | ON/OFFの効果とApplies on Generateを表示 | Aと同じ |
| Canonical state update | createDefaultFormの既存boolean。選択はcanonical requestへ一度投影 | factoryのfalseを維持。説明はstateを変更しない |
| Undo / Redo | 既存設定変更/Generate履歴。説明だけにHistoryを作らない | Aと同じ |
| Session / compatibility | 明示false/trueを保持。supported旧omissionは従来false。ResultはLoadで不変 | Aと同じ |
| Regeneration | 保存選択を使用。新factory値で上書きしない | Aと同じ |
| Export / artifact | 現在Result。fresh再生成のDefinition配置だけ変わる | 現在Result。fresh配置は現行 |
| Validation / error | 既存boolean/text_anchor受入範囲。invalidを新defaultへ修復しない | Aと同じ |
| Failure / recovery | 以前のResultと選択を維持、設定修正/Generate/Undo | Aと同じ |
| Scientific output | sequence/feature/座標は不変。Definition表示位置だけ変更 | 科学内容と配置不変 |
| Cache / provenance | 既存requestのLock値。追加cache policyなし | Aと同じ |
| Performance | 既存locked placementのみ。追加render/Worker初期化なし | 常時説明のみ |
| Architecture | 単一Web factory、既存共通Python placementを再利用。新Web位置engineなし | 既存factory/placementを維持 |
| Evidence available / missing | 既存coordinate correctionと45 testsあり。fresh/reset、legacy omission、browser gap/Sessionの変更後確認が必要 | 既存OFF baselineあり。常時説明のkeyboard/browser確認が必要 |
| Residual risk | fresh diagramのDefinitionがrowに追従しなくなり、以前の既定外観から変わる。OFFは選べる | fresh diagramの段状配置は残り、整った列には利用者のON操作が必要 |
| Route | DURABLE_AUTHORITY_REQUIRED | DURABLE_AUTHORITY_REQUIRED（常時説明の保証を含む改訂） |
| Next action | PD-OI-024のD1を限定改訂。D2-P/D3-Aを保持してauthority merge後D01–D03 | D1-Aを維持し常時説明を追加。authority merge後D02–D03と説明を確認 |

## Evidence-first option / Engineering recommendation

最終SVG上の共通左端とgapを、負/正の不均等translation、mixed shared rows、center/Similarity alignmentで測定する。保存ON/OFFとsupported旧omission、Subtitle/Replicon/text_anchorを含める。既存Lock=true修正のbrowser証拠とWeb fresh/resetの変更後証拠は別である。evidenceは整列の正否を判定するが初期値の製品選好を選ばない。evidence-onlyで既定値やauthorityは変えない。

採択結果: **A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A**。実装はこの結果の全寄与を満たす。未選択の結果は実装対象に含めない。

## Product Decision Owner response

採択Aの全文は署名済み。Bは比較用に保存した非採択回答。採択されたRationale、維持/退役範囲、riskは改変しない。正式authorityのserializationとbaseへのmergeはS00で行う。署名は外部操作の許可を含まない。

### A — 署名済み回答

```text
PRODUCT_DECISION
Concern: linear.definition-display
Scenario revision: 2
Choice: A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A
Rationale: Webで新しく作るLinear比較図ではDefinitionを共通左列にそろえ、行のalignやoffset後も名前を比較しやすくする。
Must preserve: Lock=trueの共通左端とconfigured gap、明示Lock=falseの共通幅中央とrow追従、単一/共有/混在行、既存text_anchorの受入範囲、保存Sessionの明示値と対応済み旧省略意味、読込時の保存Result、CLI/Python省略default。PD-OI-024のD2-Pの保存/手入力Subtitleと継承・ラベル区別、およびD3-AのReplicon/Organelle選択順・独立制御・既定falseをすべて維持する。Linear LayoutでON/OFFの違いとGenerate適用を常時説明する。
May retire: D1-Aのうち、Web fresh/resetがLock=falseを初期値として選ぶ部分だけ。OFFの明示操作とCLI/Python既定は退役しない。
Accepted residual risk: 新しいWeb図のDefinition外観が従来のfresh図と変わり、Definitionがrowに追従しなくなる。利用者はOFFを選べ、既存Sessionの値と保存Resultは勝手に変更しない。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### B — 非採択回答

```text
PRODUCT_DECISION
Concern: linear.definition-display
Scenario revision: 2
Choice: B / D1-A-WITH-PROMINENT-GUIDANCE; retain D2-P and D3-A
Rationale: Definitionがrowへ追従する既定配置を維持し、共通列が必要な利用者にはLinear LayoutでLockの効果を明確に伝える。
Must preserve: D1-AのWeb/CLI/Pythonの現行既定、Lock ONの共通左端とgap、OFFの共通幅中央とrow追従、全行構成、text_anchor受入範囲、Sessionの値と保存Result。D2-Pの保存/手入力Subtitleと継承・ラベル区別、D3-Aの命名選択順・独立制御・既定falseを維持する。LayoutにON/OFFの違いとGenerate適用を常時表示する。
May retire: なし。
Accepted residual risk: freshのaligned/offset図ではDefinitionの段状配置が残る。共通列には利用者がONにする必要があり、常時説明を見落とす可能性がある。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名なし。非採択。
