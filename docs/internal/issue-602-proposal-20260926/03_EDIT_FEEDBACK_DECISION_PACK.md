# Product Decision Pack — 編集・Generate・保存・Exportの適用状態を示す

状態: **SIGNED — A / DERIVED-APPLICATION-STATUS**。署名者: **satoshikawato**、署名日: **2026-09-26**。

本Packは選択肢と判断根拠を保存する。採択はA、Bは非採択。正式authorityへの反映はS00の責務で、runtime baseへのmergeは未完了。[署名済み回答と機械表現](06_SIGNED_PRODUCT_DECISIONS.md)に採択全文を保存している。

## Identity

- Concern key: `web.edit-application-feedback`
- Scenario revision: `1`
- Discovery lane: developer preflight
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Prepared for proposed branch/head: `fix/issue-602-linear-live-edit-20260926`、計画書をcommitし、実装を継続するwork branch。
- Prepared by: Codex、2026-09-26
- Related issue: [#602](https://github.com/satoshikawato/gbdraw/issues/602)
- 実装・証拠・全体依存: [00_MASTER_PLAN.md](00_MASTER_PLAN.md)

## Trigger / Authority search

BUG-13。生成設定draft、即時編集、Apply前draftが混在し、利用者がGenerateの必要性と再配置を誤解する。編集をすべて自動適用するのではなく、どの事実を画面上で保証するかを選ぶ。

| Source inspected | Result | Conflict or gap |
| --- | --- | --- |
| Product Impact map / active BD | mapはcanonical request/Resultの継続を保護。BD storeは空 | 基準SHAでこの新しい表示結果を決めるBDはない |
| Static Web runtime contract | gbdraw/web/CLAUDE.mdは右の単一編集の即時canonical commitと必要時rerenderを規定 | その契約を変更しない。全生成設定のglobal Pending表示は未規定 |
| Static Product contract | OIC-013 / PD-OI-016、PD-OI-032、PD-OI-034等は失敗時のResult保持、pending設定と対象限定操作の分離、review draftを要求 | 新しいStatusで既存の継続を変えない。全手動配置をGenerate後にも維持する権限はない |
| Domain / integrity | 正確な適用状態、source/feature identity、atomic commitmentを維持 | A/Bの表示方式を一意に選ばない |
| Released compatibility | docs/SESSION_COMPATIBILITY.mdと既存restore fixtures | supported writer/readerと保存Result/draft分離を維持 |
| Eligible exact-head decision | なし。PR-local routeはmapped AFFORDANCE_PRESERVEDに限定 | candidate authorityでruntimeを自己承認しない |
| Current code/tests | 各Packのobserved evidenceと00のbase検証 | 現在の挙動がそのまま正解になるわけではない |

Decision result: **SIGNED — A / DERIVED-APPLICATION-STATUS**。
Procedural next action: **DURABLE_AUTHORITY_REQUIRED**。選択は完了している。S00で署名全文を既存authorityへ正確にserializeし、baseへのmerge後に依存runtimeを実装する。A/Bの比較資料や署名だけでGate/test failureを免除しない。

## User journey / Current observed behavior

Actor: 両modeで図を整える利用者。Journey ID: J602-APPLICATION。Entry: sidebar、Editor、review、Generate、Save/Export。Preconditions: 既存Resultがあり、一部設定を変更。Goal: 何が今の図へ反映され、何が次回Generateで変わるか理解する。Steps: scale変更 → live色/label編集 → 状態確認 → SaveまたはExport → Generate → Undo。Checkpoints: draft編集、live適用開始/成功/失敗、review Apply/Cancel、Generate開始/終了、Session round trip、履歴復帰。Failure path: 最終成功Resultを保持し、pending設定から修正・再試行。Next actions: Generate、設定を戻す、Undo、Save、現在ResultのExport。

recordDisplayControlsとPaletteには部分的なPending説明があるが、全生成設定の状態表示はない。Palette Instant Previewは左でも即時反映する。右のAlignment reviewはApply前draft。Generate候補は既存editor mutationを継承し、成功時にzoomをresetする。幾何の再計算や一部編集の継承範囲をglobal表示が説明していない。既存session-draft-authority testはdraftとResultの分離を確認する。

## Non-waivable constraints

- Architecture: SPA/no build step、単一canonical request / Worker / Result admission / History ownerを維持。新たなrender path、active authority store、互換readerを追加しない。
- Security/privacy: genome/Resultをbrowser内に保持し、同一originの資産と既存SVG sanitizationを維持。新たな外部依存を追加しない。
- Scientific correctness: sequence、annotation、source座標、feature identity、比較根拠を変更しない。表示文言は自動決定を科学的な確実性として説明しない。
- Persisted compatibility: 現在サポートするSession/requestを維持し、読み込みだけで保存Resultを生成し直さない。新schema/branch-only migrationを作らない。
- Performance/resource safety: Status/画面配置を理由にPython/LOSATを呼び出さない。追加のgenome複製・全量hash・SVG cloneを行わない。
- Required evidence: 最終runtimeでfocused testsとbrowser acceptanceを実施。Product選択は既存Gate/test failureの免除ではない。
- 右の単一編集をGenerate/close/Saveまで延期しない。pending生成設定をtarget-only操作へ混ぜない。
- 対応済みoverrideの保持、失敗/Cancel/stale時のResultとHistory、Undo/Redoを維持。比較根拠不明やinvalidをAppliedと表示しない。

## Choice A / Choice B / Comparison matrix

表は各Choiceの完全な本文であり、共通制約も両方に適用する。実装ファイルの選択ではない。

| Dimension | Choice A（採択） | Choice B（非採択） |
| --- | --- | --- |
| Stable outcome ID | A / DERIVED-APPLICATION-STATUS | B / STATIC-APPLICATION-CONTRACTS |
| Complete normative outcome | 操作ごとのLive edit / Applies on Generate / Apply requiredを表示。Result付近に有効生成intentのPending/invalid/unknownを正確に表示。Liveのapplying/errorとPendingを独立に表示。Generate再配置、Save/Exportの意味も表示 | 同じ操作分類とGenerate/Save/Export説明を固定表示。global Pending/Appliedの判定は追加せず、現在の局所的Pendingと既存busy/errorを維持 |
| Preserved effects | 現在の適用タイミング、canonical即時commit、review draft、atomic Generate、失敗時保持 | Aと同じ |
| Added effects | 局所説明と正確なglobal Pending。値を戻した/Undoした/Loadした後も事実に合わせる | 局所説明のみ。事実に基づかないglobal clean表示を作らない |
| Lost / retired effects | 既存適用/継続の退役なし。重複表示ロジックは同じ意味で統合 | 既存適用/継続の退役なし |
| Entry / accessibility | 操作名付近とResult/Generate付近に短い説明。role=statusで重要遷移をpoliteに伝え、全keystrokeを読上げない | 同じ場所に固定説明。既存busy/errorを読上げ |
| Immediate feedback | 生成intent差とlive状態から導出。位置では分類しない | 固定説明を常時読める。個別のPending以外は差分を表示しない |
| Canonical state update | Statusは観測のみ。Generate成功/restore/履歴でartifactの比較基準を更新。live適用は実際にcommitした部分だけ反映 | Statusはstateを変更しない。global比較基準を新設しない |
| Undo / Redo | artifact/draft復帰に追従。Statusだけの履歴はゼロ | 既存局所表示が復帰、固定説明は不変 |
| Session round trip | 保存Result/requestとactive draftの差を再構成。Status enumを保存しない。settings-onlyは未生成を示す | 保存Resultとdraftを維持。固定説明だけで差分を推定しない |
| Regeneration | Pendingが次の成功Generateで適用。幾何再計算とzoom resetを事前説明。対応済みcanonical editsを維持 | Aと同じ。成功した設定との差はglobal表示しない |
| Export / artifact | 現在Resultを出力することを明示。pending/live errorを勝手に適用しない | Aと同じ |
| Validation / error | invalidと比較不能unknownをcleanと区別し、既存入力errorへ誘導 | 既存validation/errorを維持。固定文だけでAppliedと断定しない |
| Failure / recovery | 失敗・Cancel・staleで比較基準を進めず、旧Result/History/Pendingを保持。retry/Undo可能 | 既存failure/retry/local Pendingを維持 |
| Scientific output | 選択・描画・座標不変。Statusのための自動Generateをしない | Aと同じ |
| Cache / provenance | statusをcache key/科学的provenanceへ追加しない。既存resource identityのみ比較 | cache/provenance不変 |
| Performance | 既存projectionと共用の軽いintent比較。Worker/file bytes/hash/SVG clone追加ゼロ | 固定文のみ。処理量の増加は小さい |
| Compatibility | 新schemaなし。保存Result/draftをLoad時に混ぜない | Aと同じ |
| Architecture | projection/比較はcanonical owner。表示は既存ownerの事実へ依存。dirty flagの第二ownerを作らない | 新comparison owner不要。既存partial Pendingを維持 |
| Evidence available / missing | draft authorityとlive契約あり。全field coverage、同名source差替え、mixed Pending+Live、load/undo/errorのbrowser testsが必要 | 既存局所表示あり。分類の例外と説明のbrowser testsが必要 |
| Residual risk | 比較対象の漏れやlive commit基準の誤更新が誤表示を生む。invalid/unknownを明示し、生成/履歴/読込の一致テストを必須とする | 利用者が個々の変更を覚えておく必要があり、Generateが必要かの判断負担が残る |
| Route | DURABLE_AUTHORITY_REQUIRED | DURABLE_AUTHORITY_REQUIRED（新しい操作説明の保証） |
| Next action | authority merge後、projectionの出力一致を先に確認してE01–E05を実装 | authority merge後、操作分類/固定説明/例外を確認。全設定差分表示は実装しない |

## Evidence-first option / Engineering recommendation

代表scenarioはscale/crop/slots、metadata modes、file交換、inactive mode/slot、Instant Preview ON/OFF、右live編集、review draft、pendingを含むSession、Undo/Redo、失敗/Cancel/stale。差分比較はResource bytesを読まず、同名別Fileも識別することを測定する。既存保存Resultのrequestが生成意図の比較基準であり、ロード時のactive draftを基準にしてpendingを消してはいけない。evidence-onlyではcommit/default/authorityを変更しない。固定説明と動的表示のどちらを保証するかは製品判断。

採択結果: **A / DERIVED-APPLICATION-STATUS**。実装はこの結果の全寄与を満たす。未選択の結果は実装対象に含めない。

## Product Decision Owner response

採択Aの全文は署名済み。Bは比較用に保存した非採択回答。採択されたRationale、維持/退役範囲、riskは改変しない。正式authorityのserializationとbaseへのmergeはS00で行う。署名は外部操作の許可を含まない。

### A — 署名済み回答

```text
PRODUCT_DECISION
Concern: web.edit-application-feedback
Scenario revision: 1
Choice: A / DERIVED-APPLICATION-STATUS
Rationale: 利用者が編集後の図と次回Generateの変更を区別できるよう、操作の適用タイミングと生成設定の未適用状態を事実から表示する。
Must preserve: 操作単位のLive edit、Applies on Generate、Apply requiredを区別し、左Palette Instant Previewと右Alignment reviewを例外なく正しく分類する。canonical即時commitと必要時自動rerender、reviewのlocal draft、target-only操作とpending設定の分離、対応済みoverride継承、atomic Generate、失敗/Cancel/stale時の旧ResultとHistory、Undo/Redo、SessionのResult/draft分離、Exportの現在Result出力を維持する。Pendingとlive applying/errorを独立に示し、invalid/unknownをAppliedとしない。Generateによる配置再計算とzoom reset、Save/Exportの意味を事前に説明する。
May retire: 製品の適用タイミング・保存・復旧・編集機能は退役しない。全DOM座標や手動位置を再生成後にも無条件に保持する保証は新設しない。
Accepted residual risk: 生成intent比較の漏れや比較基準の誤更新は誤表示を生みうる。生成・live commit・履歴・Sessionの一致テストを必須とし、根拠不足はunknownとして表示する。Status目的のWorker呼出、genome byte読取り/hash、SVG/checkpoint cloneは受け入れない。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### B — 非採択回答

```text
PRODUCT_DECISION
Concern: web.edit-application-feedback
Scenario revision: 1
Choice: B / STATIC-APPLICATION-CONTRACTS
Rationale: 適用境界の説明を最小の変更で明確にし、全生成設定の差分比較を新たな保証にしない。
Must preserve: 操作単位のLive edit、Applies on Generate、Apply required、Palette Instant Previewとreviewの例外、canonical即時編集、自動rerender、review draft、target-only/pending分離、対応済みoverride継承、Generate失敗/Cancel/staleの旧ResultとHistory、Undo/Redo、Session/Export意味を維持する。Generate再配置とzoom reset、Save/Exportの意味を固定文で説明し、既存局所Pending/busy/errorを維持する。global Appliedと断定する表示は作らない。
May retire: なし。
Accepted residual risk: 全設定を集約したPending表示は提供せず、利用者が変更を覚えてGenerateの必要性を判断する負担が残る。説明は各操作とResult/Generate付近で常時読めることを必須とする。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名なし。非採択。
