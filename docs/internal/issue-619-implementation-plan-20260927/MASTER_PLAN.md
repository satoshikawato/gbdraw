# Circular Custom Track Slots: 数値入力と単位選択の実装計画

## 1. 目的・対象・状態

Issue: [#619](https://github.com/satoshikawato/gbdraw/issues/619)。Circular の Custom Track Slots で保存済み Width/Radius が `[object Object]` と表示される不具合を修正し、各 field を **数値入力＋明示的な単位選択**にする。

作業ブランチは `fix/issue-619-circular-track-measure-inputs`。公開準備時に取得した最新 `origin/dev` は `88028fd242d263f0fe86aaf9da57b8dc9eb082f6`。この SHA まで fast-forward してから計画書をコミットする。この計画書の公開は runtime 実装の完了や Product authority の承認を意味しない。

製品目標は数値と単位の分離。単位の表示 domain、切替動作、Auto 時の単位保持には独立した判断がある。各 [Decision Pack](DECISION_PACKS/README.md) の推奨 A は署名前の候補である。後述の実装詳細は **3 Pack の A が正式に選択された場合の具体案**。他の選択では該当部分だけを更新する。未署名の推奨案を実装者が承認済みとみなしてはならない。

対象は Circular の width/radius。Linear、純 pixel の inner/outer gap、Python renderer の geometry、CLI の scalar grammar は変更しない。

S00 の独立した検証成果は [SESSION_RESULTS/S00.md](SESSION_RESULTS/S00.md) に記録する。
既存 typed numeric-text draft は受理されたが、3 Pack の署名と
対象 authority の `origin/dev` 統合は未成立。settings-only の実 Save も既存 writer/validator 不整合で失敗。S01 の開始を許可しない。

## 2. 問題と現行契約

基準の tobacco chloroplast Gallery Session は Session 44 / request schema 8。`config.adv.circular_track_slots` に width/radius の `{value, unit}` object が保存されている。

| field | 保存値 | 現行入力 text | 目標表示（Pack 01-A） |
| --- | --- | --- | --- |
| plastome_regions width | `{value:20, unit:"px"}` | `[object Object]` | 数値 `20`、単位 `px` |
| plastome_regions radius | `{value:0.65, unit:"factor"}` | `[object Object]` | 数値 `0.65`、単位 `×R` |
| gc_content width | `{value:0.08, unit:"factor"}` | `[object Object]` | 数値 `0.08`、単位 `×R` |
| gc_content radius | `{value:0.56, unit:"factor"}` | `[object Object]` | 数値 `0.56`、単位 `×R` |
| Features width/radius | `null` | 空欄＋auto geometry | 空欄、Auto geometry は単位付きで明示 |

R は基準円半径。factor はその倍率で、無次元。width と radius の両方で px/factor が有効。現行 scalar は `1.5` を factor、`1.5px` を px、`150%` を factor 1.5 と扱う。この意味は 0.13.0 にも存在する。

request-derived projection は `projectCanonicalCircularMeasure()` で object を text 化する。一方、`restoreCurrentWriterActiveConfig()` は保存 config の編集 draft を優先する。`normalizeCircularTrackSlot()` がその object を維持し、native `v-model.trim` が text input に直接渡すため表示が壊れる。保存 draft の復元優先順を逆転してはならない。

[baseline browser observations](base-observations.json) の取得 SHA は `f5f86634459e0dcd46c1a452e9219fbba635d429`。公開基点との間で対象 Web runtime、fixture、Product/architecture policy に差分がないことを確認し、この証拠を再利用する。同 JSON は1440px/390pxの四つのobject text、focus/blurでのslot/Result/Undo不変を記録する。新しいcontrolsの成功証拠ではない。[scalar boundary observations](scalar-boundary-observations.json) はtyped numeric textの既存normalization/active-config/slot validation/payload受理を確認する。actual Save/Load、native reader、HistoryはS00で確認が必要。

再現: `node docs/internal/issue-619-implementation-plan-20260927/probe_scalar_boundary.mjs`。browser baselineは固定baseを空のsnapshotに`git archive`し、そのsnapshot内でwheelをprepareした後、`observe_base.py --source-root <snapshot> --output <new-json-path>`を実行する。snapshotに新UIを適用してbaselineのexpected failureを修正後gateとして使わない。source hashes・browser/wheel versionはJSONに記録する。

### 固定する不変条件

- committed request/Result と次回 Generate 用 draft を分ける。Load は保存 preview を保持する。
- 読み取り・panel open・focus/blur だけで scalar、Result、History を変更しない。
- canonical scalar は `{value: finite positive number, unit: "px" | "factor"}` または null。request へ raw text や `%` unit を渡さない。
- px と factor を黙示変換しない。Generate は既存の唯一の typed request/Worker/Result admission を使う。
- 不正・未完成の編集は見えるまま保持する。失敗時の直前 Result/request を維持し、訂正後に再 Generate できる。
- Session の受理 domain、disabled row、inactive stack、settings-only、config のない CLI-origin Session を保つ。不正 Session の受理は拡張しない。
- 同値の保存値は同じ canonical value/unit と geometry。入力値を表示目的で丸めない。

## 3. Product preflight と承認境界

表示バグだけの修正は既存の editable scalar 契約への適合である。数値＋単位 control は新しい編集 UX のため、developer preflight は `PRODUCT_DECISION_REQUIRED`。

| Pack | 独立した責任 | 推奨候補 |
| --- | --- | --- |
| [01](DECISION_PACKS/01_INPUT_REPRESENTATION.md) | 数値・単位の表示／入力、既存 percent と suffix 入力の扱い | A: selectors は px/×R。既存 `%` は factor 表示、単位付き入力は共有 adapter で取り込む |
| [02](DECISION_PACKS/02_UNIT_CHANGE.md) | 単位変更が表すユーザー操作 | A: 数値を保持し、選んだ単位で再解釈。図への反映は Generate |
| [03](DECISION_PACKS/03_AUTO_UNIT_LIFECYCLE.md) | Auto の単位選択の寿命と永続化 | A: Auto の unit は次の入力用の transient preference。空欄は常に null |

各 Pack に authority search、選択肢、effects、比較、route、全文記入済みの署名用 `PRODUCT_DECISION` を置く。署名は concern ごとに取得する。rationale/retirement/risk は署名前の案であり、実装者が人間の選択を推論して埋めてはならない。

S00 は明示 receipt を確認し、既存 static Product Contract にのみ承認内容をシリアライズする。新しい BD 番号を捏造せず、独自 JSON registry を作らない。machine representation を人間へ提示する。**authority-only の内容が origin/dev に入った後**に S01 以降を始める。candidate authority で同じ candidate runtime を承認しない。

authority-only PR の公開・merge は repository policy と明示 authorization に従う。計画書ブランチへの commit/push の許可だけを、別ブランチや merge の許可へ広げない。署名や authority 統合が未完なら、S00 は独立した証拠・文書作業を完了して commit/push し、依存 runtime は変更しない。

## 4. 実装アーキテクチャ（推奨 A の組合せ）

### 4.1 一つの scalar owner、一つの editor draft

`app/track-slot-validation.js` の既存 Circular scalar 検証を拡張し、canonical pair/null を返す DOM-free parse boundary にする。`circularScalarPayload()` と validation が同じ関数を呼び、重複した解釈を同じ変更で除去する。既存の有効な値・unit・percent・正数条件、typed object の扱いを保つ。Boolean/nonfinite/unsupported unit の不一致を受理拡張で解消しない。

既存 draft と typed object を読む／数値＋unit を書く小さな codec は `app/circular-track-slots/measure-editor.js` に置く。codec は既存 scalar owner を使い、canonical schema を決めない。Session projection の既存 formatter は同 codec の適切な読み取り helper に委譲し、二つの独立 formatter を残さない。歴史 reader の動作も検証する。

数値欄は `type="text" inputmode="decimal"` を使う。`type="number"` は不正文字列を空欄へ消すため、この draft 保持要件を満たすとは限らない。native Vue trim/composition behavior を維持する。

editor model の概念:

```text
既存 slot.width / slot.radius（保存値または編集 draft）
  → codec の read-only view: valueText, selectedUnit, Auto/error state
  → 数値 input＋unit select
  → owner action が元の scalar draft を1回更新
  → 同じ validation / request / Session / History
```

非空の入力は元の slot scalar を唯一の正本とする。現在の Web scalar object を使う場合、編集中は `value` に text を保持できるか S00 で実測する。例: `{value:"1.", unit:"px"}`。既存 request encoder は canonical number へ投影し、不正 text は拒否する。typed request の `value` domain を text へ広げない。

新しい Web draft shape、別の scalar array、accepted-value mirror、別 request pipeline を追加しない。既存 active-config/import/native reader がこの表現を受けない場合は、必要な representation の evidence を提示して境界を解決する。黙って schema を変えたり不正 draft を null 化したりしない。

### 4.2 専用 input component と wiring

小さな `CircularMeasureInput` を `app/circular-track-slots/measure-input.js` に置く。二つの field が同じ component/codec を使う。`app.js` で登録し、`index.html` は native object binding と重複 suffix をこの component に置き換える。top-level `createCircularTrackSlotEditor()` は既存 owner actions の入口を維持する。

component は modelValue を読み、更新を emit する。props を直接変更しない。親/editor action が元の width/radius を更新する。数値と unit の両方に slot ID＋field を含む accessible name を付ける。px／×R と R の意味を help で説明する。unit dropdown の値を見れば `1.5` の解釈が一意に分かる。

Pack 01-A では既存 typed px を数値20＋px、factor を数値0.65＋×R、`65%` を数値0.65＋×R と投影する。Load 自体では元の slot を書き換えない。完全で有効な `20px`／`65%` の入力・paste は共有 adapter で numeric text＋unit に取り込み、1回の field update にまとめる。IME 中や不正な suffix text を消さない。

Pack 02-A の unit change は numeric lexeme を維持して unit のみ変更する。`1.5 ×R → 1.5 px` であり、円半径を使う変換ではない。数値が不正でも text を保持する。同じ unit の選択は no-op。未生成の unit change は既存 Pending feedback の対象とし、Result は Generate まで変えない。

### 4.3 Auto と UI lifecycle

空欄/null は Auto。数値0は Auto と解釈しない。推定／解決 geometry は `Auto: 18.5 px` のように単位を含めて表示し、選択 unit の手入力値に見せない。

Pack 03-A では Auto に canonical unit はない。新規／再マウント時の次回入力用 unit は ×R。利用者は空欄のまま px を先に選び、次の `1.5` を1.5pxとして入力できる。Auto の unit preference だけは component の小さな transient state。非空の unit は slot が持つため mirror refs を増やさない。

空欄での selector change は request/Result と History を変えず、Session に書かない。panel 再マウント、Load、Reset で preference は既定に戻る。manual state の Undo/Redo は slot から text/unit を復元する。Auto に戻ったときの unit preference の復元は約束しない。scalar state が Auto であることは必ず復元する。この範囲は Pack 03 の署名対象。

### 4.4 History・Session・request の責任

- existing History input adapter の focus/change transaction を使う。数値編集と unit change は実際の操作単位で1 transaction。component と親の両方から transaction を開始しない。
- UI-only Auto preference を geometry intent に含めない。manual unit change は通常の draft change として保存・Undo/Redo・Pending に載る。
- History restore と Session Load の read-only projection で codec を使い、旧 unit の local state から restored scalar を上書きしない。
- `restoreCurrentWriterActiveConfig()` の保存 draft 優先を維持する。committed request から current draft を再構成して pending edits を失わない。
- canonical scalar は px/factor number/null のまま。Session 44、request 8、Circular slot schema 4 の変更を目的にしない。既存 current-reader の正当な受理を evidence で確認する。
- node helper で native import/export と writer validation を実行し、Web-only lexical draft と canonical pair を混同しない。未受理の invalid Session を新しく読み込めるようにはしない。

## 5. SOLID / KISS / DRY / YAGNI と ratchet

| 原則 | 実装で守る境界 |
| --- | --- |
| SRP | scalar parse、editor view/mutation、component presentation、History、Session、request を分離 |
| OCP | 既存 boundary を二つの field 用に拡張。汎用 units framework を作らない |
| LSP | 既存 valid scalar と canonical contracts を置換後も受理。旧科学的意味を保つ |
| ISP | component は scalar model、field identity、必要な owner action のみを使う |
| DIP | markup は component/editor に依存し、Session または Worker 内部を呼ばない |
| KISS | 二択 unit と既存 Generate。physical-unit conversion、追加 Apply、watcher 同期を作らない |
| DRY | validation/payload/projection は一つの scalar 解釈。旧 parser/formatter/直接 binding を同時削除 |
| YAGNI | unit registry、新 request schema、永続 UI preference、Linear 改修を追加しない |

source of truth は既存 slot scalar。Auto preference は scalar に unit がない場合の独立した transient UI intent で、valid scalar の複製ではない。canonical request/Worker/Result admission の owner/path は不変。codec/component は既存 owner 下の private decomposition。計画上 OE/PE/CB の増加は不要だが、実装後の source facts で確認する。例外条件が出たら完全な changed-scope sets と必要な maintainer decision を準備し、policy を弱めて通さない。

## 6. セッション分割と所有

同じ branch を複数セッションが同時に更新しない。順番は S00 → S01 → S02 → S03。各セッションは自己完結した [INSTRUCTION PROMPT](SESSION_PROMPTS/README.md) を読み、前セッションの commit SHA・受入結果を確認する。

| Session | 主な所有 | 完了条件 |
| --- | --- | --- |
| [S00](SESSION_PROMPTS/S00_AUTHORITY_AND_BOUNDARIES.md) | Product/evidence 文書、必要な既存契約の authority-only 更新、境界の disposable checks | 署名内容／codec 表現の evidence／trusted-base の成立を判定。未成立なら runtime へ進まない |
| [S01](SESSION_PROMPTS/S01_SCALAR_AND_EDITOR_MODEL.md) | scalar validation/payload、measure-editor codec、必要な既存 projection と focused Node tests | duplicate scalar 解釈を削除。typed/string/%/Auto/invalid の view・request round trip を確認 |
| [S02](SESSION_PROMPTS/S02_CONTROLS_AND_LIFECYCLE.md) | measure-input component、app wiring、markup/CSS、History/Session の必要最小接続 | 分離 controls、入力・単位変更・Auto・restore・History が承認内容どおり |
| [S03](SESSION_PROMPTS/S03_BROWSER_AND_HANDOFF.md) | browser regression、実 download/replay、既存 public technical docs、最終 review | 全受入 ID、policy gate、production/test/docs diff review、最終 commit/push |

各終了時に `SESSION_RESULTS/Sxx.md` を作る。記録は対象 commit/environment/input、commands と exits、実測結果、限界、再利用可能な evidence、次の session の前提を含む。署名、gate pass、browser verification を実測せず完了扱いにしない。

## 7. checkout・commit・push 規則

実装時は **remote の当該 branch を取得し、その既存checkoutをS00〜S03で直列に引き継ぐ**。セッションごとのclone、worktree、lock directory作成は必須にしない。別作業との隔離が必要な場合だけ、既存repositoryのGit履歴を共有するworktreeを一つ使う。他 session が使用中のcheckoutを switch/clean/reset しない。branch は origin/dev から既に作成済みなので、新しい session で main/dev から作り直さない。

詳細な開始・終了コマンドは [共通手順](SESSION_PROMPTS/README.md) に置く。各 prompt から必ず参照する。対象checkoutを他sessionが使用中でないことを確認する。server・browser・出力先は必要な検証に限って用意し、既存環境を使う場合もsourceと所有者を確認する。プロセスを停止する際は自分が起動した PID だけを対象にする。

各セッションの終了時は、成果物と結果を **当該 branch に commit して同名 remote branch へ push**する。未解決の Product/evidence boundary があっても、完了した独立文書・検証結果を commit/push して境界を引き継ぐ。runtime を途中まで公開して完成と呼ばない。署名がなく runtime が変更できないことは、未コミットのまま終わる理由にしない。

push 直前に remote state を確認する。remote が別 session によって進んでいたら停止してその内容を取り込み、共有ファイルの owner と前提を再確認する。force push、main/dev への直接 push、他 session の変更の revert は行わない。

## 8. 受入条件・検証

| ID | 必須結果 |
| --- | --- |
| C619-01 | tobacco fixture の四 field が numeric text＋px/×R selector で正しく読める。1440px /390pxで可視・操作可能、object textなし |
| C619-02 | Load/open/focus/blur は元 scalar／committed request／Result／History を変更しない。表示操作で新 Worker/helper/runを起動しない |
| C619-03 | 数値1.5の解釈は selector が決める。typed px/factor、旧 unitless／px／%、decimal/exponent、blank、precisionを検証 |
| C619-04 | 不正／未完成 textと選択unitを保持。row error、旧Result/request保持、訂正後Generate成功。0/負/nonfinite/不正unitはauto化しない |
| C619-05 | approved unit-change rule、IME、trim、keyboard、no-op、操作ごとのUndo/Redo。geometry更新はGenerate時だけ |
| C619-06 | 空欄unit選択→入力、Clear、panel remount、Load/Reset、History to AutoがPack03どおり。Auto geometryと選択unitが混同されない |
| C619-07 | valid draftのSave/Load、draft≠committed、disabled/inactive slots、settings-only、CLI-origin。canonical型/科学的意味/保存preview維持 |
| C619-08 | 未編集／同値編集の四canonical pairsとgeometry一致。変更fieldのみ意図したgeometry変更。実SVG downloadとnative Session replay |
| C619-09 | Linear、既存支持されたlegacy projection、canonical owner/path、History/Session admissionが退行しない |
| C619-10 | numeric/unit controlsが390pxで重ならず、tab順・accessible names・エラー関連付け・helpが使用可能 |

最低限の focused Node commands:

```bash
node --test tests/web/track-slot-display.test.mjs \
  tests/web/circular-track-slots.test.mjs \
  tests/web/track-slot-validation.test.mjs \
  tests/web/session-request.test.mjs \
  tests/web/session-draft-authority.test.mjs \
  tests/web/session-active-config-contract.test.mjs \
  tests/web/settings-only-session.test.mjs
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs
```

新 codec tests は `tests/web/circular-track-measure-editor.test.mjs`、browser spec は `tests/web/circular-track-measure-input.playwright.spec.js`。既存 CI の discovery に載せ、workflow や checker を本 runtime change で変更しない。

```bash
python tools/prepare_browser_wheel.py
GBDRAW_WEB_TEST_PORT=42619 npx --no-install playwright test \
  --config=playwright.functional.config.js --project=chromium --workers=1 \
  tests/web/circular-track-measure-input.playwright.spec.js
```

port は例。対象 checkout で空いている port を選び、そこから起動した server のソースを確認する。Node Playwright と Python Playwright の両方を確認する。Node がなければ Python で等価 check を行い、未実行の spec を pass としない。sandbox Chromium failure は適切な escalation で同じ check を再試行する。テスト全体は少なくとも30分の実行余地を与え、test-owned timeout を緩めない。

geometry が変わる予定はない。tracked reference を書き換えない。actual SVG comparisons は canonical geometry と実 artifact を確認し、artifact byte identity は同じ生成環境で適切な場合にのみ要求する。

## 9. documentation と公開境界

既存 `docs/REFERENCE/web-app.md` と必要なら `docs/REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md` に exact control semantics を記載する。新しい public page を増やさない。Gallery Session は入力 fixture として使う。公開 Gallery/tutorial capture がこの controls を説明している場合だけ、その owner recipe と applicable screenshot skill に従って更新する。social preview は変更しない。

プラン公開のために runtime、Gallery、wheel、reference、CI authority を変更しない。実装後のdev stagingとrelease promotionは既存workflowを使い、第二のdelivery pipelineを作らない。PR title/body 作成時は既存 write-clear-pull-request skill と language check を適用する。

実装全体の proposed commit title: `Make circular track values and units explicit`。

Summary: `Replace mixed scalar inputs with numeric fields and unit selectors while preserving scalar meaning, draft/Result separation, Session compatibility, and History.`
