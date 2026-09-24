# Issue #581 — Feature popup と Similarity Alignment 選択 UI の総合実装計画

- 状態: Product Decision ３件は承認済みで、authority-only PR #582 は 2026-09-24 に `dev` へマージ済み。S01 runtime は実装・検証済み。S02〜S04 は未実施。
- 作成日: 2026-09-24
- 対象: [Issue #581](https://github.com/satoshikawato/gbdraw/issues/581) と [2026-09-24 の設計提案コメント](https://github.com/satoshikawato/gbdraw/issues/581#issuecomment-5811037046)
- 固定実装ブランチ: `issue-581-feature-popup-alignment-ux-20260924`
- ブランチ作成時の base: `origin/dev` @ `8c54c14455fedebc55561c2c4dc952a44d8d78de`
- Product authority: PR #582、`origin/dev` @ `e52e3ea99fbbe4cdedfe9aa2606a91163e876764`（元ブランチ `issue-581-product-decisions-20260924`）

## 1. この文書の使い方

本書は過去の会話を知らない参加者向けに、問題、提案する製品結果、既存契約、
実装境界、作業順、検証、停止条件を記録する。各セッションの実行指示は同じ
ディレクトリの `SESSION_00_*.md` から `SESSION_04_*.md` に分けてある。
セッションは番号順に進め、各終了時に第10節の実施記録を更新する。

**すべての Issue #581 runtime、test、docs 作業には
`issue-581-feature-popup-alignment-ux-20260924` を使う。**
`main`、`dev`、他 Issue のブランチで作業しない。開始時にブランチ、HEAD、
upstream、作業ツリー、`origin/dev` との祖先関係を確認する。別の runtime ブランチは
作らない。新しい `origin/dev` が公開されても、既存差分を確認せずに rebase しない。
Product authority だけを先にマージするための別ブランチは第7節に記す。

計画の一次資料は上記 Issue とコメントである。実装時には更新の有無を再確認し、
本書の記述だけを最新の製品判断として扱わない。

## 2. 問題と目標

Web の feature popup は feature を開くたびに `Record actions` の回転フォームを
タブより上に展開する。一般的な feature 確認・編集が下へ押し出され、回転フォーム内の
`Cancel` は popup 全体を閉じる。

Linear Similarity Groups の `Select alignment anchors` は候補と record に内部 ID を
表示し、候補を一つ選ぶたびに Pyodide の Python Resolver を呼ぶ。複数 record の選択中、
radio が待機状態になり、比較図を見ながら判断しにくい。

提案する利用者の一連の操作は次のとおり。

1. 既存の `Align` または `Align & orient` から exact reference を指定する。
2. Python Resolver が初回に曖昧と判定した record だけ、選択パレットを開く。
3. 図の上下の対応を保ったまま候補の gene/locus、source 座標、strand を読み、
   palette の radio または図上の候補から Select／Skip する。
4. 選択は即時にローカル下書きへ反映する。全曖昧 record が解決されるまで Apply を
   無効にする。
5. Apply 時に全選択を一回だけ Python Resolver で検証し、成功した plan だけを既存の
   生成・Result・History 経路へ渡す。失敗時は修正可能な下書きを保持する。

数値としての「0 ms」は保証しない。受入条件は、個々の選択で Python/Worker 呼び出しと
diagram 再生成を行わず、radio の待機表示を生じさせないことである。Vue の局所的な
表示更新は必要な動作である。

## 3. 既存の権威と実装事実

- `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の PD-OI-026 は exact reference、
  record ごとの独立した Select／Skip、automatic resolution の優先順位を定める。
  図上の位置、代表 status、score 等を候補の自動順位に使わない。
- 同文書の PD-OI-027〜029 と PD-OI-031 は orientation、active plan、Reset／History、
  Web と Python/CLI の境界を定める。PD-OI-032 は popup からの record 回転、
  対象 record だけの更新、Cancel／失敗時の Result 保全を定める。
- 追加判断 PD-OI-033〜035 は PR #582 で `origin/dev` @ `e52e3ea9` にマージされ、
  固定実装ブランチにも取り込まれた。
- `gbdraw/layout/similarity_alignment.py` は record ごとに選択を処理する。
  Python Resolver が候補の適格性、曖昧性、plan を決める唯一の所有者である。
- `gbdraw/web/js/app/similarity-alignment.js` は現在 `answer()` から選択ごとに
  `resolveRequest()` を呼び、`applyPlan()` 失敗時に `clearDraft()` する。
- `gbdraw/web/js/app/feature-editor/svg-actions.js` には stable feature identity による
  候補ハイライトがある。`gbdraw/web/js/app/orthogroups.js` の enriched member と
  feature catalog には表示用 metadata がある。
- `gbdraw/web/js/app/app-setup.js` は popup feature が変わるたびに record rotation の
  domain draft を開く。`gbdraw/web/index.html` はその `active` だけでフォームを
  タブより上に表示し、フォーム内 `Cancel` を popup close に結び付けている。
- 現行の `tests/web/similarity-alignment-actions.test.mjs` と
  `tests/web/similarity-alignment-ui.playwright.spec.js` は Resolver と UI の既存経路を
  検証する。record rotation の単体・ブラウザテストも既存である。

これらは base の実装事実であり、製品上の未確定な選択を自動承認しない。

## 4. 提案する製品結果と境界

### 4.1 Feature popup

- rich popup の `Edit` 冒頭に、見つけやすい `Record actions` 折り畳み欄を置く。
  simple popup も同じ一つのフォームを使う。初期状態は閉じる。
- `featureRecordRotationDraft.active` は対象 feature との domain binding として維持し、
  UI 可視性は popup に属する一つの transient な `recordActionsExpanded` で表す。
  feature 切替、popup close、フォーム内 Cancel で閉じる。
- フォーム内 Cancel は draft を reset して欄を閉じる。feature popup 自体と既存 Result は
  保つ。再び開いたときは現在の popup feature へ結び直す。
- Anchor／offset／orientation／feature-end の計算、target-only Apply、History、
  Session、sidebar 側の挙動は既存 owner を使う。別の回転 engine を作らない。

専用 `Record` tab は Issue コメントに記載された別案である。`Edit` 内の折り畳み欄は
承認された `A / EDIT_DISCLOSURE` に対応する。実装に進む条件は第7節に記す。

### 4.2 候補の表示

- record の表示名は現在の Linear 表示 record の definition／accession 等から解決する。
  内部 `recordKey` は matching に残し、表示名が欠けても `Record 2` のような順番を示す。
- 候補の主行は番号、gene、locus tag、source 座標範囲、display strand とする。
  gene／locus がなければ product、最後は `CDS · 12,341..13,720 bp (+)` のような
  type／座標へ戻す。長い文字列は折り返し・省略表示して詳細を閲覧できるようにする。
- 内部 ID、role、representative、direct edge evidence は `Details` disclosure に置く。
  これらを自動推薦や信頼度の順位として扱わない。
- view model は exact identity を key に、enriched member と catalog の既存 facts から
  作る。Python helper request や保存済み plan に UI 専用 label を加えない。

### 4.3 選択と Apply

- 初回 Resolver 応答と exact request を保持し、曖昧 record の候補集合を固定する。
  選択の唯一の mutable source は recordKey ごとのローカル choice とする。
  palette と図上クリックは同一の `selectCandidate(recordKey, anchor)` を呼ぶ。
- 個々の Select／Skip は同期的に choice だけを更新する。全曖昧 record に choice があれば
  Apply を有効にする。この段階で plan が存在すると見なさない。
- Apply で choices を一回だけ request に投影し、Python Resolver 応答を既存の
  `validateResolution()` で検証する。`resolved` と plan が揃った時だけ
  `applyPlan()` を呼ぶ。曖昧さが残る／応答不整合なら生成しない。
- Python 検証失敗・生成失敗では選択と baseline を保ち、理由を表示して再試行を許す。
  成功、明示 Cancel、外部状態の無効化では既存の状態所有者が後始末する。
- パレット中に source、crop、group、committed Result が変化した場合、
  開始時の artifact と一致するか検証する。古い選択を適用せず、再開始の理由を示す。
  pan／zoom はこの検証を無効にしない。
- 曖昧 record がない場合の既存の自動適用、exact reference、orientation、
  Reset、Undo/Redo、Session 再生は維持する。

### 4.4 パレットと図上ガイド

- 選択 UI は全画面 backdrop と `aria-modal`、focus trap を外した non-modal palette とする。
  パレットは画面内で drag 可能にし、位置は transient UI state のみで管理する。
  narrow viewport では位置と高さを画面内に収め、radio、Skip、Apply が使えるようにする。
- reference feature の描画上の中心 x を通る細い縦線と、曖昧候補の番号 badge を
  **プレビュー専用 overlay** に置く。SVG Result・download・Session には混入させない。
  pan、zoom、resize、Result 置換で位置を再計算または破棄する。
- 描画 feature が見つからない、非表示、画面外、複数断片で特定できない場合も、
  候補一覧からの選択は維持する。ガイドや badge は Resolver の usable 判定を変えない。
- SVG の既存 feature identity lookup／hover owner を再利用する。候補 feature の
  click は popup open より先に一回だけ選択へ送る。非候補 click と通常の図操作は
  既存の意味を保つ。palette と feature の hover は双方向に同期する。
- keyboard だけでも全候補と Skip を選べる。non-modal なので Tab を閉じ込めない。
  Escape／Cancel と focus return は他の popup や dialog と競合しないよう扱う。

## 5. 所有者とデータフロー

`gbdraw/web/index.html` は markup／CSS、`app/app-setup.js` は接続だけを所有する。
`app/similarity-alignment.js` は draft、request、Resolver、Apply の状態遷移を所有する。
`app/feature-editor/svg-actions.js` は SVG feature lookup、hover、pointer action と
preview overlay の lifecycle を所有する。大きくなる場合だけ `feature-editor/` 配下の
private helper へ分解し、新しい状態 owner は作らない。record rotation の domain owner は
`app/record-display/feature-record-rotation.js` のままとする。

`start → Python resolve once → immutable ambiguities + local choices
→ radio/canvas selection (no Worker) → Apply → Python resolve once
→ validated plan → existing runAnalysis/Result/History`

JavaScript に Python の候補順位、crop 適格性、座標変換を再実装しない。
canonical request schema、Worker protocol、Python API、CLI、Session writer は変更しない。
UI 追加 state は保存対象にしない。

## 6. 実装セッションと受入条件

| 順序 | 指示ファイル | 所有する成果 | 開始条件 |
| --- | --- | --- | --- |
| S00 | `SESSION_00_PRODUCT_PREFLIGHT.md` | outcome 別の authority 調査、必要な Decision Pack と権威経路 | 本書と最新 base |
| S01 | `SESSION_01_POPUP_AND_DRAFT.md` | 回転欄、候補 view、ローカル選択、Apply 時 batch 検証 | 該当 runtime outcome が許可済み |
| S02 | `SESSION_02_FLOATING_PALETTE.md` | non-modal palette、drag、focus、responsive layout | palette outcome が許可済み、S01 完了 |
| S03 | `SESSION_03_CANVAS_INTEGRATION.md` | guide、badge、双方向 hover／click、cleanup | S02 完了 |
| S04 | `SESSION_04_ACCEPTANCE_AND_HANDOFF.md` | 実ブラウザ検証、必要な既存 docs 更新、全体差分監査 | S01〜S03 完了 |

全セッションの最終受入は次を満たす。

1. popup は初期状態で回転フォームを展開せず、フォーム Cancel 後も popup と Result を保つ。
2. 複数の曖昧 record で Select／Skip を何度切り替えても、各 click の Worker 呼び出しは 0。
   Apply は全選択を一回の Resolver request で検証する。
3. biological label、source 座標、strand、record 表示名が読める。欠損 metadata の fallback と
   同名 gene の識別ができ、内部 ID が主見出しにならない。
4. palette／図上の選択・hover が一致し、pan／zoom／resize 後も guide と badge が対応する。
   ガイドが出ない feature も一覧から選べる。
5. Cancel、validation failure、generation failure、stale／superseded completion は
   committed Result と History を変更しない。再試行可能な失敗では選択を保持する。
6. rich／simple popup、desktop／390 px viewport、pointer／keyboard、no-modal focus、
   sanitized Result、保存／再読込後の既存 alignment plan を確認する。

テストの主な場所は `tests/web/similarity-alignment-actions.test.mjs`、
`tests/web/similarity-alignment-ui.playwright.spec.js`、
`tests/web/feature-popup-record-rotation.test.mjs`、
`tests/web/interactive-svg-v3.playwright.spec.js` である。必要な Python resolver 回帰には
`tests/test_similarity_alignment.py` と `tests/test_similarity_alignment_web_adapter.py` を使う。
ブラウザ確認では Node と Python の Playwright を確認し、Node spec が実行できない環境でも
同等の Python Playwright による実ブラウザ確認を行う。

## 7. Product Impact とブランチ運用

Product Decision Owner `satoshikawato` は 2026-09-24、次の３件の
`PRODUCT_DECISION` 文面をそのまま承認した。文面の全フィールドは
`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` に PD-OI-033〜035 として
記録し、PR #582 で `dev` にマージした。以下の表はその選択を要約する。

| Concern | 承認された選択 | 権威の所在 |
| --- | --- | --- |
| `web.feature-popup.record-actions-presentation` | `A / EDIT_DISCLOSURE`（PD-OI-033） | `origin/dev` @ `e52e3ea9` と固定実装ブランチ |
| `web.similarity-alignment.choice-and-retry` | `A / LOCAL_BATCH_RETRY`（PD-OI-034） | 同上 |
| `web.similarity-alignment.canvas-interaction` | `A / FLOATING_GUIDE_CANVAS_PICK`（PD-OI-035） | 同上 |

既存の PD-OI-026〜032 と承認済み３件の全文を実装前に照合する。新たな material な
製品結果差が見つかった場合は `PRODUCT_IMPACT_RATCHET.md` に従って調査し、
未決の outcome に依存する runtime を開始しない。Issue コメントだけを durable authority
とみなさず、候補 authority を runtime ブランチへ直接積んで自己承認しない。

authority-only 変更は PR #582 で `dev` にマージ済みであり、固定実装ブランチにも
取り込んだ。Issue #581 runtime は固定実装ブランチに置く。権威変更、PR、merge、push は
明示された許可範囲に従う。

## 8. 設計原則と対象外

- **SOLID:** domain selection、SVG interaction、presentation、generation に一つずつ owner を
  保つ。radio と canvas は狭い同じ選択 interface に依存し、UI から Python の規則を参照しない。
  継承階層や汎用 component interface は導入しない。
- **KISS:** rich／simple の共通フォームと既存の Resolver、feature lookup、Result 経路を使う。
  一つの状態を別の可変 state へ複写しない。
- **DRY:** candidate identity、Select／Skip、Apply、hover、Cancel の判断点を各一箇所に集める。
  旧 per-click Resolve 経路を残して二つの実行経路にしない。
- **YAGNI:** score/ranking、debounce、別 Worker、Session field、汎用 overlay framework、
  CLI/Python 用の新 UI、canvas からの推測 alignment を追加しない。

architecture-bearing な変更には `ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` を適用する。
既存 owner/path を維持する通常の変更は簡潔な前後の owner/path 根拠を記録し、
例外条件に当たる場合だけ所定の OE／PE／CB 完全比較と判断を行う。
public docs は既存の該当ページに集約し、機能ごとの新規ページを安易に増やさない。

## 9. 最終チェック

- 対象 source、test、docs の diff を別々に確認し、無関係な変更を含めない。
- focused JS／Python／browser checks と必要な repository gate を実行する。
- 実ブラウザで最終 UI を読める縮尺で確認する。公開画像を更新する場合は
  元の再現 recipe から生成して目視確認する。
- プレビュー overlay が SVG download と Session に混入しないことを確認する。
- `git status`、branch/upstream、commit 履歴を確認し、作業範囲に応じた英語の
  commit title と短い summary を提示する。push／PR はそのセッションの許可に従う。

## 10. 実施記録

| Session | 状態 | HEAD／base・主な変更 | 検証・残件 |
| --- | --- | --- | --- |
| S00 | 完了 | PD-OI-033〜035 は authority-only PR #582 で `origin/dev` @ `e52e3ea9` へマージし、固定実装ブランチへ取り込んだ。 | JSON receipt ３件の構造確認、`git diff --check`、Web 変更ゲート PASS、PR の CI 成功。S01 開始可能。 |
| S01 | 完了 | 開始 HEAD `28e043d5`（固定ブランチ、`origin/dev` @ `e52e3ea9` を祖先に含む）。Edit 内の Record actions 開閉欄、Cancel 後の popup／Result 保持、生物学的候補表示、record 単位のローカル Select／Skip、Apply 時の一括 Python 検証、再試行可能な失敗後の選択保持を実装。既存 owner／canonical path と Session writer は維持。 | JS focused 33 件、Python Resolver 39 件、Chromium focused 4 件が PASS。Ruff、`git diff --check`、Web change gate PASS（通常の architecture review 要）。広めのブラウザ実行では S01 と無関係な released v40 の legacy materialization ケースで期待エラーが空となり FAIL。S02 開始可、S03／S04 は未実施。 |
| S02 | 未実施 | 変更なし | S01 完了。開始可能。 |
| S03 | 未実施 | 変更なし | S02 待ち。実装・検証なし。 |
| S04 | 未実施 | 変更なし | S01〜S03 待ち。受入検証なし。 |
