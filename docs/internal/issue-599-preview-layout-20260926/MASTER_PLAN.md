# Issue #599 — Preview 配置・操作の修正実装計画

作成日: 2026-09-26。対象: [gbdraw Issue #599](https://github.com/satoshikawato/gbdraw/issues/599)。
計画作成時の最新 `origin/dev`: `d457b7189b137185a8dec800819a312c30b969fa`。
実装ブランチ: `fix/issue-599-preview-layout-20260926`。
リポジトリ: `https://github.com/satoshikawato/gbdraw.git`。

この文書は、新規参加者が仕様、責任分担、実装順序、検証、Git 操作を理解するための総合計画である。計画作成セッションでは文書だけを追加する。実装・検証の進捗は [SESSION_STATUS.md](SESSION_STATUS.md) に記録する。

## 問題と承認済みの製品結果

gbdraw Web は、ブラウザー内の Python renderer でゲノム図を生成し、Result Preview で配置を仕上げ、SVG などへ出力する単一ページアプリである。`Result` は保存・出力する生成物、`candidate` は公開前の新しい生成物、`committed request` は現在の Result を生成した確定済みの要求を指す。生成中に画面で編集された要求とは区別する。

Issue #599 は次の3件を扱う。Product Decision Owner `satoshikawato` は2026-09-26に3件それぞれの推奨 A を明示承認した。完全な承認本文と機械表現は [承認記録](00_APPROVED_PRODUCT_DECISIONS.md) に保存する。3つの関心は独立して追跡し、新しい一括承認の仕組みは作らない。

| 関心 / Issue 内の識別子 | 実装する製品結果 | セッション |
| --- | --- | --- |
| `web.composition-decoration-continuity` / BUG-05 | 同じ図の legend、plot title、Linear scale bar の手動差分を、新しい自動配置へ1回加える。対応不能な非ゼロ差分がある場合は新候補を公開せず、旧 Result を保持して Reset または設定修正へ案内する | S01 |
| `web.layout-edit-affordance` / BUG-06 | Layout edit の明示切替を維持。OFF でも対象と有効化方法を説明し、ON で装飾をドラッグできる | S02 |
| `web.preview-search-placement` / BUG-18 | 検索を上部専用行、toolbar を下部専用行へ固定。検索の自由ドラッグと重複する位置計算を削除する | S03 |

配置継承は「新しい自動配置 + 手動差分」であり、絶対位置の固定ではない。clipping や overlap はあり得る。自動 clamp、差分縮小、viewBox 拡張は追加せず、既存 padding と Reset で回復する。

対象外: diagram 全体・個別 record・padding・legend entry order の新しい再生成継承保証、狭幅 editor の上下 dock、Similarity Group alignment review の仕様変更。これらに関係する既存機能は維持する。特に他の修正が同じ Generate 境界へ入っている場合は、装飾の差分と record の差分を混同せず、既存の単一 commit 経路へ統合する。

## 根拠と限界

[基準観測データ](base-observations.json) と [基準観測スクリプト](observe_base.py) は上記 dev SHA の修正前の証拠である。Gallery の `HmmtDNA_basic_circular` と `lambda_basic_linear` を使い、実 renderer と pointer 操作で確認した。zoom 0.6 で24×12 CSS px移動すると、legend/title/Linear scale の差分は約40×20 SVG unitsとなり、次の Generate でゼロへ戻った。Layout edit OFF の legend/title は `grab` cursor だが、装飾の drag owner は開始せず canvas pan になる。

drawer open の1024×768では検索左端14に対して Preview 左端358、900×768では検索右端1030に対して Preview 右端882だった。固定360pxの JS 退避と CSS transform 上書きが競合する。観測した初期10条件では search/toolbar 矩形交差は0であり、Issue の toolbar overlap をこの観測で再現したとは扱わない。狭幅の DOM 座標だけでは実画面上の可視性・touch・hit target は証明できない。

基準の Node 3ファイルは8/8成功した。観測スクリプトは不具合が存在することを検査するため、修正後の合格試験として実行しない。修正後には下記の反転した受入条件を独立して検査する。参照先の `GUI_AUDIT_DEV_20260926.md` は調査時の base に存在せず、未確認の監査本文を根拠にしない。観測の browser/wheel 情報は JSON に残す。

## 実装開始前の製品契約

必読: [Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md)、[Decision Packet template](../PRODUCT_DECISION_PACKET_TEMPLATE.md)、[OIPC](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)、[Architecture Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)、[Web Change Policy](../WEB_CHANGE_POLICY.md)。各セッションで repository `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md` とその時点の最新規約も読む。

計画作成時の base は OIPC revision 16、PD-OI-001–035、`tools/web-product-decisions.json` の active decisions は空である。今回の3件は developer preflight で発見した unmapped concern なので、既存 OIPC を恒久的な製品契約の記録先とする。新しい JSON registry、Product checker、decision state machine は追加しない。承認記録中の JSON はレビュー用の非実行表現であり、active authority ではない。

S00 は最新 dev の authority を再調査し、承認した3つの完全な結果を authority-only 変更として記録する。番号は最新 base で採番する。同じ結果が既に active なら再利用し、二重登録しない。mapped concern になっていれば既存 map/BD の経路を使う。競合する結果が見つかった場合は該当部分だけ停止して差分を報告し、別の製品結果を自動選択しない。既存 outcome の置換が必要なら complete replacement、scenario revision、Supersedes を規約どおり記録する。

**S01–S03 の runtime 実装は、必要な authority-only 変更が dev へマージされ、取得した `origin/dev` 上で承認内容を確認するまで開始しない。** 同じ candidate 内の契約で、その candidate の runtime を承認しない。OIPC-C05/C06 の intent 維持、PD-OI-016 の failure isolation、PD-OI-035 の mobile alignment palette coverage、および mapped contracts を維持する。選択 ID の一致だけでは不十分で、独立した全必須結果を AND-of-OR の要件として照合する。

この計画のコミット・実装ブランチへの push は依頼で承認済み。PR の作成、dev への merge、deploy は含まれない。S00 はレビュー可能な authority-only 差分をコミット・push した上で、必要な PR 作成・merge の具体的対象を示す。権限のある担当者が別途承認・実行した merge を取得して確認する。権限待ちを製品判断の再承認へ置き換えない。

## アーキテクチャと状態の責任

### 装飾差分の継承

Python は自動配置を生成する。手動差分の正本は既存 Result の clean SVG とし、UI offset refs は適用済み SVG から派生させる。

```text
Generate 開始: 旧 committed Result/request から小さい装飾 snapshot
  → canonical typed request / 既存 Worker render
  → 既存 SVG sanitizer / editor mutation
  → 対象同一性の照合 → 新 automaticTranslation + delta
  → 全 candidate Result が成功 → 既存 atomic commit
  → Preview bind/ready → 既存 History replacement
```

| 責任 | 既存 owner / 変更方針 |
| --- | --- |
| 差分取得・照合・有限値検査・適用 | `gbdraw/web/js/app/legend-layout/composition-actions.js` と既存 `legend-layout/`。必要な private helper を同じ責任範囲に置く。演算は `compositionUserDeltas()` / `applyCompositionUserDeltas()` を再利用 |
| Generate orchestration | `app/run-analysis.js`。async 前に snapshot を受け取り、既存 lock/token/Cancel/stale guard を維持 |
| 依存の接続 | `app/app-setup.js`。小さい snapshot/transform provider の wiring のみ。legend API 全体を orchestrator へ渡さない |
| 候補変換・admission | `app/candidate-render.js` → `services/svg-result-ingestion.js` の既存 `transformSvg` / `callerTransforms`。sanitizer/editor mutation 後、commit 前に適用 |
| drag と保存 | `app/legend/drag-actions.js`、`app/legend-layout/diagram-drag.js`、`canvas-actions.js` の既存経路 |
| request/resource identity | `services/session-request.js`、既存 resource/catalog。resource ID や filename だけを source 内容一致の証明にしない |
| History/Session/ready | `services/history-snapshot.js`、`services/config.js`、`app/preview-runtime.js` の既存復元・bind/ready。再加算や第二 bind を作らない |

同一性は mode/grouping、既存 validated source identity/content、biological region、record identity 集合で証明する。表示順・prefix・配列 index は identity ではない。既存 validated token/digest を使い、追加の genome read/hash を行わない。不明な旧 Session、重複・曖昧な出力は「対応不能」とする。

batch の全出力を照合し、全 candidate の変換が成功するまで公開しない。選択外は必要な旧 clean SVG だけを読み、第二の永続 delta map は作らない。scale は `#length_bar` の役割で照合し、旧 `primary[i]` を新 `primary[i]` へ丸ごとコピーしない。recordTranslations/alignment の二重適用を禁止する。

通常 Generate、committed-candidate の render、automatic reflow のすべてを同じ candidate 変換境界へ通す。Save/Load・Undo/Redo は既に保存済みの位置を復元し、継承をもう1回行わない。render 待機中に旧 Result から新たに差分を取り直さない。非ゼロ対応失敗は公開前に止め、旧 Result/request/History を保持する。ゼロ差分は追加確認を要求せず、空 transform のために EMPTY fast path を MUTATING へ変更しない。

### 操作説明

Layout edit OFF は従来の pan。supported decoration target に `help` cursor、hover 枠、`Turn on Layout edit to move this item` の説明を出す。toolbar に常設説明を置き、keyboard focus と touch でも同じ情報へ到達させる。ON は `grab`、drag 中は `grabbing`。mode を自動 ON にしない。

対象判定・feature/label/legend 個別編集・Shift/Ctrl 優先順位は既存 drag owner を再利用する。第二 pointer router は作らない。hint は Preview の派生表示であり canonical 値・History を更新しない。HTML wrapper の hint/CSS を優先し、SVG transient class が必要な場合だけ `services/svg-serialization.js::stripTransientPreviewState()` の明示 inventory へ追加する。保存 Result・plain/interactive SVG・PNG/PDF に hint を混入させない。

### 検索・toolbar の構造

```text
Result Preview
  header / notices
  search: 上部専用 auto row
  workspace: 既存 canvas と既存 editor drawer
  toolbar: 下部専用 auto row、必要時 wrap
```

geometry の owner は `gbdraw/web/index.html`。search/toolbar を workspace 外に置き、drawer の containing block を workspace とする。同じ search/canvas/editor DOM を維持する。drawer の tab/visibility/Close/Escape は `app/right-drawer.js` が引き続き所有する。

`app/feature-search/preview-actions.js` の固定 `RIGHT_DRAWER_WIDTH_PX`、dragOffset/activeDrag、位置 computed、search drag start/move/stop と関連 global listener を削除する。app-setup exports と template の style/drag binding も削除する。旧 search/toolbar absolute・translate・container-query override を同じ変更で撤去する。無関係な popup の寸法は対象にしない。

container 幅に応じた wrap、`min-width: 0`、scroll を同じ CSS owner に置く。新しい座標 ref、ResizeObserver、collision solver は不要。通常高さ740px以上の受入条件で workspace 高さ200px以上。短い画面・soft keyboard・200% zoom は Preview/card または文書の scroll で全操作へ到達させる。fixed height と overflow clip で隠す解決は認めない。query/field/regex/active match/focus を drawer/resize だけで reset しない。

## SOLID・KISS・DRY・YAGNI の実装／運用規則

| 原則 | 必須の具体化 |
| --- | --- |
| SRP | 装飾座標、検索意味、chrome geometry、candidate transaction、製品判断を各 owner に分ける |
| OCP / LSP | 既存 transform seam を拡張し、admission/ready/失敗時保持の契約を維持 |
| ISP / DIP | orchestrator へ最小 provider を注入。generic service registry や巨大 API は追加しない |
| KISS | automatic + delta と CSS の専用行で解く。設定 toggle、observer 系列、互換 fallback を増やさない |
| DRY | 差分演算、sanitize/commit、search/editor DOM、geometry、製品契約の正本を複製しない |
| YAGNI | 新 Session/request schema、全要素 position store、drag engine、Worker、全 checkpoint clone、先回りする抽象化を追加しない |
| ワークフロー | 同一 remote branch へ順番に引き継ぐ1 writer。セッション境界はレビュー可能な責任単位。各人の専用 clone で作業し、前セッションの公開済み状態を取得する |

architecture-bearing 差分は owner と canonical path の before/after、旧経路の撤去を簡潔に記録する。通常の非増加経路に静的推測の OE/PE/CB 数値を付けない。Ratchet の例外条件（OE/PE 増加、新 persisted compatibility 等）に該当した場合だけ完全な sets・算術・maintainer 判断を用意する。チェック失敗を policy/checker/baseline の緩和で回避しない。

## セッションと順序

| ID / INSTRUCTION PROMPT | 前提 | 所有する成果物 / 完了条件 |
| --- | --- | --- |
| [S00](SESSION_00_AUTHORITY_INSTRUCTION_PROMPT.md) | 承認記録と最新 dev | 3 concern の恒久記録、authority-only 差分、受入条件との対応。commit/push 後に authority dev merge の確認を引き継ぐ |
| [S01](SESSION_01_COMPOSITION_CONTINUITY_INSTRUCTION_PROMPT.md) | 必要 authority が origin/dev に存在 | 装飾 snapshot/identity/候補変換、failure/zero/batch 検査。C01–C08、R01 の該当範囲 |
| [S02](SESSION_02_LAYOUT_AFFORDANCE_INSTRUCTION_PROMPT.md) | S01 公開済み | OFF/ON の説明、既存 gesture 優先、transient 除去。A01–A02 と保存/出力の検査 |
| [S03](SESSION_03_PREVIEW_CHROME_INSTRUCTION_PROMPT.md) | S02 公開済み | search/workspace/toolbar 構造と旧位置経路撤去。P01–P03、既存検索・navigation の検査 |
| [S04](SESSION_04_INTEGRATION_INSTRUCTION_PROMPT.md) | S01–S03 公開済み | 統合受入、required gates、既存ユーザー文書、production/tests/docs/generated 各 diff review、完了判定 |

S01–S03 は個別に実装証拠を保存し、S04 は有効な既存証拠を再利用する。統合差分・環境変更・失敗・未解決の懸念がある箇所だけ再検査する。5つの branch や5つの authority store は作らない。実装担当者が前提未達なら独立してできる調査・証拠整理を行い、進捗に理由を記録して commit/push し、runtime 開始条件は維持する。

## 必須 Git 手順 — 全セッション共通

**必ず公開済みの `fix/issue-599-preview-layout-20260926` を取得して使う。** 毎回 origin/dev から別ブランチを切り直さない。共有作業ツリーの branch switch/reset/stash、他セッションの untracked ファイルの追加、既存 worktree の横取りは禁止。前セッションの終了と公開 SHA を確認してから開始する。同じ remote branch の同時 writer は許可しない。

新しい専用 clone を作る手順:

```bash
task_dir=$(mktemp -d /tmp/gbdraw-issue599-session.XXXXXX)
git clone --branch fix/issue-599-preview-layout-20260926 \
  https://github.com/satoshikawato/gbdraw.git "$task_dir/repo"
cd "$task_dir/repo"
git fetch origin
git pull --ff-only origin fix/issue-599-preview-layout-20260926
test "$(git branch --show-current)" = "fix/issue-599-preview-layout-20260926"
test "$(git rev-parse --abbrev-ref '@{upstream}')" = "origin/fix/issue-599-preview-layout-20260926"
git status --short
```

`task_dir` は取得した作業ディレクトリ。実行前に実パスを確認する。既存の専用 clone を再利用する場合も、clean tree、branch、upstream、前セッション SHA を確認して fetch/pull する。main/dev を upstream にしない。

S00 の authority が dev に merge された後、S01 は必要な記録を `git show origin/dev:docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` 等で確認し、clean な実装ブランチへ取り込む。`git merge --ff-only origin/dev` が ancestry の都合で不可能なら、両側の差分を確認して `git merge --no-edit origin/dev` で通常 merge する。公開済みブランチの rebase、reset、force push を行わない。conflict は両方の有効な変更を維持して解消する。source の同一境界に他 issue が入った場合は実装責任を増やさず統合する。

終了時は以下を必ず行う:

1. 対象検査・必要 gate・production/test/docs/generated の該当 diff review を完了する。
2. `results/S00.md`～`results/S04.md` の担当ファイルを作成し、SESSION_STATUS を更新する。結果には入力 SHA、authority base SHA/記録、変更点、実行コマンド・exit/result、artifact、制限、次の開始条件を記録する。未実施を合格と書かない。
3. branch/upstream を上記 `test` で再確認。`git fetch origin` 後に remote branch が入力 SHA から予期せず進んでいれば、その差分と担当者を確認する。force push や他者の変更を上書きしない。
4. 対象パスを列挙して `git add -- <対象パス>`。`git add .` は使わない。`git diff --cached --check` と staged diff を確認し、各セッションの英語 commit title で commit する。
5. **各セッション実装終了後は commit と同名 remote branch への push を必ず行う。** 次の公開先だけを使う。

```bash
git push origin HEAD:refs/heads/fix/issue-599-preview-layout-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-599-preview-layout-20260926
git status --short
```

remote SHA が local HEAD と一致することを確認する。push のエラー・中断時は remote の実状態を調べてから再試行する。終了報告には commit SHA、branch、検証結果、未完了条件、英語 proposed commit title と short summary を含める。変更不要だったセッションも検証結果文書を意味のある commit として公開し、空 commit は作らない。PR 作成・merge・deploy は明示された別途の権限境界を守る。

## 受入条件

| ID | 必須の結果 |
| --- | --- |
| C01 | 実 Circular/Linear で装飾 drag → 同じ Generate 2回。約40×20差分を保持し、倍増しない |
| C02 | color/font/title text/legend side を変更して新 automatic + 同じ delta。biological coordinates/comparison/scale 値は不変 |
| C03 | primary の順序に依存せず scale のみ照合。既存 recordTranslations と active alignment に二重加算なし |
| C04 | batch 複数出力を異なる差分で編集。選択・Result順・prefix変更、Save/Load後も対応 identity のみに継承 |
| C05 | target消失、source内容/crop/mode/grouping/record集合変更、不明・曖昧identityの非ゼロは公開前に停止。旧 Result/request/History保持。対象 Reset/Reset Layout または設定修正後に成功 |
| C06 | fresh/zero/対象なしは自動 Generate、追加確認・空 mutation parse なし |
| C07 | render/transform/bind失敗、Cancel/stale/supersessionで旧 Result/request/History保持。失敗 History entryなし |
| C08 | GenerateとResetのUndo/Redo、Session save/load → Generate、current Result exportで同じ配置。hint混入なし |
| A01 | OFF panとhelp/hover説明、ON grab/active grabbing。feature click、label、legend個別編集、Shift/Ctrlの既存優先維持 |
| A02 | keyboard/touchで説明と有効化へ到達。Result/load/Historyでrebind。保存Resultとplain/interactive SVG、PNG/PDFにhintなし |
| P01 | 幅1440/1024/900/768/390、高さ844/740/480、container境界前後、settings resize、drawer開閉系列、200%zoom、soft keyboardで可視性・非被覆・実hit target・scroll到達性。740px以上の通常条件はworkspace>=200px |
| P02 | query/field/regex/active match/focus保持。Prev/Next/Open/Enter、toolbar、drawer tab/Close/Escape、pan/zoomが使用可能 |
| P03 | 同じsearch/canvas/editor DOM。固定360px退避、search drag global listener、競合CSS・template bindingを撤去 |
| R01 | zero/nonzero/batchのparse/serialize/bind/Worker回数を比較。追加genome read/hash、全checkpoint clone、二重bindなし。既存responsiveness guardrail合格 |

C01/C03/C05/C07 の新しい保証には、修正前に失敗する意味のある回帰検査を設ける。base 用の disposable clone と既存観測を使い、実装ブランチを reset しない。browser は screenshot だけでなく elementFromPoint、実 pointer/keyboard、矩形・clipping を確認する。矩形交差0でも page外・clipping・drawer下に隠れていれば不合格。

拡張先は既存 `tests/web/composition-layout.test.mjs`、`candidate-render.test.mjs`、`run-analysis-simple-path.test.mjs`、`preview-feature-search.test.mjs`、`legend-layout-actions.test.mjs`、composition browser specs、`preview-navigation.playwright.spec.js`、`responsiveness-guardrail.test.mjs`。必要な契約を示す検査を選び、実装行の写しや弱い mock だけにしない。

## 検証環境・文書・完了

browser 検査時は Node `@playwright/test` と Python Playwright の両方を調べる。Node がなければ Python の等価な targeted check を使う。Chromium sandbox エラーは同じ local check を必要な escalation で再実行して診断する。wheel が必要なら `python tools/prepare_browser_wheel.py` で準備し、source と wheel の Python/TOML 一致を確認する。generated wheel を commit せず、試験準備で cache-bust を更新しない。global Worker-ready を待たず、自分が開始した operation と既存 Preview readiness receipt を待つ。

既存 [Web Change Policy](../WEB_CHANGE_POLICY.md) の required gates と trusted-base checker を使用する。candidate 側の checker・map・policy で候補を自己承認しない。全 matrix/Gallery staging が必要な gate は正確な integrated dev SHA に対して実行し、merge前の証拠と merge後の staging を区別する。リモートCI確認は5分以上の間隔。テストをタイムアウトと扱う前に30分以上確保し、途中出力を監視する。テスト自身の短い timeout assertion は変更しない。

公開説明は既存 `docs/REFERENCE/web-app.md` と Session/export 文書 owner へ集約する。新しい独立マニュアルを増やさない。実際にユーザー操作手順や screenshot を更新するセッションは `love-me-love-my-docs`、Gallery 関連なら `web-gallery-screenshot-maintenance` の適用範囲を確認して読み、再生成証拠と目視確認を残す。`examples/gbdraw_social_preview.png` は編集しない。Python geometry を変更しない計画なので tracked reference SVG を更新しない。

各関心の承認条件がすべて証拠へ対応し、required gates、owner/path review、production/tests/docs/generated review、commit/push の一致確認が揃った時点で実装完了とする。authority merge、CI、staging、権限が未達なら項目単位で未完了を明示する。署名や文書作成を runtime 実装の完了とみなさない。
