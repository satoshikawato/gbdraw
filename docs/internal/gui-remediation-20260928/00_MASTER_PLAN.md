# Web GUI 改修・計算重複防止 総合実装計画

Status: 計画。製品実装、規約改正、CI 改正の完了を表さない。
作成日: 2026-09-28。
対象: gbdraw のブラウザ版。読者は、この文書から参加する実装者・レビュー担当者。

## 1. 目的と作業場所

操作応答の遅延、Circular から Linear への設定混入、Preview の検索と Editor の配置、Similarity alignment の情報欠落と内部 ID 表示を修正する。不要な常時説明を減らし、重い計算の重複を再発させない契約と検証を設ける。

製品契約と実装を同じ PR で更新できる規約・CI 改正も対象とする。規約改正を利用した製品変更は、改正済みの基準ブランチで検査する。

**すべてのセッションは fix/gui-feedback-remediation-20260928 ブランチを使用する。**

| 項目 | 値 |
| --- | --- |
| 元リポジトリ | /mnt/c/Users/genom/GitHub/gbdraw |
| 再利用する worktree | /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928 |
| remote | origin / https://github.com/satoshikawato/gbdraw.git |
| push 先 | origin の同名 fix/gui-feedback-remediation-20260928 |
| 作成時の最新 origin/dev | 57cef3ba47f4b7790a9145f2ce6988422a00710e |
| 比較用 origin/main | 4556e04e929a4a85ad28d1833ce7304bd764881c |
| 計画入口 | このファイル |
| 進捗の記録先 | [SESSION_LOG.md](SESSION_LOG.md) |

この worktree は既存の未コミット変更を隔離するために一度だけ作成した。クローンは作成していない。以後のセッションで clone、追加 worktree、セッション別ブランチを作らない。別の機械へ移る場合も既存 checkout の有無を確認し、必要な取得は初回に一度だけ行い、その場所を SESSION_LOG に記録して再利用する。

通常の「新しい作業は dev から新規ブランチ」という一般則に対し、この計画の継続セッションでは指定済みブランチを再利用する。main/dev への直接 commit/push、公開済み履歴の force-push、無関係な変更の破棄は禁止する。

## 2. 読む順序と用語

1. リポジトリの [AGENTS.md](../../../AGENTS.md)、[CLAUDE.md](../../../CLAUDE.md)、[Web CLAUDE.md](../../../gbdraw/web/CLAUDE.md)。
2. 本計画、SESSION_LOG、実行対象セッションのプロンプト。
3. [Architecture Fitness Function Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)、[Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md)、[Web Change Policy](../WEB_CHANGE_POLICY.md)。
4. 対象機能の [Product Contract](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md) と [Session Compatibility](../../SESSION_COMPATIBILITY.md)。

Result は最後に成功した描画結果、draft は次回 Generate 用の設定、canonical request は保存・描画に共通の型付き要求を指す。owner は値の意味・妥当性・変更規則を決める唯一の責任箇所。RBH は reciprocal best hit。画面上の名前と、照合に使う安定 ID は別物である。

本計画は実行順序と受入条件を定める。未取得の人間の決定、未実行テスト、未マージの規約を承認済みとして扱わない。

## 3. 要求と観測の対応

| ID | 利用者に必要な結果 | 確認済みの原因・未確定事項 | 担当 |
| --- | --- | --- | --- |
| G01 | No comparison → LOSAT、LOSATN ↔ LOSATP、通常入力がすぐ反応する | dev の Pending 派生表示が生成 intent とラベル表を繰り返し構築。override が空でも全 feature を走査 | S03 |
| G02 | Generate の押下をすぐ受け付け、処理中・取消・失敗を正しく示す | 外側の record 準備・検証が Processing 設定より先に実行される | S03 |
| G03 | Circular Session 後の Linear にタイトル・行配置・表示設定が混入しない | title は共有。欠落した Linear layout を false に復元。保存済みの非 active Show も存在。Replicon ON は未再現 | S00、S04 |
| G04 | 検索バーを以前の大きさに戻し、Editor 開時に退避する | main は最大 39.5rem。dev の専用全幅行への変更で拡大 | S06 |
| G05 | 通常幅では Editor が Preview 上端から開き、検索バーの下から始まらない | dev の Editor containing block が検索行の下の workspace | S06 |
| G06 | livA Align で racM を正しく選択し、自動 Align する | 保存済み直接 RBH が UI catalog で省略され、resolver に渡らない | S05 |
| G07 | livE Align の方向変化を説明でき、Keep が勝手に反転しない | 現行既定は Keep。新しい反転か既存方向の保持かはブラウザで未確定 | S00、S05 |
| G08 | Review alignment options を Align の右隣へ置く | popup と drawer の markup が対象 | S06 |
| G09 | 指定された長文を削除、または必要時のヘルプにする | 対象と削除境界は第 8 節 | S03、S06 |
| G10 | Similarity group に配列・feature の名前/IDを示す | inspector が recordKey と biologicalFeatureId を直接表示 | S05 |
| R01 | Product 契約・実装・必要な通常テストを同じ PR で更新できる | 現行規約と checker の双方に同時変更の禁止がある | S01、S02 |
| R02 | 同じ操作・同じ入力の計算重複が再発したら検証で検出する | 既存の構造メトリクスを拡張可能。未計測をゼロとして扱わない検証が必要 | S01、S03、S07 |

Pending の指定表示削除、help-tip 化、検索幅と退避、Align/Review 横並びは確定した要求である。同じ UI 選択を再質問しない。旧 Session の曖昧な非 active 設定の扱い、per-mode 化する表示設定の正確な範囲、未再現症状だけを限定した判断対象とする。

## 4. 調査根拠とその限界

調査時のコード根拠は第 1 節の dev SHA。行番号は変更されるため、以下の symbol を検索する。

| 根拠 | 場所・symbol |
| --- | --- |
| Pending 再計算 | gbdraw/web/js/app/app-setup.js / generationApplicationFeedback |
| intent と生成表の重複 | gbdraw/web/js/services/config.js / getGenerationApplicationStatus、services/session-request.js / projectGenerationIntent、compareGenerationIntent、addGeneratedTableResources |
| 空 override でも全件走査 | gbdraw/web/js/app/feature-editor/label-override-table.js / buildLabelOverrideRows |
| Generate 前処理 | app/app-setup.js / runAnalysis、app/run-analysis.js の processing と afterPaint |
| mode 復元・既定値 | services/config.js / applyConfigData、services/session-request.js、services/session-active-config-contract.js、js/mode-profiles.js、app/layout-preferences.js |
| Preview 配置 | gbdraw/web/index.html / preview-feature-search、preview-editor-layout、preview-workspace、right-drawer |
| RBH 情報の省略 | gbdraw/web_support/feature_catalog.py / orthologEdges の件数化 |
| Web alignment 入力 | gbdraw/web/js/app/similarity-alignment.js / group.orthologEdges、inspectActivePlan |
| Python の一意直接 RBH 解決 | gbdraw/layout/similarity_alignment.py / resolve_similarity_alignment |
| canonical graph の取得候補 | gbdraw/web_support/similarity_alignment.py / _projection_context、gbdraw/api/record_planning.py、gbdraw/session_request_codec.py |
| 構造計測 | gbdraw/web/js/services/runtime-test-hooks.js、tests/web/helpers/session-regeneration-contract.cjs |

調査用の予備計測には大規模 Result 後の約 1 秒の long task があるが、配信 SHA、wheel、機械条件、反復統計が不足している。この値を性能保証や合格証拠にしない。S00 で基準を採取する。

BGC の純粋 resolver 検証では、同じ og_18 candidates に対し、edge 無しでは曖昧性 1 件、保存済み edge 有りでは 0 件、racM が aligned / unique_direct_rbh となった。これは domain 層の証拠であり、Worker・browser 経路の合格証拠ではない。S00/S05 で再現可能な記録を残す。

Vnig の main artifact は Session v41、dev artifact は v44。両方とも Replicon は false。Gallery を更新しただけで、利用者が保存している旧 artifact の問題まで解消したと扱わない。

## 5. アーキテクチャと設計制約

- UI は既存 action を呼び、設定意味は既存 mode/session owner、alignment の生物学的判断は Python、画面配置は CSS が決める。責務の重複を作らない。
- canonical request → 既存 diagram Worker → typed validation/planning/render → sanitizer → Result/History の経路を維持する。第二の request builder、第二の Python runtime、新しい LOSAT 経路を作らない。
- 公開 API、typed request、CLI と Web は同じ domain 意味を維持する。共有化は実在する二経路の重複を消す場合に限定し、旧経路を同じ変更で除く。
- 呼び出し元には必要な入力・結果だけを公開する。UI の status が genome graph、全 metadata、SVG checkpoint を要求する interface にしない。
- before/after、Save/Load、cold/warm、失敗/取消でも同じ契約が成立する。検証や安全処理を省く最適化をしない。
- 状態を増やす前に不要 consumer と重複処理を除去する。新 cache、汎用 scope framework、whole-app snapshot、Vue build 変更は初期スコープ外。
- SOLID は責務と依存境界、DRY は同じ意味・計算の単一 owner、KISS は少ない経路と状態、YAGNI は実測で不要な一般化を追加しないこととして適用する。短いコードだけを目的にしない。

Architecture Ratchet の owner/path evidence を残す。通常の非増加変更は簡潔な before/after owner・path と削除した旧経路を記録し、全 OE/PE/CB 集合は既存規約の例外条件に該当するときだけ作成する。

## 6. R01: Product 契約と実装の同時更新

### 改正後の通常経路

対象は静的 Product 契約 docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md の変更と、その内容を実現する runtime・関連する通常 tests/docs の同一 PR。契約は guard inventory に残し、変更時は Review REQUIRED とする。

変更の結果、保持する結果、退役範囲、既知のリスクを既存の Product Decision 形式で記録し、人間が契約と実装の一致をレビューする。明示済みの選択は利用し、不足する項目だけ確認する。新しい署名基盤、decision registry、承認 UI、Markdown 全体の独自意味 parser は作らない。

次の境界は維持する。

- 判定器・detector・workflow・authority map・maintainer allowlist は基準ブランチのものを使う。
- static 契約の変更が、機械評価対象の map/BD/architecture rule を同時に書き換えて自分を正当化する経路を認めない。
- candidate 契約の存在だけで自動 CONFORMING にしない。静的な製品判断は人間の Review、既存の機械判定は基準ルールに従う。
- base map に登録された hard contract を変更した候補テストだけで合格させない。必要な unchanged PR_GATE coverage がない場合は既存の限定手順を使う。
- 通常 runtime-only、contract-only、既存 mapped decision の経路は維持する。例外を他の governance ファイルへ一般化しない。

現在の checker には Product Contract と他ファイルの同時変更、および runtime と guard の同時変更を拒否する二つの判定がある。片方だけ外して完了としない。

### 最初の規約移行だけに必要な順序

現行 checker は checker 実装と authority 文書の同時変更も拒否する。そのため、現在の基準ルールで通る次の二段階を、同じ作業ブランチで実施する。

1. S01: 規約文書だけで限定例外を事前定義し、dev へ取り込む。
2. S02: checker とそのテストだけで例外を実装し、dev へ取り込む。
3. S03 以降: 改正済み dev を同じ作業ブランチへ取り込み、Product 契約・実装・通常テストを同一変更で更新する。

S01で移行規則を明示する。S02有効後は、新しいbase規則が旧static Contractのlifecycleや個別receiptにある「契約を先にmerge」の手続き要件だけを置き換える。製品結果・保持条件・互換性の意味は変えない。古い手続き文は最初の許可された製品変更で整合させる。

S01/S02 の merge は別途明示された外部操作の承認に従う。承認済みでなければ完成した差分・検証・PR案まで準備する。今後の製品変更ごとに契約専用 PR を要求する運用は残さない。移行自身を未マージの例外で通そうとしない。

## 7. R02: 計算重複の再発防止契約

規範は既存 Architecture Ratchet に追記し、Web CLAUDE には具体 owner と参照を置く。別の汎用計算レジストリは追加しない。

初期の自動検証対象は、生成用ラベル override 表の直列化と feature selector metadata/index、および比較切替・通常入力・Generate受付の経路に限定する。ほかの計算は設計レビューの原則として扱い、全アプリケーションの計算を追跡する仕組みは作らない。

| 契約 ID | 規則 | 検証 |
| --- | --- | --- |
| CW-01 | 操作受付、選択表示、help、status の更新だけでは canonical request、生成 TSV、全 feature metadata、SVG/checkpoint の構築・走査を開始しない | 実際の UI 操作を通した観測。比較切替の Python/LOSAT dispatch も 0 |
| CW-02 | 対象表/indexの同じ構築目的・同じ操作・同じ意味入力では、必要な構築は高々1回。実在する複数consumerはownerの結果を使う | 対象の構築境界と入力版を観測し、重複consumerの退行を検出 |
| CW-03 | ラベル・可視性 override がすべて空なら、override 出力のために feature metadata/index を作らない | metadata 入力へアクセスすると失敗する sentinel と正常な空結果を検証 |
| CW-04 | required validation、sanitize、各 Generate の描画・catalog 作成を省略しない | 異なる入力、別 Generate、retry、invalid request の負例と既存出力契約 |
| CW-05 | 表示専用 consumer を除去するとき、そのためだけの producer、subscription、保存状態も同じ変更で除去する | consumer inventory と意味のある browser regression。文字列検索だけを合格証拠にしない |
| CW-06 | 再利用の範囲・入力依存・寿命を owner が説明できる | 操作ローカル再利用を優先。長寿命 cache は実測根拠と無効化検証が必要 |

CW-02 の「同一」は関数名や配列ポインタの一致では決めない。入力データ・record/source identity・override・mode 等の意味が同じことを owner が定義する。独立した before/after 観測は、偶然同じ値でも目的が異なるため違反ではない。新しい操作、明示 retry、正当な validation と render も重複違反ではない。「一度計算したら二度と計算しない」という契約ではない。

既存 runtime-test-hooks と performance/contract helpers を使う。計測点は実処理の境界に置き、定数 0 の記録で仕事の不存在を証明しない。probe が動いたこと、必要な正の処理を観測できることを確認する。未記録・計測無効・coverage 欠落は合格にしない。

少なくとも、空 override の早期終了を除去する変異と、status から全 feature 走査を再導入する変異を使い、対象テストが失敗することを確かめる。変異は一時パッチとして実行・復元し、通常差分へ残さない。新しい metric 名が必要でも、観測のためだけの production state manager は作らない。意味入力の識別には既存のoperation/source revisionを使い、全genomeのhashや新しい汎用dependency trackerを追加しない。二つ目のconsumerが残らなければ、共有機構や仮想consumerを作らず、実在経路への重複構築変異で感度を確認する。

PR レビューでは、各重い派生値の owner、trigger、input、consumer、lifetime、廃止経路を簡潔に確認する。計算回数と応答時間の両方を検査し、高速な機械で重複が隠れることを防ぐ。

## 8. 製品変更の具体範囲

### 応答と説明文

S03 で Pending 専用 generationApplicationFeedback/getGenerationApplicationStatus と、その専用 intent baseline/bookkeeping を consumer 追跡の上で除去する。canonical request builder、compareCanonicalRenderRequests、committedCanonicalSession、artifact capture/restore は維持する。

空 override の早期終了を metadata 構築前へ置く。Generate の既存 operation owner を前処理まで広げ、最初に受付状態を公開・描画可能にする。二重クリック、前処理中 Cancel、遅着 helper、validation が review を開く場合、失敗・再試行を同じ lifecycle で処理する。

| UI 文言・領域 | 改修 |
| --- | --- |
| Result 上の Pending、unknown、Supported color…、Save stores… の常時説明 | 常時ブロックと不要な派生処理を削除。実エラーは既存 operation ref から必要時に通知 |
| Generate 上の Pending と Generate recalculates placement… | 削除。ボタン内 Processing/Canceling、実際の error/recovery は保持 |
| Align applies resolved targets immediately… | S06 で既存 help-tip へ移動 |
| ON aligns definitions in a… | S06 で既存 help-tip と統合し、常時段落を削除 |
| Current: Run LOSAT… とプログラム・閾値の要約 | S06 で削除。選択済み状態、Custom、実際のエラーは必要時に保持 |

tooltip は hover・keyboard focus・tap に対応する。説明は一つの text source から可視 tooltip と accessible description へ提供する。隠しテキストのために重い computed を存続させない。指定外の全 help/警告を一括削除しない。

### モード設定と Session

欠落、明示 false、明示 Show、Auto を区別する。Circular request に存在しない Linear フィールドを projection が発明せず、fresh Linear defaults を維持する。タイトルと合意した表示設定を既存の per-mode owner に統合する。

未訪問 Linear の受入値は title 空、Arrange in rows=true、Replicon=false、Accession/Length=auto。Auto の結果は実際の共有行 topology に従い、常に非表示とはしない。モードを往復したら、そのモードで編集した値へ戻す。

両モードを編集して保存した inactive draft の round trip は既存契約。旧 flat Session の inactive Show が生成由来か手動設定かは field の存在、既定値一致、filename、source の不存在から判定できない。S00 で旧値の保持/初期化の範囲を明確にし、S04 で採用された互換方針だけ実装する。保存済み旧 Vnig の扱いが未決定なら、G03 全体を完了にしない。

Gallery publication の未使用モード既定値を修正する。生成物は owner tool で再生成し、旧 v41/v44 fixture を修正後のデータに置き換えて検証をすり抜けない。

### Similarity alignment と名前

BGC Session の comparison-canonical-orthogroups-1 resource にある og_18 の直接 edge を利用する。livA→racM は rbh、livA→racL は coortholog。bitscore の大小や representative だけで選ばない。

既存 Worker/Python adapter が canonical resources と provenance から選択グループの edge を取り出し、既存 resolver に渡す。CLI と Web の実際の anchor/edge projection の重複だけ共有する。UI の reactive catalog へ graph 全体を戻さず、LOSAT を再実行しない。

保存時 recordIndex を現在の表示順序と同一視しない。recordKey/source/protein/feature の binding、crop、reorder、重複 source を検証する。graph が欠ける場合に推測せず、既存の一候補解決と本当の曖昧性 Review を維持する。新 ambiguity policy、rationale enum、plan schema は不要。

通常 Align は既存 Keep を維持する。livA Keep→livE Keep、および livA の明示 All right→livE Keep の両操作を検査する。

record は表示名と accession/record ID、feature は gene と locus_tag/protein ID を示す。欠落時は feature type と座標などで識別可能にし、hash は内部 binding に保持する。review と inspector の既存 display helper を共有し、表示名を内部キーにしない。

### Preview

通常幅で検索最大幅は main の 39.5rem を基準とし、利用可能幅に収める。Editor は Preview 上端から開き、検索を残り幅へ退避させる。drawer 幅を同一 CSS 変数で共有し、固定 360px の JS 移動や新座標 ref は追加しない。

コミット全体の revert はせず、検索/Editor の配置を局所修正する。検索の自由 drag 復活は要求に含めない。狭い画面での Editor/review の canvas 確保、スクロール、Close、keyboard 操作は維持する。Align と Review は popup と drawer の両方で同一行、狭い幅では切断せず折返し可能とする。

## 9. 入力、性能基準、合格条件

固定する入力:

- tests/test_inputs/MG1655.gbk と Sakai.gbk。
- tests/test_inputs/AP027131.gb と AP027132.gb。
- gbdraw/web/gallery/sessions/Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz の main v41 と dev v44。
- tests/test_inputs/GCF_000196095.1_ASM19609v1_genomic.gbff と GCF_000354175.2_ASM35417v2_genomic.gbff。
- gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json。

S00 は入力 SHA-256、配信 source SHA/dirty status、wheel hash、browser/Playwright version、CPU/機械、viewport、cross-origin isolation を記録する。main と dev で Session 互換が違う場合は、その比較を unavailable と明記し片側だけ暗黙 migration しない。baseline はsourceを変更しないsnapshotとして一つの再利用ディレクトリへgit archive等で抽出し、追加clone/worktreeを作らない。各snapshot自身のsourceからgitignored browser wheel等の生成assetを準備し、hashを記録する。作業branchのwheelをbaselineへ流用しない。

性能は生成前、生成後、保存 preview 読込後、override 無/有、cold/warm Worker を分ける。warm 操作は宣言した warm-up 後、複数 fresh page にまたがる各 20 回以上を測る。実 pointer/keyboard event から最初の可視反映までに blur/History を含め、その後の UI settlement と Generate 全体時間も記録する。handler だけの時間、spinner だけの改善で合格させない。S00でHistoryと遅延reactive処理まで含む実行可能なsettlement判定を定義し、操作が実際に選択状態を変えたことをassertする。cold/warmはfresh pageという名称だけで決めず、Workerのconstruction/initialization/run観測で証明する。

固定 runner 上の提案予算は、可視反映 p95 <=100ms/max <=250ms、比較切替 settlement p95 <=200ms、切替に起因する main-thread task >=100ms がないこと。S00 で runner と測定端点・除外条件・Generate 完了時間の許容差を候補測定前に固定する。達成不能なら根拠を記録して対策・正式な予算判断を求め、候補を通すために入力縮小・後付け除外・閾値緩和をしない。

| 受入 | 必須内容 |
| --- | --- |
| 構造と応答 | CW-01〜06、空/非空 override、mode 切替で不要 search/render がない、非空編集も残る |
| Generate | 早期受付、重複押下1操作、前処理 Cancel、遅着拒否、review 終了、失敗復旧、retry |
| Session | 旧 Vnig と現行、指定2入力、両方向 mode 切替、明示 inactive 設定、Save/fresh Load、import failure rollback |
| BGC | livA→racM direct RBH、parA→racL、genuine ambiguity Review、Keep方向不変、明示方向変更、Undo/Redo/Reset/Save/Load |
| UI | 1440×900、1024×740、768×740、390×844、390×740、200% zoom。検索/Editor非重複、canvas/Close/toolbar到達、focus/tap |
| 維持する意味 | live edit、Result/draft分離、Exportは現Result、canonical validation、sanitizer、source privacy、cancel/stale隔離 |
| 改正ルール | 許可する同時変更、禁止guard混在、Review REQUIRED、通常経路、mapped hard coverageをfixtureで確認 |

対象 unit/adapter/browser テストを先に実行する。最終段階で必要な fast Python、Web機能、performance、architecture gates を実行する。production、tests、docs、generated diff を別々にレビューし、変更がない証拠を理由なく再採取しない。

## 10. セッションと依存関係

| 順序 | プロンプト | 完了条件 |
| --- | --- | --- |
| S00 | [基準・再現・判断範囲](sessions/S00_BASELINE_AND_SCOPE.md) | 固定入力、基準計測、契約差分と未確定事項を記録 |
| S01 | [規約の事前改正](sessions/S01_GOVERNANCE_POLICY.md) | 同時更新例外とCW契約を規範へ記載、docs-onlyをdevへ反映 |
| S02 | [CIの同時更新許可](sessions/S02_GOVERNANCE_CHECKER.md) | 改正規範下でchecker-onlyを検証、devへ反映 |
| S03 | [応答改善と二重計算防止](sessions/S03_RESPONSIVENESS.md) | status除去、Generate lifecycle、CW退行検出 |
| S04 | [モード別設定とSession](sessions/S04_MODE_SESSION.md) | 合意した旧形式方針を含むG03受入 |
| S05 | [Alignmentと名前](sessions/S05_ALIGNMENT_IDENTITY.md) | domain/adapter/実browserの証拠が一致 |
| S06 | [検索・Editor・ヘルプ](sessions/S06_PREVIEW_AND_HELP.md) | 配置、説明、アクセシビリティの受入 |
| S07 | [統合検証・生成物・handoff](sessions/S07_INTEGRATION.md) | 全要求の結果表と差分レビュー、未完了を明記 |

同じ branch/worktree で順次行う。sub-agent は調査・レビューを並行可能だが、実装ファイルの責任範囲を明示する。同じ owner の同時編集を避け、他者の変更を戻さない。

S01/S02 の外部取り込み待ちでも、read-only 調査と独立した検証準備は進める。authority/checker-only PR の差分へ runtime を混ぜない。未成立の判断がある場合は依存する部分だけ止める。

## 11. 毎セッションの運用

開始時に指定 worktree で status、branch、upstream、前回の SESSION_LOG を確認し、origin を fetch する。別セッションが dev を更新していてもブランチを作り直さない。同期が必要なら既存履歴を保持して origin/dev を merge し、競合と新しい基準 SHA を記録する。

例:

    cd /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928
    git status --short
    git branch --show-current
    git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
    git fetch origin --no-tags

初回 push 前は upstream がないことが正常。設定する場合は origin/fix/gui-feedback-remediation-20260928 のみ。branch mismatch を自動 reset で解消しない。

ブラウザが必要なら Node/Python の両方を確認する。

    command -v playwright
    playwright --version
    python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
    node -e "console.log(require.resolve('@playwright/test'))"

Python/runtime 変更を browser で検証する前に、当該 checkout から python tools/prepare_browser_wheel.py を実行する。wheel は commit しない。sandbox Chromium 起動失敗は必要な権限で同じチェックを再試行する。test command の制限時間は少なくとも30分とし、進捗を定期確認する。CIのremote確認間隔は5分以上。

Gallery/tutorial の生成物を変更するときは web-gallery-screenshot-maintenance、手順文書や screenshot の再現作業は love-me-love-my-docs、PR wording は write-clear-pull-request を、その時点で読んで適用する。本計画作成だけで公開手順や全Galleryの再生成は行わない。

終了時は results/Sxx.md を作成し、開始/終了 SHA、変更、実行コマンドとexit結果、固定入力、残った判断、次の着手点を記録して SESSION_LOG を更新する。ログは一つの目次と各セッション結果を使い、別の並行進捗台帳を増やさない。未実行をpassと書かない。

各セッションは担当差分を一つの意味のあるcommitへまとめ、英語のtitle/summaryを残す。push/PR/merge/deployは、その時点の承認範囲でのみ実行する。計画書公開の承認を、将来のruntime公開・main/dev直接push・mergeの承認へ拡張しない。remote mutationの再試行前に実際のremote状態を確認する。

## 12. 完了の定義

G01〜G10とR01/R02がそれぞれ検証済み、または人間が明示した除外として記録されていること。旧Sessionや未再現の方向問題を、現行Galleryの更新だけで解決済みにしない。

規約とcheckerの同時更新許可が基準devに存在し、製品差分の契約とruntimeが同じPRでレビュー可能であること。CW違反を入れた負例が失敗し、通常処理・必要な再計算は通ること。

Result/History/Session/Exportの意味、科学的選択、privacy、安全境界を維持し、最終の図と画面を読める大きさで目視確認すること。実装していない項目があればこの計画全体を完了としない。
