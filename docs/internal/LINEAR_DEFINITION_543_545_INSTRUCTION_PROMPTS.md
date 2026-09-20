# Linear Definition・Replicon表示修正 — INSTRUCTION PROMPTS

このファイルは、gbdrawのIssue #543と#545を修正する実装セッションの開始指示である。前提となるチャット履歴はない。

表示契約、原因、受入条件、実施記録の管理先は[総合計画書](LINEAR_DEFINITION_543_545_MASTER_PLAN.md)。各プロンプトは、指定したファイルを読めるリポジトリの作業環境で使用する。保存済みの計画ファイルが新しいブランチにない場合は、この2ファイルを利用可能にしてから開始する。

Product Decision Owner `satoshikawato`が`2026-09-19`にD1-A・D2-P・D3-Aを明示的に選択済みである。回答原文、機械可読な転記、採用した表示結果は総合計画書第4節に保存されている。S0では選択を再質問せず、最新base確認とdurable authorityへの反映を扱う。

S0→S1→S2の順に使用する。各コードブロック全体をそのセッションの指示として渡せる。S1/S2は実施記録から対象ブランチと成果物を確認するため、担当者が変わっても過去の会話を必要としない。ソース・入力・環境・受入条件が変わっていない検証を無条件に繰り返さない。

## S0 — 基準確認と確定した決定のauthority反映

```text
gbdrawのLinearラベル不具合について、基準を確認し、確定したProduct Decisionを
既存のdurable authority経路へ反映する準備をしてください。

対象:
- https://github.com/satoshikawato/gbdraw/issues/543
  Show Replicon=falseでも自動推論したプラスミド名などが残る。
- https://github.com/satoshikawato/gbdraw/issues/545
  Lock Definition Column=falseの複数レコード行でDefinitionが右寄せになる。

最初に読むファイル（リポジトリルートからの相対パス）:
1. AGENTS.md、CLAUDE.md、gbdraw/web/CLAUDE.md
2. docs/internal/LINEAR_DEFINITION_543_545_MASTER_PLAN.md
3. docs/internal/PRODUCT_IMPACT_RATCHET.md
4. docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/SESSION_COMPATIBILITY.md、docs/REFERENCE/web-app.md、docs/CLI_Reference.md

総合計画書第3節の監査基準は100ac32ad5a73cae22e7b43a2dbbc0bb2126e699。
最新baseで同じ原因が残るかを調べ、既に直った部分を重複実装する計画にしない。
確認済みの原因は、WebがReplicon名をrecord_subtitleへコピーすることと、
複数レコード行のbuilderがLock=falseでtext-anchor=endを指定することである。
この監査情報を修正後の合格証拠と扱わない。

確定済みのProduct Decision:
Concern=linear.definition-display、Scenario revision=1、
Choice=D1-A, D2-P, D3-A、Owner=satoshikawato、Decision date=2026-09-19。
D1-Aは共通幅の中央揃え、D2-Pは保存済みSubtitleとプレビューの保持、
D3-AはOrganelleを含む自動名のShow Repliconによる制御を選択する。
理由、維持範囲、廃止範囲、受け入れられたリスクの原文は総合計画書第4.3節にある。
この選択を再質問せず、同じ結果だけを反映する。

作業:
1. 作業ツリー、現在のブランチ、関連する既存変更を確認して無関係な変更を保持する。
   初回の実装作業はgit fetch originの後、最新origin/devから追跡先なしの作業
   ブランチを作る。既に開始済みなら総合計画書の実施記録から対象作業を特定する。
2. 総合計画書のF1/F2を再現する。元Issueの入力ファイルを得られない場合も、
   記載されたsynthetic recordで独立して再現できる。レコード順、個数、qualifier、
   明示Subtitle、フラグ、SVGのline-kindと文字座標を区別して記録する。
3. 確定したD1-A/D2-P/D3-Aの適用について、仕様、base側の実在する決定、
   main/releaseの互換履歴、到達するWeb/CLI/Python操作を調べる。
   明示的text_anchor設定やSubtitle継承の既存の受入範囲を確認する。
   調査は選択内容の適用確認であり、別の結果を選び直す作業ではない。
4. Product Impact Ratchetに従い、最新baseのProduct Impact mapと
   docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdを確認し、正式な反映先を
   特定する。未登録のconcern用に並行したauthority storeを作らない。
5. 第4.3節の回答に従ったauthority-only変更を準備し、生成した正式表現を
   レビュー可能にする。選択、理由、維持/廃止範囲、リスクを補完・拡張しない。
   計画内のJSON転記をCIのdecision schemaやbase反映済みauthorityとみなさない。
   必要なdurable authorityをruntimeのbaseへ反映してからS1を開始する。
   外部へのpush/mergeは別途の許可に従い、決定内容の再承認と混同しない。
6. 総合計画書第10節へbase、対象版/ブランチ、再現手順、authority反映状況、
   根拠の場所、変更可能な範囲、残る作業を記録する。最新baseと決定の衝突や
   新しいユーザー可視差が判明した場合だけ追加のProduct手続きを行う。

設計上の制約:
- SOLID: Web入力、名称解決、行振り分け、配置、描画の所有者を分ける。
- DRY: 既存のReplicon行とLinear geometryへ収束させる。
- KISS/YAGNI: 既存の小さな契約を使い、新しい表示フラグ、推測migration、
  汎用ラベル基盤、並行する描画経路を計画へ追加しない。
- 異なる明示Subtitleと同名の明示Subtitleを保存する。
- 行全体に共通する文字とレコード固有文字の現在の必要な寄与を保つ。

このセッションの境界:
再現・調査・テスト準備・authority-only変更の準備を完了する。
runtimeの修正には着手しない。新しい未決事項やauthority反映待ちがある場合は、
依存する実装だけを保留し、独立した調査を終える。
通常の読み取りや再現のたびに許可を求めない。push、PR、merge、tag、deployは行わない。

完了報告:
確認した原因、確定済み決定の正式反映状況、実装可能な範囲、対象ブランチ、
証拠へのリンク、S1が開始できるかを記載する。追加の判断や外部操作の許可が
必要な場合は、対象と規則上の理由を具体的に示す。実装修正が完了したとは報告しない。
```

## S1 — 責務の集約、修正、focused tests

```text
gbdrawのIssue #543（Replicon表示スイッチの迂回）と#545（Linear Definitionの
右寄せ）について、確定した表示契約の範囲を実装し、focused検証を完了してください。

最初に読むファイル（リポジトリルートからの相対パス）:
1. AGENTS.md、CLAUDE.md、gbdraw/web/CLAUDE.md
2. docs/internal/LINEAR_DEFINITION_543_545_MASTER_PLAN.md（全体と第10節の実施記録）
3. docs/internal/PRODUCT_IMPACT_RATCHET.md
4. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
5. 総合計画書にリンクされた表示契約、Product判断、Session互換性の根拠

開始条件:
Product Decision Owner satoshikawatoが2026-09-19に選択したD1-A・D2-P・D3-Aを
総合計画書第4節で読む。同じ選択を再質問しない。S0の対象版、作業ブランチ、
必要なauthorityがruntimeのbaseに反映されていることを確認する。
ソースや入力が変わった部分の根拠だけを更新する。計画内の受領記録とコード例だけを
durable base authorityとみなさない。必要な反映が未完了ならS0の残る作業を具体化し、
独立して進められる調査・テスト準備を行う。新しいProduct差を自分で選ばない。
既存作業は由来を確認したdev起点の作業ブランチで引き継ぎ、無関係な変更を保持する。

実装:
1. 総合計画書第6節のR/A/P/X受入条件に対応する失敗例を既存test所有者へ追加する。
   バグのある版での失敗と、既存の正しい挙動を区別する。testが実装の式を
   コピーしているだけにならないよう、利用者の操作結果とSVGをassertする。
2. app-setup.jsのexpandDiscoveredLinearRecords()から、推論Replicon名を
   record_subtitleへ自動コピーする経路を取り除く。D3のOrganelle方針を満たす。
   自動名は既存のPython metadata/DefinitionGroupのReplicon行で扱う。
3. 参照検索で利用者がなくなったdiscoveryのSubtitle情報、判定フラグ、専用処理を
   同じ変更内で削除する。GenBankの軽量読み取りとGFF/FASTAのhelper経路を確認し、
   Organism/strain推論や実利用のあるコア関数は保つ。
4. 明示Subtitle、ファイル既定値、保存済みSessionはD2に従って扱う。
   Show RepliconでSubtitle全体を隠したり、文字列一致で保存値を削除したりしない。
5. 行/レコードの振り分けはassemble.pyの_split_definition_line_kinds()に保つ。
   builderに分割条件を再実装せず、先頭だけの行識別や後続の固有ラベルを欠落させない。
6. D1に従って横配置を修正する。anchorだけの変更で終わらず、基準幅、行offset、
   gap、描画境界、衝突境界を整合させる。既存のLinear geometry/layout境界の
   小さな配置結果を共有し、builders.pyと_record_collision_bands()の重複式を削る。
   単一/複数/混在行で同じ横配置契約を満たす。垂直配置の必要な挙動は保つ。
7. 総合計画書第8節のfocused checks、追加した回帰ケース、到達するrequest/
   Session/CLI/Circular testsを実行し、in-scopeの失敗を修正する。
8. productionとtestsの差分を別々に読む。責務と経路の前後、削除した旧経路、
   挙動検証、適用gate、rollbackを簡潔に記録する。例外条件に当たる場合だけ
   完全なOE/PE/CB証拠と必要な判断を準備し、checkerの権威を緩めない。

設計制約:
SOLIDは新しい継承階層を足す理由にしない。入力、意味解決、配置、描画の所有者を
明確にし、builderは必要な解決済み入力とgeometryに依存させる。
DRYは自動Repliconの二重管理と配置の重複計算をなくすことで満たす。
KISS/YAGNIに従い、別のwatcher、Showフラグ、Sessionスキーマ、policy engine、
DOM重複除去、代替生成経路を追加しない。既存のcanonical requestとWorkerを使う。

このセッションの境界:
runtimeとtestsの整合、必要なfocused checksと修正まで完了する。
公開例や参照SVGの更新はS2の生成・目視確認に結び付ける。通常テストで参照を上書き
せず、test-owned timeoutや受入条件を緩めない。無関係な全面リファクタを含めない。
push、PR、merge、tag、deployは行わない。

引継ぎ:
総合計画書第10節に対象版/ブランチ、変更ファイル、満たした受入ID、正確な検証
コマンドと結果、未実施項目、S2への注意点を記載する。未commitの差分も識別する。
current sourceのbrowser wheelによる最終検証をまだ行っていない場合は明記する。
英語のproposed commit titleと短いsummaryを示す。S2が残る間は修正全体を
完了扱いしない。
```

## S2 — 実ブラウザ、文書・図、最終受入

```text
gbdrawのIssue #543と#545の修正について、実ブラウザの生成から保存・再読込までを
検証し、必要な文書と図を整備して最終受入を完了してください。

問題の概要:
#543は自動Replicon名が通常のSubtitleにも書き込まれ、表示スイッチを迂回する問題。
#545はLinearの複数レコード行でLock Definition Column=falseが右寄せになる問題。
修正の目的は、Replicon表示を一つの経路で扱い、Definition配置の決定を描画と
衝突判定で共有することである。

最初に読むファイル（リポジトリルートからの相対パス）:
1. AGENTS.md、CLAUDE.md、gbdraw/web/CLAUDE.md
2. docs/internal/LINEAR_DEFINITION_543_545_MASTER_PLAN.md（全体と第10節の実施記録）
3. docs/internal/PRODUCT_IMPACT_RATCHET.md
4. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
5. 実施記録から参照されるS0の判断とS1の対象差分・検証結果

開始条件:
Product Decision Owner satoshikawatoが2026-09-19に選択したD1-A・D2-P・D3-Aを
総合計画書第4節で確認し、同じ選択を再質問しない。
S1の作業ブランチ、ソース、未commit変更、正式に反映されたauthorityを確認する。
承認や証拠を別candidateへ無条件に流用しない。S1未完了なら不足したローカル修正と
focused検証を先に完了する。未解決のProduct判断は依存箇所だけ停止して明示する。

作業:
1. 総合計画書第8節のPlaywright環境確認を行い、現在のPythonソースから
   python tools/prepare_browser_wheel.pyでbrowser wheelを準備する。
   wheelとソースの版・変更対象モジュールの一致を記録し、実際にそのwheelを
   使用したブラウザ生成を確認する。既存wheelのファイル名だけで同一性を判断しない。
2. 第6節の受入表を実行する。特に、複数GenBankのupload→off/on/off各Generate、
   先頭plasmid、手入力Subtitle、Line Styles、Lock両値、長短ラベル、行offset、
   Save→新しいページでLoad→Generate、対応する旧Sessionと明示クリアを確認する。
   各操作の完了を待ち、古いResultを新しい生成結果として測定しない。
3. 最終SVGのdata-definition-line-kind、文字内容と数、文字座標、gap、
   描画/衝突境界、キャンバス内の収まりを確認する。属性だけで合格にしない。
   選択された仕様の生成SVGと保存/再生成の結果を比較する。
   NodeのPlaywrightが使えなければPython Playwrightで等価の実操作を行う。
   ChromiumのOperation not permittedは同じローカルcheckを必要なsandbox権限で
   再実行して確認する。mockだけでPyodideの全経路を検証済みにしない。
4. 共通Python/JSへの変更が到達するCLI、typed Python、Circularを必要な範囲で
   確認する。有効なS1の証拠は再利用し、変更・失敗・未解決箇所だけを再検証する。
   現在のCI impact判定とProduct/architecture policyが要求するgateを完了する。
5. Web仕様、CLI help/reference、Session互換説明の該当箇所を選択された挙動へ
   更新する。保存済み自動Subtitleが残りうる制約と明示的な対処を正確に記述する。
   必要な既存文書を編集し、説明ページや仕様管理先を増やさない。
6. 影響する実用的なmulti-record図を既存recipe/Sessionから生成し、読める倍率で
   目視する。必要なラベル、legend、track、色、比較情報を保つ。
   Gallery説明/captureを変更するときだけ対応skillを読む。最小テスト図を
   公開例に使わず、無関係なGallery全体を刷新しない。
7. 参照SVGは意図した幾何変更を確認後、必要なものだけ正式な更新手順で生成する。
   SVG差分を読み、TestOutputComparisonを再実行する。browser wheelはcommitせず、
   examples/gbdraw_social_preview.pngには触れない。
8. production、tests、docs、生成物の差分を別々に確認する。検証で見つかった
   in-scopeの不具合を直し、その変更で無効になった証拠を更新する。
   計画の受入条件を、結果に合わせて弱めない。

SOLID/KISS/DRY/YAGNIの最終確認:
責務の所有者と依存方向が明確で、RepliconをSubtitleへ重複保存していないこと。
描画と衝突判定が同じ配置結果を使い、単一/複数レコードの契約が一致すること。
不要になったpayloadや判断経路を削り、将来のためだけの抽象、状態、互換分岐、
管理工程を追加していないこと。例外が必要な場合は既存policyに従うこと。

完了条件と報告:
総合計画書第9節を確認し、第10節へソース版、ブランチ、wheel対応、検証コマンド、
受入IDごとの結果、生成図、制約、rollbackを記録する。未実施・失敗を明記する。
全条件成立時だけ実装完了とし、英語のproposed commit titleと短いsummary、
変更ファイルと証拠へのリンクを示す。
push、PR、merge、tag、deployや外部への投稿は、このプロンプトでは許可されない。
別途許可されている場合はその対象と条件を確認し、該当する既存手順に従う。
```
