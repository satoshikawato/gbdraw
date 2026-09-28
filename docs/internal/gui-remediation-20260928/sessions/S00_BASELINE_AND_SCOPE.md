# S00 INSTRUCTION PROMPT: 基準、再現、判断範囲

gbdraw Web の応答遅延、モード設定混入、BGC alignment、検索/Editor配置の改修に必要な証拠を採取してください。製品runtimeは変更しないセッションです。

## 作業場所と必読

- ブランチは **fix/gui-feedback-remediation-20260928**。
- /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928 を再利用する。clone、追加worktree、別ブランチを作らない。
- [総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、checkoutのAGENTS.md、CLAUDE.md、gbdraw/web/CLAUDE.mdを読む。
- 総合計画の開始/終了、権限、検証運用を適用する。未追跡の別提案書や会話履歴を前提にしない。

## 実施

1. originをfetchし、branch/upstream/dirty stateを確認する。初期devは57cef3ba47f4b7790a9145f2ce6988422a00710e、比較mainは4556e04e929a4a85ad28d1833ce7304bd764881c。更新差分があればsymbol単位で再確認する。
2. 総合計画第9節の入力を固定し、hashと実行環境を記録する。baselineはsource不変のsnapshotとして一つの再利用場所へ抽出し、そのsourceから生成wheelを準備する。Session/wheel/sourceを混ぜた比較をしない。
3. MG1655/Sakaiを用い、生成前/生成後/override有無でNo comparison→LOSAT、LOSATN↔LOSATP、通常入力、Generateを計測する。Vnig保存preview後も測る。端点・runner・除外条件・時間予算は候補実装前に固定する。Historyと遅延reactive処理を含む実行可能なsettlement判定を用意し、各操作が本当に選択を変更したことをassertする。cold/warmはWorker construction/initialization/run countersで確認する。
4. main v41/dev v44のVnigをそれぞれLoad→Linear→GCF_000196095.1とGCF_000354175.2を読み込み→Generate。title、row grouping、Accession/Lengthのintentと実効値、Repliconを各段階で記録する。期待値のためにfixtureを修正しない。
5. BGC Sessionからcomparison-canonical-orthogroups-1をdecodeし、og_18/memberByProteinIdを使うresolver proofを再現する。同じ候補のedge無/有でambiguity=1/0、後者はracM/unique_direct_rbh。fixture/resolver blob IDとコマンドを残す。
6. BGCでfresh Load→livA Keep→livE Keep、fresh Load→livA Review/All right→livE Keepを実操作する。修正前のlivAが未選択Reviewを開く場合は、racMを明示選択→KeepまたはAll right→Applyとして次へ進む。この補助操作はbaseline限定であり、候補の通常AlignにはracM自動解決を要求する。record reverseComplement、region reverse、表示矢印、Result/Historyを記録し、新規反転と保持を区別する。
7. 契約の差分を棚卸しする。PD-OI-024/037/054、inactive draft互換、PD-OI-038/039の狭幅保護を確認する。BD番号はbaseに実在するものだけ引用する。
8. CW-01〜06の計測可能な境界、positive control、未計測失敗条件を定義する。既存runtime-test-hooks/contract helperを使う見込みを記録し、generic registryを作らない。

## 不足するProduct判断の扱い

Pending削除、help-tip、検索の縮小とEditor退避、Align右側のReview、契約/実装同時更新、二重計算防止は確定要件。同じ選択を再質問しない。

旧flat Sessionのinactive値については、意図的な保存と生成由来を区別できる証拠があるか調べる。なければ、既存明示値を保持する結果とactive-modeだけ復元する結果の保存互換・旧Vnigへの影響を具体的に比較する。後者の適法性を既存互換契約で評価し、判断を代行しない。per-modeにするtitle/font等の範囲も明示する。

必要なDecision Packは既存テンプレートを使用する。新しいdecision frameworkは作らない。未確定はその値を変更する処理だけに限定し、独立した調査を完了する。

## 成果物と終了条件

- results/S00.md: コマンド、hash、環境、観測、性能raw dataへの参照、再現限界。
- 必要な最小のreproduction scriptと計測結果。scriptは結果Markdown内に完全な実行コードとして保存できる。raw genome/巨大profileを無差別に複製しない。S01のdocs-only差分に独立した実行コードが混在する場合は、先に証拠のみの取り込みを済ませ、規約変更へ混ぜない。
- 結果文書内のProduct判断表、規約/CI改正対象、CW検証計画。
- SESSION_LOG更新と担当差分のcommit。英語title/summaryを記録する。ブラウザのconfig、test選択、wheel準備を含む実際のコマンドを残す。

S00は証拠採取の完了であり、製品の修正完了ではない。必要な承認が未取得なら不足フィールドと理由を明記する。外部公開はその時点の承認範囲で行う。
