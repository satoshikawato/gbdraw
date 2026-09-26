# INSTRUCTION PROMPT S04 — 統合検証・既存文書・完了判定

あなたは gbdraw の [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) 修正の統合担当者である。対象は装飾位置の再生成継承、Layout edit の発見、検索・toolbar 専用行の3件。個別セッションの証拠を統合し、保存・出力・性能・既存操作との組合せを確認して既存ユーザー文書を更新する。

## 取得・開始条件

S01–S03 の commit/push 完了と、同じ remote branch の writer が自分だけであることを確認する。共有作業ツリーを切り替えず、必ず公開済みの指定実装ブランチを専用 clone へ取得する。

```bash
task_dir=$(mktemp -d /tmp/gbdraw-issue599-s04.XXXXXX)
git clone --branch fix/issue-599-preview-layout-20260926 \
  https://github.com/satoshikawato/gbdraw.git "$task_dir/repo"
cd "$task_dir/repo"
git fetch origin
git pull --ff-only origin fix/issue-599-preview-layout-20260926
test "$(git branch --show-current)" = "fix/issue-599-preview-layout-20260926"
test "$(git rev-parse --abbrev-ref '@{upstream}')" = "origin/fix/issue-599-preview-layout-20260926"
git status --short
```

repository AGENTS/CLAUDE、Web CLAUDE、[総合計画](MASTER_PLAN.md)、[承認一覧と3 Pack](00_APPROVED_PRODUCT_DECISIONS.md)、[進捗](SESSION_STATUS.md)、`results/S00.md`～`S03.md` を読む。最新 dev に全必要 authority があることを確認し、入力/検査/authority SHA を記録する。足りない実装前提や失敗結果は合格とみなさない。

## 所有範囲

全受入条件の証拠対応、統合で見つかった必要な不具合修正、既存 `docs/REFERENCE/web-app.md` と Session/export owner の説明、`results/S04.md` と進捗を担当する。新しい製品結果・checker・policy・schema・public manual owner は追加しない。他 issue や前担当の有効な変更を revert しない。

## 統合タスク

1. 各 Pack の選択 outcome と must preserve/may retire/risk を総合計画 C01–C08/A01–A02/P01–P03/R01 へ照合する。独立した必須 contributions をすべて確認する。選択ID一致や一枚のscreenshotだけで完了にしない。
2. 有効な個別証拠を、対象コード・入力・環境・受入条件が不変の範囲で再利用する。統合で変わった箇所、失敗、未解決の懸念だけ再検査する。必要修正は元のowner/単一候補経路へ戻し、暫定fallbackを追加しない。
3. 実 Gallery Circular/Linear 図で装飾drag、2回Generate、search/toolbar操作、drawer開閉、Layout edit ON/OFF、色/font変更、Undo/Redo/Reset、Save/Load→Generate、current Result export の組合せを確認する。batch選択外、source/target不一致、failure/Cancel/staleのatomic保持も証拠へ対応付ける。
4. zero/nonzero/batchのparse/serialize/bind/Worker回数、追加read/hash/全checkpoint clone不在、Responsiveness guardrailを確認する。scientific coordinates/comparison/scale値と record/alignment動作を維持する。
5. CSS専用行の可視性、hit target、workspace高さ、scroll到達性、keyboard/touch/soft keyboard/zoom と drawer系列を確認する。検索state/focusと transient-free export を確認する。
6. production、tests、docs、generated のdiffを別々にreviewする。差分取得では最新origin/devとのmerge-baseと最終HEADを記録し、既にdevへ入ったauthorityをruntime新規authorityとして誤計上しない。owner/pathの収束とsuperseded paths撤去を簡潔に示す。例外条件があればratchetの完全なsets/判断を用意する。

## 必須検査と文書更新

既存 composition/candidate/run-analysis/search/browser tests と、変更に関係する保存/export tests を実行する。少なくとも既存fast Web admissionの該当検査を確認し、[Web Change Policy](../WEB_CHANGE_POLICY.md) に従い `PR / gate` と `Web base policy (trusted base)` の必要チェックを扱う。local検査とremote required statusを区別する。candidate checker/map/rules の変更で自己承認しない。

browserではNode/Python Playwrightの両方を確認し、source一致wheelと既存operation/readiness receiptを使用する。全matrix/Gallery stagingは正確なintegrated dev SHAの境界で行い、merge前の候補証拠で代用しない。CI pollingは5分以上間隔。長いtestは30分以上確保し、test自身のtimeout assertionを緩めない。

ユーザー文書には、Layout editの有効化、同じ図で手動差分が継承されること、対応不能時のReset/設定修正、絶対位置固定ではないこと、検索の専用行と自由drag退役を実画面に即して説明する。raw Python recipeだけにブラウザー手動位置再現を約束しない。既存referenceとSession/export文書に集約する。

操作手順・screenshotを更新する前に `love-me-love-my-docs` skillを読み適用をユーザーへ伝える。Gallery関連なら `web-gallery-screenshot-maintenance` skillも必要な範囲で読む。再現script・文字どおり実行した手順・regen証拠を残す。公開figureはGallery品質のrecipe/sessionを使い目視する。最小smoke図をpublic showcaseにせず、`examples/gbdraw_social_preview.png` と無関係なgenerated filesを編集しない。

Python geometryを変えない計画なのでreference SVGを更新しない。統合で科学的geometry変更が必要になれば、その根拠と新しい影響を先にレビューし、通常試験でreferenceを上書きしない。

## 完了・commit/push・権限境界

`results/S04.md` に全受入IDの証拠対応表、検査SHA/環境、コマンド・exit/result、artifact、required status、owner/pathレビュー、文書再生成結果、未完了条件を保存する。SESSION_STATUSは実際の完了範囲に合わせる。authority merge、remote gate、merge後stagingが未達なら個別に未完了を明示する。検証-onlyでも意味のある結果文書をcommitする。

branch/upstreamを再確認しfetch。予期しないremote進行は担当と差分を確認する。対象パスだけをstageして `git diff --cached --check` とstaged diff確認。**各セッションの終了条件としてcommit/pushを必ず実施する。**

```bash
git commit -m "docs: verify and document issue 599 preview layout fixes"
git push origin HEAD:refs/heads/fix/issue-599-preview-layout-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-599-preview-layout-20260926
git status --short
```

local/remote SHA一致を確認。pushエラー時は実remote stateを確認して再試行する。最終報告にbranch/commit、主要変更、検証と制限、未達の境界、英語commit titleとshort summaryを含める。PR作成・dev merge・deployはこの計画では未承認。必要なら検証済み具体的差分を提示した後で該当境界の承認を求める。PR文言を準備する場合は `write-clear-pull-request` skillとrepository validatorを適用する。

SOLID/KISS/DRY/YAGNIに従い、統合のために新しいowner/実行経路/互換branchを追加しない。失敗した境界を診断し、既存責任へ修正を集約する。
