# S07 INSTRUCTION PROMPT: 統合受入、再発防止、公開前handoff

gbdraw GUI改修と規約改正を統合検証してください。個別セッションのpassを列挙するだけでなく、実際の連続操作と公開する生成物の一致を確認します。

## 作業場所と前提

**fix/gui-feedback-remediation-20260928** と /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928を最後まで使う。clone/追加worktree/別branchを作らない。
[総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、results/S00〜S06とリポジトリ指示を読む。
originをfetchし、S01/S02の規約とcheckerがbaseにあることを確認する。必要なdev更新は同じbranchへ取り込み、変わった部分だけ再検証する。

## 実施

1. G01〜G10/R01/R02ごとに実装、契約、テスト、実画面の証拠を対応付ける。未判断の旧Sessionや未再現症状を完了へ書き換えない。
2. 固定S00入力で大規模Result後の比較切替/入力/Generate、Vnig→Linear→指定Vibrio、BGC livA→livE、Editor開閉と検索を連続実行する。
3. CWの構造検証、positive control、二つの変異負例、応答/settlement/Generate全体の予算を再確認する。必要な入力が変わらない既存passは再利用し、候補を通すための計測条件変更をしない。
4. Undo/Redo、Save/fresh Load、Export現Result、failed/canceled/stale generationの保存、source/crop/reorderとidentityを検証する。
5. 必要なGallery Session/図/tutorial/captureをowner toolで更新する。旧fixtureは保持。web-gallery-screenshot-maintenanceとlove-me-love-my-docsを実際の対象作業で適用し、完成図を目視する。
6. 公開文書の重複を増やさず既存ownerを更新する。削除されたPending/Current文字列を手順・captureが要求し続けていないことを確認する。
7. production、tests、docs、generatedを別々にレビューする。owner/path evidenceをまとめ、不要consumer、旧経路、計測専用の常設状態が残っていないことを確認する。

## 最終チェック

担当変更に適切な以下の入口を使い、実行した正確なcommand/exit/resultを記録する。

    pytest tests/ -v -m "not slow"
    ruff check gbdraw/
    npm run test:web:pr-smoke
    npm run test:web:perf-smoke

Web機能は既存functional configの対象journeyを実行する。必要な対象がPR smoke/perf profileに含まれることを確認し、未収録を実行済み扱いしない。Node unavailable時はPython Playwrightで同等journeyを検証する。実装Pythonのwheelを再準備する。

architecture/policy gateは基準側checkerでbase/headを指定して確認する。S02のfixtureだけで最終gateを代用しない。図のgeometryを意図的に変える場合だけreference更新手順を使い、期待値合わせでreferenceを再生成しない。

## 終了と外部操作

results/S07.mdに要求別結果、予算、CW感度、旧Session処遇、生成物の再現、残余リスクと除外を記載。SESSION_LOGを更新する。未完了があれば全体完了にしない。

担当差分をcommitし、英語title/summaryを提示する。push/PR/merge/deployはその時点の明示承認の範囲だけ実行する。計画書pushの承認を製品merge/deployへ転用しない。

PR文面を作成する場合はwrite-clear-pull-requestを読み、同一title/bodyに言語checkを実行する。契約・runtimeの同時更新であること、基準側例外、独立したReview責務、維持したhard coverageを説明する。remoteを再試行する前に実際の状態を確認し、CIは5分以上の間隔で確認する。
