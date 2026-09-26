# 初期状態の証拠と再現入口

このファイルは実装前の動作を説明する。実装後の受入結果は各SESSION_XX_RESULT.mdに記録する。
[総合計画書](./MASTER_PLAN.md)の対象はIssue #601のBUG-15/19である。

## Sourceと環境

- 調査source SHA: 2edc00aebc74e01003da643dfc957b513d5dcfe5。
- 計画branch base SHA: af5d942af60353dda199aa487da9152a3576b3fe。
- その間のgbdraw/tests/tools/.githubに差分がないことを確認した。static Product contractには別Issueの追加がある。
- 初期テスト環境: Node 26.8.2、Python 3.13.3、Playwright 1.61.0、Chromium。
- Python 3.10–3.12 matrix、Firefox/WebKit、公開サイト、原監査buildの受入完了ではない。

## 観測したfailureとregex

| Source | 観測 |
| --- | --- |
| services/error-normalization.js | 内部ValueError文をsummaryに表示し、traceback終端の原因も除去する |
| Worker serializeError / client deserializeWorkerError | code/stage/contextは往復しない |
| app/python-helpers.js | render exceptionをtype/str/tracebackへ文字列化 |
| runCommittedCanonicalCandidateInternal | engine-error/catchがerrorLogを設定後、status:errorのみ返す |
| similarity-alignment::applyPlan | outcome.errorがあれば表示できるが、欠落時はgeneric fallback |
| feature-editor/rule-actions.js | 任意failureをInvalid ruleと表示し、field finallyでaccepted値へ復帰 |
| app/rule-matching.js / web_support/rule_matching.py | Color/Labelは既存Python ownerへ評価委譲し、現在性確認してcommit |
| Feature Search / standalone | JS RegExp(...,'i')を使い、方言なしのInvalid regexを表示 |

[synthetic probe](./evidence/probe.mjs)と[JSON結果](./evidence/probe.json)では、
内部比較文がsummaryになり、traceback-onlyのsummaryが「stderr:」になり、
code/stage/contextがtransportで落ちた。
検索の(?i)NADH / (?P<enzyme>NADH)はInvalid regex、JSの(?<enzyme>NADH)は構文エラーなし。
検索は空catalogで構文のみ確認したため、一致件数の証拠とは区別する。

比較identity不一致を実入力から発生させた再現ではない。
Issue参照のGUI_AUDIT_DEV_20260926.mdは調査treeに存在せず、原監査の入力欄・pattern・buildは未特定。

## 初期確認の結果

| Command/範囲 | 結果 | Log |
| --- | --- | --- |
| Node error-normalization / rule-matching / similarity-alignment-actions | 41 passed | [Node](./evidence/node-baseline.log) |
| Python tests/test_web_rule_matching.py | 9 passed | [Python](./evidence/python-baseline.log) |
| Chromium Python Color rule / Label TSV の二既存シナリオ | 2 passed、21.6s | [Browser](./evidence/browser-baseline.log) |

NodeのAlign unitはrunAnalysisをstubし、errorを返すため実orchestrationの欠落を覆わない。
ChromiumはPython固有構文の適用、失敗時の旧canonical保持、Undo/Redoと再生成を確認した。
未確定field draft・Copy diagnosticsは未実装で、この成功をその受入に流用しない。

## 再現コマンド

repository rootで実行する。専用clone、環境、未使用portは総合計画の共通手順に従う。

~~~bash
node docs/internal/issue-601-bug15-bug19-implementation-20260926/evidence/probe.mjs
node --test tests/web/error-normalization.test.mjs tests/web/rule-matching.test.mjs tests/web/similarity-alignment-actions.test.mjs
python -m pytest tests/test_web_rule_matching.py -q
python tools/prepare_browser_wheel.py
GBDRAW_WEB_TEST_PORT=8791 node node_modules/@playwright/test/cli.js test -c playwright.functional.config.js python-rule-parity.playwright.spec.js --grep 'Python-only color rules|Python label TSV syntax' --workers=1
~~~

wheelはgitignored test asset。cache-bustとtracked outputは変更しない。
probeはsynthetic inputだけを使い、Worker importのためselfを用意するがmessage handlerを起動しない。
実装後の新evidenceを初期baselineへ上書きしない。
