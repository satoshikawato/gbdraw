# Issue #601 BUG-15 / BUG-19 — 総合実装計画書

Status: 製品動作の二つのA案は承認済み。runtimeの実装は未開始。
Repository: https://github.com/satoshikawato/gbdraw.git
実装ブランチ: **fix/issue-601-bug15-bug19**
作成日: 2026-09-26
作成時に取得した最新dev: af5d942af60353dda199aa487da9152a3576b3fe

新規参加者は本計画、Issue、承認記録、リポジトリから作業を開始できる。
セッション別のINSTRUCTION PROMPTを順番に実行し、結果をコミット・プッシュして引き継ぐ。

## 問題と完成後の動作

gbdrawはゲノム図を生成・編集するPythonアプリケーションで、Web版はブラウザ内で処理するSPAである。
[Issue #601](https://github.com/satoshikawato/gbdraw/issues/601)の対象は次の二点。

- BUG-15: 内部比較例外が直接表示される一方、一部のエラー伝達では原因が失われ、
  Align Applyに修正不能なgeneric errorが表示される。
- BUG-19: Python regexを使うColor/LabelとJavaScript regexを使うFeature Searchの説明が曖昧。
  rule処理の準備失敗までInvalid ruleと表示され、不正な既存rule編集値は元の値に戻される。

完成後は既知原因と次の操作を短く案内し、任意Detailsへ許可済みの有界診断を表示する。
Copy diagnosticsは安全な表示内容だけを手動コピーする。
不正な既存Color ruleのpattern編集値を未確定入力として保持し、修正・Retry・Revertを可能にする。
Python ruleの意味、正常live edit、最後の成功Result、History、Session、再生成、取消・古い応答の隔離を維持する。

PDFの日本語・中国語font対応（BUG-07）、比較identityのsource/view根本修正、
Feature SearchのPython化、全入力欄への汎用draft機構は対象外。
比較identity不一致は識別・案内の対象とし、実入力で根本原因が再現した場合は別の限定修正として追跡する。
本計画の実装完了とIssue全体の閉鎖は別の判定である。

## 承認された独立した製品判断

Product Decision Owner: satoshikawato。承認日: 2026-09-26。

| Concern | 採択済み動作 | 本文と機械表現 |
| --- | --- | --- |
| web.errors.diagnostic-disclosure / scenario 1 | A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS | [診断情報公開](./decisions/DECISION_01_ERROR_DISCLOSURE.md) |
| web.rules.rejected-pattern-edit-recovery / scenario 1 | A / KEEP_REJECTED_PATTERN_DRAFT | [pattern編集後の回復](./decisions/DECISION_02_REGEX_EDIT_RECOVERY.md) |

各記録は理由・維持条件・退役範囲・残余リスク・owner・dateを含む承認本文の保存である。
コードの都合で製品動作を再選択・拡張しない。新しいリスクを承認済みと扱わない。
このdirectoryは新しいdecision storeではなく、正式authorityは既存
[Option Integrity Product Contract](../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)に関心ごとに記録する。
S01でauthority-only変更を別ブランチへcommit/pushし、devへの統合を確認してからS02以降を開始する。
同じcandidateのauthorityでruntimeを自己承認しない。

OIPC-C07、PD-OI-016、PD-OI-031/034、Webのlive-edit契約は原因・retry・旧状態保持をすでに要求する。
番号とscenario revisionはS00/S01で最新devの原本から確認し、存在しないBD番号を付けない。

## 調査済みの事実と限界

[初期証拠](./BASELINE_EVIDENCE.md)にsynthetic probeと既存テストの再現入口を収録した。
調査SHAは2edc00aebc74e01003da643dfc957b513d5dcfe5。
計画基点af5d942aまでのgbdraw/tests/tools/.githubは未変更で、この証拠を初期状態の説明に再利用できる。

Color/Labelは[PR #538](https://github.com/satoshikawato/gbdraw/pull/538)でPython ownerに統一済み。
(?i)NADH、(?P<enzyme>NADH)は既存Node/Python/Chromium検査で受理される。
Feature SearchとstandaloneはJS RegExp(query, 'i')で同じPython専用構文を拒否する。
新しいrule評価器を作らない。

Issue参照のGUI_AUDIT_DEV_20260926.mdは調査treeにない。
原監査の入力欄・pattern・buildは未特定。S00で取得できない場合は未確認を記録する。
現在のローカルColor/Label成功を、原報告・公開サイトの修正確認と表現しない。

## 実装アーキテクチャと所有境界

~~~text
既存Python/JS producer
  → failure { code, operation, stage, context }
  → 既存Python adapter / Worker / client
  → error-normalizationの有界な表示model
  → operationに対応する既存UI
  → 既存ownerのedit / retry / export action
~~~

| 責任 | Owner（gbdraw/基点） | 方針 |
| --- | --- | --- |
| 比較identityの失敗 | exceptions.py / render/groups/linear/pairwise_match.py | 専用型で識別。ValueError捕捉・整合性検査を維持 |
| Python例外採取 | web_support内の限定adapter / web/js/app/python-helpers.js | string化前にtype/causeを採取。render/helperを収束 |
| transport | web/js/workers/diagram-generation-worker.js / web/js/services/diagram-generation.js | 同じbounded envelopeを両側同時更新 |
| 文言と公開model | web/js/services/error-normalization.js | 唯一の意味owner。有限code・allowlist、再正規化で原因維持 |
| Generate/Align outcome | web/js/app/run-analysis.js / web/js/app/similarity-alignment.js | callerへerrorを返し、global errorLogを読み戻して推測しない |
| UI composition | web/js/app/app-setup.js / web/index.html | 正しい見出し、Details/Copy、実在する回復action |
| rule field draft | web/js/app/feature-editor/rule-actions.jsと既存composition | target/revisionに結び付く一owner |
| regex検索 | web/js/app/feature-search/search-core.js / web/js/services/standalone-interactivity-assets.js | JS意味を維持して方言・不正案内を明示 |

producerは失敗の意味を所有し、利用者文言を持たない。
regexはre.errorまたは明示causeから分類し、message.includesやtraceback末尾解析を新しいclassifierにしない。
cleanup副因は主因を置換せず、unknownは実際に分かるoperation/stageとstable codeを持つ。
不明な位置・原因を捏造しない。native CLI/APIの既存例外捕捉を維持する。

normalizerの上限はsummary 1,000字、detail text 4,000字、section 8個を維持する。
contextは有限の許可fieldと型・長さを検証する。行番号・field識別子・Python character位置・有限reasonを使い、
raw pattern、sequence、file/record名、path、SVG、自由な例外文/stdout/stderr/tracebackを公開しない。
Python character位置とJS UTF-16 indexを同一扱いしない。
移行対象のknown validation修正情報を一覧化して保持し、すべてunknownへ潰さない。
置換した個別prefix・raw表示・同じfailure用特別extractorを同じ変更で除去する。

旧Result保持の文言はtransaction ownerの事実に基づく。
初回Resultなし、正常rollback、rollback自体の失敗を区別し、復旧成功を捏造しない。
canceled/stale/supersededはerrorと独立したoutcomeに保つ。

### Regexと未確定field

Color/Labelは既存preparation → EVALUATE_RULES → web_support/rule_matching.py →
rendererの既存owner → 現在性確認 → atomic commitを維持する。
構文確認と照合は同じ一回のPython評価。pendingはvalidでもnon-matchでもない。
revision/catalog/Result/modeのstale guardとprepared reuseを維持する。

Color/LabelにPython regex・case-insensitive、検索にJavaScript regex・case-insensitiveを示す。
準備失敗をInvalid ruleと呼ばない。未確認whitelist/visibility入口は実ownerを追ってから説明を変える。

draftの対象は既存Color ruleのpattern fieldだけ。
同じdocumentのdrawer close/reopen・一時mode切替では未確定textを保持する。
対象ruleのUndo/Redo置換、row削除、document/session成功置換、resetではownerが解放する。
Session置換が失敗した場合は旧documentのdraftを保つ。
非同期結果はtarget/revisionと現在documentが一致した時だけ採用する。

不正draftとRetry/Revertはartifact Historyを増やさず、成功editだけ既存一entryでcommitする。
新しいrule ID/schemaやSession fieldは不要。
Save/Generateはlast accepted rule、Exportは現在Resultを使い、Not appliedの案内にその事実を示す。
draft表示はcanonicalの第二ownerではない。template/finally/watcherにdraft遷移を散らさない。

## SOLID / KISS / DRY / YAGNI

| 原則 | コード・architecture・workflowへの適用 |
| --- | --- |
| SRP | producer、transport、文言、表示、回復action、draftを分担。二つのProduct記録も独立 |
| OCP | 既存境界へ有限codeを追加。汎用exception/plugin/draft frameworkは不要 |
| LSP | ValueError捕捉、native API、正常live edit、cancel/staleを維持 |
| ISP | UIはbounded model、rule ownerは限定evaluator/actionを受け取る |
| DIP | normalizerはVue/Pyodideを構築せず、必要なactionを既存compositionから注入 |
| KISS | 一Worker、一評価、一normalizer。一実装branchで直列に引継ぐ |
| DRY | 文言・型を収束し、取得・commit/push手順は本計画だけで定義 |
| YAGNI | regex翻訳、追加runtime、永続draft、migration、外部送信、汎用workflow基盤なし |

[Architecture Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)に従い、
能力ごとのbefore/after owner/pathと削除した旧経路を簡潔に記録する。
OEは余剰semantic owner、PEは余剰production path、CBはpersisted compatibility path。
目標はdelta(OE) <= 0、delta(PE) <= 0、delta(CB) = 0。例外条件時だけ規約の完全なpacketを用意する。
[Product policy](../PRODUCT_IMPACT_RATCHET.md)のAND-of-ORで独立要求をすべて保持し、
同じoption IDやテスト成功だけで保存・回復・accessibilityの別要求まで満たしたと扱わない。

## セッション一覧

全セッションを直列で行う。次セッションは先行push済みSHAとresultを取得してから開始する。

| Session | INSTRUCTION PROMPT | 主な所有範囲 | 開始条件 |
| --- | --- | --- | --- |
| S00 | [基準・入口・重複確認](./sessions/SESSION_00_BASELINE.md) | inventory・初期証拠。runtimeなし | 計画branch取得 |
| S01 | [正式契約化](./sessions/SESSION_01_AUTHORITY.md) | 別authority branchのOIPC・結果記録 | S00完了、二receipt確認 |
| S02 | [構造化原因とtransport](./sessions/SESSION_02_ERROR_BOUNDARY.md) | Python producer/adapter、Worker/client、normalizer、tests | S01 authorityがdev統合済み |
| S03 | [案内・Details/Copy・方言](./sessions/SESSION_03_ERROR_UI_AND_DIALECT.md) | callers、composition、index/search、tests | S02 push済み |
| S04 | [Color rule未確定入力](./sessions/SESSION_04_REGEX_DRAFT.md) | field owner、lifecycle、限定UI、tests | S03 push済み |
| S05 | [統合受入と仕様更新](./sessions/SESSION_05_ACCEPTANCE.md) | integration、既存reference docs、原因ownerの必要修正 | S02–S04 push済み |

## 他セッションとの干渉回避

実装branchのwriterは常に一セッション。隔離treeでも並行pushは許可しない。
開始前に先行writerの終了とpush済みresultを確認し、共有tree/index/editable環境を使わない。

既存fix/issue-601-export-output-20260926にもBUG-15/19作業が含まれる。
そのS03はnormalizer、Python glue、Worker/client、Generate/Align、index/searchを本計画と共有する。
Issue #598/602のAlign/compositionにも共有fileがある。
S00は各remote branchとresultを読み、BUG-15/19のowner引継ぎを確認するまで重複runtimeを編集しない。
PDF側のbranchや計画を本セッションから変更しない。

既存別候補にはweb.errors.user-facing-diagnostic-disclosureがある。
S01はその保存条件と今回のCopy/draft条件を比較し、同じ意味に競合するactive authorityを作らない。
既にdevに入った実装は再利用して残る差分だけを実装し、別error registryへ複製しない。
キー変更・製品退役・矛盾解消を実装者が推測して承認したことにしない。
実際に矛盾する条件が残れば影響concernだけを止め、根拠付き差分をProduct Decision Ownerへ提示する。

各セッション専用cloneを標準とする。cloneはrefs/index/build outputも隔離する。
環境、server port、Playwright output、wheel、temporary filesも専用にする。

## 共通取得手順

**すべての実装セッションはorigin/fix/issue-601-bug15-bug19を取得して当該ブランチで作業する。**
main/devや別作業branchから独自に実装を始めない。
SESSION_CODEだけを担当promptの値へ変更する。未使用pathを作り、既存directoryを削除しない。

~~~bash
SESSION_CODE=s00
ISSUE601_SESSION_ROOT=$(mktemp -d "/tmp/gbdraw-issue601-${SESSION_CODE}-XXXXXX")
git clone --no-checkout https://github.com/satoshikawato/gbdraw.git "$ISSUE601_SESSION_ROOT/repo"
cd "$ISSUE601_SESSION_ROOT/repo"
git fetch origin
git switch --track -c fix/issue-601-bug15-bug19 origin/fix/issue-601-bug15-bug19
git pull --ff-only
git status --short --branch
git rev-parse HEAD
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
~~~

upstreamはorigin/fix/issue-601-bug15-bug19、HEADは取得remote branchであることを確認する。
dirtyなら自分の差分か確認し、他セッションの変更をreset/clean/stashして開始しない。
先行SESSION_XX_RESULT.mdの開始条件・検証とremote headを確認する。

必読: [AGENTS.md](../../../AGENTS.md)、[CLAUDE.md](../../../CLAUDE.md)、
[Web CLAUDE](../../../gbdraw/web/CLAUDE.md)、本計画、二承認記録、担当prompt、先行result、
[Product policy](../PRODUCT_IMPACT_RATCHET.md)、[Web policy](../WEB_CHANGE_POLICY.md)。
最新の具体的規約に従い、outcome変更が必要なら黙って計画を書き換えない。

S01だけはauthority/runtime分離のため別authority branchへ移る。
S02はauthorityのdev統合SHAを確認し、実装branchへ通常mergeで取り込む。
公開履歴のrebase/force-pushはしない。各sessionでdevから別実装branchを乱造しない。

## 環境と検証の隔離

必要なセッションだけ専用clone内で環境を準備する。

~~~bash
python -m venv "$ISSUE601_SESSION_ROOT/venv"
. "$ISSUE601_SESSION_ROOT/venv/bin/activate"
python -m pip install -e ".[dev]"
npm ci
export PLAYWRIGHT_BROWSERS_PATH="$ISSUE601_SESSION_ROOT/playwright-browsers"
python -m playwright install chromium
~~~

同等な隔離環境は再利用できるがpath/versionをresultに記録する。
共有環境に別checkoutをeditable installしない。
Node @playwright/test、CLI playwright、Python importを確認し、
Node不可時はPython Playwrightで同じ受入を実行する。
Chromium sandbox起動失敗は必要なescalationで同じ確認を再試行する。

wheelは必要時にpython tools/prepare_browser_wheel.pyで生成するgitignored asset。
テスト用wheelのためにcache-bustを更新しない。
GBDRAW_WEB_TEST_PORTとPlaywright outputを専用の未使用値にする。
別checkoutの既存serverを再利用せず、自分が起動したserverだけを終了する。

## 終了時のcommit・push

**各セッションの実装終了後は検証・差分reviewを完了してコミットし、同名remote branchへ必ずプッシュする。**
通常sessionは一実装commitにまとめ、先行sessionとbase mergeの履歴は書き換えない。
S01はauthority commitと実装branch上の非規範handoff記録を分ける。

1. SESSION_XX_RESULT.mdに開始SHA、outcome、owner/path差分、削除した旧経路、commands/results、
   known failure移行一覧、限界、次開始条件を記録する。
   自分自身の最終commit SHAはcommit後のhandoffで伝え、自己参照SHAのためにamendしない。
2. production/tests/docs/generatedを別々にreviewし、変更fileを明示してgit addする。
   git add .で別作業を混ぜず、wheel/dist/egg-info/temporary/browser outputをstageしない。
3. branch/upstream、git diff --cached --check、staged statを確認する。
4. English commit titleと短いsummaryでcommitする。
5. git fetch originでremote更新を確認し、増えていたら原因を確認して通常merge・必要検証を行う。
   push失敗をforce-pushで解消しない。
6. git push origin HEAD:refs/heads/fix/issue-601-bug15-bug19でpushする。
7. git ls-remote --heads origin fix/issue-601-bug15-bug19のSHAとgit rev-parse HEADの一致を確認する。
   成功済みpushを再試行しない。失敗したmutationの再試行前にもremote実状態を確認する。
8. commit SHA、remote branch、result file、次開始条件をhandoffする。push未完了を完了扱いしない。

S01のauthority branchも同じbranch/upstream/remote SHA確認とcommit/pushを行う。
PR作成、dev merge、main promotion、release/deployは別の明示的許可が必要で、
main/devへの直接pushで代用しない。PR文面にはwrite-clear-pull-request skillを適用する。
remote CI/status監視は五分以上の間隔を守る。

## 必須統合受入

| ID | 合格条件 |
| --- | --- |
| E01 | render/helper両経路でcode/operation/stage/許可contextがWorker/client/caller/UIまで保持 |
| E02 | known causeと修正情報、idempotent normalization、二重prefixなし |
| E03 | raw/private sentinelがsummary/Details/Copy/consoleに出ず、型・長さ・section上限を維持 |
| E04 | cleanupは主因保持。初回/正常rollback/rollback失敗を真実に表示 |
| E05 | 実run-analysisのengine-error/catch→Align Applyで原因保持、旧draft/Result/orientation/History/request保護、retry成功 |
| E06 | canceled/stale/supersededにerror UIや不正commitなし。旧応答が新状態を上書きしない |
| E07 | 正しいoperation見出し、任意Details/Copy、実在action、clipboard不可時の手動選択 |
| R01 | (?i)、(?P<name>)、\\Z、\\bβ、Unicode foldのColor/Label targetがnativeと一致 |
| R02 | 不正[ / JS named groupをempty/unrelated catalogでも拒否。manual/TSV/preset/Session/History/Generate保持 |
| R03 | Feature Search/standaloneのJS意味と方言案内一致。Python runtime追加なし |
| R04 | rejected Color field、Not applied、Retry/Revert、syntax/runtime区別、成功のみHistory一entry |
| R05 | close/reopen・一時mode切替はdraft保持。row/対象History/成功Session置換/resetは解放、失敗Session置換は保持 |
| R06 | draft非永続。Save/Generateのaccepted rule、Exportの現在Resultと案内一致 |
| U01 | Circular/Linear、desktop/390 px、keyboard/focus、failure/readiness/操作到達性をbrowser確認 |
| P01 | 既存25,000-feature preparation/reuse、Worker construction/call数、追加validate-onlyなし |
| SCI01 | 有効target/SVG/科学的意味、native捕捉、Session/request schemaが不変 |
| G01 | 関連Node/Python/browser、必要なread-only SVG comparison、trusted-base Web Gate/architecture契約が合格 |

runAnalysis stubだけでE05を完了とせず、実orchestrationのfailureを接続する。
原監査build、未検証browser、未実行supported Python matrixを明示する。
test timeoutを緩めず、長いcommandは30分以上の余裕で進捗監視する。
未変更領域の有効証拠は再利用し、新変更/失敗/未解消懸念だけを追加検証する。

## 基本検証コマンド

専用cloneのrepository rootで実行し、新設した回帰検査も含める。

~~~bash
node --test tests/web/error-normalization.test.mjs tests/web/rule-matching.test.mjs tests/web/similarity-alignment-actions.test.mjs
python -m pytest tests/test_web_rule_matching.py -q
GBDRAW_WEB_TEST_PORT=8791 node node_modules/@playwright/test/cli.js test -c playwright.functional.config.js python-rule-parity.playwright.spec.js --workers=1
ruff check gbdraw/
pytest tests/test_output_comparison.py::TestOutputComparison -v
~~~

8791は例。実行時は未使用portにする。
actual headのpolicyはtrusted base checkoutのcheckerを使い、
node tools/check-web-change-budget.mjs --base BASE_SHA --head HEAD_SHAで評価する。
full SHAを指定してresultへ記録し、candidate detector/rulesでcandidate runtimeを許可しない。
新blocking受入は既存PR gate inventoryにも入れ、local-only browser成功だけで自動安全としない。
checker/authority/CI guard変更が必要なら別prerequisiteとし、runtimeに混ぜない。

公開仕様は既存REFERENCE/web-app.md、input-formats-and-tsv-schemas.md等のownerを更新し、ページを新設しない。
Galleryの操作画像が実際に変わる時だけ関連skillを適用する。
通常reference SVGはread-only、owner-maintained social previewは変更しない。

## 完了条件

S05で必須受入・該当gateが合格し、全sessionのcommit/push/remote SHA/resultが確認できたら実装完了。
authorityのdev統合、runtime PR/CI/dev staging、Issue全体の閉鎖、main/deployは別状態として報告する。
合格条件を弱めるfallbackやknown→unknown一括変換で完了扱いにしない。
計画・承認記録・基準テスト成功は、未実装runtimeの受入成功ではない。
