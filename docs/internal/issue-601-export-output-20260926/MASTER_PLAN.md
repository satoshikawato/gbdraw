# Issue 601 — PDF文字対応・エラー案内の総合実装計画

Author / Product Decision Owner: `satoshikawato`。
作成日: `2026-09-26`。
Repository: `https://github.com/satoshikawato/gbdraw.git`。
実装ブランチ: **`fix/issue-601-export-output-20260926`**。
実装ブランチの基点: fetch済みの最新`origin/dev`、`b11fe6091f179f6e251fe415ecfde269343ec764`。修正前の動作調査はCURRENT_BEHAVIOR_EVIDENCE.mdに記載したsource SHAを基準とする。

## 目的と完成後の動作

[Issue #601](https://github.com/satoshikawato/gbdraw/issues/601)のBUG-07、BUG-15、BUG-19を扱う。
日本語と中国語（簡体字・繁体字）を含む図のPDF出力に対応する。既存fontで表せるアルファベット・数字・記号とstyleを維持し、不足部分だけ対応fontへ自動切替する。PDF内に文字を保持し、検索・選択できる出力を提供する。
失敗した操作には既知の原因と次の操作を示し、許可済みの診断情報を任意Detailsで確認できるようにする。Color/LabelのPython regex評価とFeature SearchのJavaScript regex評価は維持し、それぞれの構文を入力欄で案内する。

この計画は実装担当者が会話や過去の提案を読まずに使う実行仕様である。[修正前の証拠](./CURRENT_BEHAVIOR_EVIDENCE.md)と実装後の受入を区別する。文書の保存・承認はruntimeの実装完了を意味しない。

## 承認済みの独立した製品判断

| 責任 | 承認されたoutcome | 正本 |
| --- | --- | --- |
| PDF文字対応 | SELECTIVE_JA_ZH_VECTOR_TEXT_FALLBACK / scenario 4 | [Decision Pack 01](./DECISION_PACK_01_PDF_JA_ZH.md) |
| 診断情報の公開 | ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS / scenario 1 | [Decision Pack 02](./DECISION_PACK_02_ERROR_DISCLOSURE.md) |

両Decision Packsは`satoshikawato`が`2026-09-26`に承認した。[機械表現](./APPROVED_DECISIONS.md)は本文の忠実なserializationで、独立したCI storeではない。製品outcomeを再選択せず、各Packのpreservation、retirement、riskをそのまま実現する。

[Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md)のLane Bを適用する。実装前にS00のPDF実現性証拠を取得し、S01で既存static Product Contractに二つの独立recordを記録する。PDF証拠が保留された場合は、error recordの正式契約化とS03を独立して進め、PDFの不足checkpointだけを止める。S01の各authority-only差分は`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`一ファイルに限定し、他の文書やruntimeを混ぜない。authorityのbase統合後に依存runtimeを実装する。新しいdetector、JSON registry、approval botをこの問題のために作らない。

## 現在の事実と修正方針

| 項目 | 基点の実装 | 修正 |
| --- | --- | --- |
| PDF | `pdf-fonts.js`がLiberationのfamily/styleを選び、missing glyphで拒否 | 検査を維持して不足部分に日中fontを選ぶ |
| font取得 | `READ_PDF_FONT`→`app/python-helpers.js::read_pdf_font`。Liberation 12filenameのallowlist | 検証した有限manifestから取得を許可し、選択ルールを複製しない |
| Export | `captureSvgExport`→`downloadPDF`→font準備→`doc.svg`→save、finally cleanup | 同期snapshotとgeometryを維持した一つのPDF経路 |
| Error | 既存normalizerあり。raw ValueError公開、stderrのcause消失、Align errorなしreturn | producer/transport/callerでcauseを保持し、既存normalizerでsummary/Details/actionsへ変換 |
| Regex | Color/LabelはPythonに統一済み。Feature SearchはJS | 評価器を増やさずsyntax表示と準備失敗の案内を改善 |

BUG-19は報告者の画面・pattern・buildを未特定。Color/Labelの再実装や全regex画面の「修正済み」判定は根拠にならない。現在の再現と保護する動作を区別して報告する。

## 実装アーキテクチャ

### PDFの一つの経路

```text
click時snapshotとfilename
  → text/styleを読む
  → 一次fontのglyph coverage
  → 不足する文字のまとまりに対応fontを選ぶ
  → 同font/styleの連続runへまとめる
  → 同じfontで位置・幅を計測
  → 原文Unicodeに対応したPDF textを保存
  → staging/font一時資源のcleanup
```

一次fontで表せる部分のfamily/styleを一律に変えない。Unicode code unitで分割せず、surrogate pairや結合単位を保つ。text/tspanの入れ子、anchor、baseline、dx/dy、spacing、rotation、textPathのpositionを維持する。fallback後にもcoverageを検証し、font/styleの不正代用を成功扱いにしない。

font候補は日本語・簡体中国語・繁体中国語に対応するNoto系static PDF-compatible font等をS00で比較し、版と実容量を確定する。現行JP Japanese subsetには簡体字の不足がある。明示text languageを地域字形へ利用し、情報がない場合はversioned manifestの固定順とcoverageを使う。UIの`lang=ja`、codepointだけの言語推測、OS font、font file列挙順で結果を変えない。

原本、checksum、license、style、地域・coverage、必要representationを一つの有限manifestで管理する。package準備でPDF用static fontと必要なDOM計測用representationを生成し、source配信とwheelへ同じ準備経路で含める。browserでのfont変換、runtime CDN、新しいWorker/JS build stepは設けない。必要fontだけsame-originから遅延取得し、cacheは固定asset数で有界、失敗Promiseは破棄してretry可能とする。PDF埋込みは使用glyph subsetを確認し、必要table/派生glyphを失わない。

既存libraryで扱える文字範囲・style・補助平面をS00で証拠化する。font registration成功やPDFファイルの存在だけを受入としない。追加dependency/writerが必要と判明した場合は、具体的な失敗boundaryと代案を示してarchitecture規定に従う。承認済みの製品outcomeを縮めてGateを通さない。

### Errorの責任分離

```text
producer: code / stage / safe params
  → Python adapterがPyodide string化前に識別情報を保存
  → 既存Worker serializer / client deserializer
  → callerへ明示error result
  → 既存error-normalization: summary / safe details / recovery actions
  → UIが実操作と利用可能actionを表示
```

producerはtype/causeを識別し、利用者文言は持たない。Python共通adapterが必要なら`web_support/`に限定し、既存render/helperの双方で使う。比較identity不一致は既存ValueError捕捉と互換の識別可能な例外にする。regexは`re.error`と位置を保持し、自由形式messageの文字列検索による新しいclassifierを作らない。

Pyodideが属性を失う前に構造化する。既存`serializeError`/`deserializeWorkerError`を拡張し、一時protocolを両側同時に更新する。`run-analysis.js`の全error returnを監査し、callerへcauseを明示的に返す。errorLogやglobal errorを読み戻してcauseを推測しない。cancel/stale/supersededの既存区分、Align draft、Result、canonical request、record方向とHistoryを維持する。

文言とknown codeの対応は`services/error-normalization.js`が一度だけ所有する。unknownは原因不明とactual stageを示す。再正規化でcodeを失わずprefixを重複させない。paramsはallowlist・型・長さで制限し、入力・raw exception/traceback/stdout/stderrを画面と通常consoleへ出さない。既知validationの現在の修正情報を移行一覧で保持し、すべてunknownへ落とさない。cleanupの二次失敗は主因を上書きしない。

### Regexの既存経路

Color/Labelは`createRulePreparation`→`EVALUATE_RULES`→`web_support/rule_matching.py`→既存Python ownerを維持する。入力近傍にPython regex、case-insensitive、`(?i)NADH`/`(?P<enzyme>NADH)`を案内する。Feature SearchとStandalone Interactive SVG検索にはJavaScript regex、case-insensitiveを案内する。準備失敗を「Invalid rule」と断定しない。atomic commit、History、stale guardを保持する。

## 所有者とセッションの順序

実装ブランチへの作業は直列に実施する。次のsessionは先行sessionのpush済みcommitとhandoffを取得してから開始する。異なるcloneでも同一branchへの同時pushを許可しない。

| Session | 作業 | 主な所有範囲 | 開始条件 |
| --- | --- | --- | --- |
| [S00](./SESSION_00_FONT_EVIDENCE.md) | font選定とPDF実現性証拠 | このdirectoryの`evidence/`と`SESSION_00_RESULT.md`。production変更なし | 計画branch取得 |
| [S01](./SESSION_01_AUTHORITY.md) | 承認outcomeの正式契約化 | 専用authority branch上のOIPC一ファイル | PDFはS00証拠合格。errorは独立した承認本文 |
| [S02](./SESSION_02_PDF.md) | 日中font fallback、配布、PDF受入 | `pdf-fonts.js`、export、font preparation/loader、関係tests | S01契約がorigin/devに統合済み |
| [S03](./SESSION_03_ERRORS_AND_REGEX.md) | error cause/Details、syntax案内 | Python producer/glue、Worker/client、normalizer、callers、index/search、関係tests | 通常はS02 push済み。独立進行時もS01 error契約統合済み |
| [S04](./SESSION_04_ACCEPTANCE.md) | 全経路の受入・必要修正・handoff | integration testsと`SESSION_04_RESULT.md`、原因ownerの必要修正 | S02/S03 push済み、両契約統合済み |

`app-setup.js`、`python-helpers.js`、`index.html`等の共有fileは表の順序で所有を引き渡す。同一fileを複数sessionが同時に編集しない。別Issueの実装を取り込むための再設計やコード削除をせず、変更された現行ownerへ本修正を合わせる。

## セッション共通の取得・終了手順

### 専用branchを取得して使う

各sessionは作業前にこのbranchをremoteから取得し、専用cloneで使用する。通常checkout、他sessionのworktree、main/devを切り替えて作業してはいけない。専用cloneはGit index、refs、build outputも分離する。session内でだけ使うport、browser output、temporary directoryを選ぶ。

次の`SESSION_ID`を該当promptの値へ設定し、未使用のpathで実行する。既存directoryを削除・上書きしない。

```bash
SESSION_ID=s00
SESSION_ROOT=$(mktemp -d "/tmp/gbdraw-issue601-${SESSION_ID}-XXXXXX")
git clone --no-checkout https://github.com/satoshikawato/gbdraw.git "$SESSION_ROOT/repo"
cd "$SESSION_ROOT/repo"
git fetch origin
git switch --track -c fix/issue-601-export-output-20260926 origin/fix/issue-601-export-output-20260926
git pull --ff-only
git status --short --branch
git rev-parse HEAD
```

検証用のPython environmentもsession専用にし、別checkoutへのeditable installや共有package環境の変更を避ける。必要なsessionでは次を使う。既存の同等な隔離環境がある場合はそれを再利用し、版とpathをresultへ記録する。

```bash
python -m venv "$SESSION_ROOT/venv"
. "$SESSION_ROOT/venv/bin/activate"
python -m pip install -e ".[dev]"
npm ci
export PLAYWRIGHT_BROWSERS_PATH="$SESSION_ROOT/playwright-browsers"
python -m playwright install chromium
```

node_modules、browser output、local server port、wheel、temporary font/PDF outputも専用clone内かsession専用pathへ置く。環境準備をWeb UIの新しいbuild stepとして製品へ持ち込まない。

取得後、[AGENTS.md](../../../AGENTS.md)、[CLAUDE.md](../../../CLAUDE.md)、[Web CLAUDE](../../../gbdraw/web/CLAUDE.md)、本計画、両Decision Packs、該当session prompt、先行`SESSION_*_RESULT.md`を読む。規定が変わった場合は最新の具体的規定に従い、製品outcomeは維持する。

S02開始時は`git fetch origin`でS01契約が`origin/dev`に存在することを確認し、専用実装cloneで`git merge --no-edit origin/dev`して統合済みauthorityを取り込む。公開済み履歴のrebase/force-pushをしない。競合は当該ownerの差分として解消し、不明な他sessionの変更を取り消さない。

### セッション終了時のcommitとpush

各sessionの実装終了後、検証と差分reviewを完了し、**必ずcommitし、同名remote branchへpushする**。機能実装はsessionごとに一つのcommitへまとめる。先行変更と必要なbase merge commitをsquashしたり、履歴を書き換えたりしない。

1. `SESSION_XX_RESULT.md`に対象branch、開始SHA、実装したoutcome、owner/path差分、実行commandsと結果、未確認事項、次sessionの開始条件を記録する。合格条件を緩めた記録で完了としない。
2. production、tests、docs、generated差分を別々にreviewする。対象fileを明示して`git add`し、`git add .`で別作業を混ぜない。generated browser wheel、`dist/`、egg-info、temporary/font probe outputをcommitしない。
3. `git diff --cached --check`と`git diff --cached --stat`を確認する。branchとupstreamが実装branchの同名remoteであることを確認する。
4. English commit titleと短いsummaryを付けてcommitする。
5. `git fetch origin`してremote更新を確認する。異なる更新があればpushを繰り返さず原因を確認し、通常mergeして必要検証を行う。`git push origin HEAD:refs/heads/fix/issue-601-export-output-20260926`でpushする。
6. remote head SHAとlocal SHAの一致を確認し、commit SHA、branch、結果fileをhandoffする。push失敗はremote状態を確認して未完了として報告する。

S01だけはProduct Contract規定のため専用authority branchを使用し、同様にcommit/pushする。authorityのレビュー・devへのmergeは別の公開操作で、main/devへの直接pushで代用しない。PR作成・mergeはその操作の明示的な許可に従う。PR文面を作る場合は`.agents/skills/write-clear-pull-request/SKILL.md`を適用する。

## SOLID / KISS / DRY / YAGNI

| 原則 | コード・architecture・workflowでの適用 |
| --- | --- |
| SRP | font選択、bytes取得、export orchestration、error識別、transport、文言、UIを分ける。製品判断も二つの独立record |
| OCP | 有限font manifestと既存error境界へ追加し、万能plugin/exception frameworkを作らない |
| LSP | 従来成功するPDF、PNG DPI、native ValueError捕捉、cancel/staleの結果区分を維持 |
| ISP | font ownerへ限定loaderとmanifestだけを渡し、state/Pyodide本体を渡さない |
| DIP | font準備はdoc/clone/loaderに依存し、VueやWorkerを構築しない。文言はsourceではなくnormalizerに集約 |
| KISS | 一つのPDF経路、固定fallback順、既存error経路。実装は直列の一branchで引き継ぐ |
| DRY | font原本と準備owner、正常化文言、共通session workflowを各一箇所で管理 |
| YAGNI | 全Unicode保証、font upload、外部export service、新しい永続設定・compat reader・registryを作らない |

[Architecture Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)に従い、普通の変更は簡潔なowner/path根拠と非増加のOE/PE/CBを確認する。OEは余剰の意味的所有者、PEは余剰の実行経路、CBは永続formatの互換経路を表す。旧選択・個別raw表示・同じ失敗の特別extractorなど、置き換えた経路を同じ変更で削除する。Session/request schemaを新設しない。新しいproduction dependency、権限、owner/pathの増加や例外はこのProduct承認で許可されない。Gate失敗は境界を診断して解消し、candidate authorityで同じruntimeを許可しない。

## 必須受入

| 領域 | 合格条件 |
| --- | --- |
| 日中PDF | 日本語・簡体字・繁体字・混在で原文text抽出一致、font埋込み、画像代用なし。既存Latin/Greek/mathのfamily/style維持 |
| Layout | normal/boldと既存italic/bolditalic、tspan/textPath、anchor/baseline/spacing/rotation、legendを確認。欠字/style非対応は明示failure |
| Unicode | 補助平面、結合・異体字は成功または原因を持つ明示非対応。code unit分断や文字化けを保存しない |
| Snapshot | Aの出力click後にBへ変更してもAの内容・名前。failureからSVG/PNGへ回復する場合も同snapshot。Result/History不変 |
| 資産 | license/checksum、source/wheel同一asset、必要font lazy、失敗後retry、cache/cleanup。初回/反復のpayload/parse/encode時間・メモリ |
| Error | source→Python adapter→Worker/client→caller→UIのcode/stage/safe params保持。原因・action・任意Details、summary/detail上限、raw/private sentinel非露出 |
| 整合性 | identity conflictを自動補正せず、Align draft、Result、History、record方向を保持。cancel/staleにerror表示なし。cleanupが主因を上書きしない |
| Regex | Python-only Color/Labelのatomic commit/TSV/History/staleとnative一致。JavaScript検索案内と実際の拒否が一致 |
| 操作 | desktop/390 px、keyboard/focus、エラー展開/再試行。利用可能actionが実際に機能する |
| Scientific output | 既存の図geometry、feature関係、reference SVGを維持。読みやすい実図PDFを独立rendererで確認 |

S00の小さい内部fixtureは実現性検証用、S02/S04は既存の現実的なGallery recipe/sessionで実図を確認する。public figureを作る場合はそのrecipeで再生成し、独立rendererの最終artifactを目視確認する。owner-maintained social previewと通常のreference outputsは変更しない。

## 検証commandとCIの実行境界

必要時に専用clone内で`python tools/prepare_browser_wheel.py`を実行する。deploy用のcache-bust更新は実際にdeployable bundleを準備する時だけ行う。

```bash
node --test tests/web/error-normalization.test.mjs tests/web/rule-matching.test.mjs tests/web/similarity-alignment-actions.test.mjs
pytest tests/test_web_rule_matching.py tests/test_similarity_alignment_web_adapter.py -q
npx playwright test tests/web/export-lazy-loading.playwright.spec.js tests/web/gui-audit-regressions.playwright.spec.js tests/web/python-rule-parity.playwright.spec.js tests/web/similarity-alignment-ui.playwright.spec.js --workers=1 --retries=0
ruff check gbdraw/
pytest tests/test_output_comparison.py::TestOutputComparison -v
```

追加したtestsも各担当sessionで実行する。Node Playwrightの解決を確認し、利用できなければPython Playwrightで同じ焦点のbrowser検査を実行する。Chromium sandbox errorは必要な実行権限で同じ検査を再試行する。長いtest runは30分以上の余裕を持って監視し、test-owned timeoutを短縮しない。

基点のPR gateはtop-levelのNode `*.test.mjs`と非slow Python `browser` testsを実行する。Playwright functional fullとslow offline/package testsにはdev stagingでの実行境界がある。新しいblocking受入は既存PR inventoryのNodeまたは非slow Python browser testへ入れ、local-only specだけをhard automatic safetyの根拠にしない。inventory/CI guard自体の変更が必要ならruntimeと別のchecker-only prerequisiteとして扱う。独立rendererによる目視は自動受入を補強する。

各実装sessionでbase checkerに対するWeb policyとarchitecture contractsを確認する。`node tools/check-web-change-budget.mjs --base <BASE_SHA> --head <HEAD_SHA>`をtrusted base checkoutから実行し、candidate guardを使用しない。通常PR CIと該当dev staging結果を区別して記録し、local focused passを全CI/deploy成功と表現しない。remote状態の監視は5分以上の間隔で行う。

## 完了判定とhandoff

S04で必須受入と関連gateが合格し、全sessionのcommit/pushとremote SHAが確認できたら実装完了とする。資料には変更内容、再現command、artifactの検証、未確認事項を記載する。PR/deploy/mergeの状態は実装完了と区別する。未知原因の詳しい調査が必要な場合でも、raw非公開とcause/retry/状態保持の保証を緩めない。
