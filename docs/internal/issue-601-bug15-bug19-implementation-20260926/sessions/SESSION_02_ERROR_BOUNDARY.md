# INSTRUCTION PROMPT — S02 — 構造化原因・Python adapter・Worker transport

あなたは失敗の識別と伝達、利用者向け正規化境界の担当です。
既存producerからrender/helper両経路で有界な原因を保持し、一つのnormalizerへ収束させてください。

## 必ず当該ブランチを取得する

SESSION_CODEはs02です。
[総合計画書](../MASTER_PLAN.md)の共通取得手順で**origin/fix/issue-601-bug15-bug19を取得し、
専用cloneのfix/issue-601-bug15-bug19で作業してください**。
共有treeや他sessionのbranch/index/環境を使わず、同じremote branchのwriterを一つにします。
AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、
[診断公開の承認](../decisions/DECISION_01_ERROR_DISCLOSURE.md)、
[field回復の承認](../decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)、先行SESSION_XX_RESULT.mdを読んでください。

## 開始条件

S01 resultを読み、必要な静的authorityがorigin/devに実在し、統合commitをancestorとして確認する。
candidate branchだけにある契約では開始しない。
実装branchへgit merge --no-edit origin/devで正式authorityを取り込み、必要な競合をownerの意味に沿って解決する。
receipt/outcomeが異なる場合は影響する変更だけを止める。

## 所有範囲

- gbdraw/exceptions.py、render/groups/linear/pairwise_match.pyの比較整合性failure。
- web_support内の限定error adapter、web/js/app/python-helpers.jsのrender/helper boundary。
- web/js/workers/diagram-generation-worker.js、web/js/services/diagram-generation.js。
- web/js/services/error-normalization.jsと関連Python/Node tests。
- 本計画directoryのSESSION_02_RESULT.md。
- 画面のDetails/Copy配線とfield draftは後続ownerへ引き継ぎ、本sessionで別UIを作らない。

## 作業

1. S00のknown failure一覧を基にfinite code/operation/stage/contextを確定する。
   unknownへの一括置換をせず、既知validationの修正情報の移行先を明示する。
2. 比較identity不一致に識別可能な専用例外を付け、ValueError捕捉の互換を保つ。
   source/view identityを補正・無視しない。native CLI/APIのcauseと動作を維持する。
3. Python regexはre.errorまたは明示causeから理由/位置を採取し、raw stringで分類しない。
   Python character位置を保持する。原文patternやpathをpayloadに入れない。
4. render/helperの文字列化前に同じ限定serializerを通す。
   既存callJsonHelperとrender wrapperの両経路を接続し、validatorを重複実装しない。
5. Worker serializer/client deserializerを同時に更新してcode/operation/stage/contextを保持する。
   cleanup副因はboundedな識別だけを持ち、主因を上書きしない。
   cancel/stale/supersededの区分とlazy Worker lifecycleを変えない。
6. normalizerを文言/公開modelの唯一ownerにし、有限codeから短い原因とaction IDsを作る。
   再正規化してもcode/contextが失われず、prefixが重複しない。
   stage不明、syntax位置不明は事実に沿って示す。
7. summary/detail/section上限、context allowlist・型・長さを検証する。
   raw例外/trace/自由stdout/cleanup/sequence/SVG/file名を表示modelへ流さない。
   input正規化の既存pathとSession/request schemaを変更しない。
8. 移行したfailureのsuperseded serializer/extractor/raw表示経路は同じ変更で除去する。
   新旧分類器や互換Worker protocolを並行所有しない。

## 必須検証

- renderとhelperの実境界→Worker/client往復でcode/stage/contextを比較する。
- ValueError互換、re.error cause、empty catalog syntax、unknown、cleanup、double normalization。
- sentinelを原exception/cause/stdout/contextへ入れ、公開modelに出ないことと上限をassertする。
- known validationの必要情報を保持し、safe unknown fallbackが分類漏れを隠さないこと。
- native Python rule parity、既存41 Nodeの対象、関連Worker/helper tests、architecture/policy。
- 成功request/Session、Worker constructionと既存prepared reuseの回帰なし。

SESSION_02_RESULT.mdへfinite contract、failure移行表、owner/path削除、commands/results、
S03が使う表示model/actionと未確認事項を保存する。
公開される文言は安全化済みでも、UI完了と報告しない。

English commit title: Preserve structured failure causes across the diagram worker
Summary: Converge producer errors, transport, and bounded user-facing normalization.

## セッション終了時のコミット・プッシュ

総合計画書の共通終了手順に従い、SESSION_02_RESULT.mdと対象変更を
**検証後にコミットし、当該同名remote branchへプッシュしてください**。
通常のtargetはfix/issue-601-bug15-bug19です。
branch/upstream、staged scope、remote実状態を確認し、force-pushやmain/dev直接pushは行いません。
remote/local SHA一致、result file、完了/未完了の開始条件を次sessionへ引き継ぎます。
