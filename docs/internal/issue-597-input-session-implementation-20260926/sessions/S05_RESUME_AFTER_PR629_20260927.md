# 次セッション用プロンプト — Issue #597 / S05、PR #629 統合後

gbdraw Issue #597 の S05 を再開し、残る受入条件を解消してください。
初稿・診断・計画だけで終わらず、許可されたローカル実装と関連検証を完了してください。
S06 以降には進まないでください。mapped evidence／authority／runtime の分離、
Product 判断、必要な外部承認を守り、FAIL や未測定を免除して完了扱いにしないでください。

## 1. 作業場所と actual refs

Repository: https://github.com/satoshikawato/gbdraw.git
保全済み S05 checkout: `/tmp/gbdraw-issue597-S05.ujOWyl`
S05 branch: `fix/issue-597-input-session-20260926`
Upstream: `origin/fix/issue-597-input-session-20260926`
前セッション終了時の local／actual remote S05 HEAD:
`ee6b6207845cc861ae428b4890a091509e1736b6`
確認済み actual trusted dev:
`7b321b606fc33c72827d0c3d3ef77288a69424cb`

保存 SHA を無条件に latest とせず、actual remote を確認してください。
不足・進行した refs／objects だけ取得し、latest trusted dev を専用作業環境へ
ローカル統合して、測定・実装の source を一致させてください。
S05 runtime の作業対象は指定 branch です。別の独立 delivery が必要なら、
AGENTS.md に従い latest origin/dev から新規 work branch を作ってください。
共有 checkout、他担当者の branch／server／worktree を操作しないでください。

S05 checkout の dirty state は新規 result／evidence が中心です。
未コミット evidence を未完了 runtime と決めつけず、最初に diff を分類してください。
元の tracked files 3,146個と当時の pinned evidence 183個は保全済みです。
追加された証拠も含めて開始時に保全範囲と hashes を更新し、既存 evidence、
fixture、cache、他セッションの変更を上書きしないでください。

## 2. 既に成立した前提を再び阻害条件にしない

PR #629 は evidence-only delivery として dev に MERGED:
https://github.com/satoshikawato/gbdraw/pull/629
Published head: `e85444b5e2724279bdac5b6fc9becdb5da1a0d85`
Actual merge: `7b321b606fc33c72827d0c3d3ef77288a69424cb`
Base before merge: `968d08211ce2315c274be14d879ef9ceb332bb66`
Approved file-manifest SHA:
`3e211a987c7c8fdd55dec2bb55f1558f8a7a7267a1b1897c399d02a8b5f27a6a`

必須 `Web base policy (trusted base)` と `PR / gate` はこの head で SUCCESS。
Base checker の実ログも Gate PASS / Review CLEAR を確認済みです。
24 paths は承認されたテスト1ファイル＋報告/evidence 23ファイルと一致します。
production／authority／guard／workflow／reference／fixture の変更はありません。
GitHub required approving review count は0。ユーザーはこの候補の公開・mergeを
明示承認しましたが、存在しない GitHub review を作った扱いにはしないでください。
merge 後の Dev staging は成功確認済みではありません。

更新対象は `tests/web/contracts/session-regenerate-intent.playwright.spec.js`。
trusted body SHA-256:
`732892b8d9c780cc4713fc41f4d28032e171860f5a176cb1b0bc49a1bf5693c7`
- Save helper は正規の Updating diagram busy、Result/request/overrides 保持、
  pending flag と availability を確認し、更新終了後に同じ Save を明示的に再試行する。
- 全試行が元の単一 180000 ms download deadline を共有し、deadline を再開始しない。
- divergent case は非同期 capture 後にも availability が clear で、captured SVG が
  current mounted SVG と一致するまで pre-Save snapshot を待つ。Save 後へ baseline を移さない。
- 元の142個の expect 行、case timeout、budget、Worker criterion、fixture、strict 比較を保持。
- exact path::test-name references は不変なので authority-reference amendment は不要。

この先行 mapped-evidence integration は成立済みです。旧レポートにある
「この候補の公開承認／先行 PR が必要」という阻害条件を再適用しないでください。
一方、この delivery は依存 S05 runtime を承認・admit したものではありません。
今後 mapped contract/reference を変える必要が出れば、再び所定の独立順序を守ってください。

Feature fill の旧 line-193 FAIL は PR #626 の既存 runtime 修正で解消済みです。
Runtime commit: `03621c0ecc1498738535b0ffc79d71ee160e6c4f`
元の Feature fill body を変更せず PASS。canonical rules／legend／History の
同一 transaction を使い、別の legend owner を追加していません。
利用者が明示する Legend color/stroke edit は direct-edit journey が独立して検証しています。
dev／provisional S05 各5ケース、S05 divergent 追加2回が最終候補で PASS。
同じ Save helper の実 reflow busy/retry 補助 probe も PASS。
関連 source が変わらない証拠は hashes と成立条件を示して再利用してください。

## 3. 必読資料と evidence locations

まず AGENTS.md、CLAUDE.md、gbdraw/web/CLAUDE.md と次の現行 trusted 資料を読む:
- docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
- docs/internal/PRODUCT_IMPACT_RATCHET.md
- docs/internal/WEB_CHANGE_POLICY.md
- docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md

以下は `docs/internal/issue-597-input-session-implementation-20260926/` 配下:
- MASTER_PLAN.md、SESSION_WORKFLOW.md、sessions/S05_INSTRUCTION_PROMPT.md
- results/S01_RESULT.md、S04_RESULT.md、S05_RESULT.md
- results/S05_TRUSTED_GUARD_INTEGRATION_20260927.md
- results/S05_CURRENT_EDIT_EVIDENCE_20260927.md
- results/S05_RESUME2_20260927.md
- results/S05_SAVE_RETRY_EVIDENCE_20260927.md（PR #629 で dev にある）
- results/S05_SAVE_RETRY_TRUSTED_INTEGRATION_20260927.md（最新の統合 receipt）
- evidence/S05-resume2-20260927-t65_9osr/ の validation.json、fingerprints.json、
  final-review.json、preservation-final.json、reused-evidence.json、raw-checksums.json
- evidence/S05-save-retry-publication-20260927-p4oldby6/ の integration-verification.json、
  final-review.json、poll-002.json、trusted-base-ci.json、pr-gate-ci.json、checksums.json
- evidence/S05-save-retry-evidence-20260927/ の mapped-review.json、
  assertion-budget-review.json、candidate-manifest.json（PR #629 の trusted files）

Raw measurements/comparison views:
`/tmp/issue597-S05-resume2-20260927-t65_9osr/`
Publication/merge raw evidence:
`/tmp/issue597-S05-evidence-publish-20260927-p4oldby6/`
Original fixture/cache/CLI evidence:
`/tmp/issue597-S05-resume-evidence/`

旧 raw `s05-view` は ee6b6207 + dev968d0821 の no-commit merge を保持しています。
clean checkout や最新 source とは扱わず、過去の比較 evidence として保全してください。
旧 evidence-view は公開 commit e85444b5 の clean branch で、upstream は同名 remote。
historical FAIL、候補 v1–v4 の FAIL、最終 settled-capture PASS、CI SUCCESS を区別してください。
古い authority revision／prerequisite を現行 trusted base より優先しないでください。
PR #624 merge `494091aa68ca59ffa27ecaa6c3df19da4fbf5090` の ancestry、
取り込まれた permissions／guards／mapped identities も確認してください。

## 4. 残る S05 acceptance を解消する

元の real fixture:
`/tmp/issue597-S05-resume-evidence/fixture/real-full-pairwise.gbdraw-session.json.gz`
SHA-256: `d3cafef9664bff958aec4932eea6e264ae872cce89762a7cbe617447af975388`
134471323 compressed / 466767669 expanded bytes、12 DNA records、62089 features、
144 full-pairwise protein entries。元の fixture preparation と pinned inputs を再利用する。
元の production measurement recipe は `evidence/S05_MEASURE.py` と
`tools/measure_issue597_s01.py`。各 mode の exact commands は前述の invocations にある。
測定 recipe exit 0 は観測完了であり、acceptance PASS ではありません。

1. Real gzip transport responsiveness
   最新3回の heartbeat max は 507.2 / 1053.5 / 518.9 ms。
   unchanged max≤500 ms に対して FAIL 3/3。historical 512–526 ms FAIL も保持する。
   read/decompress/decode/parse、Worker reply、native receiver、candidate preparation を
   分けて測り、一方式の transport を修正・再測定する。
   fixture 縮小、codec-only への置換、threshold 引上げ、main-thread parse fallback を禁止する。

2. Real full Load/Save responsiveness
   Load wall 123.7–129.9 s、heartbeat max 105.8–110.1 s。
   preview mount だけで約98–102 s。Save wall 20.6–22.7 s、heartbeat max 1.68–2.06 s。
   両方3/3 FAIL。projection／validation／restore／SVG admission/mount／History／
   compression／handoff を含む全経路を、元の代表入力・測定定義・budgets で検証する。
   S05 の適用可能な full-path 受入条件を「S06で対応」として免除しない。
   不要な full clone/signature/base64 再変換、第二 owner、size別二重経路を増やさない。

3. Real saved-preview Load の Python Worker count 0
   最新3回とも1。`validateConfigOverrides` が config preflight で既存 diagram Worker を
   構築する。stored unmanagedConfigOverrides がない経路を追跡する。
   必要な validation を省略して0にしない。小 fixture PASS や後続 color helper に
   許された Worker 1 と混同しない。delayed discovery と Inspect/Generate continuation を保持する。

4. Real Generate の Missing GenBank file
   safe summary は `Sequence #1: Missing GenBank file.`。
   canonical request は12 GenBank resource IDsを持つが、restored UI は入力 file のない
   placeholder 1件。Worker は init/validateConfigOverrides/evaluateRules までで render run は0。
   records-table CLI writer が cliInputs/cliTables に依存 file を保存する一方、linearSeqs は空で、
   schema-2 explicit empty bindings が復元時に優先される経路を追跡済み。
   これは diagnosis であり修正や新しい Product outcome の承認ではない。
   existing typed request/resource/inventory と source identity を調べ、正規境界で解消する。
   explicit empty draft を committed artifact から自動補充する fallback、別 renderer、argv builder、
   合成 GenBank、main-thread Python を追加しない。fixture 再生成で旧入力の FAIL を隠さない。
   成功後の catalog、編集、History、Session、export まで検証する。

5. Strict CLI SVG agreement と content equivalence
   最新3 saved filesでも CLI root baseProfile="full" が Web にないため strict FAIL。
   共通 sanitizer/admission/serializer/export を追跡し、scientific output 契約を満たす。
   既存 asymmetric verifier は expected root attribute を除くので strict PASS の根拠にしない。
   normalize／ignore／reference overwrite で差を隠さない。意図的 output変更は先に所定 review。
   real pipeline は request/resources/cache/derived cache/manifest の exact equality が true でも、
   webFiles/editorState は false。Vibrio losatCache も false。必要な semantic fields を比較し、
   狭い equality だけで全 content equivalence PASS としない。

6. 未測定 native clone bytes
   別の単回 CDP trace では main deserialize 393.596 ms、HandlePostMessage 394.173 ms、
   enclosing receiver task 394.680 ms、Worker serialize 388.160 ms を実測済み。
   trace overhead が違うので元の3回 FAIL を置換せず、元の null を勝手に埋めない。
   native wire/copy bytes は trace fields に numeric bytes がなく UNAVAILABLE。
   transfer-list buffer bytes 0 や heap delta／JSON size を native bytes の実測と扱わない。
   profiler等で元の測定定義に沿う値を得るか、測定限界と受入上の未解決を明示する。
   必須 null／UNAVAILABLE が残る状態で全測定完了や S05 COMPLETE としない。

各境界で再現・診断から必要な実装と関連チェックまで完了してください。
真の Product／authority／外部承認待ちが生じても、独立してできる調査・証拠保存を完了する。
PD-OI-044／045 等で既に選ばれた outcome を再確認させず、選択範囲を拡張しない。
materially different な Product-valid outcomes が残る部分だけ停止して Decision Pack を用意する。

## 5. 検証、保存、承認範囲

- candidate authority／guard で同じ candidate runtime を自己証明しない。
  mapped変更が新たに必要なら evidence-only → 必要な authority-reference-only → runtime の順序。
- 実装は既存 owners/paths に集約し、不要な compatibility branch／silent fallback を増やさない。
- assertions、case timeout、busy/retry semantics、p95/max budget、Worker-free criterion を弱めない。
- source contents／private genome sequence／comparison rows をログに出さない。
  tests/reference_outputs と owner-maintained social preview は変更しない。
- AGENTS.md の Node/Python Playwright 両経路を確認し、owned port/artifact directory を使う。
  必要な wheel は generated/ignored asset。既存 test-owned timeout は保持し、長い tests は
  少なくとも30分許して incremental monitoring。remote CI/status は5分以上の間隔で確認する。
- unchanged source/input/environment/acceptance の証拠は hashes と条件を示して再利用する。
  material change、failure、未解決部分だけ追加検証し、guard PASS を browser/performance PASS と混同しない。
- 新しい result/validation/fingerprints/logs/checksums は新規 paths に保存する。
  actual refs/branch/upstream/HEAD/dirty state、source/checker/guard/authority/mapped/fixture fingerprints、
  exact commands/env/exit、raw artifacts、instrumentation limits、再利用条件を残す。
  production／tests／docs-evidence／generated diffs を別々にレビューする。
- このプロンプトは必要なローカル実装・検証・独立 delivery の準備を許可する。
  PR #629 の commit/push/PR/merge承認はその候補だけで、既に実施済み。
  別の runtime／evidence／authority candidate の commit/push/PR/merge/deployment に引き継がない。
  現在の task でその具体的 target/scope の明示承認がなければ、reviewable な branch/diff/
  title/body/checks を完成させてから、その残る外部境界だけ承認を求める。
  古い S05 instruction prompt の「必ずcommit/push」を新しい外部承認と解釈しない。

全適用条件、trusted prerequisites、browser/performance/Worker/Generate/strict SVG、
必須測定、required review と最終 source/evidence 一致を満たした場合だけ S05 COMPLETE。
満たせなければ INCOMPLETE として具体的残件と次の一手を報告する。S06+ は開始しない。
最終報告に changes、actual SHAs、checks、historical FAILとの区別、evidence paths、
残る制約、English proposed commit title と short summary を含めてください。
