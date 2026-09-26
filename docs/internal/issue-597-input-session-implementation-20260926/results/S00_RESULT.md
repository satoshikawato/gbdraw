# S00 — Approved authority candidates and BUG-01 scope removal

この文書はS00実行時点の履歴である。その後のProduct authority統合・dev取り込みは [AUTHORITY_INTEGRATION_RESULT.md](./AUTHORITY_INTEGRATION_RESULT.md) を参照。

日付: 2026-09-26。**候補準備完了、devへのauthority統合は未完了**。
BUG-02 / BUG-20 の候補準備と、ユーザー追加指示によるBUG-01計画削除を実施した。
S01、runtime実装、authority用別PR作成・push・mergeは未実施。

## Checkout / base / environment

| 項目 | 実測値 |
| --- | --- |
| 独立checkout | `/tmp/gbdraw-issue597-S00.y0ZsEa` |
| Branch | `fix/issue-597-input-session-20260926` |
| Upstream | `origin/fix/issue-597-input-session-20260926` |
| 開始HEAD / 取得remote HEAD | `51a786086dc2777e5fa375e5f94d8b7ac7deeedc` |
| 計画作成base / dev由来のancestor | `d457b7189b137185a8dec800819a312c30b969fa` |
| 検証用最新dev | `302dfa1136c95ab50ddef606ad835ec4087b95b2` |
| origin/main | `4556e04e929a4a85ad28d1833ce7304bd764881c` |
| 独立dev worktree | `/tmp/gbdraw-issue597-S00-dev.5c1Lci`、上記devのdetached HEAD |
| ログ保存先 | `/tmp/gbdraw-issue597-S00-evidence.adg7t0` |
| Node / Python | `v26.8.2` / `3.13.3` |

指定のsingle-branch clone → fetch → pull --ff-only → dev/main/tags明示fetchを実行し、
clean branch/upstreamと計画baseのancestryを確認した。
共有checkout、他sessionのbranch/worktree、依存環境、serverは操作していない。
通常sandboxは起動時のbubblewrap mountエラーで失敗したため、同じ専用/tmp内の操作をsandbox escalationで実行した。

## 担当差分と縮小した計画

- BUG-01の専用承認記録 `decisions/01_CIRCULAR_SOURCES.md` と専用 `sessions/S02_INSTRUCTION_PROMPT.md` を削除。
- MASTER_PLAN、索引、関連promptsから source collection、新bindings writer、C-01〜C-05、S02依存、bindings変更を理由とした再生成・例外packetの必須化を除去。
- 残るsession IDは維持し、順序を **S00 → S01 → S03 → S04 → S05 → S06 → S07 → S08** とした。
- Circularは既存scalar入力、一ファイル内の複数records、既存source identity、保存形式・released readersを維持する。BUG-02を新parser/複数source実装に拡張しない。
- `PLAN_PREPARATION.md` は計画作成時の履歴として保持し、現在のscopeを本結果へ参照させた。
- 残る二件のreceiptはbyte単位で開始HEADと同一。再承認、rationale/retirement/riskの補充はしていない。

## Authority intake / serialization

詳細な適用・検証・外部統合手順は [authority-candidates/README.md](../authority-candidates/README.md)。
候補の唯一の適用先は次の通りで、このimplementation branchのactive filesには適用していない。

| 責任 | dev baseでの調査 | 候補 / 適用先 | 状態 |
| --- | --- | --- | --- |
| Product | 新二concernは未mapped。active BDは空。既存OIPC revision 19 / 最大PD-OI-043。 | [01-product-contract.patch](../authority-candidates/01-product-contract.patch) → `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` のみ、revision 20 / PD-OI-044・045。 | 準備・検証済み。dev merge SHA未取得。 |
| Privileged | `new Worker` は現行Diagram Worker detectorに一致。予定clientはbase allowlist外。 | [02-import-worker-permission.patch](../authority-candidates/02-import-worker-permission.patch) → `tools/web-change-policy.json` のみ。 | 一件のconstructor owner候補を検証済み。S01の最終path確認・dev merge待ち。 |
| Checker | exact-path authority enforcementとowner/importer検出は既にdevで稼働。checker/detectorのGit blobは開始HEADとdevで同一。 | 変更候補なし。 | checker-only lane不要。 |
| Mapped evidence | canonical request / saved regenerationの既存concern・contractは維持。新二concernはstatic authority経路。 | 現時点でmap/BD/architecture rules/evidence-ref変更なし。 | 実装で必要ならevidence-only → authority-ref-only → runtime。 |

二件の全9項目を human receipt → 既存JSON → OIPC候補JSONで比較し、完全一致した。
既存43 records、cross-surface/lifecycle clauses、acceptance catalog・risk末尾は不変。
`PD-OI-010` のcrop条件、`PD-OI-019` のcanvas default、`PD-OI-016` / OIPC-C07のfailure isolationを保持する。
候補のreviewed outcome digestsは承認記録の値をそのまま保存し、Markdown自体のdigestと混同しない。

| Receipt | 割当候補 | Choice | 比較 |
| --- | --- | --- | --- |
| `diagram-generation.circular-transform-discoverability` / 1 | `PD-OI-044` | `A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS` | 9/9一致、owner/date/provenance/metadata一致 |
| `web.session-operation-consistency` / 1 | `PD-OI-045` | `A / EXCLUSIVE_SEMANTIC_SESSION_OPERATION` | 9/9一致、owner/date/provenance/metadata一致 |

### Fingerprints

| 対象 | SHA-256 |
| --- | --- |
| Product patch | `8e2638b63daa9f173b9f86be0be661e5597242907e91c7f7f7fb83bb0ac2b19c` |
| Permission patch | `43de78abad3503ecdd9f6f57f1d43b37490aa1538b076edcc2bdd66b40b8c51c` |
| Receipt 02 Markdown | `9f3589dec70a8f86b4718b68cb4dd4b3ef9ad83b6ab93935a437879d64eb539a` |
| Receipt 03 Markdown | `8b9c8d7f185b2c1e65113d1a97062bf96661e2265b8e20e01c39e8583057e859` |
| dev OIPC | `ce6681683e4cf734bb44df0972162764b4dee7bcc1d02ba9c44fa54188820c00` |
| dev privileged policy | `fee8541705ea8bce25939c8aef2f727f64b2cfb1ff2bea12114471452df9b47f` |
| dev checker | `c7888ffb9c6d136284b567712f6e11a474e6b08f54d4daaee7d78e1a402d7f1c` |
| dev detector | `e05a4d4ffa394db8af69f27d82e4673e0e0c2533ba6b0cd7b43dd16e4e227274` |

## Verification

各patchはcleanな同一dev worktreeへ **単独** 適用し、検証後に `git apply --reverse` で差分を除去した。
二つをまとめて適用した状態をauthority-only候補のpass evidenceとして使っていない。
tracked Product patchの空context行はtrailing whitespaceを除去した。標準 `git apply --check` / apply / reverseが成功し、適用後contractは変更前patchとbyte単位で同一だったため、そのauthority/checker evidenceを再利用した。
[verify_candidate.py](../authority-candidates/verify_candidate.py) と候補READMEのcommandsで再現できる。

| Command / check | 結果 / 根拠 |
| --- | --- |
| 各patchの `git apply --check`、単独apply、`git diff --check` | exit 0。許可された一targetのみ変更。 |
| `python .../verify_candidate.py product <dev-worktree>` | exit 0。全18 fields、既存records/clauses、schema/owner、未mapped、provenance・metadataを検証。 |
| `python .../verify_candidate.py privileged <dev-worktree>` | exit 0。owner一件だけ追加、importer permission不変。現行detectorが一constructor・zero privileged import edgesを検出。baseは不許可、candidateは許可、別ownerは不許可。 |
| Product-onlyの `node tools/check-web-change-budget.mjs --base origin/dev` | exit 0、Gate PASS / Review REQUIRED。`product-policy.log`。 |
| Permission-onlyの同command | exit 0、Gate PASS / Review REQUIRED。operator 48 → 49、importer 54 → 54。`privileged-policy.log`。 |
| `node --test tests/web/architecture-ratchet-fixtures.test.mjs tests/web/product-impact-ratchet-fixtures.test.mjs` | exit 0、53/53 pass、skip 0。`authority-fixtures.log`。 |
| implementation checkoutで `node tools/check-web-change-budget.mjs --base 51a786086dc2777e5fa375e5f94d8b7ac7deeedc` | exit 0、Gate PASS / Review CLEAR。S00の変更scope。`s00-scope-policy.log`。 |
| implementation checkoutで `node tools/check-web-change-budget.mjs --base origin/dev` | **exit 1、Gate FAIL / Review REQUIRED**。開始時からdev未取り込みのauthority差分を検出。`implementation-latest-dev-policy.log`。 |
| `node tools/check-web-change-budget.mjs --base origin/dev --head 51a786086dc2777e5fa375e5f94d8b7ac7deeedc` | **同じexit 1**。S00前のimmutable開始HEADでもauthority separation failure。`starting-head-latest-dev-policy.log`。 |
| 最終計画の相対リンク・fences、bash blocksの `bash -n`、Python AST、BUG-01参照除去 | PASS。Markdown 19、相対リンク66、bash blocks 19、Python files 2。active手順と履歴・削除記録を区別して検査。 |

latest-dev差分FAILはOIPC revision 16のimplementation branchと、revision 19のdevを直接比較すると、
既にdevへ統合されたauthorityを候補が取り消す差分に見えるためである。checkerのseparationを弱めず、
S00開始HEADでも同failureを再現して原因を確認した。S00のscope限定PASSをbranch全体のlatest-dev gate PASSとは呼ばない。
ユーザー指定通りS00ではactive authorityを変更せず、dev mergeとその後のbranch gate再検証をruntime前の未完了条件に残す。

S00ではbrowser、性能、codec、Python scientific outputs、schema compatibilityの実装変更がない。
browser/large Session/heap/heartbeat/performance/replayは **未測定・未実施** であり、S01以後に属する。
53 fixturesはchecker mechanicsの証拠で、完成したBUG-02・BUG-20のacceptance証拠ではない。

## Diff review / architecture / acceptance

- Production: `gbdraw/` とruntime経路の変更ゼロ。既存discovery/request/validation/Result/History ownerは不変。
- Tests/checker: `tests/`、`tools/`、`.github/` のtracked変更ゼロ。candidate verifierはinert patch検証だけを担当する。
- Active authority: 実装checkoutのOIPC/map/BD/policy/architecture rules変更ゼロ。temporary devのcandidateだけを検証した。
- Docs: BUG-01専用計画の削除、二件だけの引継ぎ、候補・統合手順・本結果をreview。二件の承認本文は不変。
- Generated assets: wheel、Gallery、reference outputs、social preview、dist/egg-infoの変更ゼロ。
- OE/PE/CB: runtime owner/path/compatibilityの追加・変更なし。新bindings readerと、その例外前提を削除。現scopeにcomplete exception setsは不要。
- A-01: candidate scope/既存detector/authority分離を検証。runtime owner・cycle・性能acceptanceは未着手。
- W-01: 専用checkout、同名branch/upstream、明示stage、一commit、non-force pushとremote/local一致確認を行う。
- D-01〜D-04 / S-01〜S-05: outcomeは候補へ忠実に保存。実装・受入試験は未実施。

## Remaining integration conditions / next entry

Product-onlyとpermission-onlyの各target・PR base・適用手順・統合順序は候補READMEに固定した。
Product判断の再承認は不要だが、別authority PRのreview/mergeはこのsessionの許可範囲外。
Product dev merge SHA、privileged dev merge SHA、implementationへのdev取り込みはいずれも未完了。
permission candidateはS01のtransport/path evidenceを確認した後に統合する。

次sessionは [S01_INSTRUCTION_PROMPT.md](../sessions/S01_INSTRUCTION_PROMPT.md) を読み、独立checkoutで測定だけを行う。
その独立evidence作業はauthority未mergeでも進められるが、S03/S04/S05のruntimeは必要authorityとbase-readinessを待つ。
runtime前は `git fetch origin` → authority merge SHAのdev ancestry確認 → clean treeで `git merge --no-edit origin/dev` → required policy/gates再検証。
新二concernを無理にBD/mapへ登録せず、mapped contractを変える必要が生じれば既存のevidence-only laneを先行させる。
このsessionはS00で終了する。

## Commit / remote handoff

English commit title: `Document approved issue 597 authority integration`

English summary: Remove BUG-01 from the implementation plan and prepare separately validated Product and minimal Worker-permission patches for BUG-02 and BUG-20. Record authority-only integration steps and pending dev prerequisites without changing active authority or runtime.

結果文書を含むcommit自身のSHAはここへ自己参照として書き込まない。
同名remoteへのnon-force push後、`git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926` と
`git rev-parse HEAD` の一致、clean `git status --short --branch` をsessionの最終応答で報告する。
