# Authority integration candidates — BUG-02 / BUG-20

候補準備は完了。dev の authority 統合は未完了。
このディレクトリは inert patch と検証手順を保存する。active authority、第二の decision registry、runtime ではない。
BUG-01 はユーザーの2026-09-26の指示で今回の計画から除外した。

## 正規保存先と候補

調査した trusted base は `origin/dev` = `302dfa1136c95ab50ddef606ad835ec4087b95b2`。

| 責任 | 候補 | authority-only PR の全変更path | 内容 |
| --- | --- | --- | --- |
| Product outcome | [01-product-contract.patch](./01-product-contract.patch) | `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` のみ | OIPC revision 19 → 20。`PD-OI-044` / `PD-OI-045` に二件の完全な receipt を忠実に追加。 |
| Worker constructor の許可 | [02-import-worker-permission.patch](./02-import-worker-permission.patch) | `tools/web-change-policy.json` のみ | `allowedPrivilegedOwners["Diagram Worker"]` に `services/session-import-client.js` 一件だけ追加。 |

`tools/web-product-impact-map.json` は canonical render request と saved Session regeneration の二 concern を持つ。
今回の `diagram-generation.circular-transform-discoverability` と `web.session-operation-consistency` は未mappedで、
`tools/web-product-decisions.json` の active decisions も空である。
Product Impact Ratchet の unmapped static authority 経路に従い、既存 OIPC に追加する。
`BD-###`、map concern、scenario、detector、workflow、architecture rule を新設しない。
既存 `PD-OI-010` の一件 crop 条件、`PD-OI-019` の fresh canvas default、
`PD-OI-016` / OIPC-C07 の failure isolation は置き換えず、二件の receipt と合わせて維持する。

候補の provenance は元承認記録を含む commit `51a786086dc2777e5fa375e5f94d8b7ac7deeedc`、
承認者・日付・承認原文・承認対象 digest に固定した。判断の再承認は不要。
候補の番号と revision はこの dev base に対する割当であり、merge 時に競合していれば最新devの未使用番号へ調整し、
二件の全9項目を維持して再検証する。先行 authority を上書きしない。

## 最小 permission の根拠

既存 `tools/web-architecture-detectors.mjs` の `Diagram Worker` operator は `new Worker` を検出する。
これは detector 上の分類名であり、Session import に Python runtime を追加する意味ではない。
予定構成は一方向の `config.js → session-file.js → session-import-client.js → session-import-worker.js`。
constructor と cancellation/stale/settlement/teardown は一 client が所有し、Worker は JSON/gzip codec のみを担当する。

| Subject / edge | 現行detectorの結果 | 候補で必要な許可 |
| --- | --- | --- |
| `services/session-import-client.js` の一 `new Worker(...)` | `Diagram Worker` owner 一件 | 上記一件 |
| `session-file.js → session-import-client.js` | target は既存 privileged import target の集合外 | 追加なし |
| client の `new URL('../workers/session-import-worker.js', import.meta.url)` | privileged import edge ではない | 追加なし |
| Worker 内の File read / gzip decode / JSON parse / reply | 既存 privileged operator を使わない | 追加なし |

client/worker は `session-file.js`、config、History、canonical request、SVG admission、Pyodide helper を逆importしない。
特に codec 共用のために Worker から `session-file.js` をimportする並行ownerや循環経路を作らない。
既存の `services/session-file.js` の importer (`services/config.js`) と全54 importer permissions は不変。
operator permissions は48 → 49、追加一件・削除ゼロ。
`Session` permission、worker file の blanket permission、new importer key は追加しない。
これは予定構成の detector probe であり、完成runtimeの適合・応答性・memoryの証拠ではない。
S01がtransportと最終pathを確定した後、実際のsubjectsを再確認してpermission-only候補をreviewする。
別subjectが必要になればその具体的scopeを別deliveryで扱い、gateやdetectorを緩めない。

## 適用・検証

この session は別PR・別branchへのpush・mergeを行っていない。
maintainer は候補commitを参照し、最新 `origin/dev` から責任別の新規branchを作る。
Product例: `product/issue-597-authority-dev-20260926` → PR base `dev`。
permission例: `policy/issue-597-session-import-worker-20260926` → PR base `dev`。
implementation branch をbaseにせず、そのcommitsをcherry-pickしない。各patchの一targetだけをstageする。

再現にはこの計画branchの独立checkoutと、そのcloneに属する検証専用dev worktreeを使う。
次のcommandsはローカルの検証用で、remoteへは書き込まない。

```bash
git fetch origin dev:refs/remotes/origin/dev
ISSUE597_CANDIDATES="$PWD/docs/internal/issue-597-input-session-implementation-20260926/authority-candidates"
ISSUE597_DEV_WORKTREE=$(mktemp -d /tmp/gbdraw-issue597-authority-check.XXXXXX)
git worktree add --detach "$ISSUE597_DEV_WORKTREE" origin/dev
cd "$ISSUE597_DEV_WORKTREE"
git status --short --branch
git apply --check "$ISSUE597_CANDIDATES/01-product-contract.patch"
git apply "$ISSUE597_CANDIDATES/01-product-contract.patch"
python "$ISSUE597_CANDIDATES/verify_candidate.py" product "$PWD"
git diff --check
node tools/check-web-change-budget.mjs --base origin/dev
git apply --reverse "$ISSUE597_CANDIDATES/01-product-contract.patch"
git status --short
git apply --check "$ISSUE597_CANDIDATES/02-import-worker-permission.patch"
git apply "$ISSUE597_CANDIDATES/02-import-worker-permission.patch"
python "$ISSUE597_CANDIDATES/verify_candidate.py" privileged "$PWD"
git diff --check
node tools/check-web-change-budget.mjs --base origin/dev
git apply --reverse "$ISSUE597_CANDIDATES/02-import-worker-permission.patch"
node --test tests/web/architecture-ratchet-fixtures.test.mjs tests/web/product-impact-ratchet-fixtures.test.mjs
git status --short
cd -
git worktree remove "$ISSUE597_DEV_WORKTREE"
```

二つのpatchを同じauthority PRへまとめない。OIPC変更は一pathだけの分離がcheckerで要求される。
各candidateの実測結果は Gate PASS / Review REQUIRED。53 focused checks は全pass。
Review REQUIRED はhuman reviewの条件であり、dev mergeを済ませたという意味ではない。
最新baseのrevision/ID/schemaが変わったら、適用だけ成功しても旧検証を再利用しない。

## 統合順序と未完了条件

1. Product-only patchを最新devの別PRでreview・mergeし、`PD-OI-044` / `PD-OI-045` の実際のIDとdev merge SHAを記録する。
2. S01はauthority未mergeでも独立measurementを実施できる。transport/pathの選定と性能・heap・転送evidenceを保存する。
3. そのpathを確認してpermission-only patchを別PRでreview・mergeし、dev merge SHAを記録する。S05開始前の必須条件。
4. runtime前にclean implementation checkoutで `git fetch origin`、必要SHAのdev ancestryを確認し、`git merge --no-edit origin/dev` で取り込む。S03はdiscovery authority、S04はSession authority、S05はpermissionとS01方式に依存する。
5. mapped contract変更が実際に必要なら evidence-only merge → authority ref-only merge → runtime の順序に追加する。現時点で二つの新concernのmap/contract変更は不要で、既存canonical request/admissionのowner/edgeは維持する。

checkerのexact-path enforcement、static Product Contract、privileged owner機構は既にdevにあるため、checker-only laneは現時点で不要。
将来checker mechanics変更が必要と判明した場合は checker-only → authority-only → runtime を守る。
Product/permission merge SHAはいずれも **未取得・未統合**。
候補を作成・検証・implementation branchへpushしても、dependent runtime開始条件は満たさない。
S00でbrowser性能、codec、UI、runtimeを検証済みとは報告しない。

詳細な検証・開始SHA・差分review・dev未取り込みのboundaryは [S00_RESULT.md](../results/S00_RESULT.md) に記録する。
