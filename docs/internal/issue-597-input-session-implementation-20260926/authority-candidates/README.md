# Authority integration candidates — BUG-02 / BUG-20

Product authority 二件は PR #610 / dev merge `af5d942af60353dda199aa487da9152a3576b3fe` に統合済み。
OIPC revision 21、PD-OI-044/045 と二件の receipt は S01 で全18項目の一致を確認した。
[AUTHORITY_INTEGRATION_RESULT.md](../results/AUTHORITY_INTEGRATION_RESULT.md) と
[S01_RESULT.md](../results/S01_RESULT.md) を参照。Product 選択の再承認は不要。

このディレクトリは inert 候補であり、active authority や runtime ではない。
`01-product-contract.patch` は S00 時点の revision 19 に対する履歴 artifact。
統合済み dev へ再適用しない。旧検証・Product integration 手順は [S00_RESULT.md](../results/S00_RESULT.md) に保存している。
BUG-01、新 bindings、S02 は scope 外のままである。

## S01 で確認した最小 permission

[02-import-worker-permission.patch](./02-import-worker-permission.patch) の適用先は
`tools/web-change-policy.json` **一ファイルのみ**。constructor owner 一件と codec importer 一件を追加する。
S00 の「importer 追加なし」という予定を、実際の disposable Worker の codec 共用に合わせて更新した。

| Subject / edge | 現行 detector | 必要な候補 |
| --- | --- | --- |
| `services/session-import-client.js` の一 `new Worker` | Diagram Worker owner | allowedPrivilegedOwners に一件 |
| `workers/session-import-worker.js → services/session-file.js` | privileged import edge | session-file の allowedPrivilegedImporters に一件 |
| `config.js → session-import-client.js` | target は privileged target 集合外 | 追加なし |
| client の Worker URL | privileged import edge ではない | 追加なし |

runtime の予定は `config.js → session-import-client.js → session-import-worker.js → session-file.js`。
`session-file.js` は既存の純粋な File/gzip/fatal UTF-8/size-limit codec owner として残し、client を逆 import しない。
config は Save の既存 codec import を維持する。client が operation ID、stale rejection、error、settlement、termination を所有する。
Worker は config、History、request、SVG admission、Python helper を import しない。
codec validation を複製せず、reply は未信頼の plain candidate として既存 preflight/adoption に渡す。

現行 permissions は operator 48 → 候補49、importer 54 → 候補55。
追加二件、削除ゼロ。別 owner、Session operator、blanket Worker permission、new importer key、detector/checker の変更はない。
S01 で active policy は変更していない。方式の性能判断は S01 結果だけを参照し、permission 候補の存在を性能 PASS と扱わない。

## 候補の再検証

S01 の独立 checkout で、active policy を変更せず次を実行した。
`--inert` は一時 directory の policy **コピー**にだけ patch を適用し、既存 detector と形状を検証する。

```bash
ISSUE597_CANDIDATES="$PWD/docs/internal/issue-597-input-session-implementation-20260926/authority-candidates"
git apply --check "$ISSUE597_CANDIDATES/02-import-worker-permission.patch"
python "$ISSUE597_CANDIDATES/verify_candidate.py" privileged "$PWD" --inert
```

結果: 一 constructor owner、一 codec import edge、候補の一 target / 二 permission 以外の差分ゼロ。
別 owner は不許可。active policy SHA は測定 source / fetched dev と同一。
この局所検証は authority PR の review / merge や trusted CI gate の代わりではない。

## 別 authority PR の開始条件

S01 結果で transport の成立とこの paths を確認してから、maintainer が最新 `origin/dev` を base に
permission-only branch / PR を用意する。implementation branch の docs/tests commits を cherry-pick しない。
この session では authority PR の作成・push・merge を行っていない。

```bash
git fetch origin dev:refs/remotes/origin/dev
git switch --no-track -c policy/issue-597-session-import-worker-20260926 origin/dev
git apply --check "$ISSUE597_CANDIDATES/02-import-worker-permission.patch"
git apply "$ISSUE597_CANDIDATES/02-import-worker-permission.patch"
python "$ISSUE597_CANDIDATES/verify_candidate.py" privileged "$PWD"
git diff --check
node tools/check-web-change-budget.mjs --base origin/dev
```

この手順は**将来の別 authority delivery**向けであり、S01 implementation checkout では実行しない。
その PR は上記 policy 一ファイルだけを stage し、review / merge SHA を記録する。
最新 base の subjects / policy が変われば候補の適用成功だけで旧 evidence を再利用しない。
必要な mapped contract が変わる場合は既存 evidence-only → authority-ref-only → runtime の順序を維持する。

Worker permission の dev merge SHA は **未取得・未統合**。
S05 は S01 の方式成立、permission-only dev merge と implementation への取り込み、既存 gate PASS を必要とする。
Product authority は統合済みであり、その判断を再選択しない。
