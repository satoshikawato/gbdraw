# Issue #597 authority integration result

日付: 2026-09-26。対象はBUG-02とBUG-20。BUG-01の計画は削除済み。
S00の候補準備と、今回のdevへのProduct authority統合を区別して記録する。

## 完了した統合

Product authority二件は [PR #610](https://github.com/satoshikawato/gbdraw/pull/610) でdevへmerge済み。
統合時刻は `2026-09-26T13:22:43Z`。

| 項目 | 確定値 |
| --- | --- |
| authority work branch | `product/issue-597-authority-dev-20260926` |
| authority base | `2edc00aebc74e01003da643dfc957b513d5dcfe5` |
| authority commit | `c394924d545530086632b64eeb5d8619f7d1ce18` |
| dev merge commit | `af5d942af60353dda199aa487da9152a3576b3fe` |
| implementationへのdev取り込み | `b5c9d55039002944b6f6e3383bf4cc214ba3eb56` |
| implementation branch | `fix/issue-597-input-session-20260926` |
| 正規保存先 | `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` |
| Contract revision | `20` → `21` |
| discovery判断 | `PD-OI-044` / `diagram-generation.circular-transform-discoverability` |
| Session判断 | `PD-OI-045` / `web.session-operation-consistency` |

[discovery receipt](../decisions/02_RECORD_DISCOVERY.md) と
[Session receipt](../decisions/03_SESSION_OPERATIONS.md) のhuman/serialized/authorityの9項目は一致する。
選択、理由、must preserve、may retire、accepted residual risk、owner、日付を変更していない。
承認対象資料のSHA-256と、receiptを保存したcommit `51a786086dc2777e5fa375e5f94d8b7ac7deeedc` をauthorityに記録した。
先行devのrevision 20と既存43判断、lifecycle、cross-surface clauses、acceptance/risk文言を保持した。
S00のrevision 19向けpatchは履歴artifactとして保持し、統合済みdevには再適用しない。

PR #610の変更はContract一ファイルのみ。Product map、decision registry、checker、privileged permission、runtimeを変更していない。
implementationへのmergeは先行devの変更をそのまま取り込んだ。今回独自にruntimeを実装していない。
S00で再現した旧implementationと最新devのauthority差異は、このdev取り込みにより解消した。

## 検証

独立dev worktreeで、既存記録の保持と二件の承認内容を検証した。

```bash
python docs/internal/issue-597-input-session-implementation-20260926/authority-candidates/verify_candidate.py product "$ISSUE597_AUTHORITY_WORKTREE"
node tools/check-web-change-budget.mjs --base origin/dev
node --test tests/web/architecture-ratchet-fixtures.test.mjs tests/web/product-impact-ratchet-fixtures.test.mjs
```

Product candidateはPASS、policyはGate PASS / Review REQUIRED、focused fixturesは53件PASS・skip 0。
Review REQUIREDはauthority変更の明示的なreview要求であり、runtime admissionではない。
PR wording checkもPASS。
PR #610のtrusted Web base policy、Web change budget、CI impact plan、PR / gate、CodeQL、Workers Buildsは成功した。
対象外CI jobsはdocs-only impact planでSKIPPED。検証記録は
[CI run](https://github.com/satoshikawato/gbdraw/actions/runs/36244606758) と
[trusted policy run](https://github.com/satoshikawato/gbdraw/actions/runs/36244606745) を参照。
head SHAを固定した通常mergeを使用し、bypass/force pushを行っていない。

取り込み後のimplementationは `node tools/check-web-change-budget.mjs --base origin/dev` でGate PASS / Review CLEAR。
active Contractは取得devとbyte-identical、dev merge SHAは取得devとHEADのancestor。
取得devとの差分は計画ディレクトリのみで、runtime・checker・tests・active authority差分はない。
計画の相対リンク、fences、bash構文、Python AST、receipt保持、候補digest、`git diff --check` もPASS。
production/tests/generatedに今回の独自変更がないことと、handoff docsの差分を分けてreviewした。

## 未完了条件と次の入口

- Worker permissionは**dev未統合**。[S00候補](../authority-candidates/README.md) の一constructor owner追加は、S01でtransport成立と最終pathを確認してから別authority PRで統合する。
- この条件をS01のmeasurement開始前の追加permissionとはしない。disposable probeをapplication routingへ接続しない。
- runtime実装前にはProduct authorityだけでなく、必要なWorker permissionのmerged base適用を確認する。候補authorityで同じ候補runtimeを許可しない。
- S01以降のmeasurement/runtime実装は未実施。次は [S01 prompt](../sessions/S01_INSTRUCTION_PROMPT.md) のbaselineとtransport evidence。S03以降を連続実行しない。

English commit title: `Record issue 597 Product authority integration`

English summary: Record the merged Product decisions and dev incorporation, preserve S00 evidence, and make the remaining Worker permission condition explicit for S01.
