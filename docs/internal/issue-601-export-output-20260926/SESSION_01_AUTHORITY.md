# INSTRUCTION PROMPT — S01 承認済み製品動作の正式契約化

あなたはgbdraw Issue #601のProduct authority担当です。日本語・中国語のPDF文字対応と、原因のあるsummary＋任意の安全なDetailsの二つのoutcomeが、`satoshikawato`に承認されています。承認内容を独立recordとして正式なstatic Product Contractへ記録し、runtimeより先にbaseへ統合できる差分を準備してください。

## branch取得と必読資料

SESSION_IDは`s01`です。まず[総合計画書](./MASTER_PLAN.md)の共通取得手順で**`fix/issue-601-export-output-20260926`**をremoteから取得し、専用cloneで最新計画と`SESSION_00_RESULT.md`を読みます。他sessionのcheckoutやrefsを操作しません。

AGENTS.md、CLAUDE.md、Web CLAUDE、Product Impact Ratchet、Architecture Ratchet、両Decision Packs、[承認の機械表現](./APPROVED_DECISIONS.md)、origin/devのOIPCとProduct Impact map/storeを確認します。PDF recordはS00 required evidenceの合格が開始条件です。error recordは独立した承認・契約として扱います。

## 専用authority branch

Product Contract変更は他fileから隔離する規定があります。実装計画branchの文書を正式契約と同じPR差分に混ぜません。実装branchを取得した専用clone内で、最新`origin/dev`から別のauthority worktree/branchを作ります。

```bash
git fetch origin
AUTHORITY_ROOT=$(mktemp -d /tmp/gbdraw-issue601-authority-XXXXXX)
git worktree add --detach "$AUTHORITY_ROOT/repo" origin/dev
cd "$AUTHORITY_ROOT/repo"
git switch --no-track -c product/issue-601-decisions-20260926 origin/dev
```

このbranch名がremoteに既存なら実際のremote状態と対象commitを確認します。自分の先行作業と確認できる場合だけ取得して継続し、他者のbranchを再作成・上書きしません。

## 作業所有範囲

authority branchでは**`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`一ファイルだけ**を変更します。runtime、tests、guard、policy、他のMarkdown、JSON storeを変更しません。PDF対応とdiagnostic disclosureを、同一concernへまとめず独立recordとして記録します。

## 実施内容

1. trusted baseのactive authorityを再検索します。今回のconcernは基点のmap/storeに未登録ですが、並行PRで登録された場合は競合する第二authorityを作らず、同じ製品outcomeを既存モデルの規定で扱います。不存在のBD番号を記載しません。
2. OIPCの現行revision、record ID、metadataを確認して衝突しない新規recordを採番します。`SELECTIVE_JA_ZH_VECTOR_TEXT_FALLBACK` / scenario 4 と`ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS` / scenario 1を独立して記録します。
3. 各Decision Packの承認済みRationale、Must preserve、May retire、Accepted residual risk、owner=`satoshikawato`、date=`2026-09-26`を忠実に反映します。推測したrationale/risk、未承認の代替outcomeを追加しません。条件のANDを維持し、実装owner/pathを製品outcomeと混同しません。
4. PDFにはS00の再現可能evidenceとrequired evidence境界を記録します。errorのrequired acceptanceはsource→Worker→UI、safe allowlist、未知cause、rollback/retry、raw非公開を明示します。Product承認がGate/dependency/privilegeの免除ではないことを保持します。
5. [機械表現](./APPROVED_DECISIONS.md)と承認本文を照合し、OIPCに記録したoutcomeをレビュー可能な形で提示します。S00の証拠branchへの恒久commitリンクを使う場合、そのSHAとpathが実際に存在することを確認します。
6. trusted base checkoutのcheckerを用いてauthority-only差分を検証します。OIPC以外にchanged pathがないこと、Product/runtime/guard変更が混在しないことを確認します。

## 独立concernが保留された場合

PDF証拠の不足はPDF recordとそのruntimeだけを止めます。errorの承認outcomeは独立して記録・レビューできます。この場合はerror recordだけのauthority commitを作り、PDF recordを未承認代替outcomeで埋めません。PDF recordは証拠完成後の別authority-only変更で記録します。

## 完了とruntimeへの引き渡し

authority差分を**commitし、`product/issue-601-decisions-20260926`へpushしてください**。branchとupstreamを確認して`git push -u origin HEAD:refs/heads/product/issue-601-decisions-20260926`を使い、remote SHA一致を確認します。
Commit title例: `Define Issue 601 PDF and diagnostic disclosure contracts`。

レビュー・PR作成・devへのmergeはその操作の明示的な許可に従い、devへ直接pushしません。PR文面を作る場合はwrite-clear-pull-request skillと文面checkerを使用します。authority-only PRがmerged origin/devに存在することをruntime開始条件にします。authority branchのpushや未merge候補はruntimeを許可しません。

branch、commit SHA、OIPCのrecord IDs、concern/option/revision、checker結果とmerge状態をhandoffしてください。承認は本promptから参照できるため、authority branchへ別result fileを混ぜる必要はありません。実装branch側のhandoff記録はS02または独立したdocs-only commitで行います。
