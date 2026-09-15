# INSTRUCTION PROMPT — S02: 経路の公開・保存契約と改修設計

このセッションでは、全経路列挙の指数増加を解消するための証拠・consumer inventory・具体的設計と、必要なProduct判断資料を完成させてください。グラフ表現を採用済みと仮定してruntimeを変更しないでください。

## 必読と前提

- [総合計画書](MASTER_PLAN.md) §2〜§9。
- S01 handoff `results/S01.md` と実在する再現runner・baseline。
- repositoryガイド、`PRODUCT_IMPACT_RATCHET.md`、`PRODUCT_DECISION_PACKET_TEMPLATE.md`、`ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`。
- base branchの比較・保存契約、OIPC、関連する公開APIとrelease/main履歴。

S01が必要条件です。S03〜S05の完了は不要です。判断待ちになっても、それらの独立した改善を止める必要はありません。

## 担当範囲

- 経路のpublic / persisted / internal consumer調査。
- テスト専用prototype・小規模比較・count/ID計算の検証。
- `results/S02.md` と、必要な場合の `results/PATH_DECISION_PACK.md`。
- Product選択を受領した場合のauthority-only表現。ただし適用する受理経路・対象は既存policyに従う。
- 本番runtime・新writer・旧readerの実装はS06の担当。

## 実行内容

1. `_build_ortholog_paths()` のstart/terminal、cycle回避、dedup、sort、path連番、edgeへの最初のpath、sharedProteinIdsの規則を明文化する。
2. 到達可能なedge構築経路からDAG性を検証する。record順の単調性を根拠なく仮定しない。same-record、coortholog、重複edge、複数source/sink、逆向きevidenceを含める。
3. `OrthogroupResult` / `OrthologPath` / `OrthologEdge` の全consumerを追う。Python戻り値、dataclass型、typed encode/decode、SVG属性、block metadata、catalog縮約、popup、derived cache、Sessionを一覧化する。
4. consumerごとに必要情報を分類する: 全列挙、正確なcount、特定path、edgeのpath ID、shared情報、単なるtransport。通常renderでの隠れたlist化・再serializeも追う。
5. 実際に公開されたnamespaceとformatをmain first-parent / release tag / positive fixtureで裏付ける。branch-only形式に互換readerを追加する設計にしない。
6. `PATH-A / exhaustive-current` と `PATH-B / lossless-graph` を具体化する。後者では通常の保存とrenderが全列挙を必要としないデータ契約、既存consumerの継続、exact count、ID再現、型とerrorsを示す。
7. 小さい入力でgraph/count/順位案を旧結果と比較する。protein列によるdedupとedge列の選択も一致させる。DAG countだけのprototypeを全metadata同値の証明として使わない。
8. 全量取得の `Ω(L)`、大整数の計算・JS転送、旧tuple契約の維持コスト、循環の扱い、compatibility数を比較表へ記載する。
9. authority searchとpreflightを行う。既存authorityで完全に決まるなら根拠を示す。material outcomeが複数残るならDecision Packを作る。単に別classを選ぶような実装上の判断だけでProduct承認を要求しない。

## Product判断が必要な場合

安定choiceは総合計画書の `PATH-A` / `PATH-B` を使います。どちらも適法・実現可能かを調査結果で確認し、禁止される選択肢を承認可能として提示しないでください。

Decision Packにはユーザーjourney、現行authority、must preserve、may retire、API/保存/取得/失敗の差、compatibility、測定、rollback、engineering recommendationを含めます。推奨は選択ではありません。

Product Decision Ownerへ判断を求める場合はrepository templateに従い、以下の完全な回答を受け取ります。

```text
PRODUCT_DECISION
Concern: <concern key/title>
Scenario revision: <revision>
Choice: <stable choice code and outcome ID>
Rationale: <product-level reason>
Must preserve: <effects and affordances>
May retire: <none or explicit scope>
Accepted residual risk: <bounded risk or none>
Owner: <maintainer identity>
Decision date: <YYYY-MM-DD>
```

欠けた理由・廃止範囲・risk・owner・dateを補完しないでください。受領済みの明示的選択だけを機械表現にし、review可能な形で提示してください。必要なauthority受理・base統合はruntimeより前です。candidate authorityでcandidate runtimeを認可しないでください。

## 完了条件と引き継ぎ

- 全consumerと必要情報、保存namespace、DAG/重複/順位の検証結果が記録される。
- S06が実装できる具体的契約と、小規模の同値oracleがある。
- 判断待ちならその状態と不足事項を明記する。S02の調査完了とProduct選択完了を別々に報告する。
- PATH-Aの場合は指数出力量が残ること、PATH-Bの場合も明示的全量取得は出力量に比例することを記載する。
- `results/S02.md` にauthority source、decision revision、baseへの反映状態とS06の開始条件を記録する。
