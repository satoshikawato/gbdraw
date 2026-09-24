# Issue #581 — S00 Product preflight INSTRUCTION PROMPT

以下を一つの独立した作業指示として使用する。過去の会話は前提にしない。

## 依頼

gbdraw Web の Issue #581 は、feature popup で Record actions が常時大きく表示される問題、
Similarity Alignment の候補が内部 ID 中心で読みにくい問題、radio click ごとに Python
Resolver を呼んで選択が止まる問題を扱う。[設計提案コメント](https://github.com/satoshikawato/gbdraw/issues/581#issuecomment-5811037046)
は、Record actions の意図的な開閉、ローカル選択と Apply 時の一括検証、non-modal palette、
縦ガイドと候補番号、図上選択を提案している。

このセッションは製品結果と権威の確認だけを行う。runtime code は変更しない。

2026-09-24 に Product Decision Owner `satoshikawato` は３件の完全な
`PRODUCT_DECISION` 文面をそのまま承認した。`A / EDIT_DISCLOSURE`、
`A / LOCAL_BATCH_RETRY`、`A / FLOATING_GUIDE_CANVAS_PICK` である。
`issue-581-product-decisions-20260924` の commit `18903342` は、これらを
PD-OI-033〜035 として authority-only 文書に記録したローカル候補である。
`origin/dev` にマージされるまで依存 runtime を開始しない。承認済み３件を
再判断に戻さず、追加の material な結果差だけを新たに審査する。

## ブランチと必読資料

- 固定 Issue #581 実装ブランチは `issue-581-feature-popup-alignment-ux-20260924`。
  総合計画は `docs/internal/issue-581-feature-popup-alignment-ux/IMPLEMENTATION_PLAN.md`。
- まず `git branch --show-current`、`git rev-parse HEAD`、`git status --short --branch`、
  upstream、`origin/dev` の SHA と祖先関係を確認する。無関係な変更を保持する。
- `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、
  `docs/internal/PRODUCT_IMPACT_RATCHET.md`、
  `docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md`、
  `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`、
  `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の
  PD-OI-026〜029、031、032 と authority 候補の PD-OI-033〜035、
  `tools/web-product-impact-map.json`、
  `tools/web-product-decisions.json` を読む。
- Issue #581 の本文・全 comments と最新 `origin/dev` の該当 code/test を確認し、
  base の accepted authority と単なる提案・実装事実を分ける。

## 作業

1. 次を別々の product concern として、利用者の操作、途中状態、失敗、次の操作まで比較する。
   - popup の Record actions の配置・初期開閉・欄内 Cancel。候補は Edit 内の折り畳み欄、
     専用 Record tab、現状維持。rich と simple の到達性も比較する。
   - 候補の biological label と local Select／Skip、Apply 一括検証、検証・生成失敗後の draft。
     最終 plan の意味と既存の自動適用を維持できるか調べる。
   - modal から non-modal palette への変更、pan／zoom、guide、badge、直接図上 click、
     他の編集で開始時 artifact が変化したときの扱い。改良 modal も製品上の代替結果として比較する。
2. 各 concern を `IMPLEMENT_EXISTING_AUTHORITY`、`EVIDENCE_REQUIRED`、
   `PRODUCT_DECISION_REQUIRED`、`NOT_ALLOWED` の一つに分類し、他の分類が不適切な理由を書く。
   既存の PD-OI-026 の exact 選択と独立性、PD-OI-032 の record 回転契約を維持する。
3. 承認済み３件の外で新たに証拠や判断が必要な concern には総合計画の隣に Decision Pack を作り、
   stable な結果選択肢、維持・追加・失う効果、Decision route、contract、証拠の限界、
   Product Decision Owner 用 `PRODUCT_DECISION` response template を記す。
   Issue コメントを正式な receipt に読み替えない。rationale、may retire、risk、
   owner、date を推測して埋めない。
4. 承認済み３件を記録した authority-only commit `18903342` の内容と最新 base を
   確認する。candidate authority と依存 runtime を同じ候補で承認させない。
   authority が `origin/dev` にマージされた後にだけ、固定 runtime ブランチへ反映する。
   push、PR、merge は別途明示許可がある範囲だけ実施する。
5. 結果、根拠、残る判断、S01〜S03 の各開始可否を総合計画の実施記録へ記入する。

## 終了条件

各結果について「どの既存契約が許すか／何が未決か」が、会話なしで分かる。
未決の runtime は開始可能と報告しない。独立に実装可能な範囲は明示する。
設計上の原則は、製品結果と実装 owner を分ける SOLID、concern ごとに最小の判断を行う
KISS、既存 authority を重複させない DRY、汎用 UX policy や新たな判断基盤を作らない
YAGNI である。英語の proposed commit title と短い summary を示す。
