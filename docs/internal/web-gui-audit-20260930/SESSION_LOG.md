# 実行状況

計画: [03_IMPLEMENTATION_REFERENCE.md](03_IMPLEMENTATION_REFERENCE.md)。判断: [02_DECISION_PACK.md](02_DECISION_PACK.md)。

| 日時 | PR | merge SHA | 直した ID | テスト | 生成物 | 残った問題 |
|---|---|---|---|---|---|---|
| 2026-09-30 | — | — | — | — | — | 監査（README）、修正提案（01）、Decision Pack（02。Owner が確定）、実装リファレンス（03）を作成した。runtime、Product Contract、CI は未変更。基準は origin/dev 4c89bab1 |
| 2026-09-30 17:28 JST | #649 (P00) | e97d90fe | — | docs tests 21 passed | CHANGELOG `[Unreleased]` に PR ごとの slot（P03〜P20）を追加 | なし |
| 2026-09-30 17:47 JST | #650 (P01) | a0cb65eb | — (authority) | 02 の 30 receipt が verbatim、SHA-256 と JSON が一致することを script で照合。architecture-contracts 0 fail | Product Contract revision 29: PD-OI-056〜084、PD-OI-018 revision 4、OIC-027 | N-15 は Product Contract と同じ PR に入れられない（isolation 規則）ため、別の docs PR にした |
| 2026-09-30 18:14 JST | #651 (N-15) | bfa18b8a | N-15 | ci 229/229、docs contracts 21 passed。コード側（smoke 上限 19、ratchet の MAXIMUM_RULE_COUNT=4）が正しく、文書を合わせた | なし | 監査に無かったずれ（gallery parity command は web-pr-smoke ではなく gallery job で走る）も直した |

## 進め方と Owner-delegated の記録

表の行は PR をマージするたびに追記する。orchestrator が SESSION_LOG だけを変える docs PR でまとめて dev に入れる（並列の PR が同じ行を触って衝突しないようにするため）。

- 2026-09-30（進め方）: CHANGELOG の `[Unreleased]` に PR ごとの slot を置き、各 PR は自分の slot だけを置き換える。dev の branch protection が `strict`（最新への追従が必須）なので、マージは直列になる。
- 2026-09-30（進め方）: P02 と Wave 1（P03〜P06、P17、P18）を同時に実装する。P02 の known-defect の印は Wave 2 以降（P07〜P16）の ID に付け、Wave 1 の PR は自分の ID の失敗するテストを自分で書く（03 §4.1 の「失敗するテストを先に書く」は各 PR で満たす）。
- 2026-09-30（進め方）: N-15 は P01 から外し、別の docs PR にした。Product Contract の変更は他のパスと同じ PR に入れられない（`tools/check-web-change-budget.mjs` の isolation 規則）。
