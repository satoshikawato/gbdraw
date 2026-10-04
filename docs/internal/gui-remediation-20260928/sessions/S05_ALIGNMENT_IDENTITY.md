# S05 INSTRUCTION PROMPT: BGCの直接RBHと生物学的な表示名

保存済みの直接RBHをWeb alignmentのresolverへ渡し、livA基準でracMが自動選択されるよう修正してください。真の曖昧性のReview、既存Keep、atomic Historyを維持し、plan表示を配列・featureの名前に直します。

## 場所と入力

**fix/gui-feedback-remediation-20260928**、/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928を再利用。
[総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、results/S00.md、AGENTS/CLAUDE/Web CLAUDEを読む。
主fixtureはgbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json。

resource comparison-canonical-orthogroups-1をdecodeしたvalue.fieldsのorthologEdgesByOrthogroupId.og_18とmemberByProteinIdに根拠がある。

| gene | fixture上のrecord | protein ID | 根拠 |
| --- | --- | --- | --- |
| livA | record-1/BGC0000708 | CAG38712.1 | exact reference |
| racL | record-5/BGC0000713 | CAG34708.1 | livAからcoortholog |
| racM | record-5/BGC0000713 | CAG34707.1 | livAからrbh |

record-Nはこの固定fixtureの説明用。runtimeの照合に表示順を使わない。

## 実装責任

1. 既存diagram-generation-workerのRESOLVE_SIMILARITY_ALIGNMENTが準備するcanonical resourcesを利用する。
2. web_support/similarity_alignment.pyの_projection_contextとtyped decodeで得たorthogroup情報から、選択groupの直接edgeを取得する。
3. api/record_planning.pyのCLI側anchor/edge projectionと実際に共通化できる部分だけ共有し、既存resolve_similarity_alignmentへ渡す。
4. graph endpointのsource/record/protein/featureを現在のrecordKeyへ正確に対応付ける。materialized requestを取得しただけでrecordIndexのreorderが正しいと仮定しない。
5. group違い、cropによる候補無効、重複source/同名feature、stale responseを検証する。graph evidenceがないときは推測しない。
6. UI catalogへ全graphを戻さず、LOSATを再実行しない。新ambiguity policy、recommendation自動昇格、bitscore順位、rationale/schema変更を導入しない。
7. inspectActivePlanとreviewの表示helperを共有し、record名/accession、gene/locus_tag/protein IDを表示する。欠落はtype/座標、同名は識別情報を補う。内部hashはbinding専用。
8. index.htmlの必要最小の名前表示bindingは担当してよい。Align/Review配置・tooltip・Preview CSSはS06に任せ、共通markupの同時編集を避ける。

## 受入と失敗例

- 実BGCでlivAの通常Align→racM/unique_direct_rbh、review自動表示なし、AlignによるLOSAT job増加なし。Gallery Loadと新規LOSATP Generate直後のResultの両方で確認し、Sessionだけからedgeを取得する実装を防ぐ。
- 同じreferenceの明示Review→racM選択と根拠を表示。racLへの明示変更、Skip、方向選択は使用可能。
- parA基準→racL。representativeのracMを常に選ぶ誤実装を検出する。
- genuine multiple RBH、edge欠落、non-RBH、多段edgeは既存のambiguity/review規則を維持。
- livA Keep→livE KeepとlivA All right→livE Keepで、Keepによる新しいorientation deltaがない。
- crop/reorder/reverse/重複source/別source同名、Save/fresh Load、Generateでidentityを維持。
- failed render/retry、Cancel/stale/superseded、Undo/Redo、位置のみResetと方向を含むResetの契約を維持。
- plan inspectorは生成・再読込・並替後にも人間が識別できる名前。hashを名前にしたfallbackを使わない。

## テストと終了

tests/test_similarity_alignment.py、tests/test_similarity_alignment_web_adapter.py、tests/test_alignment_direction_projection.py、tests/web/similarity-alignment-actions.test.mjs、tests/web/alignment-direction-projection.test.mjsを対象にする。Python変更後はbrowser wheelを再準備し、similarity-alignment-ui.playwright.spec.jsと実GalleryでWorker経路を検証する。

domain proofだけ、stub JS planだけで完了にしない。results/S05.mdにresource由来のedge、adapter、browserの3段階と方向の観測を記録。実際にKeep反転があれば原因まで追い、以前の反転保持と混同しない。

production/test diffを別々にレビューし、SESSION_LOG更新と担当差分commit、英語title/summaryを残す。既存一意RBH契約の実現であり、契約変更を不要に増やさない。
