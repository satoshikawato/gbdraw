# INSTRUCTION PROMPT S00 — authority・表現境界を確定する

## 開始条件と隔離

[共通開始・終了手順](README.md)を実行する。remoteの`fix/issue-619-circular-track-measure-inputs`を取得し、対象branchの既存checkoutを前sessionから引き継ぐ。同branchのwriterは直列。cloneは不要で、worktreeは別作業との隔離が必要な場合だけ使う。他sessionのcheckout・server・browser・artifactへ干渉しない。dirty/unrelated changesをstage/revertしない。

[総合計画](../MASTER_PLAN.md)、AGENTS/CLAUDE/Web CLAUDE、architecture/Product/Web change policiesを読む。scalarの科学的意味、唯一のtyped request/Worker、draft/Result分離、Auto/invaliddraft、History/Session/privacyを守る。実装選択はSOLID/KISS/DRY/YAGNIに照らし、少数owner/pathとsuperseded pathの削除を優先する。

## 目的

Circular Width/Radiusの数値＋単位controlsについて、3つの独立Product outcomeと既存Web scalar draftの表現境界を確認する。runtime実装の前提を成立させる。現行のnumeric/string/typed/%/nullを保持できるか、deterministic disposable checksで証拠を残す。

## Ownership / 前提

所有: 本計画・Decision Packs・S00結果、必要な既存static Product Contractのauthority-only更新、runtime不変のevidence fixtures/tests。runtime、request schema、checker、workflow、Gallery/referenceは変更しない。

3 Packは未署名の候補であり、controlの目標だけを完全なreceiptとみなさない。署名があればその全文に限定して既存static Product Contractへserializeし、machine representationをレビューに提示する。候補と同じruntimeを候補authorityで承認しない。

## 作業

1. 最新origin/devでauthority searchを確認する。BD番号が実在しない場合は引用しない。署名/日付/rationale/preservation/retirement/riskの欠落を推論で埋めない。
2. disposable fixtureでwidth/radiusのtyped number、Web typed numeric text `1.`／`1e-3`、旧bare/px/%、null、invalidtext、0/負/nonfinite/unsupportedunitを確認する。
3. `normalizeCircularTrackSlot`、validation、`buildCircularTrackSlotPayload`、current active-config、Save/Load/native Session reader、History captureの各checkpointでrawtext/value/unitがどう残るか観察する。実 requestはnumber+px/factor/nullのみであることを断言する。
4. 実現可能な表現が既存契約で受理されれば、そのpositive fixtureとoutputをS01へ渡す。受理されなければ、失敗した境界と比較候補をevidence packetに記録する。schemaを黙って変えたりfallbackでacceptanceを弱めない。
5. Product receiptが未署名なら結果に依存runtimeの未成立を記載する。署名済みならauthority-onlyの内容を既存ownerへ正確にserializeする。authority-only公開/mergeに別authorizationが必要なら、具体的な文書/commitを準備してその境界だけを残す。
6. origin/devに選択outcomeのauthorityが入ったか確認する。S01開始条件は署名とtrusted-baseの成立。taskbranchのcandidateだけでは成立しない。

## 検証・受入

基準codec/requestとdraftのfocused Node checks、必要なnative session読み取りを実行。C619-03/04/07/09に対応する境界の証拠を残す。既存mapped contractを変えて唯一の安全根拠にしない。authority-only変更がある場合は既存policy/checkerのschema検証を実行する。

完了は、選択outcomeとcodec表現の受理結果を明示したS00記録。Product未署名・authority未統合・表現未受理があれば、それを未成立と記載した独立成果でありruntime着手の許可ではない。

次prompt: S01（開始条件が成立した場合のみ）。

## 完了と公開

所有範囲のproduction/test/docs diffを別々にレビューする。`SESSION_RESULTS/S00.md`にcommands/exits、基準SHA・environment/input・authority、C619結果、未実行check、残るboundary、次sessionの前提を記録する。実測していないgate/browser checkをpassと呼ばない。

**共通終了手順に従い、成果と結果文書を対象branchへcommitし、同名remote branchへpushする。** commit/pushを省略して終了しない。push直前にremote state、branch/upstreamを確認し、force push/main/devへの直接pushをしない。別writerが先へ進んでいればその内容を確認する。authority boundaryが未成立でも独立文書・evidenceはcommit/pushし、依存runtimeは変更しない。自分が起動したprocessだけを停止し、対象checkoutを次sessionへ引き継ぐ。

handoffにはremote branchと実際のcommit SHA、English commit title/summary、通ったchecks、未成立の境界、次に使うpromptを含める。
