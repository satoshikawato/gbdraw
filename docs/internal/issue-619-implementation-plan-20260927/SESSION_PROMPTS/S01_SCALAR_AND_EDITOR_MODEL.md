# INSTRUCTION PROMPT S01 — scalar ownerとeditor codecを実装する

## 開始条件と隔離

[共通開始・終了手順](README.md)を実行する。remoteの`fix/issue-619-circular-track-measure-inputs`を取得し、対象branchの既存checkoutを前sessionから引き継ぐ。同branchのwriterは直列。cloneは不要で、worktreeは別作業との隔離が必要な場合だけ使う。他sessionのcheckout・server・browser・artifactへ干渉しない。dirty/unrelated changesをstage/revertしない。

[総合計画](../MASTER_PLAN.md)、AGENTS/CLAUDE/Web CLAUDE、architecture/Product/Web change policiesを読む。scalarの科学的意味、唯一のtyped request/Worker、draft/Result分離、Auto/invaliddraft、History/Session/privacyを守る。実装選択はSOLID/KISS/DRY/YAGNIに照らし、少数owner/pathとsuperseded pathの削除を優先する。

## 目的

Circular measureの解釈を一つにし、保存scalarを数値text＋unitへ投影できるDOM-free editor codecを実装する。px/factor/%/Auto/invaliddraftの意味を保ち、次のUI sessionにlosslessなread/write boundaryを渡す。

## Ownership / 前提

所有: `app/track-slot-validation.js`のCircular scalar boundary、`app/circular-track-slots.js`のscalar payload接続、`app/circular-track-slots/measure-editor.js`、必要な`services/session-request.js` projection、対応Node tests。HTML/CSS/component/History orchestration/公開Galleryは所有しない。

S00結果の表現positive fixture、署名済み3outcomes、origin/devにあるactive authorityを確認する。署名待ちのAを実装しない。推奨以外が選ばれていればmasterの該当契約を更新してから作業する。

settings-only の前提修復と実受理証拠は
[補助 session 結果](../SESSION_RESULTS/S00_SETTINGS_SAVE_REPAIR.md)を確認する。
同結果の inactive 跨mode config 差と nonfinite History の限界も読み、
settings-only 開始条件の成立を全 runtime 受入の完了と混同しない。
この補助 session は scalar owner/editor codec を実装していない。

## 作業

1. validationとcircularScalarPayloadの重複scalar解釈を既存validation ownerへ収束させる。pair/nullを返し、payloadとrow validationが同じownerを呼ぶ。superseded parserを同じ変更で削除する。
2. typed objects、current/legacy strings、numerictext、unitによるread/write codecを実装する。Load/readはrawscalarを変更しない。typedrequestにrawtextや%unitを入れない。
3. approved legacy suffix入力adapterを同codec/ownerへ接続する。selectorで指定したplain numberのunitを勝手にfactorへ戻さない。pixelとfactorのphysical conversionは作らない。
4. invalid/partialtextを保持し、Auto/nullと区別する。zero/negative/nonfinite/unsupportedunitをautoに修復しない。displayprecisionを編集値へ適用しない。
5. session-requestのscalar projectionを共有read helperへ委譲する。width/radiusとpurepixelgap/legacyspacingの別意味を混同せず、gap zeroなど既存callerの受理を保つ。旧独立formatterを残さない。
6. typed/string/%→view→draft→canonical pairのexact equalityとno-read-mutationをNode testsで確認する。新fileは`tests/web/circular-track-measure-editor.test.mjs`。rendererやschemaを変更しない。

## 検証・受入

新codec tests、track-slot-validation、circular-track-slots、session-request、active-config/draft-authorityのうち変更経路を保護するfocused checksを実行する。C619-03/04/07/09のDOM-free部分。architecture-contractsと通常policy gateを実行し、owner/path evidenceを残す。

S02へ公開するinterfaceはread-only view、numericupdate、unitchange、legacyinputadapter、error/Autoの扱い。componentが別parserや第二draftを作る必要がない形にする。実際のUI/History/browser受入はこのsessionのpassとして報告しない。

次prompt: S02。

## 完了と公開

所有範囲のproduction/test/docs diffを別々にレビューする。`SESSION_RESULTS/S01.md`にcommands/exits、基準SHA・environment/input・authority、C619結果、未実行check、残るboundary、次sessionの前提を記録する。実測していないgate/browser checkをpassと呼ばない。

**共通終了手順に従い、成果と結果文書を対象branchへcommitし、同名remote branchへpushする。** commit/pushを省略して終了しない。push直前にremote state、branch/upstreamを確認し、force push/main/devへの直接pushをしない。別writerが先へ進んでいればその内容を確認する。authority boundaryが未成立でも独立文書・evidenceはcommit/pushし、依存runtimeは変更しない。自分が起動したprocessだけを停止し、対象checkoutを次sessionへ引き継ぐ。

handoffにはremote branchと実際のcommit SHA、English commit title/summary、通ったchecks、未成立の境界、次に使うpromptを含める。
