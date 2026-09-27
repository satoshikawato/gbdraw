# INSTRUCTION PROMPT S02 — numeric/unit controlsとlifecycleを実装する

## 開始条件と隔離

[共通開始・終了手順](README.md)を実行する。remoteの`fix/issue-619-circular-track-measure-inputs`を取得した専用作業ツリーを使用し、他sessionのcheckout・server・browser・artifactへ干渉しない。同branchのwriterは直列。dirty/unrelated changesをstage/revertしない。

[総合計画](../MASTER_PLAN.md)、AGENTS/CLAUDE/Web CLAUDE、architecture/Product/Web change policiesを読む。scalarの科学的意味、唯一のtyped request/Worker、draft/Result分離、Auto/invaliddraft、History/Session/privacyを守る。実装選択はSOLID/KISS/DRY/YAGNIに照らし、少数owner/pathとsuperseded pathの削除を優先する。

## 目的

Circular Width/Radiusをnumeric text input＋unit selectorで編集できるようにし、Auto/invaliddraft/History/Load/Generateが承認済み仕様で連続するUIを完成させる。

## Ownership / 前提

所有: `app/circular-track-slots/measure-input.js`、既存Circular editor actions、`app.js`/`app/app-setup.js` wiring、`index.html` markup/local CSS、必要最小のHistory/Session接続と対応focused tests。scalar parserを別に作らない。S01のinterface/実測結果を確認する。

S00authorityとS01完了を確認する。componentが必要な表現をS01で実証できていなければ、先にその境界を解決する。native `type=number`でinvalidtextを消すことを修正と呼ばない。

## 作業

1. 2field共通の専用componentを作る。native `type=text inputmode=decimal`、Vue trim/composition、accessible numeric/unit names、slot identityを維持する。
2. 非空modelの数値text/unitはslot scalarから読む。propsを直接mutateせず、owner actionから1回更新する。duplicate acceptedvalue/parallel scalar array/一般unit registryを作らない。
3. approved suffix/percent入力の投影とadapterをS01のcodecに委譲する。unitless数値の意味は現在のselectorに限定する。旧floating suffix/直接object bindingを同時に削除する。
4. unit changeは承認済み規則で行う。推奨Aならnumeric lexemeを保持し、unitだけ変える。自動GenerateやR-basedconversionを追加しない。sameunitとfocus/blurはno-op。
5. approved Auto lifecycleを実装する。推奨Aならlocal preferenceだけを持ち、空欄はnull、unitだけの変更はrequest/Result/History/Sessionへ入れない。manual unitはslotにあるためlocalmirrorで上書きしない。
6. existing History input adapterのtransactionを使用する。二重begin/commitをしない。Undo/Redoでrestoredscalarを読み、localunitsで上書きしない。
7. Save/Load/disabled/inactive/settings-only/CLI-originのactive draft復元を維持する。current-configuration restoreとcanonical committed artifactの責任を逆転しない。
8. 390pxでnumeric/unitが読めるlayoutにする。Auto geometryは単位付きのAuto情報として表示し、selectorのunitの入力値と混同させない。helpにRと次回Generateの意味を記載する。

## 検証・受入

focused component/state checks＋実browserでLoad、read、input、unit change、Auto、IME、Undo/Redoを確認する。C619-01/02/04/05/06/10。branchのwheelを必要時生成する。既存app-lifecycle helperを使用し、表示操作で追加Worker/helper/runがないことをoperationの差分で確認する。

History/Session関連codeを変更した場合は該当Node checksと必要browser roundtripを実行する。C619-07のDOM連携。S01の未変更parserの証拠は再利用する。runtime completionとfullbrowser/download/native replayの最終受入を区別する。

次prompt: S03。

## 完了と公開

所有範囲のproduction/test/docs diffを別々にレビューする。`SESSION_RESULTS/S02.md`にcommands/exits、基準SHA・environment/input・authority、C619結果、未実行check、残るboundary、次sessionの前提を記録する。実測していないgate/browser checkをpassと呼ばない。

**共通終了手順に従い、成果と結果文書を対象branchへcommitし、同名remote branchへpushする。** commit/pushを省略して終了しない。push直前にremote state、branch/upstreamを確認し、force push/main/devへの直接pushをしない。別writerが先へ進んでいればその内容を確認する。authority boundaryが未成立でも独立文書・evidenceはcommit/pushし、依存runtimeは変更しない。自分のlock/processを解放する。

handoffにはremote branchと実際のcommit SHA、English commit title/summary、通ったchecks、未成立の境界、次に使うpromptを含める。
