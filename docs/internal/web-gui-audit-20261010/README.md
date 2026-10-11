# Web GUI リリース前監査 — dev 11a199b6 / 0cb214fa / 8f4cffb7（2026-10-10〜11）

対象は `origin/dev` の `11a199b6`（run 1）、`0cb214fa`（run 2）、`8f4cffb7`（run 3）。release tier の `Tests`（`workflow_dispatch`、`tier=release`）run 38045602906（seed 38045602906）、run 38086450007（seed 38086450007）、run 38110730283（seed 38110730283）の結果をまとめる。run 2 と run 3 は「run 2 の結果」「run 3 の結果」に書く。

**件数**

- run 1（`11a199b6`）: 固有の問題 7 件（P2 1、P3 6）。すべて修正済みか、テストのみで扱いが決まっている（下記）。
- run 2（`0cb214fa`）: 49 ジョブのうち 40 success、4 failure（原因は 3 件）、5 skipped（PR 専用ジョブ）。見つかったのは OV-340（P3）、OV-341（P2）、OV-342（P3）、OV-343（P3）の 4 件。OV-340 と OV-343 は #996（`4f708885`）で修正した。OV-341 は Owner が A（KEEP_ANCHOR_PLAN）を選び（D-46）、OV-342 とともに #997（`8f4cffb7`）で修正した。
- run 3（`8f4cffb7`）: 49 ジョブのうち 44 success、0 failure、5 skipped。新しい問題は 0 件。recipes 40/40、journey 7 本すべて成功、sweep と random walk の finding 0。
- 合計: 固有の問題 11 件（P2 2、P3 9）。すべて修正済みか、テストのみで扱いが決まっている。run 1 と run 2 で修正した項目は run 3 でもすべて直ったままである。

**run 1 の詳細（以下の節）。** run 1 が見つけた固有の問題は 7 件（P2 1、P3 6）。OV-330 は当初 P1（暫定）としたが、原因調査で P2（製品）と P3（テストのみ）に分かれた。製品コードの修正は 3 件（OV-330、OV-334、OV-337）で、いずれも dev にマージ済み。残りはテスト oracle、CI 設定、スクリーンショット手順の修正である。

## 方法と表記

- 2026-09-30 の監査（`docs/internal/web-gui-audit-20260930/`）は領域別に人が dev を操作して再現した。今回は自動監査である。release tier の `Tests` を dev に対して dispatch すると、#974 が追加した `Promotion audit / ...` ジョブが走る。
  - `recipes`: CLI/Python のシナリオを `--check` で実行する（40 件）。
  - `sweep-parity`、`sweep-xss-viewport`、`sweep-gallery-roundtrip`、`sweep-session-load-timing`: sweep。report-only で、失敗してもジョブは落ちない。
  - `journeys-J1`..`journeys-J7`（`tools/audit/user-journey.audit.spec.js`）: ユーザー導線。正しい結果を assert し、失敗はジョブの失敗になる。画面は contact sheet（`promotion-journeys-contact-sheet`）に残る。
  - `random-walk`: seed 付きの乱歩（`tests/web/live-generate-random-walk.promotion.spec.js`）。Generate を挟んだ live 表示と Generate 結果の一致を各 step で比べる。
  - ジョブ名と手順は `docs/internal/WEB_PERIODIC_AUDIT.md` と `tools/audit/README.md` に従う。
- 見つかった項目は S2 が原因を調べた。重大度は暫定のものに「暫定」を付ける。重大度の基準は 2026-09-30 の監査と同じ（P1 誤った出力・編集の消失・有効入力で Generate 不能、P2 機能の不具合・古い状態の残存、P3 軽微）。
- ID は FINDINGS の `OV-nnn`。JS のパスは `gbdraw/web/js/` からの相対パス、それ以外はリポジトリ直下からのパス。行番号は書かない。
- 状態: 修正済み（PR 番号とマージ SHA）、未修正、テストのみ。

## run 1 の結果（ジョブ別）

`Release / gate`: success。

| ジョブ | 結果 | 項目 |
|---|---|---|
| journeys-J1, J2, J5, J6, J7 | 成功 | J7 の phone 幅の Preview が空白に見えた（OV-338） |
| journeys-J3 | 失敗 | OV-333、OV-334。J3/01 の空白は OV-338 |
| journeys-J4 | 失敗 | OV-331 |
| random-walk | 失敗 | OV-330 |
| recipes | 35/40 | OV-336（5 件、手元では bundled runtime で全件成功） |
| sweep-parity | 成功 | Circular HmmtDNA 0 DIFF/32、Linear MJNV 0 DIFF/26 |
| sweep-gallery-roundtrip | 4 件（report-only） | OV-337 |
| sweep-session-load-timing、sweep-xss-viewport | 成功 | なし |

## run 1 で見つかったもの

### P2（1 件）

| ID | 症状 | 原因 |
|---|---|---|
| OV-330(a) | random walk（Linear、Result 1 つ、Auto Reflow off）で、step 2 の `featureVisibility off`（FL2）の後に step 4 の `visibilityRuleAdd {recordId *, CDS, locus_tag ^FL1$, show}` を入れると、live 表示が Generate と異なる。Generate だけに feature `f38ba7c3f` とそのラベルが描かれ、Legend の順も異なる（live は [repeat_region, CDS]、Generate は [CDS, repeat_region]）。Generate が正しい | 古い Feature visibility rule（値は set 済み、action はまだ `off`）から描いた rerender が採用された。rule 準備時の snapshot（`app/rule-matching.js` の `drawingInputReaders`）が feature ごとの visibility 行を持ち、visibility rules を持たなかった。EU で顕在化したが、欠けていた入力はそれ以前からある。**修正済み**: #987（マージ `241dcd4b`）。snapshot に `visibilityRules()` を加えたが、これは Generate 中の rule 編集を古いまま残した（CI shard 2/8 の `gui-audit-20260930-editor`）ため、`editFeatureVisibilityRules` が「自動 rerender 中に rule が変わったら rerender をもう 1 回要求する」方式に置き換えた（`8918a928`）。`tests/web/feature-visibility-actions.test.mjs` に node test |

### P3（6 件）

| ID | 症状 | 原因 |
|---|---|---|
| OV-330(b) | 同じ random walk の step 7（`featureVisibility off` FL1、`visibilityRuleRemove 0` の後）が "did not settle within 180 s; open dialogs: [(page unreadable)]" で終わる。**テストのみ** | walk 自身の `answerDialogs` の microtask ループが Worker の返信を飢えさせ、timer が溜まって renderer が落ちた。`f09b3726`（walk spec）から。EU ではない。**修正済み**: #987（`241dcd4b`、ループが yield する）。製品は変更なし |
| OV-331 | J4 で Legend の全 step を Undo し、全 step を Redo しても Legend の順が戻らない（期待 [transfer RNA, rRNA, CDS, GC skew (+), GC skew (-)]、実際は CDS が末尾）。**テストのみ** | J4 が `__GBDRAW_APP__.sortLegendEntries('desc')` で並べ替えていた。この呼び出しには History の step がなく、記録されるのは Sort ボタンのクリックだけ（`app/legend/sort-actions.js`、`app/history-inputs.js`）。journey が "Sort Z-A" をクリックし、各 step で Undo の深さが増えることを確認するよう修正。製品は変更なし。**修正済み**: #985（`f4cc3575`）。b88ea91c でも同じで、EU ではない |
| OV-333 | J3 の v30（v30-bgc）と v44（v44-two-mode）の Session で、Generate した結果が読み込み直後の preview と異なる（Legend の位置が [0,0] から [408.8,698.8] へ、11 件以上の SVG 差）。現行 writer の Gallery Session は通る。**journey の oracle のみ** | fixture の保存済み Result が保存済み設定から描かれていなかった（v30 は rules と title の前、v44 は staged な record-display 行と main の renderer）。Load は保存済み Result を保持し、Generate は draft から描く。これは `docs/REFERENCE/session-and-request-compatibility.md` の記述どおり。journey を「保存した figure」ではなく records で比べるよう修正。製品は変更なし。**修正済み**: #985（`f4cc3575`）。b88ea91c でも同じ |
| OV-334 | J3 で Save して新しい context に Load すると、全 `config.rules[i].fromFile` が欠落から `false` に変わる（user state が一致しない）。steps 4、9、19（v30-bgc、v44-two-mode、gallery-bgc）。`ui.layoutPreferences.<他モード>` も異なる（未使用スロット） | `applyConfigData`（`services/config.js`）が `fromFile: false` を書いた。Generate は `normalizeSpecificRule` の形（`fromFile` は true のときだけ）を commit する。Load が rule を Generate と同じ形で読むよう修正。**修正済み**: #984（`7aa15f74`、`4994ee83`）。b88ea91c でも同じ |
| OV-336 | `recipes` ジョブで T-CLI-08、T-CLI-10、H-CLI-06、H-CLI-07、H-CLI-08 が "did not resolve the bundled LOSAT runtime" で失敗する（35/40）。**CI ジョブ設定（#974）** | ジョブが `gbdraw setup-losat` で managed runtime を入れる。`resolve_losat_runtime`（`gbdraw/comparisons/losat_runtime.py`）は managed を bundled より優先するが、これらの recipe は `source == bundled` を要求する（`docs/recipes/run_cli_scenarios.py`）。ソース checkout は `gbdraw/bin/linux-x86_64/losat` を持つので setup-losat は不要。probe も任意の source を「利用可」とみなしていた。H0 で bundled runtime を使った手元確認では 5 件とも成功。ジョブから `setup-losat` を外し、probe は `bundled` 以外を利用不可（`skipped-run-locally`）とした。**修正済み**: #985（`f4cc3575`、`3a795b4f`、`80c7a417`） |
| OV-337 | Circular の Gallery Session 4 件（HmmtDNA_basic_circular、HmmtDNA_ATskew、Vnig_TUMSAT-TG-2018、tobacco-chloroplast）で Load → Save すると `webFiles.bindings.linearSeqs[0].depth` が `null` から `[null]` に変わる（`sweep-gallery-roundtrip`、report-only）。当初 P3 暫定 | `reconcileDepthTrackStateAfterSessionFiles`（`services/config.js`）が、どの行にもファイルがなくても Linear の全行を Linear Depth の系列数まで埋めた。2 回目の Load が 1 回目の Save の空行を埋めた。#945（`332a993c`）から、EU より前。Circular と同じく、行にファイルがあるときだけ埋めるよう修正。failing-first は `tests/web/mode-scoped-session.test.mjs`。**修正済み**: #984（`7aa15f74`、`540ca5d0`）。`HmmtDNA_basic_circular` の sweep は成功 |
| OV-338 | J7（390x844）で Generate 後の Result Preview の下が空白、J3/01 で 0.13.0 Session 30 を Load した直後の preview が空白。いずれも SVG は描画されている（J7 の `diagram.png` は正しい）。**未修正（テストの手順のみ）** | viewport の問題で、製品の不具合ではない。J7: Result frame（372x254）がページ y=2092 にあり、撮影時は fold の下。main でも同じ（100% zoom、auto-fit なし）。J3/01: 0.13.0 の保存済み preview は上に約 520 px の余白があり（2777x1360、再生成は 2611x840）、100% の左上が空になる。Fit（19%）で全体が見える。リリース後、OV-328 のテスト PR に含めて、撮影前に Result frame を表示し Fit を押す |

## リリース準備中に見つかったもの

run 1 が見つけたものではない。WS-A、WS-C がリリース準備で見つけ、マージ済み。

| ID | 重大度 | 症状 | 原因・状態 |
|---|---|---|---|
| OV-278 | （記録なし） | Session に Web の draft がなく（CLI、Python、Gallery）、共有 GenBank の 2 番目以降の record を record ID で選び display start を持つと、Generate の request が `startCoordinate: null` になる（Vibrio の record-2 と record-4） | `projectCanonicalSessionRequest`（`services/session-request.js`）が display start を placeholder selector `'#1'` の draft に書き、行の draft は [source, '#n'] で探すため見つからない。Python の decode と保存された Session は正しい。`bindRecordDisplayDrafts`（`services/record-display-model.js`）が、行が分かった時点で draft を record ID の行に移す。**修正済み**: #977（`e10ef579`） |
| OV-325 | P3 | Circular の Session で複数 record を record ID で選び display start を持つと、draft の identity が重複する（両方 ['circular', '#1']）。Load 後の Save が "Duplicate record display draft identity" で失敗する。Linear は影響なし。CLI（`gbdraw circular --records_table` に `record_id` 行 2 つと `display_start`、`--multi_record_canvas --save_session`）がこの形の Session を書くことを確認。dev に以前からある | OV-278 と同じ placeholder。`exportSession`（`services/config.js`）が検証とスライス複製の前に両 drawing の draft を record display の行へ bind する（live の draft には書かない）。`tests/web/session-export-validation.test.mjs`（base で失敗）。**修正済み**: #977（`e10ef579`）。**未修正の範囲**: Gallery 公開（`services/gallery-session-publication.js`）は record 行を読まないので、この形の Circular Session の公開は失敗する。現在この形の Gallery Session はない。必要になったら `projectCanonicalSessionRequest` の placeholder を根本で直す |
| OV-332 | P3 | Linear 比較の Result filters（Minimum identity、Bitscore、E-value、Minimum length）が Generate 時に反映されることを、control も `docs/REFERENCE/web-app.md` の Applies on Generate の行も示していなかった。動作自体は正しい（J6 が確認: Generate 前 77 matches、identity 73 で Generate 後 37） | `gbdraw/web/index.html` の Result filters ブロックに `ui-hint` がなかった。"Applies on Generate: comparison result filters." を追加し、`web-app.md` の行に 4 つの filter を追加。**修正済み**: #977（`e10ef579`）。比較の appearance は時点を確認しておらず含めていない |

## run 1 以外で記録した OV-326..OV-329、OV-335、OV-339

S2 とオーケストレーターがレビュー中に記録した。リリースの条件ではない。

| ID | 重大度 | 症状 | 原因・状態 |
|---|---|---|---|
| OV-326 | P3 | Linear の Result で feature の stroke-width が "2.0" と描かれているとき、stroke を編集せずに Reset all strokes を実行すると、保存済みの Result が再 serialize される（`<defs/>`、`style="opacity: 1;"`、transform が変わる）。見た目の差はない | `resetAllStrokes` が "2.0" と 2 を別の値として比べる。#977 は stroke override があるときだけ呼ぶことで回避。**未修正** |
| OV-327 | P3、テストのみ | CI shard 8/8 の `tests/web/source-recipe-cli-parity.playwright.spec.js` "Circular: each Source recipe probe replays in the CLI to the GUI Result"（step "feature color rule"）が "Execution context was destroyed" で不安定に失敗する（#979、#980） | `setApp` の `page.evaluate` が長い `addSpecificRule()` を await する。Playwright が CDP の "Promise was collected" を navigation のエラーに書き換える。アプリの動作ではない。`evaluateWithRetainedPromise` を使うよう修正。**修正済み**: #982（`ed727105`） |
| OV-328 | P3、テストのみ | OV-327 と同じ形で `addSpecificRule` を `page.evaluate` 内で await する spec が約 20 箇所ある | `tests/web/playwright-long-app-promises.test.mjs` の `LONG_OPERATIONS` に `addSpecificRule` がなく、wrapper の中の body も読まない。run 2 の前に dev の push CI（shard 3、`live-generate-parity:866`）を 7aa15f74 と 4d468bce で落としたため、**#993（77a38dcd）で修正**: 17 箇所を `evaluateWithRetainedPromise` に変え、`LONG_OPERATIONS` に `addSpecificRule` を加えた（D-43） |
| OV-329 | P3、テスト・docs のみ | 現行の Session version 46 と request schema 9 を literal で書いた箇所が 2 つある（`docs/capture/flows/how_to/joint_display_placement.py`、`tests/web/losat-cache-migration.playwright.spec.js`）。今は正しいが、次の writer bump で古くなる | export された定数を import していない。**未修正**: タグ後（#993 には含めていない） |
| OV-335 | P3、潜在 | Legend の順序操作（Sort A-Z/Z-A、Sort by default、Move up/down）は、ボタンのクリックを経たときだけ History の step になる。ユーザーの操作で別経路から呼ぶものは今はない。random walk は直接呼ぶので Sort の Undo を検証していない | `app/history-inputs.js` のクリック adapter だけが記録する。**未修正**: 0.14.0 の後（D-30）。4 つの操作を `app/app-setup.js` で `undoableAction` に包み、walk はクリックにする |
| OV-339 | P3、テストのみ | `tests/web/helpers/live-generate-parity-steps.cjs` で、ラベル編集が先に終わっても 20 ms の `labelOnDialog` poll が page を閉じるまで続く。ハングはしない | `settled` フラグがない（2 行）。**未修正**: タグ後（#993 には含めていない） |

### その他の follow-up（FINDINGS の "Follow-ups noted"）

| ID | 内容・状態 |
|---|---|
| OV-384 | feature popup（z-50）が左上のヘッダボタン（z-40）に重なり、#986 の Load example ボタンがあると Undo が popup の下になる。0.14.0 には含めない。既知の問題として、popup を閉じてから操作する旨をリリースノートに書く |
| OV-401 | 逆向きに alignment された record の stable feature ID が 2 つの frame（biological と drawn）で別の owner により作られ、validator が drawn hash を biological ID とだけ比べていた。c388ccc6（main にあり、どのタグにもなし）からの回帰で 0.13.0 は影響なし。影響は Python の `enrich_svg(..., feature_catalog=)` と Gallery の asset tool のみ。**0.14.0 で修正**（biological ID か同じ行の `drawnSelector.hash` を受け付ける。PR 番号は FINDINGS に記載なし）。重複する builder の除去は identity contract の変更なので、タグの後（D-37） |
| OV-412 | Web の interactive SVG export が、crop または逆向きの record で `data-gbdraw-stable-feature-id` に biological ID でなく drawn hash を書く。**未修正**: タグの後、既知の問題（D-37） |
| OV-413 | Gallery 公開の比較が、CLI replay の書く LOSAT の真偽値セル `True`/`False` を受け付けない（保守用ツールのみ）。「majanivirus PR」で修正済み（PR 番号は FINDINGS に記載なし） |
| CI | リリース後（`test.yml` は governance）: Playwright に `actionTimeout` を設定し、覆われた click が test timeout まで待たないようにする（#986 の shard 4 が 10 分の retry 3 回で 45 分を使った） |

## 問題がなかった範囲

- `sweep-parity`: Circular HmmtDNA 0 DIFF/32、Linear MJNV 0 DIFF/26。`sweep-session-load-timing`、`sweep-xss-viewport` も成功。
- journeys J1、J2、J5、J6、J7 は成功した（J7 は OV-338 の viewport 手順の問題を除く）。`recipes` の 35/40 は成功し、残り 5 件は OV-336 の CI 設定による失敗で、手元の bundled runtime では全件成功した。
- owner graph（`origin/main` `fe6861f0` と dev `11a199b6` の比較、`WEB_PERIODIC_AUDIT.md` の手順 5）: dev で main より高い指標はない。injectionEdges 54→13、forwardClosures 6→0、stateBackdoors 2→0、wholeObjectPorts 43→10、layerImports 77→0。authority PR は不要。
- version と schema の literal: writer 定数は Python/Web とも Session version 46、request schema 9、pyproject は 0.14.0 で一致。「現行」の意味で 44/45/0.13.x を書いた literal はない（古い入力の fixture と reader の test は意図的）。名前のない 2 件が OV-329。
- deploy 専用の検査: `gallery-publication.yml`（dev への push）の run 38045467217 が 11a199b6 に対して `test:web:gallery-publication` を実行し、成功した。`deploy_web.yml` 専用の Gallery artifact 検査はローカルで実行し成功した（`gallery-session-publication.test.mjs` exit 0、`tests/test_gallery_session_semantics.py` と `tests/test_refresh_gallery_sessions.py` 87 passed、`tools/gallery_artifact_manifest.py` exit 0 でツリーはクリーン）。

- run 2（`0cb214fa`）で成功したもの: journeys J1、J2、J3、J4、J5、J7、random-walk（seed 38086450007）、`sweep-parity`（Circular HmmtDNA 0 DIFF/32、Linear MJNV 0 DIFF/26）、`sweep-gallery-roundtrip`（0 件。OV-337 は #984 で解消）、`sweep-session-load-timing`、`sweep-xss-viewport`（ともに 0 件）、recipes 39/40（LOSATP は bundled runtime）。Core 3.10–3.12、slow、Playwright 8 shard、Vibrio full generation、browser、LOSAT cache、packaging も成功。dev への push（`0cb214fa`）の `Tests` 38086337559 と Gallery publication 38086337584 も成功。
- run 2 の journey は J3 と J4 が失敗から成功に変わった。run 1 と同名の step で 1.5 倍を超えて遅くなったものはない。
- H0'（`0cb214fa`）の手元証拠:
  - owner graph（`origin/main` `fe6861f0` と比較）: run 1 と同一。dev で main より高い指標はなく、authority PR は不要（#994 の exception の CB +2 は、この凍結指標には現れない）。
  - version と schema の literal: writer 定数は変わらず Session version 46、request schema 9、pyproject 0.14.0。新しい fixture `tests/fixtures/sessions/hash-color-rule-linear-rc.v30.gbdraw-session.json.gz` は 0.13.0 tag が書いた意図的な旧 reader の positive fixture。新しい test の 44/0.13.0/version=30 は古い入力としての記述だけ。run 1 の 2 件（OV-329）は変化なし、新規の literal はない。
  - deploy 専用の Gallery 検査（`.worktrees/rel-h0b` で `heavy` 経由）: `gallery-session-publication.test.mjs` exit 0（1 pass）、`tests/test_gallery_session_semantics.py` と `tests/test_refresh_gallery_sessions.py` 90 passed（run 1 は 87）、`tools/gallery_artifact_manifest.py` exit 0 でツリーはクリーン。ただし OV-343 のとおり、後者は Python 3.10 では #988 の test が import で失敗する（#996 で修正済み）。
  - capability: run 1 と同じく `packaging`、`full` が journey に未対応（`capabilities-h0b.txt`）。

## 未検証

- どの journey も対応づけていない capability（promotion の body に列挙する。人が見るか Owner が waiver する）。
  - `packaging`: `package-lock.json`、`package.json`、`pyproject.toml`、`tools/verify_gui_offline.py`。ローカルの deploy 検査と release dispatch の packaging ジョブで確認する。
  - `full`: Python モジュール（`gbdraw/analysis/*`、`annotations/*`、`web_support/*`、`cli_utils/*`、`losat_setup.py`、bundled `gbdraw/bin/linux-x86_64/losat`）、`tools/*`、Playwright config、examples。Core pytest、40 件の recipe（LOSATP は bundled runtime）、`sweep-parity` が覆う。この証拠による waiver を勧める。
- `deploy_web.yml` の `npm run test:web:vibrio-generate` は dev では動かない。release tier の `Tests` dispatch の側で走る。
- OV-325 の Gallery 公開側（上記）。
- 画面の目視は contact sheet のレビューだけで、2026-09-30 のような領域別の手動操作は行っていない。

## 証拠と再現方法

- run 38045602906（dev `11a199b6`、2026-10-10）。summary job 114196624187。artifacts は `promotion-audit-*` と `promotion-journeys-contact-sheet`。ログと artifact の展開先は `/home/kawato/gbdraw-baselines/release-0140-20261008/`（`logs/rt1/job-<id>.log`、`evidence/rt1-review.md`、`evidence/rt1-montage`）。
- recipes（`docs/recipes/` の CLI と Python）:
  `PYTHONPATH=$PWD python docs/recipes/run_cli_scenarios.py --all --check`、
  `PYTHONPATH=$PWD python docs/recipes/run_python_scenarios.py --all --check`。
- sweep と journey: `npx playwright test -c tools/audit/playwright.audit.config.cjs <file-or-filter>`（準備は `tools/audit/README.md`）。
- 各項目の調査メモ: `ws-d/OV-330/HANDOFF.md`、`ws-d/OV-331/HANDOFF.md`、`ws-d/OV-333/HANDOFF.md`、`ws-d/OV-338/HANDOFF.md`（`/home/kawato/gbdraw-baselines/release-0140-20261008/` 以下）。全文は `FINDINGS.md`。
- 評価（run 1）: `evidence/owner-graph.md`、`evidence/version-literals.md`、`evidence/deploy-suites.md`、`evidence/capabilities.txt`。
- 評価（run 2、`0cb214fa`）: `evidence/owner-graph-h0b.md`、`evidence/version-literals-h0b.md`、`evidence/deploy-suites-h0b.md`、`evidence/capabilities-h0b.txt`。run 2 のレビューは `evidence/rt2-review.md`、`evidence/rt2-journeys.md`（montage は `evidence/rt2-montage/`）、ログは `logs/rt2/`。

## run 2 の結果

run 38086450007（`Tests`、`workflow_dispatch`、`tier=release`）を dev `0cb214fa`（H0'）に対して実行した。seed 38086450007。run 1（`11a199b6`）からの追加は PR #984–#994 の 11 件（`git log --first-parent --merges 11a199b6..0cb214fa`）。

- #984 Session Load の rule 形（OV-334、OV-337）、#985 promotion audit の first-run 修正（OV-331、OV-333、OV-336）、#987 visibility rule の walk（OV-330）: run 1 の修正。
- #986 結果表示後もヘッダに Load example ボタンを残す、#988 BGC の Gallery 例を og_6 で整列し矢印の軸幅 0.6 に（Web で保存した Session の Gallery 公開も修正）、#989 Web の alignment review dialog の再設計と alignment・popup の不具合修正、#990 CLI Session で逆向き record の比較行を search frame のまま保存（Gallery 公開で既定値の表記を受理）、#991 描画オフの record の Delete settings が編集も削除できるよう R2 を改訂（文書のみ）、#992 複数 record のファイルから描く record を選べるように、#993 Playwright spec の Specific colors rule の promise を保持（OV-328）、#994 feature ID を元の record から作り `hash=` を元の feature に向ける（OV-401、OV-412）。

49 ジョブ: 40 success、4 failure、5 skipped（PR 専用ジョブ）。`Release / gate`: failure（下の 3 原因に従う）。

| ジョブ | 結果 | 項目 |
|---|---|---|
| journeys-J1, J2, J3, J4, J5, J7 | 成功 | J3 と J4 は run 1 の失敗から成功へ（OV-331、OV-333、OV-334 の修正）。J7 の phone 幅は OV-338 のまま。J3/06 の通知は下記 |
| journeys-J6 | 失敗（step 5） | OV-341 |
| random-walk | 成功 | OV-330 の修正（#987）が効いている |
| recipes | 39/40 | OV-340（T-PY-05）。LOSATP は bundled runtime |
| gallery acceptance (Python 3.10) | 失敗 | OV-343 |
| sweep-parity | 成功 | Circular HmmtDNA 0 DIFF/32、Linear MJNV 0 DIFF/26 |
| sweep-gallery-roundtrip | 成功、0 件 | OV-337 の修正後 |
| sweep-session-load-timing、sweep-xss-viewport | 成功、0 件 | なし |
| 上記以外の release ジョブ（Core 3.10–3.12、slow、Playwright 8 shard、Vibrio full generation、browser、LOSAT cache、packaging） | 成功 | なし |

contact sheet（`evidence/rt2-journeys.md`）: 各 journey の所要時間に、run 1 と同名の step で 1.5 倍を超える遅延はない。J7 の合計は +28% で、390 px 幅の Generate step の揺れの範囲。J3/06（main が書いた Session 44 を Load した直後）にサイドバーの "Record rotation, feature placement, or tolerance has changes pending Generate" が出るが、run 1 の J3/06 と同じで、旧 Session の draft が保存済み Result と異なりうること（OV-333 の扱い）と整合する。**項目にしない。**

## run 2 で見つかったもの

### P2（1 件）

| ID | 症状 | 原因・状態 |
|---|---|---|
| OV-341 | J6 の step 5 で、Gallery の BGC Session（BGC0000708–BGC0000713）を Load して Generate し、Minimum identity を 73 に上げて Generate（77 → 37 matches、成功）、さらに編集せず Generate すると `{status: blocked, reason: stale-reference}` になる。canvas に "Showing the last successful result. Current changes were not applied." と "Alignment needs repair: The clicked feature is not a unique member of this Similarity Group."（**Reselect reference** / **Clear alignment**）が出る。以後の Generate は修復するまですべて blocked。Clear か Reselect で回復できる。journey の oracle は正しい | `validateBeforeGenerate`（`app/similarity-alignment.js`）が、plan を画面にある Result の Similarity Group（前の comparison 設定の group）に対して検証する。Minimum identity は plan を消さない（消すのは pair の source/file、program、LOSATP mode だけ。`app/app-setup.js`）。identity 73 では og_6 が {livU, parU} に分かれ、reference の neoU（record-2）が og_15 へ移る。その Generate は reference が group を離れた plan で Result を commit し、次の Generate が `buildHelperRequest` の "not a unique member" で blocked になる。#994 と #989 は無関係。**#988** が fixture の plan を og_1 から og_6/neoU に変えたことで顕在化した（4d468bce でも失敗、run 1 の `11a199b6` は og_1 で成功）。欠陥自体は alignment lifecycle の導入時からあり、Generate 時の comparison 設定（Result filters、grouping 設定）で regroup すれば起きる。重大度 P2（機能が壊れる、回避策あり。描かれた Result が科学的に誤りではない）。**状態: Owner が A を選択（D-46、2026-10-11）**。concern `diagram-generation.similarity-alignment.plan-lifecycle`（PD-OI-028 の拡張）。A KEEP_ANCHOR_PLAN（推奨）、B CLEAR_ON_FILTER_EDIT、C BLOCK_THE_REGROUPING_GENERATE、D CLEAR_AFTER_REGROUPING_GENERATE。A では、regroup しても plan を stale とみなさず、同じ anchor feature で整列を続ける（anchor feature 自体が消えたときは従来どおり修復を求める）。修正 PR は #997。リリースノートの "A comparison change clears it, ..." の一文は Result filters について誤りだったので、PR-DOCS で「Result filters は plan を保つ」と直した。詳細は `ws-d/OV-341/HANDOFF.md` |

### P3（3 件）

| ID | 症状 | 原因・状態 |
|---|---|---|
| OV-340 | `recipes` ジョブの T-PY-05 が "Published artifact is stale: docs/images/t-py-05/python_bgc_losatp_groups.svg (payload differs)" で失敗する（39/40）。docs の成果物のみ | #994 が逆向き record の SVG feature ID を変え、H-CLI-04 と T-CLI-08 の画像は更新したが、この画像を更新しなかった。1679 個の ID 属性のうち 96 個（52 種、すべて record_5 = BGC0000713）が異なり、旧 ID から新 ID に対応づけるとバイト一致する。`f104def3` では `--check` が通る。**修正**: #996（`4f708885`、画像を再生成）。guard は release tier の `recipes` ジョブ（今回これが検出した）。PR の CI は "Recipes standard" だけ |
| OV-342 | OV-341 の経路で、stale-reference の repair 通知が Align 時のメッセージ "The clicked feature is not a unique member of this Similarity Group." を再利用する（誰も feature を click していない）。期待は reference と group を名指す文言（例: "The reference is no longer a unique member of og_6."） | `app/similarity-alignment.js` の `buildHelperRequest` のエラーが `markStaleReference` 経由で出る。**修正**: #997（`8f4cffb7`）。A では regroup の経路自体がなくなり、残る stale-reference 通知は reference か整列対象の feature のどちらが消えたかを示す |
| OV-343 | `gallery acceptance (Python 3.10)` が `tests/test_refresh_gallery_sessions.py::test_refresh_resolves_web_default_colors_over_their_palette` で `ModuleNotFoundError: No module named 'tomllib'`（他の 136 件は成功）。テストのみ | #988（`8e5ce14f`）が追加した test が `tomllib` を素で import する。Python 3.10 に標準では無い。**修正**: #996（`4f708885`、他の test と同じ version 分岐の import）。他の `tomllib` import は tests/ と tools/ のすべてで分岐済み。guard は追加しない（release tier の 3.10 acceptance が検出した。PR の CI は 3.11 のみ） |

状態のまとめ: OV-340 と OV-343 は #996、OV-341（A、D-46）と OV-342 は #997 で修正済み。run 3 で確認した。

## run 3 の結果

run 38110730283（`Tests`、`workflow_dispatch`、`tier=release`）を dev `8f4cffb7` に対して実行した（D-45）。seed 38110730283。run 2（`0cb214fa`）からの追加は #995（`gbdraw/web/CLAUDE.md` の ownership 表 1 行、文書のみ）、#996（OV-340、OV-343）、#997（OV-341、OV-342）の 3 件。

- 49 ジョブのうち 44 success、0 failure、5 skipped（PR 専用ジョブ）。summary job 114388227696。
- journey J1..J7: すべて成功（7 本、失敗 0、206.6 s）。J6 は Minimum identity を上げた Generate のあと、編集なしの Generate も整列したまま成功する。
- recipes 40/40（LOSATP は bundled runtime）。gallery acceptance は Python 3.10 と 3.12 で成功。
- sweep 4 本と random walk（seed 38110730283）: finding 0。`sweep-parity` は circular HmmtDNA で 0 DIFF、0 NO_EFFECT。
- どの journey も対応づけていない capability は run 1、run 2 と同じ `packaging` と `full`（「未検証」の節）。
- dev への push（`8f4cffb7`）の `Tests`、`Gallery publication`、CodeQL も成功。
- contact sheet: run 2 から新しい異常はない。J6 の step 04（Minimum identity 73 で Generate）と 05（編集なしの Generate）に "Alignment needs repair" も "Showing the last successful result" も出ない（run 2 の 05 には両方が出た）。match 数は 77 → 37。viewport 1280×720 では canvas が切れるため（OV-338）ribbon は画面に写らず、整列の確認は match 数と step の成功による。既知の OV-338 と J3/06 の通知は run 1、2 と同じ。run 2 と比べて 1.5 倍を超えて遅くなった step（2 秒超）はない
- 新しい問題: なし。

レビューは `evidence/rt3-review.md`、`evidence/rt3-journeys.md`（montage は `evidence/rt3-montage/`）、ログは `logs/rt3/`。
