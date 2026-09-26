# INSTRUCTION PROMPT — S04 全体受入・最終修正・handoff

あなたはgbdraw Issue #601の統合受入担当です。日本語・中国語の部分font切替によるtext PDF、既知causeと安全な任意Details、Python/JavaScript regex案内が承認された製品動作として全経路で成立することを確認し、必要な修正を完了してください。

## branch取得と必読資料

SESSION_IDは`s04`です。[総合計画書](./MASTER_PLAN.md)の共通手順で**`fix/issue-601-export-output-20260926`**をremoteから取得し、専用cloneでS02/S03の成果を含む最新push済みcommitから作業してください。両正式契約が統合され、PDFとerrorの実装が完了していることを開始条件とします。AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、両Decision Packs、承認の機械表現、全先行sessionのresultとauthority merge状態を読みます。

## 所有範囲と進め方

integration tests、`SESSION_04_RESULT.md`と必要な原因ownerの修正を担当します。先行sessionの実装を尊重し、failureのownerで修正します。authority/guardの緩和、duplicated fallbackの追加、expected reference SVGの更新で受入を通しません。修正でsource/input/environmentが変わった範囲の証拠を更新し、変わらない有効evidenceは再利用します。

## 必須検証

1. merged origin/devの二つのOIPC recordと実装outcomeを照合します。各Must preserve/May retire/riskの全条件を検証表へ写し、option ID一致だけで実現としません。PDFとerrorのconcernを一つの曖昧な合格判定にしません。
2. 日本語・簡体/繁体中国語、Latin/Greek/math混在、normal/bold、既存italic/bolditalic、familyとtspan/textPath/anchor/spacing/rotationを検査します。実download PDF bytesで原文Unicode、font埋込み、画像代用なし、page寸法をassertします。supplementary/結合・異体字は字形と原文が一致するsuccessまたは明示非対応です。
3. 現実的なGallery session/recipeを読み、内容を保った実図を独立PDF rendererで表示します。可読性、glyph位置、label/legend、細線、feature関係を確認し、recipe、artifact digest、render commandsとレビュー結果を保存します。内部smoke fixtureをpublic showcaseに使いません。
4. click時snapshot/filename、図の変更後も旧click内容を保存すること、PDF failure時の同snapshot SVG/PNG recovery、PNG DPI、staging/font/object URL cleanupとretryを検査します。成功Result、Session、Historyを保持します。
5. font asset/license/checksum、sourceとwheel同一配信、same-origin/offline、必要font lazy、初回/反復/最大図のpayload/parse/encode/cache/メモリを検査します。実offline auditにはbrowser-offline-qa skillを適用し、専用cloneのwheelを準備します。
6. known/unknown/prepare/regex/identity/cleanupをsource→Worker→UIで検査します。private sentinel非露出、safe paramsと表示上限、summary/action、任意Details、実stage、主因の保持、cancel/stale、rollbackとretryを確認します。native CLI/APIと既存validationの修正情報を保持します。
7. Python-only Color/Labelのnative一致、atomic commit、TSV/preset/Session、History/staleとJS検索のsyntax案内を検査します。desktop/390 px/keyboard、Details展開と実際に使えるactionを確認します。
8. 新しいblocking assertionsがPR inventoryにあり、local-only Playwright成功だけでhard safetyを主張していないことを確認します。trusted-base Web policy、architecture contracts、関連Node/Python/browser、ruff、参照SVG comparisonを実行します。関連package/offline CIとdev stagingの実行境界を明示します。
9. production、tests、docs、generated差分を別々にreviewします。font選択、文言、transport、regex各ownerが一つで、置き換えた経路の削除、OE/PE/CB非増加、不要schema/dependency/compatibilityなしを確認します。

## 完了条件

総合計画書の必須受入表と各Decision PackのAND条件が合格し、必要な修正と関連再検証が終わった時点で実装完了とします。実行していないCI/staging/deployを成功扱いにせず、未実行と実装完了・公開状態を分けて報告してください。required local evidenceの未合格は完了扱いにしません。

## 終了時のcommit・push

`SESSION_04_RESULT.md`へ検証したhead/base、acceptance表、commands/results、artifact/recipe、owner/path review、実装状態、CI/公開状態と実際の未確認事項を保存します。最後の記録だけの変更でも、総合計画書の共通終了手順で**commitし、同名remote branchへpushしてください**。
Commit title例: `Verify Issue 601 export and error recovery behavior`。
remote SHA一致を確認し、English commit title/summary、branch、commit SHA、結果fileを最終handoffします。PR作成・merge・deployは各操作の明示的な許可に従います。
