# INSTRUCTION PROMPT — S5: 受け入れ検証、回帰確認、制作手順の確定

文書版: 2.2<br>
改訂日: 2026-09-24

## 1. 任務

`satoshikawato/gbdraw`の紹介動画「Meet gbdraw」の制作工程について、技術的な生成成功、生物学的内容と操作の正しさ、再現可能な運用、完成映像および派生資産の視認性と品質、自動視覚回帰テスト（SSIM）を網羅的に検証してください。確認で見つかった担当範囲の問題は修正し、再実行してください。成果物が存在するだけで合格とは扱わないでください。

「Meet gbdraw」は、初めてgbdrawを知る研究者へ円形図・注釈図・比較図・ラベル選択・全CDS機能別配色・D-loop領域注釈ブラケット・SVG出力を紹介する38秒動画です。
- 完成動画仕様: 1920×1080、30/1 fpsのCFR、1,140フレーム（38.0秒）、H.264/yuv420p、SAR 1:1、英語字幕（SRT焼き込み）、音声ストリームなし。
- 派生資産仕様: README用軽量ループWebP（1〜2MB、15 fps）、学会スライド用5チャプターMP4クリップ。
- 演出仕様: 2秒Afterシーン（Scene 7, 9, 11）における変更箇所のVisual Spotlight（視線誘導）、およびScene 8/10の前進予告字幕によるフリーズ感の完全解消。
- 制作アーキテクチャ: 素材取得はPython Playwright、編集・派生資産生成はPythonとFFmpeg、再編集はブラウザー不要、視覚回帰検証はFFmpeg SSIMフィルター。

## 2. 読むものと依存条件

[総合実装計画書](../01_MASTER_PLAN.md)全体、特に§12のAC-01〜AC-13を読み、現在の`AGENTS.md`、`CLAUDE.md`、適用されるアーキテクチャ・Web指示に従ってください。

`docs/videos/meet-gbdraw/handoffs/S1.md`〜`S4.md`と、その実装・成果物を確認してください。S4までに`build/render/check`（`--visual`含む）、9素材、台本、認定環境、完成MP4、派生資産（WebP, chapters）、再編集テストがあることが前提です。存在しない機能をあるものとして検査を省かないでください。記録上PASSでも実ファイルと簡単な検査で確認してください。

計画参照SHA `4556e04e929a4a85ad28d1833ce7304bd764881c` は現在のHEADや作業ブランチの指定ではありません。先行実装を保持した承認済み作業ツリーを使い、現在のブランチ規則に従ってください。既存差分を捨てたり、main/devへ直接作業したり、無断でcommit/push/PR/merge/公開したりしないでください。

## 3. 対象範囲

受け入れ検証、不足している失敗系テスト、発見した制作コードの不具合修正、派生資産の検証、SSIM視覚回帰の確認、制作README、引き継ぎを担当します。
- **対象外**: 音声追加、比較機能本体の仕様変更、製品の新挙動開発、動画編集用GUIの新設。
- **CI方針**: 初回完成の必須要件ではありません。全PRでの自動フル動画build、必須ゲート化、自動公開は行わないでください。

## 4. 検証手順

### 1. 実装とテストを別々に監査する
制作コードの差分とテストの差分を分けて確認してください。既存fixture、セレクター列、比較条件の二重管理、製品から動画モジュールへの逆依存、rendererからPlaywrightへの依存、無意味な抽象化、使われない試作を探してください。

元のGUIフローやshared helperの意味を変えていないか、既存assertionの緩和で通していないかを確認してください。現在のアーキテクチャ規則が要求する責任・依存の証拠を残し、行数だけで判断しないでください。

### 2. unit/integration/失敗系を実行する
台本のフレーム計算、SRT生成、素材pathの境界、IDとchecksum、必須証拠、フォント、FFmpeg、動画の正規化と連結、派生資産生成を確認してください。少なくとも次の失敗を意図的に作り、非ゼロ終了と診断を確認してください。

- 未知/重複素材ID、欠けた必須素材、改変したPNG/SVG/録画、壊れたmanifestやscene定義。
- `..`やsymlinkで素材ルート外へ出るpath、不在ファイル、既存の非空出力先。
- 不足FFmpeg/ffprobe/codec/字幕フォント、字幕はみ出し、FFmpegエンコード失敗。
- 長すぎる録画、途中までしかデコードできない動画、必要操作の完了を証明しないclip。
- GUI生成timeout、外部通信、download失敗などの取得エラーと、その成功扱い防止。
- 一部CDSだけproduct/geneが切り替わった状態、入力qualifierの改変、意図しないrRNA/tRNAの変更。
- CDSのカテゴリ重複・未分類、ND4/ND4Lの混同、CDSを隠して全件適用と見せる状態、凡例の色/caption不一致。
- D-loopのRegion Annotations TSV不在、Custom Track Slotsの誤設定、原点越えの直線化・分断、P3→P2→P3復元の破綻。
- SSIMスコアが閾値（0.98）を下回る場合の `--visual` 早期検出とエラー報告。

実際の故障を注入できない部分はテストdoubleで境界を検証し、その限界を報告してください。広い例外捕捉で成功を返したり、古い成功素材へフォールバックしたりしてはいけません。

### 3. 独立したフルbuildを2回行う
認定した同じアプリ実体・生成wheel・入力・撮影コード・制作環境を使い、別々の新しい出力領域へ、固定入力からのフル`build`を2回実行してください。2回目に1回目の録画・図・比較キャッシュをコピーしないでください。各GUIフローは新規context等の既存契約を守ります。

```bash
python docs/capture/build_video.py build \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/acceptance-a

python docs/capture/build_video.py build \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/acceptance-b
```

比較するのは、入力identity、生成された図と編集の意味、素材の出所、フレーム数、出力仕様、scene構成、視認性です。raw recordingやMP4のバイト列が完全一致することを要求しないでください。

### 4. 完成メディアと派生資産を実測する
各完成MP4および派生資産について、ffprobeと実デコードで検証してください。

- **完成MP4**: 解像度1920×1080、30/1 CFR、SAR 1:1、H.264、yuv420p、音声なし、実デコード1,140フレームを確認。durationが38.0秒相当であることを確認。
- **README用WebP (`final/meet-gbdraw.webp`)**: ファイルサイズが1〜2MBの範囲に収まり、無限ループ設定が有効であること。
- **学会スライド用チャプタークリップ (`final/chapters/*.mp4`)**: 全5クリップ（`01_`〜`05_`）が存在し、それぞれが独立して正常再生可能であること。

### 5. 素材を変えずに字幕・台本を再編集する
1回目の素材を使い、テスト用の台本コピーで字幕だけを変更して`render`を実行してください。入力素材とmanifestのハッシュを実行前後に計算し、不変であることを証明してください。

ブラウザー起動、ローカルアプリサーバー、LOSAT呼び出しの入口を例外を出すテストdoubleにする等、再編集が取得経路に入らないことを確認してください。

### 6. 既存機能の回帰を確認する
共有ヘルパーを変更した箇所と、その既存呼び出し元を特定してください。対象の既存GUI/図/ダウンロードテストを実行してください。

`docs/images/`、`tests/reference_outputs/`、所有者管理画像（`examples/gbdraw_social_preview.png`）、公開Gallery等について、開始前からある差分と今回発生させた差分を区別してください。意図しない差分があれば原因を直し、通常作業で基準画像を更新しないでください。

### 7. 目視と再生で映像・演出を確認する
最終版を1080pと960×540程度の表示で目視確認してください。
- 字幕の可読性とタイミング
- 完成図の識別、重要な凡例とシンテニー比較リンクの完全性
- 3組のBefore/Afterの差分明瞭性
- 全13 CDSと凡例の4色対応
- 内側トラックのD-loop弧状ブラケットの描画品質（原点越えの接続、他要素との非衝突）
- **Visual Spotlight**: 2秒Afterシーン（Scene 7, 9, 11）における視線誘導効果が自然で、図を隠していないこと
- **フリーズ感の解消**: Scene 7→8、Scene 9→10 の前進予告字幕により、停止感なく次の工程へ引き込まれること
- スライド用チャプタークリップの単体完結性

### 8. 自動視覚回帰テスト（SSIM）の実施
FFmpegの `ssim` フィルターを利用した自動比較を実行します:

```bash
python docs/capture/build_video.py check \
  --out build/videos/meet-gbdraw/acceptance-a \
  --visual \
  --baseline docs/videos/meet-gbdraw/reference-frames/
```

- 全13シーンの代表フレームにおいて、基準フレームとの構造的類似度 `SSIM >= 0.98` を達成していることを確認し、レポートを記録する。

### 9. 制作READMEを完成させる
新規参加者が次を読んで再現できるよう、`docs/videos/meet-gbdraw/README.md`を更新してください。
- 目的と成果物構成（本編MP4、WebP、チャプタークリップ）
- 認定環境の準備とブラウザーwheelの作成手順
- 3つのCLIサブコマンド（`build`, `render`, `check`）および `--visual` オプションの使い方
- 素材を変更する場合と台本・字幕だけを変更する場合の実行経路
- Visual SpotlightとRegion Annotations D-loopの設定方法
- トラブルシューティングとFAQ

## 5. 最終受け入れ条件マトリクス（AC-01〜AC-13）

計画書AC-01〜AC-13に対し、PASS/PARTIAL/BLOCKEDと証拠の場所を対応表で示してください。

| 項目ID | 項目名 | 合格判定基準 |
|---|---|---|
| **AC-01** | 素材取得の意味的完全性 | 4完成図が既存の生物学的・意味的assertionをすべて通過する |
| **AC-02** | 素材の技術的整合性 | 4完成図のSVGおよび1080p PNGが生成され、寸法・チェックサムが記録される |
| **AC-03** | 素材マニフェスト完全性 | `assets.json` が未知ID・改変ファイル・不正パスを早期拒否する |
| **AC-04** | ソースデータの不変性 | P0〜P3の全状態で元GenBankレコードのsequence・annotation情報が完全に不変である |
| **AC-05** | ラベル切替の意味的完全性 | 全13 CDSのproduct→gene切替が実SVGおよび表示文字列で照合される |
| **AC-06** | 機能別配色の意味的完全性 | 13 CDSが呼吸鎖複合体の4群（7/3/2/1）に重複・欠落なく塗り分けられ、凡例と一致する |
| **AC-07** | D-loop領域注釈ブラケットの完全性 | Region AnnotationsパネルおよびCustom Track Slotsによる内側D-loopブラケット（16024..576、原点越え）が描画され、P3→P2→P3復元往復再現性が実証される |
| **AC-08** | 本編動画の技術仕様整合性 | 1920×1080、30/1 CFR、実デコード1,140フレーム（38.0秒）、H.264/yuv420p、音声なし |
| **AC-09** | 字幕同期と視認性 | 生成SRT、焼き込み字幕、台本テキストが一致し、表示枠はみ出しや文字化けがない |
| **AC-10** | 演出・視認性と構図安定性 | 3組のBefore/Afterの構図が安定し、Visual Spotlightによる視覚誘導が適切に機能し、P3と出力録画・Outroの同一性が確認できる |
| **AC-11** | 再編集の高速性と完全分離 | 字幕・台本のみを変更した `render` が、ブラウザー・サーバー・LOSATを起動せずに成功し、素材ハッシュが不変である |
| **AC-12** | 派生資産の完全性 | README用軽量ループWebP（1〜2MB）および学会スライド用5チャプターMP4クリップが自動生成され、サイズ・再生独立性を満たす |
| **AC-13** | 自動視覚回帰テスト | `check --visual` が実行可能であり、基準フレームとの比較で全シーン `SSIM >= 0.98` を達成する |

## 6. 設計原則の監査

- **SRP**: 取得・編集・派生資産生成・検証・公開判断が分離されているか。
- **OCP**: 字幕や場面順、スポットライト座標の変更で取得コードを直す必要がないか。
- **LSP**: 既存ヘルパーや素材の契約を弱めていないか。
- **ISP/DIP**: レンダラーが製品内部状態やPlaywrightへ依存していないか。
- **DRY**: TSV正本、セレクター列、フレーム計算、字幕の正本が一元管理されているか。
- **KISS/YAGNI**: 1本の紹介動画制作に不要な過剰なフレームワークや外部依存を持ち込んでいないか。

## 7. 引き継ぎフォーマット

`docs/videos/meet-gbdraw/handoffs/S5.md`に次を記録してください。

```text
Status: PASS / PARTIAL / BLOCKED
Implementation base: branch, HEAD, dirty state, actual app/wheel identity
Acceptance matrix: AC-01 through AC-13 with evidence paths
Independent builds: two run IDs, commands, source/environment equivalence
Media measurements: codec, dimensions, fps, SAR, decoded frames, audio
Derived assets: WebP size and loop check, 5 chapter clip durations and dimensions
Visual spotlight: cue type, target coordinates, non-destructive appearance
Visual regression (SSIM): command run, min/avg SSIM scores, comparison verdict
Re-edit evidence: unchanged assets and no browser/server/LOSAT execution
Regression: selected tests, rationale, protected paths and existing diffs
Visual review: resolutions, frames, playback actually inspected, remaining review
Documentation: exact newcomer reproduction route
Artifacts: final MP4/SRT/poster/webp/chapter clips/assets/reports paths and hashes
Remaining risks: only observed or explicitly unverified items
Proposed commit title: English
Proposed commit summary: English
```

最後に、完成した範囲、成果物の場所、再実行方法、未解決事項を簡潔に報告してください。実行できなかった検証を合格扱いせず、安全に進められた作業は明確に残してください。動画やフォントのGit追加、リモートへのアップロード、公開は別の明示指示なしに行わないでください。
