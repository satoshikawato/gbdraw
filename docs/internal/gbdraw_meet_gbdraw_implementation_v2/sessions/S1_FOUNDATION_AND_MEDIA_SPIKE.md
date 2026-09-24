# INSTRUCTION PROMPT — S1: 制作基盤と最小映像経路の実証

文書版: 2.2<br>
改訂日: 2026-09-24

## 1. セッションの目的とゴール

`satoshikawato/gbdraw`に、紹介動画「Meet gbdraw」の自動制作パイプラインを構築します。
本セッションのゴールは、**制作環境を確定し、実際の図（ヒトミトゲノム円形図）からテストMP4・WebPを生成し、短いGUI操作録画、字幕描画、およびSSIM測定機能が成立することを実証すること**です。

完成動画は、1920×1080、30 fps、38秒（1,140フレーム）、英語字幕、無音のWeb紹介動画です。ヒトミトゲノム、タバコ葉緑体、Lambda–DE3比較、5件のBGC比較の完成図を紹介し、product→geneラベル切替、全CDS機能別配色、D-loop領域注釈ブラケット追加のBefore/AfterとSVG出力の実録画で構成されます。また、GitHub README用の軽量ループアニメーション（`meet-gbdraw.webp`）やスライド用クリップ（`chapters/*.mp4`）の出力経路も基礎検証します。

本セッションでは、設計説明にとどまらず、担当範囲のコード、テスト、実行証拠を確実に残してください。

## 2. 前提資料と参照基準

作業前に以下のリポジトリ内資料を確認してください。

- [総合実装計画書](../01_MASTER_PLAN.md)（文書版2.2）
- `AGENTS.md`、`CLAUDE.md`（Web関連に触れる場合は`gbdraw/web/CLAUDE.md`）
- `docs/capture/README.md`、`docs/capture/config.py`、`docs/capture/web_server.py`、`docs/capture/flows/web_capture.py`
- `docs/capture/flows/human_circular.py`、`docs/capture/flows/tutorials/gui_first_circular.py`
- `gbdraw/web/tutorial-data/manifest.json`、ブラウザーwheel準備手順

> [!IMPORTANT]
> 計画書の参照SHA（`4556e04e929a4a85ad28d1833ce7304bd764881c`）は資料の版であり、作業ブランチの指定ではありません。`AGENTS.md`の規則に従い、作業ブランチは最新の`origin/dev`から派生させ、`main`や`dev`への直接コミット・プッシュは行わないでください。

## 3. 入出力仕様とデータ契約

### 入力データ
- 固定ヒトミトゲノム: `gbdraw/web/tutorial-data/human-mitochondrion/HmmtDNA.gbk`（NC_012920.1）

### 確定・生成すべき成果物
- `docs/videos/meet-gbdraw/storyboard.json`: 総合計画書§4.2に基づく13場面、新字幕、合計1,140フレームの定義
- `docs/videos/meet-gbdraw/environment.json`: 実測した制作環境（OS、Python、Playwright、Chromium、FFmpeg、フォント等）の固定値
- `docs/capture/video/model.py`: 台本・素材契約を表現・検証する軽量データモデル
- テスト映像成果物（作業runディレクトリ内）:
  - ヒトミトゲノムの実ダウンロードSVG
  - 映像用高解像度PNG
  - 数秒（例: 3秒・90フレーム）のテストMP4（指定フォントによる英語字幕焼き込み）
  - 短いGUI操作の録画ファイル（1920×1080）
  - テストWebPアニメーション
  - SSIMフィルター動作確認レポート

## 4. 推奨実装手順

### Step 1: 実装ベースの確認と環境測定
1. 現在のブランチ、HEAD、作業ツリー差分を記録する。
2. 必要に応じて `python tools/prepare_browser_wheel.py` を実行し、生成ブラウザーwheelのハッシュを記録する。
3. Python、Pillow、Playwright、Chromium、FFmpeg、ffprobe、利用コーデック（libx264）、libass/subtitlesフィルター、およびSSIMフィルターの動作を確認する。
4. **フォントの選定**: GUI用のWOFF2とは別に、字幕・表紙用として制作環境にあるTTF/OTF書体を1つ選定し、書体名・ハッシュ・利用条件を記録する（任意のシステムフォントへの無言フォールバックは禁止）。

### Step 2: 台本モデルと台本JSONの作成
1. `docs/videos/meet-gbdraw/storyboard.json` に、総合計画書§4.2の13場面（合計1,140フレーム）と更新された字幕テキストを定義する:
   - `intro`: "Meet gbdraw — Publication-ready genome figures"
   - `circular`: "Draw circular genomes effortlessly."
   - `plastome`: "Annotate plastomes with inverted repeats."
   - `genome-comparison`: "Align and compare whole genomes."
   - `cluster-comparison`: "Visualize biosynthetic gene clusters."
   - `labels-product`: "Too crowded? Full product labels."
   - `labels-gene`: "Switch to clean, readable gene symbols."
   - `colors-before`: "Categorize features by biological role."
   - `colors-functional`: "Color CDS by respiratory complex."
   - `dloop-before`: "Add custom region annotations."
   - `dloop-bracket`: "Origin-spanning D-loop bracket."
   - `export-svg`: "Export crisp, scalable vector SVG."
   - `outro`: "Runs in your browser. Try gbdraw.app"
2. `docs/capture/video/model.py` に、台本構造のバリデーション関数を実装する。scene順序から開始フレームを累積計算し、合計フレーム数が1,140であることを確認する。未取得素材をダミー成功で偽装しない。

### Step 3: ヒトミトゲノム円形図のSVGダウンロードとPNG化
1. 既存の `human_circular.py` / `gui_first_circular.py` のフローを利用し、実UIから `NC_012920.1` のSVGをビルド作業領域へダウンロードする。既存ドキュメント画像（`docs/images/`）やreference outputsを上書きしない。
2. ダウンロードした静的SVGを、固定Chromiumと既存Webフォントを用いて映像用PNG（アスペクト比維持、必要要素を包含）へ変換する。`document.fonts.ready` 等でフォント読み込み完了を待機する。

### Step 4: テストMP4/WebPの生成と字幕・SSIM検証
1. 変換したPNGから、例えば3秒（90フレーム）のテストMP4および短尺テストWebPを生成する。
2. 選定したフォントを用い、FFmpeg（subtitlesフィルター等）により英語字幕が美しく描画されることを確認する。
3. FFmpegの `ssim` フィルターを用いて、生成フレーム間の類似度測定が正しく実行できることを確認する。
4. 出力MP4をffprobeおよび完全デコードで検査し、解像度（1920×1080）、CFR 30fps、H.264、yuv420p、音声なしを確認する。

### Step 5: Playwright実画面録画の実証
1. 実アプリの短いUI操作の録画を実証する。
2. **録画方式**: Playwright公式の `browser.new_context(record_video_dir=..., record_video_size={"width": 1920, "height": 1080})` を推奨する。またはCDP Screencast APIを用いてもよい。
3. 録画サイズが1920×1080であること、録画終了時に正しくファイルが保存されることを実測する。
4. 録画停止とコンテキスト破棄を `finally` 節で確実に実行する。

## 5. 既知の落とし穴とガードレール

- **フォントの混同**: Web UIが読み込むWOFF2フォントをFFmpeg/Pillowがそのまま読めると仮定しない。編集用TTF/OTFを環境から明示指定する。
- **既存資産の破壊防止**: `docs/images/` や `tests/reference_outputs/`、所有者管理のプレビュー画像に差分を出さない。全出力先をビルド領域に閉じる。
- **偽装の禁止**: 台本の未取得素材をダミーオブジェクトで満たして合格扱いにしない。
- **過剰実装の排除**: 汎用動画編集エンジン、プラグイン機構、別言語Playwrightラッパーなどは作成せず、最小限の関数とdataclassで実装する（KISS/YAGNI）。

## 6. 受け入れ条件（Acceptance Criteria）

1. 実アプリから得たヒトミトゲノムSVGが意味的検証を通り、映像用PNGとテストMP4へ変換できる。
2. テストMP4を実デコードでき、指定寸法（1920×1080）、30fps CFR、音声なしを確認できる。
3. Playwrightによる実画面録画が成立し、実サイズとファイル保存が確認できる。
4. 字幕書体が意図通り描画され、フォント不在・FFmpeg不足のエラーハンドリングが明確である。
5. `storyboard.json` が構造的に正しく合計1,140フレームであり、未完成素材が正しく未完成として検出される。
6. 単体テストで台本の重複ID、不正フレーム数、不正パス境界が検証されている。
7. FFmpegによるSSIM測定コマンドが正しく機能することが実証されている。
8. 既存の承認済み画像・製品挙動に変更がない。

## 7. 引き継ぎフォーマット

セッション完了後、`docs/videos/meet-gbdraw/handoffs/S1.md` に以下を実値で記録してください。

```text
Status: PASS / PARTIAL / BLOCKED
Implementation base: branch, HEAD, upstream, dirty state, preserved changes
Completed work: implemented paths and responsibility boundaries
Environment: actual versions, wheel/font hashes, recorder/container verification
Tests: exact commands, exit status, observed result, skipped/unverified items
Artifacts: paths and hashes of SVG, PNG, trial MP4, raw recording, reports
Visual review: actual frames/playback inspected and findings
Protected files: before/after differences and relevant regression results
Decisions: contracts and file locations established for S2
Remaining work: blockers and safe next steps
Proposed commit title: English
Proposed commit summary: English
```
