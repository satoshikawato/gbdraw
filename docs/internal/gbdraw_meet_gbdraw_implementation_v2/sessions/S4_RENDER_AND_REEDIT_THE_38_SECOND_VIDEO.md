# INSTRUCTION PROMPT — S4: 38秒の編集・派生資産生成・ブラウザー不要の再編集を実装する

文書版: 2.2<br>
改訂日: 2026-09-24

## 1. セッションの目的とゴール

`satoshikawato/gbdraw`の紹介動画「Meet gbdraw」について、S2/S3で取得した素材群と台本から、**38秒の完成MP4、英語字幕、表紙画像（poster.png）、README用ループWebP、学会スライド用チャプタークリップ、およびレビュー用フレーム**を生成する編集パイプラインを実装します。

本セッションのゴールは以下の通りです:
1. **完成動画の生成**: 全9素材・13 sceneから38秒（1,140フレーム、30 fps）のMP4と焼き込み字幕を生成する。
2. **Visual Spotlightの合成**: 2秒間のAfterシーン（Scene 7, 9, 11）において、変更箇所（geneラベル、機能別配色、D-loopブラケット）へ視聴者の視線を即座に誘導する控えめな視覚的合図（スポットライト）を合成する。
3. **派生資産の自動生成**:
   - `final/meet-gbdraw.webp`: READMEトップ配置用の軽量（1〜2MB）ループアニメーション。
   - `final/chapters/`: 学会発表スライド用の5つの独立チャプターMP4クリップ。
4. **高速な再編集の実証**: 素材を取得し直さず、台本（字幕・尺・スポットライト設定）の変更だけで動画や派生資産を再生成する工程（`render`）を確立する。再編集時にブラウザー、サーバー、LOSATを絶対に起動しない。
5. **CLIの統合**: `build`（全工程実行）、`render`（再編集のみ）、`check`（整合性検証、`--visual`によるSSIM回帰検証を含む）の3コマンドを完成させる。

## 2. 前提条件と依存関係

- [総合実装計画書](../01_MASTER_PLAN.md)（文書版2.2、特に§3.1、§4、§9、§10、§12）
- `AGENTS.md`、`CLAUDE.md`
- 先行セッションS1〜S3の実装、成果物、および引き継ぎ記録（`S1.md`, `S2.md`, `S3.md`）
  - 全9素材が `assets.json` に揃っていることを確認してください（特に `human.dloop-bracket`）。
- 制作環境: S1で認定されたFFmpeg、ffprobe、libass、選定された編集用TTF/OTFフォント

> [!IMPORTANT]
> S1〜S3の成果物を保持した作業ツリーから継続してください。`main`や`dev`への直接コミット・プッシュは行わず、無断のresetや古い参照コミットへの強制チェックアウトは避けてください。

## 3. 編集仕様と絵コンテ契約

完成動画仕様:
- 解像度: 1920×1080
- フレームレート: CFR 30/1 fps
- 総フレーム数: 1,140フレーム（38.0秒）
- コーデック: H.264 / yuv420p、音声ストリームなし、Web再生向けfaststart
- 字幕: 台本から生成されたUTF-8 SRTに基づく焼き込み字幕

### 13 Sceneの台本規範値（storyboard.json）

| Scene ID | フレーム数 | レイアウト | 使用素材ID | 英語字幕 | 演出・意図 | 視覚誘導 (Spotlight) |
|---|---:|---|---|---|---|---|
| `intro` | 60 | `grid2x2` | 4完成図（S2） | Meet gbdraw — Publication-ready genome figures | ツール名と「論文品質」を一瞬で認知 | なし |
| `circular` | 90 | `single` | `human.circular` | Draw circular genomes effortlessly. | 円形ゲノムを簡単に描画 | なし |
| `plastome` | 90 | `single` | `tobacco.plastome` | Annotate plastomes with inverted repeats. | 葉緑体と逆位反復の注釈 | なし |
| `genome-comparison` | 150 | `single` | `lambda-de3.comparison` | Align and compare whole genomes. | 全ゲノム比較・アライメント | なし |
| `cluster-comparison` | 150 | `single` | `bgc.comparison` | Visualize biosynthetic gene clusters. | BGC遺伝子クラスター比較 | なし |
| `labels-product` | 60 | `single` | `human.labels-product` (P0) | Too crowded? Full product labels. | Before 1: 文字が多く混雑した状態 | なし |
| `labels-gene` | 60 | `single` | `human.labels-gene` (P1) | Switch to clean, readable gene symbols. | After 1: 簡潔で可読性の高いgene名へ | **geneラベル領域の強調枠/グロー** |
| `colors-before` | 60 | `single` | `human.labels-gene` (P1) | Categorize features by biological role. | Before 2: 単一色。次の配色への前進予告 | なし（字幕でアクション予告） |
| `colors-functional` | 60 | `single` | `human.functional-colors` (P2) | Color CDS by respiratory complex. | After 2: 呼吸鎖複合体ごとの色分け | **着色CDS群のハイライト** |
| `dloop-before` | 60 | `single` | `human.functional-colors` (P2) | Add custom region annotations. | Before 3: 領域注釈追加への前進予告 | なし（字幕でアクション予告） |
| `dloop-bracket` | 60 | `single` | `human.dloop-bracket` (P3) | Highlight origin-spanning D-loop. | After 3: 原点越えD-loopブラケット表示 | **内側D-loopブラケットの強調** |
| `export-svg` | 120 | `single` | `human.svg-export` | Export crisp, scalable vector SVG. | ベクターSVG保存の実画面録画 | ダウンロードUIへの視線誘導 |
| `outro` | 120 | `single` | `human.dloop-bracket` (P3) | Runs in your browser. Try gbdraw.app | インストール不要・Web完結の最強メリット | なし |

> [!TIP]
> **フリーズ誤認の解消とVisual Spotlight**:
> 1. P1（Scene 7→8）およびP2（Scene 9→10）が連続する区間では、Scene 8（"Categorize features by biological role."）およびScene 10（"Add custom region annotations."）の字幕が次の編集工程への前進予告として機能します。
> 2. 2秒間のAfterシーン（Scene 7, 9, 11）では、変更された箇所（ラベル、配色、ブラケット）に半透明のパルスリングや薄いハイライト枠を合成することで、視聴者が2.0秒以内に迷わず変化を認知できるようにします。図自体を覆い隠さない控えめな演出としてください。

## 4. 推奨実装手順

### Step 1: 入力契約の検証（早期チェック）
1. `render.py` は、素材ファイル（`assets.json`）と台本（`storyboard.json`）のみを入力として受け取る。
2. 未知素材ID、改変ファイル、不正な相対パス（`..`等）、フレーム数不一致を検出した場合は編集開始前にエラーとする。
3. `render` の実行パスにおいて、Playwrightやブラウザー、ローカルWebサーバー、LOSATを絶対に起動しないよう依存境界を分離する（遅延インポート等を活用）。

### Step 2: 静止画シーンと表紙の合成
1. 素材のアスペクト比を維持し、字幕用の安全領域（下部マージン等）を確保した1920×1080キャンバスへ配置する。
2. 表紙（`intro`）の2×2グリッドは、S2の4完成図から自動合成する。
3. `poster.png` を同じ合成ロジックから生成する（既存の `examples/gbdraw_social_preview.png` を上書きしない）。
4. P0〜P3のBefore/Afterシーンは、S3で共通構図（同一中心・同一半径）に調整されたPNGを一様に配置する。
5. **Visual Spotlightの適用**: Scene 7, 9, 11 のフレームに対し、台本で指定された注目領域座標（または差分領域マスク）に基づき、FFmpegのオーバーレイや描画フィルターで薄いハイライト合図を合成する。

### Step 3: 字幕の生成と焼き込み
1. `storyboard.json` からシーン境界の整数フレーム累積和を計算し、UTF-8のSRTファイルを生成する。
2. S1で選定・固定したTTF/OTFフォントを用い、FFmpegの `subtitles` フィルター（libass）により字幕を正確に焼き込む。
3. 日本語パス、空白、引用符、コロン等を安全にエスケープし、コマンドインジェクションを防止する。

### Step 4: 動画素材の正規化とタイムライン連結
1. 各シーンを共通の仕様（1920×1080、30/1 fps、SAR 1:1、yuv420p、開始PTS=0）へ正規化して連結する。
2. 静止画シーンは指定フレーム数分だけ表示する。
3. 出力録画シーン（`export-svg`）は、実操作が120フレームより短い場合は末尾フレームを延長し、微細な超過時は末尾トリミングまたは微細なPTS正規化を行って正確に120フレームへ揃える。
4. シーン接続はシンプルなカット編集とし、過度なフェードやトランジションを追加しない。

### Step 5: 派生資産の自動生成（Derived Assets）
1. **README用ループWebP (`final/meet-gbdraw.webp`)**:
   - 完成MP4から、代表シーン（イントロ、円形図、比較図、機能別配色、SVG出力）を抽出し、640×360または960×540、15 fps、無限ループのWebPを生成する。
   - ファイルサイズが1〜2MBに収まるよう品質パラメータを調整する。
2. **学会スライド用チャプタークリップ (`final/chapters/`)**:
   - `01_introduction.mp4` (Scene 1: 4完成図グリッド、2.0秒)
   - `02_single_genome.mp4` (Scene 2〜3: ヒトミトコンドリア & 葉緑体、6.0秒)
   - `03_synteny_comparison.mp4` (Scene 4〜5: Lambda–DE3 & BGC比較、10.0秒)
   - `04_biological_styling.mp4` (Scene 6〜11: ラベル・配色・D-loopのBefore/After、12.0秒)
   - `05_export_and_web.mp4` (Scene 12〜13: SVGダウンロード & Web完結、8.0秒)
   - 各クリップは1920×1080、H.264、完全なスタンドアローン動画としてスライドへそのまま貼り付け可能とする。

### Step 6: CLIコマンドの統合
`docs/capture/build_video.py` に以下のサブコマンドを実装する:

```bash
# 1. 全工程の実行（素材取得＋編集＋派生資産生成＋検証）
python docs/capture/build_video.py build \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/run-001

# 2. 既存素材からの高速な再編集（ブラウザー不要、台本・字幕・派生資産更新）
python docs/capture/build_video.py render \
  --assets build/videos/meet-gbdraw/run-001/assets.json \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/render-002

# 3. 完成成果物の検証（仕様整合性）
python docs/capture/build_video.py check \
  --out build/videos/meet-gbdraw/run-001

# 4. SSIMによる画質・視覚的回帰テスト（オプション）
python docs/capture/build_video.py check \
  --out build/videos/meet-gbdraw/run-001 \
  --visual \
  --baseline docs/videos/meet-gbdraw/reference-frames/
```

- 出力先には新規または空のディレクトリを要求し、既存の成功runを無断で上書きしない。
- 失敗した場合は非ゼロ終了とし、過去の成功成果物で誤魔化さない。

### Step 7: レビュー用フレームの抽出とSSIM回帰検証
1. 完成MP4から、各シーンの代表フレーム、シーン接続境界の直前・直後フレーム、Before/Afterの切り替わりフレームを画像として抽出し、`reports/review-frames/` に保存する。
2. `--visual` 指定時は、FFmpegの `ssim` フィルターを用いて基準フレームとの類似度を計算し、全シーンで `SSIM >= 0.98` であることを自動検証してレポートに出力する。

## 5. 既知の落とし穴とガードレール

- **責任境界の破壊**: `render` 処理に Playwright の `Page` オブジェクトや gbdraw の内部状態を渡さない。CLIで `render` を指定しただけで撮影モジュールがロードされる構造を避ける。
- **素材の改変禁止**: 字幕の修正やシーン順序の変更によって、原素材（SVG、PNG、raw動画）のハッシュが変化してはいけない。
- **過剰な視覚演出の禁止**: Visual Spotlightはあくまで注目領域への視線誘導であり、ゲノム図自体の線やラベル、色を損なうような激しいアニメーションや派手なエフェクトを用いない。
- **シェルインジェクションの防止**: FFmpegの呼び出しには `subprocess.run(list_args, shell=False)` を使用し、シェル展開に依存しない。

## 6. 受け入れ条件（Acceptance Criteria）

1. **AC-08準拠**: 全13 sceneが連続し、合計1,140フレーム（38.0秒）。完成MP4を実デコードして完全一致を確認できる。
2. **AC-09準拠**: 生成されたSRT、焼き込み字幕、台本テキストが一致し、字幕の重複やはみ出しがない。
3. **AC-10準拠**: 3組のBefore/Afterの構図が安定しており、Visual Spotlightによる視覚誘導が適切に機能し、P3と出力録画・Outroの同一性が確認できる。
4. **AC-11準拠**: 字幕だけを変更した `render` コマンドが、ブラウザー・サーバー・LOSATを起動せずに成功し、素材のハッシュが不変である。
5. **AC-12準拠**: 派生資産（`final/meet-gbdraw.webp`、`final/chapters/*.mp4`）が自動生成され、サイズ要件と再生独立性を満たす。
6. **AC-13準拠**: `check --visual` が実行可能であり、基準フレームとの比較で `SSIM >= 0.98` を達成する。
7. `build`, `render`, `check` の各CLIが仕様通りに動作する。
8. 抽出されたレビュー用フレームを目視確認し、映像の欠陥がない。

## 7. 引き継ぎフォーマット

セッション完了後、`docs/videos/meet-gbdraw/handoffs/S4.md` に以下を実値で記録してください。

```text
Status: PASS / PARTIAL / BLOCKED
Implementation base: branch, HEAD, S1-S3 work confirmed
CLI: actual build/render/check commands and output rules
Implementation: renderer responsibilities, model changes, dependency boundaries
Artifacts: MP4/SRT/poster/webp/chapter clips/review frames/report paths and hashes
Measurements: real codec, dimensions, fps, SAR, decoded frames, audio absence
Derived assets: WebP size and loop check, 5 chapter clip durations and dimensions
Visual spotlight: cue type, target coordinates, non-destructive appearance
Visual regression (SSIM): command run, min/avg SSIM scores, comparison verdict
Re-edit proof: changed caption, unchanged asset hashes, no browser/server/LOSAT
Tests: exact commands, exit status, failure cases, skipped/unverified scope
Visual review: what was actually inspected and remaining issues
Regression: protected files and affected existing code
Next session: reproduction commands and S5 acceptance work remaining
Proposed commit title: English
Proposed commit summary: English
```
