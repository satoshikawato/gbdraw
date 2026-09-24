# INSTRUCTION PROMPT — S2: 4種類の検証済み完成図を素材化する

文書版: 2.2<br>
改訂日: 2026-09-24

## 1. セッションの目的とゴール

`satoshikawato/gbdraw`の紹介動画「Meet gbdraw」に使用する**4種類の完成図素材（ヒトミトゲノム、タバコ葉緑体、Lambda–DE3比較、5件のBGC比較）**を、既存のGUIフローから取得する工程を実装します。

本セッションのゴールは、**4つの図の実際のSVG出力をダウンロードし、既存の生物学的・意味的検証を通過させた上で、映像用PNGと出所追跡記録（assets.json）を完成させること**です。

これらの素材は、完成動画の冒頭18秒間で使用される4つの導入図であると同時に、後段（S4）で生成される学会発表スライド用チャプター動画（`chapters/01_introduction.mp4`, `02_single_genome.mp4`, `03_synteny_comparison.mp4`）やREADME用ループWebPの元素材となります。

> [!NOTE]
> 後段の編集実演（Before/After）や38秒動画の最終合成はS3およびS4の管轄であり、本セッションでは扱いません。

## 2. 前提条件と依存関係

- [総合実装計画書](../01_MASTER_PLAN.md)（文書版2.2）
- `AGENTS.md`、`CLAUDE.md`
- 先行セッションS1の実装および `docs/videos/meet-gbdraw/handoffs/S1.md`
  - S1で確定した認定環境、SVG→PNG変換経路、台本モデルが動作することを確認してください。
- 既存の撮影・チュートリアルフロー:
  - `docs/capture/flows/human_circular.py` / `flows/tutorials/gui_first_circular.py`
  - `docs/capture/flows/tutorials/gui_annotated_chloroplast.py`
  - `docs/capture/flows/tutorials/gui_losatn.py`
  - `docs/capture/flows/tutorials/gui_losatp_groups.py` / `flows/bgc_losatp.py`
  - `docs/capture/assertions/`

> [!IMPORTANT]
> S1の変更を保持した作業ツリーから継続してください。`main`や`dev`への直接コミット・プッシュは行わず、無断のresetや古い参照コミットへの強制チェックアウトは避けてください。

## 3. 入力データと素材契約

### 対象素材IDと入力定義

| 素材ID | 対象ゲノム・比較 | 入力データ正本 | 映像・スライド利用目的 |
|---|---|---|---|
| `human.circular` | ヒトミトゲノム（`NC_012920.1`）の完成円形図 | `gbdraw/web/tutorial-data/human-mitochondrion/HmmtDNA.gbk` | Scene 2 / 単一ゲノム紹介・スライドClip 2 |
| `tobacco.plastome` | タバコ葉緑体（`NC_001879.2`）の領域注釈付き図 | `gbdraw/web/tutorial-data/tobacco-chloroplast/NC_001879.2.gb` | Scene 3 / 逆位反復・オルガネラ紹介 |
| `lambda-de3.comparison` | Lambda（`NC_001416.1`）対 DE3（`NC_042057.1`）の全ゲノム比較図 | `gbdraw/web/tutorial-data/lambda-de3/` | Scene 4 / シンテニー比較紹介・スライドClip 3 |
| `bgc.comparison` | 5件のBGC比較図（`BGC0000708/709/711/712/713`） | `gbdraw/web/tutorial-data/bgc-five/` | Scene 5 / 多ゲノムクラスター比較紹介 |

> [!NOTE]
> BGCは生合成遺伝子クラスター領域であり、完全染色体としては扱いません。既存設定の向き・整列・比較条件を保持してください。

### キャンバス配置と映像契約
- **基準解像度**: 1920×1080 ピクセル（16:9）
- **背景色**: 純白（`#ffffff`）
- **アスペクト比維持**: 各図本来のプロポーションを維持し、上下左右に十分なセーフエリア・余白を設ける。
- **欠損防止**: タイトル、凡例、スケールバー、シンテニー比較リンクが画面端でトリミングされないよう安全マージンを確保する。

## 4. 推奨実装手順

### Step 1: 撮影オーケストレーター（capture.py）の実装
1. `docs/capture/video/capture.py` に、上記4素材を直列に取得するオーケストレーターを実装する。
2. 動画用runごとに新規出力ディレクトリ（例: `build/videos/meet-gbdraw/run-001/`）を割り当て、各フローのダウンロード先と証拠画像出力先を完全にその配下に隔離する。
3. 既存の `CaptureWebServer` によるローカル配信、ネットワーク遮断、コンテキスト初期化を利用する。

### Step 2: 既存GUIフローの実行と意味的検証
1. 4つのフローを順次実行し、それぞれの実SVGをダウンロードする。
2. 既存の `assertions/` を再利用し、フィーチャー数、比較結果のリンク数、ダウンロード成功等の意味的検証をパスすることを確認する。
3. いずれかの取得や検証が失敗した場合、処理を中断し、古い別runの成功素材で補完しない。

### Step 3: 静的SVGから映像用PNGへの変換
1. S1で実証した閉じた描画経路（固定Chromium＋ローカルフォント）を用いて、各SVGから高解像度の映像用PNGを生成する。
2. **キャンバス配置の安定性**: 図のアスペクト比を維持し、重要要素（タイトル、凡例、スケール、比較リンク）が欠けないよう適切なビューポート・余白を設定する。
3. この段階のPNGには字幕や表紙用テキストを焼き込まない（S4の編集工程に委ねる）。

### Step 4: 素材メタデータ（assets.json）の生成
1. 取得した各素材について、以下の情報を `assets.json` に記録する:
   - 素材ID、種別（`image`）
   - 相対ファイルパス、SHA-256チェックサム
   - 画像の実測寸法（width, height）
   - 元となったSVGのパスとSHA-256
   - 取得時のアプリ実体、生成wheelハッシュ、入力データハッシュ、実行環境情報
   - 意味的検証の実行証拠への参照
2. 相対パスは素材ルート内に解決されることを検証し、`..` やルート外シンボリックリンクを拒否する。

## 5. 既知の落とし穴とガードレール

- **既存資産の保護**: `docs/images/`、`tests/reference_outputs/`、`examples/gbdraw_social_preview.png` 等を上書き・更新しない。
- **素材の混同**: S2の `human.circular` は導入で見せる既存完成図です。S3で扱うproductラベル開始状態（P0）や、D-loop領域注釈追加前（P2 / `human.dloop-before`）、追加後（P3 / `human.dloop-bracket`）とは別物として管理し、同一視しない。
- **部分成果の扱い**: S2完了時点では4素材しか存在しません。完成台本が要求する9素材のうち、残る5素材（S3の4状態＋1録画）は未取得であることを明示し、全動画の完成と偽装しない。
- **スライド素材への配慮**: 生成される高解像度PNGは、S4で発表スライド用チャプタークリップにも切り出されるため、ノイズや不要なUIクロップを含めずクリーンな図として出力する。
- **過剰な共通化の回避**: 既存のprivate helperを共有する必要がある場合は、最小限のpublic helperとして抽出し、元の呼び出し元も同時に移行・テストする。無意味な巨大抽象クラスを作らない。

## 6. 受け入れ条件（Acceptance Criteria）

1. **AC-01準拠**: 4つのGUIフローがそれぞれ実行され、既存の意味的assertionsをパスする。
2. **AC-02準拠**: 4素材のSVGおよびPNGが生成され、1920×1080解像度で実寸法とSHA-256が記録される。
3. **AC-03準拠**: `assets.json` が生成され、未知ID・不正パス・改変ファイルを拒否するバリデーションが通る。入力データのチェックサム不一致時に正しく早期失敗する。
4. 既存のドキュメント画像、参照出力、ソーシャルプレビュー画像に一切差分がない。
5. 共有ヘルパーに変更がある場合、元の呼び出し元の回帰テストがパスする。
6. 生成された4画像を実際の表示サイズで目視確認し、図の欠けや凡例・ラベル・比較リンクの破綻がない。

## 7. 引き継ぎフォーマット

セッション完了後、`docs/videos/meet-gbdraw/handoffs/S2.md` に以下を実値で記録してください。

```text
Status: PASS / PARTIAL / BLOCKED
Implementation base: branch, HEAD, dirty state, S1 changes confirmed
Completed work: capture orchestration and any shared-helper changes
Contracts: real asset IDs, schema, paths, callable interfaces
Evidence: input/figure hashes, semantic reports, environment/wheel identity
Tests: exact commands, exit status, observed failures and unverified scope
Artifacts: assets.json and four SVG/PNG/evidence paths
Visual review: each figure actually inspected and findings (layout, legends, links)
Regression: protected paths and affected existing callers
Next session: how S3 attaches the edit/export flow without duplicating fixtures
Proposed commit title: English
Proposed commit summary: English
```
