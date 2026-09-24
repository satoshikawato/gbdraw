# 総合実装計画書 — gbdraw紹介動画「Meet gbdraw」

文書版: 2.2<br>
作成・改訂日: 2026-09-24<br>
対象リポジトリ: `satoshikawato/gbdraw`<br>
対象成果物: 38秒のブラウザーGUI紹介動画1本、派生アセット（README用WebP、スライド用クリップ）、およびその再生成工程<br>
文書の状態: 実装前の仕様・作業計画。コード実装、録画、動画生成、実行検証を完了したという記録ではない。

## 1. 目的と対象読者

gbdrawは、注釈付き配列から円形・線形のゲノム図や比較図を作成するツールで、ブラウザーGUIを提供する。本計画は、初めてgbdrawを知る研究者に用途を伝える紹介動画「Meet gbdraw」を自動生成するためのものである。[R2]

動画の視聴後に伝わるべき内容は、「単独のゲノムを描ける」「複数のゲノムや遺伝子クラスターを比較できる」「GUIで図を調整してSVGを出力できる」の3点とする。全操作を習得させる入門チュートリアルではない。ファイル選択からの全手順、アルゴリズム、インストール方法は扱わない。

実装担当者は、この文書、リポジトリ、担当セッションのINSTRUCTION PROMPT、実際の引き継ぎ記録だけを使って作業する。別の動画、過去の会話、未提示の設計判断を前提にしない。

**完成条件は、見栄えのよいMP4が1回出力されることだけではない。正しい素材を取得し、その出所を追跡し、字幕だけを再編集でき、失敗を検出できる小さな制作工程を完成させることである。**

## 2. 変更権限、参照基準、実装ベース

### 2.1 参照コミットと作業ブランチを混同しない

設計時の参照コミットは `4556e04e929a4a85ad28d1833ce7304bd764881c` である。このSHAは既存資産の説明用であり、実装開始時の最新コミットや作業ブランチの指定ではない。

参照時の`AGENTS.md`は、エージェントの新規作業ブランチを最新の`origin/dev`から作成し、`main`や`dev`へ直接コミット・pushしないよう指定している。また、ブラウザー用wheelの準備手順、参照画像を通常作業で更新しない規則、所有者管理の画像を変更しない規則がある。[R1]

S1で現在の`AGENTS.md`、`CLAUDE.md`、適用範囲の追加指示を読む。ブランチ作成が必要なら現在の規則に従う。S2以降は先行実装が存在する承認済み作業ツリーを確認し、セッション切り替えを理由に変更を捨てたり、履歴を巻き戻したりしない。新規ブランチが必要な場合は、先行変更の保持・統合方法を実際の状態から確定する。資料のSHAへ無条件にcheckoutしない。

本計画が許可する作業は、対象範囲の実装とローカル検証である。別の明示指示がなければ、リモートへのpush、PR作成、merge、リリース、公開サイトの変更、SNSへの投稿は行わない。コミットも自動作成せず、最後に英語のコミット案を提示する。

### 2.2 実装前に記録するもの

現在のブランチ、HEAD、upstream、作業ツリー差分、先行セッションの変更、現在の撮影環境固定値、生成ブラウザーwheelの状態を記録する。参照コミットから移動したパスや変更されたAPIは、既存コードを読んで対応表にまとめる。

Python本体のソースとブラウザーが実行するwheelが一致するとは限らない。現在のリポジトリ手順に従い、必要なら`python tools/prepare_browser_wheel.py`を実行して準備し、使用したwheelのハッシュを残す。動画のためにキャッシュバストークンを更新したり、生成wheelをコミットしたりしない。[R1]

### 2.3 設計の変更手続き

ファイル移動、同等APIへの適応、描画領域の調整など、製品の意味や動画の約束を変えない実装判断は担当者が行い、理由と証拠を引き継ぐ。生物学的題材・3つの編集目的の変更、38秒の変更、比較データの差し替え、アプリの挙動変更、検証の緩和は黙って行わない。影響する部分を未完了として切り分け、安全に進められる部分を完成させる。

## 3. 成果物と対象外

### 3.1 必須成果物および派生アセット

| 成果物 | 仕様 |
|---|---|
| 完成動画 | `meet-gbdraw.mp4`。1920×1080、30/1 fpsのCFR、1,140フレーム、38秒。H.264、yuv420p、音声ストリームなし。Web再生向けのfaststartを使用 |
| 英語字幕 | `meet-gbdraw.en.srt`。焼き込み字幕と同じ台本から生成。発話の書き起こしではなく画面説明テキスト |
| 表紙画像 | `poster.png`。今回生成した図から作成。既存の所有者管理画像を置換しない |
| README用WebP | `meet-gbdraw.webp`。GitHub READMEトップ等での自動ループ再生用軽量アニメーション（1〜2MB前後） |
| スライド用クリップ | `chapters/*.mp4`。各完成図（Circular, Plastome, Linear, Clusters, Edit）の単体ループMP4 |
| 原素材 | 完成図SVG、映像用PNG、字幕を焼き込んでいない操作録画。実際の録画コンテナー・codecを記録 |
| 素材記録 | `assets.json`。素材の出所、ハッシュ、撮影条件、検証証拠を参照 |
| 制作・検証記録 | 環境、台本ハッシュ、元素材ハッシュ、出力情報、SSIM自動測定レポート、失敗理由、確認用フレーム |
| 制作コードと手順 | Python PlaywrightとFFmpegを中心とする制作コマンド（build/render/check）、テスト、README |

MP4のバイト列の完全一致は要件にしない。再現性の対象は、生物学的内容、操作の意味、編集構成、フレーム数、表示品質、素材の出所である。

### 3.2 対象外

縦動画、多言語化、ナレーション、BGM、SNS自動投稿、紹介動画以外のチュートリアル制作、一般的な機能推定・注釈分類エンジン、単一遺伝子の恣意的な色変更を主役にする実演、汎用動画編集GUI、LLMによる毎回の自律操作、プラグイン式レンダラー、並列撮影、キャッシュ探索・増分ビルド基盤、CI必須ゲート化は作らない。別言語のPlaywright実装、Remotion、MoviePy、OBS制御を追加しない。gbdrawの通常利用にFFmpegを必須にしない。

## 4. 動画の仕様と題材

### 4.1 紹介する編集上の価値、視線誘導、視聴者体験

編集パートは、同じヒトミトゲノムについて「表示する名称を選ぶ」「全CDSを機能カテゴリ別に整理する」「原点をまたぐD-loop領域注釈ブラケットを追加する」の3つを示す。個別フィーチャーの装飾ではなく、研究用の図を論文品質へ整える作業を題材とする。

順序は、product→geneラベル、CDSの機能別配色、D-loopブラケットの追加とする。1組のBefore/Afterでは変更する意味を1つに限定し、配色・ラベル・表示集合を同時変更して差の原因を分からなくしない。productを使う図も正当な選択肢であり、変更前を意図的に低品質にしない。

**静止画フリーズ誤認の防止とテンポ感**:
同一素材（P1およびP2）を連続して用いる移行部（Scene 7→8、Scene 9→10）では、字幕を単なる状態の説明にとどめず、次のステップへの目的予告（例: 「Categorize features by biological role.」「Add custom region annotations.」）とすることで、視聴者に動画が停止したと誤認させず、スムーズな編集ステップとして知覚させます。

**視線誘導（Visual Spotlight）**:
2秒のAfterシーンにおいて、変化した領域（代表geneラベル、凡例、D-loopブラケット）への控えめな視線誘導（微細なアテンション枠やフェード強調）を施し、視聴者が瞬時に変化の核心を理解できるようにします。

**研究者ベネフィットの訴求**:
研究者にとって最大の動機である「ベクターSVGによる論文品質の出力（Publication-ready vector SVG）」および「ブラウザー完結・インストール不要（Runs in your browser）」を明確に訴求します。

### 4.2 絵コンテの規範値

38秒・1,140フレームを維持し、完成図の紹介に18秒、3組のBefore/Afterに12秒、出力録画に4秒、末尾に4秒を割り当てる。時間区間は開始を含み終了を含まない。以下の秒は説明値で、実装上の時刻は整数フレーム数の累積から算出する。

| Scene ID | 時間 | フレーム数 | 映像 | 英語字幕 | 意図・演出 |
|---|---:|---:|---|---|---|
| `intro` | 0–2秒 | 60 | 4種類の完成図を2×2で配置 | Meet gbdraw — Publication-ready genome figures | ツール名と「論文品質」を一瞬で認知 |
| `circular` | 2–5秒 | 90 | ヒトミトゲノムの円形図 | Draw circular genomes effortlessly. | 導入図1：円形ゲノムを簡単に描画 |
| `plastome` | 5–8秒 | 90 | タバコ葉緑体の領域注釈付き完成図 | Annotate plastomes with inverted repeats. | 導入図2：葉緑体と逆位反復の注釈 |
| `genome-comparison` | 8–13秒 | 150 | Lambda–DE3比較図 | Align and compare whole genomes. | 導入図3：全ゲノム比較・アライメント |
| `cluster-comparison` | 13–18秒 | 150 | 5件のBGC比較図 | Visualize biosynthetic gene clusters. | 導入図4：BGC遺伝子クラスター比較 |
| `labels-product` | 18–20秒 | 60 | P0：全CDSのproductラベル | Too crowded? Full product labels. | Before 1: 文字が多く混雑した状態 |
| `labels-gene` | 20–22秒 | 60 | P1：同じCDSのgeneラベル | Switch to clean, readable gene symbols. | After 1: 簡潔で可読性の高いgene名へ（代表geneを視線誘導） |
| `colors-before` | 22–24秒 | 60 | P1を再利用：CDSは単一色 | Categorize features by biological role. | Before 2: 単一色。次の配色への目的を予告 |
| `colors-functional` | 24–26秒 | 60 | P2：全13 CDSを4機能カテゴリで配色 | Color CDS by respiratory complex. | After 2: 呼吸鎖複合体ごとの色分け（凡例とCDSを視線誘導） |
| `dloop-before` | 26–28秒 | 60 | P2を再利用：領域注釈なし | Add custom region annotations. | Before 3: 領域注釈の追加を予告 |
| `dloop-bracket` | 28–30秒 | 60 | P3：原点をまたぐD-loopブラケット追加 | Origin-spanning D-loop bracket. | After 3: 美しい弧状ブラケットが出現（ブラケット部を視線誘導） |
| `export-svg` | 30–34秒 | 120 | P3を実際にSVGダウンロードするGUI録画 | Export crisp, scalable vector SVG. | ベクターSVG保存の実画面 |
| `outro` | 34–38秒 | 120 | その実ダウンロードSVGと利用開始案内 | Runs in your browser. Try gbdraw.app | インストール不要・Web完結の最強メリット |

本編は13 sceneで構成する。各Before/Afterは同じ構図の2枚を2秒ずつ表示する。クリック手順を4秒に圧縮するのではない。3つの編集操作は実GUIで最初から最後まで実行・録画して証拠を残すが、その全録画を本編へ入れる要件はない。本編の変更前後は実際に生成・ダウンロードしたSVGから作る静止画であり、連続操作画面に見せかけるカーソルや架空UIを合成しない。

「instant」「real-time」「one click」「38秒で解析完了」などを使わない。BGCのSimilarity groupを機能同一性やオーソロジー確定の証明として説明しない。機能別配色は既存gene注釈と明示ルールを用いる表示処理であり、gbdrawが配列から機能を新たに推定したとは説明しない。

### 4.3 素材と既存フローの対応

| 素材ID | 内容 | 根拠・担当 |
|---|---|---|
| `human.circular` | ヒトミトコンドリア`NC_012920.1`の既存完成図 | S2、`flows/human_circular.py`、`flows/tutorials/gui_first_circular.py` [R3] |
| `tobacco.plastome` | タバコ葉緑体`NC_001879.2`、領域注釈付き | S2、`flows/tutorials/gui_annotated_chloroplast.py` [R6] |
| `lambda-de3.comparison` | `NC_001416.1`対`NC_042057.1`、実ブラウザー検索 | S2、`flows/tutorials/gui_losatn.py` [R4] |
| `bgc.comparison` | `BGC0000708/0709/0711/0712/0713`の5件のBGC | S2、`flows/tutorials/gui_losatp_groups.py`、`flows/bgc_losatp.py` [R5] |
| `human.labels-product` | P0：productラベル・CDS単一色・領域注釈なし | S3。変更前SVGとPNG |
| `human.labels-gene` | P1：CDSラベルだけをgeneへ変更 | S3。ラベル変更後かつ機能別配色前の共通素材 |
| `human.functional-colors` | P2：P1の全13 CDSを4機能カテゴリに配色 | S3。配色後かつD-loopブラケット追加前の共通素材 |
| `human.dloop-bracket` | P3：P2に原点をまたぐD-loop領域注釈ブラケットを追加 | S3。出力録画で実際に得た最終SVGとPNG。末尾にも同じ素材を使う |
| `human.svg-export` | P3の実SVGダウンロードを伴う録画 | S3。4秒（120フレーム基準）、保存内容まで検証 |

完成台本には一意な9素材を要求する。P1とP2の再登場、P3の末尾への利用では同じ素材IDを再利用し、同じ状態の複製素材を増やさない。編集の全操作録画は`evidence/`から参照する証拠であり、本編用の必須動画素材を3本増やすものではない。

既存フローのパスは`docs/capture/`からの相対パス。入力の順序・チェックサム・比較条件の正本は既存manifestとフローとする。BGCはクラスター領域であり完全染色体として扱わず、既存の向き・整列・比較条件を保持する。[R5]

### 4.4 表示状態の契約

| 状態 | CDSラベル | CDS配色 | D-loop領域注釈ブラケット | 描画フィーチャー・トラック |
|---|---|---|:---:|---|
| P0 | 全13個でproduct | 1色 | なし | 37論理フィーチャー（CDS, rRNA, tRNA） |
| P1 | 全13個でgene | P0と同じ1色 | なし | 37論理フィーチャー（ラベル採用元のみ変更） |
| P2 | P1と同じgene | 全13個を4カテゴリへ | なし | 37論理フィーチャー（配色と凡例のみ変更） |
| P3 | P2と同じgene | P2と同じ4カテゴリ | **あり** | **37論理フィーチャー ＋ 1個のD-loopブラケット（内側トラック）** |

すべての状態で元GenBankレコード（配列、注釈）を維持する。P3では、公式チュートリアル（`highlight-mitochondrial-features.md`）に準拠し、`mitochondrial_regions.tsv` による原点越えブラケット（`NC_012920.1: 16024..576, wraps_origin=true`）が内側トラックに正しく描画されたことを検証する。[R7]

## 5. アーキテクチャと責任境界

### 5.1 パイプライン

```text
固定入力＋実装ベースのgbdraw＋生成ブラウザーwheel
                        │
                        ▼
既存GUIフロー＋3つの表示調整（Python Playwright）
                        │
          ┌─────────────┼─────────────┐
          ▼             ▼             ▼
        SVG/PNG       操作録画        検証証拠
          └─────────────┼─────────────┘
                        ▼
                    assets.json
                        │
             ┌──────────┴──────────┐
             │                     │
             ▼                     ▼
         編集処理            storyboard.json
      （Python＋FFmpeg）       順序・尺・字幕
             │
             ▼
         MP4 / SRT / poster / 確認用フレーム
             │
             ▼
       技術検証＋意味の検証＋実際の目視確認
```

素材取得と編集を別の工程にする。`render`は素材を読み、必要なレイアウトと字幕を作り、映像を合成するだけである。ブラウザー起動、ローカルアプリ起動、LOSAT実行、入力データ再取得、SVG再描画を行わない。SVGから映像用PNGへの変換は素材取得工程で完了させる。

### 5.2 最小のファイル構成

以下は新規ファイルの配置案。現在のリポジトリ構造に適合させる場合は、責任と依存方向を維持し、S1の記録で実配置を確定する。

```text
docs/capture/
  build_video.py           # コマンド入口・依存の組み立て
  video/
    __init__.py
    model.py               # 台本・素材のデータ型と入力検証
    capture.py             # 既存フロー呼び出し・素材取得
    render.py              # 静止画配置、字幕、FFmpeg編集
    validate.py            # 制作結果の検証とレポート
  flows/                   # 既存。実際に重複する操作だけ最小抽出
  assertions/              # 既存。意味的検証を再利用
  web_server.py            # 既存。ループバック配信

docs/videos/meet-gbdraw/
  README.md
  storyboard.json
  environment.json        # S1で実測確認した制作環境の固定値
  handoffs/                # 各セッションの人間が読める引き継ぎ

tests/                     # 現在の配置慣例に従って動画関連テストを追加
```

編集処理は`Page`やgbdrawの内部状態を受け取らない。素材取得は字幕・出力尺を操作列の正本にしない。既存GUIフローや製品コードから`video/`へ逆向きに依存させない。CLIのサブコマンド選択前にPlaywrightや重い撮影モジュールを無条件importしない。

新規の抽象基底クラス、DIコンテナー、イベントバス、ジョブキューは作らない。通常の関数、小さなdataclass、明示的な引数で実装する。外部処理の境界だけをテストで差し替えられるようにする。

### 5.3 既存基盤の再利用

ループバックサーバー、通信遮断、入力チェックサム、生成完了待ち、SVG意味検証、実ダウンロード検証は既存実装を使う。[R2][R8]

既存フローが途中画像を要求する場合は動画の作業領域へすべて出力する。不要な画像が多少増えることを理由に、全フローへ動画専用フラグを追加しない。既存の`docs/images/`へ出力してから拾い直す方法は禁止する。

共有ヘルパーを抽出するのは、実際の重複や既存APIの境界が必要になった箇所だけとする。抽出時に元の呼び出し元も同時に移行する。セレクター列のコピー、monkeypatch、ダミーのPageによる既存フロー偽装、別言語への全面移植は行わない。

## 6. データ契約と知識の所有場所

### 6.1 正本

| 知識 | 正本 |
|---|---|
| 入力データのID・ファイル・チェックサム | 既存tutorial-data manifestと対応設定 |
| GUI操作・比較条件・図の意味的検証 | 既存フロー、最小限の共有ヘルパー、既存assertions |
| geneラベル・機能別配色の規則 | 既存Priority TSV・feature-presentationの定義。動画で必要なCDS部分だけを再利用 |
| product開始状態・3つの変更の順序 | S3の明示的な実演レシピ。既存GUIの初期値とは呼ばない |
| D-loopの生物学的identity・区間 | 固定GenBank内のD-loop feature。注釈ブラケットの別TSVは追加しない |
| 動画の順序・字幕・フレーム数 | `storyboard.json` |
| 制作環境の追加固定値 | `environment.json`。既存のPlaywright固定値は重複定義せず参照 |
| 取得された素材の実体と出所 | `assets.json`と対応証拠 |
| 完成動画の生成条件と結果 | 実行レポート |

説明書の時間表は人間向けの規範仕様。コードの字幕生成と編集時刻計算は必ず`storyboard.json`から行う。レポートに同じ情報を記録することは証跡であり、設定の別の正本ではない。

### 6.2 台本契約

schema version、動画ID、出力形式、期待総フレーム数、scene配列を持つ。sceneはID、素材ID配列、`single`または`grid2x2`の配置、整数フレーム数、英語字幕を持つ。2×2配置は今回の表紙に必要な固定レイアウトであり、任意レイアウト言語へ拡張しない。動画素材は`single`のみ。

sceneの定義例:

```json
{
  "id": "colors-functional",
  "layout": "single",
  "assets": ["human.functional-colors"],
  "frames": 60,
  "caption": "Color CDS by function."
}
```

字幕にクリックセレクター、任意Pythonコード、FFmpegフィルター式、配列ハッシュを埋め込まない。scene順序から開始フレームを累積計算する。固定の受け入れ要件1,140フレームと合計を独立に照合し、単に合計値を期待値へコピーするテストにしない。

### 6.3 素材契約

素材ごとに一意ID、`image`/`video`の種別、ファイルの相対パスとSHA-256、実測サイズ、動画なら時間・実測codec・時間基準、素材取得run ID、対応する証拠のパスとハッシュを持つ。図のPNGには元SVGのパスとハッシュを関連付ける。

出所には、アプリの実際のHEAD、関連ソースの状態を特定する情報、dirty状態と差分の識別情報、生成ブラウザーwheelのハッシュ、使用入力manifest・各入力のハッシュ、撮影レシピの版、Playwright・Chromium・OS・Python・フォント・viewport情報を残す。HEADだけで未コミット変更を識別したことにしない。公開用レポートにユーザー名、秘密、不要な絶対パスを含めない。

`assets.json`は使用できる素材を明示列挙する。ディレクトリから「最新っぽいファイル」をglobで拾わない。相対パスは素材バンドルのルート内に解決されなければならず、`..`、不正な絶対パス、ルート外へ抜けるsymlinkを拒否する。ファイル名に空白があること自体はエラーにしない。

S1/S2で部分的な素材集合を作ることはできる。ただし、最終台本が必要とする全素材を満たすまでは完成動画として合格させない。同一フルbuildの素材は同じアプリ実体・入力契約から取得し、古い別runの成功素材を自動で混ぜない。

P0→P1→P2→P3の各遷移について、前後の素材ID、入力identity、操作録画のpath/hash、実際のUI操作、許容する変更、意味的検証結果を証拠へ記録する。P3→P2→P3の表示復元検証も残す。証拠は撮影側が作り、rendererは参照とハッシュだけを検査し、生物学的規則を二重実装しない。

### 6.4 検証記録

単なる`verified: true`だけで合格扱いにしない。検証名、実行コマンドまたは使用関数、入力と出力のハッシュ、終了状態、対象の意味的結果、失敗理由を保存する。編集工程は、生物学的解析を再実行せず、素材と証拠のハッシュ結合を検証する。この違いをレポートへ明記する。

## 7. 制作環境と描画方法

### 7.1 固定値の扱い

参照版の撮影基盤はPython Playwright 1.61.0、Chromium 149.0.7827.55、1440×900、device scale 1、en-US、UTC、light等を固定している。実装時には現在のコードから読み直す。動画だけを理由に固定バージョンを更新しない。[R2]

初回の認定制作環境はLinuxの1環境だけとする。Windows/macOS全環境での同一画素保証は対象外。OS、Python、Pillow、FFmpeg/ffprobeのバージョンとbuild configurationを記録し、S1の実行で成立した組合せを固定する。ランタイム中に不足依存を自動インストールしない。

FFmpegはH.264エンコード、字幕焼き込み、必要な映像フィルターを備えたビルドを1つ利用する。初回の字幕描画経路は`subtitles`/libassとし、S1で利用可能性を確かめる。コンテナーの導入自体は必須にしない。既存の再現可能な環境管理が使えればそれを使う。

### 7.2 フォントはGUI用と編集用を区別する

参照版のWeb用InterはWOFF2で配置されている。これをPillowやlibassがそのまま利用できると仮定しない。[R9]

GUIとSVGの描画には既存のブラウザー用フォントを利用する。字幕・表紙の編集には、制作環境にある利用条件を確認したTTF/OTF書体を1つ選び、実ファイルのハッシュ、書体名、供給元、利用条件をS1で記録する。自動的に任意のシステムフォントへフォールバックしない。必要なフォントがなければ制作環境の前提不足として明示する。フォントファイルを動画配布物や計画書パッケージへ同梱しない。

### 7.3 完成図のPNG化

検証済み静的SVGを固定Chromiumで映像用に描画する。低解像度の既存チュートリアルPNGを拡大して完成素材にしない。図の縦横比と全体の必要要素を維持し、フォント読み込み完了とレイアウト安定を待つ。

既存`CaptureWebServer`はパッケージ済みWebルートのループバック配信を担当する。[R8] 生成SVGは通常その外にあるため、実装では、SVG本文と許可されたローカルフォントだけを描画する閉じたページ、または同等に限定したローカル配信を用いる。生成物を製品ディレクトリへコピーしたり、リポジトリ全体を公開するサーバーへ変えたりしない。

静的SVGのscript、外部参照などは既存安全検証を通す。描画補助HTMLは素材変換専用とし、アプリの動作を偽装する画面には使わない。図のSVG構造や生物学的座標を編集ソフト側で修正しない。

### 7.4 実画面録画

Playwrightによる実画面録画は、公式のブラウザーコンテキスト録画機能 `browser.new_context(record_video_dir=..., record_video_size={"width": 1920, "height": 1080})` を標準・推奨経路とする。あるいは、低レベルのCDP Screencast API（`start_screencast`/`stop_screencast`）も利用可能である。いずれの場合も、生成された実録画ファイルから寸法、フレームレート、コーデック、コンテナー情報をffprobeで測定し検証する。[W1]

既存ドキュメント撮影の1440×900設定は変更しない。動画専用ページで1920×1080を指定し、S1でレイアウトと実録画サイズを確認する。既存のviewport依存ヘルパーを利用した際にサイズが戻らないことも検証する。実寸の不足を単に高解像度へ拡大して隠さない。

カーソルは標準の操作表示またはブラウザー既定のポインター挙動を利用し、独自JavaScriptによる装飾を合成しない。字幕・章タイトルは録画中に追加しない。標準操作表示に含まれる装飾は環境・素材記録に残す。必要なら実測した固定矩形で操作箇所を切り出せるが、追尾カメラや動的ズームは作らない。切り出し矩形は元録画の座標系と実寸に対応させ、重要操作を画面外にしない。

## 8. 生物学的に意味のある表示調整の契約

### 8.1 同じ実レコードから4状態を作る

S3では固定の`HmmtDNA.gbk`（`NC_012920.1`）を読み込み、既存の完成図レシピを土台に、§4.4のP0→P1→P2→P3を同じアプリ状態で順に作る。4枚を互いに無関係なレシピから作ってBefore/Afterとして連結しない。図の向き、スケール、トラック、rRNA/tRNAの表示・色、通常のラベル有無は、3つの変更と関係なく動かさない。

P0は実演のために設定した開始状態であり、アプリのデフォルトを再現したとは説明しない。全13 CDSのproductを正しく表示し、geneラベルや機能別配色を先に適用しない。D-loopのブラケット等、元フィーチャーとは別の領域注釈も置かない。開始状態を不自然に悪く見せる色・サイズ・ラベル衝突の演出は禁止する。

元ファイルのハッシュ、sequence、全source featureの型・区間・strand・qualifierを記録する。sourceのフィーチャー集合と表示するフィーチャー集合を区別し、D-loopを非表示にしてもソースから削除しない。

### 8.2 productからgeneへ：ラベルの採用元を変える

既存の**Labels / Priority File (TSV)**等の実UIを用いる。各CDSのproduct/geneが入力に存在することを確認し、P0ではproduct、P1ではgeneが実際に使われていることを全13 CDSで検証する。gene用TSVは既存の`CDS\tgene`を再利用する。product用の開始状態は同じ既存schemaで`CDS\tproduct`に相当するレシピを1か所に置く。[R7][R12]

13ラベルを個別に書き換える方法ではなく、qualifier優先順位の変更として実演する。入力の`product`を`gene`へ改名・削除しない。遺伝子名が存在しない別データへ暗黙に拡張する機能も作らない。

許容する変更はCDSの表示文字列、文字長に応じたラベル配置・leader・canvas bounds等の正常なレイアウト効果である。CDSのidentity、色、座標、D-loop状態、rRNA/tRNAのラベル・色を意図せず変えない。単に画像中にgene名が一つあるだけでは合格にせず、論理CDSごとの文字列対応を確認する。

### 8.3 全CDSを機能カテゴリ別に色分けする

既存feature-presentation教材の機能別配色定義を正本として、以下のCDS4群を適用する。[R7]

| カテゴリ | 対象gene | 期待CDS数 | 教材の凡例名 |
|---|---|---:|---|
| NADH dehydrogenase | ND1、ND2、ND3、ND4、ND4L、ND5、ND6 | 7 | NADH dehydrogenase |
| Cytochrome c oxidase | COX1、COX2、COX3 | 3 | Cytochrome c oxidase |
| ATP synthase | ATP6、ATP8 | 2 | ATP synthase |
| Cytochrome b | CYTB | 1 | Cytochrome b |

正規表現等の適用では完全一致を使い、ND4とND4L等を誤って重複分類しない。各CDSがちょうど1カテゴリに入り、13個すべてを覆い、未分類CDSが残らないことを検証する。既存教材のrRNA配色行は今回の変更対象から外し、CDSの4行だけを再利用する。フィルターして得たTSVと元定義のハッシュ・対応を証拠へ残す。既存定義を共有できない場合は、その定義だけを最小抽出し、既存呼び出し元も同時に移行する。

**Colors / Specific Table (-t)**等の実UIから規則を適用し、必要なGenerateを実行する。一つずつ13回クリックする、全CDSを同じ新色にする、ソースを変更する、動画編集で色を塗る、といった代用はしない。4群のfillと凡例色・captionが一致し、正常な凡例の再構成以外でrRNA/tRNA、ラベル、配列座標が変わらないことを確認する。残存CDSを隠して「全CDS」の見かけを作ってはいけない。

色とカテゴリの対応はユーザーが与える表示規則であり、配列からの機能推定ではない。新たな機能注釈サービス、LLM分類、汎用オントロジー基盤は作らない。

### 8.4 原点をまたぐD-loop領域注釈ブラケットを追加する

ヒトミトゲノム（`NC_012920.1`）のD-loopは原点（16,569 bpと1 bpの境界）をまたぐ非コード領域（`16024..576`、1,122 bp）である。本実演では、公式チュートリアル（`highlight-mitochondrial-features.md`）の正規手順に基づき、**Region Annotations（領域注釈）パネルを用いてD-loopの弧状ブラケット（bracket）を追加する操作**を実演する。[R7]

操作手順としては、Web UIの **Region Annotations** パネルから既存チュートリアルデータ `mitochondrial_regions.tsv`（`d_loop\tbracket\tNC_012920.1\t16024\t576\tsource\ttrue\tD-loop`）をインポートし、**Custom Track Slots** で内側トラックへバインドして **Generate Diagram** を実行する。

P2→P3で、37個のCDS/rRNA/tRNAは同じまま維持され、内側トラックに原点をまたぐ（`wraps_origin=true`）1個のD-loopブラケットとラベル「D-loop」が追加描画されたことを検証する。配列の回転・cropや座標変更を伴わない。

検証用には、ブラケットトラックを削除（または非表示）にしてP2へ戻し、再追加でP3へ復元できる再現性を確認する。この往復全操作を本編で繰り返す必要はない。

### 8.5 実行証拠と本編用素材を分ける

3つの変更は、必要なUI操作・Apply・Generate・結果確認を含めてそれぞれ最後まで録画する。その録画は証拠として保存し、4秒という本編の表示時間に収めるために操作を省かない。準備は通常のUIで行い、内部状態代入、直接SVG属性変更、比較キャッシュ注入、架空の成功画面は使わない。内部状態の読み取りを待機・検証に用いることは区別する。

各状態で実際のSVG出力を取得して意味的に検証し、そこから映像用PNGを作る。P3のPNGは後述の出力録画でダウンロードしたSVGから作る。本編ではP0/P1、P1/P2、P2/P3の3組を各4秒で比較する。2枚の切り替えに架空のクリック動作を合成せず、図のBefore/Afterであることを字幕とカットで示す。

4状態の映像フレームは、同じ円の中心、同じbpに対する縮尺、同じ向きで比較できるようにする。正常なlabel/legend reflowは許容するが、画像ごとに独立した最大fitをして円の大きさを変えない。素材変換時に4状態のboundsを合わせた共通キャンバス・余白を使う。SVGの科学的内容は改変しない。最初から全図が読めなければ、元の開始レイアウトを実GUIで整え、4状態を取り直す。

これは3種類の既知の操作と4状態を扱うレシピで足りる。汎用Before/Afterエンジン、動画イベント同期、複雑な自動カット、任意状態遷移システムを追加しない。

### 8.6 最終SVG出力と待機

P3を最後に復元した同じアプリ状態から、SVGボタンを押して実ダウンロードを完了させる録画を取得する。事前に通常のUIでボタンを見える位置へスクロール・配置し、録画開始からダウンロード完了までの不要なアイドル時間を排除する。`expect_download`等でファイルの非空性・静的SVG安全性・全CDSのgeneラベルと4群配色・D-loop表示・source identityを検証する。末尾とD-loop表示後のsceneは、この同じSVG由来の`human.dloop-visible`を使う。

処理完了は状態assertionと有限timeoutで待ち、演出用の停止と区別する。3つの編集証拠録画の長さに本編4秒の制限を課さない。本編へ使う出力クリップは、30 fps正規化後120フレーム（4秒）を基準とし、実操作が120フレームより短い場合は成功した末尾フレームを複製・静止して120フレームに延長する。操作完了がわずかにオーバーした（数フレーム程度）場合、ダウンロード完了後の不要な余白を末尾トリミングするか、微細なPTS正規化を行って安全に120フレームへ揃える。大幅な超過（操作完了自体が4秒に収まらない場合）は途中切断で誤魔化さず、UI操作手順の事前最適化（ボタンの可視化や不要待ちの削除）を行う。

録画停止・操作装飾解除・context終了をfinallyで保証する。失敗録画は診断用に保持できるが、検証済み素材へ登録しない。

## 9. 編集工程

### 9.1 画像・映像・字幕

静止画と短い録画の2素材種別、`single`と表紙用`grid2x2`の2配置に限定する。表紙、終了画面、字幕は編集工程で作る。表紙・posterの図は素材IDで参照し、既存のソーシャルプレビューを加工・置換しない。

図の縦横比を保って配置する。操作を見せる画面と図だけを見せる画面を区別し、字幕のための安全領域を確保する。S1で字幕用書体、文字サイズ、行間、余白、色を1組固定する。重要文字の縮小限界を決め、はみ出した字幕を黙って極小化しない。

subtitleのテキストは台本からUTF-8のSRTへ生成し、同じ生成物を焼き込みに使用する。ASS/libass固有の装飾制御が必要なら生成処理側で固定し、字幕本文に制御コードを許可しない。Pillowは表紙や文字領域の事前チェックに利用できるが、実際の字幕レンダリング結果もフレームで確認する。

### 9.2 FFmpegの呼び出し

引数配列による`subprocess.run`等を用い、`shell=True`を使わない。外部コマンドの終了コード、stderr、timeoutを記録する。字幕本文をフィルター式へ文字列連結しない。ファイルを経由してもFFmpeg独自のパス・フィルターのエスケープは必要なので、空白、日本語パス、引用符、コロン、バックスラッシュをテストする。

Before/Afterの6 sceneは通常の静止画sceneであり、P1とP2は同じ素材IDを再利用する。証拠用の長い操作録画を本編素材として誤って連結しない。

各場面を共通のfps、寸法、sample aspect ratio、pixel format、開始PTSへ正規化し、カットで連結する。複雑なトランジション、音声同期、自動的な速度調整は実装しない。出力はffprobeと実デコードで検証する。[W2][W3]

### 9.3 原素材を変更しない

`render`は入力素材とその証拠を読み取り専用で扱う。字幕差し替え・場面順序変更・割当フレームの変更で、`assets.json`、SVG、PNG、録画のハッシュが変わってはいけない。製品承認用の最終版は38秒の固定台本で検証する。

## 10. CLI、出力領域、失敗の扱い

以下のコマンドは新規に実装するインターフェースの仕様であり、既存コマンドだとは仮定しない。

```bash
# 固定入力から全素材を取得し、編集と検証を実行する
python docs/capture/build_video.py build \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/run-001

# 既存素材を使って再編集する。新しい出力先を指定する
python docs/capture/build_video.py render \
  --assets build/videos/meet-gbdraw/run-001/assets.json \
  --storyboard docs/videos/meet-gbdraw/storyboard.json \
  --out build/videos/meet-gbdraw/render-002

# 完成動画、出力ハッシュ、字幕、記録された証拠を確認する
python docs/capture/build_video.py check \
  --out build/videos/meet-gbdraw/run-001

# SSIMによる代表フレームのビジュアルリグレッション自動検査（閾値 SSIM >= 0.98）
python docs/capture/build_video.py check \
  --out build/videos/meet-gbdraw/run-001 \
  --visual
```

`build`と`render`は空または新規の出力先を要求する。既存の成功runを上書きしない。初回は`--force`や自動バックアップ管理を作らない。新しいrun名で再実行できれば十分である。

```text
build/videos/meet-gbdraw/run-001/
  assets.json
  figures/                 # SVG・映像用PNG
  raw/                     # 操作録画
  evidence/                # 既存フローの画像・意味的検証・ダウンロード
  final/
    meet-gbdraw.mp4
    meet-gbdraw.en.srt
    poster.png
    meet-gbdraw.webp       # GitHub README用軽量ループアニメーション
    chapters/              # スライド用単体チャプターMP4クリップ群
  reports/
    build-report.json
    validation.json
    ssim-report.json       # SSIM自動測定結果
    review-frames/
```

再編集runは原素材を複製する必要がない。`reports/`に、使用した元manifestのハッシュと出所・検証結果のスナップショットを残す。`check`は完成メディアとスナップショットの結合を確認するもので、ブラウザー解析を再実行した証明ではない。元素材そのものを再監査する場合は元バンドルも必要である。この区別を結果に示す。

作業中は一時領域へ出力し、必要な検証がすべて通ってから完成物として確定する。失敗したrunは非ゼロ終了とし、成功した別runの動画を返さない。部分素材・ログは原因調査用に残してよいが、最終成功の印を付けない。大容量バイナリーを通常のソース履歴にコミットしない。

## 11. 実装セッションと段階ゲート

セッションは同じ実装を順次積み上げる。先行完了は引き継ぎ文だけでなくコードと成果物で確認する。将来の機能の空実装を先に増やさない。

| セッション | 作業 | 完了ゲート |
|---|---|---|
| S1 | 環境・実装ベースの確認、台本契約、静止画→試験MP4、短いScreencast、字幕・フォントの実証 | 実環境で映像経路が成立し、認定環境と最小のデータ契約が記録される |
| S2 | 4種類の完成図取得、SVG→PNG、素材manifest、既存assertionsとの結合 | 4素材が正しく生成され、ドキュメントの承認済み画像が不変 |
| S3 | P0→P1→P2→P3の実GUI調整、全操作録画、状態SVG/PNG、最終SVG出力 | 全13 CDSのラベル・機能別配色、元D-loopの表示往復、source保存を検証し、9素材を満たす |
| S4 | 台本による表紙・字幕・場面編集、38秒のMP4、build/render/checkの統合 | 全素材から最終動画を作れ、再編集時にブラウザーやLOSATが起動しない |
| S5 | 回帰・失敗系・2回の独立build・再編集の検証、目視、READMEと引き継ぎ | 技術・意味・視認性を分けて報告し、完成の根拠を残す |

各セッションは計画説明だけで終えず、担当範囲の実装と実行検証を行う。環境制約がある場合は、実際の診断と未検証範囲を記録し、安全に実行できるテストを完了する。合格していないゲートを合格したことにしない。

## 12. 検証と受け入れ条件

### 12.1 自動検証

| ID | 受け入れ条件・失敗系 |
|---|---|
| AC-01 | 元レコード、フィーチャー集合、比較結果、領域注釈が既存の意味的検証を通る |
| AC-02 | 素材ID・path・checksum・証拠を検証。未知ID、重複、不在、改変、ルート外参照を拒否する |
| AC-03 | P0→P1→P2→P3が1意味ずつ変化。全13 CDSのproduct→gene、4カテゴリ7/3/2/1の全件配色、原点をまたぐD-loopブラケット追加を検証 |
| AC-04 | 各状態SVGが実UI出力に対応しsourceは不変。P0〜P2の37個およびP3の37個+1ブラケットの描画集合を照合。最終SVGは出力録画・本編P3・末尾と同一 |
| AC-05 | 完成動画は1920×1080、30/1 CFR、SAR 1:1、H.264、yuv420p、音声なし、実デコード1,140フレーム |
| AC-06 | scene合計1,140。PTSやdurationの丸めは1フレーム以内の説明可能な差のみ許容し、余分な映像フレームを許容しない |
| AC-07 | SRTと焼き込みテキストが同じ台本由来。フレーム境界から時刻生成し、空白・重複・はみ出しを確認する |
| AC-08 | 本編の出力録画は30 fps正規化後120フレーム（4秒）基準。短ければ成功末尾を延長し、微小超過は安全に末尾トリミングまたは正規化して120フレームに揃える。3つの編集の全録画は証拠として別保存し4秒上限を課さない |
| AC-09 | `render`はブラウザー、サーバー、LOSATを起動せず、入力素材のハッシュを変更しない |
| AC-10 | FFmpeg/codec/font不足、外部通信、GUIエラー、生成timeout、download失敗、壊れた動画で非ゼロ終了 |
| AC-11 | 既存の承認画像、reference outputs、所有者管理画像に差分がない。共有ヘルパーの既存呼び出し元も回帰確認 |
| AC-12 | 同一ソース・環境で、独立した空の出力領域と新規ブラウザーcontextによるフルbuildが2回成功する |
| AC-13 | 13 Scene代表フレームについて、承認済みリファレンス（またはRun間）のSSIM自動判定（`SSIM >= 0.98`）をパスする |

意味のない自己検証を避ける。例えばmanifestへ書いたフレーム数を読むだけで動画を測定したことにしない。ffprobeの実測、完全デコード、source/targetの意味的比較など異なる経路の証拠を使う。

テストは純粋関数のunit、FFmpegの小さなintegration、対象GUIのend-to-end、既存共有フローのregressionに分ける。通常の全テストへ高コストの動画buildを無条件に追加しない。構造・単体テストは軽く、フル制作は明示実行とする。

### 12.2 視認性と内容確認

完成動画を実際に再生し、各sceneの代表フレーム、編集の前後、接続境界の前後フレームを確認する。元の1080pに加え、960×540程度の縮小表示でも字幕と主張が読み取れるか確認する。すべての遺伝子ラベルを小画面で読めることまでは要求しないが、何を示す図か、どの操作で何が変わったかは識別できなければならない。

図やUIの文字が字幕・カーソルに隠れないこと、product/geneラベル、全CDSと凡例の機能別配色、D-loop表示の3つの前後差を追えること、凡例・重要リンクが不自然に切れないこと、比較の意味や処理速度を誤認させないことを確認する。全白・全黒検出だけを品質保証にしない。静止画中心の動画なので、静止フレームが長いこと自体を不具合扱いしない。

確認者は実際に見た範囲だけを記録する。動画再生できない環境では、フレーム確認済みと全編再生未実施を分ける。自動検証PASSと公開承認は別であり、公開は人の判断とする。

## 13. SOLID・KISS・DRY・YAGNIの適用

| 原則 | コードでの適用 | アーキテクチャ・ワークフローでの適用 |
|---|---|---|
| SRP | 取得、編集、契約検証、CLIを分離 | アプリ不具合、素材不良、編集修正、公開確認を別の問題として処理 |
| OCP | 順序・字幕・尺は台本で変更できる | 編集変更だけなら撮影しない。ただし未知の操作まで設定言語で表現しない |
| LSP | 継承を作らず、素材・ヘルパーの契約を保つ | 素材差し替えは出所・検証・寸法の契約を満たすものに限定。既存テストの意味を変えない |
| ISP | rendererにPageや巨大なアプリ状態を渡さない | 各セッションは必要な入力・成果物・証拠だけを受け取る |
| DIP | 上位台本は素材IDを参照。CLIで具体処理を結合 | 内容の企画とセレクター・外部コマンドを切り離す。不要な抽象クラスは作らない |
| KISS | 画像/動画、single/grid2x2、カット、字幕のみに限定 | 1本、1言語、1認定環境、直列の制作。新規出力先による再試行 |
| DRY | 既存fixture・操作・検証を共有し、字幕は1正本 | 出所、環境、台本の所有者を明確化。説明用の文書と実行設定を混同しない |
| YAGNI | プラグイン、動画GUI、音声、増分buildを作らない | 後続動画を先に作らず、今回の素材を再利用可能にして終了 |

原則に従っているという宣言だけでなく、変更点、依存方向、再利用先、追加しなかった機構を各セッションの引き継ぎへ示す。ファイル分割や行数の削減自体を目的にしない。

## 14. 引き継ぎと運用

各担当者は`docs/videos/meet-gbdraw/handoffs/Sn.md`へ人が読める記録を残す。大きなログや動画はbuild領域に置き、そのパスとハッシュを記録する。引き継ぎに未実行のPASSや仮のSHAを書かない。

引き継ぎの必須項目は、ステータス（PASS/PARTIAL/BLOCKED）、ブランチ・HEAD・差分、実装内容、変更ファイル、確定した契約、実行コマンドと結果、成果物、目視範囲、既存回帰確認、未完了事項、次セッションの必要条件、英語のコミットタイトル案と要約とする。詳細なテンプレートは各INSTRUCTION PROMPTに含める。

初回は手動のローカル制作までを必須とする。CIが必要になった場合に限り、S5で明示起動のworkflowを1つ追加して成果物を保存する。自動公開、全PRでのフル動画build、必須チェック化はしない。

## 15. 実装前リスクと対処

| リスク | 対処 |
|---|---|
| 参照版と現在のdevが異なる | S1で実ファイル・API対応を確認。参照SHAを最新扱いしない |
| 古いブラウザーwheelが実行される | 現行準備手順を利用し、実wheelのハッシュを記録 |
| Screencastの実寸やカーソルが想定と違う | 最小の実録画で測定。出力サイズ指定だけで合格させない |
| Web用WOFF2が字幕処理で使えない | 編集用TTF/OTFを制作環境で1書体固定。GUIフォントと別管理 |
| 表示調整に長いGenerateが必要 | 完全な操作録画を証拠へ保存し、本編は検証済みのBefore/Afterを使う。操作時間を4秒に押し込まない |
| D-loopを隠してもブラケットが残る | 入力由来featureと別annotationを混同しない。実演開始時は別D-loop annotationなしを検証 |
| 配色規則が一部CDSに当たらない | 全13個のカテゴリ被覆・排他性・凡例を検査し、未分類を隠して成功扱いしない |
| 画像ごとの自動fitで円の大きさが変わる | P0〜P3を共通尺度・向き・位置でPNG化し、差分の視認性を点検 |
| 最終出力の録画が4秒を超える | 結果確定後から撮り、必要動作を切らず未通過範囲を報告 |
| 既存フローが承認画像を上書きする | 全出力先をrun領域へ明示し、前後の差分を検証 |
| 不安定な描画IDが比較を壊す | 生物学的identityに基づく比較と、根拠のある限定的正規化 |
| 字幕修正で再解析が走る | render境界をテストし、重いモジュールを遅延import |
| メディアは出たが最後の操作が切れている | 長い録画を拒否、成功末尾のみ延長、境界フレームと全編再生を確認 |

## 16. 出典と確認範囲

以下は設計の根拠となる一次資料。リポジトリ資料は参照コミットを固定したURLで示す。実装時は作業ベースの実ファイルを優先して差分を確認する。外部ドキュメントのAPI記述は、導入済みの固定バージョンとS1の実行で再検証する。

- [R1] [gbdraw / AGENTS.md](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/AGENTS.md)。ブランチ、生成wheel、参照画像、showcase品質の規則。
- [R2] [Documentation capture README](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/README.md) と [Tutorial index](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/TUTORIALS/README.md)。既存撮影環境と対象教材。
- [R3] [human_circular.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/flows/human_circular.py) と [gui_first_circular.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/flows/tutorials/gui_first_circular.py)。ヒトミトゲノムの操作・検証。
- [R4] [gui_losatn.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/flows/tutorials/gui_losatn.py)。実ブラウザー検索を行うLambda–DE3フロー。
- [R5] [gui_losatp_groups.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/flows/tutorials/gui_losatp_groups.py) と [BGC Tutorial](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/TUTORIALS/GUI/compare-proteins-losatp.md)。BGC素材の定義と整列。
- [R6] [gui_annotated_chloroplast.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/flows/tutorials/gui_annotated_chloroplast.py)。タバコ葉緑体の完成図と検証。
- [R7] [Feature presentation Tutorial](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/TUTORIALS/GUI/highlight-mitochondrial-features.md)、[gui_feature_highlight.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/flows/tutorials/gui_feature_highlight.py)、[Feature presentation reference](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md)。既存の全CDSラベル・機能別配色定義と表示処理。教材のD-loopブラケットは、本計画の元フィーチャー表示切替には流用しない。
- [R8] [web_capture.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/flows/web_capture.py) と [web_server.py](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/docs/capture/web_server.py)。共通ブラウザー操作とループバック配信。
- [R9] [Web用Interフォントの配置](https://github.com/satoshikawato/gbdraw/tree/4556e04e929a4a85ad28d1833ce7304bd764881c/gbdraw/web/vendor/fonts/inter)。参照版はWOFF2。
- [W1] [Playwright Python — Screencast](https://playwright.dev/python/docs/api/class-screencast)。2026-09-21に公式APIを確認。実録画検証は未実施。
- [W2] [FFmpeg filters documentation](https://ffmpeg.org/ffmpeg-filters.html)。フィルター構文はS1で固定FFmpegのhelpと最小実行により確認する。この文書作成時のページ再取得はtimeoutとなったため、最新版の全構文を確認したとは扱わない。
- [W3] [ffprobe documentation](https://ffmpeg.org/ffprobe.html)。2026-09-21に情報取得・検査機能の公式資料を確認。実行検証は未実施。

- [R10] [固定ヒトミトゲノムGenBank](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/gbdraw/web/tutorial-data/human-mitochondrion/HmmtDNA.gbk)。参照ファイルでD-loop `complement(join(16024..16569,1..576))`を確認。NCBIのWebページは再取得時に確認画面となったため、新規ダウンロード版を確認したとは扱わない。
- [R11] [tutorial-data manifest](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/gbdraw/web/tutorial-data/manifest.json)。入力ID・ファイル・チェックサムの正本。
- [R12] [CDS gene label priority TSV](https://github.com/satoshikawato/gbdraw/blob/4556e04e929a4a85ad28d1833ce7304bd764881c/gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv)。`CDS\tgene`の既存入力。product開始状態は同schemaで明示設定し、実行時に全CDSの適用を検証する。

改訂2.0の追加確認範囲は、固定リポジトリの元D-loop、ラベルTSV、feature-presentation教材・仕様、およびPlaywright公式APIの読取りである。アプリ操作、4状態の生成、録画、動画生成は未実行であり、S1〜S5の実行検証を必要とする。
