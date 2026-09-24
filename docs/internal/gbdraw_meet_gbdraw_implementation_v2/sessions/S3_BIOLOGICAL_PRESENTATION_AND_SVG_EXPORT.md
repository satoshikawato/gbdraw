# INSTRUCTION PROMPT — S3: ラベル選択・CDS機能別配色・D-loop領域注釈ブラケットとSVG出力

文書版: 2.2<br>
改訂日: 2026-09-24

## 1. セッションの目的とゴール

`satoshikawato/gbdraw`の紹介動画「Meet gbdraw」に必要な、**生物学的に意味のある3つの表示調整（P0→P1→P2→P3）と最終SVG出力のGUI録画**を実装・検証します。

対象は固定ヒトミトゲノム（`NC_012920.1`、16,569 bp）です。
1. **P0 → P1**: CDSのproductラベルをgeneラベルへ切り替える（Labelsパネル）
2. **P1 → P2**: 全13 CDSを4機能カテゴリ（呼吸鎖複合体別）に配色する（Colorsパネル）
3. **P2 → P3**: 公式チュートリアルに準拠し、Region Annotationsパネルで `mitochondrial_regions.tsv` を読み込み、Custom Track Slotsで内側トラックに原点越えのD-loopブラケットを注釈表示する
4. **SVG出力録画**: P3の状態から実際にSVGダウンロードを完了する操作を録画する

本セッションのゴールは、**研究者がゲノム図を論文品質へ整える価値を実証すること**です。3つの操作は実GUIで最後まで実行・録画して検証証拠として保持し、本編で使用するBefore/After静止画（各4秒、2秒×2）とその元SVG、および本編用の最終SVG出力録画（4秒・120フレーム基準）を確定させ、S2の4素材と合わせて計9素材を揃えます。

## 2. 前提条件と依存関係

- [総合実装計画書](../01_MASTER_PLAN.md)（文書版2.2、特に§4、§6、§8、§12）
- `AGENTS.md`、`CLAUDE.md`
- 先行セッションS1/S2の実装、成果物、および引き継ぎ記録（`S1.md`, `S2.md`）
- 既存の教材・設定正本:
  - `gbdraw/web/tutorial-data/manifest.json`、固定`HmmtDNA.gbk`
  - ラベル優先度TSV: `gbdraw/web/tutorial-data/shared/cds_gene_qualifier_priority.tsv`
  - 機能別配色定義: `docs/TUTORIALS/GUI/highlight-mitochondrial-features.md` 内のCDS配色ルール（4群）
  - 領域注釈TSV: `gbdraw/web/tutorial-data/human-mitochondrion/mitochondrial_regions.tsv`（公式チュートリアル正本、D-loop: 16024..576, wraps_origin=true）
  - 既存フロー・検証: `docs/capture/flows/human_circular.py`、`assertions/gui_feature_highlight.py`

> [!IMPORTANT]
> S1/S2の成果物を保持した作業ツリーから継続してください。`main`や`dev`への直接コミット・プッシュは行わず、無断のresetや古い参照コミットへの強制チェックアウトは避けてください。

## 3. 実装する4状態（P0〜P3）のデータ契約

同じ入力データ（`NC_012920.1`、16,569 bp）を用い、同一アプリセッション内で順に状態を遷移させます。

| 状態 | CDSラベル | CDS配色 | D-loopブラケット | 描画構成（ベース＋カスタムスロット） | 前状態からの変更意味 |
|---|---|---|:---:|---|---|
| **P0** | 全13個でproduct | 単一基底色 | なし | ベース1トラック（37論理フィーチャー: CDS 13, rRNA 2, tRNA 22） | 実演開始用の準備状態（Before 1） |
| **P1** | 全13個でgene | P0と同じ単一色 | なし | ベース1トラック（37論理フィーチャー） | ラベルの採用qualifierのみ変更（After 1 / Before 2） |
| **P2** | P1と同じgene | 4機能カテゴリ配色 | なし | ベース1トラック（37論理フィーチャー） | CDSの色と対応凡例のみ変更（After 2 / Before 3: `dloop-before`） |
| **P3** | P2と同じgene | P2と同じ4色 | **内側トラックに表示** | **ベース1トラック（37個）＋内側カスタムトラック1スロット（D-loopブラケット）** | 領域注釈による原点越えD-loopの可視化（After 3: `dloop-bracket` / Outro） |

> [!NOTE]
> すべての状態で元GenBankレコードのsequenceおよび全注釈情報は不変です。P0〜P3は仕様上の論理状態名であり、汎用ステートマシンを実装する必要はありません。

## 4. 推奨実装手順

### Step 1: ソースデータの確認とスナップショット
1. `HmmtDNA.gbk` をパースし、全長16,569 bp、CDS=13, rRNA=2, tRNA=22（合計ベースフィーチャー=37）であることを確認する。
2. `mitochondrial_regions.tsv` を確認し、D-loopの区間が `16024..576`（原点越え、1,122 bp）として正しく定義されていることを確認する。
3. ソース配列と全注釈のハッシュおよびスナップショットを記録する。

### Step 2: 読める開始図（P0）の作成
1. `HmmtDNA.gbk` を読み込み、CDS・rRNA・tRNAの37論理フィーチャーを表示する。
2. 全CDSを単一基底色にし、ラベル優先順位を `CDS\tproduct` に設定して全13 CDSにproductラベルを表示する。
3. P0のSVGをダウンロードし、全13 CDSの文字列がproductであることを検証して、映像用PNGを生成する。

### Step 3: P0 → P1（ラベルをproductからgeneへ切替）
1. 操作録画を開始する。
2. Web UIの **Labels** パネルにて、既存の `cds_gene_qualifier_priority.tsv`（`CDS\tgene`）を **Priority File (TSV)** として読み込み、**Generate Diagram** を実行する。
3. 結果確認後に録画を停止し、証拠として保存する。
4. P1の実出力SVGをダウンロードし、全13 CDSが `ND1`, `COX1`, `ATP6`, `CYTB` などのgene名に更新されていることを検証する。色や表示フィーチャー集合が変わっていないことを確認する。

### Step 4: P1 → P2（全13 CDSを4機能カテゴリで配色）
1. 操作録画を開始する。
2. 既存公式チュートリアルの配色定義からCDSの4群（完全一致ルール）を適用する:
   - **NADH dehydrogenase** (7個): `ND1`, `ND2`, `ND3`, `ND4`, `ND4L`, `ND5`, `ND6`
   - **Cytochrome c oxidase** (3個): `COX1`, `COX2`, `COX3`
   - **ATP synthase** (2個): `ATP6`, `ATP8`
   - **Cytochrome b** (1個): `CYTB`
3. Web UIの **Colors** パネルにて、上記ルール（TSV）を **Specific Table (-t)** として読み込み、**Generate Diagram** を実行する。
4. 結果確認後に録画を停止し、証拠として保存する。
5. P2の実出力SVGをダウンロードし、13個のCDSが期待通りの4色に塗り分けられ、対応する凡例が正しく生成されていることを検証する。

### Step 5: P2 → P3（Region AnnotationsパネルによるD-loopブラケット追加）
1. 操作録画を開始する。
2. **UI操作手順**:
   - Web UIの **Region Annotations** パネルにて、公式チュートリアルの正本 `mitochondrial_regions.tsv` をインポートする。
   - **Custom Track Slots** 設定にて、ベーストラックの内側に新規トラック枠（inner track slot）を追加し、インポートしたD-loop領域注釈をこのスロットに割り当てる。
   - **Generate Diagram** を実行する。
3. 結果確認後に録画を停止し、証拠として保存する。
4. P3の実出力SVGをダウンロードし、ベーストラックの37フィーチャーの内側に、16024〜576の原点越え区間（1,122 bp）を示す「D-loop」弧状ブラケットが正しく描画されていることを検証する。
5. **QA検証**: 内側スロットを削除または非表示に戻してP2が再現されること（P3→P2→P3の往復再現性）を確認する。

### Step 6: 最終SVG出力録画の取得
1. P3が画面に表示されている状態から、SVGダウンロード操作の録画を開始する。
2. 事前にSVGダウンロードボタンが見える位置にスクロール・配置しておき、無駄な待機時間を排除する。
3. SVGボタンをクリックし、実ダウンロードの完了を検出して録画を停止する。
4. **尺の基準と正規化**:
   - 30 fps正規化後120フレーム（4秒）を基準とする。
   - 実録画が120フレームより短い場合は、最後の静止フレームを複製して120フレームに延長する。
   - わずかにオーバーした場合（数フレーム程度）は、ダウンロード完了後の不要な余白を末尾トリミングまたは微細なPTS正規化を行って安全に120フレームに揃える。
5. ダウンロードされた最終SVGがP3の期待値と完全一致することを検証する。

### Step 7: 共通構図でのPNG化と素材登録
1. P0〜P3の各SVGから映像用PNGを生成する。
2. **共通構図（フリーズ誤認・ガタつき防止）**:
   - 4枚の円の中心座標、半径、向きを一致させ、ラベル変化・配色変化・内側ブラケット追加による差分のみが視覚的に浮き彫りになるよう、4状態の最大boundsを包含する共通キャンバス・余白を設定する。
3. `assets.json` に以下の5素材を追加し、S2の4素材と合わせて**全9素材**を完成させる:
   - `human.labels-product` (image, P0由来)
   - `human.labels-gene` (image, P1由来、ラベルAfter & 配色Before兼用)
   - `human.functional-colors` (image, P2由来、配色After & D-loop追加前 `dloop-before` 兼用)
   - `human.dloop-bracket` (image, P3由来、D-loop After `dloop-bracket` & Outro兼用)
   - `human.svg-export` (video, 最終SVG出力録画)

## 5. 既知の落とし穴とガードレール

- **D-loopの正しい実装**: 本機能は公式チュートリアル `docs/TUTORIALS/GUI/highlight-mitochondrial-features.md` に準拠した **Region Annotations と Custom Track Slots** による内側弧状ブラケットの追加です。非推奨のnative feature toggleや手動描画の偽装を行わないでください。
- **偽装の禁止**: 13個のラベルを1個ずつ手入力で書き換えたり、全CDSを1色で塗りつぶす等の手抜きを行わない。正規のTSV適用とGenerateを実行する。
- **分類の排他性と完全性**: 正規表現照合において `ND4` と `ND4L` の完全一致を厳密に行い、誤った重複分類や未分類CDSを残さない。
- **録画尺の混同**: 120フレーム（4秒）の制限は本編用の「最終SVG出力録画」だけに適用されます。3つの編集操作の証拠録画に4秒制限を課して操作を切断してはいけません。
- **原点越えの検証**: D-loopは16024から16569を経て1から576へと環状原点を跨ぎます。直線化されたり分割されて2つの独立フィーチャーと誤認されないよう、`wraps_origin=true` の描画パス構造を検証してください。

## 6. 受け入れ条件（Acceptance Criteria）

1. **AC-04準拠**: P0/P1/P2/P3 の全状態で元レコードの配列・注釈情報が不変である。
2. **AC-05準拠**: 13 CDSすべての product→gene 切り替えが確認できる。
3. **AC-06準拠**: 4機能カテゴリの7/3/2/1被覆と排他性、凡例の対応が検証される。
4. **AC-07準拠**: Region AnnotationsパネルおよびCustom Track Slotsによる内側D-loopブラケット（16024..576、原点越え）の描画、およびP3→P2→P3復元往復再現性が検証される。
5. P3の最終ダウンロードSVGと出力録画・本編P3静止画が同一ハッシュで結合される。
6. 本編出力録画が120フレーム（4秒）に正しく正規化される。
7. 4状態のPNGが共通の中心・尺度で生成され、各Before/Afterの差分が視覚的に明瞭である。
8. S2の4素材と合わせて計9素材の完全なマニフェストが整う。

## 7. 引き継ぎフォーマット

セッション完了後、`docs/videos/meet-gbdraw/handoffs/S3.md` に以下を実値で記録してください。

```text
Status: PASS / PARTIAL / BLOCKED
Implementation base: branch, HEAD, preserved S1/S2 work
Source: input/wheel hashes, full sequence and annotation identity
States: P0/P1/P2/P3 definitions and actual SVG/PNG paths/hashes
Labels: actual priority UI/schema, all-CDS product/gene mapping
Colors: reused rule owner, CDS-only derivation, 7/3/2/1 coverage, legend evidence
D-loop bracket: Region Annotations TSV, Custom Track Slots configuration, origin-spanning proof, round trip
UI actions: actual public paths and Apply/Generate requirements
Recordings: complete edit evidence versus final export clip, real sizes/durations
Final export: actual download, P3 SVG/PNG identity, 120-frame limit
Assets: five S3 IDs plus four S2 IDs, complete nine-asset manifest
Tests: exact commands, results, injected failures, skipped/unverified scope
Visual review: all four states, common framing, raw recordings actually inspected
Regression: affected helper callers and protected paths
Next session: exact asset bundle and remaining S4 prerequisites
Proposed commit title: English
Proposed commit summary: English
```
