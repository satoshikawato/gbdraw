### Refactoring Daily Log — 2026-01-02 (Fri)

このドキュメントは、`gbdraw` のリファクタリング作業ログ（当日分）と、翌日の方針をまとめたものです。

---

### 今日やったこと（網羅）

#### **1) CLI とライブラリ層の分離（library-ization の第一段）**
- **新API追加**: `gbdraw/api/diagram.py`
  - **`assemble_linear_diagram_from_records(...) -> svgwrite.Drawing`**
  - **`assemble_circular_diagram_from_record(...) -> svgwrite.Drawing`**
  - CLI を通さず「パイプラインから直接 `Drawing` を得る」入口を用意。
- **既存フローの分解**:
  - `gbdraw/circular_diagram_components.py` に **`assemble_circular_diagram(...) -> Drawing`** を追加し、保存処理（`save_figure`）から分離。
  - `gbdraw/linear_diagram_components.py` にも同様に **`assemble_linear_diagram(...) -> Drawing`** を追加（後述でさらに分割）。
- **CLI 側は “薄く”**:
  - `gbdraw/linear.py` / `gbdraw/circular.py` は **APIで組み立て → `gbdraw.render.export.save_figure()` で保存**に整理。

#### **2) 0除算の安定性バグ修正（スモークテストで発見）**
- **`gbdraw/analysis/skew.py`**:
  - cumulative skew が全て0のケースで **0/0** になっていたので回避（全ゼロなら factor=0）。
- **`gbdraw/svg/linear_tracks.py`**:
  - GC/Skew の系列が完全にフラット（max_diff=0）なとき **0/0** になっていたので、フラット線として描画するよう回避。

#### **3) linear の “組み立て” を責務分離（巨大ファイル分割）**
- **実装を移動**: `gbdraw/diagrams/linear/*`
  - `gbdraw/diagrams/linear/assemble.py`: assemble/plot の本体
  - `gbdraw/diagrams/linear/builders.py`: group 生成＋配置の呼び出し
  - `gbdraw/diagrams/linear/positioning.py`: translate 等の座標配置
  - `gbdraw/diagrams/linear/precalc.py`: 事前計算（label height / definition width）
  - `gbdraw/diagrams/linear/__init__.py`, `gbdraw/diagrams/__init__.py`: package 整備
- **互換 façade**:
  - `gbdraw/linear_diagram_components.py` は **re-export だけ**にして、既存 import を壊さずに内部構造を整理。

#### **4) typed config の精度向上（`show_labels` の型）**
- **`gbdraw/config/models/canvas.py`**:
  - `CanvasConfig.show_labels` を `object` → **`bool | Literal["all","first","none"]`** に変更。
  - `from_dict()` で `None` を **`False` に正規化**（不定値を排除）。

#### **5) multi-track の“土台”導入（まだ renderer 未接続のデータモデル）**
- **新パッケージ**: `gbdraw/tracks/*`
  - `gbdraw/tracks/spec.py`: `TrackSpec` / `ScalarSpec` / `CircularTrackPlacement` / `LinearTrackPlacement`
  - `gbdraw/tracks/parser.py`: 最小パーサ `parse_track_spec(s, mode=...)`, `parse_track_specs([...], mode=...)`
- **API から触れる入口**:
  - `gbdraw/api/tracks.py` を追加し、`gbdraw.api` から re-export。

#### **6) circular 側で multi-track を “読み取り専用で接続”（show/hide + 一部配置）**
- **circular の track spec 接続**: `gbdraw/circular_diagram_components.py`
  - **axis/ticks/gc_content/gc_skew/legend/labels** の **show/hide** を `track_specs` で上書き可能に。
  - `@r/@w`（または `ri/ro`）で GC/Skew の **中心半径・幅**を上書き（既存 `track_dict`/`track_ratio_factors` の上に “override”）。
- **必要な注入点（小さく追加）**
  - `gbdraw/groups/circular/gc_content.py`: `norm_factor_override`（+ track width override は呼び出し側で）
  - `gbdraw/groups/circular/gc_skew.py`: `norm_factor_override`
  - `gbdraw/groups/circular/ticks.py`: `radius` override
- **API 側の受け渡し**: `gbdraw/api/diagram.py`
  - `track_specs` を `str`/`TrackSpec` で受け、`parse_track_specs(..., mode="circular")` して `assemble_circular_diagram(...)` に渡す。
  - `legend@show=false` は `legend="none"` 相当に寄せる（canvas sizing も含めて一貫させるため）。

#### **7) circular の labels を “トラック” ではなく “レイアウトレイヤー（arena）” として分離**
- **外側（非埋め込み）ラベルを record group から切り離し**
  - 新規: `gbdraw/groups/circular/labels.py` `LabelsGroup`
    - **非埋め込みラベル＋leader line** を `id="labels"` の独立Groupに。
  - `gbdraw/groups/circular/seq_record.py` は **埋め込みラベルのみ**を描くように変更（features と一体）。
  - `gbdraw/groups/circular/__init__.py` で `LabelsGroup` を export。
- **label arena の導入（最小接続）**
  - `gbdraw/labels/placement_circular.py`
    - `prepare_label_list(..., outer_arena=(inner_px, outer_px))` を追加
    - 外側ラベルの “anchor” を **arena の内側**に寄せる
    - 外側配置用の楕円（arc）を **arena の外側**に収まるようスケール（現段階は「既存アルゴリズムを arena にフィット」）
  - `gbdraw/circular_diagram_components.py`
    - `labels@ri/ro`（または `labels@r,w`）から **arena を計算**して `LabelsGroup` に渡す
    - `labels@show=false` は **外側ラベルだけ OFF**（埋め込みは `config_dict/cfg.canvas.show_labels` に従う）

---

### 現状の到達点（“今どこまで出来ているか”）

- **library API が使える**:
  - `gbdraw.api.diagram.*` で `Drawing` を返せる（パイプライン統合の芯ができた）
- **linear の組み立て責務が分割され、互換も維持**
- **typed config の重要箇所が固まってきた**
- **multi-track は “仕様モデル＋パーサ＋circularへの最小接続” まで到達**
- **labels は “トラック”ではなく“layout layer（arena）”として扱う方向が確定**

---

### 明日の方針（やる順・狙い）

#### **A) circular の labels arena を “より自由な配置” に寄せる（画像のような配置へ）**
- **狙い**: 楕円上に置くだけではなく、arena の帯の中で **半径も個別に可変**にして衝突回避できるようにする。
- **やること（候補）**
  - 初期配置（角度順）→ 衝突検出 → **半径方向に逃がす**（arena 内 clamp）→ 収束、の反復
  - leader line の折れ点（`middle_x/middle_y`）の扱いを整理（内側境界に固定 or 可変）

#### **B) “埋め込み” と “外側” の制御を track spec 側に拡張**
- **狙い**: `labels@show=false` を外側だけ/埋め込みだけ/両方、のように分けられるようにする。
- **案**
  - `labels` に `mode=external|embedded|all` を追加する
  - もしくは `embedded_labels` を別 kind として追加する（`TrackKind` 拡張）

#### **C) circular の features/definition もトラック化の土台に乗せる**
- **狙い**: 将来的な “任意トラック配置” のために、features/definition を独立Groupとして制御できるようにする。
- **やること（最小）**
  - circular も `gbdraw/diagrams/circular/*` へ分割（linear と同じ構造）
  - `TrackSpec` の `z`（重なり順）を実際の add 順に反映する方針を固める

---

### メモ（仕様上の注意点）
- 現段階の arena は **「既存の arc 配置をスケールして帯に収める」**段階で、完全な2D自由配置は明日以降。
- `labels@show=false` は **外側（non-embedded）だけ**を消す。埋め込みを消すには `cfg/config_dict` 側の `canvas.show_labels` を落とす（または明日拡張）。


