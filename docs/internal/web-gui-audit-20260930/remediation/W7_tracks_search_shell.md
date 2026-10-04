<!-- Raw design report of workstream W7 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W7 修正提案（DEV `4c89bab1`、read-only）

**前提。** パスはリポジトリ直下からの相対パスで、行番号は `4c89bab1` のものです。分類は PRODUCT_IMPACT_RATCHET の procedural classification に従い、次の略号を使います。
- IEA: IMPLEMENT_EXISTING_AUTHORITY
- PDR: PRODUCT_DECISION_REQUIRED
- ER: EVIDENCE_REQUIRED

**検証の範囲。**
- FE-06 だけは node で `search-core.js` を直接実行して再現しました（scratch の `design/w7/fe06.mjs` と `fe06b.mjs`）。
- ほかの項目は、コードを読んだ結果と監査の evidence によるものです。私自身はブラウザでは実行していません。

## PR 分割案

| PR | 内容 | 前提 |
|---|---|---|
| PR-1 | TR-02 と TR-03。managed depth 行の reconcile owner を作る（`architecture-change`） | IEA。D-TR02 は別 PR に切り出せる |
| PR-2 | TR-06、TR-08、TR-09 | IEA、Ordinary |
| PR-3 | TR-10(a) と TR-11 | IEA、Ordinary |
| PR-4 | TR-10(b)。help-tip を keyboard や tap で開けるようにする | PDR |
| PR-5 | TR-12 | IEA。production scope 外 |
| PR-6 | FE-06 と Gallery refresh | PDR |

---

## TR-02 削除・無効化・移動した depth 行が戻る

### 1. 分類

**中核の不具合は IEA。** 無関係な切り替えで行が再作成・再有効化され、別の index に付け替えられる件です。根拠は次の authority です。
- `gbdraw/web/CLAUDE.md` の「Explicit track slots are authoritative when enabled」「Capability gain does not replace a currently valid selection」「watcher execution is not an invariant mechanism」
- `docs/REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md:23-26` の「An explicit slot list is authoritative. Turning a custom stack off preserves it」
- `gbdraw/web/index.html:3264` の Reset ボタンの title「…this is the only action that regenerates it」
- `index.html:3256` の「Use the saved custom stack for rendering」

**付随する D-TR02 は PDR。** Linear でも Depth ファイルを設定したときに managed 行を追加するか、という LSP の整合の問題です。
- 現状の Circular は、watcher が行を追加しています。
- 現状の Linear は、「Add Depth TSV series」を押したときにしか追加しません。

| 案 | 内容 |
|---|---|
| **A（推奨）** | 両モードとも、論理 series が最初の source を得たときに 1 行追加する。ただし、その index を参照する行（有効・無効を問わない）がない場合に限る。source を失った series の managed 行は除去する |
| B | 自動追加をやめ、行は Reset でのみ作る。Circular の現行挙動と、GUI tutorial `docs/TUTORIALS/GUI/build-a-quantitative-genome-map.md:60-66` の前提が変わる |
| C | モード間の差を残す |

### 2. 根本原因の再確認

監査の記述は正しいです。次の点を補足します。

- (a) watcher `gbdraw/web/js/app/app-setup.js:2185-2206` は依存に `form.suppress_gc` と `form.suppress_skew` を含みます。そのため、無関係な切り替えでも `ensureCircularTrackDepthSlot` が走ります。
- (b) ensure の処理（`circular-track-slots.js:1712-1798`）には、ユーザーの編集を上書きする箇所が 3 つあります。
  - 無効な行を除外して「欠落」とみなす（`:1736-1737`）。そのうえで新しい行を作り、id が衝突すると `depth_2` になる（`:1768-1773`）。
  - 既存の行を `enabled = true` に戻す（`:1745`）。
  - 重複する managed 行を、未 claim の index へ付け替える（`:1756-1766`）。
- (c) Linear の `ensureLinearTrackDepthSlots`（`linear-track-slots.js:1356-1360`）は、managed 行を無条件に `enabled = true; params = { track_index }` にします。
  - 呼び出し元は `addLinearDepthTrack` と `removeLinearDepthTrack` だけで（`app-setup.js:1869-1871, 1944-1946`）、起動条件が Circular と異なります（LSP 違反）。
  - 監査の「legend_label が消える」はコードからの判断です。evidence（`linear-depth-reenable.json` の s1/s2）には、もともと legend_label がありません。
- (d) watcher は `semanticFileWatchersSuppressed` を見ていません。そのため、suppress_gc を戻す Undo でも同じ復活が起こりうると推定します（未再現）。

### 3. 修正案（コード）

**owner は 1 つにする。**
- `app/depth-track-state.js` に純関数 `managedDepthSlotAdditions({slots, gainedTrackIndexes, managedPredicate})` を追加する。既存の depth 行が参照している index（有効・無効を問わない）は除外する。
- 除去側は、既存の `dropInvalidManagedDepthSlots`（`:367-386`）をそのまま再利用する。

**各 editor に commit wrapper を置く。**
- Circular: `reconcileCircularDepthSlots(previousSourced)`。`makeDepthSlotForTrackIndex`（`:1585`）、`applyCircularGeometryShortcuts`、`commitManagedSlotMutation`（`:1631`）を再利用する。
- Linear: `reconcileLinearDepthSlots(previousSourced)`。既存の `defaultSlot('depth',{side:'below'})` と axis の補正を再利用する。
- activeCount には `show_depth` に依存しない値（`circular-track-slots.js:445-448`、Linear は論理幅）を使う。

**明示的な遷移からだけ呼ぶ。**
- `setCircularDepthFile`（`app-setup.js:1803-1819`）
- `removeCircularDepthTrack` の末尾（`:1899-1901`）
- `addLinearDepthTrack` と `removeLinearDepthTrack`（`:1869-1871, :1944-1946`）
- D-TR02 が A なら `setLinearDepthFiles`（`:1821-1843`）にも追加する。

**stack が非有効でも保存 stack に適用する。** 次の 2 条件を両立できるのはこの方式だけです。
- 先に depth を読み込み、あとで stack を有効化したとき、depth 行が出ること（tutorial の前提）。
- 有効化時に reconcile して、ユーザーが削除した行を復活させないこと。

**削除するコード。**
- watcher `app-setup.js:2185-2206` の全体
  - normalize と conservation の同期は、`setCircularTrackSlotsEnabled`（`circular-track-slots.js:1876-1886`）から明示的に呼ぶ。
  - conservation には別の watcher（`:2207-2217`）が既にある。
  - label の同期は `updateDepthTrackLabelFromFile` が既に行う。
- ensure の claim・再有効化・付け替え部分（`circular-track-slots.js:1733-1766`、`linear-track-slots.js:1349-1377`）
- app の export `ensureLinearTrackDepthSlots`（`app-setup.js:4528`）

**追加しないもの。**
- 汎用の slot 同期フレームワーク
- 削除した行を記録する tombstone（永続 schema が増えるため）
- 新しい watcher

### 4. 上位の対応

必要です。
- managed depth 行の reconcile は slot editor が持ち、depth 入力を変える遷移から明示的に呼ぶ。
- `gbdraw/web/CLAUDE.md` の Module ownership 表に 1 行追加する。
- D-TR02 の Decision Pack を用意する。

### 5. テスト

- **unit（`tests/web/circular-track-slots.test.mjs:175-250`）**: ensure を呼んでいる 4 件を新しい API に書き換える。挿入位置と Axis の assert は残す。次の 4 件を新たに追加する。
  - 無効な行は reconcile しても不変。
  - 削除した行は復活しない。
  - gained index には 1 行だけ追加され、既存の無効行がある index には追加されない。
  - 他の行の順序と params は不変。
- **LSP**: 同じ table を Linear editor にも流す。
- **既存テストの書き換え**: `tests/web/depth-track-session.playwright.spec.js:469-553` は ensure の重複除去を固定しているので、新しい contract に合わせて書き換える（OIPC-C08）。
- **browser（`tests/web/gui-audit-regressions.playwright.spec.js`）**:
  - 監査の 3 シナリオで、slot 列（`id:enabled:params`）と、Generate 後の `g[data-gbdraw-slot-id]` が不変であること。
  - Linear で depth 行を無効にしてから Add Depth TSV series を押しても、無効のまま params が保たれること。
  - Undo 後も不変であること。

### 6. 規模・PR

- production 4 ファイル。churn は約 200〜300、net は 0 以下の見込み。
- watcher から明示遷移へ移すのは lifecycle 責務の再定義なので、`architecture-change` label を付けます。Review REQUIRED です。

### 7. 依存・リスク

- TR-03 と同じ PR にします。
- 旧 Session で、非有効の stack に depth 行がないものは、読み込んで stack を有効化しても depth 行が自動追加されません（Reset で生成できます）。残余リスクとして web-app.md に記載します。
- Gallery `tobacco-chloroplast`（stack 有効）が不変であることを、gallery-publication で確認します。

---

## TR-03 アップローダーの Remove で depth 行が残る

### 1. 分類

**IEA。** 根拠は CLAUDE.md の「Reconcile capability loss synchronously」と「same availability predicate」です。track card の Remove（`removeCircularDepthTrack`）は既に行を落としており、同じ意図のアップローダー Remove だけが違う動作をしています。

### 2. 根本原因の再確認

監査の記述は正しいです。
- `setCircularDepthFile(idx,null)` は slot に触りません。
- availability watcher（`:1964-1970`）が `show_depth` を false にします。
- すると slot watcher が `showDepth` で gate され（`:2199`）、drop 分岐（`circular-track-slots.js:1715-1724`）に到達しません。
- さらに `desiredCircularDepthTrackCount`（`:1430`→`:450-454`）が `show_depth` に依存しています。
- 一方、stack が有効なときの描画は slot が決めています（`services/session-request.js:1597-1599`）。

### 3. 修正案（コード）

- TR-02 の reconcile を、`setCircularDepthFile` から無条件に呼びます。ファイルの追加・置換・削除で共通です。
- アップローダー専用の clear handler は追加しません。

### 4. 上位の対応

不要です（TR-02 に含まれます）。

### 5. テスト

- **browser**: stack 有効・depth ありの状態で「Depth TSV tracks」内の Remove を押すと、次のとおりになること。
  - depth 行が消え、row issue が空になる。
  - Generate が成功する。
  - Undo 1 回でファイルと行が戻る。
- **unit**: activeCount 0 のとき、managed 行は削除され、manual 行は `disabled_` と `depth_binding_error` になること。

### 6. 規模・PR

TR-02 の PR に含めます（10 行程度の増分）。

### 7. 依存・リスク

- Linear の File clear は、PD-OI-025 により論理 series を保持します。既存テスト `depth-track-session.playwright.spec.js:563-564` も、clear 後に slotIndexes が変わらないことを期待しています。そのため Linear は Circular と同じ結論になりません。
- Linear で「全 record が null の series に有効な行がある」状態で Generate が成功するかは **ER** です。
- Circular で途中の列を消すと null の穴が残り、count は減りません（compact は末尾だけ、`depth-track-state.js:103-109`）。この扱いも同じ ER に含めます。

---

## TR-06 stack 無効中に Hide GC を解除すると gc_content が無効のまま

### 1. 分類

**IEA。** 根拠は次のとおりです。
- confirm の文言（`circular-track-slots.js:2165`「Hiding … will disable those custom track slots」）と、restore 関数が存在すること。
- canonical な値は `form.suppress_gc` であり、`_suppressed_by_global` はそこから派生する印にすぎない（one canonical visibility value）。

### 2. 根本原因の再確認

監査の記述は正しいです。
- restore の先頭の `if (!…circular_track_slots_enabled) return;`（`:2187`）のせいで、印が残ります。
- 付随する問題（コードからの判断）があります。`applyGlobalSuppressToSlot`（`:2139-2145`）は、既に無効な行にも印を付けます。その結果、ユーザーが無効にした行を複製してから Hide を解除すると、その行が有効化されえます。

### 3. 修正案（コード）

- `applyCircularSuppressControlsToSlots`（`:1246-1260`）を、form から計算する双方向の純関数にします。
  - form が隠す場合: 有効な行を無効にして印を付ける。
  - 隠さず印がある場合: 有効に戻して印を消す。
- `setCircularSuppressControl`（`:2202-2230`）は、form を更新して `normalizeSlotsInPlace()` を呼ぶだけにします。
- 次を削除します。
  - `disableCircularTrackSlotsForSuppress` と `restoreCircularTrackSlotsForSuppress`（`:2172-2200`）
  - `applyGlobalSuppressToSlot` とその呼び出し（`:1908, :1930, :2135`）。どの呼び出しも直後に normalize がある。
- Generate 時の呼び出し（`app/run-analysis.js:2556`）で、旧バグ状態で保存された Session も自己修復します。

### 4. 上位の対応

不要です。

### 5. テスト

- **unit**:
  - stack 無効中に unhide すると、有効になり印が消えること。
  - ユーザーが無効にした行は、Hide/Unhide を往復しても無効のままであること。
- **browser**: 監査の suppress-A の手順で、SVG に gc_content が出ること。

### 6. 規模・PR

1 ファイル。churn は約 60、net は約 −30 です。

### 7. 依存・リスク

低いです。

---

## TR-08 無効な行に別の行の "(auto)" が出る

### 1. 分類

**IEA。** 根拠は 2 つあります。
- "(auto)" の表示は、その行自身の値であるべきです。
- 互換 fallback には release の根拠がありません。`slotIndex` と `slotId` は、同じ commit `49065afe`（0.13.0 に含まれる）で同時に導入されています。CLAUDE.md の persisted-format 規則では、根拠のない reader は保持しません。

### 2. 根本原因の再確認

監査の記述は正しく、次のように精密化できます。
- Python が出力する `slotIndex` は、出力された行（有効な行と auto underlay）の中での index です（`gbdraw/diagrams/circular/assemble.py:426`、`gbdraw/diagrams/linear/assemble.py:531`）。
- 一方、UI は全行の中での index を渡します（`circular-track-slots.js:2437-2445`、`linear-track-slots.js:1524-1531`）。
- 無効な行は出力に存在しないため、`track-slot-display.js:77-78` の index fallback が別の行に当たります。
- 無効な行と有効な行で id が重複した場合も誤表示になります。validator は有効な行どうしの重複しか検査しません。

### 3. 修正案（コード）

- `findTrackSlotGeometry`（`:53-79`）は slotId での一致だけにし、index fallback と `slotIndex` 引数を削除します。
- editor 側では、無効な行の resolved 値を引かず、既存の estimate 経路を使います。
  - Circular は `circularTrackSlotEffectiveEnabled` で判定する（`:2474-2478`）。
  - Linear は `enabled===false` で判定する（`:1555-1559`）。
- 未使用の `slotGeometryInstanceKey`（`:50-51`）を削除します。

### 4. 上位の対応

不要です。

### 5. テスト

- `tests/web/track-slot-display.test.mjs:98-118` の「older geometry without slot IDs」は、根拠がないので削除します。
- 代わりに次を追加します。
  - 無効な行、または同じ id の無効な行に対して null が返ること。
  - 両モードで、ticks を無効にして Generate した後、ticks 行に gc_content の値が出ないこと。

### 6. 規模・PR

3 ファイル。churn は約 40 です。

### 7. 依存・リスク

- Generate 後に行の id を改名すると、表示が estimate に戻ります（実態どおりの表示です）。
- PD-OI-048/050 の Auto の意味は変わりません。

---

## TR-09 Show Coordinate Scale OFF でも ticks 行が追加される

### 1. 分類

**preset Reset で showTicks が抜けている件は IEA。**
- 同じファイルの `resetCircularTrackSlotsFromSimpleControls`（`:1688-1710`）と、`services/session-request.js:1648-1656` は showTicks を渡しています。
- preset Reset も depth・GC・skew は反映しているので、showTicks の欠落だけが偶発です。

**初回の stack 有効化で ticks 行が入る件はバグではありません。**
- `setCircularTrackSlotsEnabled` は、配列が空のときしか stack を再生成しません（`:1880-1885`）。初回は、保存済みの既定 stack（`services/session-active-config-contract.js:42`）がそのまま使われます。
- これは UI 契約どおりです。
  - 「Use the saved custom stack」（`:3256`）
  - 「only action that regenerates it」（`:3264`）
  - note `index.html:4478`「Use an enabled Ticks slot to control coordinate-scale visibility」
  - 各 title の「Reset to copy this value into the stack」（例 `:4067`）

### 2. 根本原因の再確認

前半（preset Reset）は監査の記述どおりで、原因は `:1851-1858` です。

### 3. 修正案（コード）

- `resetCircularTrackSlotsFromSimpleControls` の本体を削除し、`resetCircularTrackSlotsToPreset(state.form.track_type)` に委譲します。2 つの関数の差は、preset と showTicks だけです。
- preset 版に `showTicks: state.form.show_scale !== false` を追加します。
- 副次効果として、Reset にも sessionBusy の判定が付きます。

### 4. 上位の対応

不要です。
- 初回有効化を simple controls から作り直したい場合は、別に PDR が必要です。
- 選択肢は A「未編集の既定 stack なら再生成する」と B「現状維持」で、B を推奨します。

### 5. テスト

**unit**:
- show_scale=false のとき、3 つの preset Reset のいずれでも ticks 行が入らないこと。
- Reset と「Reset to <現在の preset>」の結果が同じ slot 列になること。

### 6. 規模・PR

churn は約 30、net は約 −20 です。PR-2 に含めます。

### 7. 依存・リスク

低いです。

---

## TR-10 アクセシビリティ

### 1. 分類

**(a) 名前のない control は IEA。**
- HTML-AAM と WCAG 4.1.2 / 2.5.3（可視ラベルを名前にする）という決定的な仕様規則があります。
- 技術 owner の `docs/REFERENCE/web-app.md:714-726` も、主要な control には安定した名前があると約束しています。

**(b) help-tip が hover 専用である件は PDR。**
- `gbdraw/web/js/components.js:12-16` が「id のない tip は hover 専用の icon のままにし、周囲の label の名前を変えない」を意図として明記しています。
- したがって変更は、既存の設計を退役させる判断になります。

### 2. 根本原因の再確認

**(a) の構造的な原因。**
- `<label class="input-label">X <help-tip/></label>` が control の兄弟で、`for` を持っていません。
- help-tip 177 個のうち 147 個が label の中にあります。そのうち control を包む label は 25、`for` を持つ label は 8 で、残りの約 114 はどの control も label していません。

**監査の列挙は不完全です。** evidence（`aria-unnamed-*.json`）には 26 件あり、監査に挙がった 8 種のほかに次も名前がありません。
- depth の色 input（`:1585`）
- 色規則の input（`:3815, :3853, :3886`）
- PRESET SCHEMES の select（`:3758` 付近）
- overlay layer の select（`:3060, :3443`）
- Linear の Axis Gap（`:4052`）

**モード間の不一致もあります。**
- Linear の行（`:2996-2998`）は、title の 'Enabled'/'Disabled'・'Slot id'・'Renderer' だけが名前になっています。そのため、状態によって名前が変わり、行どうしも区別できません。
- Circular の行（`:3339-3352`）には、slot id を含む aria-label があります（LSP 違反）。
- Dinucleotide（`:4535`）は placeholder の "GC" しか名前がありません。

**(b)** `index.html:7429-7430` の v-else 分岐は、focus できない `<i aria-hidden>` に mouseenter しか付いていません。

### 3. 修正案（コード）

**(a)** 変更は `index.html` だけです。
- 既存の慣例どおり、可視ラベルと同じ `aria-label` を付けます。
  - `Window`、`Step`、`GC Content Mode`、`GC Height`、`Center Reserved Radius`、`Dinucleotide`
  - `` `Feature lane ${id}` ``、`` `Tick label layout ${id}` ``、`` `Depth track index ${id}` ``
- Linear の行には、Circular と同じ形の aria-label を付けます。既存テストが `getByTitle('Enabled')` を使っているので、title は残します。
- `for`/`id` の付与は、label 内の tip が名前に混ざるため、(b) の決定後に行います。

**(b)** 選択肢は次のとおりです。

| 案 | 内容 | 代償 |
|---|---|---|
| **A（推奨）** | 全 tip を disclosure button にする。id を自動採番し、v-else 分岐を削除し、tip を label の外へ移し、対象の control から aria-describedby で参照する | tab stop が最大 175 増える |
| B | tap での開閉だけ追加する | keyboard 利用者は開けないまま |
| C | 現状維持。残余リスクを文書化し、必須の説明は id 付きの tip に限る | ― |

A を推奨する理由です。
- 390 px の表示は PD-OI-038 で支援対象になっている。
- PD-OI-054 で常時表示の説明を tip へ移したため、tip に届くことの重要性が増した。

### 4. 上位の対応

- `gbdraw/web/CLAUDE.md` の「Changing a setting」手順 2 に、次の 1 文を追加します。「可視ラベルと一致する accessible name を付ける。placeholder や状態で変わる title は名前にしない」。
- (b) の決定後に、components.js の意図コメントと web-app.md の Accessibility を更新します。

### 5. テスト

- 新しい spec `tests/web/accessible-names.playwright.spec.js` を追加します。
  - 両モードで入力を読み込み、stack を有効にし、全 `<details>` と Custom Track Slots を開く。
  - 見えている `input`/`select`/`textarea` のすべてについて、`toHaveAccessibleName(/\S/)` を満たすこと。
  - 名前の出所が placeholder だけではないこと。
  - checkbox の名前が checked の切り替えで変わらないこと。
- (b) が A なら、次も検査します。
  - label の中に `.help-tip-trigger` がないこと。
  - すべての tip が focus できること。
- 静的な文字列テスト（`tests/test_gui_tracks_capture_contracts.py:158-176` の型）では網羅できないため、使いません。

### 6. 規模・PR

- (a): 1 ファイル、churn は約 60〜90。
- (b) を A で行う場合: churn は約 500〜800、net は +100〜200 で Size review REQUIRED。パネル単位での分割を推奨します。

### 7. 依存・リスク

- Gallery の capture が影響を受ける場合は、撮り直します。
- (b) を採ると tab の順番が長くなります。

---

## TR-11 stack 行の slot id 欄と renderer 欄が 31 px

### 1. 分類

**IEA。** control・順序・名前をすべて保ったままの reflow なので、アフォーダンスは変わりません。overflow menu にする案を取る場合は PDR になります。

### 2. 根本原因の再確認

監査の記述は正しいです。
- `grid-cols-[auto_minmax(0,1fr)_minmax(0,1fr)_auto]` の 4 列目に、24 px のボタンが 6〜8 個入ります（`index.html:3338-3367`）。
- 幅 280 px のパネルでは、id 欄と renderer 欄が 31〜44 px になります。
- Linear の行（`:2995-3008`）も同じ構造です（監査では計測されていません）。

### 3. 修正案（コード）

- `<style>` の `.track-slot-geometry-grid`（`:164-168`）の隣に、次の 2 つを追加します。
  - `.track-slot-row-head`: 3 列（checkbox / id / renderer）
  - `.track-slot-row-actions`: `grid-column:1/-1`、flex-wrap、右寄せ
- 両モードの行で同じ class を使います。JS や menu は追加しません。

### 4. 上位の対応

不要です。

### 5. テスト

幅 1280・1920・390 px のそれぞれで、両モードについて次を確認します。
- id 欄と renderer 欄の幅が 96 px 以上であること。
- 横方向のあふれがないこと。
- ボタンの数が変わらないこと。

### 6. 規模・PR

churn は約 30 です。PR-3 に含めます。

### 7. 依存・リスク

- stack パネルの tutorial media は再撮影が必要になりえます（`tools/capture_gallery_tutorial_screenshots.py`）。

---

## TR-12 Gallery tutorial のリンク切れ

### 1. 分類

**IEA。** 文書の訂正で、product outcome は変わりません。

### 2. 根本原因の再確認

監査の記述は正しく、次を補足します。
- **生成器は所有していません。**
  - `tools/prepare_interactive_gallery_assets.py:1457-1460` は、`./tutorials/<id>.json` への参照を書くだけです。
  - `tools/refresh_gallery_sessions.py` は tutorial に触れません。
  - `artifact-manifest.json` も tutorial を含みません。
  - 本文は手書きです（owner は `gallery/`、編集手順は `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md`）。
- 対象のファイルは `2badff49`（「reorganize documentation」）で削除されました。
- **相対パスの形式自体が無効です。** `../../../docs/...` は `/gallery/` から見てホストの root の外を指します。ほかの tutorial は `https://github.com/satoshikawato/gbdraw/blob/main/docs/...` を使っています（例: `lambda_basic_linear.json:459`）。

### 3. 修正案（コード）

- `gbdraw/web/gallery/tutorials/vibrio-harveyi-group-collinear.json:576-579` を次のように置き換えます。
  - リンク先: `https://github.com/satoshikawato/gbdraw/blob/main/docs/REFERENCE/web-app.md#record-selection-and-layout`（origin/main の `:106` に存在することを確認済み）
  - label: 「Read how Linear rows and record order work」
- 代わりのリンク先の候補は `docs/TUTORIALS/GUI/compare-proteins-losatp-collinear.md` です。
- Gallery の再生成は不要です。

### 4. 上位の対応

link checker のテストを追加します（下記）。

### 5. テスト

`tests/test_web_packaging.py` に `test_gallery_tutorial_links_resolve` を追加します（隣に `:624-656` があります）。全 tutorial JSON の `href`/`src` を再帰的に集め、次を検査します。
- 相対リンク: `gallery/` を基準に解決した結果が `gbdraw/web/` の中にあり、ファイルが存在すること。
- `./#id`: examples.json の id であること。
- `blob/main/<path>#anchor`: 作業木に path があり、anchor が見出しの slug と一致すること。

### 6. 規模・PR

JSON 1 行とテスト約 40 行です。

### 7. 依存・リスク

- blob/main のリンクは、main に昇格するまで 404 になりえます。

---

## FE-06 既定の All 検索で遺伝子名が配列に一致する

### 1. 分類

**PDR。**
- All が配列を含む挙動は 0.12.0 から release 済みです（`ff14c108`、`c55812db`。IUPAC 展開は `94b694f6`）。
- `docs/REFERENCE/web-app.md:651-664` は、All の構成を定めていません。
- PD-OI-054 が維持を求めているのは「全 field」、つまり field の一覧であり、All の中身ではありません。

| 案 | 内容 |
|---|---|
| **A（推奨）** | All から配列の内容を除外する。対象は Nucleotide sequence、Amino acid sequence、`/translation` qualifier の値。配列の検索は、専用の Nucleotide / Amino acid field（IUPAC）で行う |
| B | 配列を All に残し、literal 一致にする |
| C | 現状維持し、文書化する |

**node での実測。**
- `CYTB` は nt の `CCTC` に一致します。
- `dnaA` は `TCAA` や `GCAA` に一致し、ほぼすべての feature に当たります。
- **B では直りません。** All は `/translation` qualifier を literal に照合するため（`search-core.js:392-394`）、`gyrA` は `…GYRA…` を含む蛋白に一致しました。
- アミノ酸の文字で綴れる 4 文字の遺伝子名は、細菌の proteome で偶然の一致が数件出ると見込まれます。
- HmmtDNA では、literal の `GATC` が nt で 21/76 の feature に一致します。これは本物のモチーフの一致です。

### 2. 根本原因の再確認

監査の記述を精密化します。
- **原因は IUPAC 展開だけではありません。** All が配列の内容を含むこと自体が原因で、IUPAC 展開はそれを増幅しているにすぎません。
- **Interactive SVG は「同じ実装」ではありません。** 手で移植した重複です。
  - IUPAC の処理: `services/standalone-interactivity-assets.js:2206-2310`
  - All の組み立て: `:2650-2666`
  - matcher: `:2669-2742`
- Python も `gbdraw/render/interactive_svg.py:1245-1257` で同じ template literal を静的に抽出するため、CLI で出力する Interactive SVG も同じ挙動です。

### 3. 修正案（A の場合）

- `search-core.js:396-397` を削除し、`:392-394` のループで `translation` を除外します。
- `standalone-interactivity-assets.js:2663-2664` と `:2661-2662` にも同じ変更をします。
- IUPAC 経路は、専用 field のために残します。

**DRY の収束手段。** 完全な統合はできません。
- JS の build step は禁止されている。
- 埋め込み runtime は module import なしで SVG の中で動く。
- Python が literal を静的に抽出する。

そのため、次の 2 点で収束させます。
- 同じ PR で両方の実装を変更する。
- 同じ fixture で parity test を行う。既存の `embeddedFunctionSource` の型（`tests/web/orthogroups-stable-identity.test.mjs:154-170, 1616-1650`）を使う。

**追加しないもの。** query の見た目から配列かどうかを推定する heuristic、新しい field。

### 4. 上位の対応

- web-app.md に、All が検索する範囲を 1 文追記します。
- parity test を、二重実装を同時に変更することの機械的な保証とします。

### 5. テスト

- **unit**（`orthogroups-stable-identity.test.mjs:1245-1261` を拡張）:
  - All で `CYTB`、`dnaA`、`gyrA` を検索すると、名前が一致する feature だけが返ること。
  - `ATGGNN`（nucleotide）と `MPEPTIDE`（amino-acid）は従来どおり一致すること。
- **parity**: 埋め込み版と module 版で、一致する id の集合が等しいこと。

### 6. 規模・PR

- production 2 ファイル、churn は約 20 です。docs、CHANGELOG、テストも更新します。
- 埋め込み runtime のバイトが変わるため、Gallery の `gallery/examples/*.svg`（10 件すべてに埋め込まれていることを確認済み）を `tools/refresh_gallery_sessions.py` で再生成し、manifest も更新します。gallery-publication の証拠が必要です。
- `docs/images/{t-cli-11,h-cli-13,t-py-08}/*.interactive.svg` も旧 runtime を含みますが、テストはバイトを比較していないので、更新は任意です。

### 7. 依存・リスク

- 決定が出るまで、依存する実装は止まります。
- FE-11（Location 検索）が同じ `search-core.js` を触るので、担当の W と順序を調整します。

---

## 横断的な提案

- **watcher で slot を修復しない。** 入力を変える遷移が owner の reconcile を明示的に呼ぶ形にする（TR-02、TR-03、TR-06 に共通）。
  - 「depth の reconcile を watch callback から呼ばない」を unit test で固定する。
  - architecture rule にする場合は、checker と authority の分離規則に従い別 PR にする。
- **1 つの概念には 1 つの builder。** simple controls から stack を作る処理を 1 つにし（TR-09）、suppress の印は form の純関数にする（TR-06）。
- **Circular と Linear の editor に、table-driven の parity test を置く。** depth 行の削除・無効化・移動、無関係な切り替え、ファイルの追加・削除、Add series を同じ表で流す。
- **表示値は identity（slotId）で引く。** index space の異なる位置 fallback は持たない。release の根拠がない互換 reader は削除する。
- **accessible name の browser test を 1 本置く。** 全パネルを開いた状態で両モードを検査する。あわせて `gbdraw/web/CLAUDE.md` に命名規則を 1 文追加する。
- **help-tip の到達性は Decision Pack で一括して決める。** 決定までは、新しい tip は id 付き（focusable）で追加する運用にする。
- **Gallery tutorial JSON と `docs/GALLERY.md` に link checker を置く。** ホストの root 制約と、blob/main の path の存在を検査する。
- **standalone runtime を変更する PR には Gallery refresh を同梱する。** これを skill と PR テンプレートに明記する。
- **search-core と standalone の二重実装は、関数単位の parity test を増やす。** 既存の `embeddedFunctionSource` を拡張する。

## 誤分類・非バグと考える指摘

- **TR-09 の後半（初回の stack 有効化で ticks 行が入る）** はバグではありません。UI 契約どおりです。
- **TR-10 の help-tip が hover 専用である件** は、意図された設計です（`components.js:12-16`）。バグではなく PDR として扱うべきです。一方、名前のない control の列挙は不完全で、evidence には 26 件あります。
- **TR-02 の「Linear で legend_label が消える」** はコードからの判断です。evidence では観測されていません。
- **FE-06 の「Interactive SVG も同じ実装」** は正確ではありません。手で移植した重複です。また「IUPAC の文字だけでできた遺伝子名」という原因の説明も不正確です。実際の原因は All が配列と translation を含むことで、`GATC` の一致は本物です。
- **TR-12 の「生成器を特定する」** について、生成器は存在しません（手書き）。相対パスの形式自体も、ホスト上では無効です。
- **TR-03** は、Linear の同等操作が PD-OI-025 によって別の結論になります。Linear 側は ER です。
- **TR-11** は Linear の行にも同じ構造があり、影響範囲が過小評価されています。

**相互参照。** TR-04 は W3 X-02、TR-05 は W3 X-01、TR-07 は W3、TR-01 は W6 が担当です。
