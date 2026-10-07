# Web GUI 再監査（2026-10-05）の修正と UI の刷新: 実装計画

作成: 2026-10-05（監査セッション gbdraw-6f）。この計画は、`FINDINGS.md` の 59 件の修正と、Owner が採用した UI の「提案」の実装をまとめる。
次のセッションのプロンプト: `NEXT-SESSION-PROMPT.md`。

## 0. 前提と開始時点

- 監査の対象: `dev` `8e5234ac`。監査の所見: `FINDINGS.md`（このディレクトリ）。各領域の詳細: `reports/<area>.md`。probe: `probes/<area>/`。
- 計画を書いた時点の `dev`: `cdd43fcd`（2026-10-05 夕方）。
- 別セッション（override-residuals、gbdraw-51）は、まだ作業中だった。
  - 作業記録: `/home/kawato/gbdraw-baselines/override-residuals-20261005/HANDOFF.md`
  - 開いていた PR: #807（R-5）、#809（OV-30）、#810（OV-19）、#811（R-7）、#812（OV-26/27）
- その別セッションと重なる項目:
  - **TK-05** = OV-29。#805 で修正済みなので、確認だけする。
  - **FL-13** = OV-27。#812 で直る見込み。merge された後に確認する。
  - **UI-04 の一部**（アイコンだけのボタンの名前）= OV-30（#809）。残りだけを直す。
- 行番号は `8e5234ac` のもの。現在の `dev` で、関数名から場所を探し直す。

## 1. Owner の決定（2026-10-05）

| 項目 | 決定 |
|---|---|
| UI-01 | **「提案」を採用する**。案 A（部品ルールを型付き style の `@layer components` に置く）に、`probes/ui01_proposal.css` の内容を部品ルールそのものとして組み込む。比較ページ（Version 2）: https://claude.ai/artifact/2ub5cMWzCBxURVGzXy6ckE |
| UJ-03/FL-01 | live の凡例を Generate に合わせる（Generate と同じ分割と並び順） |
| UJ-09 | Load Session の前に確認を出す。History の扱いは変えない（Owner-delegated の推奨） |
| TK-13 | CLI に合わせる。record に underlay の feature がなければ、Features の行のない stack を受け付ける |
| CI-07c、FL-12 | 未決定。推奨案（今の挙動のまま、UI に説明を足す）で進め、Owner-delegated として記録する |

次の 3 件は**実装しない**。最後の報告で、選択肢と推奨を示して Owner に確認する。

- **UJ-06**（初めて開いたとき、例のデータへの道がない）
  - 推奨: 空の状態に "Load an example" を置き、Gallery の HmmtDNA Session を読み込む。
  - ほかの案: Gallery のカードに "Open in app" を付ける。
- **UI-10**（電話の幅での preview と Editor）
  - 推奨: 電話の幅では、設定と preview をタブで切り替える。
- **UI-11**（fit-to-window がない）
  - 推奨: preview のツールバーに Fit ボタンを置く。

## 2. UI-01「提案」の仕様（G01〜G03 の正）

参考の見本:
- `probes/ui01_transform.py`（案 A の差し替え）
- `probes/ui01_proposal.css`（提案）
- `screenshots/ui-01/raw/*-proposal.png`

見本は「上書きする CSS」だった。実装では、部品ルールとマークアップを直接直す。上書きの CSS を足してはいけない。

### 2.1 仕組み

- `index.html` の 18 個の `@apply` ルールを `<style type="text/tailwindcss">` の `@layer components { … }` に移す。
  - Play は型付きの style しか処理しない。普通の `<style>` に置いた `@apply` は、2025-12-17 の `1ff0c843` から一度も効いていない。
- `.form-input-compact`、`.adv-options .form-input-compact`、`.auto-value-input` は、同じ層の `.form-input` の後ろに置く。
  - こうしないと、compact の入力欄 58 個の文字が欠け、"(auto)" の表示 23 個が隠れる。
- やってはいけない方法:
  - `type` だけを変える: utility が負ける。
  - `:where()` で詳細度を 0 にする: preflight に負ける。
  - 理由は FINDINGS.md の「UI-01 案 A の実装」にある。
- 色は CSS 変数（token）で持つ。値は `ui01_proposal.css` の `:root` と同じ。
  - page `#eef2f6`、surface `#fff`、text `#1e293b`、text-2 `#334155`
  - muted `#64748b`、line `#e2e8f0`、line-strong `#cbd5e1`
  - accent `#2563eb`、accent-soft `#eff6ff`

### 2.2 見た目の規則

- **文字の階層（4 段）**
  - カードの見出し（`.card-header` と、カード直下の `details > summary`）: 14 px、600、text 色、下線なし。
  - 入れ子の disclosure の見出し: 12 px、600、text-2 色。
  - グループ見出し（設定内の `h4`）: 11 px、600、大文字、letter-spacing 0.06em、muted 色。大文字はこの段だけ。
  - 項目名（`.input-label`）: 12 px、500、大文字にしない、text-2 色、line-height 1.3。
  - 補足: 11 px 以上、muted 色、line-height 1.45。新しい class 1 つ（例 `.ui-hint`）にまとめる。
    - 今は `text-[9px]` が 87 か所、`text-[10px]` が 423 か所、slate-400 の文字がある。
    - これらのうち補足の文字だけを置き換える。compact の入力欄と badge は除く。
- **行の高さ**: `.btn` と `.input-label` は `text-base` を使わない。`text-base` は 24 px の line-height も持ち込み、要素の `text-[10px]` では行の高さが変わらないから。`.btn` の line-height は 1.25 にする。
- **アップロード欄**: 1 行、高さ約 32 px（padding 6 px 10 px、1 px の点線、角丸 8 px、地 `#f8fafc`）。
  - ホバー: accent の枠と accent-soft の地。
  - 選択後: `#86efac` の実線と `#f0fdf4` の地。
  - 黄色の地はやめる。アイコンは 15 px、line-height は 18 px。
- **カード**: 設定のカードは白（surface）で、line 色の枠。
  - Input Genomes の `border-l-4 border-l-blue-500` はやめる。
  - ページの地は page 色。
- **Result のカード**: 緑の枠と ring をやめ、1 px の line 色の枠と弱い影にする。
  - 見出しは text 色、15 px、600。緑はチェックのアイコンだけに残す。
- **ヘッダーの操作**（`.app-config-button`）: すべて surface の地と line 色の枠にする。ホバーで accent。Reset Settings の黄色はやめる。
- **出力ボタン**: PNG の `btn-primary` を `btn-secondary` にして、SVG、Interactive SVG、PNG、PDF を同じ見た目にする。
- **checkbox**: 16 px。
- **Retry Generate**: `btn btn-secondary btn-sm` にする。今は `.btn` がない。
- **`<html lang="en">`**（UI-05）。
- **マークアップの変更**（G03）:
  - Custom Track Slots の見出し行: "Applies on Generate" を見出しの下の行へ移し、見出しが切れないようにする。
  - "Reset to Tuckin / Middle / Spreadout" を、小見出し "Reset to preset" とボタン "Tuckin" "Middle" "Spreadout" に分ける。
  - Preset scheme の select は、`h-7` と `p-2` の組み合わせで文字が欠けている。直す。
  - TK-11: 行の gap と legend label の入力欄が列からはみ出している。`width:100%` と `min-width:0` を付ける。
  - UI-13: 設定パネルの中の横スクロール（3 px のはみ出し）と、Generate の上の空き（約 70 px）を直す。
- **対象外**: Tailwind Play の実行時コンパイラ（コンソールに出る警告）。事前コンパイルへの移行は、この計画に含めない。

### 2.3 退行を防ぐテスト（G01 で追加する）

Playwright の spec を足す。HmmtDNA を読み込み、すべての `details` を開き、Label Mode を Out にした状態で、1280 と 390 の幅で確かめる。

1. 部品のスタイルが効いている: `.form-input` の border が 1 px、`.card` の角丸が 0 より大きい。UI-01 の再発を検出する。
2. 見えている input と select のうち、`clientHeight - padding < fontSize*1.05` のものが 0 個（文字が欠けていない）。
3. 見えている `.auto-value-input` の背景が透明。
4. カードの見出し、項目名、補足の computed `font-size` と `text-transform` が、2.2 の段のとおり。

確かめ方の手本: `probes/ui01_typescan.py` と、FINDINGS.md の測定の probe。

## 3. 進め方の規則

- 作業場所: Linux の clone `/home/kawato/gbdraw-work`。`/mnt/c` では作業しない。
  - worktree は `.worktrees/<name>` に作る。
  - port は 4401〜4420 を使う。4301〜4319 は別セッションのもの。
  - 作業記録は `/home/kawato/gbdraw-baselines/gui-fix-20261006/` に置く（HANDOFF.md、agent-common.md、logs/、pr/）。
- 共通の規則は `/home/kawato/gbdraw-baselines/override-residuals-20261005/agent-common.md` を下敷きにする。worktree、PATH shim、テストを先に書くこと、PR の書き方、検証の予算がそこにある。
- **直す前に再現する**: 各 ID を、最新の `dev` で直す前に再現する。再現しなければ、試したことを書いて閉じる（ほかの PR で直っていることがある）。
- **根本から直す**: 同じ種類の次の control や経路でも起きないように直し、その抜けを見つけるテストを付ける。テストは、修正なしで失敗することを 1 回確かめる。
- **PR は 1 つの原因につき 1 つ**: base は `dev`。
  - PR 本文に次を書く: Product impact の分類、Owner の決定または Owner-delegated の選択、監査の ID。
  - 本文の言葉は `node tools/check-pr-language.mjs` で確かめる。
- **見た目が変わる PR には Before/After のスクリーンショットを載せる**（memory: gui-pr-before-after-screenshots）。
  - 画像は orphan branch `pr-screenshots` の `pr/<PR#>/<name>-before.png` と `-after.png` に置く。
  - 撮影の手本は `probes/ui01_screens.py`。
- **merge**: 必須の check が通ったら、`dev` への auto-merge を arm してよい（memory: auto-merge-allowed、ci-merge-practice）。
  - `main` への push と promotion はしない。promotion は Owner が決める。
  - ARCHITECTURE_EXCEPTION が要る変更は、arm せずに報告する。
- **CI**: 監視する watcher は 1 つ。エージェントに CI を待たせない。
- **モデル**: 機械的な作業は Sonnet に任せる（docs、再生成、文言、規則が明確な置き換え）。
  - Opus は同時に 2 つまで。Sonnet は 2〜3 つまで。
  - ブリーフは短くし、必要な抜粋だけを渡す。
- **検証の予算**:
  - 自分が変えたテストと、変えたコードに grep で当たるテストだけを走らせる。
  - functional の Playwright を全部ローカルで走らせない。
- **形式**: `main`（`fe6861f0`）と release 0.13.0 は Session 44 を持つ。Session 45、request schema 9、catalog 5 は `dev` にしかない。
  - 形式を変えるときは、45/9/5 のまま直し、branch が持つ生成物（Gallery Session）を作り直す。
  - 版や schema を変えたら、テストの期待値をすべて同じ PR で直す。main と deploy だけで走る suite も含める（memory: update-test-expectations-on-version-bumps）。
- **新しく見つけた不具合**: 記録して ID を付ける（`GX-01` から）。
  - 軽微なもの（数行、原因が明確、製品の選択が要らない）は、その場で直す。
  - 大きいものや、製品の選択が要るものは、選択肢と推奨を書いて報告する（memory: note-and-fix-minor-bugs）。
- **質問で止まらない**: 推奨案で進め、Owner-delegated として記録する（memory: proceed-with-recommended-options）。
  - ただし、どの決定も扱っていない製品の挙動を新しく選ぶことはしない。
- **docs**: UI の文言や挙動が変わったら、`docs/REFERENCE/web-app.md` などの持ち主の文書を同じ PR で直す。
  - docs のスクリーンショット（`docs/images/`）は、最後の G18 でまとめて作り直す。

## 4. PR の計画

列の意味: 「モデル」は担当のエージェント（O = Opus、S = Sonnet）。「画像」は Before/After のスクリーンショットが要るか。

| PR | 内容（監査の ID） | 主な場所 | モデル | 画像 | 依存 |
|---|---|---|---|---|---|
| G00 | `FINDINGS.md` と `PLAN.md` を `docs/internal/web-gui-audit-20261005/` に入れる（docs だけ） | docs | S | なし | なし |
| G01 | **UI の土台**: 2.1、2.2（マークアップの変更を除く）、2.3 のテスト。UI-01、UI-05、UI-06 | `index.html` の style、出力ボタン、Retry Generate | O | 要 | G00 |
| G02 | **補足の文字の統合**: `.ui-hint` への置き換え。UI-13（横スクロールと空き）、Preset scheme の select | `index.html` 全体 | S（O が見直す） | 要 | G01 |
| G03 | **マークアップの調整**: Custom Track Slots の見出し行、Reset to preset、TK-11 | `index.html` の track の部分 | S | 要 | G02 |
| G04 | **アクセシビリティ**: UI-02（ダイアログの共通部品）、UI-03、UI-04 の残り（+TK-16）、UI-12、UI-09 | `index.html`、`app/ui.js`、各ダイアログ | O | 要 | G03、#809 |
| G05 | **半径の指定と配置の診断**: TK-01（P1）、TK-02。CLI と Web の両方 | `gbdraw/diagrams/circular/radial_layout.py` | O | なし | なし |
| G06 | **Depth**: TK-03、TK-06、TK-07、TK-08、TK-10、TK-14（docs）。TK-05 は確認だけ | `app/app-setup.js`、`app/depth-track-state.js`、`services/session-request.js`、`services/error-normalization.js` | O | 文言のみ | なし |
| G07 | **track slot の editor**: TK-04、TK-09、TK-12、TK-13（Owner の決定）、TK-15 | `app/circular-track-slots.js`、`app/linear-track-slots.js`、`app/track-slot-validation.js` | O | 要 | G03 |
| G08 | **凡例の後処理**: FL-02、FL-03、FL-04、FL-05、FL-06、FL-07 | `app/candidate-render.js`、`app/legend/entry-actions.js`、`app/feature-editor/color-actions.js`、`rule-actions.js` | O | 要 | なし |
| G09 | **live の凡例を Generate に合わせる**: UJ-03/FL-01（Owner の決定）。live と Generate の凡例を比べる spec を付ける | 凡例の投影 | O | 要 | G08 |
| G10 | **入力の状態**: CI-03、CI-02 | `app/app-setup.js`（source watcher、`setLinearSeqPrimaryFile`）、`app/run-info.js` | O | なし | なし |
| G11 | **入力の文言**: CI-01（Python の文と Web の事前検査）、CI-08、UJ-07、UI-08、UJ-08、UI-07 | `gbdraw/crop_genbank.py`、`app/record-discovery.js`、`services/error-normalization.js`、file-uploader | S | 要 | G01 |
| G12 | **比較**: CI-04、CI-06、CI-07（a、b、d、c は説明を足す）、CI-05、UJ-05 | `gbdraw/web_support/error_adapter.py`、`app/comparison-ui.js`、`app/pairwise-match-popup.js`、`gbdraw/io/comparisons.py` | S（CI-05 は O） | 要 | なし |
| G13 | **ラベルの Undo**: UJ-01 | `app/feature-editor/label-actions.js` | O | 要 | 別セッションの label の PR が merge された後 |
| G14 | **mode、Load、Cancel**: UJ-02、UJ-09（確認を出す）、UJ-10。mode を切り替えても検索語が残る問題 | `app/record-display-options.js`、`services/config.js`、`app/run-analysis.js` | O | 要 | なし |
| G15 | **Interactive SVG**: FL-08、FL-09。popup の題を gene と product のどちらにするかそろえる（アプリと同じにする） | `standalone-interactivity-assets.js` | S | 要 | なし |
| G16 | **色、注釈、その他**: FL-10、FL-11（文言を挙動に合わせる）、FL-12（説明を足す）、FL-13（#812 の後に確認）、FL-14 | prefix の検証、presets、rules の UI、`app/annotations/` | S | 要 | G01 |
| G17 | **その他の観察**: 確かめてから直す。Linear で "Auto will show…" が出続ける、Reset alignment のダイアログの位置、footer で切れる popup、1600 px で drawer が preview に重なる、track の行の `semanticMutationAvailable`、region と逆相補の順序の docs | 各所 | S | 要 | G01 |
| G18 | **Gallery と docs の画像の作り直し**: UJ-04（Gallery Session の Output Prefix が `out`）。refresh tool を直して Session を作り直す。Gallery の SVG（G15 の runtime を埋め込む）、`docs/images/` の GUI のスクリーンショット | `tools/refresh_gallery_sessions.py`、gallery、docs/images | S | 該当なし | ほかのすべて |
| G19 | **記録**: FINDINGS の「修正の状況」、SESSION_LOG、CHANGELOG の確認 | docs | S | なし | G18 |

### 4.1 波（同時に走らせるもの）

- **波 0**
  - G00 を出す。
  - 別セッションの PR（#807、#809〜#812）と HANDOFF の状態を確かめる。
  - 監視の watcher を 1 つ起動する。
- **波 1**: G01（O）、G05（O）、G12（S）、G15（S）。
  - G01 は `index.html` の土台なので、最初に入れる。
  - `index.html` を大きく変える PR（G02、G03、G04）は、G01 の後に順番に入れる。
- **波 2**: G02（S）→ G03（S）、G06（O）、G08（O）、G10（O）、G11（S）、G16（S）。
- **波 3**: G04（O）、G07（O）、G09（O）、G13（O）、G14（O）、G17（S）。
- **波 4**: G18 → G19。

Opus は常に 2 つまで。待ち時間は、Sonnet の PR と merge の後始末に使う。

### 4.2 終わりの条件

- 59 件それぞれが、次のどれかになっている:
  - 修正済み（PR 番号）
  - 再現しない（試したこと）
  - Owner の確認待ち（UJ-06、UI-10、UI-11）
- `dev` で Tests と Gallery publication の CI が通っている。
- `main` への promotion はしない。
- 最終報告に書くもの:
  - 修正した PR の一覧
  - Owner-delegated の選択
  - 新しく見つけた不具合（GX-xx）
  - Owner に確認すること（UJ-06、UI-10、UI-11 の選択肢と推奨）
