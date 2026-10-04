<!-- Raw design report of workstream W4 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W4 修正提案: History の区切り、Session の読み込み、Generate と Worker の寿命

## 結論
- **SE-06 は監査が示した修正方法が誤っていた。** `services/imported-comparison-intent.js:191` を緩めると EDITABLE になる。その状態で Generate すると、保存された BLAST ではなく Web 既定の LOSAT protein 比較が実行された。match 要素は 8715 → 112 に減り、科学的出力が黙って置き換わる。**読み取り専用のままが正しい。** 本当の不具合は次の 2 つ。
  - record の識別子が CLI と Web で食い違う。
  - Inherit のとき、schema 7 の比較を schema 8 に昇格せずに流用する。
- **SE-08 の Qualifier Priority の無視は tobacco に限らない。** Linear の CLI Session はすべて該当する。
- **GE-06 は Product Decision が必要**（案 A を推奨）。
- **SE-08 の読み込み遅延は、証拠を用意したうえで判断が必要。**
- それ以外は既存の authority どおりに直せる。

## 検証方法
- 対象は dev `4c89bab1`。dev worktree は読み取りだけで、変更していない。
- scratch の `/tmp/claude-1000/-mnt-c-Users-genom-GitHub-gbdraw/c9538259-b3ce-4057-baa8-e783da1943ee/scratchpad/design/w4/` で実験した。
  - `site/gbdraw/web/js` は JS を複製し、使い捨ての prototype patch を当てたもの。
  - 実験スクリプトは `p01`–`p09`、`k01`–`k04`（Python Playwright）。
- **prototype で確かめたこと**
  - dev では失敗する SE-01/02/03/04、GE-06、GE-09、SE-06 の各シナリオが、prototype では期待どおりになった。
  - 既存テストも prototype で通った。
    - node: `history`、`history-inputs`、`history-canonical-owner`、`history-config-restore`、`generation-feedback`、`diagram-resource-recovery`、`session-cli-compatibility`、`imported-comparison-intent`
    - Playwright: History 系 10 件と、`session-cli-compatibility` の linear case
- **確かめていないこと**
  - SE-07 は仕組みを確認しただけで、修正を試していない。
  - SE-08 の Python 修正も試していない。

---

## SE-01（P2）Undo/Redo 後の catalog が admit されない

**1. 分類:** IMPLEMENT_EXISTING_AUTHORITY
- 根拠は `docs/REFERENCE/web-app.md:706-707`（Undo/Redo は form と editor の変更をたどる）と、`gbdraw/web/CLAUDE.md` の History 規則。
- 例外が出て Features が 0 件になるのは、どう見ても誤り。

**2. 根本原因の再確認**
- 監査の指摘は正しい。
- 補足すると、checkpoint の復元は trusted=false で行われる。そのため `normalizeEditorStateData`（`gbdraw/web/js/services/config.js:1076-1080`）と `applyEditorStateData`（`:1169-1171`）で二重に複製される。
- **同じ原因の潜在経路がもう 1 つある**（コードからの判断）。Session 読み込みの rollback は、`config.js:3579` で adopted catalog を参照のまま取得する。しかし `:3636` の復元で複製され、admit されない catalog になる。
- `services/history-snapshot.js:1334-1344`（export は `:1588`）の `applyGeneratedArtifactSnapshot` は、gbdraw と tests のどこからも呼ばれていない。
- 性能面でも損がある。checkpoint ごとに、catalog 全体を JSON で複製して署名している。HmmtDNA の色変更 1 件の checkpoint は 1,397,698 → 1,067,930 bytes（-24%）に減った（prototype で計測）。

**3. 修正案**
- 役割の分担
  - 取得は `services/history-snapshot.js` が持つ。
  - catalog を state に入れる窓口は `config.js` の `applyEditorStateData` 1 か所に絞る。
- 手順
  1. `buildArtifactCheckpoint`（`:1514-1529`）
     - catalog の参照は、`captureGeneratedArtifactOwnerSet`（`:710`）と同じ `artifactOwnedValue(getGeneratedArtifactRef(state.featureCatalog, null))` で取る。
     - 本体の JSON からは `editorState.featureCatalog` を除く。`buildEditorStateData({ preserveAdoptedCatalog: true })` を使い、複製せずに取り除く。
     - 参照は、checkpoint object を key とするモジュール内の WeakMap に保持する。History は `record.checkpoint` の同一性を保つので、この方法で足りる。署名（`:1584`）は変えなくてよい。
  2. `applyArtifactCheckpoint`（`:1570` の後）で、参照を `setGeneratedArtifactRef(state.featureCatalog, ref)` で戻す。`installGeneratedArtifactOwnerSet`（`:800`）と同じやり方。
  3. `config.js` の `applyEditorStateData` と `normalizeEditorStateData` では catalog を複製しない。`isAdoptedFeatureCatalog` が真のものだけを state に入れ、それ以外は null にする。
     - `adoptCatalog` option と複製の分岐（`:1169-1171`）、`:4646` の引数を削除する。
     - これで import rollback の経路も同時に直る。
  4. 使われていない `applyGeneratedArtifactSnapshot` を削除する。
- **採らない方法:** Undo のたびに `admitFeatureCatalog` で admit し直すこと。catalog 全体を毎回たどるうえ、二重の複製が残り、復元の途中で失敗する道も増える。
- **prototype の結果**（`p04_se01.py`）: Undo、Redo、Undo の後に Linear → Circular で 37 features のまま、エラーなし、adopted は常に true。dev は adopted=false になり、例外が出た。

**4. 上位の対応:** `gbdraw/web/CLAUDE.md` の History 段落に 1 文を足す。「History の履歴項目と復元経路は、Generate が所有する artifact（feature catalog）を参照で持つ。`state.featureCatalog` は null か admit 済みの catalog だけ」。

**5. テスト**
- node の `tests/web/history.test.mjs`（`catalogA` を使う既存の仕組み、`~1735-2099`）
  - checkpoint → undo → redo の後、`state.featureCatalog.value === catalogA` であること。
  - checkpoint の JSON に catalog が入っていないこと。
- `applyEditorStateData` に admit されていない object を渡すと null になること。
- Playwright（`history-generated-authority` を拡張）
  - Gallery の HmmtDNA で、凡例項目が増える色変更と Reset Settings を行う。
  - それぞれ Undo/Redo してから Linear → Circular に切り替え、pageerror がなく 37 features であること。
  - 性能の証拠として `getDiagnostics().checkpointEstimatedBytes` が減ること。

**6. 規模・PR:** production 2 ファイル、churn 約 60、net は負。Ordinary、単独の PR。

**7. 依存・リスク**
- 古い catalog の参照が checkpoint に残り、その分のメモリが byteSize に数えられない（許容できる）。
- catalog をその場で書き換えるコードはない（grep で確認）。

---

## SE-02（P2）ラベルクリックとキーボード操作が Undo できない

**1. 分類:** IMPLEMENT_EXISTING_AUTHORITY（`web-app.md:706-707`）

**2. 根本原因の再確認**
- 監査より範囲が広い。キーボードの Space（checkbox）と ArrowRight（radio）も、dev で 0 step だった（`k01`）。
- 原因は、checkbox/radio の transaction を pointerdown（`gbdraw/web/js/app/history-inputs.js:108-112`）でしか始めないこと。ラベルクリック、Space、矢印キーはどれも、v-model が値を反映した後の bubble の `change`（`:131-137`）にしか届かない。
- `k02` での確認: root の capture-phase の `click` と `change` の時点では、reactive state はまだ変更前だった（v-model の listener は input 要素上の target phase で動くため）。

**3. 修正案**
- 担当は `app/history-inputs.js`。
- capture-phase の `change` listener を追加し、checkbox/radio なら `beginForElement(target)` を呼ぶ。commit は既存の bubble 側の `onChange` のままにする。
- pointerdown の分岐（`:110`）から checkbox と radio を外す。
- label を control に対応づける処理や、keydown 用の特例は足さない。
- **prototype の結果**（`k04`、`p01`）: ラベルクリック、Space、Arrow のそれぞれで 1 step が記録され、Undo の順序も正しかった。

**4. 上位の対応:** 規則として明文化する。「discrete control の transaction は、値を確定するイベントの capture phase で始める（checkbox/radio は `change`、button は `click`）。text 系は focus で始める」。

**5. テスト（入力イベントの組み合わせ）**
- `tests/web/history-inputs.playwright.spec.js` に追加する組み合わせ:
  - 対象: checkbox、radio、select、text、mode button
  - 操作: control 本体のクリック、ラベル文字のクリック、キーボード（Space、Arrow、Enter）
  - 前提: フォーカスなし、text 欄にフォーカス、text 欄に入力済み
- 確認すること:
  - 変更した control ごとに Undo step がちょうど 1 つ増える。
  - Undo で変更前の値に戻り、戻る順序が後入れ先出しになっている。
- node の `history-inputs.test.mjs` で、capture の change から bubble の change までが 1 step になることも確認する。

**6. 規模・PR:** PR-C（SE-02/03/04 をまとめる）

**7. 依存・リスク**
- ユーザー操作で発生したイベントでは、listener と listener の間に microtask が処理される（HTML 仕様）。修正はこれに依存している。
- script の `dispatchEvent` ではこの処理が起きないため、unit test では模擬が必要。

---

## SE-03（P2）text 欄にフォーカスがあるときのクリックが落ちる・まとめられる

**1. 分類:** IMPLEMENT_EXISTING_AUTHORITY（同上）

**2. 根本原因の再確認**
- 監査の指摘は正しい。`services/history.js:368` の `begin()` は、誰が開いたかに関係なく、開いている transaction を返す。
- drag も同じ型の潜在不具合を持つ（コードからの判断）。`app/legend/drag-actions.js:87` と `app/legend-layout/diagram-drag.js:88` は source だけを渡して `begin` するため、text 欄がフォーカスを保っていると同じことが起きうる。

**3. 修正案（History の区切りを管理する `services/history.js` が担当）**
1. `begin(label, { source, owner })` に所有者の概念を入れる。
   - 別の owner の transaction が開いていれば、先にそれを確定してから新しく始める。
   - 同じ owner か owner 未指定（`runUndoable` の合流）なら、今の transaction を返す。
2. 「開いている intent を先に確定する」処理が 3 か所に重複している。これを `settlePendingIntent()` 1 つにまとめ、重複を削除する。
   - `beginCheckpoint`（`:507-516`）
   - `beginArtifactReplacement`（`:672-681`）
   - `runUndoableCommand`（`:918-920`）
3. `history-inputs.js`
   - `:61` で `owner: element` を渡す。
   - button の transaction は click の capture（`:166-173`）だけで始める。pointerdown 側の分岐（`:96-102`）は削除する。こうすると、前の欄の @blur による値の正規化（例: `index.html:2274`）が、その欄の step に入る。
4. drag の 2 ファイルでも owner を渡す。
- `pendingCommit`（`:48, 58, 82-87`）は不要になるかもしれない。組み合わせテストが通った場合だけ削除する。
- **prototype の結果**（`p01`）
  - text 入力 → checkbox は 2 step に分かれた（dev では 1 step にまとめられた）。
  - text 入力 → Linear も 2 step に分かれ、Undo すると先に Circular に戻った（dev では Linear のままだった）。

**4. 上位の対応:** `gbdraw/web/CLAUDE.md` に不変条件として書く。「transaction は 1 つの owner（control または gesture）に属する。別の owner の transaction を始めるときは、開いているものを先に確定する」。

**5. テスト**
- SE-02 と同じ組み合わせテスト。
- node: owner A の後に owner B を始めると、A が確定して B が開くこと。owner なしの `runUndoable` は合流すること。

**6. 規模・PR:** PR-C

**7. 依存・リスク**
- Undo の step が増える（意図どおりで、ユーザーから見える変化）。
- Safari ではクリックしてもフォーカスが移らないが、owner の規則でこの場合も区切られる。

---

## SE-04（P3）select にフォーカスがあると Ctrl+Z などが効かない

**1. 分類:** IMPLEMENT_EXISTING_AUTHORITY
- 根拠は `web-app.md:706`（Undo/Redo は form の変更をたどる）と、CLAUDE.md の「有効な操作が状態を変えずに黙って戻ってはならない」。
- select には、守るべきブラウザ標準の undo がない。
- 監査のいう range は `index.html` に存在しない（0 件）。

**2. 根本原因の再確認**
- 原因は監査のとおり（`app/history-shortcuts.js:7, 18`）。
- **修正すると表に出る問題が 1 つある。** ショートカットが効くようになると、Ctrl の keydown で adapter（`history-inputs.js:123-129`）が select の新しい transaction を開く。この transaction は「変更後の値」を前状態として持つので、Undo の後に focusout すると Undo が新しい「Change setting」として記録され、Redo が消える。

**3. 修正案**
- `history-shortcuts.js:7` で select を text 系から外す。
- `history-inputs.js` の `onKeyDown` では、Ctrl/Meta/Alt を押したままのキーを無視する。
- `history.js` の `undo()` と `redo()` の先頭で `settlePendingIntent()` を呼ぶ。
- **prototype の結果**（`p02`）: dev は無反応だった。prototype では Undo と Redo が効き、Undo してから blur しても redo=1 のままだった。

**4. 上位の対応:** 不要

**5. テスト**
- 組み合わせテストに行を追加する: select にフォーカスして ArrowDown、Ctrl+Z、Ctrl+Shift+Z、Ctrl+Y、blur の順に操作し、回数と値を確認する。
- text 欄では Ctrl+Z がブラウザ標準の undo のままであることも確認する。

**6. 規模・PR:** PR-C（3 ファイルに drag の 2 ファイルを加えて計 5 ファイル。churn 約 100、net は +10〜30。Ordinary）

**7. 依存・リスク:** 低い。

---

## GE-06（P2）Generate 中の Ctrl+Z

**1. 分類:** PRODUCT_DECISION_REQUIRED（現在の挙動は NOT_ALLOWED）
- `web-app.md:105-106` と `:130-137` は、Generate 中に使えないものとして Save/Load しか挙げていない。Generate 中の Undo をどうするかは決まっていない。
- 一方、今の挙動は CLAUDE.md の History 規則（transaction 中は before が current のままでなければならない）と、PD-OI-016（UI が状況を示す）に反する。
- 選択肢
  - **A（推奨）:** History 自身の artifact replacement か checkpoint の transaction が開いている間は、Undo/Redo を busy として拒否する。ボタンにも同じ判定を使い、Generate はそのまま確定させる。overlay にはすでに理由が表示されている。
  - **B:** Undo/Redo で先に Cancel を実行し、直前の Result を保ったうえで Undo する。
- **別に判断が必要な問い**（横断的な提案 8 を参照）: overlay は inert ではないので、キーボードなら後ろの設定を編集できる。Generate 中の編集全般を止めるかどうか。

**2. 根本原因の再確認**
- 監査は `state.js:793-802` を原因としたが、これはずれている。Generate 中の編集は仕様上許されており、`'mutation'` が processing を見ないのは意図どおり。
- 破られているのは History 側の不変条件。
  - `undo` と `redo`（`services/history.js:1000-1049`）は、開いている artifact transaction を確認しない。
  - `runUndoableArtifactReplacement`（`:827-877`）には、実行中を示す印がない。
  - `activeCheckpoint`（`:575-593`）も、`undo` からは参照されない。
- dev で再現: Rendering 中の Ctrl+Z で draft の編集が戻り、Generate C は古い draft を確定した（`p03`: draft 5000、確定済み 1000）。

**3. 修正案（A の場合）**
- `history.js` に `activeReplacement` を足す。
- 判定は `historyAvailability()` 1 つにまとめる: `mutationAvailability() ||`（`activeCheckpoint` か `activeReplacement` があれば busy）。
- この判定を `undo`、`redo`、`canUndo`、`canRedo`（`:1052-1053`）で共有する。
- **reactive にする必要がある。** 開始時と解除時に `touchTransaction()` を呼び、判定の中で `transactionRevision` を読む。prototype では、これがないと Undo ボタンが無効のまま戻らなかった。
- 文言は state.js の既存 busy 文言を使うか、W3（X-01）に任せる。
- shortcut 側に Generate 専用の分岐は足さない。
- **prototype の結果:** 拒否され、C は確定し、undo の数が 1 増えた。

**4. 上位の対応:** 決定後に、`web-app.md:130-137` へ「Generate 中は Undo/Redo も使えない」と追記する。

**5. テスト**
- node: `runUndoableArtifactReplacement` が保留中なら `undo()` が busy を返し、`canUndo` が false、保留が解けると true に戻ること。
- Playwright: `__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse`（`services/diagram-generation.js:361-372`）で Generate を止めて Ctrl+Z/Ctrl+Y を押す。止めている間は確定済み request と undo/redo の数が変わらず、再開後に C が確定して履歴が 1 つ増えること。

**6. 規模・PR:** `history.js` に約 20 行。決定後に PR-C へまとめてよい。

**7. 依存・リスク:** A では、Generate 内部で LOSAT が長く動く間もキーボードの Undo が止まる（意図どおり）。

---

## GE-09（P3）JS 側の準備中の Cancel で Worker が終了する

**1. 分類:** IMPLEMENT_EXISTING_AUTHORITY（CLAUDE.md の lazy Worker の規則「後の操作は Worker を再利用する」）

**2. 根本原因の再確認**
- 監査の指摘は正しい。`services/diagram-generation.js:605-616` は `hadActiveRequest` を計算しているのに、`:614` で無条件に終了させる。
- Cancel の後は、`runAnalysisInternal` の `throwIfGenerationCanceled`（`app/run-analysis.js:2012-2016`）が、`executeCanonicalCandidate`（`:4680`）より前に働く。そのため、終了させなくても後から Worker に要求が送られることはない。

**3. 修正案:** 終了処理の前に `if (!hadActiveRequest) return false;` を入れる。要求ごとの所有者の区別は足さない（YAGNI）。

**4. 上位の対応:** 不要

**5. テスト**
- `tests/web/generation-feedback.test.mjs` にケースを追加する: Worker が起動済みで要求がないとき、false が返り、Worker は終了せず、次の実行で同じ instance を使う。
- このテストは dev で失敗し（terminated=true）、prototype で通ることを確認した。既存の 3 件と 9 件も通った。

**6. 規模・PR:** 1 ファイル 2 行、単独の PR-A。

**7. 依存・リスク:** なし。ほかの helper が動いている最中の Cancel が、それも終了させる点は今と同じ。

---

## SE-06（P2）CLI の Linear+BLAST Session で Generate できない

**1. 分類:** Inherit の不具合は IMPLEMENT_EXISTING_AUTHORITY
- 根拠は PD-OI-008 の「実行可能だが投影できない比較は、明示的な inheritance で使い続けられる」。
- **読み取り専用（PRESERVED_READ_ONLY）の扱いは正しい。**
- EDITABLE にするのは新しい挙動になるので、望むなら PRODUCT_DECISION が必要。今はやらないことを推奨する。

**2. 根本原因の再確認（監査は一部誤り）**
- **監査は `imported-comparison-intent.js:191` を不具合としたが、これは誤り。**
  - ここを緩めると EDITABLE になり、Generate は Web 既定の LOSAT を実行した。
  - comparisons は `[nucleotideBlast]` から `[precomputedProteinComparison, orthogroupResult]` に変わり、match は 8715 → 112 になった。
  - これは PD-OI-008 の「比較を黙って置き換えない」に反する。
  - CLI の sidecar は `webFiles.bindings.linearComparisons: []` を明示的に書き、Web は CLI 由来の Session に比較の draft を作らない（`docs/SESSION_COMPATIBILITY.md:162-163`）。
- **本当の不具合 A（record の識別子）**
  - CLI は bindings の uid に `cli-seq-N` を書く（`gbdraw/cli_utils/session.py:1321`）。一方、同じ Session の request では `record-N` や `record-N:k` を使う。複数 record のファイルでは `record-1:1`、`record-1:2`、`record-2` になることを確認した。
  - `services/session-request.js:2906-2925` は binding の uid をそのまま draft の recordKey にする（`:910`）。
  - そのため Inherit の対応づけ（`imported-comparison-intent.js:437-447`）で例外になる。
- **本当の不具合 B（新しく見つけたもの）**
  - Inherit は、確定済みの comparisons を schema 8 に昇格しないまま候補にコピーする（`:449`）。
  - main で作った sidecar（v42、schema 7）には `settings.alignOrthogroupFeature` が残っており、Python が拒否して VALIDATION_UNCLASSIFIED になる。`origin/main` の `4556e04e` の CLI 出力で再現した。
- **互換性の証拠**
  - `cli-seq-N` は 0.12.0 タグと main に存在する。
  - 0.13.0 の CLI は v30 を書くため、Web では読めない。
  - main の CLI は v42 を書き、Web で読める。形も同じで、nucleotideBlast と無効な generatedProteinComparison を持ち、uid は `cli-seq-N`。
  - したがって、**Web の読み込み側を直す必要がある。** CLI の書き出し側は変えない。書き出しを変えると uid の形が 2 種類になるが、読み込み側は `cli-seq` を今後も受け付け続けなければならない。生成済みの成果物も書き換えることになる。

**3. 修正案**
1. `session-request.js` の CLI sidecar 用の分岐（`:3785-3799`、条件は `initializeCliInputs && storedConfig == null`）。Linear のときは、各 binding の uid を確定済みのファイル単位の recordKey に置き換える。
   - ファイル単位の key は、recordKey の末尾の `:<n>` を除いた値を、順序を保って重複なく並べたもの。
   - 置き換えるのは、key の数が bindings の数と一致するときだけ。
   - Web の draft はこの分岐を通らず、自分の uid を保つ。
2. `app/run-analysis.js:4590-4593` では、`promoteCanonicalRenderRequestToCurrent(committed.renderRequest, { featureCatalog })` を通して渡す（`session-request.js:4784`。`:4826` で legacy の欄を削除する）。昇格の経路をこれ 1 つにする。
3. 比較を選ばずに Generate したときの文言は W3（X-01）が担当する。
- **prototype の結果**
  - main の v42 と dev の v44 の sidecar のどちらも、Inherit → Generate が ok になった。
  - dev の v44 では SVG が byte 単位で一致した。main の v42 では match 数が 8715 のままだった。

**4. 上位の対応**
- `docs/REFERENCE/session-and-request-compatibility.md:50-53` に 1 文を足す。「CLI sidecar の record の識別子は `renderRequest.records[].recordKey`。binding の uid は初期値にすぎない」。
- 互換性の規則に従い、main で作った小さな v42 の CLI Linear+BLAST sidecar と、その来歴を fixture として固定する（`single.v41-bindings1` と同じ方式）。

**5. テスト（CLI → Web の往復）**
- `tests/web/session-cli-compatibility.playwright.spec.js` の cases に `linear blast` を追加する。小さな入力と数行の outfmt 6 を使う。
- 確認すること:
  - 読み込み後の扱いが PRESERVED_READ_ONLY であること（LOSAT への置き換えを防ぐ）。
  - draft の uid がファイル単位の key と一致すること。
  - Inherit → Generate が ok で、comparisons が確定済みのものと一致し（昇格の分を除く）、svgSemantics も一致すること。
- 固定した v42 の fixture も同じ手順で確認する。
- unit では、schema 7 の pipeline を Inherit できることを確認する。

**6. 規模・PR**
- production 2 ファイル、約 30 行。fixture、spec、文書を 1 文追加。
- 互換経路の変更なので、Review REQUIRED。

**7. 依存・リスク:** `session-request.js` の CLI 分岐を SE-07 と共有するので、順に進めるか 1 つの PR にまとめる。

---

## SE-07（P2）CLI Session で凡例の位置が失われる

**1. 分類:** IMPLEMENT_EXISTING_AUTHORITY
- 根拠は `session-and-request-compatibility.md:50-52`（Web は設定を renderRequest から初期化する）、PD-OI-019（保存された layout の選択が優先される）、OIPC-C05。

**2. 根本原因の再確認（監査は正しい。`k03` で仕組みを確認）**
- 投影された form では、`legend`（`session-request.js:4410`）が `multi_record_canvas`（`:4412`）より前に並ぶ。
- `safeDeepMerge`（`config.js:1902`）は、アクセサ（`state.js:235-245`）を通して legend を書く。この時点では既定の MRC=true なので、値は circular.multi の枠に入る。
- 次に MRC が false になる。
- 最後に `restoreLayoutPreferences(ui, { preserveActive: true })`（`config.js:1845-1884`、呼び出しは `:4554`）が、merge 後の値 'left' から全部の枠を作り直す。その結果、upper_left はどの枠にも残らない（実測: single も multi も 'left'）。
- 同じ種類の例: `rendered-v27.v40-schema6` の fixture。保存された form が projection を上書きするため、`none` が `bottom` になる。
- Linear の通常の CLI Session は影響を受けない。

**3. 修正案（`ui.layoutPreferences` の表現を持つ `session-request.js` が担当）**
1. projection は、確定済みの (mode, grouping) の layoutPreferences を返す。作り方は `app/layout-preferences.js` の `createDefaultLayoutPreferences` と `updateActiveLayoutPreference` を使い、値は `options.output.legend` と `plotTitlePosition` から取る。
   - アクセサ経由の `legend`（`:4410`）と `plot_title_position`（`:4497`）は projection から削除する。
2. `config.js` の `restoreLayoutPreferences` は、保存された値がないときは projection の値を使う。
   - `activeBeforeRestore` と `preserveActive` の分岐（`:1845-1880`）は削除する。
3. `config.js:4405` の `legendSide` は、保存済み Result の側を表すので、draft ではなく確定済みの値を使う。
- merge の順序で解決する細工や、新しい session の欄は入れない。

**4. 上位の対応:** 不要。現在の所有者の規定をそのまま実装するだけ。

**5. テスト**
- CLI 往復の cases に次を追加する。
  - circular `--legend upper_left`
  - circular `--multi_record_canvas --legend upper_right`
  - linear `--legend left`
- 確認すること: 読み込み後に `form.legend === request.output.legend` であり、最初の Generate の後も request と SVG の凡例位置が変わらない。
- projection の unit test も足す。rendered-v27 の期待値の変更は、レビューのうえで受け入れる。

**6. 規模・PR:** 2 ファイル、churn 約 60。Review REQUIRED。

**7. 依存・リスク**
- Gallery の publication（`services/gallery-session-publication.js:40-61`）と `CURRENT_WRITER_FORM_FIELDS`（`services/session-active-config-contract.js:65`）が影響を受ける。
- Gallery を refresh したとき、Session が一致したままか確認が必要。

---

## SE-08（P2）

### (a) Qualifier Priority の編集が無視される

**1. 分類:** IMPLEMENT_EXISTING_AUTHORITY（PD-OI-009 の「editor のない leaf だけを保持し、管理される隣の欄は守る」、OIPC-C03）

**2. 根本原因の再確認（監査より範囲が広い）**
- CLI の encoder は、処理済みの filtering をそのまま `diagramOptions.config` に書く（`gbdraw/session_request_codec.py:2238-2243`、`gbdraw/circular.py:1051-1054`、`gbdraw/linear.py:1428`）。
  - この中には、`preprocess_label_filtering`（`gbdraw/labels/filtering.py:406-438`）が作る派生値 `whitelist_map`、`priority_map`、`label_override_rules` が含まれる。
- Web で保持する設定を導く処理（`gbdraw/web_support/config_overrides.py:101-116`）は、これらを unmanaged とみなし、`labels.filtering.raw` として残す。
- この raw が、Generate のたびに request に入る（`session-request.js:1272-1280`）。
- Python の `_resolve_diagram_options_config`（`gbdraw/api/diagram.py:283-291`）は、新しい table を付けても古い `priority_map` と `whitelist_map` を残す。消すのは `label_override_rules` だけ。
- その結果、「処理済み」と判定されて table の処理が飛ばされる（`filtering.py:412-414`）。label override の table がなければ、Qualifier Priority と Whitelist の両方が無視される。
- **影響範囲**
  - dev と main の CLI が書いた Linear の Session はすべて該当する。`lin1` で 'CDS: gene' に変えても反映されないことを確認した。
  - Circular では、table を使った Session（tobacco など）が該当する。通常の Circular CLI Session は影響を受けない（確認済み）。

**3. 修正案（Python のみ）**
- `config_overrides.py` では、派生した key を、比較の対象からも投影する raw からも除く。
  - 残りの変更がすべて GUI で管理される欄なら、raw は捨てる（既存の `:114-115` の規則）。
  - 保存済みの明示的な raw（`requireUnmanagedOnly` の経路）にも同じ処理をかける。main の Web には同じコードがあり、main で保存された Session のための互換 fixture が必要。
- 派生した key の集合は、`labels/filtering.py` で 1 か所に定義する。
- 追加で防御も入れる: `diagram.py:283-291` で、table を付けるときに古い map を捨てる（`:290-291` と同じ扱い）。
- CLI の書き出し側は、この PR では変えない。

**4. 上位の対応:** 不要

**5. テスト**
- pytest
  - tobacco の config と通常の CLI Linear の config は、どちらも保持する設定が空（`{}`）になる。
  - 派生した key を含む明示的な raw は、取り除かれる。
  - 本当に unmanaged な filtering は保持される。
  - 古い `priority_map` がある config に table を付けて描くと、table のとおりのラベルになる。
- ブラウザ: CLI Linear の Session で Priority を変えるとラベルが変わり、preserved の欄に raw が出ない。

**6. 規模・PR:** Python のみ（Web の production の規模には入らない）。browser wheel は生成される。

**7. 依存・リスク:** 低い。

### (b) 読み込みに約 10 秒かかる

**1. 分類:** EVIDENCE_REQUIRED の後に PRODUCT_DECISION_REQUIRED
- 「lazy Worker の規則に反する」とは言えない。preview は Python を要しないが、PD-OI-009 の型付き検証は、config があると Python を必要とする（`config.js:1106-1136` が早期に抜けるのは config が null のときだけ）。
- 影響はすべての CLI sidecar に及ぶ。tobacco は全体 9.9 秒で、そのうち Worker の起動が 9.2 秒、helper は 20 ミリ秒。
- 選択肢
  - **A:** 読み込み時の検証を続ける。
  - **B:** 最初に Python を使う操作まで検証を遅らせる。Load の一括性（PD-OI-045）と、設定の開示のタイミング（PD-OI-009）が変わる。
  - **C:** CLI が full config ではなく、差分だけの configOverrides を書く。replay で結果が変わらないことの証明が必要で、既存のファイルには効かない。
- **推奨:** 既存のファイルには A を適用する。C は、損失のない replay を証明したうえで、別の architecture PR として進め、Gallery の tobacco を refresh する。B は推奨しない。

---

## PR の分け方と順序
PR-A（GE-09）→ PR-B（SE-01）→ PR-C（SE-02/03/04）→ PR-D（GE-06、決定後。PR-C に含めてもよい）→ PR-E（SE-06）→ PR-F（SE-07）→ PR-G（SE-08a、Python）。SE-08b と GE-06 は Decision Pack にする。

## 横断的な提案
1. History に不変条件を入れる。transaction は 1 つの owner に属し、別の owner を始めるときは先に確定する。`settlePendingIntent` を `begin`、`beginCheckpoint`、`beginArtifactReplacement`、`runUndoableCommand`、`undo`、`redo` で共有し、重複を削除する。
2. discrete control は、値を確定するイベントの capture phase で始める。操作の手段に左右されない。
3. History の利用可否は、reactive な判定 1 つにまとめる。`undo`、`redo`、`canUndo`、`canRedo` はこれを共有する。
4. Generate が所有する artifact（catalog など）は、どの履歴項目でも復元経路でも参照で持つ。JSON で複製したり署名したりしない。
5. CLI sidecar の初期化は、`session-request.js` の 1 つの分岐が担当する。確定済みの request を正とし、uid や layout はそこから決める。
6. CLI → Web の往復の組み合わせテストを、受け入れの基準にする。最初の Generate の request が確定済みのものと等しいこと（凡例、comparisons、record key、filtering）を確認し、main で作った Session を固定の証拠として置く。
7. 処理の途中で作られる派生値は、設定として永続化の境界を越えさせない（定義は 1 か所）。
8. Generate 中の操作を 1 つの Product Decision にまとめる。対象は Undo、inert でない overlay の後ろでのキーボード編集。推奨は、一貫して modal にすること。
9. Worker を終了させるのは、Worker の中の処理を取り消すときだけにする。要求ごとの所有者の区別は今は不要。
10. 確定済み request を使い回すときは、必ず `promoteCanonicalRenderRequestToCurrent` を通す。

## 監査の分類の誤り・不具合ではないもの
- **SE-06:** `imported-comparison-intent.js:191` は不具合ではない（上記）。本当の原因は識別子と Inherit の昇格漏れ。
- **SE-08b:** lazy Worker の規則には反していない。tobacco だけでなく、CLI Session 全般の設計上の緊張関係。
- **SE-08a:** tobacco だけの問題ではない。Linear の CLI Session 全般に及ぶ。
- **SE-02:** ラベルの対応づけだけの問題ではない。キーボード操作も同じ原因。
- **SE-04:** range の input は存在しない。
- **GE-06:** state.js を原因とするのはずれている。担当は History。
- **SE-01:** 同じ原因の潜在経路（import rollback）が、コードから見て残っている。
- **drag:** `drag-actions` と `diagram-drag` に、SE-03 と同じ型の潜在不具合があるとコードから判断した。
