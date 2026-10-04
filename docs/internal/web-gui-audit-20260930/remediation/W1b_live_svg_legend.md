<!-- Raw design report of workstream W1b (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W1b 修正提案: live SVG の Result 書き込み、凡例編集、装飾配置の継続（dev 4c89bab1）

コードは変更していない。実行した実験は1つだけで、scratch 上で dev CLI を動かした（PV-04 の新規所見）。ブラウザでの再現はしていない。行番号は repo 直下からの `path:line`（4c89bab1 時点）。

---

## IN-01（= GE-01, P1）

**1. 分類: PRODUCT_DECISION_REQUIRED**
- 現状は OIC-009/OIC-010 に違反し、科学的出力を黙って壊す。Residual-risk boundary はこれを許さないので、維持は NOT_ALLOWED。
- **A（推奨）**: Circular の Species / Strain / Plot title / Title position / Title font / Default font size / Keep Full Definition を Applies on Generate にし、live で定義を書き換える経路を削除する。
- **B**: live を残す。committed canonical Session にこれらのフィールドだけを投影し、正規 renderer で再描画する。
- **C（不採用）**: helper に committed の selector / region / reverse / label / 順序を渡す。record planner の複製になる。定義の寸法による track の再配置（GE-03 で def font 30 → CANNOT_FIT）を再現できない。
- A を推奨する理由:
  - Linear は 15f0d3e2（2026-04-05）で同じ watcher を外し、Applies on Generate になっている。`gbdraw/web/js/app/results.js:371-475` の Linear 分岐は到達不能。A で両モードが揃う（LSP）。
  - `docs/REFERENCE/web-app.md:700` に "Generate again after changing color, font, title text, or legend side" とある。
  - 削除だけで P1 を即時に止められる。
  - B には非ブロッキングの committed 投影再描画が必要だが、まだない。性能の計測も要る。

**2. 根本原因の再確認**
- 監査のとおり。`app/results.js:254-296` が draft の `files.c_gb` の全 record を `python-helpers.js:1734` で parse し直し、record planner を通さずに `DefinitionGroup(gb_record=record)` を組む。
- 書き戻し（`:354-361`）は committed request を更新しない。そのため Session では Result と request も食い違う。

**3. 修正案（A）**
- owner: `app/results.js`（createResultsManager）。
- 削除するもの:
  - `results.js` の `:14-55`、`:132-138`、`:177-232`、`:234-492` と関連 import
  - `app/watchers.js:186-196, :574-585` と cancel 呼び出し（`:191, :258, :376`）
  - `app/app-setup.js:2458, :2487-2489, :1403, :2670, :2882, :3230`
  - `app/python-helpers.js:1693-1827, :1866`
  - `services/diagram-worker-protocol.js:15`
  - `workers/diagram-generation-worker.js:745-773`
  - operation 名（`services/error-normalization.js:4` と `gbdraw/web_support/error_adapter.py:20`）
  - 呼び出し元がなくなる `legend-layout/reposition-actions.js:155-173`（refreshCompositionGeometry）と `legend-layout/composition-actions.js:850-894`（reconcileCompositionTitle）
- 追加するもの: `index.html` の Species/Strain（`:2413-2424`）と Titles & Record Labels に、既存の「Applies on Generate」表記（`index.html:1472` と同じ span）。
- 追加しないもの: guard、新しい state、helper の引数。
- B を選ぶ場合: `services/session-request.js:4942` の projectCommittedRecordTransform と同じ型の committed projector を1つ置く。投影するフィールドは既存の CONFIG_OVERRIDE_PATHS を使って data として宣言する。削除対象は A と同じ。

**4. 上位の対応**
- Decision Pack を作り、承認後に OIPC へ static co-change で PD を記録する。
- `web-app.md:89-92` の表を更新する。
- web/CLAUDE.md に横断ルールを追記する（後述の3）。
- **guard-only PR を先に出す**: `tests/web/architecture-contracts.test.mjs:307-314`（`app/results.js` を importer とする厳密な列挙）と `:593`（results.js の件数 4）が、runtime での縮小を阻む。
- `tools/web-change-policy.json` から `app/results.js` を外す縮小は、runtime と同じ PR でよい（narrow contraction の規定内）。

**5. テスト**
- Playwright（`tests/web/gui-audit-regressions.playwright.spec.js` に追加）:
  - Region 1000–9000 と Record label を付けて Generate する。
  - 定義系フィールドを編集して 1 秒待つ。`results[0].content` と Save Session の Result がバイト一致すること。
  - Generate 後は `8,001 bp`、label を保持し、Species を含むこと。
  - grid の並べ替えと、Generate 前のファイル差し替えでも `data-gbdraw-record-id` が変わらないこと。
- 削除: `tests/test_web_feature_metadata.py:137-230` と `tests/web/definition-layout-completion.test.mjs` の2本目。
- 元の不具合を検出できた assertion: 「Generate 以外の経路で定義の text / length / GC% / record-id が変わらない」。

**6. 規模と PR**
- production は 10 files、gross 約 650、net 約 −600。
- 変更ファイル数が Ordinary の上限 8 を超えるので、Architecture profile + `architecture-change` label（lifecycle 経路と Worker の helper op を削除するため）。
- Review REQUIRED: size profile、architecture-bearing な削除、PD の co-change。
- 単独 PR。GE-01 も同時に close する。

**7. 依存とリスク**
- 先に guard-only PR が要る。
- B を選ぶなら、横断7の reflow 仮説の解消が前提。
- 既存 Session に保存された誤った Result は、再 Generate すれば正しくなる。migrator は不要。

---

## GE-02（P2）

**1. 分類**
- (a) 空値や不正値を Result に書かない: **IMPLEMENT_EXISTING_AUTHORITY**（OIPC-C01/C03。これで OIC-013 の「Result を保持」表示が実態どおりになる）。
- (b) live 経路の扱い: **PRODUCT_DECISION_REQUIRED**
  - A: 経路を維持して修正する。
  - B: committed 投影で自動再描画する。
  - **C（推奨）**: global stroke の live 反映をやめ、Applies on Generate にする。
- C を推奨する理由（コードで確認）:
  - `app/svg-styles.js:420-432` は、per-feature の stroke override を持つ block も上書きする。Generate では `candidate-render.js:147-170` の override が勝つ。
  - `:496-540` の凡例 swatch 判定は path 形状によるヒューリスティック。dev CLI の実 SVG では `#feature_legend` がなく、GC content の swatch（stroke="none"）も条件に一致する。Generate ではこの swatch に block stroke は付かない。
  - Auto の値を JS は持っていない（`auto-value-display.js` は表示用）。
  - batch では表示中の Result だけが変わる。
- 退役には Decision が必要: `index.html:4329/4458` がこれを Live edit と明示しており、`tests/web/history-generated-authority.playwright.spec.js:180` もテストしている。

**2. 根本原因の再確認**
- 監査のとおり。空文字は `!== null` を素通りする。null を受けても Generate 時の値に戻さない。`v-model.number`（`index.html:4082, 4098, 4352, 4372, 4493`）が '' を入れる。
- 「失敗した Generate が Result を保持したと表示する」点は、表示としては正しい（live 編集を含む最後の成功 Result が保持されている）。欠陥は不正値を受け付けることだけ。

**3. 修正案**
- C: `svg-styles.js:409-548`（applyStylesToSvg）、`:650-667`（watcher）、`:681`（export）を削除する。`index.html:4329/4458` の文言を Applies on Generate に直す。
- A を選ぶ場合: 同じ owner で、null と不正値は書かない。対象は `[data-gbdraw-feature-part="block"]` のうち override のないものに限る。それでも凡例 swatch の判定に正本がなく、近似が残る。

**4. 上位の対応**: 同じ Decision Pack に concern として加える。`web-app.md:91-92` を更新する。

**5. テスト**
- C: Generate 後に stroke を 5 → 空 → −1 と変えても Result がバイト一致すること。−1 のまま Generate すると validation error になり、Result もバイト一致のままであること。
- `history-generated-authority.playwright.spec.js:180-200` を書き換える。
- A: live 後と fresh Generate で、block、GC swatch、override 付き feature の stroke 属性を比べる。

**6. 規模と PR**: C は 2 files、gross 約 150、net 約 −140。Ordinary。

**7. 依存とリスク**
- 不正な数値の扱いとインラインのメッセージは X-02（W3）。batch への伝播は FE-02/PV-09（W1a）。
- `originalSvgStroke`（`canvas-actions.js:87-97`）は per-feature stroke の Auto 復元に使われているので残す。

---

## PV-01（P2）

**1. 分類: 受け入れ済みの残余リスク**
- PD-OI-052 の Accepted residual risk に「clipping/overlap、自動 clamp なし」とある。照合キーに side は含まれない。`web-app.md:700-706` も、side を変えた後も offset を保持すると書いている。したがって欠陥ではない。
- 変えるなら PD-OI-052 revision 2 の Decision が要る。
  - **A（推奨）**: 現状維持。
  - B: side を変えたら delta を 0 にする。無言の削除になるので、開示が必須。
  - C: bbox が viewBox と交わらないときは、既存の DECORATION_CONTINUITY で候補を公開する前に止め、Reset または padding を案内する。metadata からの算術で判定できる。対象は `decoration-continuity.js:72-107`。
- 注意: PV-07 のため、052 が挙げる緩和策「padding で調整」は毎回の Generate で消える。改善するなら PV-07 の B が先。

**2. 根本原因の再確認**: 監査の説明は正確。ただし決定済みの挙動である。

**3. 修正案**: A なら変更なし。C を選ぶ場合は `decoration-continuity.js:72-107` に判定を加える（1 file、約 25 行）。

**4〜7**: A なら不要。C を選ぶ場合は OIPC の supersession、unit テスト（画面外で停止、部分的な clip は通過）、Ordinary。リスクは Generate が止まる負担。

---

## PV-02（P2）

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**
- `web-app.md:101-102` の "Supported color, label, visibility … edits are carried forward"。
- `index.html:6512` の「legend text」の Live edit 表示。
- `candidate-render.js:228-242` の `legendRenames` の契約（`tests/web/svg-result-ingestion.test.mjs:394-406` が検証）。
- feature のある項目の rename は保持される（LSP）。

**2. 根本原因の再確認**
- `feature-editor/color-actions.js:661`（syncOriginalLegendMetadataRename、`:421-440`）が originalLegendOrder と originalLegendColors を書き換える。`:666` が originalCaption を上書きする。
- その結果 caption === originalCaption となり、Generate 時に rename の op が生成されない。

**3. 修正案**
- owner: `renameLegendEntryInSvg`（`color-actions.js:623-676`）。
- renderer が生成した行（originalLegendOrder に含まれる行）では、caption だけを更新し、original 系は変えない。
- syncOriginalLegendMetadataRename は、renderer が新しい caption を出す rule 経路（`:701`）だけで使う。
- 新しい state や field は追加しない。

**4. 上位の対応**: 不要。

**5. テスト**
- node（`feature-color-actions.test.mjs`）: rename 後に originalCaption が 'GC content' のままで、`compileDirectEditorMutationPlan` が `legendRenames` を出すこと。
- Playwright: Circular と Linear で rename → Generate、Save → Load → Generate、batch で保持されること。

**6. 規模と PR**: 1 file、約 15 行。PV-04、PV-12 と同じ PR にする（2 files、約 60 行、Ordinary）。

**7. 依存とリスク**
- `legendRenames` は xPos/yPos を再生する（`svg-result-ingestion.js:455-460, :402-425`）。Sort や side 変更の後は古い anchor に戻って重なる可能性がある（未検証）。rename 後の幅も Generate 時には reflow されない。horizontal legend + sort + rename + Generate の Playwright を必須にする。
- 既存 Session で上書き済みの identity は復元できない。

---

## PV-03（P2）

**1. 分類**
- 現行の authority では欠陥ではない。PD-OI-052 が「legend順の新しい継承保証は含めない」とし、`web-app.md:105` もそう書いている。
- 保持するなら **PRODUCT_DECISION_REQUIRED**:
  - A: 現状維持にして、その旨を表示する。
  - **B（推奨）**: 既存の candidate mutation plan に順序の op を data として追加する。順序の意図は Session の `editorState.legend.entries`（`services/config.js:1027`）に既に保存されているので、新しい schema は要らない。
  - C: Python に request field として渡す（YAGNI）。
- B を推奨する理由: 052 と同じ「やり直しの負担」の journey であり、既存の境界に乗る。

**2. 根本原因の再確認**: 監査のとおり。順序は Result の DOM にしかなく、`legend/entry-actions.js:861-887` の抽出で作り直される。

**3. 修正案（B）**
- `compilePlanBundle`（`candidate-render.js:103-316`）で、順序が originalLegendOrder と違うときだけ `legendOrder` op を出す（全 Result が対象）。
- 実行は `svg-result-ingestion.js:427-` の applyLegendOperations。slot を割り当てるアルゴリズムは `legend/sort-actions.js:71-128` から移し、live の sort も同じ関数を呼ぶようにする（DRY）。未知の caption は後ろに回す。適用順は最後にする。

**4〜7**
- 上位: PD の記録と `web-app.md` の更新。
- テスト: node で op の生成、zero fast path、未知 caption。Playwright で Sort Z-A → Generate、Move → Save → Load → Generate、batch。
- 規模: 3 files、net 約 +30、Ordinary。PV-07 と同じ PR にできる。
- リスク: PV-02 の anchor 再生との相互作用。

---

## PV-04（P2）

**1. 分類: PRODUCT_DECISION_REQUIRED**
- UNKNOWN と表示される点は PD-OI-046 違反。これは X-01（W3）の mapping で扱う。
- **A（推奨）**: feature のない項目で使っている既存の target dialog（Merge / Suffix / Cancel、`color-actions.js:725-735, :801-829`）を、feature のある項目にも使う。
- B: PD-OI-042 を拡張し、palette や生成された凡例 key との衝突も hex suffix で自動的に区別して通知する。
- C: 明示的な validation エラーとして拒否する。
- 例外: 衝突先が specific-color rule の caption なら、dialog を出さず既存の PD-OI-042 の正規化に任せる。

**2. 根本原因の再確認**
- 監査のとおり。`continueLegendRenameRequest`（`:737-841`）の feature/rule 分岐（`:780-785`）が、衝突判定より前に return する。`legend/entry-actions.js:731-736` は plain Error を投げる。
- **新規所見（監査外、dev CLI で確認）**: rule の caption が既定の type 名と同じで色が違うと、CLI の凡例が rule の色を落とす。
  - `tRNA product .* #ff0000 rRNA` + HmmtDNA で試した。凡例は `[rRNA, CDS, GC content, GC skew (+), GC skew (-)]` で、rRNA は既定色 `#71ee7d`。赤い tRNA 22 個に対応する凡例行がない。
  - 原因は `gbdraw/legend/table.py:175-230`。`gbdraw/features/colors.py:39-51` は rule の caption しか予約しない。
  - PD-OI-042 の「既存 legend key と衝突した場合」がこのケースを含むかは解釈が必要。Python の owner 向けで、別の workstream。

**3. 修正案（A）**
- 衝突判定（`:787-798`）を分岐の前に移す。
- 衝突先が rule 所有でなく、色が違えば dialog を開く。merge / suffix は rule 経路にも通す。`:780-785` の早期 return は削除する。
- `entry-actions.js` の throw は防御として残す。

**4〜7**
- 上位: Decision（B を選ぶなら PD-OI-042 の改訂）。
- テスト: node で feature ありの tRNA→rRNA が dialog を開き commit しないこと、merge / suffix / rule 所有のケース。Playwright で errorLog が null であること。
- 規模: 1 file、約 40 行。PV-02 と同じ PR。
- リスク: S03 の PD-OI-042 実装との境界。

---

## PV-07（P3）

**1. 分類**
- 欠陥ではない。PD-OI-052 は padding を継承保証の対象外としている。
- **PRODUCT_DECISION_REQUIRED**:
  - A: 現状維持（Generate ごとに 0）。
  - **B（推奨）**: Session の `ui.canvasPadding`（`config.js:3728`）の単一の値を、候補を公開する前に全出力へ適用する。
  - C: request/CLI の field にする（YAGNI）。
- B を推奨する理由: 052 が示す緩和策を持続させられる。既存の seam を再利用するだけで済む。

**2. 根本原因の再確認**: `app/app-setup.js:2292-2302` が 0 に戻す。padding は `canvas-actions.js:42-70` で表示中の Result にだけ後から適用される。

**3. 修正案（B）**
- `app-setup.js:2298-2301` の 0 への初期化を削除する。
- decoration continuity と同じ候補 transform に、padding を1つ加える。
- viewBox の計算は、composition の `editorPadding`（`composition-actions.js:833`）と共通の純関数にまとめる。

**4〜7**
- 上位: PD の記録と文言の更新。
- テスト: Right 150 → Generate で幅が 1641。Save → Load → Generate、batch、Reset。二重に適用されないこと。
- 規模: 3 files、約 50 行、Ordinary。
- リスク: Load（trusted restore）時の二重適用。

---

## PV-10（P3）

**1. 分類: EVIDENCE_REQUIRED**
- 結果に応じて、IMPLEMENT_EXISTING_AUTHORITY か PRODUCT_DECISION_REQUIRED に分かれる。
- 根拠:
  - live の side 変更は長く続いている、テスト済みの機能（`tests/web/composition-layout-real.playwright.spec.js:306-325`）。
  - JS の composition は Python oracle と parity テストされている（`tests/web/composition-runtime-parity.test.mjs`、`composition-python-oracle.py`）。
- 差分の出所は未特定。仮説（未検証）は、legend reflow の入力（`legend/layout-actions.js:478-`、widthHint は `reposition-actions.js:127-134`）の寸法。
- 推奨: まず evidence を取り、parity を直す。一致させられない場合は、A（退役、GE-07 も消える）か B（再描画）の Decision。

**2〜7**
- 観測: Linear と Circular の各 side で、live と fresh Generate の legend bbox と viewBox を比べる。
- 差が reflow にあれば、oracle に reflow のケースを加え、`layout-actions.js` を一致させる。補正係数や第2の計算経路は足さない。
- 規模は evidence の後に決める（parity 修正なら 1〜2 files、Ordinary）。GE-07 と同じ PR。

---

## PV-11（P3）

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**
- OIC-006、OIC-025、PD-OI-038 の must-preserve（keyboard、visibility のみ）。
- 先例: `web-app.md:591` の review は、閉じると起動元へ focus を戻す。
- 「Escape でどこからでも閉じる」範囲そのものは規定がないので、変更しない。

**2. 根本原因の再確認**: `app/ui.js:350-358` と Close ボタン（`index.html:6477`）が `app/right-drawer.js:92-99` を呼ぶが、focus を移さない。

**3. 修正案**
- owner は `createRightDrawerController`（`architecture-contracts.test.mjs:474` で唯一の owner と固定）。
- `closeRightDrawer` で、focus が drawer の中にあれば Editor toggle（`index.html:6459`）へ戻す。drawer 内かどうかと戻し先は、狭い callback で受ける。
- focus trap や watcher は追加しない。

**4〜7**
- テスト: `right-drawer.playwright.spec.js`。1600 と 390 の幅で、Escape と Close の後に toggle に focus があること、drawer の外に focus があれば動かないこと。
- 規模: 2 files、約 25 行、Ordinary、単独 PR。
- 上位の対応: 任意で web/CLAUDE.md に1行追記。

---

## PV-12（P3）

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**（OIC-024 と PD-OI-037 が求める表示の真実性）。色の操作を追加するなら、新しい機能として別の Decision。

**2. 根本原因の再確認**: 監査のとおり。色の入力欄は一度も存在していない。文言は 82e9756d（2026-09-26）で入った。`updateLegendEntryColor` は UI に接続されておらず、テストの app API からしか使われていない。

**3. 修正案**: `index.html:6512` を実際の操作に合わせる（例: "legend text, order, stroke, and removal"）。関数はテストが使っているので今回は残す。

**4〜7**: 上位の対応は不要。テストは文言の assertion。1 行で、PV-02 と同じ PR。リスクなし。

---

## GE-07（P3）

**1. 分類: IMPLEMENT_EXISTING_AUTHORITY**（web/CLAUDE.md の「利用できない要求は決定的な fallback に解決する」、OIPC-C07）。

**2. 根本原因の再確認**
- 監査のとおり。`reposition-actions.js:116-119` の throw は、`watchers.js:234-253` の nextTick の中でしか到達しない。
- Circular でも同じ分岐に入る可能性がある（未検証）。

**3. 修正案**
- throw を `return false` に置き換える。Result は変えず、次の Generate で draft が適用される。
- PV-10 が A になれば、分岐ごと削除する。
- watcher に try/catch は足さない。

**4〜7**
- テスト: Linear と Circular で、None → Generate → Left の後に pageerror が出ず Result がバイト一致すること、Generate で左に出ること。
- 規模: 1 file、−3 行。PV-10 と同じ PR。ただし PV-10 の evidence が長引くなら先に出す。

---

## 横断的な提案

1. **共通原因 A: Result の「第二の書き手」。** 正規の renderer 以外が近似で Result を書き換え、しかも入力に draft を使う（IN-01、GE-02、PV-10）。どれも表示中の Result だけを直接書き、committed request は更新しない。
2. **共通原因 B: live 編集の記録が、Generate 時の再適用コンパイラ（`candidate-render.js::compilePlanBundle` → `svg-result-ingestion.js`）の前提とずれている**（PV-02）。あるいは再適用の op がそもそもない（PV-03、PV-07）。live 側は executor と別に DOM を実装していて、DRY に反する。
3. **ルールを1つにする（web/CLAUDE.md の live-edit invariants に追記）。** live 編集が current Result を変えてよいのは次の3つだけとする。
   - (a) Generate のコンパイラが同じ executor で再適用する editor intent
   - (b) Python oracle の parity で守られた composition owner の編集
   - (c) committed Session + 宣言済みの投影フィールド + intent による自動再描画

   draft や helper による部分的な再描画は Result に書かない。どれにも当たらない設定は Applies on Generate と表示する。PD-OI-052 の「同じ候補境界」を一般化したもの。
4. **強制の手段 1（ratchet）。** `tools/web-change-policy.json` の "Mounted SVG/Result replacement" の allowlist を、縮小方向にしか動かない ratchet として使う。そのため先に guard-only PR を出し、`architecture-contracts.test.mjs:580-600` の件数を厳密な値から上限（`:258-266` legendReplacementCount と同じ形）に変える。S05 では、この厳密な件数のために runtime の縮小 PR が進められなかった。
5. **強制の手段 2（振る舞い）。** table-driven の Playwright で「live→Generate parity」を確かめる。各 live 操作について、live 後の Result と、同じ状態から fresh Generate した Result の対象要素を比べる。Applies on Generate の項目は「Generate 前は Result がバイト不変」を assert する。IN-01、GE-02、PV-02、PV-10 はすべてこの1本で検出できた。
6. **Python helper の登録表**（`python-helpers.js:1855-1872`）。SVG 断片を返す helper は正規 renderer の外にある描画経路なので、追加には architecture review を求める。
7. **検証が必要な仮説（W1a と共有、未検証）。**
   - `runLabelReflow` は draft（`run-analysis.js:1898, :4579`）から request を組み、adopt しないまま Result を置き換える（`:4770`）。そのため、Applies on Generate のはずの draft がラベルの reflow で Result に漏れる可能性がある。
   - 正しい入力は、`session-request.js:4942` の「No live form state participates」と同じく committed Session。
   - IN-01、GE-02、PV-10 で B を選ぶなら、これが前提になる。
8. **表示の真実性（OIC-024）。** Result に影響するすべての項目に、Live edit か Applies on Generate のどちらかを表示する。今は Species/Strain、Titles & Record Labels、Legend Position に表示がない。
9. **Decision は1つの Pack にまとめる。** IN-01、GE-02、PV-04、PV-03、PV-07、任意で PV-01 を、それぞれ独立した concern として1つの Decision Pack に入れ、推奨案を一括で承認できる形にする。IN-01 を先頭に置く。
10. **batch。** 規則 (a) の executor を live でも使えば、Generate と同じく全 Result に適用される。W1a の FE-02/PV-09 と owner（`svg-result-ingestion.js` の executor）を共有できる。

---

## 誤分類、またはバグではないもの

- **PV-01**: 受け入れ済みの残余リスク（PD-OI-052、`web-app.md:700-706`）。「P1 に近い」は過大。
- **PV-03、PV-07**: 052 が継承保証の対象外と明記している。機能要望で、Decision が要る。
- **GE-02 の「Result を保持したと表示する」**: 表示は正しい。欠陥は不正値を受け付けることだけ。
- **PV-10 の「Generate 前の export が live の配置になる」**: live 編集の定義どおり。欠陥は live と Generate の parity の差だけ。
- **IN-01**: Linear には同じ不具合はない（dead code）。GE-01 と重複している。
- **GE-07**: Linear に限らない可能性がある。
- **PV-12**: 機能の欠落ではなく、表示の不一致。
- **新規（監査外）**: CLI の凡例衝突（PV-04 の新規所見）。
