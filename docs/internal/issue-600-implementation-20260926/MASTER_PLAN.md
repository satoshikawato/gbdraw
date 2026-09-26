# Issue #600 — 総合実装計画

Status: **４つの Product outcome は承認済み、runtime は未実装**。
この文書は実装範囲・所有境界・実行順・検証・ブランチ運用の唯一の計画 owner。
Product authority は既存 `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` が保持する。
この計画の記載だけでは、未登録の仕様を同じ候補 runtime に適用する根拠にはならない。

- Repository: `https://github.com/satoshikawato/gbdraw.git`
- Issue: [#600](https://github.com/satoshikawato/gbdraw/issues/600)
- Implementation branch: `fix/issue-600-annotations-styles-20260926`
- Initial base: `origin/dev@d457b7189b137185a8dec800819a312c30b969fa`（2026-09-26 に fetch 後確認）
- Approval: [承認済み Product Decisions](APPROVED_PRODUCT_DECISIONS.md)
- Initial observations: [EVIDENCE.md](EVIDENCE.md)

## 1. 問題と完成時の挙動

gbdraw は Python の genome 図描画エンジンと、ブラウザー内で動く Web UI を持つ。
Circular は single/grid/batch、Linear は複数 record を並べる１つの図を生成する。
annotation は領域の line/bracket/band/highlight、specific-color rule は feature の
色と凡例 caption、track slot は図上の track 配置と寸法を指定する。

| 問題 | 調査した変更前の挙動 | 承認済みの完成挙動 |
| --- | --- | --- |
| BUG-08 | Annotation TSV に notes 等の未知列があると Web/Python とも import 拒否 | 未知列を捨て、列名と非保存を通知。既知値・必須列・重複・不正 row は検証 |
| BUG-09 | １ selector の未一致で図全体が失敗 | binding 成功済みの該当 annotation 行だけを warning 付きで skip。partial geometry は作らない |
| BUG-10 | Web は同名異色を拒否。Python 凡例は最後の色が先の色を上書き | canonical rule の caption に hex suffix を付けて区別し、各使用色を既存 solid 凡例行で説明 |
| BUG-14 | validator/draft/payload の px 解釈が不同意。不正値が null/0 になる経路もある | 純 pixel の listed slot fields で optional px を受理し、共有 parser と現行 typed 型へ収束 |

対象の仕様は [承認記録](APPROVED_PRODUCT_DECISIONS.md) に完全な条件で記載する。
４件の scope は独立であり、他の TSV、record broadcast、multi-swatch 凡例、任意 CSS
unit は対象外。既存有効入力の図、最後の成功 Result、保存 preview、private data を保持する。

用語: **draft** は編集中の user state、**canonical request** は Python が検証する
描画要求、**Result** は最後に成功して取り込まれた図、**Session** は入力・設定・
保存 preview を持つファイル、**owner** はその意味を決める単一の責任箇所。

## 2. Architecture / SOLID・KISS・DRY・YAGNI

既存 canonical path を維持する。

```text
Web editor / Session
  → session-request.js の canonical projection
  → run-analysis.js の orchestration
  → diagram-generation.js / diagram Worker
  → Python typed request / shared preparation / resolver
  → drawing + structured warnings
  → sanitized Result admission / live editing / export
```

- **SRP:** codec は表の conversion、resolver は annotation target の解決、feature-input
  owner は caption 正規化、UI は owner action と通知を担当する。
- **OCP / LSP:** canonical contract の明示的な extension で実装し、valid-input の
  戻り値・規則順序・selector・座標・transform を保持する。catch-all fallback を足さない。
- **ISP / DIP:** 既存 request/result、immutable bundle、typed Worker helper の小さい
  境界を使う。別 generation pipeline や class 階層を増やさない。
- **KISS:** strict/lenient toggle、fuzzy header correction、broadcast、multi-swatch
  renderer、汎用単位変換を追加しない。文書とセッションも１つの実装ブランチへ収束する。
- **DRY:** Python の caption normalizer を native/Web で共用する。Web の pixel parser
  を validator/normalizer/payload で共用し、旧同義 parser を同じ変更で除去する。
  生物学的 matching を JS に複製しない。承認本文は別ファイル、acceptance は本書が owner。
- **YAGNI:** unknown metadata 保存 schema、speculative compatibility reader、final-result
  cache、新しい decision UI/evaluator/registry を作らない。

`CLAUDE.md` と `gbdraw/web/CLAUDE.md`、architecture/Product ratchet に従う。
本書は Product に新しい outcome を選び直させる手続きではなく、承認済み outcome を
既存 durable authority と runtime に正しく反映する手順である。

## 3. セッションの順序・所有範囲

セッションは **S00 → base authority 確認 → S01 → S02 → S03 → S04 → S05** の
順で実行する。同じ remote branch に対する実装 writer は常に１セッションだけ。
独立した４つの関心を、共有 branch の競合回避のために逐次実装する。

| Session / Instruction prompt | 開始条件 | 担当と終了条件 |
| --- | --- | --- |
| [S00 — authority](sessions/S00_AUTHORITY.md) | 本書と承認記録を取得 | ４ concern を別 record として静的契約へ登録。authority-only commit を push。dev 統合は別の maintainer 境界 |
| [S01 — TSV](sessions/S01_TSV_IMPORT.md) | ４ outcome の durable authority が origin/dev に存在 | BUG-08、TSV-01～03、focused/browser tests と session report を commit/push |
| [S02 — selector](sessions/S02_SELECTOR_RESOLUTION.md) | S01 の pushed completion | BUG-09、SEL-01～05、resolver/reporting を通して commit/push |
| [S03 — colors](sessions/S03_COLOR_CAPTIONS.md) | S02 の pushed completion | BUG-10、CLR-01～05、normalization/admission/live/native を一体実装し commit/push |
| [S04 — pixels](sessions/S04_PIXEL_INPUT.md) | S03 の pushed completion | BUG-14、PX-01～03、旧 parser を削除し commit/push |
| [S05 — integration](sessions/S05_INTEGRATION.md) | S01～04 の pushed completion | INT-01～02、regression/policy/browser/output verification と必要修正を commit/push |

Product concern を４つの独立 record として維持することと、authority-only delivery を
１セッションで行うことは別である。policy の証拠置換順序を無視して１ runtime PR に
authority を混ぜない。S00 の diff は計画/evidence/authority だけであり、runtime を含めない。
この時点の work branch を authority-only 候補として dev に統合した後、同じ名前の
work branch を再取得して runtime commits を進められる。S00 の authority-only 統合時は
この remote work branch を残す。自動削除された場合は maintainer が同名 branch を
統合済み dev から復元し、remote に計画と authority があることを確認してから S01 を開始する。

S00 終了時に authority が dev に未統合なら、S01～04 を開始しない。
S00 の candidate file は自身の runtime を承認できない。
PR 作成や dev merge は、この計画の commit/push 指示だけから許可を推測しない。
既存の有効な publication 承認がある場合は、その target/条件を維持して実施する。

### 3.1 他セッションに干渉しない checkout

**毎回、remote の当該実装ブランチを取得して使用する。元の共有 checkout や他の
worktree の branch を切り替えず、独立 clone に同名 branch を作る。**
同一 Git repository の複数 worktree に同名 branch を checkout する制約を避けるため、
実装セッションの標準手順は独立 clone とする。検証のために repo-global Git 設定や
他セッションの server/process/port を変更しない。browser fixture/server は当該 clone
と専用 ephemeral port、tmp directory だけを使う。

```bash
SESSION_DIR="$(mktemp -d /tmp/gbdraw-issue600-session-XXXXXX)"
git clone --no-checkout https://github.com/satoshikawato/gbdraw.git "$SESSION_DIR"
cd "$SESSION_DIR"
git fetch origin
git switch --no-track -c fix/issue-600-annotations-styles-20260926 origin/fix/issue-600-annotations-styles-20260926
git status --short --branch
git branch --show-current
git for-each-ref --format='%(upstream)' refs/heads/fix/issue-600-annotations-styles-20260926
```

開始時に remote work branch SHA を記録し、前 session の report と実コードを読む。
base dev の取得は freshness/authority 確認のためであり、work branch を新しい dev tip
から作り直して計画や前 session の commits を失ってはいけない。
S01 は authority を含む `origin/dev` を取り込み、以後も material base change のみ
適切に統合して必要な checks を行う。許可のない force-push/reset/cherry-pick による
他 session の履歴置換は行わない。

### 3.2 各セッション終了時の commit / push

**担当実装・focused verification・必要な fixes・session report を完成させ、同名の
remote work branch に commit と push を行ってから終了する。** 未検証の draft を
完成扱いにしない。S00 の authority merge 待ちや、新しい Product outcome が本当に
必要な場合は未完了理由を report に明記し、許可された独立成果だけを push する。

1. `SESSION_RESULTS/<session>.md` を作成。start work SHA、authority base SHA、変更 owner/
   path、実際の command と exit/result、未実行範囲、artifact、残る問題を記す。
2. production/tests/docs/generated diff を別々に review。`git diff --check` と、既存
   acceptance を弱めていないことを確認。session completion の要件に結び付ける。
3. commit 前に `git fetch origin` し、remote work branch が start SHA と同じか確認。
   他 writer が進めていたら作業を停止して実際の remote state を調べる。force-push で
   解消せず、担当者が統合と再検証の範囲を決める。
4. 現 branch が `fix/issue-600-annotations-styles-20260926`、upstream が未設定または
   `origin/fix/issue-600-annotations-styles-20260926` であることを確認。
5. `git add -- <担当変更の明示 path>` で scope だけを stage。無関係な変更は含めない。
   session を１ commit にまとめ、英語 title/summary を使う。
6. 次の明示 refspec で push し、`git ls-remote` の SHA が local `HEAD` と一致することを
   確認。push エラー後は先に実際の remote state を確認し、成功済み操作を繰り返さない。

```bash
git push --set-upstream origin HEAD:refs/heads/fix/issue-600-annotations-styles-20260926
git rev-parse HEAD
git ls-remote --heads origin refs/heads/fix/issue-600-annotations-styles-20260926
```

既存 branch の advance は normal fast-forward push のみ。main/dev への直接 commit/push、
他の remote branch への push、merge、deploy、release/tag はこの指示の対象外。
次 session に pushed SHA、検証結果、残る具体的境界を伝える。会話ログを再現する
handoff ではなく、当該 branch の report とコードだけで次の担当者が開始できるようにする。

## 4. BUG-08: table admission

対象: `gbdraw/web/js/app/annotations/table-codec.js`、`annotations.js`、
`gbdraw/annotations/io.py`、必要最小限の Web markup。

1. header を正規化し、必須列と重複列を検証する。BOM/空白は現在の表入力契約と
   一致させる。必須列の綴り違いは missing required のまま。
2. annotation の既知列だけを射影する。未知列を `metadata` や request params に
   転写しない。未知列名を１ import 当たり１件の集約通知として返す。
3. Web の既存 `parseAnnotationTable(text)` の配列戻り値を維持する。
   codec 内部の１処理を、既存戻り値用 wrapper と UI 用の通知付き結果から共有する。
   別々の parser にしない。Python の tuple 戻り値も維持し、未知列の summary を
   logger に配送する。notification のために render worker を起動しない。
4. 完全に parse/validate した後、`replaceSets` を一度だけ呼ぶ。known enum、数値、
   target、ID、style が不正なら全 import を拒否し、直前 draft/Result を保つ。
5. `fill_colour` も「未知列」として列名を通知する。類似語から `fill` と推測しない。
   TSV 再出力・Session に付加列が残らないことを通知する。

Python の DataFrame 列を set 化する前に、正規化後の重複を検証する。
ファイル経由は pandas が duplicate header を `id.1` 等に変える前に header を
検証する。未知列許可によって重複が未知列扱いで通過する穴を作らない。
行の列数が header と違う場合は拒否し、ヘッダにある未知列と malformed row を
区別する。これらは受入範囲拡大に必要な境界保護であり、汎用 TSV parser の
全置換を目的にしない。

## 5. BUG-09: annotation resolution and reporting

対象: `gbdraw/annotations/resolve.py`、`models.py`、`planning.py`、
`gbdraw/api/request_render.py` と builder 引数の既存境界、
`gbdraw/web_support/request_render.py`、Web Result admission/通知。

- `_feature_segments` は matched geometry と未一致の事実を返す private result に
  分解する。selector syntax/record binding の例外を catch-all で握りつぶさない。
- `resolve_annotation_set` が set/row/record の文脈を付け、１件でも未一致ならその
  行の geometry を生成せず `feature_selector_unmatched` を１件返す。
  同じ selector 重複は warning 内で重複排除する。未一致 selector の値や genome
  contents を console に dump せず、UI は record/set/annotation ID と欠落件数を示す。
- すべて未一致・一部未一致を同じ row-skip にする。完全一致だが drawable location が
  空なら既存 `empty_span`。coordinate の clip/skip/error は既存 policy のまま。
- record ID の欠落、重複、index 範囲外、multi-record の record 省略は fatal を維持。
  `record` 空欄を「全 record」へ変更しない。複数 genome への同じ注釈は明示的に
  record を指定した複数行で表す。broadcast は別の将来要件。
- `ResolvedAnnotationBundle` を解決された値として drawing と reporting が共有する。
  typed plan の record materialization/transform 確定後に一度だけ解決し、single/grid/
  batch と Linear の既存 builder 引数へ渡す。grid の既存事前解決を置換し、
  assembler が受け取った bundle を再解決しない。batch は record index の元の意味を
  保って分割する。CLI/API も同じ resolver を使い、API が解決結果の warning を
  取得できる immutable result field と CLI の表示 adapter を用意する。
- 新しい logger 捕捉・global warning collector・Drawing への private side-channel
  を作らない。`PreparedDiagramRequest` → `RequestRenderResult` の additive field で
  warning を明示的に運び、Web response の metadata に必要な識別情報だけ射影する。
  metadata schema が固定キーならその既存 schema を正しく更新し、勝手に未知 field を
  追加して通す案にはしない。persisted schema の変更が必要と判明した場合は、
  schema namespace・released fixture を確認し、承認済み仕様の scope を越える migration は
  別 Pack とする。
- warning は成功 Result に紐付けて一度表示する。失敗/cancel/stale の候補の warning
  は current Result に混ぜない。全行スキップでも正常な genome diagram を返し、
  「注釈 N 行をスキップ」を隠さない。描画されない mark の legend entry は作らない。
  request に含まれた明示 slot と Web の既存自動 projection の位置・gap は維持し、
  skip を理由に削除・縮小しない。Python が resolved marks から新規の自動 slot を
  作る場合は、空の set に新しい slot を作らない。Web の `automaticCircularAnnotationSlots`
  を生物学的 matching owner にせず、この違いを SEL-03 で確認する。

## 6. BUG-10: canonical caption normalization

対象: `gbdraw/features/colors.py`、`gbdraw/api/prepared.py::resolve_feature_inputs`、
`gbdraw/web_support/rule_matching.py` と既存 typed Worker helper dispatch、
Web `specific-color-rules.js`、`rule-matching.js`、`watchers.js`、
`feature-editor/rule-actions.js`・`color-actions.js` の既存 action 境界。

採用仕様は caption の自動区別。例えば次の入力を受け入れる。

```text
CDS  gene  abc  #112233  Transporter
CDS  gene  mfs  #445566  Transporter
```

正規化した canonical rule/出力 TSV の caption はそれぞれ
`Transporter [#112233]`、`Transporter [#445566]`。各使用色を別の solid 凡例行にする。
同名異色を１つの caption と複数 swatch にまとめる機能は、この修正の対象外。

1. Python に I/O のない `normalize_specific_color_captions` 相当の private helper を
   置き、`resolve_feature_inputs` で **compile より前**に一度適用する。canonical
   table と compiled rule に同じ caption を渡す。legend builder が raw table を
   参照し直さない。Web helper も同じ関数を呼ぶ。
2. 同 caption・同正規化色は同名のまま。異色がある caption は **全色**に normalized
   lowercase hex を付け、first/last-wins を廃止する。空 caption は現在の「凡例を
   作らない」を維持する。rule の match precedence と行順は変えない。
3. 生の caption 全体を先に予約し、自動生成 label が既存の文字列 caption と衝突する
   場合は `(2)` 等を付ける。allocation は caption/color の安定順で決めるが、出力の
   rule 行順は元の順序を維持する。unique rule caption→color を postcondition とし、
   二度適用しても変わらないことを検証する。新しく生成した suffix と他の legend key の衝突を黙って上書きしない。
   既存 `_unique_legend_key` 相当の allocator を再利用し、確定 caption を rule に
   返す。無関係な manual legend との衝突は既存の atomic validation で扱い、勝手に
   置換しない。従来からある type/numeric legend 全般の命名体系はこの修正では変更しない。
4. Web は full candidate rules（file と retained manual の両方）を既存 Worker 経由で
   正規化してから legend intents を作る。正規化後は単色 caption の既存 editor
   モデルを使える。Python が返した normalized rules/intents を使い、JS に caption
   suffix 生成を複製しない。`reject`/`last-wins` の旧分岐は削除する。
5. matching の pattern cache が hit しても caption/color の正規化を skip しない。
   `createRulePreparation` のキャッシュは matching 用に留める。手動 rule の
   追加・変更・feature recolor も同じ preparation action を通す。template/watchers
   に新しい「後で帳尻を合わせる」state mirror を設けない。
6. asynchronous preparation は source file/rules/catalog/Result の snapshot を
   検証し、stale/cancel/error なら mutation しない。current の候補だけを既存 History
   transaction で rules・provenance・legend・Result に一度 commit する。非同期処理前に
   legend だけ変えることは避ける。無関係な manual legend entry は置換しない。
7. Save は admitted canonical rules の文字列を保存し、TSV download も同じ caption を
   出す。原 input の filename/content 自体は変更しない。Session Load は既存 preview と
   draft を保ち自動 Generate しない。過去の曖昧 rule がある draft は新しい rule edit/
   Generate の通常 preparation で明示的に正規化する。migration や dual writer は足さない。
   native RenderResult の request/provenance も実際に使用した normalized rules と一致
   させ、raw caption の request を「実際に使用した値」として保存しない。
8. legend rename/recolor/remove/sort は正規化した該当行だけに作用する。多色化後に
   suffix を文字列から逆解析して source rule を探さない。追加 rule・色変更で新たな
   衝突が起きたら同じ preparation に戻す。

これは suffix を Web だけに付ける修正ではない。native/browser、fresh/live、保存後の
再生成が同じ規則集合を使うまでを１つの色規則修正として扱う。

## 7. BUG-14: pure pixel input contract

対象: Web `track-slot-validation.js`、`linear-track-slots.js`、`circular-track-slots.js`。
Python `gbdraw/tracks/scalars.py` の pixel-text helper と `tracks/circular.py`、
Linear の該当 pixel parse 呼出し。

- scope は Linear `height`/`spacing` と Circular `inner_gap_px`/`outer_gap_px` の
  text/draft/CLI slot/track-table 入口。radius/width の unitless=factor、`%`、typed
  `{value,unit}` の意味を変えない。一般の font-size 等へ規則を広げない。
- 省略/null/trim 後の空文字は auto。数値は有限値だけ。文字列は decimal/exponent と
  optional case-insensitive `px`（前後の空白可）。JS の `Number('')`、`0x10`、bool、
  array、`px`、`10%`、`10em`、Infinity、NaN は許可しない。
- gap/spacing は `>=0`、height は `>0`。受入後の canonical Circular gaps は数値、
  Linear は `{value:number,unit:'px'}`。既存 schema を変えない。
- shared Web helper を既存 `track-slot-validation.js` の純関数として置く。
  両 slot module がすでにこの module を使う方向を維持し、import cycle を作らない。
  `parseOptionalCircularGap`、`parseOptionalLinearPx`、`normalizePxNumberText`、payload
  用の同義 parser を置換する。normalizer は不正値を null/0 にしない。
  編集途中の不正文字列は draft に保持し、row error を出す。submission は同じ
  parser で必ず拒否する。
- Python の text grammar も同じ fixture を使う。typed JSON reader を string 許可に
  緩めず、CLI/TSV text adapter で数値へ変換してから渡す。retired keys の reader は
  増やさない。


## 8. 受け入れ条件

| ID | 受け入れ条件 |
| --- | --- |
| TSV-01 | JS と Python が notes/gene_desc/pmid を持つ同じ TSV を受理し、annotation 値が付加列なしの fixture と同じ |
| TSV-02 | missing required、trim 後 duplicate/BOM、row-width mismatch、known enum/数値/style 不正は fatal。draft/Result 無変更 |
| TSV-03 | unknown の列名と保存されない旨を import notice に表示。render request/Session/TSV export に unknown がない |
| SEL-01 | record ごとに matched/missed の行が混在しても両 mode は成功。他の注釈の geometry は一致 |
| SEL-02 | 複数 selector の部分未一致は行全体を skip、完全一致は元の envelope/segments/circular_path、無効 record は fatal |
| SEL-03 | 全注釈 missing でも通常図が成功し、成功 Result に warning。empty mark の legend を作らず、request に含まれた slot/gap は保つ。resolved mark から新規 auto slot を作る native path は empty slot を作らない |
| SEL-04 | single/grid/batch と Linear の記録 identity、crop/reverse/rotation、index binding、同名 record を検証。余計な resolve と warning 重複なし |
| SEL-05 | CLI/API/Worker warning が同じ skip outcome を示す。cancel/stale/fail、History、Save/Load で別 Result の warning を混ぜない |
| CLR-01 | 同名同色の dedupe、同名異色、named/hex 同色、blank caption、literal suffix 衝突、default/numeric legend 衝突を検証 |
| CLR-02 | normalization は idempotent、行順と match precedence を不変にし、caption-color 対応は rule reorder でも同じ |
| CLR-03 | Web live import／manual edit と native fresh Generate が同じ normalized rules・使用色・caption を持つ。unused-rule 凡例は生成しない |
| CLR-04 | rename/recolor/remove/sort、再 import、Undo/Redo、Save/Load/Generate、dual legend、file/manual 同名と stale/error rollback を検証 |
| CLR-05 | 図の feature colors と各凡例 swatch が一致。元ファイルは保持し、canonical rules の export/provenance に実際の caption を使う |
| PX-01 | 10 / 10px / 10PX / 10 px / .5px / 1e1px は同値。omitted/blank、0、負数、非有限、px、0x10、10%、10em、bool/object を境界別に検証 |
| PX-02 | raw validation・normalized draft・payload・CLI slot/TSV が同じ semantic value。不正入力は null/0 化されない |
| PX-03 | typed gaps の numeric-only、Linear ScalarSpec、Circular factor/%、disabled draft と Save/Load の現行形式を維持 |
| INT-01 | 従来有効入力の SVG が変わらない。新しく受理する多色凡例・skipped annotation のみ期待差分。reference の上書きを回避 |
| INT-02 | browser keyboard/390px、status 通知、実図・download を実機で確認。失敗時は最後の Result を保つ（PD-OI-016） |

対象の既存 Python tests: `test_annotations.py`、`test_annotation_planning.py`、
`test_circular_annotation_tracks.py`、`test_linear_annotation_tracks.py`、
`test_feature_visibility.py`、`test_web_rule_matching.py`、track-slot/session/request tests。
Web: `annotations.test.mjs`、`file-imports.test.mjs`、`track-slot-validation.test.mjs`、
`circular-track-slots.test.mjs`、rule/legend/History/session tests、および該当 Playwright journeys。


## 9. architecture ratchet / rollback

| concern | before → after の owner/path | 除去する重複判断 |
| --- | --- | --- |
| Annotation table | 各 surface codec → 同じ codec、known column projection + notice | strict unknown branch を許可範囲が確定した同じ場所で置換 |
| Selector resolution | Python resolver → 同じ resolver、bundle を plan/drawing/report が共有 | warning 用の別 resolve を作らず、既存 grid 事前解決を plan に収束 |
| Caption normalization | JS reject/last-wins・Python caption overwrite → Python input 正規化を共用 | JS conflict policy と last-wins、raw caption による別 legend semantics |
| px parsing | validator/normalizer/payload の個別 parser → shared lexical helper + 各 typed adapter | 同義 parser と silent invalid→auto/zero |

通常の非増加変更として owner/path evidence を示す。共有 helper は２つ以上の実在経路を
統合するためだけに追加し、旧経路を同じ修正で消す。OE/PE/CB の全 repo 数値を推測しない。
新 compatibility reader、複数 owner の意図的維持、hard invariant waiver が必要になったら
通常ルートを続けず、該当 ratchet の完全 exception evidence と maintainer 判断を用意する。
Product approvalで architecture Gate を waiver できない。

Rollback は関心ごとの runtime revert。承認済み authority を無言で元に戻さず、
未修正状態を回帰として明示し、恒久仕様を戻すなら新 receipt と authority supersession を使う。
既存 canonical schema を変えず、新しい caption も通常の文字列なので speculative reverse
migration は設けない。

## 10. 検証と実装完了の定義

各 session は担当 acceptance に対する positive/negative/parity を focused checks で
確認する。S05 で以下の共有 gates を１回実行し、S05 で code を直した場合だけ影響した
check と必要な regression を再実行する。証拠は source/input/environment/criterion が
変わらない限り再利用する。文書の links/Markdown 変更だけで無変更 runtime の test を
再実行しない。

```bash
ruff check gbdraw/
pytest tests/ -v -m "not slow"
pytest tests/test_output_comparison.py::TestOutputComparison -v
node --test tests/web/architecture-contracts.test.mjs tests/web/architecture-ratchet-fixtures.test.mjs tests/web/product-impact-ratchet-fixtures.test.mjs
```

Web fast contracts / Playwright の正式な scope は `.github/workflows/test.yml` と base CI
impact plan を読む。`.mjs` を production ESM として直接実行する仮定をせず、既存 tests
の disposable-copy loader を使う。Node dependencies は lockfile から当該 clone 内へ
インストールする。新 production dependency を追加しない。

Browser verification が必要なら CLI/ Python Playwright と Node `@playwright/test` の
両方を確認する。Node runner がなければ Python Playwright で同じ targeted assertion。
Chromium sandbox permission failure は同じ check を必要な escalation で再実行する。
Python-only import probe を画面 acceptance の代わりにしない。既存 Gallery-quality
session を使う場合も、その inputs の local copy で実行し、Generator-owned public
artifacts と `examples/gbdraw_social_preview.png` を変更しない。

browser wheel が test に必要な場合だけ `python tools/prepare_browser_wheel.py` で作る。
cache-bust は deployable bundle 準備以外で更新しない。wheel は commit しない。
tracked reference は通常 test で read-only。既存 valid-input の geometry 変更が必要に
なったら原因を説明し、reference を先に更新して差分を隠さない。

trusted-base policy gate は **base の checker** を使い、`--base <実際のauthority済みdev
SHA> --head <実装commit SHA>` を指定する。未 commit diff や candidate checker が
candidate authority を通した結果を、base の検証と取り違えない。mapped behavior
contract の置換が必要なら evidence-only→authority-only→runtime の repository 順序を
適用する。単なる test 更新を authority として扱わない。

テストは30分未満で任意 timeout と判断せず、長い run を増分で監視する。remote CI
の polling は５分より短くしない。gate の failure を修正し、criterion を弱めて通さない。

実装完了は、４ outcome と全 acceptance が満たされ、旧同義 parser/last-wins 等の
superseded paths が削除され、checks と artifacts の根拠が report にあり、全 session
の commit が指定 remote branch に存在する状態。レビュー用 PR の本文が必要になった
場合は `.agents/skills/write-clear-pull-request/SKILL.md` と language check に従う。
public procedural docs/screenshots を編集する場合だけ、対応する docs/Gallery skill を読む。

## 11. References / report convention

- [Product decision records](APPROVED_PRODUCT_DECISIONS.md)
- [Baseline evidence](EVIDENCE.md)
- [Repository guidance](../../../CLAUDE.md)
- [Web guidance](../../../gbdraw/web/CLAUDE.md)
- [Product Impact Ratchet](../PRODUCT_IMPACT_RATCHET.md)
- [Architecture Ratchet](../ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md)
- [Web Change Policy](../WEB_CHANGE_POLICY.md)

`SESSION_RESULTS/` は各 session が自分の完了・境界を記録する場所。最初から空の
completion reports を生成しない。field は Start work SHA / Authority base SHA / Scope /
Owner-path evidence / Commands and results / Artifacts / Remaining work。完成後の own
commit SHA は committed report に自己参照で埋めず、最終 handoff と remote 確認で報告する。
