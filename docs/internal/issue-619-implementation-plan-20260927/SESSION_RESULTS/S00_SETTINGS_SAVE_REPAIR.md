# S00 prerequisite repair — settings-only Session Save

今回の範囲は settings-only Save writer の metadata 投影と受理証拠のみ。
S01 の scalar owner/editor codec、numeric/unit controls、renderer は未着手。

## Checkout・authority・scope

- Branch/upstream: `fix/issue-619-circular-track-measure-inputs` /
  `origin/fix/issue-619-circular-track-measure-inputs`。
- Handoff SHA: `74f9db766708d90e867dbd47fb970cf731b70fb0`。
- Fetched/merged dev: `252986d096011fcf1a0f5564e940480d3b92844d`。
- Merge commit / verified source base: `c75abf72d227443cbc01e0d892e4bc37851b2558`。
  既存公開 branch を merge で継続し、rebase/force push なし。
- Active authority: Product Contract revision 23、`PD-OI-048`〜`050`、
  scenario revision 1、3 Pack の signed A 全文。Owner `satoshikawato`、
  Decision date `2026-09-27`。authority を変更していない。
- 対象範囲の dirty edits はなく、他 checkout の writer/server と対象外の
  8 untracked proposal directories を保持した。検証は専用 snapshot/server/browser。

Developer preflight: `IMPLEMENT_EXISTING_AUTHORITY`。
merged `docs/SESSION_COMPATIBILITY.md` は biological input 前の実 Save と
settings-only Load/native reader を明示的に支持し、inactive biological inputs と
committed artifacts の混入を禁止する。`PD-OI-045` の supported settings-only、
旧 request/resources/Result/History と atomic Load の維持も適用する。
新しい選択肢、型、default、retirement を選ぶ変更ではないため、
`EVIDENCE_REQUIRED` / `PRODUCT_DECISION_REQUIRED` / `NOT_ALLOWED` に該当しない。
実受理証拠はその既存 authority への適合を検証するものであり、authority を代替しない。

## 修復と owner/path

before: `config.js::exportSession` は source-free 状態にも
`runMetadata.annotationWarnings: []` を出力する。
既存 `session-authority.js::validateSettingsOnlyDocument` は非空 metadata を拒否する。
最新 dev merge 後も writer regression で同じ実エラーを再現した。

after: 既存 `settingsOnly` 判定が真の場合だけ `runMetadata: {}` を書く。
他の Session は既存の annotationWarnings / trackSlotGeometry 投影をそのまま使う。
state、biological inputs、committed artifacts は削除せず、validator を変更していない。
Session 44 / request 8 / Circular slot 4 / draft shape と scalar grammar は不変。

Save の owner は `services/config.js`、settings-only 判定/inventory/admission は
`services/session-authority.js`、resource assembly は既存 resource owner、
scalar projection は既存 `app/circular-track-slots.js` / validation owner のまま。
Save→既存 assembly→既存 admission→gzip/download、Load→既存 import transaction
の canonical path は一つのまま。新 owner/path/compatibility branch/registry はない。
changed-scope の OE/PE/CB delta は各 0。superseded path はない。
rollback は writer の1行を revert する（元の Save 不整合は再発する）。

## 実受理

詳細: [settings-save-repair-observations.json](settings-save-repair-observations.json)。
commands/provenance: [settings-save-repair-manifest.json](settings-save-repair-manifest.json)。

- 無編集 fresh control: source/Result/committed request なし。Save Session ボタンの
  gzip download、fresh app の file upload、同じ download の
  `load_session_document(path).to_dict()` が成功。metadata は空。
  空の title を prompt で維持し、無編集条件と History 不変を保った。
- S00 の positive 13 fixtures を width/radius 両方へ適用。
  typed px/factor number、`1.` / `1e-3` 両 unit、precision、legacy bare number/text、
  px suffix、percent、null、blank が実 Save/download/fresh Load/native reader で成功。
  raw number/text/unit と config を完全比較し、既存 canonical payload を
  finite positive number＋px/factor または null と完全比較した。
  `65%` は factor 0.65、bare 1.5 は factor 1.5、Auto は既存 null 意味。
- Rendered control: tobacco Gallery の既存 SVG/geometry と schema-valid な
  synthetic nonempty annotation warning を使用。
  Save/fresh Load/native reader で warning と trackSlotGeometry を維持。
  request、Result bytes の SHA-256、config、biological resources を保持した。
  warning を実 Generate で発生させた証拠ではない。
- inactive biological control: 上記の Circular source を残して Linear へ切替。
  settings-only に分類されず、render request/resources、Circular slot draft、metadata を保存・復元した。
  source resource は公開 seed と byte 一致を別途確認した。
  Load は既存 committed request の Circular mode へ戻るため、inactive 保存時の
  `modeProfiles.activeMode` は linear→circular、adv の axis_stroke_color/evalue/identity/
  pairwise_match_style も対応 profile に変わる。**この control の config 全体一致は未成立**。
  required checkpoint の inactive source 誤分類防止とは区別し、mode/profile 復元を修復・承認しない。
  unchanged inventory Node tests は未生成の inactive source も保護する。
- invalid/unfinished text 5 fixtures: rendered writer 拒否、実 file upload 拒否、
  元の draft/request/Result/History を維持。settings-only の typed incomplete
  text も writer 拒否。null 化や silent fallback はない。
- 全成功例で read/Save は config/request/Result hashes/metadata/biological presence/
  Undo/Redo counts/Worker count を変更しない。
  settings-only fresh Load の Worker construction は 0。
  rendered Load の既存 Worker baseline は保存前と同じ。

既存 settings-only browser spec は非default Circular/Linear profile、auxiliary bytes、
実 Generate→Undo/Redo、failed Load rollback を保護する。既存 Save lifecycle spec は
single-flight、pending、download と各 settlement を保護する。両方を変更せず実行した。

## Commands・environment・証拠

専用 root: `/tmp/gbdraw-619-settings-save-20260927/`。
実 download と negative candidates: `downloads-complete/`。
過去の診断出力は `downloads/` と `logs/` に残す。
旧 S00 の seven evidence hashes は全て unchanged。旧 failure 観測を上書きしていない。
scalar normalization/validation/payload の S00 source hashes も unchanged。
変更された Save 境界の成功は全て今回の source で取得し、旧成功を流用していない。

環境は Python 3.13.3、Node v26.8.2、Python Playwright 1.61.0。
Node `@playwright/test` と Python Playwright の両方を確認した。
Chromium、wheel SHA、runtime/fixture hashes、native import path は観測 JSON に記録。
HEAD の archive に writer の1行差分を適用した dedicated snapshot が対象 checkout と
byte 一致することを照合し、その snapshot 内で wheel を prepare した。
wheel の全 tracked Python files を source と byte 照合した。
他 session の root wheel/server/browser は操作せず、cache-bust を更新していない。

| Command / evidence | Exit / result |
| --- | --- |
| `git fetch origin` / target `git pull --ff-only` / handoff ancestry | 0、up to date / ancestor |
| explicit dev refspec fetch / `git merge --no-edit origin/dev` | 0、上記 merge SHA |
| unchanged-writer new regression | 1、既知の committed render artifacts 拒否を再現 |
| focused Node selection (manifest に全8 files) | 0、20 passed |
| `node --test tests/web/architecture-contracts.test.mjs` | 0、139 passed |
| `pytest -q tests/test_settings_only_session.py` | 0、12 passed |
| snapshot `python tools/prepare_browser_wheel.py --no-build-isolation` | 0 |
| unchanged settings-only / Save lifecycle Playwright specs、専用 config | 0、5 passed |
| `observe_settings_save_repair.py` (下記) | 0、16 positive round trips / 6 negative controls |
| probe `ruff check` / new regression `node --check` | 0 |
| `git diff --check` / ordinary policy gate | 0、Gate PASS / Review REQUIRED（registered session/compatibility path の変更による manual review） |

再現（出力は未使用の専用 directory を指定する）:

```bash
PYTHONPATH=. python docs/internal/issue-619-implementation-plan-20260927/SESSION_RESULTS/observe_settings_save_repair.py \
  --source-root /tmp/gbdraw-619-settings-save-20260927/source \
  --artifacts /tmp/gbdraw-619-settings-save-20260927/downloads-complete \
  --output /tmp/gbdraw-619-settings-save-20260927/recheck.json
```

snapshot 作成手順は manifest、実行した Node browser config は専用 root に保持。
公開後の実 HEAD/remote SHA の照合は `/tmp/gbdraw-619-settings-save-20260927/publication.json`
と最終 handoff に記録する。
probe は過去 progress の自動 reuse を行わない。

診断 failure は pass に含めない。sandbox は host mount error で command 開始前に失敗し、
同じローカル検証を escalation で実行した。new regression / 初回 probe の config 比較は
既存 undefined omission を JSON representation として比較するよう修正した。
初回 probe の新 title 入力は既存 History entry を作るため、無編集 control は空 title を
維持する入力へ修正した。rendered Load に settings-only と同じ Worker 0 条件を誤適用した
probe は、同一 admitted rendered baseline と比較するよう修正した。
別 attempt は fixture 設定中に execution context が切断したため、過去出力を保持して
同じ条件を新 artifact directory で再実行した。runtime/admission/schema/timeout は変更していない。inactive control に settings-only と
同じ全 config 一致を当初要求した probe は、その差を隠さず観測 JSON と本結果に
残し、要求された source 誤分類防止と scalar/request/Result/metadata の保持を独立確認した。

## Review・限界・S01 handoff

production diff は writer の1行。test diff は新 writer regression のみ、
既存 mapped contracts は変更なし。docs diff は本結果と現在の前提条件への追記。
generated diff は新観測 JSON / manifest のみ。旧観測、Gallery/reference outputs、
wheel、dist/egg-info、checker/workflow/Product authority を stage しない。
各区分を別々にレビューし、明示 path だけを stage する。

full pytest、slow/performance、全 browser regression、全 reference comparison、
native render replay、remote CI/staging/deployment は未実行。
1行の metadata projection 修復に対し、変更経路の focused checks と既存 required gates を
実行した。新 UI/IME/unit-change/Auto lifetime/geometry の受入は S01〜S03 の範囲。
S00 の既知の nonfinite-number History loss と typed Boolean coercion は未修復で、
それらを valid numeric draft domain の証拠として扱わない。native invalid-config 受理も
Web admission の根拠にしない。今回の finite valid typed numeric-text の開始条件とは区別する。

S01 の開始条件（3 signed outcomes の base 統合、valid typed numeric-text の
settings-only Save/download/fresh Load/native reader 受理）は成立。
inactive 跨modeの全 config 一致、nonfinite 内部値の History 復元は未成立の観測として
残すが、今回の valid settings-only 表現の開始条件を不成立にするものではない。
S01以降は既存契約どおり、担当変更経路と残る全受入を個別に検証する。
次は [SESSION_PROMPTS/S01_SCALAR_AND_EDITOR_MODEL.md](../SESSION_PROMPTS/S01_SCALAR_AND_EDITOR_MODEL.md)。
この session は S01 を開始しない。

English commit title: `Fix settings-only Session Save metadata`

Summary: `Omit render metadata from settings-only Session exports while preserving rendered metadata; verify real downloads, fresh imports, native readers, and strict draft rejection.`
