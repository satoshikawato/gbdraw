# Resume Issue #597 / S08 after S07 (2026-09-28)

gbdraw Issue #597 の S08（Integrated acceptance and final head handoff）を実装・検証してください。計画だけで終えず、許可されたローカル作業を完了してください。BUG-01 には進まないでください。commit、push まで行い、既存の PR #641 を更新してください。remote merge は、下記の gate がすべて満たされた場合に限り行ってください。

## 作業場所と保存先

- 作業 checkout: `/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/issue597-S05-recovered-20260928`
  （branch `fix/issue-597-input-session-20260926`、upstream は同名 remote）。S07 終了時は clean で、local と remote はともに `6c8baf3e8bd50cbab883cc2b96665770e80ac173` でした。別の独立 checkout を使う場合も `/tmp` 以外に置いてください。
- 永続 raw evidence: `/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/issue597-S05-recovery-evidence-20260928/S08/`（新規作成）。S05〜S07 の evidence は同じ親 directory にあります。
- `/tmp` に作業・証拠・生成物を保存しないでください。browser profile、`TMPDIR`、pytest `--basetemp` など大文字小文字を区別する一時領域は、native Linux の専用 directory（例 `/home/kawato/gbdraw-issue597-s08-scratch`）を使ってください。browser profile を `/mnt/c` に置くと、Chromium headless shell は `ILL_ILLOPN`、full Chromium は `SIGTRAP` で起動に失敗します。native Linux に置けば Python Playwright と Node `@playwright/test`（親の shared checkout の `node_modules` から解決、1.61.1）の両方が起動します。
- 共有 checkout と Issue #619 の変更は保全してください。dirty な共有 checkout を reset・clean・switch しないでください。

## 必読資料

最初に次を全文確認してください。
- checkout の `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`
- `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`、`PRODUCT_IMPACT_RATCHET.md`、`WEB_CHANGE_POLICY.md`、`OPTION_INTEGRITY_PRODUCT_CONTRACT.md`（PD-OI-044/045/046 を含む）
- 本 directory の `MASTER_PLAN.md`、`SESSION_WORKFLOW.md`、`sessions/S08_INSTRUCTION_PROMPT.md`
- `results/S03_RESULT.md`、`S04_RESULT.md`、`S05_RESULT.md`、`S05_RECOVERY_20260928.md`、`S05_RESUME_RESULT_20260928.md`、`S05_COMPUTATION_AUDIT_20260928.md`、`S06_RESULT.md`、`S06_FOLLOWUP_RESULT_20260928.md`、`S07_RESULT.md`

承認済み Product outcome（PD-OI-044/045）は再選択を求めずに従ってください。ユーザーが明示的に許容した例外は二つだけです。
- heartbeat max ≤500 ms の未達。
- 実データ Linear の saved preview に非既定の saved config がある場合、その validation のために Python Worker が 1 件起動すること（command-line Session 全般で同じ）。

どちらも閾値は変えず、測定結果は FAIL のまま記録してください。fresh upload の Worker 0 や他の条件へ拡張しないでください。

## intake

- 状態を staged／unstaged／untracked に分類して `S08/intake/` に保存し、actual remote refs と PR #641 の `mergeable`／`mergeStateStatus` を読み取ってください。
- S07 終了時の値は次のとおりですが、latest と決めつけないでください。
  - HEAD／remote: `6c8baf3e`（親 `3f74d507` は `ee6b620` と `0073885` の merge）
  - `origin/dev`: `57cef3ba47f4b7790a9145f2ce6988422a00710e`
  - PR #641: `CONFLICTING`／`DIRTY`
- dev が進んでいる場合は、その差分を読んで統合範囲を判断してください。

## 実施内容

1. **dev 統合（S07 で未実施）**
   - clean tree で `git merge --no-edit origin/dev` を実行してください。統合 merge は実装 commit と分けて記録します。
   - S07 の `git merge-tree` 試算では、次の 6 files が conflict しました（`S07/merge-sim/`）。
     - `gbdraw/web/index.html`
     - `gbdraw/web/js/app/run-analysis.js`
     - `gbdraw/web/js/services/config.js`
     - `gbdraw/web/js/services/error-normalization.js`
     - `tests/web/feature-selection.test.mjs`
     - `tests/web/right-drawer.playwright.spec.js`
   - S03〜S06 の挙動（discovery／disclosure、Session 排他、import Worker、bounded decoding、memory 修正）と、dev 側の Issue #599／#601／#619 の変更を両方保持してください。#619 には `57cef3b` の「current Web Session の Load で保存 active mode を復元」も含まれます。
   - 解消内容と根拠を file ごとに記録してください。mapped contract や guard／authority file を、自己承認のために変更しないでください。
2. **S07 の F1**
   - Circular 探索失敗時に、`Source records` status へ正規化済み error object が raw JSON で表示されます（`S07/observe-followup/discovery-error-{1440,390}.png`）。
   - 正規化された summary を表示するよう修正してください。
   - ファイル名の表示は authority に従ってください：PD-OI-044 と MASTER_PLAN の「安全な filename/stage/error」、PD-OI-046 の「file/record 名を自動公開しない」。両者が両立しない場合は Product preflight とし、outcome を自分で選ばないでください。
   - `tests/web/record-display-discovery.playwright.spec.js` の該当 assertion は、弱めずに authority と一致させてください。
3. **S07 の F4**
   - `session-save-lifecycle.playwright.spec.js` の single-flight test が、圧縮失敗時の戻り値に `error` object がないため fail します。
   - 実装と test のどちらが PD-OI-045/046 に反しているかを判断し、assertion を弱めずに修正してください。
4. **S07 の F5（public docs capture の残り）**
   - T-GUI-05 は stale、T-GUI-06／T-GUI-10 は閉じた disclosure、T-GUI-12 は `Custom Track Slots`、H-GUI-16 は `Import TSV` で止まります。
   - 共有 helper `open_ancestor_details` を使って flow を修正し、`python docs/capture/run_all.py --scenario … --tier extended` で再生成し、`--check` を通してください。
   - 完成図は、label・凡例が preview toolbar に隠れないことを目視で確認してください（S07 では T-GUI-01／09 を 50% にした前例あり。register の規則に従う）。
   - text が変わらない画像でも、手編集はしないでください。source bytes、labels、legend、tracks、comparison context を保持してください。
5. **F2／F3**（dev でも再現）
   - F2：Session 容量制限や未対応 browser の error が汎用 `UNKNOWN` になる。
   - F3：Result がない page の error 後に、status が `Invalid settings · Canonical resource … missing` になる。
   - Issue #597 の authority が修正を要求するかを判断してください。範囲外なら、再現条件付きで別 issue 候補として記録してください。
6. **S08 acceptance matrix**：D-01〜D-04、S-01〜S-05、A-01、W-01 を、実際の tests と results に結び付けてください。未実施・条件一致での reuse・failed を区別してください。次を確認してください。
   - Circular／Linear
   - single／grid／batch
   - native／restored／inactive／settings-only
   - draft と committed の不一致
   - History／rollback
   - full replay
7. **必須 gate**
   - focused Node と full Web Node（30 分以上を許容し、incremental に監視）
   - `pytest tests/ -m "not slow"`（native basetemp）
   - 通常 PR の browser contracts（Node `@playwright/test`、専用 port、`--workers=1`）
   - `TestOutputComparison`（read-only）
   - `ruff check gbdraw/`
   - Vibrio の perf spec と、再構成済み実データ（`reconstructed-fixture/real-full-pairwise.gbdraw-session.json.gz`、SHA-256 `1a89693457e8…`）の Load／Save／Generate
   - architecture checker（`--base origin/dev`、および commit 後の `--head HEAD`）
   - 必要な場合だけ browser wheel を prepare してください。cache-bust は deployable bundle を準備する場合だけ更新してください。

## 守ること

- 閾値、mapped assertions、tests／timeouts、budgets を弱めないでください。
- unchanged evidence の reuse は、source／input／environment／条件が一致するときに限ってください。
- 未測定を PASS と呼ばないでください。heartbeat max 未達、native structured-clone bytes UNAVAILABLE、元 gzip の喪失は未解決のまま記録してください。
- 変更しないもの：`examples/gbdraw_social_preview.png`、`tests/reference_outputs/`、`dist/`、`gbdraw.egg-info/`、browser wheel。
- Gallery session と生成 artifact は owner generator だけで更新してください。
- production・tests・docs・generated の差分は別々にレビューしてください。
- architecture 例外条件が実際に生じた場合だけ、final exact head SHA を入れた review packet を repo 外（`S08/`）の Markdown で作成してください。approval comment は投稿しないでください。
- `pgrep`／`pkill` のパターンが自分の shell command に一致しないよう注意してください。Playwright の route handler は `(route, request)` で呼ばれます。

## commit、push、PR、merge

- 統合 merge と S08 の実装・tests・docs・結果の commit を、同名 remote へ non-force で push し、local と remote の SHA 一致を確認してください。English commit title は `Complete input and Session regression acceptance for issue 597` です。
- PR #641 の本文を更新する場合は `.agents/skills/write-clear-pull-request/SKILL.md` に従い、`node tools/check-pr-language.mjs` を 1 回実行してから `gh pr edit` してください。
- CI の polling は 5 分以上の間隔で行ってください。
- PR #641 の remote merge は、次をすべて満たす場合に限り行ってください。
  - conflict がない
  - required checks（`Web base policy (trusted base)`、`PR / gate`）がすべて成功している
  - 必要な human architecture review が揃っている
  - branch protection を迂回しない

  満たさない場合は merge せず、状態を報告してください。

## 結果

`results/S08_RESULT.md` と `S08/` 配下の raw evidence に、次を保存してください。
- commands と results
- source／input／environment の SHA
- owners と paths
- authority／base references
- acceptance ID
- conflict 解消の記録
- 残る FAIL と未測定、その reuse 条件
- 次の作業の入口

最後に次を報告してください。
- 変更内容、actual SHAs、分類
- 検証結果と evidence paths
- PR／merge の状態
- English proposed commit title と short summary
