# CIテスト信頼性・実行時間改善計画

- Status: in progress (Phase 1-3 complete; Phase 4 corrected after two remote failures; Phase 0 and Phase 5 remote evidence pending)
- Date: 2026-08-10
- Scope: PR #319で確認した成果物比較、Python recipe contract、重複する高コストテスト、GitHub ActionsのPR必須チェック
- Execution: 2026-08-10にPhase 1-3とPhase 4のCI partitionを実装・ローカル検証した。最初のremote run後にT-PY-07上限を300秒へ変更したが、次のremote runでも300秒でtimeoutした。2026-08-11のユーザー指示により、自動testとして固有の価値を持たないT-PY-05/T-PY-07の全量再生成と専用CI jobを削除した。修正後のremote evidenceは未実施。

**訂正 (2026-08-10):** 最新のremote logで、runnerを`ubuntu-24.04`へ変更しても同じrunner familyが選ばれ、`fonts-liberation`も導入済みで0 packageだったため、両変更はno-opだったことが確認された。ここでのH-CLI-13 PNG対策は[H-CLI-13 bundled-font raster implementation plan](HCLI13_BUNDLED_FONT_RASTER_IMPLEMENTATION_PLAN_2026-08-10.md)に置き換えられた。test partitionとruntime改善は引き続き有効だが、runner/font installをremedyまたは成功根拠として扱わない。

**訂正 (2026-08-11):** T-PY-05/T-PY-07は小型のLOSATP、orthogroup、collinearity testが所有する実装契約に、公開成果物の全量再生成を重ねていた。両scenarioは`docs/recipes/run_python_scenarios.py`による明示的な成果物検証として残すが、pytest collectionとPR workflowから外す。timeout延長や専用heavy marker/jobは維持しない。

関連文書:

- [PR 319 CI, Cloudflare, and CodeQL remediation plan](PR_319_CI_CLOUDFLARE_CODEQL_REMEDIATION_PLAN_2026-08-10.md)
- [Documentation simplification implementation plan](DOCUMENTATION_SIMPLIFICATION_IMPLEMENTATION_PLAN_2026-08-09.md)
- [Repository guidance](../../AGENTS.md)
- [Project guidance](../../CLAUDE.md)
- [Failing GitHub Actions run](https://github.com/satoshikawato/gbdraw/actions/runs/31349212637)

## 1. 目的

次の三点を同時に満たす。

1. native Cairo差を公開成果物のstale判定と誤認しない。
2. 各test moduleが固有に所有する契約だけを実行し、全量Tutorial成果物の再生成は明示的な手動検証へ分離する。
3. PR必須チェックの実行開始から最後の必須チェック完了までを、通常時5分以内に収める。

ここでいうPR feedback timeは、最初の必須jobの`run_started_at`から最後の必須jobの`completed_at`までとする。queue待ちは除外し、checkout、依存導入、browser/runtime準備は含める。

目標値:

- 連続5回のPR runの中央値: 5分以内
- 同じ5回の最大値: 7分以内
- rerunで成功したrunは達成根拠に数えない
- `not slow` testのcoverageを落とさず、main向け`slow` gateも維持する

5分を一度だけ下回ることではなく、再実行なしで安定して短いことを完了条件にする。

## 2. 確認済みベースライン

対象run: `31349212637`、head `561ef0d7efd12c0a40ac44d6cdd4d0b420035440`

| Python | pytest時間 | 主な結果 |
| --- | ---: | --- |
| 3.10 | 1,399.23秒（23分19秒） | H-CLI-13失敗、Python recipesが180秒でtimeout |
| 3.11 | 1,180.43秒（19分40秒） | 同上 |
| 3.12 | 1,333.95秒（22分13秒） | 同上 |

workflow全体のwall-clockは約24分40秒だった。三つのPython jobは並列だが、各jobが同じ`not slow` suiteを直列実行するため、最も遅いjobがPR feedback timeを決めている。

Python 3.12 jobの主なfile別時間:

| Test file | 時間 |
| --- | ---: |
| `tests/test_circular_label_placement.py` | 439.2秒 |
| `tests/test_python_howto_recipe_contracts.py` | 180.4秒でtimeout |
| `tests/test_refresh_gallery_sessions.py` | 125.5秒 |
| `tests/test_gallery_session_semantics.py` | 88.3秒 |
| `tests/test_reproduce_examples.py` | 69.7秒 |
| `tests/test_linear_track_layout.py` | 69.0秒 |

現在の`not slow`は軽量gateではない。`slow`で除外されるのは6件だけで、recipe、Gallery再構築、dense label corpus、browser assertionが通常のmatrixに入っている。

### 2.1 H-CLI-13の原因

`docs/recipes/_scenario_support.py::publish_output()`は、`--check`時に生成物と公開物をpayload完全一致で比較する。H-CLI-13のPNGには環境非依存の正規化がない。

制御比較で次を確認した。

- CI: native Cairo 1.18.0
- 公開PNGを再現するlocal環境: native Cairo 1.18.4
- CairoSVG 2.8.2/2.9.0を入れ替えても、native Cairoが同じなら結果は同じ
- 1168×973 RGBAのうち異なるのは71 pixel
- 差分は`x=885..897`、`y=622..632`の13×11領域
- alphaは同一、最大RGB差は21/48/62

したがって原因はPython version、hash seed、CairoSVG version、生成漏れではなく、native Cairoのアンチエイリアス差である。

PNG判定を通した場合、PDFの`/Producer`とEPS/PSの`%%Creator`に残るCairo version差で次に失敗する。現行normalizerはCreationDateだけを除去している。

### 2.2 Python recipe timeoutの原因

`tests/test_python_howto_recipe_contracts.py`が宣言するownerは`H-PY-01`から`H-PY-05`だが、実行commandは`run_python_scenarios.py --all --check`である。

`--all`の実体:

- 15 scenarios
- 21 outputs
- `T-PY-01..09`（10なし）、`T-PY-11`、`H-PY-01..05`

testが結果をassertするのは`T-PY-01`と`H-PY-01..05`だけで、残り9 Tutorialは実行するがassertしない。`T-PY-01`は`tests/test_onboarding_recipe_contracts.py`でも独立に実行される。

同一interpreterでのlocal計測:

| Scenario | 時間 | 主因 |
| --- | ---: | --- |
| `T-PY-07` | 66.53秒 | 13回の逐次LOSATP subprocessが58.71秒 |
| `T-PY-05` | 13.86秒 | 25回のLOSATP search |
| 残り13 scenarios | 約4.2秒 | 通常のparse/render/check |

`H-PY-01..05`だけを既存のper-scenario subprocess形式で実行した場合は20.55秒で完了した。180秒を300秒へ増やすだけでは、owner外の処理と三versionでの重複を残す。

### 2.3 dense label testの原因

`tests/test_circular_label_placement.py`は、同じMjeNMV recordと同じinner-label設定から作る高コストなlabel layoutを複数testで再計算する。Python 3.12では同型のtestが一件約40秒を占める。さらに同じhelperを呼ぶ同値testが二件ある。

production algorithmをCI専用に変える前に、test fixtureのownerと重複assertionを整理する必要がある。

## 3. 守る契約と非目標

### 守る契約

- PRではcore compatibilityをPython 3.10、3.11、3.12で実行する。
- recipe、Gallery、browserの必要なacceptanceはPRではcanonical Python 3.11で一度実行し、mainまたはscheduled gateで3.10/3.12も実行してsupported-version coverageを維持する。
- H-CLI-13はSVG、Interactive SVG、PNG、PDF、EPS、PSの六形式を実際に生成する。
- 公開PNGの寸法、alpha、視覚内容に大きな差があれば失敗する。
- Tutorial recipeは公開されたliteral codeをclean directoryで実行する。
- T-PY-05/T-PY-07の公開成果物を明示的に再検証する場合はprecomputed resultやmockで済ませず、手動scenario runnerで実LOSATPを実行する。
- browser、offline、session、Gallery、reference outputの既存契約を維持する。
- `tests/reference_outputs/`は通常実装で更新しない。

### 非目標

- timeout値を上げて遅さを隠さない。
- 重いtestへ`slow` markerを付けるだけでPR coverageから外さない。
- Cairo 1.18.0でPNGを再生成し、Cairo 1.18.4側を失敗させない。
- 5分達成のためにproductionのLOSATP回数や描画結果を変えない。
- 新しいCI timing service、専用scheduler、独自test runnerを作らない。
- 計測前にGallery/session cacheの汎用層を追加しない。

## 4. 完了条件

### 正しさ

- Cairo 1.18.0と1.18.4で確認した狭いPNG差は許容し、寸法、alpha、広域差、大きなchannel差は拒否する。
- PDF/EPS/PSはCreationDateとrenderer version以外の差を拒否する。
- H-CLI-13の六形式すべてが最後まで検証される。
- H-PY contractは`H-PY-01..05`だけを実行する。
- 自動化するPython Tutorial recipeは別ownerで一度ずつ検証され、T-PY-05/T-PY-07は自動collectionへ含まれない。

### Coverage

- 現在の`pytest tests/ -m "not slow" --collect-only` node ID集合と、PR job全partitionの和集合が一致する。
- partition間の意図しない重複を記録し、同じscenarioの重い再生成は一ownerにする。
- PRで3.11だけを使うacceptance nodeは、mainまたはscheduled gateで3.10/3.12にも割り当てる。
- main branchの`slow` testは従来どおり残す。

### 時間

- `test_python_howto_recipe_contracts.py`のrecipe regenerationはCIで60秒以内。
- MjeNMV inner-label layoutの同一設定は一process内で一度だけ構築する。
- 各PR必須jobはcold setup込み7分未満。
- 連続5回のPR feedback time中央値が5分以内、最大が7分以内。

## 5. Phase 0: baselineとpartition台帳を固定する

Status: pending

### 作業

1. 実装開始時のHEAD、worktree差分、対象run IDを記録する。既存のdirty fileを変更対象へ混ぜない。
2. 現在の全node ID、`not slow` node ID、file別時間を保存する。
3. GitHub Actionsの各step開始・終了時刻、Python/Cairo/CairoSVG/Pillow versionを記録する。
4. CI commandへ`--durations=30`を付け、遅いtest名を毎runで確認できるようにする。新しいreport scriptは作らない。
5. 次のtest ownerを台帳化する。

| Owner | 対象 |
| --- | --- |
| compatibility/core | recipe、Gallery、browser、slow以外 |
| Python how-to recipes | `H-PY-01..05` |
| Python Tutorial recipes | `T-PY-02..04`、`T-PY-06`、`T-PY-08..09`、`T-PY-11`。T-PY-05/T-PY-07は手動成果物検証 |
| onboarding | `T-PY-01`と対応する初回CLI recipe |
| CLI recipes/export | `T-CLI-*`、`H-CLI-*` |
| Gallery/session | bundled Galleryのparse、rebuild、semantic parity |
| browser | Chromiumを必要とするPython/Node/Playwright test |

### Baseline commands

```bash
python -m pytest tests/ -m "not slow" --collect-only -q
python -m pytest tests/ -m "not slow" -q --durations=30
python docs/recipes/run_cli_scenarios.py --scenario H-CLI-13 --check
python docs/recipes/run_python_scenarios.py --all --check
```

### 完了条件

- test node IDの基準集合と、現在の重複ownerが記録されている。
- 以後の時間短縮をtest件数の減少だけで説明できない状態になっている。

## 6. Phase 1: H-CLI-13の環境非依存比較を修正する

Status: completed

Evidence: bounded PNG比較、renderer version正規化、RGBA alpha分離の回帰testを含む対象42件が16.60秒で通過し、H-CLI-13の六形式`--check`が4.49秒で通過した。

### 対象

- `docs/recipes/_scenario_support.py`
- `docs/recipes/run_cli_scenarios.py`
- `docs/capture/run_all.py`
- `tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py`
- `tests/test_documentation_capture_contracts.py`

### 作業

1. `docs/capture/run_all.py::_images_match()`のbounded raster比較を、recipeとcaptureが共有できる小さなhelperへ移す。二つの別実装を作らない。
   Pillowは既存dev dependencyのまま、PNGを`--check`する時だけimportする。package runtime dependencyや`export` extraへ新しい依存を追加しない。
2. RGBA画像はRGBとalphaを分けて比較する。`ImageChops.difference()`をRGBAのまま`getbbox()`へ渡すと、RGB差があってもdifference alphaが0のため見逃すケースを回帰testにする。
3. H-CLI-13 PNGには次の全条件を要求する。
   - width、height、mode一致
   - alpha完全一致
   - changed RGB pixelsが100以下
   - 最大channel差が64以下
   - 差分bounding box面積が256以下
4. 上記thresholdはCairo 1.18.0/1.18.4の観測値だけを覆う初期値とする。別環境で超えた場合は値を無条件に広げず、差分画像とmetricsを確認する。
5. `publish_output()`の比較callbackを、payload変換だけでなくexpected/actualのpairを判定できる形へ置き換える。既存callbackを残したまま第二経路を追加しない。
6. PDFの`/Producer (cairo X.Y.Z...)`、EPS/PSの`%%Creator: cairo X.Y.Z`を正規化する。CreationDate以外のobject/stream/PostScript本文は従来どおり完全一致にする。
7. 比較失敗時はchanged pixel数、最大channel差、bounding boxをerrorへ含める。

### Test matrix

- 同一RGB、同一RGBA: pass
- 71 pixel、13×11、観測済み最大差以内、alpha同一: pass
- 101 pixel: fail
- channel差65: fail
- 小数pixelが広域へ散る: fail
- alpha差: fail
- PDF/EPS/PSのrenderer version/dateだけが異なる: pass
- PDF streamまたはPostScript本文が異なる: fail

### 検証

```bash
python -m pytest \
  tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py \
  tests/test_documentation_capture_contracts.py -q
python docs/recipes/run_cli_scenarios.py --scenario H-CLI-13 --check
git diff --exit-code -- tests/reference_outputs/
```

### 中止条件

- bounded比較がlabel、geometry、色、alphaの実変更を許容する場合は中止する。
- 複数のnative Cairo差を安全に分離できない場合は、fuzzy thresholdを広げず、成果物生成用containerでnative stackを固定する別案へ切り替える。
- PNGの再生成だけで完了にしない。

## 7. Phase 2: Python recipeのownerを分離する

Status: completed, then corrected to remove non-essential full-data pytest owners

Evidence: H-PY、onboarding、Tutorial ownerの対象21件が87.49秒で通過した。H-PY regenerationは5.34秒、全15 scenarioの手動`--all --check`も通過した。

Remote correction: 最初のpartitioned runでは`T-PY-07`がGitHub-hosted runner上で180秒を超え、subprocess上限を300秒へ変更した。しかし次のremote runも300秒でtimeoutし、32 CPUのlocal 60.22秒を4 CPUのhosted runnerへ外挿した前提が誤りだった。T-PY-05/T-PY-07の実装契約は既存のfocused testが所有するため、両scenarioの全量再生成pytest ownerを削除し、timeoutは180秒へ戻した。

### 対象

- `tests/test_python_howto_recipe_contracts.py`
- `tests/test_onboarding_recipe_contracts.py`
- Python Tutorial recipeの新しい単一owner test module
- `docs/recipes/run_python_scenarios.py`（runner変更が必要な場合のみ）

### 作業

1. `test_python_evidence_recipes_regenerate_from_a_clean_external_context`は既存`SCENARIO_IDS`を使い、`H-PY-01..05`を`--scenario ... --check`で個別実行する。
2. `T-PY-01`は既存onboarding testをownerとし、別moduleで再実行しない。
3. T-PY-05/T-PY-07を除くPython Tutorialをscenario単位にparametrizeする。失敗したscenarioがpytest node IDから分かる形にする。
4. T-PY-05/T-PY-07はpytest markerや専用jobを持たず、`--scenario`または`--all`による明示的な成果物検証だけにする。
5. 自動化するTutorialのliteral code、clean-directory実行、成果物freshness、workdir cleanupを弱めない。

### 検証

```bash
python -m pytest \
  tests/test_python_howto_recipe_contracts.py \
  tests/test_onboarding_recipe_contracts.py \
  tests/test_python_tutorial_recipe_contracts.py \
  -q --durations=30
# 公開成果物を明示的に再検証する場合だけ実行
python docs/recipes/run_python_scenarios.py --all --check
```

### 完了条件

- H-PY regeneration testがCIで60秒以内に終わる。
- 13 Python scenariosが自動testで一ownerずつ実行される。
- T-PY-05/T-PY-07はpytest collectionとPR workflowへ含まれない。
- 手動runnerはT-PY-05/T-PY-07の実LOSATP成果物検証能力を維持する。

## 8. Phase 3: 高コストtest fixtureの重複を除く

Status: completed

Evidence: MjeNMV inner-labelの同一設定をmodule-scoped fixtureで一度だけ構築し、同じ構築を使う8件のassertionと`y_overlap` counterで共有した。同値のinner-order testは一ownerへ統合した。対象moduleは68件が83.75秒で通過し、CI baselineの439.2秒から短縮、180秒以内の条件を満たした。

ユーザー指示により30秒閾値は適用せず、Gallery/sessionを含む全suiteを重複監査した。static session decode、SVG parse、figure spec、同一arrow reproductionを既存fixture/cacheへ集約し、Gallery partition全113件が216.05秒で通過した。CIで直接実行するNode specを再実行していた、dataset-specific pytest wrapper 4件を削除した。Node spec 62件は6.98秒、追加のproject-session 3件は3.98秒で直接ownerが通過した。

### 対象

- `tests/test_circular_label_placement.py`
- Phase 3後のprofileで30秒以上の重複が確認されたGallery/session testだけ

### 作業

1. MjeNMVの同一record、config、inner-label layoutをmodule-scoped fixtureで一度だけ作る。
2. `y_overlap` call数も同じfixture構築時に取得し、counter用に同じlayoutをもう一度作らない。
3. consumer testがlabel dictを変更する場合は、fixtureの共有値を変更せず必要なtestだけcopyする。
4. 同じhelperを呼び、同じassertionを行う二つのinner-order testは一ownerへ統合する。
5. outer-only、inner-enabled、middle/resolve-overlapsなどconfigが異なるlayoutは別fixtureのままにする。
6. Phase 3後に再profileし、Gallery/sessionで同一sessionのparse/renderが30秒以上重複することが確認できた場合だけ既存pytest fixtureへ寄せる。汎用cache layerは作らない。

### 検証

```bash
python -m pytest tests/test_circular_label_placement.py -q --durations=30
python -m pytest \
  tests/test_gallery_session_semantics.py \
  tests/test_refresh_gallery_sessions.py \
  tests/test_reproduce_examples.py \
  -q --durations=30
```

### 完了条件

- inner-enabled MjeNMV layoutの同一構築が一process一回である。
- label-placement moduleがCIで180秒以内に終わる。
- assertion coverageをfixture共有だけで減らしていない。

## 9. Phase 4: PR CIを用途別に分割する

Status: locally corrected after two remote failures; remote runtime verification pending in Phase 5

### Marker方針

既存`slow`に加え、次のsuite owner markerだけを追加する。

- `recipe`（公開documentationとCLI/Python recipe acceptance）
- `gallery`
- `browser`

markerなしの`not slow` testはcoreとする。各testはcore、recipe、gallery、browserのうち一ownerだけを持つ。

### 初期job topology

| Job | Python | Selection | Runtime準備 | 目標 |
| --- | --- | --- | --- | ---: |
| `core` matrix | 3.10/3.11/3.12 | core、`not slow` | project/test depsのみ | 各4分以内 |
| `recipes-standard` | 3.11 | `recipe and not slow` | CairoSVG/native Cairo | 3分以内 |
| `gallery` | 3.11 | `gallery` | Galleryに必要なdeps | 4分以内 |
| `browser` | 3.11 | 通常の`browser`とNode/Playwright | Node、Chromium、browser wheel | 5分以内 |
| `losat-cache-browser-acceptance` | 3.11 | 既存の専用browser acceptance | Node、Chromium、browser wheel | 5分以内 |
| `lint` | 3.11 | Ruff | Ruffのみ | 2分以内 |
| `acceptance-supported-main` | 3.10/3.12 | `recipe or gallery or browser` | surfaceごとに必要なruntime | PR必須外、main/scheduled |
| `slow-main` | supported matrix | `slow` | 必要なdeps | PR必須外、mainのみ |

### 作業

1. full `not slow` suiteを三versionで繰り返す現行`test` matrixを、上記partitionへ置き換える。
2. Chromium、Node、browser wheel準備はbrowser jobだけで行う。
3. `libcairo2-dev`はcompileが必要なjobへ限定する。CairoSVG実行だけならrunner既存のruntime libraryを使えるか確認し、不要な`apt-get update`を外す。
4. `actions/setup-python`のpip cacheを有効にする。cache missでも7分上限を満たすことを確認する。
5. 各pytest commandへ`--durations=30`を付ける。
6. まずnative job partitionだけで計測する。core test時間が4分を超える場合に限り、`pytest-xdist`をdev dependencyへ追加し、core jobだけ`-n auto --dist loadfile`を試す。
7. xdistでorder dependencyやresource contentionが出た場合は、そのfileをserial jobへ戻す。timeoutを広げない。
8. PRで3.11だけを使うacceptance selectionを、mainまたはscheduled workflowでは3.10/3.12にも実行する。PR latency改善をsupported-version coverage削減として実装しない。

### Coverage audit

実装後にnode IDを保存し、次を機械的に比較する。

```bash
python -m pytest tests/ -m "not slow" --collect-only -q
python -m pytest tests/ -m "not slow and not (recipe or gallery or browser)" --collect-only -q
python -m pytest tests/ -m "recipe and not slow" --collect-only -q
python -m pytest tests/ -m "gallery" --collect-only -q
python -m pytest tests/ -m "browser" --collect-only -q
```

partitionの和集合が最初の集合と一致しない場合、workflow変更は未完了とする。

2026-08-10のdeduplication後のlocal collectionは次のとおり。owner partition間の重複はなく、件数の和は`not slow`集合と一致した。

| Selection | Node count |
| --- | ---: |
| all | 2,936 |
| `not slow` | 2,930 |
| core | 2,620 |
| recipe standard | 184 |
| Gallery | 113 |
| browser | 13 |

当初の「変更前pytest node ID集合を不変にする」条件は、実装中のユーザー指示「重複testをすべて削除する」を優先した。削除した7 ownerの内訳は、直接Node suiteと同じspecを再実行するpytest wrapper 4件、同じhelper・assertionを持つMjeNMV test 1件、focused実装testへ全量成果物再生成を重ねていたT-PY-05/T-PY-07の2件である。対応behaviorは直接Node owner、focused pytest owner、または明示的なmanual recipe runnerに残る。2026-08-11に再収集した2,930 non-slow nodesは上表のPR partitionで過不足なく所有される。

coreはserial計測が4分を超えたため計画どおり`pytest-xdist`をcoreだけへ追加した。Python compatibilityと無関係なdocumentation、capture、tutorial fixture contract 128件はrecipe acceptanceへ移し、PRではcanonical Pythonで一度だけ実行する。GitHub-hosted runnerのworker数を想定したlocal 4 worker coreは2,603 passed、17 skipped、213.09秒でgreen、Galleryは216.05秒、browser pytestは17.92秒、recipe standardは61.09秒で、test実行部分は各目標内だった。recipe standardのlocal結果は133 passed、46 failedで、失敗はすべて実装開始前から存在するworktreeのCRLF化によるfixture size/checksum差だった。clean checkoutでのrecipe green、setup込み時間、remote timingはPhase 5に残す。

最初のremote runでは、`recipe`へ移したcapture contractが`docs.capture.run_all`経由でPlaywrightをimportする一方、`recipes-standard`が`.[export]`しか導入していないrouting不整合が6件を失敗させた。同じrunのH-CLI-13は、runnerにLiberation Sansがなく別fontへfallbackしたため、PNGが3,541 pixel、最大channel差129、広域bounding box、alpha差ありとなった。`recipes-standard`は既存`.[dev]` extraを使い、`fonts-liberation`を明示導入し、成果物比較環境を`ubuntu-24.04`へ固定した。main向けrecipe acceptanceにも同じfont準備とrunner固定を適用した。

修正後は、以前失敗したPlaywright import 6件が2.60秒で通過した。clean checkout相当ではrecipe standard 179件が65.65秒で通過した。一方、次のremote runでT-PY-07は300秒でも完走しなかったため、全量再生成2件と専用jobを削除した。CIと同じCairoSVG 2.9.0、Pillow 12.3.0、native Cairo 1.18.0とLiberation Sansを使ったH-CLI-13は六形式すべて通過した。workflow YAMLもparse済みである。削除後のremote setup時間と必須check結果はPhase 5に残す。

### 中止条件

- focused testが所有していない固有behaviorをPRから落とす必要がある場合は中止する。
- xdist導入後にflaky failure、port競合、同一fixtureの破壊が出る場合は、plugin設定を増やさずnative shardへ戻す。
- browser testをNode package不在だけでskipしない。必要ならPython Playwrightで同じbehaviorを検証する。

## 10. Phase 5: broad gateとremote実測

Status: pending after two failed remote runs and local remediation

### Local/focused gates

```bash
python -m pytest \
  tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py \
  tests/test_documentation_capture_contracts.py \
  tests/test_python_howto_recipe_contracts.py \
  tests/test_onboarding_recipe_contracts.py \
  tests/test_python_tutorial_recipe_contracts.py \
  tests/test_circular_label_placement.py \
  -q --durations=30

python docs/recipes/run_cli_scenarios.py --scenario H-CLI-13 --check
python -m pytest tests/ -v --tb=short -m "not slow" --durations=30
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
ruff check gbdraw/
git diff --exit-code -- tests/reference_outputs/
```

### Remote evidence

1. PRへpushした最初のrunを記録する。失敗してもrerunせず原因を直す。
2. 全必須jobがgreenになった後、意味のないcommitを追加せず、同じworkflowを`workflow_dispatch`で合計5回計測する。
3. 各runについて次を記録する。
   - run URLとhead SHA
   - 最初のrequired job開始時刻
   - 最後のrequired job完了時刻
   - 各jobのsetup/test時間
   - cache hit/miss
   - top 30 durations
4. 中央値5分、最大7分を満たさない場合、最遅jobの上位testかsetupだけを次の対象にする。job数を無条件に増やさない。

pushと`workflow_dispatch`は外部状態を変更するため、実装セッションで明示的な許可を得てから行う。許可がない場合、Phase 5はremote evidence pendingのままにし、local結果だけで計画をcompletedにしない。

### 完了条件

- focused、full serial、partitioned CIの全gateがgreen。
- reference SVG、公開recipe成果物、Gallery生成物に意図しないdiffがない。
- 連続5回のremote evidenceが時間目標を満たす。
- 計画のStatusと各Phase statusを観測結果に合わせて更新している。

## 11. 変更範囲の見込み

| 区分 | 主なfile | 意図 |
| --- | --- | --- |
| Documentation evidence support | `docs/recipes/_scenario_support.py`, `docs/recipes/run_cli_scenarios.py`, bounded image helper | native renderer差とstale artifactを分離 |
| Recipe tests | Python how-to/onboarding/Tutorial recipe modules、CLI export contracts | scenario ownerを一意化 |
| Performance tests | `tests/test_circular_label_placement.py` | 同一fixture再計算と同値testを除去 |
| Test config | `pyproject.toml` | marker、必要な場合だけxdist |
| CI | `.github/workflows/test.yml` | runtime別partitionと準備の限定 |
| Generated artifacts | 原則なし | semantic output変更時だけauthoritative runnerで再生成 |

production package `gbdraw/`の描画・LOSATP・session実装は、profileで別のproduction問題が証明されない限り変更対象にしない。

## 12. リスクと対策

| リスク | 対策 |
| --- | --- |
| PNG toleranceが実変更を隠す | 寸法・alpha完全一致、pixel数・channel差・bboxの三重上限、reject test |
| Cairo metadata normalizerが本文差を消す | version/dateの狭いregexだけを対象にし、stream/本文差のnegative testを残す |
| scenario分割でcoverageを落とす | collect-only node ID集合の和集合比較 |
| fixture共有でtest間mutationが漏れる | shared valueをread-onlyにし、mutationするconsumerだけcopy |
| xdistでorder dependencyが表面化する | core限定・計測後導入、問題fileはserial shardへ戻す |
| full-data LOSATPのCPU oversubscription | T-PY-05/T-PY-07は自動testにせず、公開成果物を検証する明示的なmanual runnerだけで実行する |
| cache hit時だけ5分を満たす | cache missを含む5回の最大値7分も条件にする |
| job分割がworkflow保守を悪化させる | owner marker三種に限定し、独自sharding scriptを作らない |

## 13. 実装時の進捗チェックリスト

- [ ] Phase 0: baseline、node ID集合、owner台帳を記録
- [x] Phase 1: bounded PNG比較とPDF/EPS/PS normalizerを修正
- [x] Phase 2: H-PYとTutorial recipe ownerを分離
- [x] Phase 3: MjeNMV、Gallery、Node wrapperの重複構築・ownerを除去
- [x] Phase 4: CI partitionとruntime準備を限定（remote timingはPhase 5）
- [ ] Phase 5: focused/full/remote evidenceを取得
- [ ] 連続5 runの中央値5分以内、最大7分以内
- [x] reference/generated artifact差分を個別レビュー（EOL差を除くsemantic diffなし。clean checkoutのexact gateはPhase 5）
- [x] 計画statusとevidence logを更新

## 14. Evidence log template

| Phase | Command/run | Result | Duration | Remaining risk |
| --- | --- | --- | ---: | --- |
| 0 | baseline collection | pending | - | - |
| 1 | `pytest tests/test_cli_tables_tracks_sessions_exports_recipe_contracts.py tests/test_documentation_capture_contracts.py -q`; `run_cli_scenarios.py --scenario H-CLI-13 --check` | 42 passed; six formats verified | 16.60秒; 4.49秒 | remote Cairo 1.18.0実走はPhase 5で確認 |
| 2 | `pytest tests/test_python_howto_recipe_contracts.py tests/test_onboarding_recipe_contracts.py tests/test_python_tutorial_recipe_contracts.py -q --durations=30`; historical `run_python_scenarios.py --all --check`; clean focused `T-PY-07` | historical all 15 scenarios verified;後続remoteでT-PY-07が300秒timeoutしたため全量pytest owner 2件を削除; 残したTutorial 7件はclean cloneで通過 | historical 87.49秒; local T-PY-07 60.22秒; remote T-PY-07 >300秒; retained 7件 8.42秒 | full-data 2件は公開成果物更新時のmanual runnerだけで検証する |
| 3 | `pytest tests/test_circular_label_placement.py -q --durations=30`; `pytest -m "gallery and not slow"` | 68 passed; 113 passed | 83.75秒; 216.05秒 | clean CIでのPython 3.11再測定はPhase 5 |
| 4 | collect-only owner partition; core xdist; recipe standard; browser pytest; Node direct suite; workflow YAML parse | 2,930 non-slow nodesをcore 2,620、recipe 184、Gallery 113、browser 13へ重複・欠落なしで分割; heavy marker/jobを削除; workflow YAML parse通過; clean editable installでrecipe 184件通過 | recipe 92.92秒; focused LOSATP/orthogroup/collinearity 213 passed, 1 skipped in 3.52秒; prior core 213.09秒（4 worker）; browser 17.92秒; Node 6.98秒 | 削除後のremote checkはPhase 5 |
| 5 | remote run 1–5 | pending | - | - |

この計画を実装するセッションは`$execute-plan-with-evidence`を使用し、実測なしでPhaseをcompletedへ変更しない。
