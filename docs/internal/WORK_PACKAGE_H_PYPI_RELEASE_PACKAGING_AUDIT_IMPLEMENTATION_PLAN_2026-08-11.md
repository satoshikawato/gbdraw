# Work package H: PyPI release and packaging audit 実装計画

- Date: 2026-08-11
- Status: planned; implementation not started
- Planning baseline: `docs_renovation` / `1a1780d3a9a7406e74475ffff64761c0fd9cdf2c` と文書作成時点の inventoried working tree
- Source: [`gbdraw_v0.14.0_codex_roadmap.md`](./gbdraw_v0.14.0_codex_roadmap.md) の Work package H
- Target release: `v0.14.0`
- Execution points: H1 before Gate 0; H-final before J-RC and again before J-Final

関連 contract と計画:

- [Repository guidance](../../AGENTS.md)
- [Project guidance](../../CLAUDE.md)
- [Work package A1 plan](WORK_PACKAGE_A1_FINAL_RELEASE_SYNCHRONIZATION_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package I plan](WORK_PACKAGE_I_FEATURE_ANALYTICS_PRIVACY_CONSENT_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package J plan](WORK_PACKAGE_J_QA_COMPATIBILITY_RELEASE_ENGINEERING_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Current installation guide](../INSTALL.md)
- [Current CLI reference](../CLI_Reference.md)
- [Web packaging tests](../../tests/test_web_packaging.py)
- [Browser/offline verifier](../../tools/verify_gui_offline.py)
- [Browser wheel builder](../../tools/prepare_browser_wheel.py)
- [Canonical hosted-bundle candidate](../../tools/prepare_cloudflare_pages.py)

この文書は Work package H を実装可能な変更単位と gate に分ける。文書作成時点では
package code、workflow、generated browser wheel、release artifact、PyPI project、
tag、外部 deployment を変更していない。実装開始時には working tree を再確認し、
無関係な変更を巻き戻さない。

## 1. 結論

`v0.14.0` は一つの `gbdraw` distribution を維持する。ただし、現在の
`py3-none-any` wheel に Linux x86-64 専用 LOSAT executable を含める構成は
universal wheel contract と両立しない。H1 では tracked native binary と bundled
runtime branch を package contract から外し、CLI は explicit executable、`PATH` 上の
`losat`、NCBI BLAST+ `blastp` の順で外部 runtime を解決する。Web は LOSAT Wasm を
引き続き使う。

release build は prepared browser wheel を一度だけ生成して outer wheel と sdist を
build する。その後、accepted outer wheel の nested member を bytes を変えずに staging
へ抽出し、これだけを hosted-bundle builder へ渡す。prepared、outer、sdist、hosted の
四者に同名の異なる wheel を置くことを禁止する。candidate ごとに wheel、sdist、Web
bundle を一度だけ作り、SHA-256 と content manifest で downstream test と publication
を拘束する。

H は artifact を作って検査する。J は candidate を認定する。K だけが明示的な
maintainer authorization の後に tag、TestPyPI/PyPI upload、deploy、archive を行う。
H の test や workflow dry-run は外部 publish の許可にならない。

## 2. 目標、権限境界、対象外

### 2.1 目標

- `pip install gbdraw` を source checkout や editable install に依存しない形で検証する。
- base install と `[export]` を別 environment で検証する。
- universal wheel、sdist、hosted bundle の required/forbidden inventory を機械判定する。
- generated browser wheel が未準備、古い、再帰的、version 不一致なら build を失敗させる。
- wheel、sdist、hosted bundle 間で browser-wheel bytes が一致することを証明する。
- version、metadata、entry point、README rendering、license/notice、size、hash を検証する。
- supported OS/Python matrix で同じ outer-wheel bytes を clean-install する。
- installed package resources から local GUI を起動し、external network blocked で
  Circular、Linear、session replay、export を確認する。
- build/test job と Trusted Publishing job の権限を分離する。
- RC と final の tag/version/index state machine を fail-closed にする。
- H-AC1 から H-AC5 の証拠を Work package J の ledger に渡す。

### 2.2 権限境界

この計画の実装権限に含まれるのは local code、tests、workflow definition、
candidate build、read-only artifact inspection、local dry-run までである。次は
Work package K の個別 authorization が必要になる。

- tag 作成または push;
- TestPyPI/PyPI project 作成、pending publisher 設定、upload;
- GitHub Environment の production approval rule 変更;
- hosted Web deployment;
- GitHub Release、Zenodo、Bioconda、preprint、journal への外部 action。

H1 は external state が未設定でも workflow static contract と non-publishing
candidate run まで完了できる。外部 project ownership は K の entry condition として
記録し、token fallback や個人 credential を repository に追加しない。

### 2.3 対象外

- `gbdraw-core`、`gbdraw-gui`、platform wheel family の新設;
- native LOSAT を platform wheel で再配布する設計;
- formal CycloneDX/SPDX SBOM を新設すること;
- byte-for-byte reproducible double build、`SOURCE_DATE_EPOCH` 導入;
- custom SLSA、cosign、signed checksum、signed tag infrastructure;
- runtime dependency lockfile または全 dependency の upper-bound campaign;
- PyPy、musl/Alpine、BSD、ARM matrix を release blocker にすること;
- full Python suite、full Playwright suite、documentation 全再現を各 OS/Python cell で繰り返すこと;
- Gallery を wheel/sdist に戻すこと;
- `examples/gbdraw_social_preview.png` または reference SVG を packaging の都合で更新すること;
- publish 失敗を `skip-existing` や同一 version の rebuild で隠すこと。

## 3. 確認済み planning baseline

実装開始時に同じ inventory を current HEAD で取り直す。ignored `dist/` や
`gbdraw/web/*.whl` は過去 build を含み得るため、baseline evidence に使わない。

| 領域 | 現状 | 必要な変更 |
|---|---|---|
| Project metadata | `pyproject.toml` は `gbdraw 0.14.0b0`、`requires-python >=3.10`、`py3-none-any` 相当、OS Independent、Python 3.10–3.12 classifiers | tested support と minimum Python を分け、version/index gate を追加する |
| Native LOSAT | `gbdraw/_build_support.py` と `MANIFEST.in` が `gbdraw/bin/linux-x86_64/losat` を wheel/sdist に含める | tracked binary、package-data pattern、bundled resolution branch を除去する |
| Browser wheel | `gbdraw/web/gbdraw-<version>-py3-none-any.whl` は generated かつ gitignored | outer build 前に一度だけ prepare し、未準備なら build 自体を失敗させる |
| Build behavior | `setup.py` は browser wheel が存在する場合だけ package data に含め、存在しなくても build を継続できる | outer `bdist_wheel` と `sdist` に fail-closed validation を置く |
| Built-artifact tests | `tests/test_web_packaging.py` は member inspection を持つが、outer wheel を clean install しない | exact artifact install と capability smoke を追加する |
| Offline GUI | `tools/verify_gui_offline.py smoke-test` は source `WEB_ROOT` を serve する | installed/extracted package resource を明示入力にできるようにする |
| Hosted bundle | Cloudflare builder は browser wheel を自分で再生成し、GitHub Pages workflow は outer wheel を同名で site root に置く | canonical path を一つにし、一度だけ作った browser wheel を明示入力にする |
| Gallery | local package-data と current tests は `gbdraw/web/gallery` を除外する | wheel/sdist では forbidden、hosted bundle では別 inventory とする |
| Version owners | `pyproject.toml` の version/classifier、`gbdraw/version.py`、Web wheel config、`recipe/meta.yaml`、`CITATION.cff`、CHANGELOG/release-note filename/link、notices、tests に分散 | single validator で concordance を強制し、test の hard-coded beta filename を除く |
| CI | `.github/workflows/test.yml` は主に editable install、Linux Python 3.10–3.12 | build-once artifact と 3 OS clean-install matrix を追加する |
| Publishing | TestPyPI/PyPI Trusted Publishing workflow は存在しない | protected manual publisher workflow と static contract tests を追加する |
| Attribution | generated open-source notices は vendored wheel/license inventory を持つ | candidate artifact 内の notice、license、component inventory を検査する |
| Local path leak | `gbdraw/web/wasm/losat/README.md` に maintainer の `/mnt/c/Users/...` path があり、package data に含まれる | repository-relative command に直し、artifact text scan を追加する |
| Public docs | README は Bioconda/source 中心、INSTALL は PyPI route を current のように扱う箇所がある | A1 が candidate state に合わせて同期し、H は artifact evidence を渡す |

## 4. 固定する実装判断

### D1. Distribution identity は一つ

`distribution=gbdraw`、`import=gbdraw`、`command=gbdraw` を維持する。GUI asset の
実測値が Gate 0 budget を満たす限り base wheel に local Web application を残す。
`[export]` は CairoSVG-based output だけを追加する。

### D2. Native LOSAT bundle は repository と package contract から除去する

`gbdraw/bin/linux-x86_64/losat` を削除し、`_NATIVE_RUNTIME_PACKAGE_DATA`、
`recursive-include gbdraw/bin *`、bundled-resource lookup とその専用 test/docs claim を
除く。source checkout と PyPI install の runtime precedence を同じにする。

Native resolution order は次で固定する。

1. explicit `--losatp_bin`;
2. explicit `--ncbi_blastp_bin`;
3. executable `losat` on `PATH`;
4. executable NCBI `blastp` on `PATH`;
5. existing typed error with checked locations and remediation。

LOSAT Wasm は browser artifact であり、native executable 禁止には該当しない。

### D3. Browser wheel は一度だけ生成する

H-final の browser wheel generation は一回である。Prepared file の SHA-256 を記録して
outer wheel と sdist を build し、accepted outer wheel の
`gbdraw/web/<wheel-name>` member を bytes を変えずに staging へ抽出する。Hosted builder
は source-tree 側の file を再利用または再buildせず、この extracted member だけを
validate/copy する。

次の四者の SHA-256 を比較する。

- prepared browser wheel;
- outer wheel の nested member と、その byte-identical extracted file;
- sdist の `gbdraw/web/<wheel-name>` member;
- hosted bundle root の `<wheel-name>`。

同名でも bytes が違えば P0 packaging failure とする。Outer wheel build 後は extracted
nested member が canonical hosted-bundle input であり、mutable source-tree copy ではない。

### D4. Outer build は missing browser wheel を許さない

outer `bdist_wheel` と `sdist` command は
`validate_browser_wheel_prepared()` を最初に実行する。browser-only wheel build は既存の
`GBDRAW_BUILDING_BROWSER_WHEEL=1` でこの outer validation を通らない。通常の source
checkout inspection と editable development を不必要に壊さず、release artifact command
だけを fail-closed にする。

### D5. Artifact policy は machine-readable owner を一つ持つ

`tools/release_artifact_policy.json` を required/forbidden path、magic bytes、
text scan、artifact ごとの compressed/unpacked limit の owner とする。数値は H1 の
fresh baseline measurement 後に review して埋める。`0`、`TBD`、無制限のまま Gate 0 を
通さない。category は必ず測定するが、個別 threshold は実測上 material な category
だけに置く。

`tools/audit_release_artifacts.py` が source/version check、wheel/sdist/Web archive check、
hash、size、browser-wheel identity、evidence JSON を一つの実装で扱う。GUI 固有の
asset semantics は `verify_gui_offline.py` に残し、同じ検査を auditor に複製しない。

### D6. Gallery owner を artifact ごとに分ける

wheel と sdist では `gbdraw/web/gallery/` を一件でも検出したら失敗する。hosted bundle
では Gallery shell、session、example、media、`gallery/remote-assets.json` を required
inventory として検査する。remote URL は full commit SHA に固定し、mutable branch ref を
禁止する。

### D7. Build once、test many

candidate build job は outer wheel、sdist、deployable Web bundle を一度だけ作る。
downstream matrix job は workflow artifact を download し、manifest hash を照合してから
install/test する。sdist qualification のために local derived wheel を作ることは認めるが、
その file は test output であり、release wheel として保存、認定、publish しない。

### D8. Version owner は検証で統合する

`v0.14.0` では dynamic-version refactor を行わない。source owner を手動で更新する既存
構成を維持し、一つの validator で不一致を拒否する。少なくとも次を比較する。

- `pyproject.toml [project].version`;
- `gbdraw.version.__version__` と `__version_display__`;
- `gbdraw/web/js/config.js` の browser-wheel name;
- prepared/nested wheel filename と `METADATA Version`;
- outer wheel/sdist filename と metadata;
- `recipe/meta.yaml` の source version;
- `pyproject.toml` の final/prerelease `Development Status` classifier;
- `CITATION.cff`;
- CHANGELOG heading、`docs/RELEASE_NOTES_<version>.md` filename/heading、
  `docs/DOCS.md` と `docs/ABOUT.md` の release-note links;
- generated open-source notice の project version;
- `gbdraw --version` と `importlib.metadata.version("gbdraw")`;
- release manifest の intended tag。

J より前は actual tag を要求しない。K が tag 作成後に actual
tag/commit/version/hash を検証する。H は `recipe/meta.yaml` を含む source concordance を
検証するが、Bioconda publication action は K が所有する。

### D9. Tested support と minimum Python を区別する

initial tested matrix は既存 CI/classifier と同じ Python 3.10、3.11、3.12、
Linux/macOS/Windows とする。Gate 0 前にそれより新しい stable CPython を probe し、
passing lane を追加するか、release-tested claim の外に置くかを明示する。
CI に合わせるためだけの artificial upper bound は追加しない。

full Python suite は J が Linux の frozen Python matrix で所有する。H は全 OS/Python cell
の lightweight exact-wheel smoke、reference lane の H-AC4、Linux min/max の sdist、
representative lane の H-AC5 を所有する。

### D10. Publisher job は build capability を持たない

publisher workflow は accepted candidate run ID と manifest を入力に exact artifacts を
download する。checkout、build frontend、compiler、test dependency は入れない。
`id-token: write` は literal `testpypi` または `pypi` environment の publisher job
だけに与える。branch push と pull request から upload job へ到達できないようにする。

### D11. PyPI artifact は置換不能

TestPyPI/PyPI へ一度受理された filename/version を再uploadしない。
pre-publication failure は new candidate、post-publication failure は yank assessment と
patch release に戻す。`skip-existing` は workflow と documentation の両方で禁止する。

## 5. Target artifact flow

```text
frozen candidate commit + intended tag + support manifest
    |
    +--> prepare browser wheel once
    |       |
    |       +--> validate name / metadata / contents / non-recursion / hash
    |
    +--> build outer wheel + sdist in empty staging directory
            |
            +--> both embed the exact prepared browser-wheel bytes
            |
            +--> extract nested member from accepted outer wheel
                    |
                    +--> build canonical hosted bundle in empty staging directory
                            |
                            +--> copies the extracted member byte for byte
                            +--> stamps candidate commit/profile in staged output only
                            +--> owns Gallery and immutable remote-assets manifest

staged outer wheel + sdist + hosted tree
    |
    +--> strict metadata + inventory + license + text + size audit
    +--> sorted Web tree manifest + canonical Web archive created once
    +--> release-evidence.json + SHA256SUMS
    +--> immutable workflow artifact
            |
            +--> exact-wheel OS/Python matrix
            +--> exact-sdist qualification
            +--> installed-resource offline GUI
            +--> J-RC or J-Final
                    |
                    +--> K-only authorised publication/deployment
```

`SHA256SUMS` は evidence 用の plain digest list であり、この計画では signature を要求しない。
publisher は list の文字列だけでなく `release-evidence.json` の artifact role、version、
target index、source commit も検証する。

## 6. Artifact contract

### 6.1 Outer wheel

Required:

- `gbdraw` Python packages and declared runtime data;
- `gbdraw` console entry point and correct dependencies/extras metadata;
- config, palettes, required fonts;
- local Web shell, JS modules, Pyodide/vendor assets, tutorial data;
- LOSAT serial/threaded Wasm and cleaned README;
- current open-source notices;
- exactly one non-recursive browser wheel matching project version.

Forbidden:

- ELF, PE, or Mach-O executable payload;
- `gbdraw/bin/`;
- `gbdraw/web/gallery/`;
- `tests/`, `dist/`, `build/`, `*.egg-info`, cache/coverage/editor files;
- private key, token, credential file or matching secret sentinel;
- absolute maintainer path such as `/mnt/c/Users/`, `/home/<user>/`, or
  `C:\Users\...`;
- a second/stale `gbdraw-*.whl` or recursive embedded wheel;
- source map/debug dump not named in policy.

### 6.2 Source distribution

Required:

- build configuration and source needed for isolated wheel construction;
- project license/readme and package data inventory;
- exactly one prepared browser wheel;
- only the tests/tools approved in `release_artifact_policy.json`.

Forbidden:

- native LOSAT or another native executable;
- Gallery;
- previous `dist/`, `build/`, egg-info, caches, credentials, local paths;
- a generated artifact not needed to construct or verify the distribution.

An isolated `pip install <sdist>` must create a functioning local wheel. The derived wheel hash is
recorded under `qualificationOutputs`, never under `releaseArtifacts`.

### 6.3 Hosted Web bundle

Required:

- browser wheel bytes identical to the nested outer-wheel member;
- index/runtime/vendor/Wasm/font/tutorial assets;
- deployment headers and build/profile stamp;
- Gallery shell and declared assets;
- full-SHA remote-assets manifest;
- analytics-disabled or enabled profile exactly matching Work package I Gate 0.

Forbidden:

- outer CPython wheel substituted under the browser-wheel filename;
- stale beta/final wheel beside the intended one;
- mutable Gallery remote ref;
- source/local profile with Google-capable config;
- credential, cache, local path, debug dump, undeclared external application dependency.

### 6.4 Size policy

H1 records, for wheel, sdist, Web archive, and unpacked Web tree:

- total compressed and unpacked bytes;
- Pyodide core/runtime wheels;
- other browser wheels;
- LOSAT Wasm;
- fonts;
- tutorial data;
- Gallery local and remote-declared bytes;
- Python code/data and remaining bytes.

The policy contains a hard compressed/unpacked maximum per artifact. Category measurements are
always recorded; only a category large enough to explain a material artifact regression needs its
own threshold. Any service upload limit remains an outer hard constraint. A budget increase
requires a reviewed policy diff and invalidates prior H-AC1/H-AC3 evidence; the auditor never
warns and continues.

## 7. Planned file ownership

| Owner | Planned files | Responsibility |
|---|---|---|
| Distribution contents | `gbdraw/_build_support.py`, `MANIFEST.in`, `setup.py`, `gbdraw/bin/linux-x86_64/losat` | remove native payload; fail-closed outer build; exact package-data contract |
| Native runtime resolution | `gbdraw/analysis/protein_colinearity.py` | remove bundled branch; keep explicit/PATH fallback and typed error |
| Browser wheel | `tools/prepare_browser_wheel.py`, `gbdraw/web/js/config.js` | prepare once; return/record exact path and hash; no recursive package |
| Hosted bundle | `tools/prepare_cloudflare_pages.py`, canonical deployment workflow selected with Work package I | accept the browser wheel extracted from the accepted outer wheel plus explicit safe output/profile/commit; do not rebuild it |
| Artifact policy/audit | new `tools/release_artifact_policy.json` and `tools/audit_release_artifacts.py` | exact inventory, metadata, hash, size, license, text scan, evidence JSON |
| Installed smoke | new `tools/smoke_installed_distribution.py` and extended `tools/verify_gui_offline.py` | no source import; H-AC4/H-AC5 semantic oracles |
| Package tests | `tests/test_web_packaging.py`, `tests/test_protein_colinearity.py`, new `tests/test_release_artifacts.py` | fail/pass boundaries and archive inspection |
| Workflow tests | new `tests/test_release_workflow_contract.py` | permission, trigger, environment, routing, no-rebuild assertions |
| Metadata | `pyproject.toml` and current version owners | release tool dependency/classifier/version agreement |
| Publishing automation | new `.github/workflows/package_candidate.yml` and `.github/workflows/publish_pypi.yml` | build/test without OIDC; K-only protected publishers |
| Public wording | README, INSTALL, CLI Reference, release notes, CHANGELOG, CITATION through A1 | actual route/support/native-runtime/version claims |
| Evidence | CI `release-evidence.json`, `SHA256SUMS`, this plan ledger, J ledger | candidate identity and gate results; generated evidence is not silently committed |

Do not create a second artifact-inspection library under `gbdraw/`. Release tooling is not runtime
API. Shared browser asset semantics stay in the existing `gbdraw._build_support` and
`verify_gui_offline.py` owners.

## 8. H1 implementation phases

### H1-0. Re-audit and freeze local decisions

Status: pending.

Tasks:

1. Record current branch, full commit SHA, dirty/untracked inventory, Python/pip/build versions.
2. Re-read `pyproject.toml`, `setup.py`, `MANIFEST.in`, package-data helpers,
   release notes, install docs, release workflows, Web package tests, and native runtime resolver.
3. Create a new empty staging root and build a provisional browser wheel, outer wheel, sdist,
   and hosted bundle only to measure baseline. Mark all ignored pre-existing artifacts stale.
4. Record compressed/unpacked/category sizes and fill concrete policy budgets.
5. Probe newer stable CPython interpreters available to CI, then record the proposed frozen
   tested matrix without changing claims prematurely.
6. Record whether the `gbdraw` PyPI/TestPyPI project and Trusted Publisher configuration are
   available. This is a read-only/precondition record; do not create or change external state.
7. Confirm H1 and I1 use one canonical hosted-bundle/deployment path.

Gate:

- baseline commands, tool versions, sizes, path inventories, and unresolved external prerequisites
  are recorded;
- no ignored `dist/` file is treated as candidate evidence;
- policy budgets are numeric and reviewed before H1-2 is accepted;
- any package-name ownership uncertainty is a K prerequisite, not a reason to use a different
  distribution name silently.

### H1-1. Remove the native LOSAT packaging branch

Status: pending.

Production work:

1. Delete `gbdraw/bin/linux-x86_64/losat`.
2. Remove `_NATIVE_RUNTIME_PACKAGE_DATA` and the `gbdraw/bin` manifest include.
3. Remove bundled-resource lookup and bundled-first runtime precedence.
4. Preserve explicit LOSAT, explicit NCBI blastp, PATH LOSAT, PATH blastp, and typed failure.
5. Remove source/PyPI behavioral divergence from CLI documentation through the A1 owner.
6. Keep browser LOSAT Wasm and its notices; rewrite the packaged Wasm README with
   repository-relative build/copy commands.

Focused tests:

```bash
python -m pytest tests/test_protein_colinearity.py -q
python -m pytest tests/test_web_packaging.py -q -k "package_data or wheel or sdist or losat"
```

Required negative oracles:

- package-data patterns and sdist manifest contain no `gbdraw/bin`;
- archive magic-byte scan finds no ELF/PE/Mach-O file;
- source code has no successful bundled-resource branch;
- explicit/PATH fallback order is deterministic on mocked Linux, macOS, and Windows;
- missing runtime error names both supported external remedies;
- browser Wasm tests remain green.

Gate:

- H-AC1 native-executable and fallback portion passes;
- no docs claim still promises a bundled native LOSAT in the PyPI distribution;
- reference SVGs are unchanged.

### H1-2. Make preparation/build/audit fail closed

Status: pending.

Production work:

1. Add outer `bdist_wheel` and `sdist` preflight validation without breaking the
   `GBDRAW_BUILDING_BROWSER_WHEEL` inner build.
2. Make the missing, wrong-name, wrong-version, recursive, or oversized browser wheel fail with
   an actionable message that names `python tools/prepare_browser_wheel.py`.
3. Extend the hosted-bundle command with explicit output, commit, profile, and extracted
   browser-wheel inputs. H-final mode must not rebuild that wheel or read a mutable source-tree
   browser wheel after the outer build.
4. Require the hosted output to be a new or empty child of the recorded staging root. Reject a
   filesystem root, repository root or source-tree path, existing non-empty path, symlink, and any
   path that escapes staging after resolution. Do not apply an unconditional `rmtree` to an
   arbitrary caller-supplied path. Apply the same staging-containment and no-clobber contract to
   auditor extraction, archive, tree-manifest, evidence, and checksum output paths.
5. Add `release_artifact_policy.json` with reviewed numeric budgets and path rules.
6. Implement `audit_release_artifacts.py` with exact-path arguments. It must:
   - require exactly one outer wheel and one sdist;
   - extract exactly one nested browser-wheel member from the accepted outer wheel to a safe
     staging child without recompression or byte changes;
   - call or compose the existing browser-wheel verifier rather than duplicate it;
   - compare prepared/outer-extracted/sdist/hosted browser-wheel hashes;
   - parse wheel/sdist metadata, tags, entry point, requirements, extras, license, README;
   - run required/forbidden path and executable-magic checks;
   - scan text members for local paths, credentials, cache/debug sentinels;
   - calculate compressed/unpacked/category sizes;
   - write the sorted Web tree manifest, create the canonical Web archive exactly once, reopen it
     to verify the member tree against that manifest, hash the literal archive, and forbid any
     downstream re-archive;
   - emit deterministic, sorted `release-evidence.json` and `SHA256SUMS`.
7. Make `twine check --strict` a release-tool dependency and workflow gate.
8. Detect tracked-source changes before/after build. Generated ignored browser-wheel/build output
   and staged bundle stamping are the only allowed changes.

Focused tests:

```bash
python -m pytest tests/test_web_packaging.py tests/test_release_artifacts.py -q
python tools/generate_open_source_notices.py --check
```

Fixture set for auditor tests:

- valid minimal wheel/sdist/Web inventories;
- missing and two browser wheels;
- browser wheel with recursive wheel;
- wrong `Version` or wheel tag;
- ELF/PE/Mach-O sentinel under an innocuous filename;
- Gallery in wheel/sdist;
- outer wheel substituted into hosted root under the expected filename;
- local Windows/WSL/Linux path sentinel;
- credential/cache/debug sentinel;
- budget boundary at exactly limit and one byte above;
- stale second outer wheel in staging;
- mismatched prepared/outer/sdist/hosted browser-wheel hashes;
- any hosted/auditor output set to filesystem root, repository root, source tree, existing
  non-empty path, symlink, or a path that resolves outside the staging root.

Gate:

- direct outer wheel/sdist build fails before preparation and succeeds after preparation;
- exact prepared browser-wheel bytes match outer, extracted, sdist, and hosted destinations;
- auditor failures are hard failures, not warning-only output;
- valid evidence JSON is stable under repeated inspection of the same bytes;
- H-AC1 and the inspection portion of H-AC3 pass.

### H1-3. Add artifact-only clean-install smokes

Status: pending.

Production work:

1. Implement `tools/smoke_installed_distribution.py` without adding the repository root to
   `sys.path`. Give it explicit fixture/output/profile arguments.
2. Run it from a temporary directory with `PYTHONPATH` unset and verify
   `gbdraw.__file__` and `importlib.metadata` resolve inside the candidate environment.
3. Add `basic` and `full` profiles. Basic is fast enough for every OS/Python cell; full owns
   the complete H-AC4 surface on the reference lane.
4. Extend `verify_gui_offline.py` so the Web root can come from installed package resources or
   an extracted artifact. Preserve source-tree mode for development tests.
5. Ensure GUI browser routing aborts and records every external request. Local source/wheel mode
   must make zero external attempts, including forged analytics consent storage.
6. Keep base and `[export]` virtual environments separate.

The smoke script may receive fixture paths from a checked-out test-data directory, but its current
working directory and import path must remain outside the checkout. It must not call an editable
install or set `PYTHONPATH` to the repository.

Focused tests:

```bash
python -m pytest tests/test_release_artifacts.py tests/test_web_packaging.py -q -k "install or smoke or offline"
```

Gate:

- test fails when source-tree import leakage is injected;
- base environment contains no CairoSVG;
- basic profile is stable on Linux, macOS, and Windows runners;
- full profile passes H-AC4 on the declared reference lane;
- separate export/local-GUI profile passes H-AC5 on its declared lanes.

### H1-4. Version contract and release workflows

Status: pending.

Production work:

1. Replace hard-coded `0.14.0b0` artifact filenames in packaging tests with the shared version
   parser unless the test is intentionally an old-version fixture.
2. Add source/candidate version concordance to the artifact auditor, including
   `recipe/meta.yaml`, the final/prerelease `Development Status` classifier, and the release-note
   filename/heading plus links from `docs/DOCS.md` and `docs/ABOUT.md`. H validates these source
   owners; K still owns Bioconda publication.
3. Add `package_candidate.yml`:
   - read-only permissions and no OIDC;
   - clean checkout of the requested commit;
   - prepare browser wheel once;
   - build/audit once;
   - upload exact artifact/evidence bundle;
   - fan out exact-wheel matrix, sdist qualification, H-AC4, H-AC5;
   - aggregate without rebuilding.
4. Add `publish_pypi.yml` with two literal publisher jobs:
   - `publish-testpypi` uses protected `testpypi` environment;
   - `publish-pypi` uses protected `pypi` environment;
   - each is reachable only by explicit K workflow dispatch and matching target/version state;
   - each downloads accepted artifacts by candidate run ID, verifies manifest/hashes, and
     invokes a reviewed, commit-pinned Trusted Publishing action;
   - neither checks out source nor installs build tooling.
5. Add static workflow contract tests for triggers, permissions, environments, conditions,
   artifact download, hash check, absence of build/checkout, and absence of branch/PR upload.
6. Keep external project/environment creation as a documented K prerequisite.

State-machine tests:

| Input | Expected |
|---|---|
| branch/beta version | build/audit only; no publish job |
| intended `vX.Y.ZrcN` + matching prerelease metadata + J-RC manifest | TestPyPI job can await protected approval; PyPI job unreachable |
| intended `vX.Y.Z` + matching final metadata + J-Final manifest | PyPI job can await protected approval; TestPyPI job unreachable |
| tag/version mismatch | fail before OIDC/publish action |
| prerelease routed to PyPI | fail |
| final routed to TestPyPI | fail |
| missing/expired candidate artifact | fail; require new H-final candidate |
| hash mismatch or rebuilt wheel | fail |
| `skip-existing` present | static test fails |

Focused tests:

```bash
python -m pytest tests/test_release_workflow_contract.py tests/test_release_artifacts.py -q
```

Gate:

- H-AC2 passes without contacting an index;
- publisher jobs alone request OIDC and protected environments;
- no workflow path can build and publish in the publisher job;
- intended-tag checks occur before K tag creation, actual-tag checks remain in K.

### H1-5. Non-publishing rehearsal

Status: pending.

Run the complete H flow against the current beta candidate without tag or upload.
Only this rehearsal may use a dirty source tree with an explicit included/excluded inventory.
Such evidence remains `stage=rehearsal` and cannot be promoted to RC/final. Prefer the same clean
detached-checkout path used by H-final whenever practical.

Required rehearsal:

1. Create an empty staging directory.
2. Prepare one browser wheel.
3. Build one outer wheel and one sdist.
4. Extract the browser wheel from the accepted outer wheel and build one canonical hosted bundle
   from that exact staged member.
5. Run strict audit and generate evidence.
6. Download the workflow artifact in every matrix job and verify its hash.
7. Run basic matrix, full H-AC4, sdist, export, and installed-GUI gates.
8. Inspect production/test/workflow/generated-artifact diffs separately.

Gate:

- H-AC1, H-AC2, H-AC4, and H-AC5 pass;
- H-AC3 structure passes but evidence is explicitly labelled rehearsal, not RC/final;
- no external tag, upload, deploy, or environment change occurred;
- H1 is ready for Gate 0 handoff.

## 9. H-final execution

H-final is repeated work, not a one-time implementation phase. Run it once for every candidate
entering J-RC and again after final version/A1 closeout for J-Final.

### 9.1 Entry conditions

H-final is artifact-only. All package/version metadata, A1-owned prose, and generated shipped
files already match the intended tag on the recorded candidate commit; H-final makes no
shipped-source edit.

- Work packages B–I have frozen shipped behavior and files for this candidate.
- Gate 0 records feature scope, schemas, Python/OS/browser support, analytics profile,
  native-runtime policy, artifact inventories/budgets, waiver owner, and intended tag.
- A1 has synchronized candidate-facing public docs and generated assets.
- package/version metadata, the intended tag, A1-owned prose, and generated shipped files already
  agree on the recorded candidate commit.
- execution starts from a clean detached checkout or verified source archive of the recorded
  candidate commit; a shared dirty worktree is never an RC/final build source, even when its
  changes appear unrelated.
- browser wheel, wheel, sdist, and hosted bundle staging directories are new and empty.
- H1 workflow/static/focused gates pass on the same source line.

If an entry condition fails, do not create a release-looking artifact set.

### 9.2 Exact candidate commands

The implemented CLI may refine argument spelling, but it must retain this order and exact-path
discipline:

```bash
python tools/prepare_browser_wheel.py
python tools/verify_gui_offline.py check-assets
python -m build --outdir <staging>/python
python tools/audit_release_artifacts.py extract-browser-wheel \
  --wheel <staging>/python/<wheel> \
  --output <staging>/browser/<expected-browser-wheel>
python tools/prepare_cloudflare_pages.py \
  --output <staging>/web \
  --browser-wheel <staging>/browser/<expected-browser-wheel> \
  --commit-sha <candidate-sha> \
  --profile <disabled-or-enabled>
python -m twine check --strict <staging>/python/<wheel> <staging>/python/<sdist>
python tools/verify_gui_offline.py inspect-wheel <staging>/python/<wheel>
python tools/audit_release_artifacts.py audit \
  --wheel <staging>/python/<wheel> \
  --sdist <staging>/python/<sdist> \
  --prepared-browser-wheel gbdraw/web/<expected-browser-wheel> \
  --extracted-browser-wheel <staging>/browser/<expected-browser-wheel> \
  --web-root <staging>/web \
  --web-archive-out <staging>/release/gbdraw-web-<version>.tar.gz \
  --web-tree-manifest-out <staging>/release/web-tree-manifest.json \
  --intended-tag <tag> \
  --evidence-out <staging>/release-evidence.json \
  --checksums-out <staging>/SHA256SUMS
```

Do not express `<wheel>` or `<sdist>` as `dist/*` in the implementation. Resolve the one expected
filename from the project version, verify the directory has exactly the expected pair, then pass
literal paths to every check.

The auditor validates the Web root, writes its sorted per-file tree manifest, creates the canonical
archive once, then includes both literal outputs in `release-evidence.json` and `SHA256SUMS`. No
later job may re-archive the tree. The deployment job extracts this accepted archive and verifies
that its tree matches the recorded manifest.

### 9.3 H-final RC

Status: pending; blocked on H1, Gate 0, feature freeze, and A1 candidate synchronization.

RC-specific entry prerequisites are prerelease metadata matching intended `v0.14.0rcN`, the
recorded intended tag, and A1 synchronization performed after those owners were updated.

1. Verify those prerequisites on the recorded candidate commit without editing source.
2. Create a new empty staging directory, run H-final, and complete H-AC1 through H-AC5.
3. Hand exact artifacts and evidence to J-RC.
4. Any accepted fix or shipped-file change discards the artifact set, returns to the applicable
   metadata/A1 synchronization owner, and requires a new H-final run.
5. Only after J-RC passes may K request authorization for the actual RC tag, TestPyPI upload,
   and staged deployment.

### 9.4 K-owned TestPyPI/staging observation

This is a mandatory release gate but not an H execution authority.

After an authorised upload, K must:

1. create a clean environment outside the checkout;
2. install the accepted local wheel once with production PyPI as the dependency index;
3. uninstall only `gbdraw`, retaining those production-index dependencies;
4. inspect TestPyPI release metadata and require the accepted wheel and sdist filenames and
   SHA-256 values, with neither file missing or replaced;
5. download `gbdraw==<rc-version>` wheel from TestPyPI with
   `--index-url https://test.pypi.org/simple/ --only-binary=:all: --no-deps --no-cache-dir`
   and no extra index, and retrieve the exact sdist URL reported by the same release metadata;
6. compare both retrieved SHA-256 values, filenames, metadata version, and index response with the
   accepted H-final manifest;
7. install the retrieved wheel with `--no-deps`, run `pip check` and H-AC4;
8. compare staged hosted build stamp/tree hash/profile and run the Work package I privacy smoke;
9. record results in the J/K ledger.

Do not configure TestPyPI as an extra index for dependency resolution. A mixed-index install can
select an unintended package. Any unresolved retrieved-hash, install, H-AC4, hosted-version, or
privacy failure blocks final-version work.

### 9.5 H-final Final

Status: pending; blocked on successful RC observation and final-version/A1 closeout.

Final-specific entry prerequisites are that upstream version owners already contain final
`0.14.0`, the intended tag is already recorded as `v0.14.0`, and A1 has synchronized every
shipped version or release-state owner on that same candidate commit. The actual tag is not yet
required.

1. Verify those prerequisites on the recorded candidate commit without editing source.
2. Create a new empty staging directory and rerun all H-final commands.
3. Produce new wheel/sdist/Web hashes; RC hashes are invalid for final.
4. Rerun H-AC1 through H-AC5 and hand the final set to J-Final.
5. Any change after this point returns to the applicable source owner and then H-final;
   documentation-only is not exempt when it changes
   sdist or hosted bytes.

J-Final pass permits K to request separate authorization for tag, PyPI publication, and deployment.
The publisher/deployer must use the accepted final artifacts rather than rebuilding from the tag.

## 10. Clean-install matrix

Gate 0 may expand this table after the newer-stable-Python probe. A reduction requires a written
support decision before `rc1` and matching metadata/docs. Silent missing jobs are failures.

| Lane | Default scope | Install source | Required oracle | Owner |
|---|---|---|---|---|
| Build/audit | Ubuntu, one frozen build Python | clean candidate | one browser wheel, outer wheel, sdist, Web bundle; H-AC1/H-AC3 | H |
| Base wheel basic | Python 3.10/3.11/3.12 × Ubuntu/macOS/Windows | same exact outer wheel | `pip check`, import/location/version, CLI version/help, small Circular/Linear static SVG | H |
| Base wheel full | Ubuntu, highest frozen Python | same exact outer wheel | all H-AC4 capability oracles | H |
| Sdist min | Ubuntu, Python 3.10 | same exact sdist | isolated derived-wheel build, basic smoke, derived hash recorded | H |
| Sdist max | Ubuntu, highest frozen Python | same exact sdist | isolated derived-wheel build, full or basic smoke per manifest | H |
| Export | Ubuntu, highest frozen Python | exact outer wheel + `[export]` in a new environment | PNG/PDF signatures and semantic dimensions | H |
| Installed GUI | Ubuntu, highest frozen Python, Chromium | exact outer wheel base install | H-AC5 installed-resource/offline flow | H |
| Full Python/core | every frozen Python on Ubuntu | candidate environment selected by J | complete suite | J |
| Browser/offline/security | frozen browser matrix | accepted source/local/hosted artifacts | complete Work package I/J matrix | J |
| TestPyPI observation | Ubuntu, highest frozen Python | retrieved RC wheel `--no-deps` after production-index dependencies | retrieved hash/version/index + H-AC4 | K |
| PyPI post-publication | frozen post-release smoke lane | retrieved final wheel | retrieved hash/version + basic/H-AC4 as release checklist declares | K |

If public documentation promises non-SVG export on macOS or Windows, Gate 0 adds one latest-frozen
`[export]` lane for each claimed OS. H does not infer the claim from the base wheel's
OS-independent tag.

Every artifact-install job must record:

- runner image and architecture;
- Python and pip version;
- requested and resolved dependency list from `pip freeze`;
- candidate artifact filename and pre-install SHA-256;
- `gbdraw.__file__`, package metadata version, CLI version;
- output inventory, assertions, duration, exit status;
- source checkout path and proof that it was absent from `sys.path`/`PYTHONPATH`.

## 11. H-AC4 and H-AC5 smoke specification

### 11.1 Basic profile for every wheel matrix cell

Use small repository-owned fixtures copied into a temporary input directory.

1. Install only the exact base wheel from a literal path.
2. Assert `pip check` succeeds.
3. Assert `import gbdraw`, metadata version, module version, and CLI version agree.
4. Assert the imported module path is under the active environment and not the checkout.
5. Run `gbdraw --help` and both subcommand help paths.
6. Render one small Circular static SVG and one small Linear static SVG.
7. Parse each SVG as XML and assert the expected root, non-empty viewBox, named diagram groups,
   and record identity rather than checking only file existence.
8. Assert output files are confined to the temporary output directory.

### 11.2 Full H-AC4 reference profile

Reuse current repository fixtures instead of adding large new biological data.

| Capability | Candidate fixture | Required oracle |
|---|---|---|
| Package-root Python API | `tests/test_inputs/HmmtDNA.gbk` | `read_genbank` and `draw_circular` return public objects; `to_svg()` parses and contains expected record semantics |
| DDBJ/GenBank | one `AP...` DDBJ accession fixture and `HmmtDNA.gbk` | parser preserves accession/length and renders |
| GFF3 + FASTA | `NC_013668.gff3` + `NC_013668.fasta` | records/features load and both a Python/CLI output assertion pass |
| Circular CLI | `HmmtDNA.gbk` | static SVG with labels and GC/skew semantic groups |
| Linear CLI | existing two-record GenBank pair | both records and expected layout groups |
| Interactive SVG | same small circular record with `interactive_svg` | offline controls/schema are embedded, IDs resolve, no external asset URL |
| Typed request | smallest existing `RenderRequest` recipe | typed render result names/formats and semantic SVG agree |
| Session replay | current small session fixture | load/render/save/reload succeeds; replayed semantic output matches declared oracle |
| Annotation | existing region-annotation TSV/session fixture | expected annotation IDs/marks/legend appear |
| Precomputed comparison | existing `tblastx`/comparison fixture | expected comparison group/count appears and LOSAT subprocess is not called |

Fixture names may be changed only to a smaller existing fixture with equivalent assertions. Record
the exact selected paths and SHA-256 in the smoke manifest. A documentation recipe may be reused,
but the runner must have an installed-artifact mode that does not inject the repository into
`PYTHONPATH`.

### 11.3 H-AC5 export profile

Create a new environment from the exact wheel with `[export]`. Assert:

- base-only negative control cannot claim CairoSVG-backed formats;
- `pip check` succeeds after installing the extra;
- PNG starts with the PNG signature and has non-zero expected dimensions;
- PDF starts with `%PDF-` and has a non-empty page;
- the corresponding SVG semantic source is the declared record;
- no development extra is installed merely to make the test pass.

### 11.4 H-AC5 installed local GUI profile

The browser verifier must prove its served root comes from
`importlib.resources.files("gbdraw") / "web"` inside the installed environment. Then:

1. launch the isolated local server without opening a real user browser;
2. assert COOP/COEP/CORP and CSP headers;
3. block every non-loopback request and retain the request log;
4. set a forged analytics-allowed storage record and prove it cannot enable a local analytics path;
5. wait for Pyodide/browser-wheel readiness;
6. generate one Circular result and one Linear result;
7. save and replay a session;
8. export static SVG and one browser-supported raster/PDF output;
9. assert zero external request attempts, not merely zero successful responses;
10. terminate server/workers and assert no orphan process or leaked object URL remains.

Source-tree smoke remains useful development evidence but cannot substitute for this installed
artifact gate.

## 12. Evidence contract

### 12.1 Machine-readable manifest

`release-evidence.json` has a versioned schema. Minimum fields:

```json
{
  "schemaVersion": 1,
  "candidate": {
    "commit": "<full-sha>",
    "version": "<pep440-version>",
    "intendedTag": "<tag>",
    "stage": "rehearsal|rc|final",
    "worktree": "clean|declared"
  },
  "support": {
    "python": [],
    "operatingSystems": [],
    "browserProfiles": [],
    "analyticsProfile": "disabled|enabled"
  },
  "releaseArtifacts": [],
  "browserWheelIdentity": {
    "preparedSha256": "<sha256>",
    "outerMemberSha256": "<sha256>",
    "extractedSha256": "<sha256>",
    "sdistMemberSha256": "<sha256>",
    "hostedSha256": "<sha256>"
  },
  "webBundleIdentity": {
    "archiveSha256": "<sha256>",
    "treeManifestSha256": "<sha256>",
    "treeHash": "<sha256>"
  },
  "checks": [],
  "matrix": [],
  "qualificationOutputs": [],
  "sourceState": {
    "before": "<inventory-hash>",
    "after": "<inventory-hash>",
    "allowedGeneratedPaths": []
  }
}
```

Each release artifact entry records role, literal filename, media/container type, target, version,
compressed/unpacked size, SHA-256, tree hash where applicable, build command identifier, and policy
version. Every check records stable acceptance ID, command/job, environment, oracle, status, and
evidence path.

`worktree=declared` is valid only when `stage=rehearsal`. The schema validator rejects it for
`stage=rc` or `stage=final`, which require a clean detached checkout or verified source archive
matching `candidate.commit`.

### 12.2 Human-readable ledger

This plan keeps phase status and links to CI runs or local logs. Work package J imports the stable
acceptance IDs and exact artifact hashes into its release ledger. Do not paste multi-megabyte logs
into Markdown; retain them as CI artifacts and record the URL/hash.

Valid statuses:

- `pending`: not run against the required candidate;
- `in progress`: implementation or verification active;
- `passed`: all listed gate evidence exists for the named candidate;
- `failed`: a gate ran and failed;
- `invalidated`: later bytes/source/support policy changed;
- `blocked`: an entry condition outside H is unresolved.

`passed` never means an upload was authorised.

### 12.3 Evidence invalidation

Invalidate at least H-AC1/H-AC3 and every downstream artifact test when any of these changes:

- package version or intended tag;
- package source/data, package-data pattern, build backend/config;
- generated browser wheel;
- Web runtime/vendor/Wasm/font/tutorial data;
- Gallery or hosted-bundle builder/profile/stamp;
- artifact policy or size budget;
- dependency or extra metadata;
- license/readme/notice content;
- workflow artifact packaging.

A test-only change may retain artifact hashes only when it cannot affect build bytes, but the
changed test must rerun against those hashes. Record the justification; do not infer it silently.

## 13. Verification strategy

### 13.1 Focused implementation loop

Run after the owning phase changes:

```bash
python -m pytest tests/test_protein_colinearity.py -q
python -m pytest tests/test_web_packaging.py tests/test_release_artifacts.py -q
python -m pytest tests/test_release_workflow_contract.py -q
python tools/generate_open_source_notices.py --check
ruff check gbdraw/ tools/
```

If `ruff` is not configured for every tool file, use the repository's actual lint scope and record
the release-tool syntax/import checks separately. Do not broaden lint configuration casually in H.

### 13.2 Packaging rehearsal gate

Run only in a fresh staging directory after focused tests:

```bash
python tools/prepare_browser_wheel.py
python tools/verify_gui_offline.py check-assets
python -m build --outdir <fresh-staging>/python
python tools/audit_release_artifacts.py extract-browser-wheel \
  --wheel <fresh-staging>/python/<exact-wheel> \
  --output <fresh-staging>/browser/<expected-browser-wheel>
python tools/prepare_cloudflare_pages.py \
  --output <fresh-staging>/web \
  --browser-wheel <fresh-staging>/browser/<expected-browser-wheel> \
  --commit-sha <candidate-sha> \
  --profile <disabled-or-enabled>
python -m twine check --strict \
  <fresh-staging>/python/<exact-wheel> \
  <fresh-staging>/python/<exact-sdist>
python tools/verify_gui_offline.py inspect-wheel \
  <fresh-staging>/python/<exact-wheel>
python tools/audit_release_artifacts.py audit <exact-arguments>
```

Then run the declared local smoke lanes or the non-publishing CI workflow. Keep staging outside
tracked reference-output directories.

### 13.3 Broader gates before H1 handoff

```bash
python -m pytest tests/ -v -m "not slow"
python -m pytest tests/test_web_packaging.py -v
```

H1 does not claim the full J suite, full browser matrix, or documentation nightly gate. It records
which broader gates J/A1 must still run.

### 13.4 Diff review

Review separately:

1. production/package resolution and deleted binary;
2. build/release tooling;
3. workflow permissions/triggers;
4. tests/fixtures;
5. public documentation/version metadata;
6. generated browser wheel and candidate artifacts.

The generated browser wheel and `dist/` outputs remain ignored and uncommitted. Candidate evidence
must name their hashes even though Git does not track them.

## 14. Risks, stop conditions, and recovery

| Risk | Detection | Response |
|---|---|---|
| Removing bundled LOSAT breaks an advertised native workflow | fallback unit tests or A1 recipe fails | fix explicit/PATH discovery and docs; do not reinsert Linux-only binary into universal wheel |
| Platform wheels appear necessary | required native functionality cannot meet contract externally | stop H, propose a separately approved platform-wheel plan; do not hide platform behavior in `py3-none-any` |
| Browser wheel silently omitted | fail-closed build preflight | run preparation; if generated bytes are wrong, fix preparation rather than weakening preflight |
| Hosted root contains outer wheel under browser filename | cross-artifact hash/inner-inventory mismatch | fail candidate and remove the competing bundle recipe |
| Bundle recipe conflict with Work package I | two production builders or profile mismatch | stop and select one H1/I1 canonical owner before continuing |
| Artifact exceeds reviewed budget or service limit | auditor size failure | remove unintended/duplicate data first; review a GUI split only if measured size remains blocking |
| RC/final starts from a dirty or mismatched tree | source-state/commit mismatch | stop and rebuild from a clean detached checkout or verified source archive; declared dirty state is rehearsal-only |
| Build mutates tracked notices/config | before/after source-state mismatch | regenerate and review the owner before candidate build, commit it through A1, then create a new candidate |
| Matrix imports checkout | module-path/PYTHONPATH oracle | fail job and fix harness; do not accept output parity as proof |
| TestPyPI dependency confusion | mixed index is used or retrieved hash differs | use production-index dependency bootstrap plus TestPyPI `--no-deps` retrieval |
| Workflow artifact expired | publisher cannot retrieve accepted run | create a new H-final candidate and rerun J; never rebuild only in publisher |
| PyPI project/Trusted Publisher unavailable | K prerequisite check | stop external release preparation; do not rename distribution or add API token silently |
| A published RC/final is bad | post-upload smoke | RC blocks final; final triggers documented yank/patch assessment, not overwrite |
| New Python version fails | Gate 0 probe | keep it outside release-tested claim and record issue; do not invent an upper bound without compatibility rationale |
| License/notice gap in vendored asset | strict inventory check | block release until provenance/notice is corrected or the asset is removed |

Before publication, recovery is a new candidate and new hashes. After PyPI publication, package
files are immutable. The release owner decides yank versus patch release and records user impact.
H tooling must preserve the failed artifact/evidence long enough to diagnose it.

## 15. Acceptance mapping

| Acceptance ID | H implementation evidence | J/K follow-up |
|---|---|---|
| H-AC1 distribution contract | H1-1/H1-2 tests, artifact policy report, native fallback tests, size inventory | J checks exact candidate report and docs parity |
| H-AC2 protected workflow | workflow static tests, non-publishing dispatch, permission/trigger report | K configures protected environments and invokes only after authorization |
| H-AC3 exact artifact set | H-final evidence JSON, SHA256SUMS, browser-wheel identity, source-state check | J stores/certifies hashes; K publishes/deploys only those hashes |
| H-AC4 base artifact capability | exact-wheel matrix, full reference profile, sdist qualification | K repeats against TestPyPI RC and post-PyPI route as declared |
| H-AC5 export and installed GUI | separate export env, installed-resource offline browser trace/results | J incorporates browser/security matrix and accepted hosted bundle |

H is not complete merely because `python -m build` succeeds. All five acceptance IDs need named
evidence or an explicit `invalidated`/`blocked` status.

## 16. Status ledger

Initial state at plan creation:

| Phase | Status | Candidate/evidence | Exit condition |
|---|---|---|---|
| H1-0 re-audit and budgets | pending | none | fresh inventory, numeric policy, support/topology decisions |
| H1-1 native LOSAT removal | pending | none | H-AC1 native/fallback checks pass |
| H1-2 fail-closed build/auditor | pending | none | H-AC1 plus H-AC3 inspection pass |
| H1-3 artifact-only smokes | pending | none | H-AC4/H-AC5 pass on declared lanes |
| H1-4 workflows/version state | pending | none | H-AC2 and concordance checks pass |
| H1-5 non-publishing rehearsal | pending | none | complete rehearsal artifact/evidence set |
| Gate 0 handoff | blocked | B–I/support freeze not yet recorded | J support/scope manifest accepted |
| H-final RC | blocked | feature/A1 freeze absent | exact RC artifact set passes H-AC1..5 |
| K TestPyPI/staging observation | blocked | external authorization/J-RC required | retrieved/staged gates recorded |
| H-final Final | blocked | RC observation and final/A1 closeout absent | exact final artifact set passes H-AC1..5 |
| J-Final handoff | blocked | final artifacts absent | J receives final hashes and evidence |

Update this table only from observed evidence. A partial matrix, skipped browser test, stale artifact,
or local source smoke cannot be marked complete.

## 17. Definition of done and handoff

H1 is done when:

1. universal wheel/sdist contain no native executable and native fallback is documented/tested;
2. outer build fails without a valid prepared browser wheel;
3. prepared browser-wheel bytes are identical across outer wheel, its staged extraction, sdist,
   and hosted bundle;
4. machine-readable inventories and numeric budgets reject all declared forbidden cases;
5. exact-wheel clean-install, H-AC4, H-AC5, and sdist qualification work without source leakage;
6. version/tag/index concordance and protected workflow static gates pass;
7. a complete non-publishing rehearsal passes;
8. no tag, upload, deployment, or external state change was performed.

Each H-final run is done when:

1. one exact wheel, one exact sdist, one exact Web archive, its sorted tree manifest,
   `release-evidence.json`, and `SHA256SUMS` exist in a fresh staging set;
2. H-AC1 through H-AC5 pass against those exact bytes;
3. candidate commit, intended tag, support/profile manifest, sizes, hashes, build commands,
   environment, and source-state check are recorded;
4. no later shipped-file or policy change has invalidated the set;
5. the RC/final source was a clean detached checkout or verified source archive matching the
   recorded candidate commit;
6. H-final made no shipped-source edit;
7. J receives the immutable artifacts and evidence without a rebuild or re-archive.

The final implementation handoff must:

- list changed/deleted production, test, workflow, documentation, and generated-artifact owners
  separately;
- include focused and broader command results, matrix job links, skips/retries, and artifact hashes;
- state every remaining K external prerequisite;
- propose an English commit title and short summary;
- make no claim that PyPI, TestPyPI, hosted Web, or a tag is live until K has observed it.
