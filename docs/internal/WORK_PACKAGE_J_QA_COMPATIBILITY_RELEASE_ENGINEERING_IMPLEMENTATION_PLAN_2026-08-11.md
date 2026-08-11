# Work package J: QA・互換性・release engineering実装計画

- Date: 2026-08-11
- Status: planned; implementation not started
- Planning audit baseline: `docs_renovation` / `6f14e2c4fd2a`; the inventory below
  is historical planning input and Phase 0 must bind execution to a new clean
  integration baseline
- Source: [`gbdraw v0.14.0 Release Roadmap and Codex Implementation Brief`](./gbdraw_v0.14.0_codex_roadmap.md), Work package J
- Target release: `v0.14.0`

関連するcontractと計画:

- [Repository guidance](../../AGENTS.md)
- [Project guidance](../../CLAUDE.md)
- [Web application contract](../../gbdraw/web/CLAUDE.md)
- [Session compatibility history](../SESSION_COMPATIBILITY.md)
- [Current session/request compatibility reference](../REFERENCE/session-and-request-compatibility.md)
- [Work package A1 implementation plan](WORK_PACKAGE_A1_FINAL_RELEASE_SYNCHRONIZATION_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package B implementation plan](WORK_PACKAGE_B_LOSAT_COMPARISON_UI_IMPLEMENTATION_PLAN_2026-08-10.md)
- [Work package C implementation plan](WORK_PACKAGE_C_RECORD_COORDINATE_ANCHORING_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package D implementation plan](WORK_PACKAGE_D_MANUAL_FEATURE_PLACEMENT_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package E implementation plan](RESPONSIVE_PREVIEW_PHASE_1_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package G implementation plan](WORK_PACKAGE_G_ANNOTATION_UX_FINISH_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Work package I implementation plan](WORK_PACKAGE_I_FEATURE_ANALYTICS_PRIVACY_CONSENT_IMPLEMENTATION_PLAN_2026-08-11.md)
- [Web performance remediation plan](WEBAPP_PERFORMANCE_REMEDIATION_IMPLEMENTATION_PLAN_2026-08-10.md)
- [CI reliability plan](CI_TEST_RELIABILITY_AND_RUNTIME_IMPLEMENTATION_PLAN_2026-08-10.md)

Roadmapが参照するA1 planはinitial audit時点では存在せず、J planのdraft中に追加された。
Phase 0でその内容、entry gate、accepted baselineを再確認し、fileの存在だけをA1完了や
J-RCの根拠に数えない。

[Work package H implementation plan](WORK_PACKAGE_H_PYPI_RELEASE_PACKAGING_AUDIT_IMPLEMENTATION_PLAN_2026-08-11.md)
は packaging contract、candidate/publisher workflow、artifact auditor、H-final artifact
production の実装 owner である。J は Phase 0 でその current baseline と evidence contract を
確認し、H が作った exact artifact を consumer/certifier として扱う。J は並行する release
workflow または artifact auditor を新設しない。

この文書はWork package Jを実行可能なrelease gateへ分解する。文書作成時点では
production code、test、workflow、fixture、generated artifact、version、tag、外部release
状態を変更していない。公開、tag、push、deploy、archive、submissionは本計画の実行権限
に含まれない。

## 1. 結論

Jは「最後に不足testを書くwork package」ではない。B–Iが実装と同時に所有testを完成
させ、Jは次の三つを行う。

1. Gate 0でscope、schema、support、artifact、privacy、waiverの境界を固定する。
2. J-RCで一つのcandidate commitと一度だけbuildしたartifact setを、互換性、Web、
   offline、security、platform、documentation、reproducibilityまで含めて認定する。
3. RC後の変更をすべてinvalidate/retestし、final-version candidateをJ-Finalで再認定して
   Work package Kへ渡す。

gate passはreadinessの記録であり、外部actionの許可ではない。Kがtag、TestPyPI、PyPI、
hosted deployment、Bioconda、Zenodo、preprint、journalを行うたびに、明示的なmaintainer
権限を別に得る。

## 2. 目標と対象外

### 2.1 目標

- roadmapのmandatory contractごとに一つ以上のtest/evidence ownerを置く。
- candidate SHA、version、scope manifest、support matrix、artifact hash、test result、
  visual evidence、defect dispositionを一つのledgerで追跡する。
- current writerとsupported readerをnamespaceごとに固定し、first-parent/tag由来fixtureで
  replayを証明する。
- Python 3.10–3.12のcore coverageと、宣言したOS/browserのbuilt-artifact smokeを分離する。
- wheel、sdist、deployable Web bundleをcandidateごとに一度だけbuildし、downstream jobが
  同じbytesを検査・installする。
- base installを`[export]`やdev dependenciesから独立に検証する。
- complete Node Playwright suite、packaged offline cold start、network/privacy canary、
  repeated/cancelled worker lifecycleをrelease gateにする。
- public figureとmanuscript candidate figureをclean archive、declared input、recipe、hash、
  readable-scale visual reviewへ結び付ける。
- P0/P1/P2を実際のstop/promotion/waiver ruleとして運用する。

### 2.2 対象外

- B–Iのproduction featureをJ独自のfallbackで完成させること。
- branch-only schemaを互換formatとして残すこと。
- public Python/CLIにsilent compatibility aliasを追加すること。
- platform wheel family、desktop package、new cloud serviceをv0.14のために新設すること。
- releaseのためだけにreference SVGの比較許容を広げること。
- `examples/gbdraw_social_preview.png`を更新すること。
- unavailable inputを合成fixtureへ置換して、同じpublic figureを再現したと主張すること。
- test timeoutを延ばすだけでrace/flakinessを隠すこと。
- tag、publish、push、deploy、archive、external message、submissionを行うこと。

## 3. 確認済みplanning baseline

Phase 0では同じinventoryを実装開始時のHEADで取り直す。

| 領域 | 現状 | Jへの影響 |
|---|---|---|
| Python test | `tests/test_*.py`は116 module。CIはLinuxでPython 3.10/3.11/3.12のcore、recipe、Gallery、browser、slow partitionを持つ | 強い基盤だが、release artifactではなく主にeditable checkoutをtestしている |
| Node/Web test | 63 Node test moduleと11 Playwright spec | 通常browser jobはgeneral Playwright suiteを実行していない。LOSAT cache specだけ専用runnerがある |
| Browser matrix | `playwright.config.js`はChromium/Desktop Chromeだけ | 対応browserをfreezeするか、Firefox/WebKit smokeを追加する必要がある |
| Rendering reference | 16 tracked SVG referenceとsemantic comparator | default parityの基盤。public/manuscript figureのvisual reviewとは別ownerにする |
| CI trigger | `.github/workflows/test.yml`はbranch push/PR中心でtag gateを持たない | candidate/tag向けrelease workflowが必要 |
| Packaging test | `tests/test_web_packaging.py`はwheel/sdist内容とoffline assetを検査する | built wheel/sdistをclean environmentへinstallしていない。複数箇所が`0.14.0b0`をhard-codeしている |
| Publishing | Trusted Publishing/TestPyPI/PyPI workflowは存在しない | H1でprotected workflowを作り、Jはdry-run/artifact gateだけを認定する |
| Hosted deploy | `.github/workflows/deploy_web.yml`は`main` pushをbuild/deployし、prepared browser wheelではなくouter CPython wheelをsite rootへcopyする | tagged releaseと無関係な後続`main` buildやwrong wheel substitutionをpublication baselineにできない |
| Package truth | metadataはOS Independent、outer wheelは`py3-none-any`、wheel/sdist package dataはLinux x86-64 LOSAT executableを含む | universal distribution contractをGate 0前に解決するrelease blocker |
| Compatibility | session writer 40、reader 27–33/39–40、request writer 5、reader 1/2/5 | C/Dのcomplete next writer前にJはfreezeできない。authentic fixture coverageもmigration familyごとに補う |
| Public contract | package-root APIとCLI snapshotがある | documented `gbdraw.api` export/signature/dataclass/default/error snapshotを追加する |
| Web session UI | current sessionは`ui.zoom`と`ui.canvasPan`を保存する | transaction-only viewport snapshotとsaved navigation policyを混同しない |
| Offline tooling | `tools/verify_gui_offline.py`はasset/wheel inspectionとsource-tree smokeを持つ | extracted installed artifact、mobile、cold/warm repeat、cancelを拡張する必要がある |
| Analytics | roadmapはconsent前no-requestを要求する一方、Cloudflare preparation/testはGA injectionを前提にする | retained analyticsまたはcomplete cutを一つ選び、CSP/offline/local contractを一致させる |
| Version owners | `pyproject.toml`、`gbdraw/version.py`、Web wheel config、`recipe/meta.yaml`、`CITATION.cff`、release notes、notices、testsに分散 | single sourceから検証するconcordance gateが必要。generated file名を手でrenameしない |
| Reproducibility | executable recipe/capture/Gallery harnessがある | planning auditでは4 public Gallery figureがmissing input。claimを維持するならinput/provenanceを解決する |
| Existing artifacts | `dist/`やgitignored browser wheelは過去buildを含み得る | clean isolated output directoryだけをcandidate sourceにする |

## 4. 固定するrelease判断

### D1. Gate passとexternal authorityを分離する

Jのstatusは`not run`、`failed`、`passed`、`invalidated`のいずれかとする。
`passed`はKへのhandoff requestを作れるだけであり、外部actionを開始しない。

### D2. Feature testはshift-leftする

B–Iのacceptance criterionは各work packageがproduction変更と同時にtestする。Jはmissing
ownerを見つけた場合、そのupstream packageをreopenする。J-only mockやmanual checkで穴を
埋めない。

### D3. Release scope manifestを先にfreezeする

少なくとも次をbinary decisionとして記録する。

- `feature_overlap_tolerance_bp` retained/cut;
- `Annotate selection` shortcut retained/cut;
- exact navigation control set;
- slow-preview fallback tuning retained/cut;
- hosted analytics retained/cut;
- supported Python/OS/browser set;
- native LOSAT bundled/external policy;
- public/release/manuscript figure inventory。

retained itemは全surfaceとreplay gateを必須にする。cut itemはUI/API/session/docs/event/CSPから
除去し、absence testを持つ。

### D4. Candidate artifactは一度だけbuildする

build jobがbrowser wheelをrepository toolで準備し、isolated directoryへwheelとsdistを
一度だけ生成する。同じworkflowでdeployable Web bundleも作る。matrix jobはdownloadした
同一artifactを使い、再buildしない。

Prepared browser wheelはouter CPython wheelとは別artifact identityを持つ。Outer wheel内の
nested memberとhosted bundle内のbrowser-wheel fileがprepared browser wheelと同じSHA-256で
あることを証明する。Outer wheelをbrowser-wheel filenameへcopyしてはならない。

### D5. Support matrixはclaimとtestを一致させる

- Linux: Python 3.10、3.11、3.12でcomplete core suite。
- Linux/macOS/Windows: frozen Python matrixでexact built-artifact clean-install smoke。
- Browser: Chromium complete suite。Firefox/WebKitをsupport claimへ含める場合はfocused
  cold-start、worker/Wasm、gzip session、download、privacy smokeを追加する。
- untested later Pythonやbrowserは「supported」と記載しない。Pyodide runtime Pythonと
  native Python supportを別contractとして記録する。

### D6. C/D writerを一回だけadvanceする

request schema 6/session 41がfirst-parent/release evidenceに基づくnext writerのままなら、
CとDのcomplete fieldsを同時にwriteする。実装開始前に別versionがrelease-backedになった
場合は、compatibility evidenceからnext versionを再割当てし、既存versionを再定義しない。

### D7. v0.14 universal distributionからsingle-platform native executableを外す

default implementationは、`py3-none-any`/OS Independent wheelとsdistからLinux-only LOSAT
binaryを除外し、native CLIではdocumented external discoveryを使い、WebではWasmを使う。外部
discoveryがcurrent contractを満たさない場合は作業を停止し、platform wheel設計を別plan
として承認する。single universal distributionの中でplatform branchを隠さない。

### D8. Analyticsはoptional integrationでありruntime dependencyではない

retained時もlocal GUIにはanalytics code/configを注入しない。hosted collector/scriptが
blocked、offline、errorでもapp ready、generate、session、exportは成功する。pre-consent、
reject、revoke後のrequestはP0。cut時はGA injection、CSP allowance、UI、storage、event code、
docs claimをすべて除く。

Enabled profileのreal-tag transport captureとGA property approvalはI3が明示的なstaging権限の
下で所有する。Jはその既存evidenceをexact hosted bundleへ照合するだけで、live endpointへ
test dataを送らない。Gate 0までにI3 evidenceを用意できなければanalyticsをcutする。

### D9. Retry/skipを成功根拠にしない

unexpected skip/xfail、missing job、retry-only pass、manual rerun-only passはrelease evidence
に数えない。Playwright trace、JUnit、skip一覧、retry countをartifactとして残す。

### D10. Candidate byteが変わればevidenceをinvalidateする

RC後のすべてのfixはregression test、affected-evidence list、新candidate hashを持つ。
version-only/docs-only変更もsdist/Web bundleを変えるためJ-Final artifact gateを再実行する。

### D11. A1とH-finalの循環を二passで解く

A1はfeature freeze後にprovisional wheelでdocumentationを同期する。H-finalはその変更後の
exact artifactをbuildする。J-RCがexact artifactに対してrecipe/documentation gateを再実行
する。J-RCのprerequisiteはA1 Phase 0–5のcandidate-synchronization milestoneであり、
post-publication Phase 7を含むA1全体のcompleteではない。Kの外部action後、A1はlive
URL/tag/DOIを検証し、inventory済みのpublication-only ownerだけを更新して
K-Publicationへ渡す。shipped ownerの変更が必要ならpatch candidateとしてA1/H-final/Jを
再実行する。

### D12. Final software releaseとscholarly publicationを分ける

J-Finalはsoftware artifactを認定する。K-Publicationはlive PyPI/Web/Bioconda/Zenodo状態と
preprint/manuscript bundleを認定する。journal submissionをsoftware releaseのopen condition
にしない。

## 5. Evidence contract

### 5.1 Human-readable ledger

この文書のSection 12をstatus ledgerとして更新する。各phaseに次を記録する。

- candidate commit/tree state;
- exact commandまたはhosted job URL;
- environment/tool versions;
- result、duration、skip/retry count;
- inspected artifactとSHA-256;
- visual review対象;
- deviation、P2 waiver、invalidated evidence。

### 5.2 Machine-readable workflow artifact

H1で一つのrelease evidence generatorを実装し、candidateごとに
`release-evidence.json`をCI artifactとして保存する。H plan Section 12.1のversioned schemaが
唯一のownerであり、Jはfieldをrenameした別schemaを作らない。J固有のtraceability fieldは
同じschema versionのreviewed extensionとして追加する。最低限のshapeは次とする。

```json
{
  "schemaVersion": 1,
  "candidate": {"commit": "...", "version": "...", "intendedTag": "...", "scopeId": "..."},
  "support": {},
  "schemas": {},
  "releaseArtifacts": [],
  "browserWheelIdentity": {
    "preparedSha256": "...",
    "outerMemberSha256": "...",
    "extractedSha256": "...",
    "sdistMemberSha256": "...",
    "hostedSha256": "..."
  },
  "webBundleIdentity": {
    "archiveSha256": "...",
    "treeManifestSha256": "...",
    "treeHash": "..."
  },
  "checks": [],
  "matrix": [],
  "qualificationOutputs": [],
  "test_reports": [],
  "reproduction": [],
  "defects": [],
  "sourceState": {}
}
```

artifact entryはfilename、target、version、size、SHA-256、wheel tag、source SHA、build job、
inspection resultを持つ。Repository-tracked public fixtureはstable fixture ID、relative path、hashで
参照してよいが、runtimeで読み込んだsequence、file content、user由来identifierはevidenceへ
入れない。Browser-wheel identityはprepared file、outer nested member、そのstaged extraction、
sdist member、hosted bundle memberのexpected/observed hash relationを記録する。

### 5.3 Requirements traceability

| Source | Jが要求するowner evidence |
|---|---|
| Release contracts 4.1–4.5 | atomic Result、source coordinate、cross-interface parity、no-analysis、local-first privacy |
| Work package B | progressive disclosure、keyboard/focus、inactive draft、no-comparison no-work、viewport |
| Work packages C/D | complete display transform、identity、placement/tolerance、joint writer、all surfaces |
| Work package E | effect registry、copy-on-write projection、single-flight、atomic commit、history、viewport、performance |
| Work packages F/G | frozen navigation set、annotation TSV lossless round-trip |
| Work package H | truthful artifact topology、metadata、clean installs、Trusted Publishing readiness |
| Work package I | retained analytics consent/privacy or complete absence |
| A1 | candidate docs、recipes、captures、Gallery、figures、compatibility text |

### 5.4 K external-action handoff

Jはexternal actionを実行しないが、Section 12.2にK action ledgerのtemplateを用意する。J-Final
時点でexpected commit、intended tag、accepted artifact/bundle hash、required post-action smokeを
埋め、observed URL、published hash、action authority、job/run、resultはWork package Kが実行後に
記録する。K-Publicationはこのledgerが埋まるまでpassしない。

## 6. Target files and ownership

| Area | Primary owners | Planned change |
|---|---|---|
| Roadmap/ledger | `docs/internal/gbdraw_v0.14.0_codex_roadmap.md`, this plan | Gate status, evidence, invalidation, J-to-K handoff and K action template |
| General CI | `.github/workflows/test.yml` | Keep PR/main partitions; expose complete evidence and correct timeouts/skip reporting |
| Release workflow | H-owned new `.github/workflows/package_candidate.yml` and `.github/workflows/publish_pypi.yml` | J consumes the build-once evidence, certifies the exact artifacts, and does not create a parallel workflow |
| Production Web deploy | `.github/workflows/deploy_web.yml`, `tools/stamp_web_build.py`, `tools/prepare_cloudflare_pages.py` | Consume/identify accepted release artifact; no unrelated `main` build as release evidence |
| Browser matrix | `playwright.config.js`, `package.json`, `package-lock.json` | General suite gate, reports/traces, frozen browser projects |
| Package metadata | `pyproject.toml`, `setup.py`, `MANIFEST.in`, `gbdraw/_build_support.py` | Truthful universal wheel and isolated artifact inventory |
| Version owners | `pyproject.toml`, `gbdraw/version.py`, `gbdraw/web/js/config.js`, `recipe/meta.yaml`, `CITATION.cff`, `CHANGELOG.md`, release notes, CLI/reference version claims, notices, intended tag in the manifest | Concordance gate; derive expected wheel names instead of hard-coding beta strings and declare historical-document exceptions |
| Packaging tests | `tests/test_web_packaging.py`, `tests/test_web_runtime_capabilities.py` | Exact content/version/offline assertions; no source-only substitute |
| Release contract test | new `tests/test_release_contract.py` if existing owners cannot hold the checks without duplication | Version/support/workflow/artifact invariants only |
| Compatibility | `gbdraw/session_io.py`, `gbdraw/session_request_codec.py`, Web `config.js`/`session-request.js`, `tests/test_session_*.py` | Final writer/readers, migration/negative/atomic failure |
| Historical fixtures | `tests/fixtures/sessions/`, its `README.md` | Authentic fixture per migration branch with commit/path/hash/oracle |
| Public contracts | `tests/test_public_contract.py`, `tests/fixtures/public_contract.json` | Extend one owner to `gbdraw.api`; do not add a parallel snapshot system |
| Web acceptance | existing focused Node/Playwright specs plus one release acceptance spec only if cross-package flow has no owner | Avoid one giant duplicate UI test |
| Offline/security | `tools/verify_gui_offline.py`, `tests/test_web_packaging.py`, Web sanitization/session tests | Installed/extracted artifact, blocked network, repeat/cancel, mobile, canary |
| Artifact evidence | H-owned new `tools/audit_release_artifacts.py` | J consumes exact-path metadata/content/hash/JSON evidence and adds no second verifier |
| Reproducibility | recipe runners, `docs/capture/run_all.py`, Gallery tools, `tools/reproduce_examples.py` | Clean artifact execution, missing-input disposition, visual evidence |

新規helperは上表のownerを一つにする場合だけ追加する。release logicをworkflow、test、scriptへ
三重実装しない。

## 7. Execution phases

各phaseはproduction、test、workflow、documentation、fixture、generated artifactのdiffを
分けてreviewする。前phaseのgateが通るまで後phaseをcompleteにしない。

### Phase 0: baseline、dependency、scope templateを固定する

Status: pending

#### 作業

1. branch、HEAD、dirty tree、relevant refs、latest release tagを記録する。
2. B–I、performance、CI reliability、A0/A1のstatus/evidenceを再監査する。
3. A1 implementation planのpath、entry gate、accepted baselineを検証する。
4. Python/Node/Playwright test inventory、marker partition、duration、skip/xfail/retryを保存する。
5. current version/schema/support/package/deploy/figure inventoryを保存する。
6. Gate 0 scope/support/compatibility manifest templateをSection 12へ追加する。

#### Baseline commands

```bash
git status --short
git rev-parse HEAD
git describe --tags --always
python -m pytest tests/ --collect-only -q
node --test tests/web/*.test.mjs
npx playwright test --list --project=chromium
python tools/reproduce_examples.py --group all --list
```

#### Completion gate

- baselineは一つのcommit/tree stateに結び付く。
- dependency statusはobserved evidenceと一致し、planned itemをcompleteと数えない。
- scope manifestの未決定欄とdecision ownerが明示される。

#### Stop condition

A1 plan、C/D joint writer、H1 owner、I retain/cut decisionのいずれもowner不在なら、J-RCの
日程を約束しない。automationの独立作業だけを進める。

### Phase 1: upstream requirements-to-test traceabilityを完成する

Status: pending

#### 作業

1. Roadmap 4.1–4.5とB–Iのmandatory acceptance criterionへstable IDを付ける。
2. 各IDへproduction owner、test node/job、fixture、oracle、surface、severityを割り当てる。
3. missing coverageをupstream planへ返し、owner packageで実装する。
4. conditional featureはscope decisionからtest matrixを生成する。
5. public claim、docs scenario、Gallery example、manuscript figureを同じcapability IDへ結ぶ。

#### Completion gate

- mandatory IDにownerなしの行がない。
- retained featureはunit→request/session→render/Web→E2Eの必要層を持つ。
- cut featureはabsence testとclaim removalを持つ。

### Phase 2: compatibilityとpublic contractをfreezeする

Status: pending

#### Owners

- `gbdraw/session_io.py`
- `gbdraw/session_request_codec.py`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-request.js`
- `tests/fixtures/sessions/`
- session/request/public/documentation contract tests

#### 作業

1. C/D complete writerのversionをfirst-parent/release evidenceから確定する。
2. session、request、raw/derived cache、identity manifest、editor metadata、Web file binding、
   table、Gallery、interactive SVG namespaceごとのwriter/reader表を作る。
3. distinct migration branchごとにauthentic fixtureとsemantic oracleを追加する。
4. current load/save/reload、supported migration、future/retired/malformed/truncated/hash mismatch、
   live-state rollbackをPython/Webでtestする。
5. package-rootと`gbdraw.api`、CLI help/default/error snapshotをfreezeする。
6. `ui.zoom`/`ui.canvasPan`のsaved navigation policyとtransaction-only snapshotを文書化する。

#### Focused gate

```bash
python -m pytest \
  tests/test_session_io.py \
  tests/test_session_compat.py \
  tests/test_session_request_codec.py \
  tests/test_api_session.py \
  tests/test_public_contract.py \
  tests/test_documentation_contracts.py \
  tests/test_documentation_reference_contracts.py -v

node --test \
  tests/web/session-request.test.mjs \
  tests/web/gallery-session-migration.test.mjs \
  tests/web/session-authority.test.mjs \
  tests/web/session-resources.test.mjs
```

#### Completion gate

- current Python/Web/docs constantsが一致する。
- positive fixtureはsource commit/path/hashを持つ。
- synthetic version mutationだけでsupportを主張するbranchがない。
- failed importはResult、resources、catalogue、editor、selection、viewportを部分変更しない。

#### Stop condition

partial C-only/D-only writerが一度でもmain/release-backedになった場合、chain migrationを追加
せず、history evidenceに基づくformat remediation planを先に承認する。

### Phase 3: CI topologyとcomplete Playwright gateを作る

Status: pending

#### 作業

1. PR/mainの既存partitionをCI reliability planと照合し、全node IDのownerと重複を確認する。
2. general Node Playwright suiteを少なくともmain/release gateで実行する。
3. Chromium suiteは`workers=1`でrace acceptanceを安定させ、report/trace/retryを保存する。
4. frozen browser supportに応じFirefox/WebKit smoke projectを追加する。
5. job timeoutはrepository guidanceに従い、長いtestを15分で誤殺しない。重いtestはowner別に
   shardし、timeout延長だけで解決しない。
6. unexpected skip/xfail/retry-only passをrelease jobで失敗または明示review対象にする。
7. LOSAT cache専用acceptanceとgeneral suiteの重複ownerを整理する。

#### Gates

```bash
node --test tests/web/*.test.mjs
python tools/prepare_browser_wheel.py
npx playwright test --project=chromium --workers=1
python -m pytest tests/ -m "browser and not slow" --durations=30
python tests/run_losat_cache_browser_acceptance.py
```

#### Completion gate

- 11 existing Playwright specsとnew retained-feature specsがhosted gateでdiscovered/executedされる。
- retry、trace、failure screenshotがrelease evidenceに残る。
- one flaky rerunをpass evidenceに数えない。

### Phase 4: package truth、version concordance、native LOSATを修正する

Status: pending

このphaseのproduction/automation変更はH1が所有する。Jはfocused gateとartifact evidenceを
受け取り、同じrelease invariantを別実装しない。

#### 作業

1. universal wheelとsdistからsingle-platform native LOSATを除外し、external
   discovery/failure messageをLinux/macOS/Windowsでtestする。
2. base install、`[export]`、browser wheel、hosted bundleのcontent boundaryを固定する。
3. `pyproject.toml`をversion sourceとし、import version、CLI、Web wheel、recipe、citation、
   notices、`CHANGELOG.md`、release notes、current CLI/reference docs、wheel/sdist METADATA、
   manifestのintended tagを一つのconcordance testで比較する。Historical release記録の旧versionは
   owner-approved exceptionとしてpath/field単位で列挙する。
4. hard-coded `gbdraw-0.14.0b0-py3-none-any.whl` testをderived expectationへ置換する。
5. package dataのlicense/notices、font、Wasm、browser wheel、tutorial data、Gallery exclusionを
   auditする。
6. Python/OS support metadataとdocumented matrixを一致させる。

#### Focused gate

```bash
python -m pytest \
  tests/test_web_packaging.py \
  tests/test_web_runtime_capabilities.py \
  tests/test_public_contract.py \
  tests/test_release_contract.py -v
```

`tests/test_release_contract.py`を追加しなかった場合は、同じassertionを既存owner testへ
割り当て、存在しないfileをcommandへ残さない。

#### Completion gate

- OS Independent/`py3-none-any` claimとwheel/sdist bytesが一致する。
- base importはCairoSVGやdev dependencyを要求しない。
- version bumpはowner sourceを更新し、generated wheel名を手でrenameしない。
- intended tag、current public version claim、artifact metadataにundeclared mismatchがない。

### Phase 5: build-once release workflowとclean-install matrixを実装する

Status: pending

このphaseのworkflow実装はH1が所有し、Jはreadiness dry runとexact-artifact flowを認定する。
J は H の workflow または auditor と同じ release logic を別名で実装しない。

#### Workflow design

H-owned `.github/workflows/package_candidate.yml` と
`.github/workflows/publish_pypi.yml` は、責務と権限を分離しながら次の一つの logical job
graph を実現する。

```text
source + scope manifest
  -> build wheel/sdist/Web bundle once
  -> inspect + evidence + artifact upload
  -> Linux/macOS/Windows clean-install matrix downloads exact artifacts
  -> offline/browser/reproduction gates download exact artifacts
  -> J decision
  -> protected publish/deploy jobs (K authority only)
```

#### 作業

1. buildはclean isolated output directoryを使い、既存`dist/*`をglobしない。Output countは
   wheel 1、sdist 1、canonical Web bundle 1でなければfailする。
2. `python tools/prepare_browser_wheel.py`を通常artifact build前に一度だけ実行し、そのexact
   pathとSHA-256をmanifestへ記録する。Outer wheel build後、canonical staged-bundle builderは
   accepted outer wheelのnested browser-wheel memberをextractして使い、browser wheelを再build
   しない。Extracted member、prepared input、sdist member、Web bundle memberのhashを照合し、
   outer wheel substitutionをfailする。
3. H-finalではtracked sourceへ`--refresh-cache-bust`を実行しない。H1でcanonical bundle
   builderをstaging-only stampへ変更し、accepted wheel hashまたはcommitから決まる
   deterministic tokenをstaged copyだけへ書く。tracked tokenをownerに残す場合はcandidate
   commit前に更新し、その変更後に全artifactをbuildし直す。
4. exact pathをartifact manifestへ渡し、compressed/unpacked/category size、budget、hash、
   content、METADATA、license/notice resultを記録する。BudgetはGate 0でfreezeする。
5. clean matrixはsource root外のvenvへwheelまたはsdistをinstallし、`pip check`、import、CLI、
   Circular/Linear SVG、session replay、GUI launchを実行する。
6. baseと`[export]`を別environmentにする。
7. publish jobだけにOIDCを与え、TestPyPI/PyPI environment approvalを必須にする。
8. tag/publish conditionをRCとfinalで明示し、workflow_dispatchのreadiness runはpublishしない。

Exact sdistのinstall smokeではinstallerが一時的なderived wheelをbuildしてよい。そのwheelは
qualification outputとしてhashを記録するが、certified release wheelの代替やpublish対象に
しない。source sdistやcertified wheel/Web bundleをdownstream jobで再buildしてはならない。
Build前後のtracked treeを比較し、staging以外の変更があればcandidateをfailする。

#### Local artifact preparation gate

```bash
python tools/prepare_browser_wheel.py
python tools/verify_gui_offline.py check-assets
python -m build --outdir <clean-artifact-directory>
python -m twine check --strict <exact-wheel> <exact-sdist>
unzip -l <exact-wheel>
tar tf <exact-sdist>
python tools/verify_gui_offline.py inspect-wheel <exact-wheel>
```

このgateはexact inventory、forbidden-file scan、prepared/outer/sdist/hosted の
browser-wheel four-way hash、category budgetも実行する。検査は H-owned
`tools/audit_release_artifacts.py` へ集約し、J用の並行 verifierを作らない。
Checksum commandのplatform差はrelease evidence ownerで吸収し、shellごとに別manifestを作らない。

#### Completion gate

- matrix job logが同じwheel/sdist hashを示す。
- prepared browser wheel、outer nested member、そのstaged extraction、sdist member、hosted
  bundle memberのhashが一致し、outer wheelをbrowser assetとしてcopyしていない。
- wheel/sdist/Web bundleのinventory、license/notice、forbidden-file、frozen size budgetがpassする。
- build job以外がcertified wheel、sdist、Web bundleを再buildしない。Exact sdist installが
  作るderived temporary wheelは別hashでqualification outputとして記録する。
- readiness runはexternal upload/deployを一切行わない。
- publish jobはaccepted J evidenceとprotected environmentなしに到達できない。

### Phase 6: packaged Web、offline、security、privacy acceptanceを完成する

Status: pending

このphaseを実行するときは`browser-offline-qa` skillを使い、そのoffline/security/lifecycle
contractに従う。

#### 作業

1. installed/extracted wheelのWeb rootをserveするsmoke modeをoffline verifierへ追加する。
2. Playwright web serverへexplicit artifact rootを渡すrelease modeを追加し、source
   `gbdraw/web` fallbackをrelease jobではfailする。Exact hosted-bundle flowも同じroot contractを使う。
3. cold start、warm repeat、generation、session save/reopen、SVG/PNG/PDF download、cancel、
   supersession、worker disposal、object-URL cleanupをtestする。
4. desktopとfrozen mobile/narrow profileを含める。
5. 全external requestをblock/logし、app-owned assetがsame-originでreadyになることを証明する。
6. unique biological canaryをinput label/sequence/annotation/errorへ配置し、request URL/header/bodyに
   現れないことをassertする。
7. generated/imported/session SVGがshared sanitizerを通り、script/event/foreign/unsafe URLを
   live DOMへ入れないことをtestする。
8. retained analyticsはallow/reject/revoke/storage cleanup/collector failureをintercepted endpointで
   deterministicにtestする。cut時はCloudflare/GitHub bundle、CSP、DOM、storageにabsenceをassertする。
9. source-tree smokeをartifact smokeの代用にしない。

Retained時のreal-tag payload/property evidenceは、I3が別途authorised staging captureで取得した
ものをJが照合する。Jのautomated gateはcollectorをinterceptし、live GA propertyへ送信しない。

#### Gates

Phase 6実装前に存在するcommand/configはinspectionとsource-tree serverだけであり、artifact
acceptanceの代用にしない。Phase 6で`smoke-test --wheel <exact-wheel>`とexplicit Playwright
artifact root、または同等のfail-closed inputを追加した後、target gateは次とする。

```bash
python tools/verify_gui_offline.py inspect-wheel <exact-wheel>
python tools/verify_gui_offline.py smoke-test --wheel <exact-wheel>
GBDRAW_WEB_ROOT=<exact-extracted-bundle-root> npx playwright test --project=chromium --workers=1
```

CLI/env名を実装時に変更した場合はplan、`--help`/config、workflowを同時に更新する。現行の
引数なし`smoke-test`やsource-rooted PlaywrightだけがpassしてもPhase 6はcompleteにしない。

#### Completion gate

- external network 0でもcore workflowが成功する。
- retained analyticsだけがconsent後にbounded requestを行い、failureはappへ影響しない。
- local GUIではanalytics injection/requestが0。
- repeated/cancelled run後も前runのResult/resource/worker stateが混入しない。

### Phase 7: documentation、Gallery、public/manuscript figure reproducibilityを閉じる

Status: pending

#### 作業

1. A1でpublic docsをfeature-complete codeへ同期する。
2. exact H-final wheelをclean directoryへinstallし、literal CLI/Python recipeを再実行する。
3. docs capture、Gallery session regeneration、tutorial captureをgenerator ownerから実行する。
4. `tools/reproduce_examples.py --group all --list`のmissing inputを解消するか、public/release claimを
   正直に変更する。
5. release/manuscript figure manifestにrecipe/session、input hash/license、tool versions、expected
   semantic/hash evidence、visual reviewを記録する。
6. `git archive`相当のclean candidate treeから再現し、local-only path/cacheへ依存しないことを確認する。
7. final artifactをreadable scaleで目視する。minimal smoke diagramをpublic figureに昇格しない。
8. protected social previewはhash確認だけ行い、再生成しない。

Exact-artifact runはsource checkoutをimport pathにしないclean venvで行い、`gbdraw.__file__`と
distribution versionがinstalled wheelを指すことを最初にassertする。Recipe/capture ownerが
repository metadataを必要とする場合はclean archive copyをfixture/runnerとして使うが、package
importはそのcopyではなくinstalled artifactから解決する。Source checkout上のpassだけを
exact-artifact evidenceに数えない。

#### Gates

```bash
python docs/recipes/run_cli_scenarios.py --all --check
python docs/recipes/run_python_scenarios.py --all --check
python docs/capture/run_all.py --tier core --check
python docs/capture/run_all.py --tier extended --check
python tools/capture_gallery_tutorial_screenshots.py --all --check --strict
python tools/reproduce_examples.py --group all --list
python tools/reproduce_examples.py --group all \
  --output-root <empty-reproduction-root> \
  --missing-report <missing-input-report.json>
```

T-PY-05/T-PY-07のfull LOSATP reproductionは通常PR testへ戻さず、release evidenceとして
明示実行する。実装開始時に各CLI optionの存在を`--help`で再確認し、stale plan commandを
無条件に実行しない。`--list`はinventory/readiness確認にすぎず、render evidenceに数えない。
Rendered outputsはfigure manifestのsemantic/hash oracleとreadable-scale visual reviewへ渡す。

#### Completion gate

- candidate docsがunreleased routeやmissing DOIをliveと主張しない。
- exact built artifactでmanifest-declared recipesが成功する。
- public/release/manuscript figureは再現可能、またはnon-reproducibleと明示されclaimから外れる。
- changed figureにはartifact path、command、review evidenceがある。

### Phase 8: J-RC candidate gateを実行する

Status: pending

#### Entry gate

- Gate 0 manifest frozen;
- B–I evidence、A1 Phase 0–5 candidate-synchronization milestone、H-final
  exact-artifact evidence complete;
- exact artifact hashes fixed;
- no unowned mandatory requirement;
- no unresolved P0/P1。

#### Full release gates

```bash
python -m pytest tests/ \
  -m "not slow and not (recipe or gallery or browser)" \
  -n auto --dist loadfile --durations=30
python -m pytest tests/ -m "recipe and not slow" --durations=30
python -m pytest tests/ -m "gallery and not slow" --durations=30
python -m pytest tests/ -m "browser and not slow" --durations=30
python -m pytest tests/ -m "slow" --timeout=600 --durations=30
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
node --test tests/web/*.test.mjs
node tests/web/session-request.test.mjs \
  --project-session gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz
node tests/web/session-request.test.mjs \
  --project-session gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json
node tests/web/session-request.test.mjs \
  --project-session gbdraw/web/gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz
npx playwright test --project=chromium --workers=1
python tests/run_losat_cache_browser_acceptance.py
ruff check gbdraw/
git diff --exit-code -- tests/reference_outputs/
```

OS/Python/browser matrix、exact artifact smoke、offline、docs/reproductionはhosted workflowの
job URLとhashを同じledgerへ結ぶ。

上のpytest commandはcandidate sourceのregression gateである。
`tests/test_web_packaging.py`などがtest用distributionをbuildしても、そのbytesをH-final artifact
evidenceへ流用しない。Exact-artifact certificationはPhase 5–7のuploaded path/hashをconsumeする
別jobだけが所有し、source regressionと二重計上しない。

#### Decision

- P0 = 0、P1 = 0。
- P2はissue、impact、workaround、owner approval、release-note decisionを持つ。
- unexpected skip/xfail/missing job/retry-only pass = 0。
- evidence manifestはcandidate SHA/artifact hashと一致する。

J-RC pass後はKへRC tag/TestPyPI/staged smokeのauthorisation requestだけを渡す。

### Phase 9: RC fix invalidationとJ-Finalを実行する

Status: pending

#### 作業

1. RC/TestPyPI/staged smokeで見つかったdefectをP0/P1/P2へ分類する。
2. accepted fixごとにregression testとinvalidated evidenceを記録する。
3. new RCが必要ならPhase 8をnew candidate/hashで再実行する。
4. A1 D13のsingle-owner transactionを要求する。H ownerがfinal version/development
   classifierを、A1がrelease notes/citation/public proseを、declared generatorが
   recipe/Web wheel/noticesなどのderived fileを更新する。Jはそのtransactionを記録・認定し、
   別のfinal-doc rewriteを所有しない。
5. final-version commitからartifactを一度buildし、Phase 2–8のrequired gatesを新hashで再実行する。
6. final release manifestとrollback/yank/patch procedureをKへ渡す。
7. Section 12.2のK action ledgerへexpected commit/tag/hashとpost-action smokeを記入する。

#### Completion gate

- final commit、intended tag、wheel、sdist、Web bundle、session/request versions、docs recipes、
  figure evidenceが一つのmanifestにbindされる。
- RC artifact hashをfinal artifact evidenceとして流用しない。
- Jはtag/publish/deployせず、explicit K action待ちで停止する。
- K action ledgerはexpected欄がcompleteで、observed external-action欄はK所有の`pending`である。

## 8. CI and support matrix

### 8.1 Required release jobs

| Job | Matrix | Source or artifact | Scope |
|---|---|---|---|
| Core | Linux × Python 3.10/3.11/3.12 | source + locked candidate deps | complete non-recipe/Gallery/browser core |
| Recipes | canonical Linux Python plus declared supported-version acceptance | exact wheel | literal CLI/Python recipes |
| Gallery/docs | canonical Linux Python | exact wheel/Web bundle | sessions, regenerated assets, captures |
| Browser unit | canonical Linux/Node | frozen candidate source modules | all `*.test.mjs`; do not mislabel source imports as artifact evidence |
| Playwright full | Chromium | exact extracted Web bundle through explicit artifact root | all `*.playwright.spec.js` |
| Browser smoke | each additional supported engine | exact Web bundle | ready, worker/Wasm, gzip session, export, privacy |
| Clean install base | Linux/macOS/Windows × frozen Python matrix | exact wheel and exact sdist; record the sdist-derived temporary wheel separately | import, pip check, CLI, SVG, session, GUI |
| Clean install export | representative supported OS/Python cells | exact wheel | PNG/PDF and optional dependency boundary |
| Offline/security | Chromium desktop+narrow | installed/extracted artifact | blocked network, sanitizer, canary, repeated/cancelled run |
| Slow regression | canonical release environment | frozen candidate source | slow tests; any test-built distribution is disposable regression output |
| Reproduction | canonical release environment | exact installed wheel plus clean archive/bundle inputs | LOSATP/manual scenarios and public figures |

### 8.2 Matrix reduction rule

runtimeが過大な場合、full logic suiteをLinux全Python、OS matrixをmin/max Python smokeへ縮小
できる。ただしGate 0でsupport manifestと理由をfreezeし、metadata/docsを一致させる。実行後に
red jobを消すためのmatrix reductionは認めない。

## 9. Severity and invalidation

| Severity | Release effect | After publication |
|---|---|---|
| P0 | candidate testingとrelease actionを停止 | immediate rollback/yank assessment |
| P1 | J-RC/J-Final promotionをblock | patch-release decision |
| P2 | maintainer-approved waiverだけ許可 | tracked follow-up |

P0例: data loss、incorrect biological coordinate、corrupted session、forbidden network content、
pre-consent/reject/revoke analytics request、install/import failure、materially incorrect export。

P1例: core crash、render-only analysis access、concurrent/stale commit、partial Result transaction、
undo mismatch、last-good loss、supported platform/browser core failure、required control inaccessible、
documented command failure。

次の変更は該当evidenceを必ずinvalidateする。

- productionまたはdependency change: focused + full + artifact hash;
- schema/compatibility change: all migration/public/session/docs gates;
- Web asset/CSP/analytics change: Web bundle、offline、privacy、deploy hash;
- docs/recipe/figure change: A1、reproduction、sdist/Web bundle where shipped;
- version/metadata change: build、concordance、clean install、hash、citation;
- workflow-only change: hosted dry runとpermission/artifact-flow review。

## 10. Stop conditions

次の場合はJをcompleteにせず、owner packageまたはmaintainer判断へ戻す。

- upstream mandatory criterionにtest/evidence ownerがない;
- partial/unreleased schemaをsupportする必要が生じる;
- universal wheelがsingle-platform native binaryを含んだままになる;
- exact artifactがmatrix jobで再buildされる;
- supported platform/browser jobがskipまたはretryだけでpassする;
- analytics/offline/CSP contractが両立しない;
- packaged Webがexternal networkなしにreadyにならない;
- public/manuscript figureがmissing/unlicensed inputへ依存したままreproducibleと主張される;
- candidate fix後にartifact hashまたはaffected evidenceが更新されない;
- canonical hosted bundle/deployment pathがKのpost-smoke用にaccepted release commit/hashを
  識別できない;
- readiness workflowがtag、publish、deploy、archive、submissionを実行しなければpassできない設計になっている。

## 11. Risk register

| Risk | Impact | Mitigation |
|---|---|---|
| Playwright suite runtime/flakiness | release gateが不安定 | one-worker race shard、trace、deterministic fixture、retry-only pass拒否 |
| OS/Python matrix cost | feedback遅延 | logic suiteとartifact smokeを分離し、Gate 0でmatrixをfreeze |
| Large Web wheel/bundle | upload/install timeout | size inventory/budget、build once、cacheではなくartifact reuse |
| Native LOSAT removal | Linux convenience regression | external discovery smoke、clear diagnostics、Web Wasm parity、docs |
| Historical fixture不足 | compatibility claimが弱い | first-parent/tag fixture、hash/provenance、namespace別support |
| GA/CSP conflict | privacy/offline P0 | retained/cut binary decision、intercepted network tests、local absence |
| Version duplication | wrong wheel/Web/session/docs | concordance gate、derived expected names、generated file手編集禁止 |
| A1/H-final cycle | stale docsまたはartifact | provisional sync → exact build → J exact-artifact rerun |
| Missing figure inputs | paper reproducibility failure | obtain licensed hash-pinned inputs or remove/label claim |
| main auto-deploy drift | tagged releaseとhosted app不一致 | accepted bundle deployment、build stamp/hash post-smoke |

## 12. Evidence ledger

### 12.1 J phase ledger

| Phase | Status | Required evidence |
|---|---|---|
| 0 Baseline/scope template | pending | HEAD/tree, dependency ledger, inventories, scope/support undecided list |
| 1 Traceability | pending | mandatory requirement map with no unowned row |
| 2 Compatibility/public contract | pending | writer/readers, fixture provenance/hashes, focused Python/Node results |
| 3 CI/Playwright | pending | hosted job URLs, discovery count, reports/traces, skip/retry summary |
| 4 Package truth/version | pending | universal-wheel decision, concordance tests, content audit |
| 5 Build once/matrix | pending | exact artifact names/sizes/hashes, install matrix, permission dry run |
| 6 Offline/security/privacy | pending | blocked-network log, canary result, sanitizer/lifecycle evidence |
| 7 Docs/reproducibility | pending | recipe/capture/Gallery results, figure manifest, visual review |
| 8 J-RC | pending | complete gate table, P0/P1=0, approved P2, candidate manifest |
| 9 J-Final | pending | final-version hashes, full rerun, K handoff, no external action by J |

Phase statusは実行結果を観測した時だけ更新する。plan作成、test collection、near-pass、budget不足、
外部action待ちを`complete`として扱わない。

### 12.2 K action ledger template

この表のexpected欄はJ-Finalが所有し、observed欄とaction authorityはWork package Kが所有する。
`pending`はJの失敗ではないが、K-Publicationはpassしていない。

| External action | Expected baseline from J | Action authority | Observed evidence and post-smoke | Status |
|---|---|---|---|---|
| RC tag | intended RC tag and accepted commit | pending maintainer approval | pending | pending |
| TestPyPI upload/exact-download smoke | accepted RC wheel/sdist hashes and H-AC4 oracle | pending maintainer approval | pending | pending |
| Staged hosted deploy/smoke | accepted RC Web bundle hash/build stamp and privacy profile | pending maintainer approval | pending | pending |
| Final tag | intended final tag and accepted commit | pending maintainer approval | pending | pending |
| PyPI publication/post-smoke | accepted final wheel/sdist hashes and clean-install oracle | pending maintainer approval | pending | pending |
| Hosted production deploy | accepted Web bundle hash/build stamp | pending maintainer approval | pending | pending |
| Bioconda | final tag and source/artifact hash | pending maintainer approval | pending | pending |
| Zenodo/version DOI | final tag, archive inventory/hash | pending maintainer approval | pending | pending |
| A1 live-identifier closeout | observed package/deploy/archive identifiers | pending maintainer approval | pending | pending |
| Preprint revision | K-Publication bundle and approved claims | pending maintainer approval | pending | pending |
| Journal submission | final manuscript/cover-letter bundle | pending maintainer approval | pending | pending |

## 13. Definition of done

Work package Jの実装は次をすべて満たしたときcompleteである。

- Gate 0 scope/schema/support/compatibility manifestがfrozenである。
- mandatory B–I/A1/H criterionにnamed evidence ownerがある。
- current writerとsupported migrationsがauthentic fixtureでPython/Web/CLI replayされる。
- package-root、typed API、CLI、version metadataがapproved contractと一致する。
- complete core、rendering、Node、Playwright、offline/security/privacy、slow/reproduction gateがpassする。
- wheel/sdist/Web bundleが一度だけbuildされ、同じhashでcontent auditとsupported matrix installをpassする。
- base installと`[export]`が分離され、universal wheelのplatform claimがtruthfulである。
- public/release/manuscript candidate figuresがreproducibleまたは正直にscope外と記録される。
- no unresolved P0/P1がなく、全P2にapproved dispositionがある。
- J-RCとJ-Final ledgerがexact candidate/final commitとartifact hashへ結び付く。
- K action ledgerのexpected baselineとpost-action smoke requirementsが記入されている。
- Jはexternal actionを実行せず、Kへexplicit-authorisation handoffを返す。
