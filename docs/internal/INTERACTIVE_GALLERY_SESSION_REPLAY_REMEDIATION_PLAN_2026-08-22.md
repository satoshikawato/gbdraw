# Interactive Gallery session replay remediation plan

- 作成日: 2026-08-22
- 改訂日: 2026-08-22
- 文書区分: 内部実装計画
- 状態: Ready for implementation — audit revisions incorporated; implementation evidence is still pending
- 対象: Web session import、次回 Generate、Gallery session publication、公開 Gallery 11 session、refresh-owned test fixture 2 session、Gallery 生成成果物、browser regression、tutorial media
- 非対象: 一般ユーザー session の未生成 draft 契約変更、session version 41、canonical request schema 6、renderer geometry の仕様変更

## 1. 結論

Interactive Gallery の不具合は、単一の JSON 欠損ではない。次の三つを、それぞれの既存 owner に沿って修正する。

1. Gallery の保存済み preview と、次回 Generate に使う active Web state が分離している。
2. WSSV の遅延 session resource を、Generate 側が File ではないと誤判定している。
3. tobacco chloroplast の current v40 session に、廃止済み Web field が混入している。

Gallery には次の公開不変条件を導入する。

> 公開 Gallery session は、current reader が例外なく受理する clean committed session である。import 直後の active state から行う最初の Generate は、保存済み committed render intent と意味的に同じ request と図を生成しなければならない。

一般ユーザーの current session は、最後に生成した `renderRequest` と未生成の編集 draft を同時に保存できる。この契約は変えない。通常 importer の `config` を常に committed request で上書きせず、整合性の強制は Gallery publication pipeline に限定する。

監査後の実装単位は3 PRとする。PR 1とPR 2は並行開発できるが、PR 2の最終 head はPR 1を取り込み、WSSVを含む全公開例のgateを通したdeploy可能な状態でのみmergeする。

| PR | 目的 | 指示書 |
| --- | --- | --- |
| 1 | WSSV lazy resource の読込境界を修正し、専用browser ownerを最初から最終配置に置く | [PR 1](./INTERACTIVE_GALLERY_SESSION_REPLAY_PR_01_WSSV_LAZY_RESOURCE_READ.md) |
| 2 | strict current-v40契約、Gallery publication、全13 refresh、公開asset、replay/CI/deploy gateを一つのdeployable vertical sliceとして実装する | [PR 2](./INTERACTIVE_GALLERY_SESSION_REPLAY_PR_02_PUBLICATION_VERTICAL_SLICE.md) |
| 3 | tutorial、capture metadata、必要なoperation mediaを修復済みsessionに合わせる | [PR 3](./INTERACTIVE_GALLERY_SESSION_REPLAY_PR_03_TUTORIAL_MEDIA_AUDIT.md) |

この統合は単なるPR数削減ではない。旧PR 2のschema変更、旧PR 3のpublication path、旧PR 4のartifact/gateは、どれか一つだけをmergeするとinvalid artifactまたは未検査artifactをdeployし得るため、同じdeployable unitにする。production、test、generated artifact、workflowのdiffはPR 2内で別commitまたは別review sectionとして分離する。

## 2. 観測済みの基準状態

### 2.1 ユーザー可視の再現結果

| Example | import | 最初の Generate | 観測された差 |
| --- | --- | --- | --- |
| `HmmtDNA_ATskew` | 成功 | 成功 | SVG 130,397文字から129,177文字へ変化。AT skew slotは残るが、slot表現とcolor/priority resource bindingがround tripしない |
| `majanivirus_orthogroup` | 成功 | 成功 | bitscore 100から50、identity 20から70、arrow shaft ratio 0.5から1.0。pairwise matches 627から262。`WSSV-like proteins` captionも消える |
| `vibrio-harveyi-group-collinear` | 成功 | 成功 | 保存済みpreviewにはcomparisonがあるが、active `linearComparisonPlan` は `none`。Generate後はcomparison group、match、legendが0 |
| `WSSV_genome_comparison` | 成功 | 失敗 | 20 FASTAは復元済みだが、`Input file is not available for browser FASTA extraction.` で停止 |
| `tobacco-chloroplast` | 失敗 | 未到達 | `config.adv.cli_circular_track_order` と `config.adv.cli_circular_track_slots` がcurrent v40 allowlistにない |

文字数や現在のmatch数は診断時のevidenceであり、raw byte固定値を将来の唯一の契約にはしない。完成図の科学的・視覚的意味、committed request、初期preview、first Generateの関係を契約にする。

### 2.2 現行テストの盲点

診断時には次がすべて成功した。

```bash
python -m pytest \
  tests/test_gallery_session_semantics.py \
  tests/test_wssv_gallery_fastas.py \
  -q

node tests/web/wssv-gallery-fastas.test.mjs
node tests/web/gallery-session-migration.test.mjs
node tests/web/session-draft-authority.test.mjs
```

Python 25件と静的JavaScript検査が通っても、実browserの不具合は再現した。現在の検査は主にcommitted `renderRequest`、保存済み`results`、source/example SVG間を比較し、import後のactive stateから次回requestを構築していない。

既存のVibrio full-generation specは、active plan `none`、comparison count 0、生成後comparison 0を期待しており、不具合を正常動作として固定している。また専用npm scriptがCIから呼ばれていない。この期待値とCI wiringを反転修正する。

## 3. 監査後の設計判断

### 3.1 Session version契約はcurrentとlegacyを混ぜない

| 入力 | 許可する処理 | 結果 |
| --- | --- | --- |
| canonical `renderRequest`を持つsupported v31–33/39 | version-specific Web/Gallery migration後にcurrent validatorを通す | valid v40、または明示的失敗 |
| supported v27–30 | 既存のCLI replay inputとしてのみ受理する | 通常CLI互換性を維持し、synthetic Web/Gallery draftを作らない |
| development-only v34–38 | 拒否する | outputなし |
| valid v40 | migration前にcurrent validatorを通し、identity pathで扱う | 内容をsilent cleanupしない |
| invalid v40 | migrationやpublication preparationの前に拒否する | path付きvalidation error |

ここでidentityとはadmission/version migration段階の契約である。valid v40をcurrent-version migrationでcleanupしない。admission後に明示的なGallery generatorが行うpublication rewriteは別段階であり、ordinary importから到達不能、committed request semanticsを保存、canonical equivalence成功時だけwrite、という条件でactive draftをwriter-native stateへ揃えられる。

したがって、current v40のpromoterに廃止fieldのsilent cleanupを追加しない。`cli_circular_track_order` と `cli_circular_track_slots` を含むtobacco chloroplast artifactは、PR 2の作業中に一度だけproducer-side data migrationを行う。その変換は次を満たす。

- production importer、promoter、refresh ownerに恒久的な例外分岐を残さない。
- 一時変換コードは最終headから除去する。
- 変換前後のsemantic diff、実行command、current validator通過結果をPR evidenceに残す。
- invalid tracked sourceを先に置換せず、staged overrideとして全13 refresh/asset生成と同じ外側transactionへ渡す。
- staged bootstrap後のfailure injectionで、final replacement前ならtracked sourceが元のままであることを証明する。
- 変換後は通常のpublication/refresh owner commandだけで再生成できる状態にする。

PythonとWebのcurrent validatorは、廃止fieldを同じく拒否する。ただしPythonへWeb全state schemaを複製せず、この既知producer driftに対するcontract testだけを共有fixtureで固定する。Gallery publicationのhistorical admissionはv31–33/39に限定し、v27–30の通常CLI replay compatibilityを本計画の都合で変更しない。

### 3.2 Publication semanticsには専用ownerを置く

707行規模の`gallery-session-migration.js`をpublication coordinatorとしてさらに肥大化させない。責務を次のように分ける。

| 責務 | 最終owner |
| --- | --- |
| supported legacyからcurrentへのversion migration | `gbdraw/web/js/services/gallery-session-migration.js` |
| current writer active-configのfield/shape admission contract | `config.js`から抽出するDOM-free shared contract（仮称`session-active-config-contract.js`） |
| current requestのprojection/buildとsemantic canonical form | `gbdraw/web/js/services/session-request.js`と既存codec owner |
| Gallery active state、metadata policy、resource binding、prepare/finalize validation | 新規の小さな`gbdraw/web/js/services/gallery-session-publication.js` |
| bulky result/cacheのstaging、subprocess、file I/O、restore | `tools/refresh_gallery_sessions.py` |
| session content capability判定 | `readFileText()`を持つ既存file-content owner |

shared active-config contractはreactive `state`やDOMをimportせず、current writerの許可field/shapeを一度だけ定義する。通常import/restoreとGallery publication admissionが同じvalidatorを呼ぶ。`gallery-session-migration.js`へcurrent validatorを戻して循環dependencyを作らず、migration moduleはhistorical-onlyを保つ。

Python helperは、request/resource precedence、active state、metadata retentionの意味を決めない。Pythonは大きな`results`やcacheをstagingし、JavaScript ownerへ必要なsemantic projectionとresource descriptorを渡し、戻り値を機械的に組み立てる。

特にrefreshの既存silent/semantic pathsは名前付きで処理する。

- `_load_gallery_refresh_source()`を削除し、通常`load_session()`のfail-closed admissionへ統一する。invalid currentのResult/catalogを消して再validationしない。
- refresh前の`normalize_current_session_artifacts(session)`呼出しは削除する。core writerで必要な同名normalizerまで無関係に削除せず、Gallery source repairとして使わない。
- `_omit_regenerable_gallery_derived_cache()`のVibrio固有policyは削除するか、一般化されたJavaScript finalizer policyへ移す。
- `_preserve_gallery_cli_invocation()`のraw argv rewriteは削除し、opaque provenanceの保持をJavaScript finalizerの一つのpolicyにする。

invalid catalog、cache、resource bindingを欠落情報の削除で通すfallbackを残さず、negative testでfail closedを固定する。

current v40とGallery publicationでは、`cliInvocation.args`をsemantic authorityにしない。active planは、まずcommitted topology、canonical request、analysis metadataから決定する。それでも表現不能なGallery固有情報が実在する場合だけ、versioned structured Gallery recipeを最小追加する。raw argvはprovenanceと整合性確認に留め、argv parserを第二のcurrent-request ownerにしない。

ただし、supported pre-v40のversion-specific migrationで、当時のpersisted contractにargv由来情報が含まれていたことをpositive fixtureとrelease/commit evidenceで証明できる既存pathは維持できる。その解釈をcurrent/publicationへ流用・拡張せず、legacy migration owner内に閉じ込める。今回のremediationを理由に既存の実証済み互換性まで削除しない。

### 3.3 Request equivalenceは一つのcanonical formを再利用する

publication専用のfield-by-field comparatorを別実装しない。`session-request.js`またはそのcodec ownerに、render requestを比較用canonical formへ写す関数を置き、次の両方が利用する。

- publication preparation後のcommitted request対rebuilt next request
- staged refresh後のcommitted request対rebuilt next request

canonical formは少なくとも次を保持する。

- mode、grouping、record数・順序・source kind・selector・region・presentation
- complete `diagramOptions`のrender-affecting fields
- track slot ID、renderer、enabled、side、radius、width、gap、z、params、axis index
- annotation set、feature visibility、label/filter/priority rule、color、legend、title
- layout、record position、record gap、comparison height
- comparison kind、endpoint、program/mode/filter、pair、canonical analysis metadata
- resource kindとdecoded payload SHA-256
- saved preview/editor semanticsに影響するinteractive metadata policy

resource IDは`{kind, payloadSha256}`へ置換してからstable serializeする。除外できるのは、明示された非semantic fieldだけである。初期allowlistは`createdAt`、runtime telemetry、result index、output prefix、overwrite、artifact file name、publication toolが強制するoutput-format representationに限定する。除外fieldを増やす場合は、SVG、feature catalog、editor、cache identityのいずれにも影響しないevidenceとnegative testを同じ変更に置く。

SVGについても別言語のfingerprint policyを増やさない。`tests/utils/svg_compare.py`の既存canonicalization/comparisonを唯一のsemantic comparatorとして拡張し、必要ならdigestとstructured diffを返す薄いAPIを加える。Playwright側は生成物を渡すだけのwrapperにし、独立した`gallery-svg-fingerprint.cjs`は作らない。既存`session-regenerate-intent.playwright.spec.js::svgEquivalence()`のDOM semantic normalization/OR判定もこのshared comparatorへ収束させる。raster/direct-edit evidenceを残す場合はsemantic equivalence ownerではない別目的を明記する。curated per-example anchorは「両方が同じように壊れる」場合を検出する別目的なので維持する。

### 3.4 Inventory ownerを一つにする

公開11例のauthoritative inventoryは`tools/prepare_interactive_gallery_assets.py::EXAMPLES`とする。

- `GALLERY_SESSION_FILES`の独立列挙は削除するか、`EXAMPLES`から機械的に導出する。
- `examples.json`はgenerated projectionであり、generator入力や第二のauthorityにしない。
- JS/browser testが`examples.json`を列挙する場合は、`EXAMPLES`から生成したmanifestとの集合・順序一致をread-only guardで検証する。
- `TEST_INPUT_SESSION_FILES`の2件は公開inventoryへ混ぜず、unfiltered refreshが所有するfixture inventoryとして明示的に維持する。

これによりpublic 11、refresh-owned 13、browser表示inventoryの関係を、重複手書きlistではなく1 owner + 1 separate fixture setで表す。

### 3.5 Transaction保証は実装より強く書かない

既存`_gallery_file_transaction()`は、全validation後にstaged fileを逐次置換し、捕捉した例外ではbackupからrestoreする。process kill、host crash、runner停止を含むcrash-atomic transactionではない。

本計画が要求する保証は次のとおりである。

> 全13 inputと全generator-owned outputをstagingで検証してから置換を開始し、置換中に捕捉した例外では変更済みfileをrestoreする。

このcaught-exception transactional replacementをfailure-injection testで証明する。crash atomicityが本当に必要になった場合は、directory-level swapやjournalを別設計として扱い、本PRへ暗黙に追加しない。

### 3.6 最大sessionの性能を観測値だけにしない

Vibrio等の最大sessionについて、unfiltered refreshとpublication projectionを同じrunner/imageでbaselineとheadを各3回測る。

baseline/headは記録したSHAの別clean worktreeで実行し、全出力をdisposable staging pathへ向けて破棄する。測定runでreview worktreeのGallery targetを置換しない。review対象artifactを作る一回のowner commandとは分離する。

- wall timeは3回のmedianを比較する。
- peak RSSは3回のmaximumを比較する。
- headはwall timeとpeak RSSの双方でbaselineの1.25倍以下をmerge gateとする。
- 閾値を超える必要がある場合は、原因、実測値、代案をarchitecture review packetで明示し、単に記録してpassさせない。

実装は`results`やlarge non-projection stateをclone/hash前に分離する。resource digestはchunk処理し、JSON全体の重複cloneや全buffer化を新たに増やさない。

### 3.7 Hosted buildはread-onlyに統一する

GitHub PagesとCloudflare Pagesの双方で、hosted buildはGallery成果物を変更しない。

- Galleryを変更できるowner commandは明示的な`python tools/refresh_gallery_sessions.py`だけにする。
- `deploy_web.yml`からdeploy時refreshを除去する。
- `tools/prepare_cloudflare_pages.py --refresh-gallery`は削除または非mutating verificationへ再分類し、build時の別producer pathを残さない。
- 共通のread-only verifierがauthoritative inventory、physical file、manifest/hash、bundle byteを照合する。
- hosted buildの前後でGallery targetにworktree diffがあれば失敗する。

## 4. 成功条件と不変条件

### 4.1 機能上の不変条件

- 全11 Gallery sessionがcurrent Web importerで成功する。
- importだけではdiagram Workerを起動せず、保存済みpreviewを表示できる。
- import直後のfirst Generateが`status: ok`で完了する。
- committed request、active-stateからbuiltしたrequest、初期preview、first Generateがrender-affecting semanticsで一致する。
- cache/lifecycleが重要な例ではsecond Generateも同じintentと図を生成する。
- WSSVは外部FASTAを要求せず、埋め込み20 FASTAと保存済みcacheだけで生成できる。
- Vibrioは11 records、record order、5 rows、48 px gap、collinear comparisons、legendを保持する。
- Majanivirusは9 records、committed filter値、627 matches、curated captions/colorsを保持する。
- AT skewは`features,gc_content,gc_skew,a_skew_2,ticks`の順序、AT skew colors/caption、labels、paletteを保持する。
- chloroplastは`plastome_regions`の4 annotation、3 slots、upper-left legendを保持する。
- 通常ユーザーsessionのcommitted resultと未生成draftの分離を維持する。

### 4.2 Offline、privacy、failure不変条件

- runtime network dependencyを追加しない。
- Gallery sessionは入力、比較資源、再生成に必要なprepared dataを自己完結して保持する。
- WSSV、Majanivirus、VibrioのGenerateがcache missを外部実行で隠さない。
- sequence本文やfull comparison rowsをlogへ出さない。
- resource failureは安全なstage/resource labelを含む構造化errorとして返す。

### 4.3 成果物不変条件

- `gbdraw/web/gallery/sessions/`、`sources/`、`examples/`、`thumbnails/`、`examples.json`は一つのgenerator pathから更新する。
- 無引数refreshが所有する2 test fixtureも、公開11件と同じstaging/validation transactionで扱う。
- 生成物を手編集しない。
- `tests/reference_outputs/`はrenderer geometryの意図的変更がない限り変更しない。
- `examples/gbdraw_social_preview.png`は変更しない。
- tutorial mediaは実UIと修復済みexact sessionからのみcaptureする。

## 5. Target architecture

### 5.1 Gallery publication flow

```text
authoritative EXAMPLES entry + source session
  -> shared DOM-free admission contract
     -> validate current v40 without cleanup
     -> or migrate canonical legacy v31–33/39 then validate v40
     -> keep v27–30 outside Gallery publication as CLI replay-only compatibility
  -> JavaScript prepareGallerySessionForPublication
     -> project committed request through session-request owner
     -> rebuild writer-native active config/files/comparison/resource bindings
     -> build next request through the same session-request owner
     -> compare the two canonical request forms
  -> Python stages bulky result/cache and invokes canonical CLI replay
  -> JavaScript finalizeGallerySessionPublication
     -> bind refreshed resources/results without changing semantic policy
     -> rebuild and compare next request again
  -> validate every staged input/output
     -> public 11 sessions
     -> refresh-owned 2 test fixtures
     -> source/example/thumbnail/examples.json manifest and hashes
  -> caught-exception transactional replacement
  -> read-only deploy verifier
```

`prepareGallerySessionForPublication()`と`finalizeGallerySessionPublication()`は専用moduleの小さなpublic surfaceとする。legacy migration、request codec、resource content readerを再実装せず、それぞれのownerを呼ぶ。

### 5.2 Complete owner/path declaration

#### Capability: current Web active-config admission

Owners before:

```text
{
  config.js::validateCurrentWriterActiveConfig coupled to reactive state,
  refresh_gallery_sessions.py::_load_gallery_refresh_source cleanup/retry
}
```

Owners after:

```text
{
  DOM-free session-active-config contract consumed by ordinary restore and Gallery publication
}
```

`config.js` remains the restore/save orchestrator but consumes the extracted contract. `gallery-session-migration.js` does not regain current admission. Python keeps its narrow envelope/current-retired-field checks and does not clone the full Web allowlist. The cleanup/retry path is removed rather than classified as compatibility.

#### Capability: Gallery active render-state publication

Complete scope: Gallery refresh時のactive Web config、canonical request、resource binding、metadata retention、next Generate一致を決める全production path。

Owners before:

```text
{
  gbdraw/web/js/services/gallery-session-migration.js legacy promotion,
  tools/refresh_gallery_sessions.py::_sync_legacy_legend_control_with_render_request,
  tools/refresh_gallery_sessions.py::_sync_circular_track_draft_with_render_request,
  tools/refresh_gallery_sessions.py::_restore_rendered_palette_file_binding,
  tools/refresh_gallery_sessions.py metadata/resource/request merge helpers
}
```

Owners after:

```text
{
  gallery-session-migration.js: supported version migration only,
  session-active-config contract: DOM-free current admission,
  session-request.js + codec: request projection/build/equivalence,
  gallery-session-publication.js: Gallery prepare/finalize policy,
  refresh_gallery_sessions.py: process, staging, I/O, restore only
}
```

Canonical paths before:

```text
{
  current v40 -> Python ad-hoc sync -> CLI replay -> stale config merge,
  legacy -> Node promotion -> CLI replay -> Python semantic merge
}
```

Canonical paths after:

```text
{
  every refresh-owned session -> validate/migrate -> JS prepare
  -> CLI replay -> JS finalize -> staged validation -> replace
}
```

Compatibility paths after: supported pre-v40 migrationだけを維持する。invalid current-v40 artifact用のreader fallback、silent cleanup、恒久repair modeは0件。

#### Capability: browser file-content availability

Owners before:

```text
{
  services/file-content-cache.js::readFileText,
  app/run-analysis.js FASTA extraction .text capability prechecks
}
```

Owners after:

```text
{
  services/file-content-cache.js::readFileText
}
```

FASTA extractionはconsumerとして残るが、File-like resourceのreadabilityを独自判定しない。

#### Capability: public Gallery inventory and hosted bytes

Owners before:

```text
{
  prepare_interactive_gallery_assets.py::EXAMPLES,
  refresh_gallery_sessions.py::GALLERY_SESSION_FILES,
  gallery/examples.json,
  deploy-time refresh paths
}
```

Owners after:

```text
{
  prepare_interactive_gallery_assets.py::EXAMPLES as source inventory,
  refresh_gallery_sessions.py::TEST_INPUT_SESSION_FILES as separate fixture inventory,
  refresh_gallery_sessions.py as the only mutating owner,
  shared verifier as read-only consumer
}
```

PR 2のarchitecture packetでは、実装後の正確なsymbol名を使ってbefore/after set、削除したhelper/path、追加moduleの責務、compatibility path数、performance evidenceを更新する。agentはarchitecture approvalを自分で投稿しない。

## 6. PR dependency and merge strategy

| Order | Branch suggestion | Dependency | Merge condition |
| --- | --- | --- | --- |
| PR 1 | `fix/wssv-lazy-session-resource` | latest `origin/dev` | shared-reader境界、focused test、専用WSSV real Generate成功 |
| PR 2 | `fix/gallery-publication-vertical-slice` | 開発はPR 1と並行可。最終head/mergeはPR 1必須 | strict v40、publication、全13 refresh、全公開asset、browser/CI/deploy gate、performance gateが一体で成功 |
| PR 3 | `docs/gallery-replay-tutorial-media` | PR 2 merge後 | strict capture audit、必要mediaだけ再撮影、desktop/mobile Gallery確認 |

各branchは最新`origin/dev`をfetchした後、repository policyに従い`git switch --no-track -c <name> origin/dev`で作る。PR 1のWSSV browser specをPR 2で移動・複製しない。PR 2はその最終ownerをCIから呼ぶ。

PR 2内では少なくとも次のreview unitを分ける。

1. production contract/owner convergence
2. focused/static/browser testsとperformance harness
3. one-time chloroplast data migration evidence
4. 全13 sessionと公開Gallery generated artifact
5. `gbdraw/web/CLAUDE.md` owner tableと、必要最小限の`tools/web-change-policy.json` architecture registry
6. CI、GitHub Pages、Cloudflare read-only deploy wiring

どの中間commitも単独merge対象にはしない。最終diffがdeploy可能であることをPR 2のacceptanceとする。

## 7. Test and CI strategy

### 7.1 Layer A: owner-boundary tests

- `readFileText(SessionResourceFileView)`が`.text()`なしで動作する。
- FASTA extractionがshared readerの成功・失敗semanticsをそのまま使う。
- canonical legacy v31–33/39 fixtureはGallery publication用currentへmigrationでき、retired fieldを出力しない。
- v27–30は通常CLI replay互換性を維持するが、synthetic Gallery/Web publication draftへ入らない。
- valid current v40はidentity semanticsを保つ。
- invalid current v40はJavaScript/Pythonの双方でmigration前に拒否される。
- ordinary restoreとGallery publicationが同じDOM-free active-config contractを使い、migration moduleとのcycleを作らない。
- normal current-session draft authority testはdivergent draftを引き続き保持する。
- publication coordinatorがPythonのsemantic mergeなしでprepare/finalizeできる。
- equivalence canonical formはresource ID差を許すが、kind/payload、track、filter、comparison、caption差を拒否する。

### 7.2 Layer B: DOM-free all11 publication gate

authoritative `EXAMPLES`から公開11件を列挙し、一件ごとに別processで検査する。大きなgzip/session dataを次ケースへ保持しない。

各sessionで次を確認する。

1. current Web config validation成功。
2. resource table adoption成功。
3. active config/files restore成功。
4. first Generate canonical request構築成功。
5. committed requestとのsemantic equivalence成功。
6. 全resource bindingが存在し、payloadが非空で、hashが一致。
7. generated `examples.json` inventoryがauthoritative inventoryと一致。

このgateはGallery publication専用であり、一般sessionへ適用しない。

### 7.3 Layer C: actual browser replay

全11を実Web appで検査する。

```text
fresh browser context
  -> import exact published session
  -> assert saved preview visible and Worker idle
  -> capture initial SVG for the shared Python comparator
  -> invoke Generate
  -> assert status ok and no unsafe fallback/network
  -> capture generated SVG
  -> compare initial and generated semantics
  -> run second Generate where cache/lifecycle behavior matters
```

- WSSVはPR 1で置く専用serial spec/config/npm scriptを唯一のfull-render ownerとし、移動しない。
- Vibrioは既存`vibrio-full-generation.serial.spec.js`を唯一のfull-render ownerとして期待値を反転する。
- 残り9例はGallery replay specでcircular/linearを分離する。
- WSSVとVibrioは別process、`workers: 1`で実行し、cache missやmemory failureを明確に帰属させる。
- failure時だけinitial/generated SVG、structured diff、trace、safe diagnosticsをartifactとして保存する。

### 7.4 Curated per-example anchors

dynamic parityだけでは両方が同じように壊れる可能性がある。既存`test_gallery_session_semantics.py`の完成図contractを維持する。

| Example | Curated anchors |
| --- | --- |
| HmmtDNA basic | 1 record、GC content/skew、ticks、right legend |
| Lambda | 1 record、comparisonなし、ruler、left legend |
| AT skew | 5 slots、AT skew positive/negative、ajisai、gene labels |
| tobacco chloroplast | 4 region annotations、3 slots、upper-left legend |
| Vnig | 6 records、multi-record positions、left legend、bottom title |
| Hepatoplasmataceae collinear | 5 records、collinear blocks、record-local tracks |
| Vibrio | 11 records、5 rows、48 px gap、nonzero committed-equivalent collinear blocks |
| Hepatoplasmataceae orthogroup | 5 records、similarity links、captions |
| BGC | 5 records、curated color categories、comparison links |
| Majanivirus | 9 records、627 matches、3 custom captions、committed filters |
| WSSV | 20 ordered conservation rings、labels/colors、cache-only replay |

### 7.5 Refresh、transaction、package tests

- unfiltered refreshのcomplete input setはpublic 11 + fixture 2であることをassertする。
- 全13のprepare、CLI replay、finalize、staged validationが成功するまでtracked targetを置換しない。
- invalid Result/catalog、cache、resource bindingは`_load_gallery_refresh_source()`型の削除/retryで通らず、source admissionで失敗する。
- refreshがcurrent sourceへ`normalize_current_session_artifacts()`をpre-cleanupとして適用せず、Vibrio cacheや`cliInvocation`のsemantic policyをPython helperが決めない。
- chloroplast staged bootstrap後、final replacement前にfailureを注入し、invalid tracked originalを含む全targetが未変更であることを確認する。
- validation failureと置換途中のcaught exceptionを注入し、変更済みtargetがrestoreされることを確認する。
- public assetsはauthoritative inventoryから一回だけ生成され、orphan/missing fileを許さない。
- GitHub Pages/Cloudflare packageのmanifest/hashがchecked-in bytesと一致する。
- hosted builder実行後にGallery targetのworktree diffがない。

### 7.6 CI placement

高コストbrowser setupを無関係なPRで4重に実行しない。

- Gallery/session/publication/resource/generator/workflowの関連path変更時は、DOM-free all11 gateを必須にする。
- 同じ関連path変更時は、影響するcircular/linear browser groupと専用WSSV/Vibrio jobをchange scopeに従って実行する。
- PR 2最終head、`dev`/main integration、predeployではpublic all11の実browser replayを必須にする。
- predeployで既存結果を再利用する場合は、deploy workflowがrequired job名とexact deployed SHAの成功を機械的に検証する。該当結果がなければall11を実行するかdeployをblockし、手書きevidenceやbranch-level greenで代替しない。
- Cloudflare production triggerにも同じexact-SHA prerequisiteを適用し、hosted buildがgateより先に走る自動経路を残さない。
- prepared browser wheelをjob間で再利用できる場合はartifact化し、各jobが同じbuildを検査する。
- WSSV/Vibrioはresource isolationのため専用processを維持するが、同じassertionを別specへ複製しない。

## 8. Generated artifact and media policy

PR 2で、明示的なowner commandから次をstaging更新する。

- `gbdraw/web/gallery/sessions/`
- `gbdraw/web/gallery/sources/`
- `gbdraw/web/gallery/examples/`
- `gbdraw/web/gallery/thumbnails/`
- `gbdraw/web/gallery/examples.json`
- `tests/test_inputs/AP027280_comparison.gbdraw-session.json`
- `tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json`

production、test、workflow、generated artifactのdiffを別々にreviewする。gzip sessionは解凍後のsemantic diff、resource inventory、request equivalenceも確認する。source/result/example/thumbnailの視覚差はreadable scaleで確認する。

PR 3では対象5例を優先しつつ、全11のstrict tutorial metadata checkを実行する。

- `HmmtDNA_ATskew`
- `majanivirus_orthogroup`
- `vibrio-harveyi-group-collinear`
- `WSSV_genome_comparison`
- `tobacco-chloroplast`

Vibrio tutorial/captureがsessionの`none` planを補正actionで隠す状態を廃止する。published session自体がcomparison intentを復元し、capture metadataはそのexact stateをassertする。

画像はstate、geometry、captionが変わったoperationだけ再captureする。同じ表示寸法でold/newを比較し、単なるsession JSON byte変更を理由に全mediaを無差別再撮影しない。

## 9. Risks and mitigations

| Risk | Mitigation |
| --- | --- |
| 一般sessionの未生成draftを消す | Gallery専用coordinatorに限定し、`session-draft-authority.test.mjs`をnegative guardにする |
| invalid currentをcleanupしてcontractを曖昧にする | currentはvalidate-first/reject。chloroplastだけone-time producer migrationし、最終codeからrepair pathを除去する |
| JavaScript/Pythonが再び意味決定を分担する | JavaScriptがprepare/finalizeとrequest canonical formを所有し、Pythonはprocess/I/O/staging/restoreに限定する |
| comparatorが重要fieldを除外する | exclusion allowlistを小さく固定し、追加時はnegative testを要求する |
| resource ID差だけで誤検出する | `{kind,payloadSha256}`へ正規化する |
| comparison topologyをraw argvから誤推測する | committed topology/analysis metadataをauthorityにし、表現不能時だけstructured recipeを設計する |
| WSSV/VibrioがCI memoryを圧迫する | 別process、worker 1、caseごとのfresh context、failure時だけtraceを保存する |
| publicationが巨大sessionを複製する | bulky resultをprojection前に分離し、chunk hashと明示的25% regression gateを置く |
| refreshが部分更新する | staging後に置換し、捕捉例外時restoreをfailure-injection testで証明する。crash atomicとは呼ばない |
| inventoryが再びdriftする | `EXAMPLES`をsourceにし、generated manifest/physical filesをread-only verifierで照合する |
| hosted buildだけ別artifactを作る | 両deploy pathをread-only化し、build前後hashとworktree diffを検査する |
| tutorial capture actionが壊れたsessionを隠す | capture前にexact session stateをassertし、補正actionを除去する |
| 過度な単純化でreal-browser回帰を失う | DOM-free all11、actual browser all11、WSSV/Vibrio isolation、curated anchors、tutorial auditをそれぞれ独立した目的で維持する |

## 10. Global acceptance criteria

- [ ] v27–30 CLI-only、v31–33/39 canonical legacy、valid current、invalid currentのversion matrixがfocused testで固定される。
- [ ] invalid current v40にsilent cleanup pathがなく、chloroplastのone-time migration helperが最終headに残らない。
- [ ] valid v40のidentityはadmission/version migrationに限定され、明示的Gallery publication rewriteとordinary importの境界がtestで固定される。
- [ ] ordinary restoreとGallery publicationが同じDOM-free current-active-config contractを使い、migration moduleはhistorical-onlyである。
- [ ] `gbdraw/web/CLAUDE.md`と必要最小限の`tools/web-change-policy.json`がfinal owner/import graphを表し、architecture checkerが成功する。
- [ ] 全11 sessionがreal current Web importerを通る。
- [ ] 全11 sessionのfirst Generateが成功する。
- [ ] 全11 sessionのcommitted/next-request publication gateが成功する。
- [ ] 全11 sessionのinitial/generated SVG semanticsが共有Python comparatorで一致する。
- [ ] `session-regenerate-intent`を含むbrowser testに独立したSVG semantic OR/fallbackが残らず、raster/direct-edit evidenceは別目的として区別される。
- [ ] request equivalence ownerが一つで、resource ID以外のsemantic差を拒否する。
- [ ] WSSVは20 FASTAをlazy resource viewのまま読み、LOSAT fallbackなしで20 ringsを再生成する。
- [ ] tobacco chloroplast current artifactに廃止2 fieldがない。
- [ ] AT skew、Majanivirus、Vibrioの観測済み差が解消する。
- [ ] Vibrio browser testがcomparison消失ではなくparityを期待する。
- [ ] normal session draft divergence contractが維持される。
- [ ] Python refresh helperがactive state/resource/requestのsemantic policyを所有しない。
- [ ] `_load_gallery_refresh_source()` cleanup/retryがなく、named Python cleanup/cache/provenance pathsが削除またはJavaScript ownerへ移管される。
- [ ] public inventory ownerが`EXAMPLES`一つで、fixture 2件は別setとして明示される。
- [ ] 無引数refreshの全13 sessionがstaging、semantic review、caught-exception restore testに含まれる。
- [ ] chloroplast bootstrapはtracked sourceを先行置換せず、同じ外側transactionのstaged overrideとして失敗時未変更を証明する。
- [ ] performanceはSHA-pinned clean worktreeとdisposable outputで測り、median wall timeとmaximum RSSが同一runner baselineの1.25倍以下である。
- [ ] Gallery generator-owned assetsが一回の明示owner commandから再現する。
- [ ] GitHub PagesとCloudflareの双方がGalleryを再生成せず、検証済みchecked-in bytesを配布する。
- [ ] 両production deploy経路がexact deployed SHAのall11 browser成功を機械的にrequireするか、packaging前に同gateを実行する。
- [ ] 対象tutorial/capture metadataが修復済みexact sessionを表す。
- [ ] `tests/reference_outputs/`と`examples/gbdraw_social_preview.png`は変更されない。
- [ ] production/test/docs/workflow/generated-artifact diffを別々にreviewする。

## 11. Plan status and evidence ledger

この表は実装中に更新し、観測していない項目を完了にしない。本文とPR指示書の改修は実装完了evidenceではない。

| Work item | Status | Required evidence |
| --- | --- | --- |
| PR 1 WSSV owner-boundary fix | Pending | focused JS + retained dedicated real-browser WSSV Generate |
| PR 2 strict v40 contract | Pending | legacy/current version matrix + one-time migration evidence + no retained repair path |
| PR 2 publication owner convergence | Pending | all11 DOM-free parity + owner/path declaration + negative equivalence tests |
| PR 2 refresh/artifact transaction | Pending | bootstrap-staged all13 refresh + pre-replace failure + caught-exception restore + generated semantic/visual review |
| PR 2 browser/CI/deploy gates | Pending | public all11 replay + isolated WSSV/Vibrio + exact-SHA deploy prerequisite + both hosted builders read-only |
| PR 2 performance gate | Pending | clean-worktree/disposable-output same-runner 3-trial baseline/head measurements within threshold |
| PR 3 tutorial/media audit | Pending | strict all11 metadata checks + affected old/new visual review + desktop/mobile |

## 12. Final verification command set

以下は計画時点のcommand contractであり、未実装のtest名を成功済みとは扱わない。各PRのfocused commandを先に実行し、PR 2とPR 3の指示書で確定した実名へ更新する。

```bash
node --test \
  tests/web/gallery-session-publication.test.mjs \
  tests/web/wssv-gallery-fastas.test.mjs \
  tests/web/session-draft-authority.test.mjs \
  tests/web/gallery-session-migration.test.mjs

python -m pytest \
  tests/test_gallery_session_semantics.py \
  tests/test_refresh_gallery_sessions.py \
  tests/test_wssv_gallery_fastas.py \
  tests/test_gallery_capture_contracts.py \
  -q

python tools/prepare_browser_wheel.py
npm run test:web:gallery-replay
npm run test:web:wssv-generate
npm run test:web:vibrio-generate

python tools/capture_gallery_tutorial_screenshots.py --all --check --strict
npx playwright test tests/web/gallery-tutorial.playwright.spec.js --project=chromium

node tools/check-web-change-budget.mjs
pytest tests/test_output_comparison.py::TestOutputComparison -v
pytest tests/ -v -m "not slow"
```

長時間testは30分以上の余裕を持って監視し、timeout assertionを短縮しない。Chromium sandbox failureは同じcommandを必要な権限で再実行する。
