# Session 03 Instruction Prompt — Record state、canonical request、Session compatibility

## Prompt

あなたはgbdraw Issue #563のSession 03を担当する。過去の会話は前提にしない。

### Branchとprecondition

必ず`issue-563-feature-popup-record-rotation-20260922`を使用する。別branchを作らない。
Session 01のpreflightがruntime実装を許可し、Session 02のsource profile helperとpure JS resolverが
focused testsをpassしていることを確認する。不足があれば同じownerで補い、parallel resolverを作らない。

branch、HEAD、upstream、statusを確認し、無関係な差分を保持する。`dev/main`へ直接commitしない。

### 最初に全文を読む

1. `AGENTS.md`
2. `CLAUDE.md`
3. `gbdraw/web/CLAUDE.md`
4. `docs/internal/issue-563-feature-popup-record-rotation/IMPLEMENTATION_PLAN.md`
5. `BASELINE_AND_PREFLIGHT.md`
6. Session 02で変更されたproduction/tests
7. `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`
8. `docs/internal/PRODUCT_IMPACT_RATCHET.md`

### Session goal

per-record origin/orientation/provenanceを一つのrecord display stateへ統合し、canonical request schema 7へ
正しくprojectする。feature catalog schema 4とSession v44を一度に導入し、v43を安全に読む。
popup UIとcandidate execution/historyはまだ実装しない。

### State contract

`recordDisplayDrafts`を次の意味へ拡張する。

```text
scope, sourceUid, selector, recordId,
topologyOverride, startCoordinate,
reverseComplementOverride, anchorIntent
```

- `reverseComplementOverride`はcomplete recordのabsolute booleanまたは継承を示すnull。
- cropped recordのreverseはregion ownerのままにし、このoverrideを適用しない。
- `anchorIntent`はschema 1のprovenanceで、rendering authorityではない。
- rendering authorityはcanonical requestのabsolute start/orientationだけ。
- manual startまたはorientation editはanchorIntentをclearする。
- effective orientation precedenceは一関数/ownerだけで定義する。

### 実装作業

1. Session 02のPython source profileをfeature metadata/catalogへ接続し、feature catalog schemaを3から4へ
   上げる。JS admission、compaction、tests、fixturesを同時に更新する。
2. Session versionを43から44へ上げる。v43/catalog 3のpositive fixtureを保持する。
   - 保存済みResultを壊さない。
   - single exactなど既存dataから安全に証明できるcapabilityだけをmigrationする。
   - それ以外は再Generateが必要というexplicit disabled stateにする。
   - future/invalid schemaは引き続きrejectする。
3. `record-display-options.js`のvalidation、draft creation、reconciliation、manual edits、effective stateを更新する。
4. `session-request.js`のcanonical builderでrecord別start/orientationをprojectする。
5. `cardinality: all`はrecord間のdisplayまたはorientation差があるときだけexact-one recordsへmaterializeする。
6. materialization時にrecord order、selector、gridRow、source input index、tracks、depth binding、comparison
   endpoint、resource IDsを保持する。
7. current Session save/load projectionへabsolute transformとanchorIntentを追加する。
8. existing source-level reverse controlsはeffective orientation ownerへのadapterにする。互いに競合する
   second ownerを作らない。cropped regionだけは既存region-bound reverseを維持する。
9. sidebarの既存feature shortcutをSession 02 resolverへ置き換え、global selectionから明示targetへ変換する
   adapterだけを残す。

### Migration/architecture evidence

- canonical request schemaは7のまま。
- Worker protocol、Python request codec、renderer pathは変更しない。
- feature catalog 3 -> 4とSession 43 -> 44のnamespace、reader、fixture、removal conditionを
  Session 01のarchitecture ledgerへ反映する。
- schema numberを上げずにfieldの意味を変えない。
- compatibility branchを複数fileへ散らさず、既存migration/admission ownerへ置く。

### Required tests

- exact draft key/validationとanchorIntent schema
- manual start/orientationでprovenance clear
- complete record versus cropped reverse ownership
- multi-record one-sourceでtargetだけorientation/startが異なるmaterialization
- no difference時は不要なmaterializationをしない
- duplicate record IDsとsourceUid/selector isolation
- tracks/comparison/resource bindingsの保存
- schema 4 current catalog、catalog 3 legacy safe admission、invalid/future rejection
- v44 save/fresh Load/regenerate
- v43 positive fixture migrationと保存済みResult保持
- request schemaが7であること

### Verification

```bash
python -m pytest tests/test_web_feature_metadata.py tests/test_web_feature_catalog.py tests/test_session_io.py -v
node --test tests/web/feature-catalog.test.mjs
node --test tests/web/record-display-options.test.mjs
node --test tests/web/session-request.test.mjs tests/web/session-export-validation.test.mjs
node --test tests/web/session-authority.test.mjs tests/web/session-draft-authority.test.mjs
node tests/web/architecture-contracts.test.mjs
ruff check gbdraw/
git diff --check
```

### 禁止事項

- popup component、Apply action、Worker execution pathの追加
- anchorIntentからrender時に座標を再計算すること
- legacy catalogのambiguous partsをguessすること
- request schema 8、別request builder、別record orientation collection
- target以外のrecordをmaterialize/変更する無条件処理

### Handoff

current/legacy schema matrix、effective orientation precedence、materialization前後の例、tests、architecture
ledger deltaを報告する。commit/push/PRは呼出時のuser指示がある場合だけ行う。
