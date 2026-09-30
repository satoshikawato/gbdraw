<!-- Raw design report of workstream W8 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

## Tier map (applies to every finding below)
- **PR gate** (`.github/workflows/test.yml`, only when `ci-impact` selects the job):
  - `web-contracts-pr` (:374) runs `node --test` over `tests/web/*.test.mjs` at maxdepth 1 (:410-415).
  - It also runs `pytest -m "browser and not slow"` (:417-420). That includes `tests/test_linear_comparison_browser_contracts.py:13-20`, which runs `npm run test:web:comparison-contracts`. So every `@comparison-contract` Playwright test is in the PR gate.
  - `core-pr` runs pytest with `-m "not slow and not (recipe or gallery or browser)"`.
  - `@pr-smoke` runs via `npm run test:web:pr-smoke` (:458-459).
- **After merge only** (push or dispatch on dev):
  - Playwright functional-full (:470-478) runs every other `*.playwright.spec.js`, including `tests/web/contracts/*.playwright.spec.js`.
  - `tests/run_losat_cache_browser_acceptance.py` (:696-734).

## Class 5

**CO-02: LOSATP reuse ignores Feature visibility**
- **How the code behaves:** at run-analysis.js:3069-3101, when `canReuseResolvedProteinArtifacts` is true, `useLosat` becomes false and extraction is skipped entirely. The protein extraction key does include `featureVisibility` (run-analysis.js:3330), but that layer never runs.
- **Tests that exist:**
  - `tests/web/run-analysis-simple-path.test.mjs:1771-1825`: only a positive reuse check, "resolved protein artifacts must bypass further LOSAT execution". The other `losatCalls` assertions (1482, 1651, 1760, 1769, 1881, 1930, 2201) also only check that no new call happened. There is no negative case for `canReuseResolvedProteinArtifacts` and no direct unit test of it. Tier: PR.
  - `tests/web/run-analysis-derived-cache.test.mjs`: derived identity is checked for memberMaxHits plus 10 Collinear params (:160-189), mode scoping (:191-224), Pairwise maxHits (:226-233) and pathRepresentation (:241-244). Raw identity is checked to ignore appearance-only settings (:50-62). Not covered: thresholds, `recordPayloads` (proteinCacheKey, displayBindingHash, viewTransform) and visibility. Tier: PR.
  - `tests/web/losat-cache.test.mjs:328-335`: raw-entry rejection for schema, program, outfmt, args, searchContext, protein-set hashes, record-instance keys, runtime-binding hashes and text integrity. Visibility is not a dimension. Tier: PR.
  - `tests/web/linear-multi-record.playwright.spec.js:2848` (the OIC-019 owner, `@comparison-contract`, PR): varies memberMaxHits (2983), candidateLimit (2996, 3004, 3012), clearLosatCache (3006), primary-file content (3014-3021), Collinear params (3058-3076), member limit (3078), Pairwise maxHits (3091-3107), and maxHits=0 rejection (3109-3129). No visibility change. Its LOSAT is a stub.
  - `tests/web/orthogroup-computation-cache.test.mjs:12`: a UI group-lookup memo, not related to the LOSAT cache.
  - Feature visibility tests (`feature-visibility*.test.mjs`, `source-visibility-reconciliation.playwright.spec.js`, `test_web_render_preparation_cache.py:393`) never touch LOSAT or protein.
- **Answer to your question:** no test has the form "reuse result == fresh, cache-cleared run result" after changing an input. The nearest are `explicitDefaults.svg === omittedDefaults.svg` (2852:3054), which checks default equivalence, and `restored.geometry === reversed.geometry` (3926), where both runs are reuse runs.

**CO-03: Orient forward + LOSATN ribbons wrong**
- `feature-popup-record-rotation.test.mjs:110-135`: draft defaults, with `orientForward` false (:116). No comparisons. Tier: PR.
- `feature-record-rotation.test.mjs:58-116`: `orientForward: true` projects to `presentation.reverseComplement` in the request. The runner is stubbed and there are no comparisons (the file mentions LOSAT 0 times). Tier: PR.
- `record-display-options.test.mjs:83-97`: the override has one owner in the request transform only. Tier: PR.
- `session-request.test.mjs:2011-2029`: the override becomes `presentation.reverseComplement` on the request. No comparisons checked. Tier: PR.
- `run-analysis-simple-path.test.mjs:729-762`: sets `reverseComplement: true` with a stub executor returning `[]`, and asserts `targetLosatExecutorJobs === 0`. This enshrines "no re-search" but never checks where the comparison lands. Tier: PR.
- `linear-multi-record.playwright.spec.js:3812-3935`:
  - LOSATP only, with a stub that returns identical rows.
  - Popup rotation never toggles Orient forward.
  - Reversal is done through `seq.region_reverse` (3906), not `reverseComplementOverride`.
  - Geometry is only checked with `not.toBe` (3889, 3911).
  - Tier: `@comparison-contract`, PR.
- `interactive-svg-v3.playwright.spec.js:60, 339, 430`: rotation UI only, no Orient checkbox, no LOSAT. Tier: after merge only.
- Python tests use supplied source-coordinate rows and are not the Web LOSAT pipeline:
  - `test_record_display_comparisons.py`: :155-169, :280-293, :335-372, :424
  - `test_alignment_direction_projection.py`: all tests
  - `test_reverse_complement_feature_identity.py:152`: LOSATP ids
- **Answer to your question:** no test rotates or orients a record, runs a nucleotide comparison, and checks that ribbons meet the homologous region.

**CO-04: one file = one LOSAT batch (packaging dependence)**
- `linear-sources.test.mjs:190-238`: batch counts (4 jobs for 8 records), result split, searchContext changes with database contents (:218-219), reorder invariance (:220-222), conflicting gencodes (:224-230), no-self (:232-238). All use stub TSV. Tier: PR.
- `linear-multi-record.playwright.spec.js`: 3504 (real Wasm, but all 5 records share one random sequence and thresholds are set to 0), 3538 and 3596 (job shapes and SVG pair presence). Tier: `@comparison-contract`, PR.
- `linear-comparisons.test.mjs`: plan resolution only.
- `tests/test_linear_multi_record_comparisons.py`: render-level only.
- **Answer to your question:** no packaging-invariance test (same records in one file vs separate files). No test combines `--max-target-seqs` with multi-record batches and checks hit parity.

**CO-06: uploaded tables assigned by order, not ID**
- `tests/test_linear_multi_record_comparisons.py:35-37`: the fixture rows use `q`/`s`, while the records are `r1..r4` (:40-44). Tests at :455-463, :465-505 and :507-540 therefore enshrine assignment that ignores IDs.
- `test_comparisons.py` (:130-498): no record-ID checks.
- `match-sequences.test.mjs:511`: uses ID-based matching, but only for Circular conservation companions.
- No test uploads a table whose qseqid/sseqid disagree with the slot's records. Tier: core, PR.

**CO-10: cropped coords in match popup and FASTA vs original coords in feature popup**
- `match-sequences.test.mjs:375-391`: header format, uncropped sources. Tier: PR.
- `match-sequences.test.mjs:393-432`: for a cropped + reverse-complemented region, it asserts the restored source is the view sequence `'CGGT'`, which enshrines view-frame sources. Tier: PR.
- `test_record_display_interactive_browser.py:69-205`: covers reverse and start rotation, and compares FASTA bodies and rotation invariance. No crop, and no comparison with feature-popup coordinates. Tier: browser, PR.
- `test_pairwise_match_popup.py:17`: display IDs and section titles only. Tier: browser, PR.
- `interactive-svg-search-performance.playwright.spec.js:452-461`: header format for Circular, uncropped.
- `Query interval` appears in no test.

**CO-01: threaded LOSAT without cross-origin isolation fails as UNKNOWN**
- The Playwright webServer is plain `python3 -m http.server` (playwright.config.js:18-22) and sends no COOP/COEP headers. `gbdraw/cli.py:38-39` and the Cloudflare `_headers` (test_web_packaging.py:979-980) do send them.
- Threaded tests add the headers themselves with `page.route` and assert `crossOriginIsolated === true`: `losat-thread-fault.playwright.spec.js:41-70` (after merge only) and `linear-multi-record:3937-3974` (`@comparison-contract`).
- Unit tests force `globalThis.crossOriginIsolated = true`: `losat-diagnostics.test.mjs:35`, `losat-thread-fault.test.mjs:118`.
- Every other test that runs real LOSAT avoids the default path in one of three ways:
  - It sets `executionMode = 'serial'`: `linear-multi-record` 1805, 2914, 3177; `generation-feedback.playwright` 190, 230; `similarity-alignment-ui` 2044, 2181; `run_protein_local_browser_acceptance.py:317`.
  - It stubs `__GBDRAW_LOSAT_EXECUTOR__`.
  - It loads a gallery session saved with `"auto"` (9 of 10 sessions).
- `linear-multi-record:3763` asserts the fresh default is `'threaded'`. `losat-diagnostics.test.mjs:106` asserts that explicit threaded fails with no serial fallback. Since the default is threaded, the default behaves as explicit.
- **Answer to your question:** no test runs the default mode with `crossOriginIsolated=false`. I checked that `normalizeUserFacingError(new Error('Cross-origin isolation is not enabled.'))` returns `UNKNOWN`.

## Class 6

**X-01 and related: known messages shown as UNKNOWN**
- `error-normalization.test.mjs` is table-driven over mapped messages only: :108-135 (18 `requireCurrent*` validators), :195-211, :216-232, :235-244, :252-261.
  - :29-34 asserts unmapped input maps to UNKNOWN, for privacy.
  - The comment at :105-106 claims "changing a validator's failure template must be detected" but it only applies to the listed validators.
  - Tier: PR.
- `test_web_error_adapter.py:71-99, 160-208`: parametrized with hand-picked examples. Tier: core, PR.
- `error-boundary.playwright.spec.js:77`: one injected REGEX_SYNTAX case plus INPUT_REQUIRED. `helpers/operation-error.cjs` and `helpers/structured-error-oracle.py` also only use regex. Tier: `@pr-smoke`.
- Tests assert the raw message but never pass it through normalization:
  - `session-active-files.test.mjs:213-225` ("choose a Record…"; normalizes to UNKNOWN)
  - `session-request.test.mjs:3915-3924` (Feature Width, from circular-track-slots.js:772; UNKNOWN)
- "Invalid session file." (config.js:1189, 1588): no test; normalizes to UNKNOWN.
- Identity 150: `mode_profiles.py:43` raises "identity must be a finite number in [0, 100]."; I checked it serializes as VALIDATION_UNCLASSIFIED. `test_mode_profiles.py:307-345` (identity 101 at :320) only asserts that ValidationError is raised, not how the adapter classifies it.
- SE-05 (legacy v39 session save):
  - `contracts/current-session-lazy-materialization.playwright.spec.js:1124-1180` runs Generate (:1142-1146) before Save (:1149). Save without a prior Generate is only tested for the synthetic current schema (:257-311). Tier: after merge only.
  - `composition-layout.playwright.spec.js:260` only loads.
- PV-04 (legend rename collision): the message comes from entry-actions.js:186, 203, 735 and maps to UNKNOWN. No test uses the text "different color". `feature-color-actions.test.mjs:739` is a plain rename.
- SE-06 and IN-04: not assessed. I had no description, and they appear in no doc outside the audit folders.

**Counts (gbdraw/web/js, excluding vendor)**
- 744 `throw new Error(` sites in 68 files.
  - 407 distinct literal messages; 72 are mapped by `nativeValidation` and 335 are not.
  - 259 distinct interpolated (template) messages.
  - 24 throws with non-literal arguments.
- `NATIVE_VALIDATIONS` has 102 exact entries (counted at runtime), plus 35 regex templates in `nativeValidation` (:264-334).
- Python: about 1,632 `raise ValidationError/ValueError/ParseError/TypeError` sites, with 873 distinct plain literals. Only 47 of those avoid VALIDATION_UNCLASSIFIED. The adapter has 33 `_EXACT`, 26 `_TEMPLATES` and 24 `_CONSTRAINTS` entries.
- **Answer to your key question:** no test lists every throw or raise message and checks that it maps to a known code.

**X-02: invalid numbers silently replaced by defaults**
- `normalizeBlastThresholdNumber` (run-analysis.js:885-891; applied at :2716-2731 and :2997-3010) replaces invalid values in place with defaults. Evalue only goes through `normalizeBlastThresholdText`, so `Number('1e-50x')` becomes NaN (session-request.js:2463). `optionalPositiveInteger` (:561-564) turns window/step 0 into null (Auto). No test covers any of these.
- Invalid-value tables exist only for single fields:
  - `optional-positive-number.test.mjs` + `fixtures/optional-positive-number.json` (comparison_height only)
  - `track-slot-validation.test.mjs:723-745` (pixel-track-inputs.json, 4 slot pixel fields)
  - `current-option-values.test.mjs` (enums and integers)
  - `run-analysis-derived-cache.test.mjs:14-19` (candidate limit)
- `test_mode_profiles.py:307-345` rejects invalid thresholds at the Python model level. The Web layer coerces them before they reach Python.
- The dormant-invalid fixtures (`session-request.test.mjs:4513-4519`, `run-analysis-simple-path.test.mjs:2022`) cover only the "no comparisons" path, where the values are ignored.
- `run-info.test.mjs` never checks `--evalue`, `--identity` or `-w`.
- `failed-generate-bindings.playwright.spec.js` and `generation-feedback.test.mjs` do not touch numeric inputs.
- **Answer to your question:** no table-driven test covers every numeric option.

## Contract (OPTION_INTEGRITY_PRODUCT_CONTRACT.md)
- OIPC-C01 (:285-291): "Invalid explicit values are rejected, not silently coerced." It names no regression owner.
- OIPC-C03 (:299-302): every accepted public value must reach its consumer or be rejected before execution. It names no regression owner.
- OIC-021 (:2343, AC-01..20 at :2372-2395) names observations but no test files. The evidence claims are in `docs/internal/issue-563-feature-popup-record-rotation/FINAL_ACCEPTANCE.md:110-111`: AC-14 is "LOSATP 4->4" and AC-15 is "existing transform contracts". No LOSATN or Orient-forward comparison geometry test exists, so AC-05 and AC-15 are missing for the audited case.
- Only OIC-015 (:2397-2446), OIC-016 (:2457), OIC-017..019 (:2471-2482, file owners) and OIC-020 (:2355) name regressions. OIC-001..014 and OIC-022..026 have none.
- Required regressions that are missing or do not cover the audited case:
  - **OIC-019** ("raw-setting/input changes … prevent incompatible reuse"): its owner test has no visibility input and does not exercise the resolved-artifact bypass (CO-02).
  - **OIC-015**:
    - The display-reverse clause is tested with LOSATP and `region_reverse` only (CO-03).
    - The contract *requires* source-file batching ("5 records require four source jobs") and says nothing about packaging invariance (CO-04).
    - It requires fresh Execution to be `threaded` but says nothing about fallback when the page is not cross-origin isolated (CO-01).
  - **OIC-005**: cache identity has no named owner.
  - **PD-OI-046** (:1833-1875): "keep correction info for every known validation", with no named regression.

## Systemic patterns
1. Tests check that a cache is reused (call count unchanged, zero jobs). They rarely check that reuse is invalidated, and never compare reuse output with a fresh run. LOSAT is almost always stubbed with fixed rows, so geometry correctness is invisible.
2. Tests use test-only side paths (`region_reverse` instead of the draft override, `executionMode='serial'`, forced `crossOriginIsolated`, `route`-injected COOP/COEP), which skip the production default path.
3. Tests assert raw `Error.message` text at the producer, while normalization is tested only against its own allowlist. The two are never joined by a census of all messages.
4. Numeric validation is tested one field at a time. The generic coercion helpers have no invalid-value corpus.
5. Several contract clauses enshrine the audited behaviour (batching, default threaded, zero jobs on transform) without the matching correctness requirement.
6. Rotation, v39 Save and thread-fault browser tests only run after merge.
