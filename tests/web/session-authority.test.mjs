import assert from 'node:assert/strict';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
// R6 (R1-M3): the Save checks of a Result set, a resource and a Web binding refuse with the typed Session-fields error.
const sessionFields = (error) => error?.code === 'INPUT_INVALID' && error.context?.reason === 'FIELDS';
const {
  SESSION_TOP_LEVEL_AUTHORITY,
  isArtifactBand,
  isAuthorityBand,
  isCanonicalBand,
  isModeScopedBand,
  isDrawingsBand,
  projectArtifactState,
  projectDocumentMetadata,
  projectWebOnlyEditorMetadata,
  readSessionInput,
  validateSessionAuthorityInventory
} = await import(pathToFileURL(join(repoRoot, 'gbdraw/web/js/services/session-authority.js')));

assert.equal(SESSION_TOP_LEVEL_AUTHORITY.modes.authority, 'editor-metadata');
assert.equal(SESSION_TOP_LEVEL_AUTHORITY.cliOptions.authority, 'provenance');
assert.equal(SESSION_TOP_LEVEL_AUTHORITY.renderRequest.authority, 'canonical-render');
assert.equal(SESSION_TOP_LEVEL_AUTHORITY.resources.authority, 'resource');
assert.equal(SESSION_TOP_LEVEL_AUTHORITY.webFiles.authority, 'resource-binding');
// R4: the content part each key belongs to (null: no content part).
assert.equal(SESSION_TOP_LEVEL_AUTHORITY.results.part, 'results');
assert.equal(SESSION_TOP_LEVEL_AUTHORITY.losatCache.part, 'cache');
assert.equal(SESSION_TOP_LEVEL_AUTHORITY.title.part, null);

const session = {
  format: 'gbdraw-session', version: 39, createdAt: 'now', title: 'Canonical',
  renderRequest: {}, resources: {}, webFiles: {}, config: {}, files: {},
  ui: {
    cInputType: 'gff', lInputType: 'gb',
    mode: 'linear', legend: 'left', linearPlotTitlePosition: 'top', zoom: 1.5,
    canvasPan: { x: 3, y: 4 }, generatedLegendPosition: 'right', downloadDpi: 300,
    linearTypographyLinked: false,
    appliedPaletteName: 'orchid', appliedPaletteColors: { CDS: '#123456' },
    pendingPaletteName: 'mint', pendingPaletteColors: { CDS: '#abcdef' }
  },
  features: {
    extractedFeatures: [{ id: 'f1' }], featureColorOverrides: { f1: '#ffffff' },
    featureVisibilityManualRules: [{ action: 'off' }], labelTextFeatureOverrides: { f1: 'stored' }
  },
  editorState: { legend: { entries: [] } }, results: [{ name: 'preview', content: '<svg/>' }],
  orthogroupState: {}, losatCache: {}, losatDerivedCache: {},
  proteinIdentityManifest: {
    schema: 2,
    proteinSets: {},
    recordAnalyses: {},
    recordInstances: {}
  },
  legacyArtifacts: { proteinRawCandidates: { schema: 1, entries: [] } },
  runMetadata: { trackSlotGeometry: { schema: 1, records: [] } },
  cliInvocation: null
};
validateSessionAuthorityInventory(session, 39);
const currentSession = {
  ...session,
  version: 40,
  features: {
    featureColorOverrides: { f1: '#ffffff' },
    featureVisibilityManualRules: [{ action: 'off' }],
    labelTextFeatureOverrides: { f1: 'stored' }
  },
  editorState: {
    legend: { entries: [] },
    featureCatalog: { schema: 3, items: [] }
  }
};
delete currentSession.files;
assert.doesNotThrow(() => validateSessionAuthorityInventory(currentSession, 40));
for (const field of ['cInputType', 'lInputType']) {
  assert.throws(
    () => validateSessionAuthorityInventory({
      ...currentSession,
      ui: { ...currentSession.ui, [field]: 'unknown' }
    }, 40),
    new RegExp(`Session ui.${field} must be gb or gff`)
  );
}
const currentWebDraft = {
  ...currentSession,
  resources: {
    'comparison-resource': { kind: 'web-file' }
  },
  config: {
    linearRecordLayout: { enabled: false, recordGap: 24, rows: [] },
    linearComparisonPlan: {
      mode: 'selected',
      defaultSource: 'losat',
      edges: [{
        id: 'edge-a-b',
        queryUid: 'a',
        subjectUid: 'b',
        included: true,
        fileActive: true,
        losatFilenameActive: false,
        source: 'upload',
        losatFilename: ''
      }]
    }
  },
  webFiles: {
    bindings: {
      schema: 1,
      linearComparisons: [{
        id: 'edge-a-b',
        file: { resourceId: 'comparison-resource' }
      }]
    }
  }
};
assert.doesNotThrow(() => validateSessionAuthorityInventory(currentWebDraft, 40));
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentWebDraft,
    webFiles: {
      bindings: {
        ...currentWebDraft.webFiles.bindings,
        schema: 99
      }
    }
  }, 40),
  sessionFields
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentWebDraft,
    webFiles: {
      bindings: {
        ...currentWebDraft.webFiles.bindings,
        linearComparisons: [{ id: 'edge-a-b', file: null }]
      }
    }
  }, 40),
  /requires a file resource binding/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentWebDraft,
    webFiles: {
      bindings: {
        ...currentWebDraft.webFiles.bindings,
        linearComparisons: [{
          id: 'edge-a-b',
          file: { resourceId: 'missing-comparison-resource' }
        }]
      }
    }
  }, 40),
  sessionFields
);
// OV-636: a `resources` that holds no resource is a missing resource whatever it
// is; `null` used to end in a TypeError. Python refuses each of them
// (gbdraw/session_io.py:632-635).
for (const resources of [null, undefined, [], 5, 'x']) assert.throws(
  () => validateSessionAuthorityInventory({ ...currentWebDraft, resources }, 40),
  sessionFields,
  String(JSON.stringify(resources))
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentWebDraft,
    webFiles: {
      bindings: {
        ...currentWebDraft.webFiles.bindings,
        linearComparisons: []
      }
    }
  }, 40),
  /Active comparison file is missing its Web file binding/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentWebDraft,
    config: { linearRecordLayout: { enabled: false, recordGap: 24, rows: [] } }
  }, 40),
  /requires config\.linearComparisonPlan/
);
for (const corrupt of [
  {
    ...currentWebDraft,
    config: { ...currentWebDraft.config, blastSource: 'losat' }
  },
  {
    ...currentWebDraft,
    ui: { ...currentWebDraft.ui, blastSource: 'upload' }
  },
  {
    ...currentWebDraft,
    config: {
      ...currentWebDraft.config,
      linearRecordLayout: {
        ...currentWebDraft.config.linearRecordLayout,
        comparisons: []
      }
    }
  }
]) {
  assert.throws(
    () => validateSessionAuthorityInventory(corrupt, 40),
    /cannot contain/
  );
}
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentWebDraft,
    webFiles: {
      bindings: {
        ...currentWebDraft.webFiles.bindings,
        linearComparisons: [{
          ...currentWebDraft.webFiles.bindings.linearComparisons[0],
          queryUid: 'a',
          source: 'upload'
        }]
      }
    }
  }, 40),
  /duplicates plan metadata/
);
assert.doesNotThrow(() => validateSessionAuthorityInventory({
  ...currentSession,
  results: [],
  editorState: { featureCatalog: null }
}, 40));
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    results: 'not-an-array',
    editorState: { featureCatalog: null }
  }, 40),
  sessionFields
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    results: [],
    editorState: {}
  }, 40),
  sessionFields
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    features: []
  }, 40),
  /Session features must be an object/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    orthogroupState: []
  }, 40),
  /Session orthogroupState must be an object/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    results: [{
      name: 'preview.interactive.svg',
      content: '<svg><metadata id="gbdraw-interactive-feature-metadata"/></svg>'
    }],
    editorState: {
      featureCatalog: { schema: 3, items: [] }
    }
  }, 40),
  sessionFields
);
assert.throws(
  () => validateSessionAuthorityInventory({ ...currentSession, files: {} }, 40),
  /version 40 cannot contain legacy files/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    features: {
      ...currentSession.features,
      featureCatalog: {
        schema: 1,
        encoding: 'biological-authority-v1',
        extracted: []
      }
    }
  }, 40),
  /branch-only feature field.*featureCatalog/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    features: {
      ...currentSession.features,
      biologicalFeatures: []
    }
  }, 40),
  /branch-only feature field.*biologicalFeatures/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    orthogroupState: { groups: [] }
  }, 40),
  /cannot contain duplicated orthogroup groups/
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    editorState: { featureCatalog: { schema: 1 } }
  }, 40),
  sessionFields
);
assert.throws(
  () => validateSessionAuthorityInventory({
    ...currentSession,
    editorState: { featureCatalog: null }
  }, 40),
  sessionFields
);
assert.deepEqual(projectWebOnlyEditorMetadata(session).ui, {
  cInputType: 'gff',
  lInputType: 'gb',
  legend: 'left',
  linearPlotTitlePosition: 'top',
  zoom: 1.5,
  canvasPan: { x: 3, y: 4 },
  downloadDpi: 300,
  linearTypographyLinked: false,
  appliedPaletteName: 'orchid',
  appliedPaletteColors: { CDS: '#123456' },
  pendingPaletteName: 'mint',
  pendingPaletteColors: { CDS: '#abcdef' }
});
const layoutPreferences = {
  circular: {
    single: { legend: 'left', plotTitlePosition: 'none' },
    multi: { legend: 'right', plotTitlePosition: 'bottom' }
  },
  linear: { legend: 'bottom', plotTitlePosition: 'top' }
};
assert.deepEqual(
  projectWebOnlyEditorMetadata({ ui: { layoutPreferences } }).ui,
  { layoutPreferences }
);
assert.deepEqual(projectArtifactState(session).features, {
  extractedFeatures: [{ id: 'f1' }],
  featureColorOverrides: { f1: '#ffffff' },
  featureVisibilityManualRules: [{ action: 'off' }],
  labelTextFeatureOverrides: { f1: 'stored' }
});
assert.deepEqual(projectArtifactState(session).ui, { generatedLegendPosition: 'right' });
assert.deepEqual(
  projectArtifactState(session).proteinIdentityManifest,
  session.proteinIdentityManifest
);
assert.deepEqual(projectArtifactState(session).legacyArtifacts, session.legacyArtifacts);
assert.deepEqual(projectArtifactState(session).runMetadata, session.runMetadata);
assert.equal(projectDocumentMetadata(session).title, 'Canonical');
assert.throws(
  () => validateSessionAuthorityInventory({ ...session, unknownField: true }, 39),
  /unclassified top-level field.*unknownField/
);
assert.throws(
  () => validateSessionAuthorityInventory({ ...currentSession, unknownField: true }, 40),
  /unclassified top-level field.*unknownField/
);
assert.doesNotThrow(() => validateSessionAuthorityInventory({ ...session, unknownField: true }, 30));

const { readFile } = await import('node:fs/promises');
const frozen = JSON.parse(await readFile(join(repoRoot, 'tests/fixtures/sessions/single.v41-bindings1.json')));
assert.doesNotThrow(() => validateSessionAuthorityInventory(frozen, 41));
const compositeSession = () => {
  const session = structuredClone(frozen);
  const leaf = session.webFiles.bindings.c_gb;
  session.webFiles.bindings = { schema: 2, c_gb: { kind: 'composite',
    components: [leaf, { ...leaf, name: 'repeat.gb' }], name: '', type: '', lastModified: 0.5 } };
  return session;
};
assert.doesNotThrow(() => validateSessionAuthorityInventory(compositeSession(), 41));
for (const change of [
  s => { s.webFiles.bindings.schema = 99; },
  s => { s.webFiles.bindings.schema = 1; },
  s => { s.webFiles.bindings.c_gb.kind = 'multipart'; },
  s => { s.webFiles.bindings.c_gb.components = []; },
  s => { s.webFiles.bindings.c_gb.components.length = 1; },
  s => { s.webFiles.bindings.c_gb.components[0] = structuredClone(s.webFiles.bindings.c_gb); },
  s => { delete s.webFiles.bindings.c_gb; },
  s => { s.webFiles.bindings.c_gb.components[0] = null; },
  s => { s.webFiles.bindings.c_gb.components[0].resourceId = 'missing'; },
  s => { s.webFiles.bindings.c_gb.components[0].type = false; },
  s => { delete s.webFiles.bindings.c_gb.components[0].name; },
  s => { s.webFiles.bindings.c_gb.components[0].extra = true; },
  s => { s.webFiles.bindings.c_gb.resourceId = 'mixed'; },
  s => { s.webFiles.bindings.c_gb.size = 0; },
  s => { s.webFiles.bindings.c_gb.lastModified = -1; },
  s => { s.webFiles.bindings.c_gb.lastModified = Infinity; },
  s => { s.webFiles.bindings.c_fasta = s.webFiles.bindings.c_gb; },
  s => { Object.values(s.resources)[0].encoding = 'raw'; },
  s => { Object.values(s.resources)[0].data = null; }
]) {
  const session = compositeSession();
  change(session);
  const before = structuredClone(session);
  assert.throws(() => validateSessionAuthorityInventory(session, 41));
  assert.deepEqual(session, before);
}
for (const version of [27, 33, 39, 40, 43]) {
  assert.throws(() => validateSessionAuthorityInventory(compositeSession(), version), /requires session version 41/);
}

const { validateAnnotationWarnings } = await import(pathToFileURL(join(repoRoot, 'gbdraw/web/js/services/session-feature-metadata.js')));
const warning = { code: 'feature_selector_unmatched', setId: 's', annotationId: 'a',
  recordId: 'same', recordIndex: 1, missingCount: 2, message: 'Skipped annotation: 2 feature selector(s) unmatched.',
  resultIndex: 0, resultName: 'out.svg' };
const warningResults = [{ name: 'out.svg', content: '<svg />' }];
const warningList = [warning];
assert.equal(validateAnnotationWarnings(warningList, warningResults), warningList);
assert.deepEqual(validateAnnotationWarnings(undefined, warningResults), []);
for (const invalid of [
  null, {}, [null], [{ ...warning, unknown: true }], [{ ...warning, missingCount: 0 }],
  [{ ...warning, missingCount: 1.5 }], [{ ...warning, recordIndex: -1 }],
  [{ ...warning, recordIndex: true }], [{ ...warning, code: 'unknown' }],
  [{ ...warning, resultName: 'stale.svg' }], [{ ...warning, resultIndex: 1 }],
  [{ ...warning, annotationId: '' }], [{ ...warning, message: {} }]
]) assert.throws(() => validateAnnotationWarnings(invalid, warningResults), /Result metadata schema/);
assert.doesNotThrow(() => validateAnnotationWarnings([{ ...warning, code: 'empty_span', missingCount: 0 }], warningResults));

// PD-OI-074: unmatched comparison table IDs reach the Web as Result metadata.
const { validateComparisonWarnings } = await import(pathToFileURL(join(repoRoot, 'gbdraw/web/js/services/comparison-warnings.js')));
const comparisonWarning = { code: 'comparison_record_id_unmatched', queryRecordIndex: 0, subjectRecordIndex: 1,
  queryRecordId: 'R2', subjectRecordId: 'R3', rowCount: 2, exampleIds: ['contig_A', 'contig_B'],
  message: 'Comparison between query record #1 \'R2\' and subject record #2 \'R3\': 2 row(s) use sequence IDs that match no displayed record (\'contig_A\', \'contig_B\'); these rows are drawn on the records assigned by position.',
  resultIndex: 0, resultName: 'out.svg' };
const comparisonWarningList = [comparisonWarning];
assert.equal(validateComparisonWarnings(comparisonWarningList, warningResults), comparisonWarningList);
assert.deepEqual(validateComparisonWarnings(undefined, warningResults), []);
assert.deepEqual(validateComparisonWarnings([], warningResults), []);
for (const invalid of [
  null, {}, [null], [{ ...comparisonWarning, unknown: true }], [{ ...comparisonWarning, rowCount: 0 }],
  [{ ...comparisonWarning, rowCount: 1.5 }], [{ ...comparisonWarning, queryRecordIndex: -1 }],
  [{ ...comparisonWarning, subjectRecordIndex: true }], [{ ...comparisonWarning, code: 'unknown' }],
  [{ ...comparisonWarning, resultName: 'stale.svg' }], [{ ...comparisonWarning, resultIndex: 1 }],
  [{ ...comparisonWarning, exampleIds: 'contig_A' }], [{ ...comparisonWarning, exampleIds: [1] }],
  [{ ...comparisonWarning, message: '' }], [{ ...comparisonWarning, message: {} }]
]) assert.throws(() => validateComparisonWarnings(invalid, warningResults), /Result metadata schema/);
// The Session validator checks a saved warning against the saved Results.
assert.throws(() => validateSessionAuthorityInventory({
  ...frozen, runMetadata: { comparisonWarnings: [comparisonWarning] }
}, 41), /Result metadata schema/);

// E1: a Session keeps the other diagram mode's Result set in
// `otherModeResult`, whose fields mirror one committed top-level set
// (Session 46; the Gallery Sessions are current-writer Sessions).
{
  const { readFileSync } = await import('node:fs');
  const read = (name) => JSON.parse(readFileSync(join(repoRoot, 'gbdraw/web/gallery/sessions', name), 'utf8'));
  const circularSession = read('HmmtDNA_basic_circular.gbdraw-session.json');
  const linearSession = read('lambda_basic_linear.gbdraw-session.json');
  const otherModeResult = {
    renderRequest: linearSession.renderRequest,
    results: linearSession.results,
    editorState: {
      featureCatalog: linearSession.editorState.featureCatalog,
      alignmentResetReceipt: null,
      legend: { originalOrder: [], originalColors: {} },
      originalSvgStroke: { color: null, width: null }
    },
    ui: { selectedResultIndex: 0 },
    runMetadata: linearSession.runMetadata
  };
  const both = { ...circularSession, otherModeResult };
  assert.doesNotThrow(() => validateSessionAuthorityInventory(both, 46));
  assert.throws(() => validateSessionAuthorityInventory({
    ...both, version: 44, editorState: { ...both.editorState, featureCatalog: { ...both.editorState.featureCatalog, schema: 4 } }
  }, 44), /cannot contain otherModeResult/);
  assert.throws(() => validateSessionAuthorityInventory({
    ...both, otherModeResult: { ...otherModeResult, renderRequest: circularSession.renderRequest }
  }, 46), sessionFields);
  assert.throws(() => validateSessionAuthorityInventory({
    ...both, renderRequest: null, results: [], editorState: { ...both.editorState, featureCatalog: null }
  }, 46), sessionFields);
  assert.throws(() => validateSessionAuthorityInventory({
    ...both, otherModeResult: { ...otherModeResult, results: [], editorState: { ...otherModeResult.editorState, featureCatalog: null } }
  }, 46), sessionFields);
  assert.throws(() => validateSessionAuthorityInventory({
    ...both, otherModeResult: { ...otherModeResult, config: {} }
  }, 46), sessionFields);
  const unboundRequest = structuredClone(linearSession.renderRequest);
  unboundRequest.records[0].source.resourceId = 'record-9-genbank';
  assert.throws(() => validateSessionAuthorityInventory({
    ...both, otherModeResult: { ...otherModeResult, renderRequest: unboundRequest }
  }, 46), sessionFields);
}

// `readSessionInput` is the one reader of a parsed Session file. Its refusals
// keep their order: not an object, a version that is no integer, a version
// this app does not read, then a format that is not `gbdraw-session`.
{
  const refusal = (value) => {
    try {
      readSessionInput(value);
    } catch (error) {
      return `${error.code}:${error.context?.reason}`;
    }
    return 'accepted';
  };
  const format = 'INPUT_INVALID:SESSION_FORMAT';
  const fields = 'INPUT_INVALID:FIELDS';
  for (const value of [null, undefined, 5, 'x', true, [], [{ format: 'gbdraw-session', version: 46 }]]) {
    assert.equal(refusal(value), format, JSON.stringify(value));
  }
  for (const version of [undefined, null, '46', 46.5, true, [46], {}]) {
    assert.equal(refusal({ format: 'gbdraw-session', version }), format, `version ${JSON.stringify(version)}`);
  }
  for (const version of [50, 45, 43, 38, 34, 26, 0, -1]) {
    // B51 (Q2): a Session older than 27 is refused as an old version; OV-1253:
    // a newer one, and one in a gap, are refused as such.
    const reason = version === 26 ? 'SESSION_TOO_OLD' : version > 49 ? 'SESSION_TOO_NEW' : 'SESSION_VERSION_UNSUPPORTED';
    assert.equal(refusal({ format: 'gbdraw-session', version }), `INPUT_INVALID:${reason}`, `version ${version}`);
    if (reason === 'SESSION_TOO_OLD') continue;
    // The refusal names the versions this reader takes.
    assert.throws(() => readSessionInput({ format: 'gbdraw-session', version }), (error) => {
      assert.deepEqual(error.context, reason === 'SESSION_TOO_NEW'
        ? { field: 'schema', reason, sessionVersion: version, latestSessionVersion: 49 }
        : { field: 'schema', reason, sessionVersion: version, readSessionVersions: '27–33, 39–42, 44, and 46–49' });
      return true;
    });
    // The version is read before the format.
    assert.equal(refusal({ format: 'other', version }), fields, `version ${version}, other format`);
  }
  for (const formatValue of [undefined, null, 'other', 'GBDRAW-SESSION', ['gbdraw-session'], 46]) {
    assert.equal(refusal({ format: formatValue, version: 46 }), format, `format ${JSON.stringify(formatValue)}`);
  }
  for (const version of [27, 28, 29, 30, 31, 32, 33, 39, 40, 41, 42, 44, 46, 47, 48, 49]) {
    const document = { format: 'gbdraw-session', version, anything: { unvalidated: true } };
    assert.equal(readSessionInput(document), document, `version ${version}`);
  }
}

{
  // K16: each band predicate narrows to the bands its type names, for every
  // supported version, as the generated admission description places it.
  const { SESSION_ADMISSION_DESCRIPTION } = await import(pathToFileURL(join(repoRoot, 'gbdraw/web/js/contracts/session-shape.generated.js')));
  const { SUPPORTED_SESSION_VERSIONS } = await import(pathToFileURL(join(repoRoot, 'gbdraw/web/js/contracts/schema.generated.js')));
  const predicates = [
    [isCanonicalBand, ['SessionChecked31', 'SessionChecked39', 'SessionChecked40', 'SessionChecked46']],
    [isArtifactBand, ['SessionChecked39', 'SessionChecked40', 'SessionChecked46']],
    [isAuthorityBand, ['SessionChecked40', 'SessionChecked46']],
    [isModeScopedBand, ['SessionChecked46']],
    // A Session 49 is in none of the four (DESIGN §6: no `>=` branch admits it).
    [isDrawingsBand, ['SessionChecked49']]
  ];
  const { bands } = SESSION_ADMISSION_DESCRIPTION;
  assert.deepEqual(bands.flatMap(({ versions }) => versions), [...SUPPORTED_SESSION_VERSIONS].sort((a, b) => a - b));
  assert.deepEqual(bands.map(({ root }) => root), ['SessionChecked27', ...predicates[0][1], ...predicates[4][1]]);
  for (const version of SUPPORTED_SESSION_VERSIONS) {
    const band = bands.find(({ versions }) => versions.includes(version));
    assert.equal(band.minVersion, Math.min(...band.versions));
    for (const [predicate, roots] of predicates) {
      assert.equal(predicate({ format: 'gbdraw-session', version }), roots.includes(band.root), `${predicate.name} ${version}`);
    }
  }
}

{
  // DESIGN §6: the writer flip (SESSION_VERSION 49) leaves the settings-only
  // form and the schema-2 bindings of Sessions 46-48 as they are, and gives a
  // Session 49 no settings-only form: its drawing without a request is read
  // by that drawing (`projectsWithoutRequest`). The child reads
  // session-authority.js with the generated SESSION_VERSION set to 49.
  const { spawnSync } = await import('node:child_process');
  const run = spawnSync(process.execPath, ['--experimental-test-module-mocks', '--no-warnings',
    join(repoRoot, 'tests/web/helpers/session-authority-at-flip.mjs')], { encoding: 'utf8' });
  assert.equal(run.status, 0, run.stderr);
  const atFlip = JSON.parse(run.stdout.trim().split('\n').at(-1));
  assert.equal(atFlip.sessionVersion, 49);
  assert.deepEqual(atFlip.settingsOnly, { 46: true, 47: true, 48: true, 49: false });
  assert.deepEqual(atFlip.bindings, { 46: 'read', 47: 'read', 48: 'read', 49: 'read' });
  assert.deepEqual(atFlip.projectsWithoutRequest, { 30: false, 48: true, 49: true });
}
