// Load of a Session 27-44 (DESIGN 4.3): the stored figure shows at once,
// Python moves the file to a Session 49 in one Worker call, and Load reads
// Python's document as it reads a Session 49 file, with the stored Results
// a Session 49 cannot hold shown in memory. A Worker failure installs nothing.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { evaluatePython, pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { body: { appendChild: () => {} }, createElement: () => ({ addEventListener: () => {}, click: () => {}, parentNode: null }) };
globalThis.alert = () => {};
const { installFakeSvgDom } = await import('./fake-svg-dom.mjs');
installFakeSvgDom();
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { state } = await import('../../gbdraw/web/js/state.js');
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}.gbdraw-session.json.gz`, import.meta.url))));
const pythonUpgrade = async (document) => JSON.parse(await evaluatePython(`
from gbdraw.web_support.legacy_session import legacy_session_to_drawings_json
result = legacy_session_to_drawings_json(payload)
`, JSON.stringify(document)));
const workerCalls = [];
setLegacySessionUpgrader(pythonLegacySessionUpgrader(workerCalls));
const load = (session) => importSession({
  target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
}, {
  confirmLoadRepairs: async () => true,
  // The composition root's port; these Results need no change.
  transformLegacyResultSvg: Object.assign(() => false, { appliesToContent: () => false })
});
const drawingsShown = () => state.drawingList().map(({ id, mode }) => ({ id, mode }));

test('a Session 44 Load makes one upgradeLegacySession call and installs the drawings of Python\'s Session 49', async () => {
  const session = fixture('HmmtDNA_basic_circular.v44-schema8');
  const expected = await pythonUpgrade(session);
  workerCalls.length = 0;
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(workerCalls.length, 1);
  assert.deepEqual(drawingsShown(), expected.document.drawings.map(({ id, mode }) => ({ id, mode })));
  assert.equal(state.activeDrawingId.value, expected.document.ui.activeDrawingId);
  assert.equal(state.results.value.length, 1);
});

test('a Session 39 Load shows the stored Results Python returns beside its document, in the drawing that opens', async () => {
  const session = fixture('BGC0000708-BGC0000713.v39');
  const expected = await pythonUpgrade(session);
  workerCalls.length = 0;
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(workerCalls.length, 1);
  assert.equal(state.activeDrawingId.value, expected.document.ui.activeDrawingId);
  assert.deepEqual(Object.keys(expected.storedResults), [expected.document.ui.activeDrawingId]);
  assert.equal(state.results.value.length, expected.storedResults[expected.document.ui.activeDrawingId].length);
});

// OV-2072: Python's conversion counts the hash= annotation targets it moved
// to their source feature (feature_1 of the Session main saved), and the
// notice after Load names them, as the Load before the Worker route did.
test('a Session 44 Load names the hash= annotation targets Python moved in its notice', async () => {
  const { ANNOTATION_TARGET_MIGRATION_NOTICE } = await import('../../gbdraw/web/js/services/feature-edit-migration.js');
  const loaded = await load(fixture('selected-feature-annotations.v44'));
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(loaded.notices, [ANNOTATION_TARGET_MIGRATION_NOTICE(1)]);
});

// OV-2000: a Session 27-39's Results are read again (feature recovery); the
// drawing that opens keeps the per-feature edits Python's conversion gave it.
test('a Session 33 Load keeps the feature edits of the drawing that opens through its feature recovery', async () => {
  const session = fixture('feature-edits-circular.v33');
  const expected = await pythonUpgrade(session);
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const opened = expected.document.drawings.find(({ id }) => id === state.activeDrawingId.value);
  const rows = (overrides) => Object.values(overrides).map(({ recordKey, biologicalFeatureId, featureVisibility, labelVisibility, labelText }) => (
    [recordKey, biologicalFeatureId, featureVisibility, labelVisibility, labelText]));
  assert.equal(state.featureExtractionPending.value, false);
  assert.ok(state.extractedFeatures.value.length > 0, 'the plan installed the file\'s feature metadata');
  assert.deepEqual(rows(state.activeDrawing().featureOverrides), rows(opened.features.featureOverrides));
  assert.equal(rows(state.activeDrawing().featureOverrides).length, 3);
});

// OV-2052: a CLI-written Session 33 with a schema-2 request keeps the record
// labels and subtitles of its `cliInvocation.args`, which Python's conversion
// gives its request, as Load before the Worker route (2d71b1611^) showed them.
test('a Session 33 schema-2 Load shows the CLI record labels and subtitles on the Linear records', async () => {
  const loaded = await load(fixture('BGC0000708-BGC0000713.schema-v2'));
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const linear = state.drawingList().find(({ mode }) => mode === 'linear');
  assert.deepEqual(linear.linearSeqs.map(({ definition, record_subtitle }) => [definition, record_subtitle]), [
    ['<i>Streptomyces lividus</i> CBS 844.73', 'Lividomycin biosynthetic gene cluster'],
    ['<i>Streptomyces fradiae</i> ATCC 10745', 'Neomycin biosynthetic gene cluster'],
    ['<i>Streptomyces fradiae</i> MCIMB 8233', 'Neomycin biosynthetic gene cluster'],
    ['<i>Streptomyces rimosus</i> subsp. <i>paromomycinus</i> NRRL 2455', 'Paromomycin biosynthetic gene cluster'],
    ['<i>Streptomyces ribosidificus</i> ATCC 21294', 'Ribostamycin biosynthetic gene']
  ]);
});

// OV-2001: a Linear Session 30's fills and Legend stay on the Linear drawing
// that opens, as Load before the Worker route (2d71b1611^) left them: the four
// caption rules color the features, the 66 saved per-feature fills are not
// restored, and the Circular drawing has the rules but no Legend entries.
test('a Linear Session 30 Load keeps its fills and Legend entries on the drawing that opens', async () => {
  const loaded = await load(fixture('BGC0000708-BGC0000713.v30'));
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(state.featureExtractionPending.value, false);
  const shown = Object.fromEntries(state.drawingList().map((drawing) => [drawing.mode, {
    active: drawing.id === state.activeDrawingId.value,
    fills: drawing.featureColorOverrides,
    rules: drawing.manualSpecificRules.map(({ cap, color, qual, val }) => [cap, color, qual, val]),
    legend: drawing.legendEntries.value.map(({ caption, color }) => [caption, color])
  }]));
  const rules = [
    ['Core biosynthetic genes', '#d03535', 'gene_kind', 'biosynthetic$'],
    ['Additional biosynthetic genes', '#f787a9', 'gene_kind', 'biosynthetic-additional'],
    ['Transport-related genes', '#577edb', 'gene_kind', 'transport'],
    ['Regulatory genes', '#57b767', 'gene_kind', 'regulatory']
  ];
  assert.deepEqual(shown, {
    linear: {
      active: true, fills: {}, rules,
      legend: [...rules.map(([caption, color]) => [caption, color]), ['other proteins', '#dddddd']]
    },
    circular: { active: false, fills: {}, rules, legend: [] }
  });
});

test('a request-less Session 30 shows its alignment target and saves a file Python reads (OV-2006, OV-2047)', async () => {
  // The stored figure is left out: a figure without a catalog is not saved
  // before a Generate (D-25).
  const session = { ...fixture('BGC0000708-BGC0000713.v30'), results: [] };
  assert.equal(session.orthogroupState.selectedOrthogroupAlignmentFeature, 'og_1');
  assert.equal((await load(session)).status, 'ok');
  assert.equal(state.selectedOrthogroupAlignmentFeature.value, 'og_1');

  const saved = await exportSession('alignment target', { askChoice: async () => 'download' });
  assert.equal(saved.status, 'saved', JSON.stringify(saved.error));
  const written = JSON.parse(gunzipSync(Buffer.from(await saved.blob.arrayBuffer())));
  // A Session 49 drawing's target is its typed request's (D-15-149): a drawing
  // that holds none writes none beside the request it is saved with.
  assert.deepEqual(written.drawings.map(({ id, renderRequest, orthogroupState }) => [
    id, renderRequest && (renderRequest.layout.similarityAlignment ?? null), orthogroupState ?? null
  ]), [['linear', null, { selectedOrthogroupId: '' }], ['circular', null, null]]);
  assert.deepEqual(JSON.parse(await evaluatePython(`
import json, tempfile
from gbdraw.session import load_session_document, materialize_session, session_to_request
document = load_session_document(payload)
with tempfile.TemporaryDirectory() as tmp, materialize_session(document, output_directory=tmp) as materialized:
    result = json.dumps([[spec.id, type(session_to_request(materialized, drawing=spec.id)).__name__]
                         for spec in document.drawings if spec.has_canonical_request])
`, written)), [['linear', 'LinearDiagramRequest']]);
  // The reloaded drawing has no target: the saved file holds none.
  assert.equal((await load(written)).status, 'ok');
  assert.equal(state.selectedOrthogroupAlignmentFeature.value, '');
});

test('a Worker failure of a Session 27-44 Load installs nothing, withdraws the stored figure and reports one diagnostic', async () => {
  const before = await load(fixture('HmmtDNA_basic_circular.v44-schema8'));
  assert.equal(before.status, 'ok');
  const drawings = drawingsShown();
  const results = state.results.value;
  setLegacySessionUpgrader(async () => { throw { code: 'HELPER_PROTOCOL', stage: 'helper', operation: 'upgradeLegacySession' }; });
  try {
    const loaded = await load(fixture('BGC0000708-BGC0000713.v39'));
    assert.equal(loaded.status, 'error');
    assert.equal(loaded.error.code, 'HELPER_PROTOCOL');
    assert.equal(loaded.error.operation, 'upgradeLegacySession');
    assert.equal(state.errorLog.value, loaded.error);
    assert.deepEqual(drawingsShown(), drawings);
    assert.equal(state.results.value, results);
  } finally {
    setLegacySessionUpgrader(pythonLegacySessionUpgrader(workerCalls));
  }
});
