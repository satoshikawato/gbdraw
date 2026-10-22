// The lift of a stored render request of schema 1 or 2 (Sessions 31-33) into
// the schema-5 request that means the same (PR-D): the arm cases and the
// request corpus of `main` history (both in the shared vector), and the
// full-config projection it feeds.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';

import { caseRequest, readVectors, valueAt } from './helpers/legacy-request-vectors.mjs';

import { liftLegacyRenderRequest } from '../../gbdraw/web/js/services/legacy-render-request.js';
import {
  projectCanonicalSessionRequest,
  readStoredRenderRequest
} from '../../gbdraw/web/js/services/session-request.js';
import { loadRepair } from '../../gbdraw/web/js/utils/load-repairs.js';

const SESSIONS = new URL('../fixtures/sessions/', import.meta.url);
const gunzipped = (name) => JSON.parse(gunzipSync(readFileSync(new URL(name, SESSIONS))));
// The shared vector: the corpus and the arm cases.
const vectors = readVectors(new URL('../fixtures/legacy-render-request-vectors.json', import.meta.url));
const corpus = vectors;
const bgcSession = gunzipped('BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz');
// Empty resources for every resource a request names.
const resourcesOf = (renderRequest) => {
  /** @type {Record<string, object>} */
  const resources = {};
  JSON.stringify(renderRequest, (key, value) => {
    if (/resourceId$/i.test(key) && typeof value === 'string') {
      resources[value] = { kind: 'canonical-tsv', name: `${value}.txt`, type: 'text/plain', size: 0, lastModified: 0, encoding: 'base64', data: '' };
    }
    return value;
  });
  return resources;
};
// The one reader of a stored request, as Load calls it (step (2) lifts).
const readStored = (renderRequest, context = {}) => readStoredRenderRequest(renderRequest, {
  featureCatalog: null, legacyOrthogroupState: null,
  resourceIds: new Set(Object.keys(resourcesOf(renderRequest))), path: 'renderRequest', ...context
});

const request = ({ schema = 2, mode = 'linear', options = {}, ...rest } = {}) => ({
  schema,
  mode,
  records: [],
  layout: {},
  comparisons: [],
  diagramOptions: { config: { canvas: {} }, configOverrides: {}, ...options },
  ...rest
});
const lifted = (value) => liftLegacyRenderRequest(value).value;

test('a value that is not a schema 1 or 2 request passes through unchanged', () => {
  for (const value of [
    null, 'text', [], {}, { schema: 5, mode: 'linear', diagramOptions: {} }, { schema: 3, mode: 'circular', diagramOptions: {} },
    { schema: 2, mode: 'grid', diagramOptions: {} }, { schema: 1, mode: 'linear', diagramOptions: null }
  ]) {
    const result = liftLegacyRenderRequest(value);
    assert.equal(result.value, value);
    assert.deepEqual(result.repairs, []);
  }
});

test('the lift writes schema 5 with an explicit grouping and leaves its input alone', () => {
  const input = request({ mode: 'circular', records: [{ recordKey: 'a' }, { recordKey: 'b' }] });
  const before = JSON.stringify(input);
  const value = lifted(input);
  assert.equal(JSON.stringify(input), before);
  assert.deepEqual([value.schema, value.grouping], [5, 'grid']);
  assert.equal(lifted(request({ mode: 'circular', records: [{ recordKey: 'a' }] })).grouping, 'single');
  assert.equal(lifted(request({ mode: 'circular', records: [{}], layout: { columns: 2 } })).grouping, 'grid');
  assert.equal(lifted(request({ mode: 'linear', records: [{}, {}] })).grouping, 'single');
});

test('forms the lift leaves to the reader stay as stored, with no row', () => {
  const slots = lifted(request({
    mode: 'circular',
    options: { tracks: { circularTrackSlots: [{ kind: 'circularTrackSlot', id: 'odd', renderer: 'features', spacing: '3 parsecs' }] } }
  })).diagramOptions.tracks.circularTrackSlots;
  assert.equal(slots[0].spacing, '3 parsecs', 'a spacing that cannot be read is left for the reader to refuse');
  const partial = lifted(request({
    comparisons: [{ kind: 'generatedProteinComparison', pairs: [{ q: 0 }], settings: { collinearityParams: { kind: 'standard', parameters: { minAnchors: 1 } } } }]
  })).comparisons[0];
  assert.deepEqual([partial.pairs, partial.settings.collinearityParams.kind], [[{ q: 0 }], 'standard']);
  assert.deepEqual(liftLegacyRenderRequest(request({ options: { configOverrides: { comparison_height: '' } } })).repairs, [],
    'a blank comparison height is unset, not invalid');
  const records = lifted(request({ schema: 1, records: [{}, { recordKey: 'kept' }] })).records;
  assert.deepEqual(records.map(({ recordKey }) => recordKey), ['record-1', 'kept']);
});

// Each case of the shared vector (tests/fixtures/legacy-render-request-vectors.json):
// a corpus request with edits, the lifted values and rows it expects. The
// cross-runtime pytest decodes the same cases in Python.
test('each arm case of the shared vector lifts to its expected values and rows', () => {
  assert.ok(vectors.cases.length >= 30, String(vectors.cases.length));
  for (const vectorCase of vectors.cases) {
    if (vectorCase.expect.refused) {
      // A form no writer stored: refused as the base refused it (D-OD13 F3, D-OD16 F2).
      assert.throws(() => liftLegacyRenderRequest(caseRequest(vectors, vectorCase)), vectorCase.expect.refused, vectorCase.id);
      continue;
    }
    const { value, repairs } = liftLegacyRenderRequest(caseRequest(vectors, vectorCase));
    const { values = [], absent = [], repairs: rows } = vectorCase.expect;
    for (const [path, expected] of values) assert.deepEqual(valueAt(value, path), expected, `${vectorCase.id}: ${path.join('.')}`);
    for (const path of absent) assert.equal(valueAt(value, path), undefined, `${vectorCase.id}: ${path.join('.')} is absent`);
    assert.deepEqual(repairs, rows.map((row) => loadRepair(row)), vectorCase.id);
    assert.equal(value.schema, 5, vectorCase.id);
  }
});

test('each request of the corpus lifts and reads without repairs to the current schema', () => {
  assert.equal(corpus.requests.length, 23);
  const promotion = {
    featureCatalog: bgcSession.editorState?.featureCatalog ?? null,
    legacyOrthogroupState: bgcSession.orthogroupState
  };
  for (const { id, renderRequest } of corpus.requests) {
    const { value, repairs } = liftLegacyRenderRequest(renderRequest);
    assert.deepEqual(repairs, [], id);
    assert.equal(value.schema, 5, id);
    assert.equal(liftLegacyRenderRequest(value).value, value, `${id}: a lifted request is not lifted again`);
    const read = readStored(renderRequest, /^BGC/.test(id) ? promotion : {});
    assert.ok(!('kind' in read), `${id}: ${JSON.stringify(read)}`);
    assert.deepEqual(read.repairs, [], id);
    assert.equal(read.value.schema, 9, id);
  }
});

test('a full config whose small tick interval Python reads as unset projects as null (OV-791)', () => {
  const { renderRequest } = corpus.requests.find(({ id }) => id === 'lambda_basic_linear@10d3a3d2');
  assert.equal(renderRequest.diagramOptions.config.objects.gc_content.small_tick_interval, 'none');
  const read = readStored(renderRequest);
  assert.ok(!('kind' in read));
  const { config } = projectCanonicalSessionRequest({
    renderRequest: read.value, resources: resourcesOf(renderRequest), deferResourceContent: true
  });
  assert.deepEqual([config.adv.gc_content_small_tick_interval, config.adv.depth_small_tick_interval], [null, null]);
});

// A resolved config always names a font list. The packaged list is the Font
// default (null), not a chosen font; another list is the font chosen (OV-1014).
test('a full config with the packaged font list projects as the default Font', () => {
  const { renderRequest } = corpus.requests.find(({ id }) => id === 'lambda_basic_linear@10d3a3d2');
  const fontOf = (request) => {
    const read = readStored(request);
    assert.ok(!('kind' in read));
    return projectCanonicalSessionRequest({
      renderRequest: read.value, resources: resourcesOf(request), deferResourceContent: true
    }).config.adv.font_family;
  };
  assert.equal(
    renderRequest.diagramOptions.config.objects.text.font_family,
    "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"
  );
  assert.equal(fontOf(renderRequest), null);
  const chosen = structuredClone(renderRequest);
  chosen.diagramOptions.config.objects.text.font_family = "'DejaVu Serif', serif";
  assert.equal(fontOf(chosen), "'DejaVu Serif', serif");
});
