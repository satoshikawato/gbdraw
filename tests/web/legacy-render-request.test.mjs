// The lift of a stored render request of schema 1 or 2 (Sessions 31-33) into
// the schema-5 request that means the same (PR-D): each arm, the request
// corpus of `main` history, and the full-config projection it feeds.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';

import { liftLegacyRenderRequest } from '../../gbdraw/web/js/services/legacy-render-request.js';
import {
  projectCanonicalSessionRequest,
  promoteCanonicalRenderRequestToCurrent
} from '../../gbdraw/web/js/services/session-request.js';
import { loadRepair } from '../../gbdraw/web/js/utils/load-repairs.js';

const SESSIONS = new URL('../fixtures/sessions/', import.meta.url);
const gunzipped = (name) => JSON.parse(gunzipSync(readFileSync(new URL(name, SESSIONS))));
const corpus = gunzipped('legacy-render-requests.vectors.json.gz');
const bgcSession = gunzipped('BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz');

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

test('the retired output prefix is dropped (A1)', () => {
  const value = lifted(request({ options: { output: { outputPrefix: 'x', legend: 'right' } } }));
  assert.deepEqual(value.diagramOptions.output, { legend: 'right' });
});

test('the feature visibility table takes its current names; a current value wins (A15, OV-792)', () => {
  const old = lifted(request({ options: { featureTable: 'T', featureTableFile: 'F' } })).diagramOptions;
  assert.deepEqual([old.featureVisibilityTable, old.featureVisibilityTableFile, 'featureTable' in old, 'featureTableFile' in old],
    ['T', 'F', false, false]);
  const both = lifted(request({ options: { featureTable: 'old', featureVisibilityTable: 'new' } })).diagramOptions;
  assert.deepEqual([both.featureVisibilityTable, 'featureTable' in both], ['new', false]);
  const unset = lifted(request({ options: { featureTable: 'old', featureVisibilityTable: null } })).diagramOptions;
  assert.equal(unset.featureVisibilityTable, 'old');
});

test('the sparse defaults are written out (A4, A5, A6)', () => {
  const options = lifted(request({ options: { config: null, configOverrides: { 'canvas.dpi': null, label_placement: 'auto' } } })).diagramOptions;
  assert.ok(Array.isArray(options.selectedFeaturesSet) && options.selectedFeaturesSet.length > 0);
  assert.equal(options.featureShapes.repeat_region, 'rectangle');
  assert.equal('canvas.dpi' in options.configOverrides, false, 'a null override is dropped');
  assert.ok(Object.keys(options.configOverrides).length > 1, 'the mode defaults are added');
  const kept = lifted(request({ options: { selectedFeaturesSet: ['CDS'], featureShapes: { repeat_region: 'arrow' } } })).diagramOptions;
  assert.deepEqual([kept.selectedFeaturesSet, kept.featureShapes.repeat_region], [['CDS'], 'arrow']);
  const full = lifted(request({ options: { config: { canvas: {} }, configOverrides: { 'canvas.dpi': null } } })).diagramOptions;
  assert.deepEqual(full.configOverrides, { 'canvas.dpi': null }, 'with a full config the overrides stay as stored');
});

test('retired override names and Linear tokens take their current forms (A7, OV-795)', () => {
  const linear = lifted(request({
    options: {
      configOverrides: {
        depth_tick_interval: 5, label_placement: ' On_Feature ', 'labels.linear.placement': 'on_feature',
        linear_track_layout: 'SpreadOut', 'canvas.linear.track_layout': 'tuckin'
      },
      config: { labels: { linear: { placement: 'on_feature' } }, canvas: { linear: { track_layout: 'tuckin' } } }
    }
  })).diagramOptions;
  assert.deepEqual(linear.configOverrides, {
    depth_large_tick_interval: 5, label_placement: 'above_feature', 'labels.linear.placement': 'above_feature',
    linear_track_layout: 'above', 'canvas.linear.track_layout': 'below'
  });
  assert.deepEqual([linear.config.labels.linear.placement, linear.config.canvas.linear.track_layout], ['above_feature', 'below']);
  const both = lifted(request({ options: { configOverrides: { depth_tick_interval: 5, depth_large_tick_interval: 7 } } }));
  assert.deepEqual(both.diagramOptions.configOverrides, { depth_large_tick_interval: 7 });
  const circular = lifted(request({ mode: 'circular', options: { configOverrides: { label_placement: 'on_feature' } } })).diagramOptions;
  assert.equal(circular.configOverrides.label_placement, 'on_feature', 'Circular keeps its tokens');
});

test('an invalid Linear comparison height is dropped with one repair row (B17)', () => {
  const result = liftLegacyRenderRequest(request({ options: { configOverrides: { comparison_height: 'tall', 'canvas.linear.comparison_height': 40 } } }));
  assert.deepEqual(result.repairs, [loadRepair({
    code: 'setting-replaced', path: 'renderRequest.diagramOptions.configOverrides.comparison_height', effect: 'dropped', stored: 'tall'
  })]);
  assert.deepEqual(result.value.diagramOptions.configOverrides, { 'canvas.linear.comparison_height': 40 });
  assert.deepEqual(liftLegacyRenderRequest(request({ options: { configOverrides: { comparison_height: '' } } })).repairs, []);
});

test('Circular slots lose the legacy spacing and tokens (A8, A9)', () => {
  const slots = lifted(request({
    mode: 'circular',
    options: {
      tracks: {
        circularTrackSlots: [
          { kind: 'circularTrackSlot', id: 'gc', renderer: 'dinucleotide_content', side: 'inside', spacing: { value: 3, unit: 'px' } },
          { kind: 'circularTrackSlot', id: 'odd', renderer: 'features', spacing: '3 parsecs' },
          'gc:dinucleotide_content@side=inside,r=0.5,w=0.1,spacing=2px,strict'
        ]
      }
    }
  })).diagramOptions.tracks.circularTrackSlots;
  assert.deepEqual([slots[0].innerGapPx, slots[0].outerGapPx, 'spacing' in slots[0], 'inner_gap_px' in slots[0]], [3, 3, false, false]);
  assert.equal(slots[1].spacing, '3 parsecs', 'a spacing that cannot be read is left for the reader to refuse');
  assert.equal(slots[2], 'gc:dinucleotide_content@side=inside,r=0.5,w=0.1,inner_gap_px=2px,outer_gap_px=2px');
});

test('a Circular sqrt size mode reads as auto (A3)', () => {
  assert.equal(lifted(request({ mode: 'circular', layout: { multiRecordSizeMode: ' SQRT' } })).layout.multiRecordSizeMode, 'auto');
  assert.equal(lifted(request({ mode: 'circular', layout: { multiRecordSizeMode: 'fixed' } })).layout.multiRecordSizeMode, 'fixed');
});

test('a schema-1 Linear layout is ignored, as the Web read it (OV-793)', () => {
  assert.deepEqual(lifted(request({ schema: 1, layout: { rows: [[0]] } })).layout, {});
  assert.deepEqual(lifted(request({ schema: 2, layout: { rows: [[0]] } })).layout, { rows: [[0]] });
});

test('schema-1 records get the record keys the Web gave them', () => {
  const records = lifted(request({ schema: 1, records: [{}, { recordKey: 'kept' }, {}] })).records;
  assert.deepEqual(records.map(({ recordKey }) => recordKey), ['record-1', 'kept', 'record-3']);
  assert.deepEqual(lifted(request({ schema: 2, records: [{}] })).records, [{}]);
});

test('comparisons get explicit record indices and lossless collinearity (A12, A13, A14, OV-794, OV-795)', () => {
  const standard = {
    minAnchors: 3, maxGeneGap: 4, blockMergeGap: 9, singletonMergeGap: 9, maxDiagonalDrift: 5, maxConflictsInMergeGap: 6,
    maxParalogLinksPerOrthogroup: 7, gapPenalty: 9, nearbyDuplicateWindow: 9, scoreMode: 'x', constantAnchorScore: 9,
    minBlockScore: 9, blockEvalue: 9
  };
  const comparisons = (schema, mode = 'collinear', links = 2) => lifted(request({
    schema,
    comparisons: [
      { kind: 'nucleotideBlast', queryRecordIndex: 'x' },
      { kind: 'nucleotideBlast', queryRecordIndex: 4, subjectRecordIndex: '5' },
      { kind: 'precomputedProteinComparison' },
      { kind: 'precomputedProteinComparison' },
      {
        kind: 'generatedProteinComparison', mode,
        settings: { collinearityParams: { kind: 'standard', parameters: standard }, collinearMaxParalogLinksPerOrthogroup: links }
      }
    ]
  })).comparisons;
  const [first, second, protein1, protein2, generated] = comparisons(1);
  assert.deepEqual([first.queryRecordIndex, first.subjectRecordIndex], [0, 1]);
  assert.deepEqual([second.queryRecordIndex, second.subjectRecordIndex], [4, 5]);
  assert.deepEqual([protein1.queryRecordIndex, protein1.subjectRecordIndex, protein2.queryRecordIndex, protein2.subjectRecordIndex], [0, 1, 1, 2]);
  assert.deepEqual(generated.pairs, []);
  assert.deepEqual(generated.settings.collinearityParams, {
    kind: 'lossless',
    parameters: { minAnchors: 3, maxUnitGap: 4, maxDiagonalDrift: 5, maxConflicts: 6, mergeOrientation: 'either' }
  });
  assert.equal(generated.settings.collinearMaxParalogLinksPerOrthogroup, 7, 'the standard value replaces the default');
  assert.equal(comparisons(1, 'collinear', 3)[4].settings.collinearMaxParalogLinksPerOrthogroup, 3, 'a set value is kept');
  assert.equal(comparisons(1, 'orthogroup')[4].settings.collinearMaxParalogLinksPerOrthogroup, 2);
  assert.equal('queryRecordIndex' in comparisons(2)[2], false, 'schema 2 stored the protein indices');
  const partial = lifted(request({
    comparisons: [{ kind: 'generatedProteinComparison', pairs: [{ q: 0 }], settings: { collinearityParams: { kind: 'standard', parameters: { minAnchors: 1 } } } }]
  })).comparisons[0];
  assert.deepEqual([partial.pairs, partial.settings.collinearityParams.kind], [[{ q: 0 }], 'standard']);
});

test('each request of the corpus lifts without repairs to a request promoted to the current schema', () => {
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
    assert.equal(promoteCanonicalRenderRequestToCurrent(value, /^BGC/.test(id) ? promotion : {}).schema, 9, id);
  }
});

test('a full config whose small tick interval Python reads as unset projects as null (OV-791)', () => {
  const { renderRequest } = corpus.requests.find(({ id }) => id === 'lambda_basic_linear@10d3a3d2');
  assert.equal(renderRequest.diagramOptions.config.objects.gc_content.small_tick_interval, 'none');
  /** @type {Record<string, object>} */
  const resources = {};
  JSON.stringify(renderRequest, (key, value) => {
    if (/resourceId$/i.test(key) && typeof value === 'string') {
      resources[value] = { kind: 'canonical-tsv', name: `${value}.txt`, type: 'text/plain', size: 0, lastModified: 0, encoding: 'base64', data: '' };
    }
    return value;
  });
  const { config } = projectCanonicalSessionRequest({ renderRequest, resources, deferResourceContent: true });
  assert.deepEqual([config.adv.gc_content_small_tick_interval, config.adv.depth_small_tick_interval], [null, null]);
});
