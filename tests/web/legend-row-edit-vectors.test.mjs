// The Legend row edits Python draws equal what the Web's executor makes of
// Python's un-edited Legend (vector:tests/fixtures/legend_row_edit_vectors.json,
// tools/twins.json). Each case runs Generate admission (`applyLegendOperations`
// and `nameLegendRows`) over `unedited` with the case's operations: the feature
// Legend rows it leaves equal `expected` (Python's `edited` without geometry),
// and a second pass over Python's `edited` changes nothing.
// tools/generate_legend_row_edit_vectors.py writes the cases;
// tests/test_legend_row_edit_vectors.py redraws them with Python.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { normalizeGenerationResponse } from '../../gbdraw/web/js/services/diagram-generation.js';
import { FEATURE_CATALOG_SCHEMA, admitFeatureCatalog } from '../../gbdraw/web/js/services/feature-catalog.js';
import { getAllFeatureLegendGroups } from '../../gbdraw/web/js/services/legend-svg.js';
import { admitCurrentGeneratedResults, createEmptySvgMutationPlan } from '../../gbdraw/web/js/services/svg-result-ingestion.js';

installFakeSvgDom();

const VECTORS = JSON.parse(readFileSync(new URL('../fixtures/legend_row_edit_vectors.json', import.meta.url), 'utf8'));
// The executor keeps the rows' places until the Web lays the Legend out; Python lays it out again.
const GEOMETRY = new Set(['transform']);

const decode = (value) => String(value)
  .replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&quot;/g, '"').replace(/&apos;/g, "'")
  .replace(/&#(\d+);/g, (_match, code) => String.fromCodePoint(Number(code)))
  .replace(/&amp;/g, '&');

const attributesOf = (element) => Object.fromEntries([...element.attributes.entries()]
  .filter(([name]) => !GEOMETRY.has(name))
  .map(([name, value]) => [name, decode(value)]));

const nodeOf = (element) => ({
  tag: element.localName,
  attributes: attributesOf(element),
  text: element.children.length > 0 ? null : element.textContent,
  children: element.children.map(nodeOf)
});

// The feature Legend groups as the vectors hold them: each group and its rows, in document order.
const legendGroupsOf = (svg) => getAllFeatureLegendGroups(svg).map((group) => {
  const direct = group.children.filter((child) => child.localName === 'g' && child.hasAttribute('data-legend-key'));
  const rows = direct.length > 0 ? direct : group.querySelectorAll('g[data-legend-key]');
  return { attributes: attributesOf(group), rows: rows.map(nodeOf) };
});

const catalog = () => ({
  schema: FEATURE_CATALOG_SCHEMA,
  items: [{
    resultIndex: 0,
    resultName: 'diagram.svg',
    recordKeys: [],
    features: [],
    biologicalFeatures: [],
    orthogroups: [],
    annotations: [],
    comparisonMatches: []
  }]
});

// Generate admission of one Result with the case's Legend operations. Returns
// the committed content, the feature Legend groups the executor left, and how
// many times admission serialized the Result (none when nothing changed).
const admit = (content, operations) => {
  const featureCatalog = catalog();
  const response = normalizeGenerationResponse({ results: [{ name: 'diagram.svg', content }], metadata: {}, featureCatalog });
  const catalogAdmission = admitFeatureCatalog(featureCatalog, response.results, { adopt: true, mode: 'linear' });
  let groups = null;
  const plan = createEmptySvgMutationPlan(1);
  const mutationPlan = {
    ...plan,
    kind: 'MUTATING',
    operationsByResult: [{
      ...plan.operationsByResult[0],
      ...operations,
      // Reads the executor's rows; returning false, it changes nothing itself.
      callerTransforms: [(svg) => { groups = legendGroupsOf(svg); return false; }]
    }]
  };
  const metrics = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
  try {
    const [result] = admitCurrentGeneratedResults(response, {
      catalogAdmission,
      mutationPlan,
      sanitizer: { sanitize: (value) => String(value) },
      parser: globalThis.DOMParser
    });
    const serialized = metrics.filter(({ name }) => name === 'svgSerializationCount').reduce((total, { value }) => total + value, 0);
    return { content: result.content, groups, serialized };
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
};

test('the vectors cover every edit kind in both modes, with a comparison block in each', () => {
  const ids = VECTORS.cases.map(({ id }) => id);
  for (const kind of ['rename-plain', 'rename-merges', 'rename-numbers', 'delete-one', 'delete-all-but-one',
    'order-reverse', 'order-sort-a-z', 'combined']) {
    assert.ok(ids.includes(`circular-${kind}`) && ids.includes(`linear-${kind}`), kind);
  }
  for (const id of ['circular-delete-merged', 'circular-order-move', 'circular-order-partial', 'circular-missing-key',
    'circular-comparison-ring', 'linear-pairwise-comparison', 'linear-combined-top', 'circular-delete-all',
    'linear-delete-all', 'circular-rule-rows', 'linear-rule-rows', 'linear-missing-key']) {
    assert.ok(ids.includes(id), id);
  }
});

for (const vector of VECTORS.cases) {
  test(`${vector.id}: the executor over Python's un-edited Legend draws Python's edited rows`, () => {
    const { groups } = admit(vector.unedited, vector.operations);
    assert.deepEqual(groups, vector.expected);
    const none = { legendRenames: [], legendDeletes: [], legendOrder: [] };
    assert.notDeepEqual(admit(vector.unedited, none).groups, vector.expected, 'the edits change the rows');
  });

  test(`${vector.id}: a second pass over Python's edited Legend changes nothing`, () => {
    const { content, groups, serialized } = admit(vector.edited, vector.operations);
    assert.equal(serialized, 0, 'admission keeps the bytes of a Result it leaves as it is');
    assert.equal(content, vector.edited);
    assert.deepEqual(groups, vector.expected);
  });
}
