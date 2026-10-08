import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';

import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
import {
  requireUniqueEditableLabelBindings
} from '../../gbdraw/web/js/app/feature-editor/label-actions.js';
import { normalizeGenerationResponse } from '../../gbdraw/web/js/services/diagram-generation.js';
import {
  admitFeatureCatalog,
  biologicalFeatureKey
} from '../../gbdraw/web/js/services/feature-catalog.js';
import { normalizeLogicalResults } from '../../gbdraw/web/js/services/result-normalization.js';
import {
  admitCurrentGeneratedResults,
  admitCurrentSessionResults,
  admitLegacyImportedResults,
  createCurrentSessionResultSource,
  createEmptySvgMutationPlan,
  createLegacyImportResultSource,
  createSavedResultPlan,
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  isCommittedSvgResult,
  markCommittedSvgResultMounted,
  markCommittedSvgResultUnmounted,
  reconcileMountedResult
} from '../../gbdraw/web/js/services/svg-result-ingestion.js';
import { stripResultBaseAttributes } from '../../gbdraw/web/js/services/result-paint-bases.js';
import { recordRuleMatches, ruleKey } from '../../gbdraw/web/js/services/rule-matchers.js';
import { displayedFeatureAddressing, featureOverrideKey } from '../../gbdraw/web/js/services/feature-override-identity.js';

class FakeElement {
  constructor(tagName, attributes = {}, children = []) {
    this.tagName = tagName;
    this.localName = tagName;
    this.attributes = new Map(Object.entries(attributes));
    this.children = [];
    this.parentElement = null;
    this.textContent = '';
    this.style = { removeProperty() {} };
    children.forEach((child) => this.appendChild(child));
  }

  get id() { return this.getAttribute('id') || ''; }
  getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }
  setAttribute(name, value) { this.attributes.set(name, String(value)); }
  removeAttribute(name) { this.attributes.delete(name); }
  hasAttribute(name) { return this.attributes.has(name); }
  appendChild(child) { child.parentElement = this; this.children.push(child); return child; }
  remove() {
    if (!this.parentElement) return;
    this.parentElement.children = this.parentElement.children.filter((child) => child !== this);
    this.parentElement = null;
  }
  cloneNode(deep = false) {
    const clone = new FakeElement(this.tagName, Object.fromEntries(this.attributes));
    clone.textContent = this.textContent;
    if (deep) this.children.forEach((child) => clone.appendChild(child.cloneNode(true)));
    return clone;
  }
  getElementById(id) {
    return this.walk().find((element) => element.id === id) || null;
  }
  walk() {
    return [this, ...this.children.flatMap((child) => child.walk())];
  }
  matches(selector) {
    if (/^\[[\w-]+\](?:,\s*\[[\w-]+\])*$/.test(selector)) {
      return selector.split(',').some((part) => this.hasAttribute(part.trim().slice(1, -1)));
    }
    if (selector.startsWith('.')) return false;
    if (selector === 'path') return this.tagName === 'path';
    if (selector === 'text') return this.tagName === 'text';
    if (selector === 'textPath') return this.tagName === 'textPath';
    if (selector === 'g[data-legend-key]') {
      return this.tagName === 'g' && this.hasAttribute('data-legend-key');
    }
    if (selector === 'text[data-label-feature-id]') {
      return this.tagName === 'text' && this.hasAttribute('data-label-feature-id');
    }
    if (selector.includes('data-gbdraw-feature-id') || selector.includes('id^="f"')) {
      return ['path', 'polygon', 'rect'].includes(this.tagName)
        && (this.hasAttribute('data-gbdraw-feature-id') || this.id.startsWith('f'));
    }
    if (selector.startsWith('#')) return this.id === selector.slice(1);
    return false;
  }
  querySelectorAll(selector) {
    return this.walk().slice(1).filter((element) => element.matches(selector));
  }
  querySelector(selector) { return this.querySelectorAll(selector)[0] || null; }
}

const serializeNode = (node) => {
  const attributes = Array.from(node.attributes.entries())
    .map(([name, value]) => ` ${name}="${value}"`).join('');
  const content = `${node.textContent || ''}${node.children.map(serializeNode).join('')}`;
  return `<${node.tagName}${attributes}>${content}</${node.tagName}>`;
};

globalThis.XMLSerializer = class {
  serializeToString(node) { return serializeNode(node); }
};

const buildSvgRoot = ({
  missingFeature = false,
  missingLabel = false,
  missingLegend = false
} = {}) => {
  const root = new FakeElement('svg', { xmlns: 'http://www.w3.org/2000/svg' });
  if (!missingFeature) {
    root.appendChild(new FakeElement('path', {
      id: 'f0001',
      'data-gbdraw-feature-id': 'f0001',
      'data-gbdraw-feature-part': 'block',
      'data-gbdraw-stable-feature-id': 'stable-a',
      'data-gbdraw-record-index': '0',
      'data-gbdraw-record-id': 'record-a',
      fill: '#aaaaaa'
    }));
  }
  if (!missingLabel) {
    const label = new FakeElement('text', { 'data-label-feature-id': 'f0001' });
    label.textContent = 'old label';
    root.appendChild(label);
  }
  if (!missingLegend) {
    const legend = new FakeElement('g', { id: 'legend' });
    const featureLegend = new FakeElement('g', { id: 'feature_legend' });
    const entry = new FakeElement('g', { 'data-legend-key': 'CDS' });
    entry.appendChild(new FakeElement('path', { fill: '#aaaaaa' }));
    const caption = new FakeElement('text');
    caption.textContent = 'CDS';
    entry.appendChild(caption);
    featureLegend.appendChild(entry);
    legend.appendChild(featureLegend);
    root.appendChild(legend);
  }
  return root;
};

class FakeDomParser {
  static calls = 0;

  parseFromString(content, mediaType) {
    FakeDomParser.calls += 1;
    assert.equal(mediaType, 'image/svg+xml');
    const isSvg = content.includes('<svg');
    return {
      documentElement: isSvg
        ? buildSvgRoot({
            missingFeature: content.includes('missing-feature'),
            missingLabel: content.includes('missing-label'),
            missingLegend: content.includes('missing-legend')
          })
        : new FakeElement('html'),
      querySelector: () => null
    };
  }
}

const catalog = () => ({
  schema: 5,
  items: [{
    resultIndex: 0,
    resultName: 'diagram.svg',
    recordKeys: ['record-a'],
    features: [{
      svgId: 'f0001',
      recordKey: 'record-a',
      biologicalFeatureId: 'feature-a',
      fillColor: '#aaaaaa',
      drawnSelector: { hash: 'f0001', location: '0..3', recordLocation: 'record-a:0..3:+' }
    }],
    biologicalFeatures: [{
      recordKey: 'record-a',
      biologicalFeatureId: 'feature-a',
      stableFeatureId: 'stable-a',
      record_id: 'record-a',
      type: 'CDS',
      start: 0,
      end: 3,
      strand: 1,
      anchorProfile: {
        precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+'
      }
    }],
    orthogroups: [],
    annotations: [],
    comparisonMatches: [],
    sequenceSources: []
  }]
});

const sanitizer = (counter) => ({
  sanitize(content) {
    counter.calls += 1;
    return String(content).replace(/<script>[\s\S]*?<\/script>/gi, '');
  }
});

const currentFixture = (content = '<svg><path id="f0001"/></svg>') => {
  const featureCatalog = catalog();
  const response = normalizeGenerationResponse({
    results: [{ name: 'diagram.svg', content }],
    metadata: { featureCatalog }
  });
  const admission = admitFeatureCatalog(featureCatalog, response.results, {
    adopt: true,
    mode: 'linear'
  });
  return { response, admission };
};

const metricTotal = (metrics, name) => metrics
  .filter((metric) => metric.name === name)
  .reduce((total, metric) => total + metric.value, 0);

const captureMetrics = (run) => {
  const metrics = [];
  const events = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => metrics.push(metric),
    onSessionLifecycleEvent: (event) => events.push(event)
  };
  try {
    return { value: run(), metrics, events };
  } finally {
    delete globalThis.__GBDRAW_TEST_HOOKS__;
  }
};

test('current-worker EMPTY admission sanitizes once and performs zero application SVG work', () => {
  FakeDomParser.calls = 0;
  const calls = { calls: 0 };
  const { response, admission } = currentFixture(
    '<svg><script>bad()</script><path id="f0001"/></svg>'
  );
  const { value: results, metrics, events } = captureMetrics(() => admitCurrentGeneratedResults(
    response,
    {
      catalogAdmission: admission,
      mutationPlan: createEmptySvgMutationPlan(1),
      sanitizer: sanitizer(calls),
      parser: FakeDomParser,
      selectedFeatureTypes: ['CDS', 'tRNA']
    }
  ));
  const committed = results[0];
  // R13: a Result carries the feature types of the request that drew it.
  assert.deepEqual(getCommittedSvgResultMetadata(committed).selectedFeatureTypes, ['CDS', 'tRNA']);
  assert.equal(Object.isFrozen(getCommittedSvgResultMetadata(committed).selectedFeatureTypes), true);

  assert.equal(calls.calls, 1);
  assert.equal(FakeDomParser.calls, 0);
  assert.equal(committed.content, '<svg><path id="f0001"/></svg>');
  assert.equal(isCommittedSvgResult(committed), true);
  assert.equal(
    getCommittedSvgResultMetadata(committed).renderedFeatureIdentities
      .byRenderedId.get('f0001').stableId,
    'stable-a'
  );
  assert.equal(metricTotal(metrics, 'svgSanitizationCount'), 1);
  for (const name of [
    'applicationSvgParseCount',
    'svgMutationIndexBuildCount',
    'featureDomFullScanCount',
    'legendDomFullScanCount',
    'svgIdentityScanCount',
    'svgSerializationCount',
    'currentLegacyNormalizationCount'
  ]) assert.equal(metricTotal(metrics, name), 0, name);
  const completedCandidate = events.find(({ name }) => name === 'artifact.candidate-completed');
  assert.deepEqual(completedCandidate.catalogFootprint, admission.scalarMetrics);

  assert.equal(markCommittedSvgResultMounted(committed), true);
  const persistedEdit = { ...committed, content: '<svg><path fill="#abcdef"/></svg>' };
  assert.equal(getCommittedSvgResultMetadata(persistedEdit), getCommittedSvgResultMetadata(committed));
  assert.equal(getCommittedSvgResultMetadata(normalizeLogicalResults([committed])[0]), getCommittedSvgResultMetadata(committed));
  assert.equal(getCommittedSvgContent(persistedEdit), committed.content);
  assert.equal(markCommittedSvgResultUnmounted(persistedEdit), true);
  assert.equal(getCommittedSvgContent(persistedEdit), persistedEdit.content);
  const persistedRoundTrip = JSON.parse(JSON.stringify(committed));
  assert.equal(isCommittedSvgResult(persistedRoundTrip), false);
  assert.equal(getCommittedSvgResultMetadata(persistedRoundTrip), null);
});

test('JSON cannot forge current-worker provenance and malformed envelopes fail closed', () => {
  const { response, admission } = currentFixture();
  assert.throws(
    () => admitCurrentGeneratedResults(JSON.parse(JSON.stringify(response)), {
      catalogAdmission: admission,
      mutationPlan: createEmptySvgMutationPlan(1),
      sanitizer: { sanitize: (value) => value }
    }),
    /runtime Worker provenance/
  );

  const malformed = currentFixture('not svg');
  assert.throws(
    () => admitCurrentGeneratedResults(malformed.response, {
      catalogAdmission: malformed.admission,
      mutationPlan: createEmptySvgMutationPlan(1),
      sanitizer: { sanitize: (value) => value }
    }),
    /malformed SVG/
  );

  const unrelated = currentFixture();
  assert.throws(
    () => admitCurrentGeneratedResults(response, {
      catalogAdmission: unrelated.admission,
      mutationPlan: createEmptySvgMutationPlan(1),
      sanitizer: { sanitize: (value) => value }
    }),
    /do not own the admitted feature catalog/
  );
});

// Per-feature edits are identity rows of the Result's mode (design Q4, R2).
// A row of the Linear drawing (PR-1: the drawing is the mode).
const identityKey = JSON.stringify(['record-a', 'feature-a']);
const identityRow = (fields) => ({
  recordKey: 'record-a',
  biologicalFeatureId: 'feature-a',
  featureVisibility: null,
  labelVisibility: null,
  labelText: null,
  labelSourceText: null,
  ...fields
});

const planOptions = {
  fill: (admission) => ({
    catalogAdmission: admission,
    featureColorOverrides: {
      [biologicalFeatureKey('record-a', 'feature-a')]: '#112233'
    }
  }),
  stroke: (admission) => ({
    catalogAdmission: admission,
    featureStrokeOverrides: {
      [biologicalFeatureKey('record-a', 'feature-a')]: {
        strokeColor: '#223344',
        strokeWidth: 2
      }
    }
  }),
  visibility: (admission) => ({
    catalogAdmission: admission,
    featureOverrides: { [identityKey]: identityRow({ featureVisibility: 'off' }) }
  }),
  Label: (admission) => ({
    catalogAdmission: admission,
    featureOverrides: { [identityKey]: identityRow({ labelText: 'new label', labelVisibility: 'off' }) }
  }),
  Legend: (admission) => ({
    catalogAdmission: admission,
    legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#334455' }],
    originalLegendOrder: ['CDS'],
    legendColorOverrides: { CDS: '#334455' },
    legendStrokeOverrides: { CDS: { strokeColor: '#445566', strokeWidth: 3 } }
  })
};

for (const [domain, makeOptions] of Object.entries(planOptions)) {
  test(`${domain}-only current mutation uses one admitted root and shared lazy index`, () => {
    FakeDomParser.calls = 0;
    const { response, admission } = currentFixture();
    const plan = compileDirectEditorMutationPlan(makeOptions(admission));
    assert.equal(plan.kind, 'MUTATING');
    const { metrics } = captureMetrics(() => admitCurrentGeneratedResults(response, {
      catalogAdmission: admission,
      mutationPlan: plan,
      sanitizer: { sanitize: (value) => value },
      parser: FakeDomParser
    }));
    const detachedSvgWork = domain === 'Label' ? 0 : 1;
    assert.equal(FakeDomParser.calls, detachedSvgWork);
    assert.equal(metricTotal(metrics, 'applicationSvgParseCount'), detachedSvgWork);
    assert.equal(metricTotal(metrics, 'svgMutationIndexBuildCount'), detachedSvgWork);
    assert.equal(metricTotal(metrics, 'svgSerializationCount'), detachedSvgWork);
    assert.equal(metricTotal(metrics, 'svgIdentityScanCount'), 0);
    assert.equal(metricTotal(metrics, 'currentLegacyNormalizationCount'), 0);
  });
}

// R2 (R-1): the Circular drawing's row of the same record key and feature
// waits for a Circular Result; the plan of this Linear Result reads the Linear
// drawing only (PR-1), so it draws none of it live.
test('a per-feature edit of the other mode does not mutate the Result', () => {
  const { admission } = currentFixture();
  const drawings = {
    circular: { featureOverrides: { [identityKey]: identityRow({ featureVisibility: 'off', labelText: 'other mode' }) } },
    linear: { featureOverrides: {} }
  };
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    featureOverrides: drawings.linear.featureOverrides
  });
  assert.equal(plan.kind, 'EMPTY');
  assert.equal(compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    featureOverrides: drawings.circular.featureOverrides
  }).kind, 'MUTATING', 'the same row in the drawing of the Result mutates it');
});

test('combined current mutations share one root/index and serialize once', () => {
  FakeDomParser.calls = 0;
  const { response, admission } = currentFixture();
  const plan = compileDirectEditorMutationPlan({
    ...planOptions.fill(admission),
    ...planOptions.stroke(admission),
    ...planOptions.Legend(admission),
    // One identity row carries the feature and label edits.
    featureOverrides: {
      [identityKey]: identityRow({ featureVisibility: 'off', labelText: 'new label', labelVisibility: 'off' })
    }
  });
  const { value: results, metrics } = captureMetrics(() => admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: plan,
    sanitizer: { sanitize: (value) => value },
    parser: FakeDomParser
  }));
  assert.equal(FakeDomParser.calls, 1);
  assert.equal(metricTotal(metrics, 'applicationSvgParseCount'), 1);
  assert.equal(metricTotal(metrics, 'svgMutationIndexBuildCount'), 1);
  assert.equal(metricTotal(metrics, 'featureDomFullScanCount'), 1);
  assert.equal(metricTotal(metrics, 'legendDomFullScanCount'), 1);
  assert.equal(metricTotal(metrics, 'svgSerializationCount'), 1);
  assert.match(results[0].content, /fill="#112233"/);
  assert.match(results[0].content, /display="none"/);
  assert.doesNotMatch(results[0].content, /data-gbdraw-label-visibility-preview="off"/);
  assert.match(results[0].content, />old label<\/text>/);
  assert.match(results[0].content, /fill="#334455"/);
});

test('direct Legend rename, deletion, and addition are applied through catalog-backed admission', () => {
  const cases = [
    {
      editor: {
        legendEntries: [{
          caption: 'Genes', originalCaption: 'CDS', color: '#aaaaaa', xPos: 20, yPos: 30
        }],
        originalLegendOrder: ['CDS']
      },
      // A rename keeps the row where the renderer drew it; the Legend layout
      // places it (OV-156).
      assertContent: (content) => {
        assert.match(content, /data-legend-key="Genes"/);
        assert.doesNotMatch(content, /transform=/);
      }
    },
    {
      editor: {
        legendEntries: [],
        deletedLegendEntries: [{ caption: 'CDS', originalCaption: 'CDS' }],
        originalLegendOrder: ['CDS']
      },
      assertContent: (content) => assert.doesNotMatch(content, /data-legend-key="CDS"/)
    },
    {
      editor: {
        legendEntries: [{
          caption: 'New', originalCaption: 'New', color: '#556677', xPos: 40, yPos: 50
        }],
        originalLegendOrder: ['CDS']
      },
      assertContent: (content) => {
        assert.match(content, /data-legend-key="New"/);
        assert.match(content, /data-legend-owner="direct-editor"/);
        assert.match(content, /transform="translate\(40, 50\)"/);
      }
    }
  ];

  cases.forEach(({ editor, assertContent }) => {
    const { response, admission } = currentFixture();
    const plan = compileDirectEditorMutationPlan({ catalogAdmission: admission, ...editor });
    const results = admitCurrentGeneratedResults(response, {
      catalogAdmission: admission,
      mutationPlan: plan,
      sanitizer: { sanitize: (value) => value },
      parser: FakeDomParser
    });
    assertContent(results[0].content);
  });
});

// OV-86: Python never draws a row the Legend editor added. Admission adds the row
// before it styles rows, so the row's own fill and stroke reach it, and an added row
// takes the renderer's row shape without another row's style.
test('a fill and stroke on a Legend editor added row style the row the plan adds', () => {
  const { response, admission } = currentFixture();
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries: [
      { caption: 'CDS', originalCaption: 'CDS', color: '#aaaaaa' },
      { caption: 'Manual row', originalCaption: 'Manual row', color: '#7b2cbf' },
      { caption: 'Plain row', originalCaption: 'Plain row', color: '#118833' }
    ],
    originalLegendOrder: ['CDS'],
    legendColorOverrides: { 'Manual row': '#7b2cbf' },
    legendStrokeOverrides: {
      'Manual row': { strokeColor: '#e63946', strokeWidth: 2 },
      CDS: { strokeColor: '#445566', strokeWidth: 3 }
    }
  });
  // The style still requires its row; the same plan adds it.
  assert.deepEqual(plan.operationsByResult[0].legendFills.map(({ caption, allowMissing }) => [caption, allowMissing]), [
    ['Manual row', false]
  ]);
  const swatches = (content) => Object.fromEntries(
    [...withoutBases(content).matchAll(/<g data-legend-key="([^"]+)"[^>]*><path ([^>]*)>/g)].map(([, caption, path]) => [caption, path])
  );
  const expected = {
    CDS: 'fill="#aaaaaa" stroke="#445566" stroke-width="3"',
    'Manual row': 'fill="#7b2cbf" stroke="#e63946" stroke-width="2"',
    'Plain row': 'fill="#118833"'
  };
  const [result] = admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: plan,
    sanitizer: { sanitize: (value) => value },
    parser: FakeDomParser
  });
  assert.deepEqual(swatches(result.content), expected);
  // A displayed Result receives the same operations through the same executor (D-07).
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, plan.operationsByResult[0]);
  assert.deepEqual(swatches(serializeNode(mounted)), expected);
});

const withoutBases = (content) => content.replace(/ data-gbdraw-base-[\w-]+="[^"]*"/g, '');

// D-07 (PD-OI-062): a batch Result shows the editor operations through the
// executor when it is displayed. Showing them again changes nothing, and an
// operation for a feature or Legend row the Result does not draw is skipped.
const mountedPaint = (admission) => compileDirectEditorMutationPlan({
  ...planOptions.fill(admission),
  ...planOptions.stroke(admission),
  ...planOptions.Legend(admission),
  featureOverrides: { [identityKey]: identityRow({ featureVisibility: 'off' }) }
}).operationsByResult[0];

test('the mounted executor is idempotent', () => {
  const { admission } = currentFixture();
  const operations = mountedPaint(admission);
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, operations);
  const once = serializeNode(mounted);
  assert.notEqual(once, serializeNode(buildSvgRoot()));
  reconcileMountedResult(mounted, operations);
  assert.equal(serializeNode(mounted), once);
});

test('a mounted batch Result skips the features and Legend rows it does not draw', () => {
  const { admission } = currentFixture();
  const renamed = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries: [{ caption: 'Genes', originalCaption: 'CDS', color: '#aaaaaa' }],
    originalLegendOrder: ['CDS']
  }).operationsByResult[0];
  [mountedPaint(admission), renamed].forEach((operations) => {
    const mounted = buildSvgRoot({ missingFeature: true });
    mounted.querySelector('g[data-legend-key]').setAttribute('data-legend-key', 'tRNA');
    const drawn = serializeNode(mounted);
    reconcileMountedResult(mounted, operations);
    assert.equal(serializeNode(mounted), drawn);
  });
});

// OV-144: the executor records Python's value of each paint attribute it
// changes, so a reconcile returns what no operation names to Python's bytes.
test('a mounted Result shown again after its edits were removed returns to what Python drew', () => {
  const { admission } = currentFixture();
  const mounted = buildSvgRoot();
  const python = serializeNode(mounted);
  reconcileMountedResult(mounted, mountedPaint(admission));
  assert.match(serializeNode(mounted), /data-gbdraw-base-stroke="" stroke="#223344"/);
  reconcileMountedResult(mounted, createEmptySvgMutationPlan(1).operationsByResult[0]);
  assert.equal(serializeNode(mounted), python);
});

test('a reconcile returns only the listed paint domains to Python and applies every operation', () => {
  const { admission } = currentFixture();
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, mountedPaint(admission));
  const painted = serializeNode(mounted);
  const strokesOnly = compileDirectEditorMutationPlan({
    ...planOptions.stroke(admission),
    legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#334455' }],
    originalLegendOrder: ['CDS']
  }).operationsByResult[0];
  // Fills, visibility, and the Legend fill keep the paint no operation names.
  reconcileMountedResult(mounted, strokesOnly, { domains: ['featureStrokes', 'legendStrokes'] });
  const content = serializeNode(mounted);
  assert.match(content, /fill="#112233"/);
  assert.match(content, /stroke="#223344"/);
  assert.match(content, /display="none"/);
  assert.match(content, /<path fill="#334455" data-gbdraw-base-fill="#aaaaaa"><\/path>/);
  // The Legend row stroke is gone; the feature's own stroke is applied again.
  assert.doesNotMatch(content, /#445566/);
  reconcileMountedResult(mounted, mountedPaint(admission), { domains: [] });
  assert.equal(serializeNode(mounted), painted, 'no listed domain: the operations apply as they are');
});

// Generate admission records Python's values too, so the first display of a
// Result Generate drew with edits can revert them. Exports strip the records,
// which leaves the bytes the executor wrote before it kept them.
test('admission records the paint Python drew once; an export strips the records', () => {
  const { response, admission } = currentFixture();
  const [result] = admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: compileDirectEditorMutationPlan({
      ...planOptions.fill(admission),
      ...planOptions.stroke(admission),
      ...planOptions.Legend(admission)
    }),
    sanitizer: { sanitize: (value) => value },
    parser: FakeDomParser
  });
  assert.match(result.content, /data-gbdraw-base-fill="#aaaaaa" fill="#112233"|fill="#112233"[^>]*data-gbdraw-base-fill="#aaaaaa"/);
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, mountedPaint(admission));
  reconcileMountedResult(mounted, mountedPaint(admission));
  assert.equal((serializeNode(mounted).match(/data-gbdraw-base-fill="#aaaaaa"/g) || []).length, 2,
    'a second change keeps the value Python drew');
  stripResultBaseAttributes(mounted);
  assert.equal(serializeNode(mounted), withoutBases(serializeNode(mounted)));
  assert.doesNotMatch(serializeNode(mounted), /data-gbdraw-base-/);
  assert.match(serializeNode(mounted), /<path fill="#334455" stroke="#445566" stroke-width="3"><\/path>/);
});

// A row the Legend editor adds copies Python's first row as drawn, without the
// records of that row's edits; its own color is its drawn fill.
test('an added Legend row carries no record of the row it copies', () => {
  const { admission } = currentFixture();
  const compile = (legendEntries) => compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries,
    originalLegendOrder: ['CDS'],
    legendStrokeOverrides: { CDS: { strokeColor: '#445566', strokeWidth: 3 } }
  }).operationsByResult[0];
  const cds = { caption: 'CDS', originalCaption: 'CDS', color: '#aaaaaa' };
  const operations = compile([cds, { caption: 'Manual row', originalCaption: 'Manual row', color: '#7b2cbf' }]);
  // The first row shows its stroke edit before the row is added.
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, compile([cds]));
  reconcileMountedResult(mounted, operations);
  reconcileMountedResult(mounted, operations);
  const rows = Object.fromEntries([...serializeNode(mounted).matchAll(/<g data-legend-key="([^"]+)"[^>]*><path ([^>]*)>/g)]
    .map(([, caption, path]) => [caption, path]));
  assert.equal(rows['Manual row'], 'fill="#7b2cbf"');
  assert.match(rows.CDS, /data-gbdraw-base-stroke="" stroke="#445566"/);
});

test('a requested Legend addition reuses an exact renderer-produced caption', () => {
  const { response, admission } = currentFixture();
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries: [{
      caption: 'CDS',
      originalCaption: 'CDS',
      color: '#556677',
      xPos: 40,
      yPos: 50
    }],
    originalLegendOrder: ['other proteins']
  });
  const results = admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: plan,
    sanitizer: { sanitize: (value) => value },
    parser: FakeDomParser
  });
  assert.equal((results[0].content.match(/data-legend-key="CDS"/g) || []).length, 1);
  assert.match(results[0].content, /fill="#556677"/);
  assert.match(results[0].content, /data-legend-owner="direct-editor"/);
});

test('a renderer-derived Legend style may be absent when no current binding remains', () => {
  const { response, admission } = currentFixture('<svg>missing-legend</svg>');
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    featureColorOverrides: {
      [biologicalFeatureKey('record-a', 'retired-feature')]: {
        color: '#334455',
        caption: 'CDS'
      }
    },
    legendEntries: [{
      caption: 'CDS',
      originalCaption: 'CDS',
      color: '#334455'
    }],
    originalLegendOrder: ['CDS'],
    addedLegendCaptions: new Set(['CDS']),
    legendColorOverrides: { CDS: '#334455' },
    legendStrokeOverrides: { CDS: { strokeColor: '#445566', strokeWidth: 3 } }
  });
  assert.equal(plan.operationsByResult[0].legendFills[0].allowMissing, true);
  assert.equal(plan.operationsByResult[0].legendStrokes[0].allowMissing, true);
  assert.doesNotThrow(() => admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: plan,
    sanitizer: { sanitize: (value) => value },
    parser: FakeDomParser
  }));
});

// R6: the admission owner classifies the failure; Generate never reports UNKNOWN.
test('an unexplained missing Legend binding remains rejected', () => {
  const { response, admission } = currentFixture('<svg>missing-legend</svg>');
  const plan = compileDirectEditorMutationPlan(planOptions.Legend(admission));
  assert.throws(
    () => admitCurrentGeneratedResults(response, {
      catalogAdmission: admission,
      mutationPlan: plan,
      sanitizer: { sanitize: (value) => value },
      parser: FakeDomParser
    }),
    { code: 'RESULT_INVALID', stage: 'result-admission' }
  );
});

// OV-46, OV-63: the plan requires a generated Legend row (fill, stroke, rename) in
// each Result. Python's per-Result Legend row facts (`metadata.legendRows`) say which
// absences are legitimate: a row the draft removed (suppressed), or a row another
// Result of the batch draws or has removed. Any other absence is a stale operation.
const legendRowsFixture = ({ contents, legendRows, withFacts = true, rename = false }) => {
  const featureCatalog = catalog();
  const [first] = featureCatalog.items;
  if (contents.length > 1) {
    featureCatalog.items.push({
      ...first,
      resultIndex: 1,
      resultName: 'diagram-2.svg',
      recordKeys: ['record-b'],
      features: first.features.map((feature) => ({ ...feature, recordKey: 'record-b' })),
      biologicalFeatures: first.biologicalFeatures.map((feature) => ({
        ...feature, recordKey: 'record-b', record_id: 'record-b', stableFeatureId: 'stable-b'
      }))
    });
  }
  const names = contents.map((_, index) => (index === 0 ? 'diagram.svg' : 'diagram-2.svg'));
  const response = normalizeGenerationResponse({
    results: contents.map((content, index) => ({ name: names[index], content })),
    metadata: {
      featureCatalog,
      ...(withFacts ? {
        legendRows: legendRows.map(({ drawn = [], suppressed = [] }, index) => ({
          resultIndex: index, resultName: names[index], drawn, suppressed
        }))
      } : {})
    }
  });
  const admission = admitFeatureCatalog(featureCatalog, response.results, { adopt: true, mode: 'linear' });
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries: [{ caption: rename ? 'Genes' : 'CDS', originalCaption: 'CDS', color: '#00aa00' }],
    originalLegendOrder: ['CDS'],
    ...(rename ? {} : {
      legendColorOverrides: { CDS: '#00aa00' },
      legendStrokeOverrides: { CDS: { strokeColor: '#445566', strokeWidth: 3 } }
    })
  });
  const admit = () => admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: plan,
    sanitizer: { sanitize: (value) => value },
    parser: FakeDomParser
  });
  return { plan, admit };
};

const WITH_ROW = '<svg><path id="f0001"/></svg>';
const NO_ROW = '<svg>missing-legend</svg>';

test('a row Python reports as suppressed by the draft may be missing from its Result', () => {
  const { plan, admit } = legendRowsFixture({ contents: [NO_ROW], legendRows: [{ suppressed: ['CDS'] }] });
  assert.equal(plan.operationsByResult[0].legendFills[0].allowMissing, false);
  assert.doesNotThrow(admit);
  const rename = legendRowsFixture({ contents: [NO_ROW], legendRows: [{ suppressed: ['CDS'] }], rename: true });
  assert.doesNotThrow(rename.admit);
});

test('a row neither drawn nor suppressed is a stale operation and is rejected', () => {
  for (const rename of [false, true]) {
    const { admit } = legendRowsFixture({
      contents: [NO_ROW], legendRows: [{ drawn: ['GC content'], suppressed: ['repeat_region'] }], rename
    });
    assert.throws(admit, { code: 'RESULT_INVALID', stage: 'result-admission' });
  }
});

test('a row Python says it drew but the SVG lacks is rejected', () => {
  const { admit } = legendRowsFixture({ contents: [NO_ROW], legendRows: [{ drawn: ['CDS'] }] });
  assert.throws(admit, { code: 'RESULT_INVALID', stage: 'result-admission' });
});

test('without Legend row facts a required row stays required', () => {
  const { admit } = legendRowsFixture({ contents: [NO_ROW], legendRows: [], withFacts: false });
  assert.throws(admit, { code: 'RESULT_INVALID', stage: 'result-admission' });
});

test('a row another Result draws may be missing from a Result that does not name it', () => {
  const { admit } = legendRowsFixture({
    contents: [WITH_ROW, NO_ROW], legendRows: [{ drawn: ['CDS'] }, {}]
  });
  const results = admit();
  assert.match(results[0].content, /fill="#00aa00"/);
  assert.doesNotMatch(results[1].content, /#00aa00/);
});

test('a row only another Result suppresses may be missing; a row no Result reports is rejected', () => {
  assert.doesNotThrow(legendRowsFixture({
    contents: [NO_ROW, NO_ROW], legendRows: [{ suppressed: ['CDS'] }, {}]
  }).admit);
  const { admit } = legendRowsFixture({ contents: [NO_ROW, NO_ROW], legendRows: [{}, {}] });
  assert.throws(admit, { code: 'RESULT_INVALID', stage: 'result-admission' });
});

test('malformed Legend row facts are rejected at admission', () => {
  const base = { contents: [NO_ROW] };
  for (const legendRows of [[], [{ drawn: 'CDS' }], [{ suppressed: [1] }]]) {
    const { admit } = legendRowsFixture({ ...base, legendRows });
    assert.throws(admit, { code: 'RESULT_INVALID', stage: 'result-admission' });
  }
});

test('source replacement may retire styled, renamed, or deleted generated Legend categories', () => {
  const { response, admission } = currentFixture('<svg>missing-legend</svg>');
  const plan = compileDirectEditorMutationPlan({
    ...planOptions.Legend(admission),
    sourceReplaced: true,
    legendEntries: [{ caption: 'Genes', originalCaption: 'CDS', color: '#334455' }],
    originalLegendOrder: ['CDS', 'rRNA'],
    deletedLegendEntries: [{ caption: 'rRNA', originalCaption: 'rRNA' }],
    legendColorOverrides: { Genes: '#334455' },
    legendStrokeOverrides: { Genes: { strokeColor: '#445566', strokeWidth: 3 } }
  });
  assert.doesNotThrow(() => admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: plan,
    sanitizer: { sanitize: value => value },
    parser: FakeDomParser
  }));
});

test('explicit category deletion stays idempotent while the category is absent', () => {
  const { response, admission } = currentFixture('<svg>missing-legend</svg>');
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    originalLegendOrder: ['rRNA'],
    deletedLegendEntries: [{ caption: 'rRNA', originalCaption: 'rRNA' }]
  });
  assert.doesNotThrow(() => admitCurrentGeneratedResults(response, {
    catalogAdmission: admission, mutationPlan: plan,
    sanitizer: { sanitize: value => value }, parser: FakeDomParser
  }));
});

test('a dormant category style applies on its first returning Result without creating a manual row', () => {
  const { response, admission } = currentFixture();
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries: [], originalLegendOrder: [],
    legendColorOverrides: { CDS: '#224466' }
  });
  assert.equal(plan.operationsByResult[0].legendAdds.length, 0);
  const results = admitCurrentGeneratedResults(response, {
    catalogAdmission: admission, mutationPlan: plan,
    sanitizer: { sanitize: value => value }, parser: FakeDomParser
  });
  assert.match(results[0].content, /fill="#224466"/);
});

test('category absence capabilities never admit ambiguous Legend bindings', () => {
  class DuplicateLegendParser extends FakeDomParser {
    parseFromString(...args) {
      const document = super.parseFromString(...args);
      const group = document.documentElement.querySelector('#feature_legend');
      group.appendChild(group.children[0].cloneNode(true));
      return document;
    }
  }
  for (const sourceReplaced of [false, true]) {
    const { response, admission } = currentFixture();
    const plan = compileDirectEditorMutationPlan({
      ...planOptions.Legend(admission), sourceReplaced
    });
    assert.throws(() => admitCurrentGeneratedResults(response, {
      catalogAdmission: admission, mutationPlan: plan,
      sanitizer: { sanitize: value => value }, parser: DuplicateLegendParser
    }), /ambiguous Legend binding/);
  }
});

test('Label replay is deferred until mounted identity binding', () => {
  const { response, admission } = currentFixture('<svg>missing-label</svg>');
  const plan = compileDirectEditorMutationPlan(planOptions.Label(admission));
  const admitted = admitCurrentGeneratedResults(response, {
    catalogAdmission: admission,
    mutationPlan: plan,
    sanitizer: { sanitize: (value) => value },
    parser: FakeDomParser
  });
  assert.equal(isCommittedSvgResult(admitted[0]), true);
  assert.equal(isCommittedSvgResult(response.results[0]), false);
});

test('mounted Label binding accepts one target and rejects missing or ambiguous identity', () => {
  const label = (featureId) => new FakeElement('text', {
    'data-label-feature-id': featureId
  });
  assert.doesNotThrow(() => requireUniqueEditableLabelBindings(
    [label('f0001')],
    ['f0001']
  ));
  assert.throws(
    () => requireUniqueEditableLabelBindings([], ['f0001']),
    { code: 'LABEL_NOT_DRAWN', stage: 'render', context: { reason: 'FORCED_LABEL', featureId: 'f0001' } }
  );
  assert.throws(
    () => requireUniqueEditableLabelBindings(
      [label('f0001'), label('f0001')],
      ['f0001']
    ),
    { code: 'RENDER_FAILED', stage: 'render', context: { featureId: 'f0001' } }
  );
  assert.doesNotThrow(() => requireUniqueEditableLabelBindings(
    [],
    ['f0001'],
    { allowMissing: true }
  ));
});

test('current-session and legacy-import remain explicit, incompatible boundaries', () => {
  FakeDomParser.calls = 0;
  const calls = { calls: 0 };
  const persistedResults = [{
    name: 'diagram.svg',
    content: '<svg><script>bad()</script><path id="f0001"/></svg>'
  }];
  const admission = admitFeatureCatalog(catalog(), persistedResults, { mode: 'linear' });
  const currentSource = createCurrentSessionResultSource(persistedResults, admission);
  const current = admitCurrentSessionResults(currentSource, {
    mutationPlan: createEmptySvgMutationPlan(1),
    sanitizer: sanitizer(calls),
    parser: FakeDomParser
  });
  assert.equal(FakeDomParser.calls, 0);
  assert.equal(current[0].content.includes('<script>'), false);
  assert.equal(getCommittedSvgResultMetadata(current[0]).selectedFeatureTypes, null);
  assert.throws(
    () => admitLegacyImportedResults(currentSource),
    /legacy-import/
  );

  const legacy = admitLegacyImportedResults(
    createLegacyImportResultSource([{ name: 'legacy.svg', content: '<svg></svg>' }]),
    { sanitizer: { sanitize: (value) => value }, parser: FakeDomParser, selectedFeatureTypes: ['CDS'] }
  );
  assert.equal(FakeDomParser.calls, 1);
  assert.equal(getCommittedSvgResultMetadata(legacy[0]).sourceClass, 'legacy-import');
  assert.deepEqual(getCommittedSvgResultMetadata(legacy[0]).selectedFeatureTypes, ['CDS']);
});

// The elements and attributes of an SVG text as FakeElements (the saved
// Result of a fixture Session).
const parseSvgText = (content) => {
  const decode = (value) => value.replace(/&quot;/g, '"').replace(/&apos;/g, "'")
    .replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&amp;/g, '&');
  const root = new FakeElement('#document');
  let parent = root;
  const tokens = /<!--[\s\S]*?-->|<!\[CDATA\[[\s\S]*?\]\]>|<[?!][^>]*>|<\/([\w:-]+)\s*>|<([\w:-]+)((?:\s+[\w:-]+="[^"]*")*)\s*(\/?)>/g;
  for (const [token, closing, name, attributes, selfClosing] of content.matchAll(tokens)) {
    if (closing) { parent = parent.parentElement; continue; }
    if (!name) continue;
    const element = parent.appendChild(new FakeElement(name, Object.fromEntries(
      [...attributes.matchAll(/([\w:-]+)="([^"]*)"/g)].map(([, key, value]) => [key, decode(value)])
    )));
    if (!selfClosing && !token.endsWith('/>')) parent = element;
  }
  return root.children[0];
};

// EU U2a: a Session 46 Result saved before the executor recorded Python's
// paint shows its edits without records. Load records Python's paint from the
// Session (`originalStroke*`, the Legend's `originalColors`, and the catalog
// fills), so a reconcile without the edits returns each part and swatch to
// it (tests/fixtures/sessions/forced-label-underlay-strokes.provenance.json).
test('Load records Python\'s paint on a Session 46 Result saved without records', () => {
  const session = JSON.parse(gunzipSync(readFileSync(new URL(
    '../fixtures/sessions/forced-label-underlay-strokes.v46.gbdraw-session.json.gz', import.meta.url
  ))).toString('utf8'));
  const results = normalizeLogicalResults(session.results);
  const admission = admitFeatureCatalog(session.editorState.featureCatalog, results, { adopt: true, mode: 'circular' });
  const { editorState, features } = session.modes.circular;
  const edits = {
    featureColorOverrides: features.featureColorOverrides || {},
    featureStrokeOverrides: editorState.featureStrokes.overrides,
    legendEntries: editorState.legend.entries,
    legendColorOverrides: editorState.legend.colorOverrides,
    legendStrokeOverrides: editorState.legend.strokeOverrides,
    originalLegendColors: session.editorState.legend.originalColors,
    originalSvgStroke: session.editorState.originalSvgStroke
  };
  const svg = parseSvgText(results[0].content);
  const paint = () => Object.fromEntries([
    ...svg.querySelectorAll('[data-gbdraw-feature-id]').filter((element) => element.getAttribute('data-gbdraw-feature-part') === 'block')
      .map((element) => [element.getAttribute('data-gbdraw-rendered-feature-id') || element.id, element]),
    ...svg.querySelectorAll('g[data-legend-key]').map((row) => [row.getAttribute('data-legend-key'), row.querySelector('path')])
  ].map(([name, element]) => [name, ['fill', 'stroke', 'stroke-width'].map((attribute) => element.getAttribute(attribute)).join(' ')]));
  const saved = paint();
  assert.equal(saved.f38ba7c3f, '#54bcf8 #2a9d8f 2.0');
  assert.equal(saved.CDS, '#54bcf8 #e63946 3');
  assert.equal(saved.repeat_region, '#f4a261 gray 2.0');

  const plan = createSavedResultPlan(results, admission, edits);
  assert.equal(plan.kind, 'MUTATING');
  plan.operationsByResult[0].callerTransforms.forEach((transform) => transform(svg));
  assert.deepEqual(paint(), saved, 'the records leave the paint as saved');
  const recorded = [{ ...results[0], content: serializeNode(svg) }];
  assert.equal(createSavedResultPlan(recorded, admission, edits).kind, createEmptySvgMutationPlan(1).kind,
    'a Result with records is admitted as saved');
  reconcileMountedResult(svg, createEmptySvgMutationPlan(1).operationsByResult[0]);
  const python = paint();
  // A width the Session kept as a number returns as that number.
  assert.equal(python.f38ba7c3f, '#54bcf8 gray 2.0');
  ['f841fb8a8', 'f2f7a48cf__instance_4_4b227777d4dd1fc6', 'f2f7a48cf__instance_5_ef2d127de37b942b']
    .forEach((renderedId) => assert.equal(python[renderedId], '#54bcf8 gray 2', renderedId));
  assert.equal(python.CDS, '#54bcf8 gray 2');
  assert.equal(python.repeat_region, '#d3d3d3 gray 2.0');
  assert.equal(python.f9cf91913, saved.f9cf91913, 'an unedited feature keeps its paint');
});

// U2b: the live preview compiles the palette, the specific-color rules (their
// prepared matches), and Feature visibility with Generate's precedence, so every
// feature and Legend row fill the displayed Result shows is an operation.
const previewFeature = (admission) => admission.renderedFeaturesByResult[0].get('f0001');
const previewPlan = (admission, { paletteColors = { CDS: '#aaaaaa' }, drawnContext = null, ...options } = {}) => (
  compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#aaaaaa' }],
    originalLegendOrder: ['CDS'],
    livePreview: { paletteColors, drawnContext },
    ...options
  }).operationsByResult[0]
);
const swatchFill = (svg) => svg.querySelector('g[data-legend-key]').querySelector('path').getAttribute('fill');
const featureFill = (svg) => svg.querySelector('[data-gbdraw-feature-id]').getAttribute('fill');

test('live fills follow override > rule > palette > what Python drew', () => {
  const { admission } = currentFixture();
  const rule = { feat: 'CDS', qual: 'gene', val: 'x', cap: 'Rule row', color: '#ff0000' };
  const fills = (options) => previewPlan(admission, options).featureFills.map(({ color }) => color);
  assert.deepEqual(fills(), [], 'the palette Python drew needs no operation');
  assert.deepEqual(fills({ paletteColors: { CDS: '#00ff00' } }), ['#00ff00']);
  assert.deepEqual(fills({ paletteColors: { CDS: '#00ff00' }, manualSpecificRules: [rule] }), [],
    'a feature whose rule match is not known keeps what Python drew');
  recordRuleMatches([previewFeature(admission)], [ruleKey(rule)], () => ({ matched: [0], priorities: [0], declined: [] }));
  assert.deepEqual(fills({ paletteColors: { CDS: '#00ff00' }, manualSpecificRules: [rule] }), ['#ff0000']);
  assert.deepEqual(fills({
    paletteColors: { CDS: '#00ff00' }, manualSpecificRules: [rule],
    featureColorOverrides: { [biologicalFeatureKey('record-a', 'feature-a')]: '#123456' }
  }), ['#123456'], 'a fill edit wins');
  assert.deepEqual(compileDirectEditorMutationPlan({
    catalogAdmission: admission, manualSpecificRules: [rule]
  }).operationsByResult[0].featureFills, [], 'Generate leaves the rules and palette to Python');
});

test('live visibility hides the rendered features the resolver does not draw', () => {
  const { admission } = currentFixture();
  const context = (selected) => ({ featureOverrides: {}, rules: [], selectedTypes: new Set(selected), colorRules: [] });
  assert.deepEqual(previewPlan(admission, { drawnContext: context(['tRNA']) }).featureVisibility,
    [{ renderedId: 'f0001', mode: 'off' }]);
  assert.deepEqual(previewPlan(admission, { drawnContext: context(['CDS']) }).featureVisibility, []);
});

// OV-146: a swatch whose Legend-only color was retired shows its rule's or
// palette's color on the next Legend fill reconcile, not the fill Python drew.
test('a Legend fill reconcile shows the rule or palette color of a row whose own color left', () => {
  const { admission } = currentFixture();
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, previewPlan(admission, { legendColorOverrides: { CDS: '#334455' } }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#334455');
  const rule = { feat: 'CDS', qual: 'gene', val: 'x', cap: 'CDS', color: '#ff0000' };
  reconcileMountedResult(mounted, previewPlan(admission, { manualSpecificRules: [rule] }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#ff0000');
  reconcileMountedResult(mounted, previewPlan(admission, { paletteColors: { CDS: '#00ff00' } }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#00ff00');
  reconcileMountedResult(mounted, previewPlan(admission), { domains: ['featureFills', 'legendFills'] });
  assert.equal(swatchFill(mounted), '#aaaaaa');
  assert.equal(featureFill(mounted), '#aaaaaa');
});

// A Result without a feature catalog (a Session older than 40) is reached
// through the features read from it, so live edits still show.
test('a Result without a feature catalog receives the editor operations of its features', () => {
  const feature = { svg_id: 'f0001', type: 'CDS', id: 'f0001' };
  const addressing = displayedFeatureAddressing([feature], ['diagram.svg'], 0);
  const operations = compileDirectEditorMutationPlan({
    catalogAdmission: addressing,
    featureStrokeOverrides: { [featureOverrideKey(feature)]: { strokeColor: '#e63946', strokeWidth: 3 } },
    livePreview: { paletteColors: { CDS: '#00ff00' }, drawnContext: null }
  }).operationsByResult[0];
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, operations, { domains: ['featureStrokes', 'featureFills'] });
  const element = mounted.querySelector('[data-gbdraw-feature-id]');
  assert.deepEqual(['stroke', 'stroke-width', 'fill'].map((name) => element.getAttribute(name)), ['#e63946', '3', '#00ff00']);
  reconcileMountedResult(mounted, createEmptySvgMutationPlan(1).operationsByResult[0], { domains: ['featureStrokes', 'featureFills'] });
  assert.deepEqual(['stroke', 'stroke-width', 'fill'].map((name) => element.getAttribute(name)), [null, null, '#aaaaaa']);
});

// Without a catalog, Python's fill of each feature is read from the displayed
// Result: a palette color equal to it is no operation, and a Legend row
// stroke reaches the features drawn with the row's color.
test('a Result without a feature catalog takes Python\'s fills from the displayed Result', () => {
  const feature = { svg_id: 'f0001', type: 'CDS', id: 'f0001' };
  const mounted = buildSvgRoot();
  const operations = compileDirectEditorMutationPlan({
    catalogAdmission: displayedFeatureAddressing([feature], ['diagram.svg'], 0, mounted),
    legendEntries: [{ caption: 'CDS', color: '#aaaaaa' }],
    legendStrokeOverrides: { CDS: { strokeColor: '#e63946' } },
    livePreview: { paletteColors: { CDS: '#AAAAAA' }, drawnContext: null }
  }).operationsByResult[0];
  assert.deepEqual(operations.featureFills, []);
  assert.deepEqual(operations.legendStrokes.map(({ renderedIds }) => renderedIds), [['f0001']]);
});

// A catalog row that does not say which fill Python drew gives no record, so
// a later fill reconcile cannot remove the fill.
test('the Load normalizer records no fill Python is not known to have drawn', () => {
  const featureCatalog = catalog();
  featureCatalog.items[0].features[0].fillColor = '';
  const response = normalizeGenerationResponse({ results: [{ name: 'diagram.svg', content: '<svg/>' }], metadata: { featureCatalog } });
  const admission = admitFeatureCatalog(featureCatalog, response.results, { mode: 'linear' });
  const svg = buildSvgRoot();
  svg.querySelector('[data-gbdraw-feature-id]').setAttribute('fill', '#112233');
  const plan = createSavedResultPlan(response.results, admission, {
    featureColorOverrides: { [biologicalFeatureKey('record-a', 'feature-a')]: '#112233' },
    featureStrokeOverrides: {}, legendEntries: [], legendColorOverrides: {}, legendStrokeOverrides: {},
    originalLegendColors: {}, originalSvgStroke: { color: null, width: null }
  });
  plan.operationsByResult[0].callerTransforms.forEach((transform) => transform(svg));
  assert.doesNotMatch(serializeNode(svg), /data-gbdraw-base-fill/);
});
