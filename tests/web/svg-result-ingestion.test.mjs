import assert from 'node:assert/strict';
import test from 'node:test';

import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
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
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  isCommittedSvgResult,
  markCommittedSvgResultMounted,
  markCommittedSvgResultUnmounted
} from '../../gbdraw/web/js/services/svg-result-ingestion.js';

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
    if (selector === '[style]') return this.hasAttribute('style');
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

const buildSvgRoot = ({ missingFeature = false, missingLabel = false } = {}) => {
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
            missingLabel: content.includes('missing-label')
          })
        : new FakeElement('html'),
      querySelector: () => null
    };
  }
}

const catalog = () => ({
  schema: 3,
  items: [{
    resultIndex: 0,
    resultName: 'diagram.svg',
    recordKeys: ['record-a'],
    features: [{
      svgId: 'f0001',
      recordKey: 'record-a',
      biologicalFeatureId: 'feature-a',
      fillColor: '#aaaaaa'
    }],
    biologicalFeatures: [{
      recordKey: 'record-a',
      biologicalFeatureId: 'feature-a',
      stableFeatureId: 'stable-a',
      record_id: 'record-a',
      type: 'CDS',
      start: 0,
      end: 3,
      strand: 1
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
      parser: FakeDomParser
    }
  ));
  const committed = results[0];

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
    featureVisibilityOverrides: { f0001: 'off' }
  }),
  Label: (admission) => ({
    catalogAdmission: admission,
    labelTextFeatureOverrides: { f0001: 'new label' },
    labelVisibilityOverrides: { f0001: 'off' }
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
    assert.equal(FakeDomParser.calls, 1);
    assert.equal(metricTotal(metrics, 'applicationSvgParseCount'), 1);
    assert.equal(metricTotal(metrics, 'svgMutationIndexBuildCount'), 1);
    assert.equal(metricTotal(metrics, 'svgSerializationCount'), 1);
    assert.equal(metricTotal(metrics, 'svgIdentityScanCount'), 0);
    assert.equal(metricTotal(metrics, 'currentLegacyNormalizationCount'), 0);
  });
}

test('combined current mutations share one root/index and serialize once', () => {
  FakeDomParser.calls = 0;
  const { response, admission } = currentFixture();
  const plan = compileDirectEditorMutationPlan({
    ...planOptions.fill(admission),
    ...planOptions.stroke(admission),
    ...planOptions.visibility(admission),
    ...planOptions.Label(admission),
    ...planOptions.Legend(admission)
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
  assert.match(results[0].content, /data-gbdraw-label-visibility-preview="off"/);
  assert.match(results[0].content, />new label<\/text>/);
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
      assertContent: (content) => {
        assert.match(content, /data-legend-key="Genes"/);
        assert.match(content, /transform="translate\(20, 30\)"/);
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

test('missing requested target rejects current admission before commit', () => {
  const { response, admission } = currentFixture('<svg>missing-label</svg>');
  const plan = compileDirectEditorMutationPlan(planOptions.Label(admission));
  assert.throws(
    () => admitCurrentGeneratedResults(response, {
      catalogAdmission: admission,
      mutationPlan: plan,
      sanitizer: { sanitize: (value) => value },
      parser: FakeDomParser
    }),
    /missing or ambiguously binds an editable Label/
  );
  assert.equal(isCommittedSvgResult(response.results[0]), false);
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
  assert.throws(
    () => admitLegacyImportedResults(currentSource),
    /legacy-import/
  );

  const legacy = admitLegacyImportedResults(
    createLegacyImportResultSource([{ name: 'legacy.svg', content: '<svg></svg>' }]),
    { sanitizer: { sanitize: (value) => value }, parser: FakeDomParser }
  );
  assert.equal(FakeDomParser.calls, 1);
  assert.equal(getCommittedSvgResultMetadata(legacy[0]).sourceClass, 'legacy-import');
});
