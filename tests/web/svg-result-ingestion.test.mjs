import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';

import { compileDirectEditorMutationPlan, displayedLegendRowColors, editorPaintDomains, LIVE_EDIT_DOMAINS } from '../../gbdraw/web/js/app/candidate-render.js';
import {
  requireUniqueEditableLabelBindings
} from '../../gbdraw/web/js/app/feature-editor/label-actions.js';
import { buildMatchPopupPayload } from '../../gbdraw/web/js/app/pairwise-match-popup.js';
import { normalizeGenerationResponse, selectHeldRenderState } from '../../gbdraw/web/js/services/diagram-generation.js';
import { SEGMENTS } from '../../gbdraw/web/js/services/current-worker-result-source.js';
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
import { PAINT_SLICE } from '../../gbdraw/web/js/services/drawing/paint.js';
import { LEGEND_ORDER_RECORD, stripResultBaseAttributes } from '../../gbdraw/web/js/services/result-paint-bases.js';
import {
  getLegendEntrySwatch as legendSwatch, legendRowKeys, legendRowTakesPlace, shownPythonLegendRow
} from '../../gbdraw/web/js/services/legend-svg.js';
import { recordRuleMatches, ruleKey } from '../../gbdraw/web/js/services/rule-matchers.js';
import { displayedFeatureAddressing, featureOverrideKey } from '../../gbdraw/web/js/services/feature-override-identity.js';

// Session Load (services/config.js) reads the editor state through Vue state.
globalThis.window = { Vue: {
  ref: (value) => ({ value }), reactive: (value) => value, toRaw: (value) => value,
  computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const { savedResultEdits } = await import('../../gbdraw/web/js/services/config.js');

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
  appendChild(child) { child.remove?.(); child.parentElement = this; this.children.push(child); return child; }
  insertBefore(child, reference) {
    child.remove?.();
    child.parentElement = this;
    const at = this.children.indexOf(reference);
    this.children.splice(at < 0 ? this.children.length : at, 0, child);
    return child;
  }
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
    if (selector === 'g[data-legend-key][display="none"]') {
      return this.tagName === 'g' && this.hasAttribute('data-legend-key') && this.getAttribute('display') === 'none';
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
    metadata: {},
    featureCatalog
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

// Step 3b S: a Generate Result arrives as document pieces; the page sanitizes
// each children text once per revision and keeps it on the Result it commits.
const encoder = new TextEncoder();
const revisionOf = (text) => createHash('sha256').update(encoder.encode(text)).digest('hex');
const DOCUMENT_HEAD = '<svg xmlns="http://www.w3.org/2000/svg" width="10">';
const documentReply = (children, { held = new Set(), starts = [] } = {}) => {
  const chunks = [];
  let offset = 0;
  const pieces = [{ kind: 'raw', text: '<defs><clipPath id="c"/></defs>' }, ...children.map((text, index) => {
    const bytes = encoder.encode(text);
    const piece = {
      kind: 'content', start: starts[index] ?? `<g id="g${index}">`, revision: revisionOf(text), bytes: bytes.byteLength
    };
    if (held.has(piece.revision)) return piece;
    chunks.push(bytes);
    const sent = { ...piece, offset, length: bytes.byteLength };
    offset += bytes.byteLength;
    return sent;
  })];
  const segments = new Uint8Array(offset);
  chunks.reduce((at, bytes) => { segments.set(bytes, at); return at + bytes.byteLength; }, 0);
  return { name: 'diagram.svg', document: { head: DOCUMENT_HEAD, tail: '</svg>', pieces }, ...(chunks.length ? { segments } : {}) };
};
const joinedDocument = (children, starts = []) => `${DOCUMENT_HEAD}<defs><clipPath id="c"/></defs>${
  children.map((text, index) => `${starts[index] ?? `<g id="g${index}">`}${text}</g>`).join('')}</svg>`;
const admitDocument = (reply, held, fakeSanitizer, mutationPlan = createEmptySvgMutationPlan(1)) => {
  const featureCatalog = catalog();
  const response = normalizeGenerationResponse({ results: [reply], metadata: {}, featureCatalog }, held);
  const catalogAdmission = admitFeatureCatalog(featureCatalog, response.results, { adopt: true, mode: 'linear' });
  return captureMetrics(() => admitCurrentGeneratedResults(response, {
    catalogAdmission, mutationPlan, sanitizer: fakeSanitizer, parser: FakeDomParser, selectedFeatureTypes: null
  }));
};
const COLD_CHILDREN = ['<path id="f0001" d="M0,0"/><script>bad()</script>', '<text>title µ</text>'];

test('current-worker EMPTY admission of a document reply sanitizes each new piece once and performs zero application SVG work', () => {
  FakeDomParser.calls = 0;
  const calls = { calls: 0 };
  const { value: results, metrics } = admitDocument(documentReply(COLD_CHILDREN), undefined, sanitizer(calls));
  const committed = results[0];

  // Once to learn the wrapper's sanitized prefix, then the head, the raw
  // piece, two start tags and two children texts: each piece once.
  assert.equal(calls.calls, 7);
  assert.equal(FakeDomParser.calls, 0);
  // The joined pieces are what one sanitizer call on the whole document gives.
  assert.equal(committed.content, sanitizer({ calls: 0 }).sanitize(joinedDocument(COLD_CHILDREN)));
  assert.equal(committed.content.includes('<script>'), false);
  assert.equal(isCommittedSvgResult(committed), true);
  assert.equal(metricTotal(metrics, 'svgSanitizationCount'), 1);
  assert.equal(metricTotal(metrics, 'svgSanitizedPieceCount'), 2);
  assert.equal(metricTotal(metrics, 'svgMemoHitCount'), 0);
  for (const name of [
    'applicationSvgParseCount',
    'svgMutationIndexBuildCount',
    'featureDomFullScanCount',
    'legendDomFullScanCount',
    'svgIdentityScanCount',
    'svgSerializationCount',
    'currentLegacyNormalizationCount'
  ]) assert.equal(metricTotal(metrics, name), 0, name);

  // The committed Result keeps today's members plus SEGMENTS: no decoded
  // document, no raw children; each segment holds sanitized text only.
  assert.deepEqual(Object.keys(committed).sort(), ['content', 'name']);
  assert.deepEqual(Object.getOwnPropertySymbols(committed).filter((key) => key === SEGMENTS), [SEGMENTS]);
  const segments = committed[SEGMENTS];
  assert.deepEqual([...segments.keys()], COLD_CHILDREN.map(revisionOf));
  for (const segment of segments.values()) {
    assert.deepEqual(Object.keys(segment).sort(), ['bytes', 'revision', 'sanitized']);
    assert.equal(Object.isFrozen(segment), true);
  }
  assert.equal(segments.get(revisionOf(COLD_CHILDREN[0])).sanitized, '<path id="f0001" d="M0,0"/>');
  assert.equal(segments.get(revisionOf(COLD_CHILDREN[1])).bytes, encoder.encode(COLD_CHILDREN[1]).byteLength);
  // Session save and History copy the members, never the symbol key.
  assert.deepEqual(JSON.parse(JSON.stringify(committed)), { name: 'diagram.svg', content: committed.content });
});

test('a warm Generate sanitizes only the children texts the page does not hold', () => {
  const calls = { calls: 0 };
  const fake = sanitizer(calls);
  const cold = admitDocument(documentReply(COLD_CHILDREN), undefined, fake).value[0];
  const held = selectHeldRenderState(null, [cold], null);
  assert.deepEqual([...held.segments.keys()].sort(), COLD_CHILDREN.map(revisionOf).sort());

  // A title edit: the first group's children are the same, the second's change.
  const retitled = [COLD_CHILDREN[0], '<text>new title</text>'];
  calls.calls = 0;
  const { value: [warm], metrics } = admitDocument(
    documentReply(retitled, { held: new Set(held.segments.keys()) }), held, fake
  );
  // Head, raw piece, two start tags, one changed children text.
  assert.equal(calls.calls, 5);
  assert.equal(metricTotal(metrics, 'svgSanitizedPieceCount'), 1);
  assert.equal(metricTotal(metrics, 'svgMemoHitCount'), 1);
  assert.equal(warm.content, fake.sanitize(joinedDocument(retitled)));
  // The unchanged segment is the same object in both Results; the old title's is gone.
  assert.equal(warm[SEGMENTS].get(revisionOf(COLD_CHILDREN[0])), cold[SEGMENTS].get(revisionOf(COLD_CHILDREN[0])));
  assert.deepEqual([...warm[SEGMENTS].keys()], retitled.map(revisionOf));

  // A start-tag-only edit (a moved group) holds every children text.
  const moved = ['<g id="g0" transform="translate(5,0)">', '<g id="g1">'];
  calls.calls = 0;
  const { value: [shifted], metrics: movedMetrics } = admitDocument(
    documentReply(retitled, { held: new Set(retitled.map(revisionOf)), starts: moved }),
    selectHeldRenderState(null, [warm], null),
    fake
  );
  assert.equal(calls.calls, 4);
  assert.equal(metricTotal(movedMetrics, 'svgSanitizedPieceCount'), 0);
  assert.equal(metricTotal(movedMetrics, 'svgMemoHitCount'), 2);
  assert.equal(shifted.content, fake.sanitize(joinedDocument(retitled, moved)));
});

test('a mutated document Result keeps the segments of what Python wrote', () => {
  const fake = sanitizer({ calls: 0 });
  const { value: [committed] } = admitDocument(
    documentReply(COLD_CHILDREN),
    undefined,
    fake,
    { kind: 'MUTATING', operationsByResult: [{ ...createEmptySvgMutationPlan(1).operationsByResult[0], callerTransforms: [() => true] }] }
  );
  assert.equal(committed.content.startsWith('<svg'), true);
  assert.equal(committed[SEGMENTS].get(revisionOf(COLD_CHILDREN[1])).sanitized, '<text>title µ</text>');
});

test('a document reply that does not match the held snapshot or its wrapper fails at result admission', () => {
  const resultInvalid = (error) => error.code === 'RESULT_INVALID' && error.stage === 'result-admission';
  const fake = sanitizer({ calls: 0 });
  const cold = admitDocument(documentReply(COLD_CHILDREN), undefined, fake).value[0];
  const held = selectHeldRenderState(null, [cold], null);
  const allHeld = new Set(COLD_CHILDREN.map(revisionOf));
  // A held piece the request's snapshot does not hold, or names with another length.
  assert.throws(() => admitDocument(documentReply(COLD_CHILDREN, { held: allHeld }), undefined, fake), resultInvalid);
  const wrongLength = documentReply(COLD_CHILDREN, { held: allHeld });
  wrongLength.document.pieces[1] = { ...wrongLength.document.pieces[1], bytes: 1 };
  assert.throws(() => admitDocument(wrongLength, held, fake), resultInvalid);
  // A sanitizer whose output loses the wrapper's prefix, or a start tag's end tag.
  const dropsPrefix = { sanitize: (content) => (String(content).includes('title') ? String(content).slice(5) : String(content)) };
  assert.throws(() => admitDocument(documentReply(COLD_CHILDREN), undefined, dropsPrefix), resultInvalid);
  const dropsEndTag = { sanitize: (content) => String(content).replace('<g id="g1"></g>', '<g id="g1">') };
  assert.throws(() => admitDocument(documentReply(COLD_CHILDREN), undefined, dropsEndTag), resultInvalid);
});

test('the held snapshot reads Results through toRaw and is fixed when the request starts', () => {
  const fake = sanitizer({ calls: 0 });
  const committed = admitDocument(documentReply(COLD_CHILDREN), undefined, fake).value[0];
  // A reactive Result hands out a proxy for its symbol key; toRaw gives the raw Map.
  const reactiveMap = new Proxy(committed[SEGMENTS], {});
  const reactive = new Proxy(committed, { get: (target, key) => (key === SEGMENTS ? reactiveMap : target[key]) });
  const toRaw = globalThis.window.Vue.toRaw;
  globalThis.window.Vue.toRaw = (value) => (value === reactive ? committed : value);
  try {
    const held = selectHeldRenderState(null, [reactive], null);
    assert.equal(held.segments.get(revisionOf(COLD_CHILDREN[0])), committed[SEGMENTS].get(revisionOf(COLD_CHILDREN[0])));
    // A Result committed after the snapshot was taken (the race) does not change it.
    const later = admitDocument(documentReply(['<text>later</text>']), undefined, fake).value[0];
    assert.equal(held.segments.has(revisionOf('<text>later</text>')), false);
    assert.equal(later[SEGMENTS].size, 1);
    const { value: [warm] } = admitDocument(documentReply(COLD_CHILDREN, { held: new Set(held.segments.keys()) }), held, fake);
    assert.equal(warm.content, committed.content);
  } finally {
    globalThis.window.Vue.toRaw = toRaw;
  }
  // A Result admitted from its whole text (a Session, a legacy import) holds nothing.
  const { response, admission } = currentFixture();
  const [whole] = admitCurrentGeneratedResults(response, {
    catalogAdmission: admission, mutationPlan: createEmptySvgMutationPlan(1), sanitizer: fake, parser: FakeDomParser
  });
  assert.equal(Object.getOwnPropertySymbols(whole).includes(SEGMENTS), false);
  assert.equal(selectHeldRenderState(null, [whole, null, { content: '<svg/>' }], null).segments.size, 0);
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
      // The deleted row stays hidden with Python's key (exports strip it).
      assertContent: (content) => assert.match(content, /<g data-legend-key="CDS" data-gbdraw-base-display="" display="none">/)
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

// PD-OI-042 scenario revision 2 (D-09): the executor names the rows the Legend
// shows by the naming rule, from each row's name and the swatch it shows. A
// row's key (`data-legend-key`, here `Shared [#112233]`) is never its text.
test('the executor draws Legend names by the naming rule and never a row key', () => {
  const texts = (root) => root.querySelectorAll('g[data-legend-key]')
    .filter((row) => row.getAttribute('display') !== 'none').map((row) => row.querySelector('text').textContent);
  const operations = (legend) => ({ ...createEmptySvgMutationPlan(1).operationsByResult[0], ...legend });
  const rule = (color) => ({ caption: `Shared [${color}]`, name: 'Shared', color, xPos: null, yPos: null, ifAbsent: true });
  const mounted = buildSvgRoot();
  // One caption in two colors: two rows, numbered in legend order.
  const two = operations({ legendAdds: [rule('#112233'), rule('#445566')] });
  assert.deepEqual(reconcileMountedResult(mounted, two), { changed: true, legendChanged: true });
  assert.deepEqual(texts(mounted), ['CDS', 'Shared (1)', 'Shared (2)']);
  assert.deepEqual(mounted.querySelectorAll('g[data-legend-key]').map((row) => row.getAttribute('data-legend-name')), [null, 'Shared', 'Shared']);
  assert.deepEqual(reconcileMountedResult(mounted, two), { changed: false, legendChanged: false }, 'named once');
  // A rename onto the shared name joins the numbering; its own key stays its record.
  const renamed = operations({ ...two, legendRenames: [{ from: 'CDS', to: 'Shared', allowMissing: false }] });
  reconcileMountedResult(mounted, renamed);
  assert.deepEqual(texts(mounted), ['Shared (1)', 'Shared (2)', 'Shared (3)']);
  // A literal `Shared (1)` keeps its text, and the numbers of the two rows skip it.
  reconcileMountedResult(mounted, operations({
    ...two, legendRenames: [{ from: 'CDS', to: 'Shared (1)', allowMissing: false }]
  }));
  assert.deepEqual(texts(mounted), ['Shared (1)', 'Shared (2)', 'Shared (3)']);
  // A name left with one row has no number again.
  reconcileMountedResult(mounted, operations({
    ...two, legendDeletes: [{ caption: 'Shared [#445566]', allowMissing: true, retire: true }]
  }));
  assert.deepEqual(texts(mounted), ['CDS', 'Shared']);
});

// The Legend tab, the popups, and the dialogs show a row by the name the
// Legend draws for it, never by its key (PD-OI-042 scenario revision 2), and
// typing a row's own name or drawn name into its name field is no rename.
test('the Legend owner shows a row by its drawn name and reads which typed names are the row itself', async () => {
  const { createLegendEntryActions } = await import('../../gbdraw/web/js/app/legend/entry-actions.js');
  const mounted = buildSvgRoot();
  const rule = (color) => ({ caption: `Shared [${color}]`, name: 'Shared', color, xPos: null, yPos: null, ifAbsent: true });
  const two = { ...createEmptySvgMutationPlan(1).operationsByResult[0], legendAdds: [rule('#112233'), rule('#445566')] };
  const legendRowNames = { value: new Map() };
  const owner = createLegendEntryActions({
    state: {
      results: { value: [] }, svgContainer: { value: { querySelector: () => mounted } },
      originalLegendOrder: { value: [] }, originalLegendColors: { value: {} },
      activeDrawing: () => ({ manualSpecificRules: ['#112233', '#445566'].map((color) => ({ feat: 'CDS', qual: 'gene', val: color, color, cap: 'Shared' })) })
    },
    readLegendRowColors: (drawing) => displayedLegendRowColors({}, drawing, []),
    legendRowNames
  });
  reconcileMountedResult(mounted, two);
  // A row the displayed Result does not hold (not read yet, deleted, or dormant)
  // is named as the color rules that draw it are captioned, never by its key (D-15).
  assert.deepEqual(['Shared [#112233]', 'Shared [#445566]', 'No such row'].map(owner.shownLegendName), ['Shared', 'Shared', 'No such row']);
  assert.deepEqual(['Shared [#112233]', 'Typed name'].map(owner.legendRowName), ['Shared', 'Typed name']);
  owner.readLegendRowNames();
  assert.deepEqual(['CDS', 'Shared [#112233]', 'Shared [#445566]', 'No such row'].map(owner.shownLegendName),
    ['CDS', 'Shared (1)', 'Shared (2)', 'No such row']);
  assert.deepEqual(['Shared [#112233]', 'Shared', 'Shared (1)', 'Shared (2)', 'Other'].map((typed) => owner.legendRowIsNamed('Shared [#112233]', typed)),
    [true, true, true, false, false]);
  assert.equal(owner.legendRowIsNamed('No such row', 'No such row'), true);
  // A row the Legend no longer draws is shown by its name, without a number.
  reconcileMountedResult(mounted, { ...two, legendDeletes: [{ caption: 'Shared [#445566]', allowMissing: true, retire: true }] });
  owner.readLegendRowNames();
  assert.deepEqual(['Shared [#112233]', 'Shared [#445566]'].map(owner.shownLegendName), ['Shared', 'Shared']);
  assert.deepEqual(['Shared [#112233]', 'Shared [#445566]'].map(owner.legendRowName), ['Shared', 'Shared']);
});

// PD-OI-061 scenario revision 5: the Legend owner names each row by its key,
// also a renamed one, and finds the row a typed name names: a drawn name
// (a shown numbered name names that row), else a name, the row in the given
// color first, else the first in legend order.
test('the Legend owner finds the row a typed name names, by its key', async () => {
  const { createLegendEntryActions } = await import('../../gbdraw/web/js/app/legend/entry-actions.js');
  const mounted = buildSvgRoot();
  const rule = (color) => ({ caption: `Shared [${color}]`, name: 'Shared', color, xPos: null, yPos: null, ifAbsent: true });
  const two = { ...createEmptySvgMutationPlan(1).operationsByResult[0], legendAdds: [rule('#112233'), rule('#445566')] };
  const legendEntries = { value: ['CDS', 'Shared [#112233]', 'Shared [#445566]'].map((caption) => ({ caption, originalCaption: caption, color: caption === 'CDS' ? '#54bcf8' : caption.slice(8, 15) })) };
  const owner = createLegendEntryActions({
    state: {
      results: { value: [] }, svgContainer: { value: { querySelector: () => mounted } },
      originalLegendOrder: { value: [] }, originalLegendColors: { value: {} },
      activeDrawing: () => ({ manualSpecificRules: [], legendEntries, deletedLegendEntries: { value: [] } })
    },
    readLegendRowColors: (drawing) => displayedLegendRowColors({}, drawing, []),
    legendRowNames: { value: new Map() }
  });
  reconcileMountedResult(mounted, two);
  owner.readLegendRowNames();
  assert.deepEqual(['Shared (2)', 'shared (1)', 'Shared', 'SHARED', 'Nothing'].map((typed) => owner.findLegendRowNamed(typed, '#445566')),
    ['Shared [#445566]', 'Shared [#112233]', 'Shared [#445566]', 'Shared [#445566]', null]);
  assert.equal(owner.findLegendRowNamed('Shared', '#000000'), 'Shared [#112233]', 'no row in the color: the first in legend order');
  // A renamed row is shown and found by its key.
  reconcileMountedResult(mounted, { ...two, legendRenames: [{ from: 'CDS', to: 'Gamma', allowMissing: false }] });
  legendEntries.value[0].caption = 'Gamma';
  owner.readLegendRowNames();
  assert.deepEqual(['CDS', 'Gamma'].map(owner.shownLegendName), ['Gamma', 'Gamma']);
  assert.equal(owner.legendRowIsNamed('CDS', 'Gamma'), true);
  assert.equal(owner.findLegendRowNamed('gamma', '#000000'), 'CDS');
  assert.deepEqual(legendEntries.value.map(owner.shownLegendEntryName), ['Gamma', 'Shared (1)', 'Shared (2)']);
});

// PR 3 (D-25): the Legend owner stores a row's color, stroke, and open Stroke
// options under the row's key, so two rows of one name keep their own edits;
// a deleted row returns in the color stored under its key.
test('the Legend owner keeps the edits of two rows of one name under their keys', async () => {
  const { createLegendEntryActions } = await import('../../gbdraw/web/js/app/legend/entry-actions.js');
  const { createLegendStrokeActions } = await import('../../gbdraw/web/js/app/legend/stroke-actions.js');
  const mounted = buildSvgRoot();
  const added = { caption: 'GC content', color: '#aaaaaa', xPos: null, yPos: null };
  reconcileMountedResult(mounted, { ...createEmptySvgMutationPlan(1).operationsByResult[0], legendAdds: [added] });
  const paint = PAINT_SLICE.create();
  PAINT_SLICE.install(paint.handle, { legendEntries: [
    { caption: 'CDS', originalCaption: 'CDS', color: '#aaaaaa' },
    { caption: 'CDS', originalCaption: 'GC content', color: '#aaaaaa' }
  ] });
  const drawing = { manualSpecificRules: [], ...paint.view, slices: { paint: paint.handle } };
  const strokeOptionsOpen = new Set();
  const state = {
    results: { value: [] }, svgContainer: { value: { querySelector: () => mounted } },
    originalLegendOrder: { value: ['CDS', 'GC content'] }, originalLegendColors: { value: {} },
    legendStrokeOptionsOpen: strokeOptionsOpen,
    activeDrawing: () => drawing
  };
  const owner = createLegendEntryActions({ state, readLegendRowColors: (drawing) => displayedLegendRowColors({}, drawing, []) });
  const strokes = createLegendStrokeActions({ state });
  assert.equal(owner.updateLegendEntryColor(1, '#ff0000'), true);
  assert.equal(owner.updateLegendEntryColor(0, '#0000ff'), true);
  assert.deepEqual(drawing.legendColorOverrides, { 'GC content': '#ff0000', CDS: '#0000ff' });
  assert.equal(strokes.updateLegendEntryStrokeColor(1, '#123456'), true);
  assert.deepEqual(Object.keys(drawing.legendStrokeOverrides), ['GC content']);
  assert.deepEqual([strokes.getLegendEntryStrokeColor(0), strokes.getLegendEntryStrokeColor(1)], ['', '#123456']);
  strokes.toggleLegendStrokeOptions(drawing.legendEntries.value[1]);
  assert.deepEqual([...strokeOptionsOpen], ['GC content']);
  assert.deepEqual(drawing.legendEntries.value.map((entry) => strokes.isLegendStrokeOptionsOpen(entry)), [false, true]);
  owner.deleteLegendEntry(1);
  drawing.deletedLegendEntries.value[0].color = '#aaaaaa';
  owner.restoreDeletedLegendEntries();
  assert.deepEqual(drawing.legendEntries.value.map(({ originalCaption, color }) => [originalCaption, color]),
    [['CDS', '#0000ff'], ['GC content', '#ff0000']], 'a row returns in the color stored under its key');
});

// OV-345: a reconcile says whether it changed the Result, so a display,
// History restore, or Reset that changes nothing keeps the Result's bytes.
test('a reconcile reports whether it changed the Result and whether the Legend rows changed', () => {
  const { admission } = currentFixture();
  const fill = compileDirectEditorMutationPlan(planOptions.fill(admission)).operationsByResult[0];
  const reset = createEmptySvgMutationPlan(1).operationsByResult[0];
  const mounted = buildSvgRoot();
  const python = serializeNode(mounted);
  assert.deepEqual(reconcileMountedResult(mounted, reset), { changed: false, legendChanged: false }, 'nothing to show');
  assert.deepEqual(reconcileMountedResult(mounted, fill), { changed: true, legendChanged: false });
  assert.deepEqual(reconcileMountedResult(mounted, fill), { changed: false, legendChanged: false }, 'shown again');
  assert.deepEqual(reconcileMountedResult(mounted, reset), { changed: true, legendChanged: false }, 'returned to Python');
  assert.equal(serializeNode(mounted), python);
});

// OV-345: Load keeps a saved Result's bytes when its records change nothing
// (an EU Result whose edits Python drew at Generate), and admits the records
// of a Result that shows an edit without them.
test('Load keeps a saved Result\'s bytes when its records change nothing', () => {
  const content = '<svg data-python="bytes"><path id="f0001" fill="#aaaaaa"></path></svg>';
  const admit = (fillColor) => {
    const featureCatalog = catalog();
    featureCatalog.items[0].features[0].fillColor = fillColor;
    const results = [{ name: 'diagram.svg', content }];
    const admission = admitFeatureCatalog(featureCatalog, results, { mode: 'linear' });
    const plan = createSavedResultPlan(results, admission, {
      featureColorOverrides: { [biologicalFeatureKey('record-a', 'feature-a')]: '#aaaaaa' },
      featureStrokeOverrides: {}, legendEntries: [], legendColorOverrides: {}, legendStrokeOverrides: {},
      originalLegendColors: {}, originalLegendOrder: [], rules: [], mode: 'linear', sessionVersion: 44
    });
    assert.equal(plan.kind, 'MUTATING');
    return captureMetrics(() => admitCurrentSessionResults(createCurrentSessionResultSource(results, admission), {
      mutationPlan: plan, sanitizer: sanitizer({ calls: 0 }), parser: FakeDomParser
    }));
  };
  const drawn = admit('#aaaaaa');
  assert.equal(drawn.value[0].content, content, 'Python drew the edit: the bytes are kept');
  assert.equal(drawn.metrics.filter(({ name }) => name === 'svgSerializationCount').length, 0);
  const recorded = admit('#bbbbbb');
  assert.match(recorded.value[0].content, /data-gbdraw-base-fill="#bbbbbb"/);
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
    // The Result's one row is another category's, drawn with its own name.
    const row = mounted.querySelector('g[data-legend-key]');
    row.setAttribute('data-legend-key', 'tRNA');
    row.querySelector('text').textContent = 'tRNA';
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
    featureCatalog,
    metadata: {
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

// The Session fixtures Load reads below, and Load's records on one saved
// Result set of `mode` (`patch` changes a saved Result before Load reads it).
const sessionFixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(
  `../fixtures/sessions/${name}`, import.meta.url
))).toString('utf8'));
const FORCED_STROKES = 'forced-label-underlay-strokes.v44.gbdraw-session.json.gz';

// `main` from 2026-06-18 to 2026-07-19 (Sessions 27 to 32), which includes
// gbdraw 0.12.0, 0.12.1 and 0.13.0 (the Session 29 and 30 fixtures below), named
// a comparison match only by `data-gbdraw-pairwise-match-id`. The legacy admission (here
// with a sanitizer that keeps its input) renames it once: the admitted Result
// has the stored IDs under `data-gbdraw-match-id` only, and the popup payload
// reads them from there. The real sanitizer and selection on a legacy Result
// are covered by "Pairwise pointer feedback" in linear-multi-record.playwright.spec.js.
test('a Session 29 or 30 Result is admitted with one match attribute', () => {
  const textParser = class { parseFromString(content) { return { documentElement: parseSvgText(content) }; } };
  [
    ['BGC0000708-BGC0000713.v29.gbdraw-session.json.gz', 77],
    ['BGC0000708-BGC0000713.v30.gbdraw-session.json.gz', 77],
    ['cli-linear-protein.v30.gbdraw-session.json.gz', 2]
  ].forEach(([name, count]) => {
    const session = sessionFixture(name);
    const stored = parseSvgText(session.results[0].content);
    const legacyIds = stored.querySelectorAll('[data-gbdraw-pairwise-match-id]')
      .map((element) => element.getAttribute('data-gbdraw-pairwise-match-id'));
    assert.equal(legacyIds.length, count, name);
    assert.equal(stored.querySelectorAll('[data-gbdraw-match-id]').length, 0, name);

    const [admitted] = admitLegacyImportedResults(
      createLegacyImportResultSource(normalizeLogicalResults(session.results)),
      { sanitizer: { sanitize: (value) => value }, parser: textParser }
    );
    const content = getCommittedSvgContent(admitted);
    assert.equal(content.includes('data-gbdraw-pairwise-match-id'), false, name);
    const matches = parseSvgText(content).querySelectorAll('[data-gbdraw-match-id]');
    assert.deepEqual(matches.map((element) => element.getAttribute('data-gbdraw-match-id')), legacyIds, name);
    assert.deepEqual(matches.map((element) => buildMatchPopupPayload(element)?.id), legacyIds, name);
  });
});
const loadSavedResults = (session, {
  mode = 'circular', results = normalizeLogicalResults(session.results),
  featureCatalog = session.editorState.featureCatalog, patch = () => {}
} = {}) => {
  const admission = admitFeatureCatalog(featureCatalog, results, { adopt: true, mode });
  const plan = createSavedResultPlan(results, admission, savedResultEdits(session, mode, session.editorState, session.version));
  return results.map((result, index) => {
    const svg = parseSvgText(result.content);
    patch(svg, index);
    (plan.operationsByResult[index].callerTransforms || []).forEach((transform) => transform(svg));
    return svg;
  });
};
// What Reset shows: a reconcile without edits.
const resetPaint = (svg) => {
  reconcileMountedResult(svg, createEmptySvgMutationPlan(1).operationsByResult[0]);
  return svg;
};
const partStrokes = (svg) => Object.fromEntries(svg.querySelectorAll('[data-gbdraw-feature-id]')
  .map((element) => [element.id, `${element.getAttribute('stroke')} ${element.getAttribute('stroke-width')}`]));
const connectorPart = (featureId, stroke, width) => new FakeElement('path', {
  id: `${featureId}__line1`, 'data-gbdraw-feature-id': featureId, 'data-gbdraw-feature-part': 'connector',
  'data-gbdraw-record-index': '0', fill: 'none', stroke, 'stroke-width': width
});

// EU U2a: main's writer (Session 44) saved a Result that shows its stroke
// and color edits without records of Python's paint. Load records Python's
// paint (the strokes of the parts no stroke edit reached, the catalog fills,
// the Legend's `originalColors`), so a reconcile without the edits returns
// each part and swatch to it
// (tests/fixtures/sessions/forced-label-underlay-strokes.provenance.json).
test('Load records Python\'s paint on a Session 44 Result saved without records', () => {
  const session = sessionFixture(FORCED_STROKES);
  assert.equal(session.version, 44);
  const results = normalizeLogicalResults(session.results);
  const admission = admitFeatureCatalog(session.editorState.featureCatalog, results, { adopt: true, mode: 'circular' });
  const edits = savedResultEdits(session, 'circular', session.editorState, session.version);
  const svg = parseSvgText(results[0].content);
  const paint = () => Object.fromEntries([
    ...svg.querySelectorAll('[data-gbdraw-feature-id]').filter((element) => element.getAttribute('data-gbdraw-feature-part') === 'block')
      .map((element) => [element.getAttribute('data-gbdraw-rendered-feature-id') || element.id, element]),
    ...svg.querySelectorAll('g[data-legend-key]').map((row) => [row.getAttribute('data-legend-key'), row.querySelector('path')])
  ].map(([name, element]) => [name, ['fill', 'stroke', 'stroke-width'].map((attribute) => element.getAttribute(attribute)).join(' ')]));
  const saved = paint();
  assert.equal(saved.f38ba7c3f, '#54bcf8 #2a9d8f 2.0');
  assert.equal(saved.repeat_region, '#f4a261 #e63946 3');

  const plan = createSavedResultPlan(results, admission, edits);
  assert.equal(plan.kind, 'MUTATING');
  const { callerTransforms } = plan.operationsByResult[0];
  assert.deepEqual(callerTransforms.map((transform) => transform(svg)), [true], 'the records');
  assert.deepEqual(paint(), saved, 'the records leave the paint as saved');
  assert.deepEqual(callerTransforms.map((transform) => transform(svg)), [false], 'a Result with records gets none again');
  resetPaint(svg);
  const python = paint();
  assert.equal(python.f38ba7c3f, '#54bcf8 gray 2.0', 'the stroke of the blocks no edit reached');
  ['f841fb8a8', 'f970830ef__instance_4_4b227777d4dd1fc6', 'f970830ef__instance_5_ef2d127de37b942b', 'CDS']
    .forEach((name) => assert.equal(python[name], saved[name], `${name} has no edit`));
  assert.equal(python.repeat_region, '#d3d3d3 gray 2.0');
});

// OV-982: a Legend row whose key is an Object.prototype member name
// (`constructor`) has no Python color unless the Session saved one; Load
// records none for its swatch, while a row with a saved color gets it.
test('Load records no Python color for a Legend row named like an Object member', () => {
  const content = '<svg><path id="f0001" fill="#aaaaaa"></path><g id="legend"><g id="feature_legend">'
    + '<g data-legend-key="constructor"><path fill="#ff0000"></path><text>constructor</text></g>'
    + '<g data-legend-key="CDS"><path fill="#00ff00"></path><text>CDS</text></g></g></g></svg>';
  const results = [{ name: 'diagram.svg', content }];
  const admission = admitFeatureCatalog(catalog(), results, { mode: 'linear' });
  const plan = createSavedResultPlan(results, admission, {
    featureColorOverrides: {}, featureStrokeOverrides: {}, legendEntries: [],
    legendColorOverrides: { constructor: '#ff0000', CDS: '#00ff00' }, legendStrokeOverrides: {},
    originalLegendColors: { CDS: '#123456' }, originalLegendOrder: [], rules: [], mode: 'linear', sessionVersion: 44
  });
  const svg = parseSvgText(content);
  (plan.operationsByResult[0].callerTransforms || []).forEach((transform) => transform(svg));
  const swatchBase = (key) => svg.querySelectorAll('g[data-legend-key]')
    .find((row) => row.getAttribute('data-legend-key') === key).querySelector('path').getAttribute('data-gbdraw-base-fill');
  assert.equal(swatchBase('CDS'), '#123456');
  assert.equal(swatchBase('constructor'), null);
});

// Python's paint is recorded only for a Session older than 46: a Session 46
// Result carries its own records, and its bytes are admitted as saved.
test('Load records Python\'s paint only on the Result of a Session older than 46', () => {
  const session = sessionFixture(FORCED_STROKES);
  const results = normalizeLogicalResults(session.results);
  const admission = admitFeatureCatalog(session.editorState.featureCatalog, results, { adopt: true, mode: 'circular' });
  const edits = savedResultEdits(session, 'circular', session.editorState, session.version);
  assert.equal(createSavedResultPlan(results, admission, edits).kind, 'MUTATING');
  assert.deepEqual(createSavedResultPlan(results, admission, { ...edits, sessionVersion: 46 }), createEmptySvgMutationPlan(1));
});

// D-19: a build of first-parent main drew a color rule's row under its key
// (`Complex IV [#ffe700]`) and wrote no `data-legend-name`
// (tests/fixtures/sessions/rule-captions-with-color-id.provenance.json). Load
// names such a row once, by the caption of the loaded rule whose row has that
// key, and draws the names; the keys, and the edits stored under them, stay.
const RULE_CAPTIONS = 'rule-captions-with-color-id.v46.gbdraw-session.json.gz';
const savedLegendRows = (svg) => svg.querySelectorAll('g[data-legend-key]').map((row) => [
  row.getAttribute('data-legend-key'), row.getAttribute('data-legend-name'), row.querySelector('text').textContent
]);
test('Load names the color-rule rows of a Result an older build saved, by the loaded rules (D-19)', () => {
  const session = sessionFixture(RULE_CAPTIONS);
  const results = normalizeLogicalResults(session.results);
  const admission = admitFeatureCatalog(session.editorState.featureCatalog, results, { adopt: true, mode: 'circular' });
  const edits = savedResultEdits(session, 'circular', session.editorState, session.version, session.modes.circular.config.rules);
  const svg = parseSvgText(results[0].content);
  const ruleRows = () => savedLegendRows(svg).filter(([key]) => key.startsWith('Complex'));
  const stroke = () => svg.querySelectorAll('g[data-legend-key]')
    .find((row) => row.getAttribute('data-legend-key') === 'Complex IV [#ffe700]').querySelector('path').getAttribute('stroke');
  // `parseSvgText` reads no text node: the saved texts are the keys (the Playwright case reads them).
  const keysAndNames = (rows) => rows.map(([key, name]) => [key, name]);
  assert.deepEqual(keysAndNames(ruleRows()), [['Complex IV [#ffe700]', null], ['Complex IV [#1f77b4]', null], ['Complex I', null]], 'as saved');
  const saved = keysAndNames(savedLegendRows(svg));
  const { callerTransforms } = createSavedResultPlan(results, admission, edits).operationsByResult[0];
  assert.deepEqual(callerTransforms.map((transform) => transform(svg)), [true]);
  assert.deepEqual(ruleRows(), [
    ['Complex IV [#ffe700]', 'Complex IV', 'Complex IV (1)'], ['Complex IV [#1f77b4]', 'Complex IV', 'Complex IV (2)'],
    ['Complex I', null, 'Complex I']
  ]);
  assert.deepEqual(keysAndNames(savedLegendRows(svg)).filter(([key]) => !key.startsWith('Complex IV')), saved.filter(([key]) => !key.startsWith('Complex IV')));
  assert.deepEqual(savedLegendRows(svg).filter(([key, , drawn]) => !key.startsWith('Complex IV') && drawn !== key), [], 'every other row is drawn by its key');
  assert.equal(stroke(), '#d62728', 'the stroke saved under the key stays on its row');
  assert.deepEqual(callerTransforms.map((transform) => transform(svg)), [false], 'a named Result is named once');
  assert.equal(reconcileMountedResult(svg, createEmptySvgMutationPlan(1).operationsByResult[0]).legendChanged, false, 'the naming pass keeps the names');
  // Without rules, and for a Result that shows no color ID, Load reads no Result.
  assert.deepEqual(createSavedResultPlan(results, admission, { ...edits, rules: [] }), createEmptySvgMutationPlan(1));
  assert.deepEqual(createSavedResultPlan([{ ...results[0], content: '<svg/>' }], admission, edits), createEmptySvgMutationPlan(1));
});

// The other allocated key: a rule named like a renderer row of another swatch
// (`tRNA` in red beside the tRNA type row). A row this build named, a row a
// rename shows under another key, and a row no loaded rule draws stay as saved.
test('Load names a saved rule row keyed for a renderer row of another swatch, and no other row (D-19)', () => {
  const ruleRow = (key, y, attributes = {}) => {
    const row = legendRow(key, y, attributes);
    row.querySelector('path').setAttribute('fill', '#ff0000');
    return row;
  };
  const savedFigure = (...rows) => {
    const svg = legendSvg(rows);
    // The block stroke of the figure: a rule row is drawn with it.
    svg.appendChild(new FakeElement('path', { id: 'f1', stroke: '#000000', 'stroke-width': '0.5' }));
    return svg;
  };
  const edits = {
    featureColorOverrides: {}, featureStrokeOverrides: {}, legendEntries: [], legendColorOverrides: {}, legendStrokeOverrides: {},
    originalLegendColors: {}, originalLegendOrder: ['tRNA', 'tRNA [#ff0000]', 'rRNA [#ff0000]'], mode: 'circular', sessionVersion: 46,
    rules: [{ feat: 'CDS', qual: 'gene', val: '^a$', color: '#ff0000', cap: 'tRNA' }]
  };
  const name = (svg) => {
    const results = [{ name: 'figure.svg', content: serializeNode(svg) }];
    const { callerTransforms } = createSavedResultPlan(results, null, edits).operationsByResult[0];
    return callerTransforms.map((transform) => transform(svg));
  };
  const older = savedFigure(legendRow('tRNA', 7), ruleRow('tRNA [#ff0000]', 31), ruleRow('rRNA [#ff0000]', 55));
  assert.deepEqual(name(older), [true]);
  assert.deepEqual(savedLegendRows(older), [
    ['tRNA', null, 'tRNA (1)'], ['tRNA [#ff0000]', 'tRNA', 'tRNA (2)'], ['rRNA [#ff0000]', null, 'rRNA [#ff0000]']
  ], 'a key no loaded rule is allocated stays a name');
  const named = savedFigure(legendRow('tRNA', 7), ruleRow('tRNA [#ff0000]', 31, { 'data-legend-name': 'Mine' }));
  const renamed = savedFigure(legendRow('tRNA', 7), ruleRow('tRNA [#ff0000]', 31, { 'data-gbdraw-base-data-legend-key': 'CDS' }));
  for (const svg of [named, renamed]) {
    const before = serializeNode(svg);
    assert.deepEqual(name(svg), [false]);
    assert.equal(serializeNode(svg), before);
  }
});

// F2: a caption saved with its own color ID whose rule is the only rule of
// that caption (the older build kept `Complex IV [#ffe700]` once the other
// color's rule was gone; here that rule is `Complex III`) is read as
// `Complex IV`, so its row's key changes from the saved text to the caption.
// Load moves the row of the saved figure and every Legend record saved under
// the old key to the new key, so the row shows its name and keeps its edits;
// a per-feature color edit saved with the old caption is read the same way.
// A caption that keeps several colors keeps its key.
const ONE_COLOR_CAPTION = (session) => JSON.parse(JSON.stringify(session).split('Complex IV [#1f77b4]').join('Complex III'));
test('Load moves a one-color rule row saved under its color ID, with its stored Legend edits, to the key of the read caption (F2)', () => {
  const full = sessionFixture(RULE_CAPTIONS);
  assert.equal(savedResultEdits(full, 'circular', full.editorState, full.version, full.modes.circular.config.rules).legendRekeys.size, 0,
    'a caption in two colors keeps its keys');
  const session = ONE_COLOR_CAPTION(full);
  const results = normalizeLogicalResults(session.results);
  const admission = admitFeatureCatalog(session.editorState.featureCatalog, results, { adopt: true, mode: 'circular' });
  const edits = savedResultEdits(session, 'circular', session.editorState, session.version, session.modes.circular.config.rules);
  assert.deepEqual([...edits.legendRekeys], [['Complex IV [#ffe700]', 'Complex IV']]);
  assert.deepEqual(edits.rules.map((rule) => rule.cap), ['Complex IV', 'Complex III', 'Complex I']);
  assert.deepEqual(Object.keys(edits.legendStrokeOverrides), ['Complex IV'], 'the stroke moves to the new key');
  assert.deepEqual(edits.originalLegendOrder.slice(1, 4), ['rRNA', 'Complex IV', 'Complex III'], 'in its place');
  assert.deepEqual(edits.legendEntries.slice(2, 3).map(({ caption, originalCaption }) => [caption, originalCaption]), [['Complex IV', 'Complex IV']]);
  assert.equal(edits.featureColorOverrides['record-1\u0000f48b7bf2f'].caption, 'Complex IV', 'item 3: the feature edit is read like its rule');
  const svg = parseSvgText(results[0].content);
  const { callerTransforms } = createSavedResultPlan(results, admission, edits).operationsByResult[0];
  callerTransforms.forEach((transform) => transform(svg));
  assert.deepEqual(savedLegendRows(svg).filter(([key]) => key.startsWith('Complex')), [
    ['Complex IV', null, 'Complex IV'], ['Complex III', null, 'Complex III'], ['Complex I', null, 'Complex I']
  ]);
  const row = svg.querySelectorAll('g[data-legend-key]').find((each) => each.getAttribute('data-legend-key') === 'Complex IV');
  assert.equal(row.querySelector('path').getAttribute('stroke'), '#d62728', 'the saved stroke stays on its row');
  assert.equal(reconcileMountedResult(svg, createEmptySvgMutationPlan(1).operationsByResult[0]).legendChanged, false);
});

// OV-288 (R15-4): the stroked repeat_region row has a Legend color of its
// own. Load finds the row's features by Python's row color, not the
// swatch's, so Reset returns them to Python's stroke.
test('Load records Python\'s stroke on the features of a stroked Legend row with a Legend color (OV-288)', () => {
  const [svg] = loadSavedResults(sessionFixture(FORCED_STROKES));
  assert.equal(partStrokes(svg).f9cf91913, '#e63946 3', 'the saved row stroke');
  assert.equal(partStrokes(resetPaint(svg)).f9cf91913, 'none 0');
});

// U2a review #2 (OV-195): a writer kept as a feature edit's
// `originalStroke*` the stroke the feature showed when its popup opened, so a
// feature stroked after its Legend row kept the row's stroke there. Load does
// not read Python's stroke from it.
test('Load does not take Python\'s stroke from a feature edit\'s originalStroke*', () => {
  const session = sessionFixture(FORCED_STROKES);
  Object.assign(Object.values(session.editorState.featureStrokes.overrides)[0], {
    originalStrokeColor: '#e63946', originalStrokeWidth: 3
  });
  const [svg] = loadSavedResults(session);
  assert.equal(partStrokes(resetPaint(svg)).f38ba7c3f, 'gray 2.0');
});

// U2a review #3 (OV-195): Python strokes a connector with the line stroke, a
// block with the block stroke of its record's size class. Load reads each
// kind from the parts of the Result that no stroke edit reached, not from one
// stroke the Session kept.
test('Load gives each kind of feature part the stroke Python drew for it in this Result', () => {
  const [svg] = loadSavedResults(sessionFixture(FORCED_STROKES), {
    patch: (root) => {
      root.appendChild(connectorPart('f38ba7c3f', '#2a9d8f', '2'));
      root.appendChild(connectorPart('f841fb8a8', 'lightgray', '5.0'));
    }
  });
  assert.equal(partStrokes(resetPaint(svg)).f38ba7c3f__line1, 'lightgray 5.0');
});

// A kind of part every one of which shows a stroke edit has no known Python
// stroke: its parts get no record, so Reset and Undo leave the saved stroke
// on them until Generate.
test('Load records no stroke Python is not known to have drawn', () => {
  const [svg] = loadSavedResults(sessionFixture(FORCED_STROKES), {
    patch: (root) => root.appendChild(connectorPart('f38ba7c3f', '#2a9d8f', '2'))
  });
  assert.equal(svg.getElementById('f38ba7c3f__line1').hasAttribute('data-gbdraw-base-stroke'), false);
  assert.equal(partStrokes(resetPaint(svg)).f38ba7c3f__line1, '#2a9d8f 2');
});

// Sessions 40-45 (main writes 44) keep one flat draft, whose edits the
// committed Result shows (tests/fixtures/sessions/two-mode-project.provenance.json:
// a Linear Legend row stroke, a row color, feature fills). Load records
// Python's paint on it the same way.
test('Load records Python\'s paint on a Session 44 Result', () => {
  const session = sessionFixture('two-mode-project.v44.gbdraw-session.json.gz');
  const edits = savedResultEdits(session, 'linear', session.editorState, session.version);
  assert.deepEqual(Object.keys(edits.legendStrokeOverrides), ['tRNA_RENAMED']);
  assert.deepEqual(edits.legendColorOverrides, { CDS: '#112233' });
  assert.equal(Object.keys(edits.featureColorOverrides).length, 2);
  const [svg] = loadSavedResults(session, { mode: 'linear' });
  assert.match(serializeNode(svg), /#333333/);
  Object.entries(partStrokes(resetPaint(svg))).forEach(([id, stroke]) => (
    assert.equal(stroke, id.includes('__line') ? 'lightgray 5.0' : 'gray 2.0', id)
  ));
  svg.querySelectorAll('g[data-legend-key]').forEach((row) => {
    const swatch = row.querySelector('path');
    const caption = row.getAttribute('data-legend-key');
    assert.equal(`${swatch.getAttribute('stroke')} ${swatch.getAttribute('stroke-width')}`, 'gray 2.0', caption);
    if (caption === 'CDS') assert.equal(swatch.getAttribute('fill'), '#54bcf8');
  });
});

// U2b: the live preview compiles the palette, the specific-color rules (their
// prepared matches), and Feature visibility with Generate's precedence, so every
// feature and Legend row fill the displayed Result shows is an operation.
const previewFeature = (admission) => admission.renderedFeaturesByResult[0].get('f0001');
const PREVIEW_DOMAINS = ['featureFills', 'featureVisibility', 'legendFills'];
// Python's rows of the displayed Result (`pythonLegendRows`): the CDS row as
// `buildSvgRoot` draws it; the features it draws (`displayedDrawnFeatures`).
const pythonRowsOf = (color) => new Map([['CDS', { key: 'CDS', color }]]);
const previewPlan = (admission, {
  paletteColors = { CDS: '#aaaaaa' }, drawnContext = null, domains = PREVIEW_DOMAINS, pythonRows = pythonRowsOf('#aaaaaa'), ...options
} = {}) => (
  compileDirectEditorMutationPlan({
    catalogAdmission: admission,
    legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#aaaaaa' }],
    originalLegendOrder: ['CDS'],
    livePreview: { drawnTracks: [], domains, paletteColors, drawnContext, pythonRows, features: [...admission.renderedFeaturesByResult[0].values()] },
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
  assert.deepEqual(fills({ paletteColors: { CDS: '#00ff00' }, manualSpecificRules: [rule] }), [null],
    'a feature whose rule match is not known keeps the fill the Result shows');
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
  // The rule draws the CDS row (Python drew the row in its color); a rule
  // whose caption names a row Python drew in another color draws its own row
  // (N-06, OV-291), whatever color the row's swatch shows (review H1).
  const rule = { feat: 'CDS', qual: 'gene', val: 'x', cap: 'CDS', color: '#ff0000' };
  const ruleRow = [{ caption: 'CDS', originalCaption: 'CDS', color: '#ff0000' }];
  reconcileMountedResult(mounted, previewPlan(admission, {
    manualSpecificRules: [rule], legendEntries: ruleRow, pythonRows: pythonRowsOf('#ff0000')
  }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#ff0000');
  reconcileMountedResult(mounted, previewPlan(admission, { manualSpecificRules: [rule] }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#aaaaaa', 'the row the caption names keeps its palette color');
  reconcileMountedResult(mounted, previewPlan(admission, { manualSpecificRules: [rule], legendEntries: ruleRow }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#aaaaaa', 'a swatch in the rule\'s color leaves the rule in its own row (review H1)');
  reconcileMountedResult(mounted, previewPlan(admission, { paletteColors: { CDS: '#00ff00' } }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#00ff00');
  reconcileMountedResult(mounted, previewPlan(admission), { domains: ['featureFills', 'legendFills'] });
  assert.equal(swatchFill(mounted), '#aaaaaa');
  assert.equal(featureFill(mounted), '#aaaaaa');
});

// OV-306: Python draws no row under a feature type's name once a captioned
// rule of that type colors a feature (`_generated_legend_rows`), so a rule
// captioned with its type's name that colors a feature draws that row in its
// color, whatever color Python drew the row in before the rule.
test('a rule captioned with its type that colors a feature draws the type row', () => {
  const { admission } = currentFixture();
  const mounted = buildSvgRoot();
  const rule = { feat: 'CDS', qual: 'gene', val: 'y', cap: 'CDS', color: '#ff0000' };
  recordRuleMatches([previewFeature(admission)], [ruleKey(rule)], () => ({ matched: [0], priorities: [0], declined: [] }));
  reconcileMountedResult(mounted, previewPlan(admission, { manualSpecificRules: [rule] }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#ff0000');
  reconcileMountedResult(mounted, previewPlan(admission, { manualSpecificRules: [rule], domains: ['legendFills'] }), { domains: ['legendFills'] });
  assert.equal(swatchFill(mounted), '#ff0000', 'also when the compile shows only the Legend fills');
});

// Review U2b #1, U3a A2a: a live rename is the executor's, which keeps
// Python's key for the row, so a palette or rule reconcile addressing that key
// keeps the row's Legend color or shows the new palette color, as Generate
// draws the renamed row.
const RENAMED_ROW = { legendEntries: [{ caption: 'Proteins', originalCaption: 'CDS', color: '#aaaaaa' }] };
const RENAME_DOMAINS = ['legendRenames', 'legendDeletes', 'legendAdds', 'legendOrder'];
test('a palette reconcile keeps the Legend color of a row renamed live', () => {
  const { admission } = currentFixture();
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, previewPlan(admission, { legendColorOverrides: { CDS: '#ff0000' } }), { domains: ['legendFills'] });
  reconcileMountedResult(mounted, previewPlan(admission, {
    ...RENAMED_ROW, legendColorOverrides: { CDS: '#ff0000' }, domains: RENAME_DOMAINS
  }), { domains: RENAME_DOMAINS });
  assert.equal(mounted.querySelector('g[data-legend-key]').getAttribute('data-legend-key'), 'Proteins');
  reconcileMountedResult(mounted, previewPlan(admission, {
    ...RENAMED_ROW, legendColorOverrides: { CDS: '#ff0000' }, paletteColors: { CDS: '#00ff00' }
  }), { domains: ['featureFills', 'legendFills'] });
  assert.equal(swatchFill(mounted), '#ff0000', 'the renamed row keeps its Legend color');
});
test('a palette reconcile shows the current palette color on a row renamed live', () => {
  const { admission } = currentFixture();
  const palette = buildSvgRoot();
  reconcileMountedResult(palette, previewPlan(admission, { paletteColors: { CDS: '#00ff00' } }), { domains: ['legendFills'] });
  reconcileMountedResult(palette, previewPlan(admission, { ...RENAMED_ROW, domains: RENAME_DOMAINS }), { domains: RENAME_DOMAINS });
  reconcileMountedResult(palette, previewPlan(admission, { ...RENAMED_ROW, paletteColors: { CDS: '#0000ff' } }), { domains: ['legendFills'] });
  assert.equal(swatchFill(palette), '#0000ff', 'the renamed row shows the current palette color');
});

// Review U2b #5: a rule whose match is not known (pending or declined) keeps
// the palette paint the Result shows, as the painter before U2b did.
test('a feature fill reconcile keeps the shown palette paint while a rule match is unknown', () => {
  const { admission } = currentFixture();
  const mounted = buildSvgRoot();
  reconcileMountedResult(mounted, previewPlan(admission, { paletteColors: { CDS: '#00ff00' } }), { domains: ['featureFills'] });
  assert.equal(featureFill(mounted), '#00ff00');
  const declined = { feat: 'CDS', qual: 'location', val: '1..9', cap: 'Declined', color: '#ff0000' };
  recordRuleMatches([previewFeature(admission)], [ruleKey(declined)], () => ({ matched: [], priorities: [], declined: [0] }));
  reconcileMountedResult(mounted, previewPlan(admission, { paletteColors: { CDS: '#00ff00' }, manualSpecificRules: [declined] }), { domains: ['featureFills'] });
  assert.equal(featureFill(mounted), '#00ff00', 'the palette paint stays');
});

// The work guard (allowlist): each edit kind compiles exactly the stages its
// shown domains need, read through the compile's structural metric. Any stage
// not listed fails, and a stroke edit reads no specific-color rule match. The
// domains are the ones app/app-setup.js shows for each edit kind
// (`LIVE_EDIT_DOMAINS`, and `editorPaintDomains` for a History step's changes).
const {
  strokes: STROKES, fills: FILLS, legendFills: LEGEND_FILLS, visibility: VISIBILITY, legendStructure: LEGEND_STRUCTURE,
  deletedRows: DELETED_ROWS
} = LIVE_EDIT_DOMAINS;
const EDIT_KIND_STAGES = [
  ['feature stroke (popup)', STROKES, ['strokes']],
  ['Legend row stroke', STROKES, ['strokes']],
  ['Legend row color (no rule)', LEGEND_FILLS, ['legendFills']],
  ['palette change', FILLS, ['fills', 'rules', 'legendFills']],
  ['color rule commit', FILLS, ['fills', 'rules', 'legendFills']],
  // A rule commit that adds or retires Legend rows shows them in its compile.
  ['color rule commit with Legend rows', [...FILLS, ...LEGEND_STRUCTURE], ['fills', 'rules', 'legendFills', 'legend']],
  ['Legend add or sort', LEGEND_STRUCTURE, ['legend']],
  // A deleted row's stroke leaves its features; a Restore returns it (OV-293).
  ['Legend delete', [...LEGEND_STRUCTURE, ...STROKES], ['legend', 'strokes']],
  ['Legend rename of a row without rules or features', [...LEGEND_STRUCTURE, ...LEGEND_FILLS, ...STROKES],
    ['legend', 'legendFills', 'strokes']],
  ['Legend Restore', DELETED_ROWS, ['legend', 'legendFills', 'strokes']],
  ['History step of a Legend delete or Restore', editorPaintDomains([{ path: ['editorState', 'legend', 'deletedEntries'] }]),
    ['legend', 'legendFills', 'strokes']],
  ['History step of a Legend row edit', editorPaintDomains([{ path: ['editorState', 'legend', 'entries'] }]),
    ['legend', 'legendFills']],
  ['Feature visibility edit', VISIBILITY, ['visibility']],
  ['History step of a Legend row stroke and color', editorPaintDomains([
    { path: ['editorState', 'legend', 'strokeOverrides', 'Proteins'] },
    { path: ['editorState', 'legend', 'colorOverrides'] }
  ]), ['legendFills', 'strokes']],
  ['History step of a feature stroke', editorPaintDomains([{ path: ['editorState', 'featureStrokes'] }]), ['strokes']],
  ['Result display, every paint domain changed', [...LEGEND_STRUCTURE, ...FILLS, ...VISIBILITY, ...STROKES],
    ['legend', 'fills', 'rules', 'legendFills', 'visibility', 'strokes']],
  ['Result display, nothing changed', LEGEND_STRUCTURE, ['legend']],
  ['Generate', null, ['fills', 'visibility', 'labels', 'legend', 'legendFills', 'strokes']]
];
test('each edit kind runs exactly the compile stages of its allowlist', () => {
  const { admission } = currentFixture();
  const feature = previewFeature(admission);
  const rule = { feat: 'CDS', qual: 'gene', val: 'x', cap: 'Rule row', color: '#ff0000' };
  recordRuleMatches([feature], [ruleKey(rule)], () => ({ matched: [], priorities: [], declined: [] }));
  // A feature's rule match results are read through its raw object
  // (services/rule-matchers.js), so a read of this feature's results counts.
  let ruleReads = 0;
  const vue = globalThis.window.Vue;
  const toRaw = vue.toRaw;
  vue.toRaw = (value) => { if (value === feature) ruleReads += 1; return toRaw(value); };
  const key = biologicalFeatureKey('record-a', 'feature-a');
  const context = { featureOverrides: {}, rules: [], selectedTypes: new Set(['CDS']), colorRules: [rule] };
  const metrics = [];
  const hooks = globalThis.__GBDRAW_TEST_HOOKS__;
  globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
  try {
    for (const [kind, domains, allowed] of EDIT_KIND_STAGES) {
      metrics.length = 0;
      ruleReads = 0;
      compileDirectEditorMutationPlan({
        catalogAdmission: admission,
        featureStrokeOverrides: { [key]: { strokeColor: '#e63946' } },
        featureOverrides: { [key]: { recordKey: 'record-a', biologicalFeatureId: 'feature-a', featureVisibility: 'off', labelText: 'x' } },
        legendEntries: [{ caption: 'Proteins', originalCaption: 'CDS', color: '#aaaaaa' }],
        originalLegendOrder: ['CDS'],
        legendColorOverrides: { Proteins: '#334455' },
        legendStrokeOverrides: { Proteins: { strokeColor: '#e63946' } },
        manualSpecificRules: [rule],
        livePreview: domains && { drawnTracks: [], domains, paletteColors: { CDS: '#00ff00' }, drawnContext: context }
      });
      assert.equal(ruleReads > 0, allowed.includes('rules'), `${kind}: a specific-color rule match is read only by the rules stage`);
      const compiles = metrics.filter(({ name }) => name === 'editorPlanCompile');
      assert.equal(compiles.length, 1, `${kind}: one compile`);
      assert.deepEqual(new Set(compiles[0].stages), new Set(allowed), `${kind}: the stages of its allowlist`);
    }
  } finally {
    globalThis.__GBDRAW_TEST_HOOKS__ = hooks;
    vue.toRaw = toRaw;
  }
});

// A Result without a feature catalog (a Session older than 40) is reached
// through the features read from it, so live edits still show.
test('a Result without a feature catalog receives the editor operations of its features', () => {
  const feature = { svg_id: 'f0001', type: 'CDS', id: 'f0001' };
  const addressing = displayedFeatureAddressing([feature], ['diagram.svg'], 0);
  const operations = compileDirectEditorMutationPlan({
    catalogAdmission: addressing,
    featureStrokeOverrides: { [featureOverrideKey(feature)]: { strokeColor: '#e63946', strokeWidth: 3 } },
    livePreview: { drawnTracks: [], domains: ['featureStrokes', 'featureFills'], paletteColors: { CDS: '#00ff00' }, drawnContext: null }
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
  const domains = ['featureFills', 'legendStrokes'];
  const operations = compileDirectEditorMutationPlan({
    catalogAdmission: displayedFeatureAddressing([feature], ['diagram.svg'], 0, mounted),
    legendEntries: [{ caption: 'CDS', color: '#aaaaaa' }],
    legendStrokeOverrides: { CDS: { strokeColor: '#e63946' } },
    livePreview: { drawnTracks: [], domains, paletteColors: { CDS: '#AAAAAA' }, drawnContext: null }
  }).operationsByResult[0];
  assert.deepEqual(operations.featureFills, []);
  reconcileMountedResult(mounted, operations, { domains });
  assert.equal(mounted.querySelector('[data-gbdraw-feature-id]').getAttribute('stroke'), '#e63946');
});

// OV-288 (R15-4): the executor reaches a Legend row's features by the color
// the renderer gives them, Python's row color (or the draft's where the
// Result predates it), never the color the row's swatch shows; a feature drawn
// in another color by this pass is not reached.
test('a Legend row stroke reaches the features in Python\'s row color, not in the swatch color', () => {
  const mounted = buildSvgRoot();
  const feature = mounted.querySelector('[data-gbdraw-feature-id]');
  const swatch = mounted.querySelector('g[data-legend-key]').querySelector('path');
  const domains = ['featureFills', 'legendFills', 'legendStrokes'];
  const operations = ({ draftColor = null, featureFills = [] } = {}) => ({
    ...createEmptySvgMutationPlan(1).operationsByResult[0],
    featureFills,
    legendFills: [{ caption: 'CDS', color: '#7b2cbf' }],
    legendStrokes: [{
      caption: 'CDS', strokeColor: '#e63946', strokeWidth: null,
      reach: { listedIds: [], namedIds: [], ownStrokeIds: [], draftColor }
    }]
  });
  reconcileMountedResult(mounted, operations(), { domains });
  assert.deepEqual([swatch.getAttribute('fill'), feature.getAttribute('stroke')], ['#7b2cbf', '#e63946']);
  reconcileMountedResult(mounted, operations({ draftColor: '#123456' }), { domains });
  assert.equal(feature.getAttribute('stroke'), null, 'a draft color the Result predates decides');
  reconcileMountedResult(mounted, operations({ featureFills: [{ renderedId: 'f0001', color: '#123456' }] }), { domains });
  assert.equal(feature.getAttribute('stroke'), null, 'a feature this pass draws in another color');
  reconcileMountedResult(mounted, operations({ draftColor: '#123456', featureFills: [{ renderedId: 'f0001', color: '#123456' }] }), { domains });
  assert.equal(feature.getAttribute('stroke'), '#e63946', 'the draft recolors the row and its feature alike');
});

// Review L3 of OV-288: one compile gives every Result of a batch the row's
// stroke, and the executor reaches in each Result the features drawn in the
// color Python gave the row there. Result 2 draws the row in another color
// (a rule's row of that caption): the same fills are reached in one Result
// and not in the other; a feature with its own stroke keeps it.
test('a Legend row stroke reaches in each Result of a batch the features of Python\'s row there', () => {
  const resultRoot = (rowColor, features) => {
    const root = new FakeElement('svg', { xmlns: 'http://www.w3.org/2000/svg' });
    features.forEach(([id, fill]) => root.appendChild(new FakeElement('path', {
      id, 'data-gbdraw-feature-id': id, 'data-gbdraw-feature-part': 'block', fill
    })));
    const legend = new FakeElement('g', { id: 'legend' });
    const featureLegend = new FakeElement('g', { id: 'feature_legend' });
    const row = new FakeElement('g', { 'data-legend-key': 'CDS' });
    row.appendChild(new FakeElement('path', { fill: rowColor }));
    featureLegend.appendChild(row);
    legend.appendChild(featureLegend);
    root.appendChild(legend);
    return root;
  };
  const roots = [
    resultRoot('#54bcf8', [['f0001', '#54bcf8'], ['f0003', '#54bcf8'], ['f0004', '#d3d3d3']]),
    resultRoot('#d3d3d3', [['f0002', '#d3d3d3'], ['f0005', '#54bcf8']])
  ];
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: {
      resultNames: ['record-a.svg', 'record-b.svg'],
      renderedTargetsByOverrideKey: new Map([['stable-a', [{ resultIndex: 0, renderedId: 'f0001' }]]]),
      resultIndexesByRenderedId: new Map([['f0001', new Set([0])]])
    },
    featureStrokeOverrides: { 'stable-a': { strokeColor: '#2a9d8f' } },
    legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#7b2cbf' }],
    originalLegendOrder: ['CDS'],
    legendStrokeOverrides: { CDS: { strokeColor: '#e63946', strokeWidth: 3 } }
  });
  const domains = ['featureStrokes', 'legendStrokes'];
  const stroked = roots.map((root, index) => {
    reconcileMountedResult(root, plan.operationsByResult[index], { domains });
    return root.querySelectorAll('[data-gbdraw-feature-id]').map((feature) => [feature.getAttribute('id'), feature.getAttribute('stroke')]);
  });
  assert.deepEqual(stroked, [
    [['f0001', '#2a9d8f'], ['f0003', '#e63946'], ['f0004', null]],
    [['f0002', '#e63946'], ['f0005', null]]
  ]);
});

// A catalog row that does not say which fill Python drew gives no record, so
// a later fill reconcile cannot remove the fill.
test('the Load normalizer records no fill Python is not known to have drawn', () => {
  const featureCatalog = catalog();
  featureCatalog.items[0].features[0].fillColor = '';
  const response = normalizeGenerationResponse({ results: [{ name: 'diagram.svg', content: '<svg/>' }], metadata: {}, featureCatalog });
  const admission = admitFeatureCatalog(featureCatalog, response.results, { mode: 'linear' });
  const svg = buildSvgRoot();
  svg.querySelector('[data-gbdraw-feature-id]').setAttribute('fill', '#112233');
  const plan = createSavedResultPlan(response.results, admission, {
    featureColorOverrides: { [biologicalFeatureKey('record-a', 'feature-a')]: '#112233' },
    featureStrokeOverrides: {}, legendEntries: [], legendColorOverrides: {}, legendStrokeOverrides: {},
    originalLegendColors: {}, originalLegendOrder: [], rules: [], mode: 'linear', sessionVersion: 44
  });
  plan.operationsByResult[0].callerTransforms.forEach((transform) => transform(svg));
  assert.doesNotMatch(serializeNode(svg), /data-gbdraw-base-fill/);
});

/** @param {Parameters<typeof reconcileMountedResult>} args */
const legendRowsChanged = (...args) => reconcileMountedResult(...args).legendChanged;

// U3a: the executor reconciles the Legend structure. A row of Python's keeps
// Python's key when renamed and stays hidden when deleted, so a reconcile
// without the edit returns the row as Python drew it.
/** @param {string} caption @param {number} y @param {Record<string, string>} [attributes] */
const legendRow = (caption, y, attributes = {}) => {
  const row = new FakeElement('g', { 'data-legend-key': caption, ...attributes });
  row.appendChild(new FakeElement('path', {
    fill: '#aaaaaa', stroke: attributes['data-test-stroke'] || '#000000', 'stroke-width': '0.5', transform: `translate(0, ${y})`
  }));
  const label = new FakeElement('text', { transform: `translate(22, ${y})` });
  label.textContent = caption;
  row.appendChild(label);
  return row;
};
const legendSvg = (rows = [legendRow('CDS', 7), legendRow('tRNA', 31), legendRow('rRNA', 55)]) => {
  const root = new FakeElement('svg', { xmlns: 'http://www.w3.org/2000/svg' });
  const legend = root.appendChild(new FakeElement('g', { id: 'legend' }));
  const featureLegend = legend.appendChild(new FakeElement('g', { id: 'feature_legend' }));
  rows.forEach((row) => featureLegend.appendChild(row));
  return root;
};
const legendOperations = (operations = {}) => ({ ...createEmptySvgMutationPlan(1).operationsByResult[0], ...operations });
const legendRows = (svg) => svg.querySelectorAll('g[data-legend-key]').map((row) => [
  row.getAttribute('data-legend-key'), row.querySelector('text').textContent, row.getAttribute('display')
]);

// D-28 (PD-OI-061 scenario revision 5): rows of one name and one swatch are
// one row. The naming pass hides the later of them and shows it again once its
// swatch differs; a delete of the merged row keeps it hidden, and a reconcile
// without the delete returns it to the merge.
test('rows of one name and one swatch are one row: the later row is hidden until it differs (D-28)', () => {
  const svg = legendSvg();
  const merged = () => svg.querySelectorAll('g[data-legend-key]').filter((row) => row.hasAttribute('data-legend-merged'))
    .map((row) => row.getAttribute('data-gbdraw-base-data-legend-key'));
  const renamed = legendOperations({ legendRenames: [{ from: 'tRNA', to: 'CDS', allowMissing: false }] });
  assert.deepEqual(reconcileMountedResult(svg, renamed), { changed: true, legendChanged: true });
  // A hidden row keeps the text it last showed.
  assert.deepEqual(legendRows(svg), [['CDS', 'CDS', null], ['CDS', 'tRNA', 'none'], ['rRNA', 'rRNA', null]]);
  assert.deepEqual(merged(), ['tRNA']);
  // D-35: the merged row records the key of the row it is drawn as, so an
  // edit of that row reaches both keys.
  assert.deepEqual(legendRowKeys(svg, 'CDS'), ['CDS', 'tRNA']);
  assert.deepEqual([legendRowKeys(svg, 'tRNA'), legendRowKeys(svg, 'rRNA')], [['tRNA'], ['rRNA']]);
  assert.deepEqual(reconcileMountedResult(svg, renamed), { changed: false, legendChanged: false }, 'merged once');
  const recolored = legendOperations({ ...renamed, legendFills: [{ caption: 'tRNA', color: '#ff0000' }] });
  reconcileMountedResult(svg, recolored);
  assert.deepEqual(legendRows(svg), [['CDS', 'CDS (1)', null], ['CDS', 'CDS (2)', null], ['rRNA', 'rRNA', null]], 'another swatch: two rows');
  assert.deepEqual(merged(), []);
  reconcileMountedResult(svg, renamed);
  assert.deepEqual(merged(), ['tRNA'], 'one swatch again: one row');
  assert.deepEqual(legendRowKeys(svg, 'CDS'), ['CDS', 'tRNA']);
  reconcileMountedResult(svg, legendOperations({ ...renamed, legendDeletes: [{ caption: 'tRNA', allowMissing: false }] }));
  assert.deepEqual([legendRows(svg)[1][2], merged()], ['none', []], 'a deleted merged row stays hidden as deleted');
  reconcileMountedResult(svg, renamed);
  assert.deepEqual([legendRows(svg)[1][2], merged()], ['none', ['tRNA']], 'without the delete it is merged again');
  reconcileMountedResult(svg, legendOperations());
  assert.deepEqual(legendRows(svg), [['CDS', 'CDS', null], ['tRNA', 'tRNA', null], ['rRNA', 'rRNA', null]], 'without the rename, Python\'s rows');
});

// U3a H1: main's writer renamed a generated row (GC content to "GC %") by
// rewriting its key without a record; Load records Python's key on that row
// once, so a live fill addressed by Python's key reaches it and a reconcile
// without the rename returns Python's row
// (tests/fixtures/sessions/forced-label-underlay-legend-rows.provenance.json).
test('Load records Python\'s key on a row a Session saved before U3a renamed', () => {
  const session = sessionFixture('forced-label-underlay-legend-rows.v44.gbdraw-session.json.gz');
  const keyRecord = 'data-gbdraw-base-data-legend-key';
  const recordsOf = (svg) => svg.querySelectorAll('g[data-legend-key]')
    .map((row) => [row.getAttribute('data-legend-key'), row.getAttribute(keyRecord)]);
  const [svg] = loadSavedResults(session);
  assert.deepEqual(recordsOf(svg), [
    ['CDS', null], ['GC %', 'GC content'], ['GC skew (+)', null], ['GC skew (-)', null], ['Added', null]
  ]);
  const fill = { caption: 'GC content', color: '#123456', allowMissing: true };
  reconcileMountedResult(svg, legendOperations({ legendRenames: [{ from: 'GC content', to: 'GC %' }], legendFills: [fill] }), {
    domains: [...LEGEND_STRUCTURE, 'legendFills']
  });
  const row = svg.querySelectorAll('g[data-legend-key]')[1];
  assert.equal(legendSwatch(row).getAttribute('fill'), '#123456', 'a fill addressed by Python\'s key reaches the row');
  reconcileMountedResult(svg, legendOperations({ legendAdds: [{ caption: 'Added', color: '#123456' }] }), { domains: LEGEND_STRUCTURE });
  assert.deepEqual(legendRows(svg)[1], ['GC content', 'GC content', null], 'without the rename the row is Python\'s');

  // Not provably a renamed row of Python's: an editor row, a second row of the
  // caption, or a Result that still has Python's row.
  const unmarked = (patch) => recordsOf(loadSavedResults(session, { patch })[0]).every(([, record]) => record === null);
  assert.ok(unmarked((saved) => saved.querySelectorAll('g[data-legend-key]')[1].setAttribute('data-legend-owner', 'direct-editor')));
  assert.ok(unmarked((saved) => {
    const rows = saved.querySelectorAll('g[data-legend-key]');
    rows[0].parentElement.appendChild(legendRow('GC %', 300));
  }));
  assert.ok(unmarked((saved) => {
    const rows = saved.querySelectorAll('g[data-legend-key]');
    rows[0].parentElement.appendChild(legendRow('GC content', 300));
  }));
});

// U3a review M1, L1: what a History restore reads of the restored Result.
test('a restored Result shows a relabeled row in place only where Python\'s order puts it', () => {
  const ruleRow = (caption, y) => legendRow(caption, y, { 'data-legend-owner': 'specific-color-file' });
  const placed = legendSvg([ruleRow('Alpha row', 7), legendRow('alpha', 7, { display: 'none' }), legendRow('tRNA', 31)]);
  assert.equal(legendRowTakesPlace(placed, 'alpha', 'Alpha row'), true, 'the commit\'s placement');
  const appended = legendSvg([legendRow('alpha', 7, { display: 'none' }), legendRow('tRNA', 31), ruleRow('Alpha row', 55)]);
  assert.equal(legendRowTakesPlace(appended, 'alpha', 'Alpha row'), false, 'appended before the rerender');
  // The Legend order the editor replays moves the hidden row last and records Python's order first.
  const reordered = (rows) => {
    const svg = legendSvg(rows);
    svg.querySelector('#feature_legend').setAttribute(LEGEND_ORDER_RECORD, JSON.stringify(['alpha', 'tRNA']));
    return svg;
  };
  assert.equal(legendRowTakesPlace(reordered([ruleRow('Alpha row', 7), legendRow('tRNA', 31), legendRow('alpha', 55, { display: 'none' })]),
    'alpha', 'Alpha row'), true, 'placed, the hidden row moved last');
  assert.equal(legendRowTakesPlace(reordered([legendRow('tRNA', 7), ruleRow('Alpha row', 31), legendRow('alpha', 55, { display: 'none' })]),
    'alpha', 'Alpha row'), false, 'appended, the hidden row moved last');
  assert.equal(legendRowTakesPlace(legendSvg([legendRow('Alpha row', 7), legendRow('tRNA', 31)]), 'alpha', 'Alpha row'), true,
    'Python drew the row');
  assert.equal(shownPythonLegendRow(placed, 'Alpha row'), false, 'a row a commit showed before Python drew it');
  assert.equal(shownPythonLegendRow(placed, 'alpha'), false, 'a hidden row');
  assert.equal(shownPythonLegendRow(placed, 'tRNA'), true);
});

test('a reconcile that hides a Legend row reports a Legend change once', () => {
  const svg = legendSvg();
  const deleted = legendOperations({ legendDeletes: [{ caption: 'tRNA' }] });
  assert.deepEqual(reconcileMountedResult(svg, deleted, { domains: LEGEND_STRUCTURE }), { changed: true, legendChanged: true });
  assert.deepEqual(reconcileMountedResult(svg, deleted, { domains: LEGEND_STRUCTURE }), { changed: false, legendChanged: false });
});

// EU review A-L2: removing the Legend order record is a change even when no
// row moves, so the stored Result loses the record with the mounted SVG.
test('a reconcile that removes the Legend order record reports a change though no row moves', () => {
  const svg = legendSvg();
  const group = svg.querySelector('#feature_legend');
  group.setAttribute(LEGEND_ORDER_RECORD, JSON.stringify(['CDS', 'tRNA', 'rRNA']));
  assert.equal(reconcileMountedResult(svg, legendOperations(), { domains: LEGEND_STRUCTURE }).changed, true);
  assert.equal(group.getAttribute(LEGEND_ORDER_RECORD), null);
  assert.deepEqual(legendRows(svg).map(([key]) => key), ['CDS', 'tRNA', 'rRNA']);
  assert.deepEqual(reconcileMountedResult(svg, legendOperations(), { domains: LEGEND_STRUCTURE }), { changed: false, legendChanged: false });
});

test('a deleted Legend row is hidden with Python\'s key and a reconcile without the delete shows it in its slot', () => {
  const svg = legendSvg();
  const drawn = serializeNode(svg);
  assert.equal(legendRowsChanged(svg, legendOperations({ legendDeletes: [{ caption: 'tRNA' }] }), { domains: LEGEND_STRUCTURE }), true);
  assert.deepEqual(legendRows(svg), [['CDS', 'CDS', null], ['tRNA', 'tRNA', 'none'], ['rRNA', 'rRNA', null]]);
  assert.equal(legendRowsChanged(svg, legendOperations(), { domains: LEGEND_STRUCTURE }), true);
  assert.equal(serializeNode(svg), drawn);
});

test('a fill addressed by Python\'s key reaches a renamed row; a reconcile without the rename restores Python\'s text', () => {
  const svg = legendSvg();
  const drawn = serializeNode(svg);
  const fill = { caption: 'CDS', color: '#123456', allowMissing: true };
  reconcileMountedResult(svg, legendOperations({ legendRenames: [{ from: 'CDS', to: 'Genes' }] }), { domains: LEGEND_STRUCTURE });
  reconcileMountedResult(svg, legendOperations({ legendRenames: [{ from: 'CDS', to: 'Genes' }], legendFills: [fill] }), {
    domains: [...LEGEND_STRUCTURE, 'legendFills']
  });
  assert.deepEqual(legendRows(svg)[0], ['Genes', 'Genes', null]);
  assert.equal(legendSwatch(svg.querySelectorAll('g[data-legend-key]')[0]).getAttribute('fill'), '#123456');
  assert.equal(legendRowsChanged(svg, legendOperations({ legendFills: [fill] }), {
    domains: [...LEGEND_STRUCTURE, 'legendFills']
  }), true);
  assert.deepEqual(legendRows(svg)[0], ['CDS', 'CDS', null]);
  assert.equal(legendSwatch(svg.querySelectorAll('g[data-legend-key]')[0]).getAttribute('fill'), '#123456');
  reconcileMountedResult(svg, legendOperations(), { domains: [...LEGEND_STRUCTURE, 'legendFills'] });
  assert.equal(serializeNode(svg), drawn);
});

// OV-1385: a Legend row's swatch is the path Python draws in the row, whatever
// fill it shows now. Python draws a row of color `none` with a swatch of fill
// `none`, so a live fill of `none` shows what Generate draws, and the next
// live edits still find the swatch. A row Python draws with a paint server
// (a gradient bar or a hatch pattern) has no solid swatch, as before.
test('a Legend swatch of fill none is still the row\'s swatch, and live fills after none reach it (OV-1385)', () => {
  const svg = legendSvg();
  const drawn = serializeNode(svg);
  const cds = () => svg.querySelectorAll('g[data-legend-key]')[0];
  const domains = ['legendFills', 'legendStrokes'];
  const fill = (color) => legendOperations({ legendFills: [{ caption: 'CDS', color }] });
  reconcileMountedResult(svg, fill('none'), { domains });
  assert.equal(legendSwatch(cds())?.getAttribute('fill'), 'none');
  reconcileMountedResult(svg, fill('#123456'), { domains });
  assert.equal(legendSwatch(cds()).getAttribute('fill'), '#123456', 'a live fill after none reaches the swatch');
  reconcileMountedResult(svg, fill('none'), { domains });
  reconcileMountedResult(svg, legendOperations({
    legendStrokes: [{
      caption: 'CDS', strokeColor: '#e63946', strokeWidth: 2,
      reach: { listedIds: [], namedIds: [], ownStrokeIds: [], draftColor: null }
    }]
  }), { domains });
  assert.equal(legendSwatch(cds()).getAttribute('stroke'), '#e63946', 'a live stroke after none reaches the swatch');
  reconcileMountedResult(svg, legendOperations(), { domains });
  assert.equal(serializeNode(svg), drawn, 'a reconcile without the edits returns Python\'s row');

  const pythonNone = legendRow('GC content', 7);
  pythonNone.children[0].setAttribute('fill', 'none');
  assert.equal(legendSwatch(pythonNone), pythonNone.children[0], 'Python\'s own swatch of fill none');
  const gradient = new FakeElement('g', { 'data-legend-key': 'Identity' }, [
    new FakeElement('linearGradient', { id: 'g1' }), new FakeElement('text'), new FakeElement('path', { fill: 'url(#g1)' })
  ]);
  const hatch = new FakeElement('g', { 'data-legend-key': 'Region' }, [
    new FakeElement('path', { fill: 'url(#hatch)' }), new FakeElement('text')
  ]);
  assert.deepEqual([legendSwatch(gradient), legendSwatch(hatch)], [null, null], 'a paint-server row has no solid swatch');
});

test('an editor row is cloned from Python\'s first row and removed by a reconcile without its add', () => {
  const editorRow = legendRow('Mine', 7, { 'data-legend-owner': 'direct-editor', 'data-test-stroke': '#ff00ff' });
  const svg = legendSvg([editorRow, legendRow('CDS', 31), legendRow('tRNA', 55)]);
  const add = { caption: 'New', color: '#556677', xPos: null, yPos: null };
  const mine = { caption: 'Mine', color: '#aaaaaa', xPos: null, yPos: null };
  // Python's first row, renamed and hidden, is copied as Python drew it.
  reconcileMountedResult(svg, legendOperations({
    legendRenames: [{ from: 'CDS', to: 'Genes' }], legendDeletes: [{ caption: 'CDS' }], legendAdds: [mine]
  }), { domains: LEGEND_STRUCTURE });
  assert.equal(legendRowsChanged(svg, legendOperations({
    legendRenames: [{ from: 'CDS', to: 'Genes' }], legendDeletes: [{ caption: 'CDS' }], legendAdds: [mine, add]
  }), { domains: LEGEND_STRUCTURE }), true);
  const added = svg.querySelectorAll('g[data-legend-key]').find((row) => row.getAttribute('data-legend-key') === 'New');
  assert.equal(serializeNode(added), serializeNode(legendRow('New', 31))
    .replace('<g data-legend-key="New">', '<g data-legend-key="New" data-legend-owner="direct-editor">')
    .replace('fill="#aaaaaa"', 'fill="#556677"'));
  assert.equal(legendRowsChanged(svg, legendOperations({
    legendRenames: [{ from: 'CDS', to: 'Genes' }], legendDeletes: [{ caption: 'CDS' }], legendAdds: [mine]
  }), { domains: LEGEND_STRUCTURE }), true);
  assert.deepEqual(legendRows(svg).map(([key]) => key), ['Mine', 'Genes', 'tRNA']);
});

test('a reconcile without an order returns Python\'s order', () => {
  const svg = legendSvg();
  const drawn = serializeNode(svg);
  // The order names rows by their keys (PR 3, D-25).
  const order = { captions: ['rRNA', 'CDS', 'tRNA'] };
  const renames = [{ from: 'CDS', to: 'Genes' }];
  assert.equal(legendRowsChanged(svg, legendOperations({ legendRenames: renames, legendOrder: [order] }), { domains: LEGEND_STRUCTURE }), true);
  assert.deepEqual(legendRows(svg).map(([key]) => key), ['rRNA', 'Genes', 'tRNA']);
  assert.equal(legendRowsChanged(svg, legendOperations({ legendRenames: renames }), { domains: LEGEND_STRUCTURE }), true);
  assert.deepEqual(legendRows(svg).map(([key]) => key), ['Genes', 'tRNA', 'rRNA']);
  reconcileMountedResult(svg, legendOperations(), { domains: LEGEND_STRUCTURE });
  assert.equal(serializeNode(svg), drawn);
});

test('the visibility domain does not un-hide a deleted row', () => {
  const svg = legendSvg();
  reconcileMountedResult(svg, legendOperations({ legendDeletes: [{ caption: 'tRNA' }] }), { domains: LEGEND_STRUCTURE });
  reconcileMountedResult(svg, legendOperations({ legendDeletes: [{ caption: 'tRNA' }] }), { domains: ['featureVisibility'] });
  reconcileMountedResult(svg, legendOperations(), { domains: ['featureVisibility', 'legendFills', 'legendStrokes'] });
  assert.deepEqual(legendRows(svg)[1], ['tRNA', 'tRNA', 'none']);
});

test('a second Legend structure reconcile changes nothing', () => {
  const svg = legendSvg();
  const operations = legendOperations({
    legendRenames: [{ from: 'CDS', to: 'Genes' }],
    legendDeletes: [{ caption: 'rRNA' }],
    legendAdds: [{ caption: 'New', color: '#556677', xPos: null, yPos: null }],
    legendOrder: [{ captions: ['tRNA', 'CDS', 'New'] }],
    legendFills: [{ caption: 'CDS', color: '#123456' }]
  });
  const domains = [...LEGEND_STRUCTURE, 'legendFills'];
  assert.equal(legendRowsChanged(svg, operations, { domains }), true);
  const once = serializeNode(svg);
  assert.equal(legendRowsChanged(svg, operations, { domains }), false);
  assert.equal(serializeNode(svg), once);
});

test('an export strips the Legend structure records and the rows a delete hid', () => {
  const svg = legendSvg();
  reconcileMountedResult(svg, legendOperations({
    legendRenames: [{ from: 'CDS', to: 'Genes' }],
    legendDeletes: [{ caption: 'rRNA' }],
    legendOrder: [{ captions: ['tRNA', 'CDS'] }]
  }), { domains: LEGEND_STRUCTURE });
  assert.match(serializeNode(svg), /data-gbdraw-base-legend-order=/);
  stripResultBaseAttributes(svg);
  assert.doesNotMatch(serializeNode(svg), /data-gbdraw-base-|display="none"/);
  assert.deepEqual(legendRows(svg), [['tRNA', 'tRNA', null], ['Genes', 'Genes', null]]);
});

// U3a A2a (gaps 1 and 2): a rule commit's row is added once where the Result
// lacks it, before the row it replaces, and owned by the rules; a retired row
// stays hidden without a record. A later reconcile keeps both until Python
// draws the rows again; a Result that draws the rule's row keeps Python's.
test('a rule commit\'s row is added where the Result lacks it and a retired row stays hidden', () => {
  const svg = legendSvg();
  const operations = legendOperations({
    legendAdds: [{ caption: 'Rule', color: '#123456', xPos: null, yPos: null, ifAbsent: true, before: 'tRNA' }],
    legendDeletes: [{ caption: 'tRNA', allowMissing: true, retire: true }]
  });
  assert.equal(legendRowsChanged(svg, operations, { domains: LEGEND_STRUCTURE }), true);
  const shown = [['CDS', 'CDS', null], ['Rule', 'Rule', null], ['tRNA', 'tRNA', 'none'], ['rRNA', 'rRNA', null]];
  assert.deepEqual(legendRows(svg), shown);
  const rule = svg.querySelectorAll('g[data-legend-key]')[1];
  assert.equal(rule.getAttribute('data-legend-owner'), 'specific-color-file');
  assert.equal(legendSwatch(rule).getAttribute('fill'), '#123456');
  assert.equal(legendRowsChanged(svg, operations, { domains: LEGEND_STRUCTURE }), false);
  assert.equal(legendRowsChanged(svg, legendOperations(), { domains: LEGEND_STRUCTURE }), false);
  assert.deepEqual(legendRows(svg), shown);
  const drawn = legendSvg([legendRow('CDS', 7), legendRow('Rule', 31)]);
  const bytes = serializeNode(drawn);
  assert.equal(legendRowsChanged(drawn, operations, { domains: LEGEND_STRUCTURE }), false);
  assert.equal(serializeNode(drawn), bytes);
});

// U3a A2a (gap 3): editor rows follow the order of their adds, as Generate
// appends them, so a renamed editor row keeps its place among them.
test('editor rows follow the order of their additions', () => {
  const svg = legendSvg();
  const add = (caption) => ({ caption, color: '#556677', xPos: null, yPos: null });
  reconcileMountedResult(svg, legendOperations({ legendAdds: [add('A'), add('B')] }), { domains: LEGEND_STRUCTURE });
  assert.equal(legendRowsChanged(svg, legendOperations({ legendAdds: [add('A2'), add('B')] }), { domains: LEGEND_STRUCTURE }), true);
  assert.deepEqual(legendRows(svg).map(([key]) => key), ['CDS', 'tRNA', 'rRNA', 'A2', 'B']);
});
