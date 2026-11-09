import assert from 'node:assert/strict';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import { createHash } from 'node:crypto';
import { mkdtemp, mkdir, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import { test } from 'node:test';
import { pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { readSavedSession } from './helpers/session-documents.mjs';

const repoRoot = new URL('../..', import.meta.url);
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-pairwise-payload-'));

const modulePaths = [
  'gbdraw/web/js/services/pairwise-match-popup.js',
  'gbdraw/web/js/services/feature-utils.js',
  'gbdraw/web/js/services/feature-sequence-fasta.js',
  // The popup names the reads of member and span sequences (D5); Python reads them when a sequence action is pressed.
  'gbdraw/web/js/services/catalog-sequences.js',
  // feature-catalog.js reads a catalog's drawn units as their records (`catalogUnitRecords`).
  'gbdraw/web/js/services/linear-units.js',
  'gbdraw/web/js/services/linear-comparisons.js',
  'gbdraw/web/js/services/linear-comparison-retained.js',
  'gbdraw/web/js/services/input-types.js',
  'gbdraw/web/js/services/record-draw-selection.js',
  'gbdraw/web/js/services/drawing/feature-edits.js',
  'gbdraw/web/js/services/drawing/paint.js',
  'gbdraw/web/js/services/drawing/slice-handle.js',
  'gbdraw/web/js/services/json-clone.js',
  'gbdraw/web/js/utils/reactivity.js',
  'gbdraw/web/js/utils/download-names.js',
  'gbdraw/web/js/services/record-source-coordinates.js',
  'gbdraw/web/js/utils/vue-raw.js',
  'gbdraw/web/js/services/losat-normalization.js',
  'gbdraw/web/js/contracts/request-shape.generated.js',
  'gbdraw/web/js/utils/contract-shape.js',
  'gbdraw/web/js/utils/optional-positive-number.js',
  'gbdraw/web/js/utils/error-normalization.js',
  // error-normalization.js reads Python's figure length units from the generated mode profiles.
  'gbdraw/web/js/mode-profiles.generated.js',
  'gbdraw/web/js/services/feature-catalog.js',
  'gbdraw/web/js/contracts/schema.generated.js',
  'gbdraw/web/js/services/feature-identity.js',
  'gbdraw/web/js/services/orthogroup-feature-metadata.js',
  'gbdraw/web/js/services/runtime-test-hooks.js',
  'gbdraw/web/js/utils/i18n.js'
];

await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');
for (const relativePath of modulePaths) {
  const target = join(tempRoot, relativePath);
  await mkdir(dirname(target), { recursive: true });
  await writeFile(target, await readFile(new URL(relativePath, repoRoot), 'utf8'), 'utf8');
}

const popupPath = join(tempRoot, 'gbdraw/web/js/services/pairwise-match-popup.js');
const instrumentedPopupPath = join(
  tempRoot,
  'gbdraw/web/js/services/pairwise-match-popup-instrumented.js'
);
const popupSource = await readFile(popupPath, 'utf8');
// Each instrumented specifier must be in the source: a replace that matches nothing would leave a call uncounted.
const replaceOnce = (source, from, to) => {
  assert.equal(source.split(from).length, 2, `instrumentation target not found once: ${from}`);
  return source.replace(from, to);
};
const instrumentedPopupSource = [
  [
    "import { formatFastaEntry } from './feature-sequence-fasta.js';",
    "import { formatFastaEntry as formatFastaEntryBase } from './feature-sequence-fasta.js';"
  ], [
    "import { featureSequenceReads } from './catalog-sequences.js';",
    "import { featureSequenceReads as featureSequenceReadsBase } from './catalog-sequences.js';"
  ], [
    '  featureIdentity,\n  identityMatches,\n  renderedFeatureIdentity,',
    '  featureIdentity as featureIdentityBase,\n  identityMatches,\n  renderedFeatureIdentity as renderedFeatureIdentityBase,'
  ], [
    'export const PAIRWISE_MATCH_SELECTOR = [',
    `const pairwiseMetrics = () => globalThis.__GBDRAW_PAIRWISE_TEST_METRICS__;
const formatFastaEntry = (...args) => {
  if (pairwiseMetrics()) pairwiseMetrics().fastaCalls += 1;
  return formatFastaEntryBase(...args);
};
const featureSequenceReads = (features) => {
  if (pairwiseMetrics()) pairwiseMetrics().sequenceReadNames += features.length;
  return featureSequenceReadsBase(features);
};
const featureIdentity = (source, options = {}) => {
  if (pairwiseMetrics()) pairwiseMetrics().identityCalls.push({
    source,
    mode: options.allowLegacySvgStable ? 'legacy' : 'standard'
  });
  return featureIdentityBase(source, options);
};
const renderedFeatureIdentity = (source) => {
  if (pairwiseMetrics()) pairwiseMetrics().identityCalls.push({ source, mode: 'rendered' });
  return renderedFeatureIdentityBase(source);
};

export const PAIRWISE_MATCH_SELECTOR = [`
  ], [
    'const createPairwisePayloadContext = ({',
    `const createPairwisePayloadContext = (options = {}) => {
  if (pairwiseMetrics()) pairwiseMetrics().contextCreations += 1;
  return createPairwisePayloadContextImplementation(options);
};

const createPairwisePayloadContextImplementation = ({`
  ]
].reduce((source, [from, to]) => replaceOnce(source, from, to), popupSource);
assert.notEqual(instrumentedPopupSource, popupSource, 'pairwise instrumentation was not applied');
await writeFile(instrumentedPopupPath, instrumentedPopupSource, 'utf8');

const popupModule = await import(pathToFileURL(join(
  tempRoot,
  'gbdraw/web/js/services/pairwise-match-popup.js'
)));
const instrumentedPopupModule = await import(pathToFileURL(instrumentedPopupPath));
const { admitFeatureCatalog } = await import(pathToFileURL(join(
  tempRoot,
  'gbdraw/web/js/services/feature-catalog.js'
)));
const { matchSequenceRecordsOf } = await import(pathToFileURL(join(
  tempRoot,
  'gbdraw/web/js/services/catalog-sequences.js'
)));

const elementFrom = (attributes) => ({
  style: { fill: attributes.fill || '' },
  getAttribute: (name) => attributes[name] || ''
});

// The rows a popup payload exposes: the feature rows of the row sections or of
// a collinear block's gene pairs, and the member rows of every group view.
const payloadFeatureRows = (payload) => [
  ...payload.sections.flatMap((entry) => entry.featureRows || []),
  ...(payload.block?.genePairs || []).flatMap((pair) => [pair.query, pair.subject].filter(Boolean))
];
const payloadMemberRows = (payload) => [
  ...(payload.group?.memberRows || []),
  ...(payload.blockOrthogroups || []).flatMap((group) => group.memberRows || [])
];

const sha256 = (value) => createHash('sha256')
  .update(JSON.stringify(value))
  .digest('hex');

const sortedCanonical = (value, featureTokens = new WeakMap()) => {
  if (value === null || typeof value !== 'object') return value;
  const token = featureTokens.get(value);
  if (token) return { __featureReference: token };
  if (Array.isArray(value)) {
    return value.map((entry) => sortedCanonical(entry, featureTokens));
  }
  return Object.fromEntries(
    Object.keys(value)
      .sort()
      .map((key) => [key, sortedCanonical(value[key], featureTokens)])
  );
};

const makeSyntheticState = () => {
  const renderedQuery = {
    fileIdx: 0,
    sourceFeatureIndex: 3,
    stable_feature_id: 'stable-query',
    svg_id: 'rendered-query',
    record_id: 'query-record',
    start: 9,
    end: 18,
    strand: '+',
    protein_id: 'QUERY.1',
    product: 'query enzyme'
  };
  const renderedSubject = {
    fileIdx: 1,
    sourceFeatureIndex: 7,
    stable_feature_id: 'stable-subject',
    svg_id: 'rendered-subject',
    record_id: 'subject-record',
    start: 29,
    end: 38,
    strand: '-',
    protein_id: 'SUBJECT.1',
    product: 'subject enzyme'
  };
  const sourceQuery = { ...renderedQuery, svg_id: 'stable-query' };
  const sourceSubject = { ...renderedSubject, svg_id: 'stable-subject' };
  const queryMember = {
    recordIndex: 0,
    featureIndex: 3,
    stableFeatureSvgId: 'stable-query',
    renderedFeatureSvgId: 'rendered-query',
    recordId: 'query-record',
    start: 9,
    end: 18,
    strand: '+',
    proteinId: 'QUERY.1',
    role: 'anchor',
    confidence: 'high',
    product: 'query enzyme'
  };
  const subjectMember = {
    recordIndex: 1,
    featureIndex: 7,
    stableFeatureSvgId: 'stable-subject',
    renderedFeatureSvgId: 'rendered-subject',
    recordId: 'subject-record',
    start: 29,
    end: 38,
    strand: '-',
    proteinId: 'SUBJECT.1',
    role: 'member',
    confidence: 'medium',
    product: 'subject enzyme'
  };
  const globalGroup = {
    id: 'og-main',
    name: 'Main enzymes',
    description: 'Synthetic similarity group',
    scope: 'global',
    member_count: 2,
    record_coverage_count: 2,
    members: [queryMember, subjectMember]
  };
  const localGroup = {
    ...globalGroup,
    name: 'Local enzymes',
    scope: 'cross_record',
    presentationScope: 'adjacent_local',
    collinearGroupScope: 'adjacent_local',
    groupKind: 'collinear_gene_group'
  };
  return {
    renderedQuery,
    renderedSubject,
    sourceQuery,
    sourceSubject,
    queryMember,
    subjectMember,
    globalGroup,
    localGroup,
    featureLookup: new Map([
      ['rendered-query', renderedQuery],
      ['rendered-subject', renderedSubject]
    ]),
    sourceFeatures: [sourceQuery, sourceSubject]
  };
};

const commonAttributes = {
  'data-gbdraw-match-id': 'match-synthetic',
  'data-query-record-id': 'query-record',
  'data-subject-record-id': 'subject-record',
  'data-query-record-index': '0',
  'data-subject-record-index': '1',
  'data-query-feature-index': '3',
  'data-subject-feature-index': '7',
  'data-query-feature-svg-id': 'rendered-query',
  'data-subject-feature-svg-id': 'rendered-subject',
  'data-query-stable-feature-svg-id': 'stable-query',
  'data-subject-stable-feature-svg-id': 'stable-subject',
  'data-query-protein-id': 'QUERY.1',
  'data-subject-protein-id': 'SUBJECT.1',
  'data-query-locus-id': 'query-locus',
  'data-subject-locus-id': 'subject-locus',
  'data-query-display-name': 'query gene',
  'data-subject-display-name': 'subject gene',
  'data-qstart': '10',
  'data-qend': '18',
  'data-sstart': '30',
  'data-send': '38',
  'data-identity': '88.5',
  'data-alignment-length': '9',
  'data-evalue': '1e-10',
  'data-bitscore': '42',
  fill: '#123456'
};

// The records the synthetic Result draws (`matchSequenceRecordsOf`): two drawn units, one per record ID.
const syntheticSequenceRecords = ({ comparisonSupplied = true } = {}) => ({
  unitKeys: ['query-unit', 'subject-unit'],
  unitKeysOfRecordId: (recordId) => ({ 'query-record': ['query-unit'], 'subject-record': ['subject-unit'] })[recordId] ?? [],
  comparisonSupplied: () => comparisonSupplied
});

const buildSyntheticCases = () => {
  const state = makeSyntheticState();
  const fullOptions = {
    featureLookup: state.featureLookup,
    sourceFeatures: state.sourceFeatures,
    orthogroups: [state.globalGroup, state.localGroup],
    matchSequenceRecords: syntheticSequenceRecords()
  };
  const duplicateRendered = {
    ...state.renderedQuery,
    svg_id: 'rendered-query'
  };
  return {
    ordinary: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'pairwise',
        'data-orthogroup-id': 'og-main'
      }),
      options: fullOptions,
      featureReferences: [state.renderedQuery, state.renderedSubject]
    },
    orthogroup: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'orthogroup',
        'data-orthogroup-id': 'og-main'
      }),
      options: fullOptions,
      featureReferences: [state.renderedQuery, state.renderedSubject]
    },
    collinear: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'collinear',
        'data-orthogroup-id': 'og-main',
        'data-collinearity-block-id': 'block-synthetic',
        'data-collinear-group-scope': 'adjacent_local',
        'data-collinearity-block-kind': 'cluster',
        'data-collinearity-orientation': 'plus',
        'data-collinearity-color-mode': 'orientation_identity',
        'data-collinearity-block-score': '42',
        'data-collinearity-anchor-index': '1',
        'data-collinearity-anchor-count': '1'
      }),
      options: fullOptions,
      featureReferences: [state.renderedQuery, state.renderedSubject]
    },
    homologyUnavailable: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'homology',
        'data-source-index': '4',
        'data-reference-side': 'query',
        'data-track-label': 'Comparison ring',
        'data-reference-record-id': 'query-record'
      }),
      // The comparison source supplied no FASTA: its span says so; the reference reads its drawn unit.
      options: { matchSequenceRecords: syntheticSequenceRecords({ comparisonSupplied: false }) },
      featureReferences: []
    },
    ambiguousFeature: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'pairwise'
      }),
      options: {
        featureLookup: new Map([
          ['rendered-query', state.renderedQuery],
          ['duplicate-key', duplicateRendered],
          ['rendered-subject', state.renderedSubject]
        ]),
        sourceFeatures: state.sourceFeatures
      },
      featureReferences: [state.renderedSubject]
    },
    conflictingEndpoint: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'pairwise',
        'data-query-stable-feature-svg-id': 'wrong-stable-id'
      }),
      options: fullOptions,
      featureReferences: [state.renderedSubject]
    },
    duplicateGroup: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'orthogroup',
        'data-orthogroup-id': 'og-duplicate'
      }),
      options: {
        featureLookup: new Map(),
        sourceFeatures: [],
        orthogroups: [
          { id: 'og-duplicate', scope: 'global', members: [] },
          { orthogroupId: 'og-duplicate', scope: 'global', members: [] }
        ]
      },
      featureReferences: []
    },
    scopedGroup: {
      element: elementFrom({
        ...commonAttributes,
        'data-match-kind': 'orthogroup',
        'data-orthogroup-id': 'og-main',
        'data-group-scope': 'adjacent_local'
      }),
      options: fullOptions,
      featureReferences: [state.renderedQuery, state.renderedSubject]
    }
  };
};

const SYNTHETIC_BASELINES_GZIP = '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';
const SYNTHETIC_BASELINES = SYNTHETIC_BASELINES_GZIP
  ? JSON.parse(gunzipSync(Buffer.from(SYNTHETIC_BASELINES_GZIP, 'base64')).toString('utf8'))
  : {};

test('compact pairwise payloads and hover rows match the accepted semantic baselines', () => {
  const captured = {};
  for (const [name, fixture] of Object.entries(buildSyntheticCases())) {
    const payload = popupModule.buildMatchPopupPayload(fixture.element, fixture.options);
    const featureTokens = new WeakMap(
      fixture.featureReferences.map((feature, index) => [feature, `${name}:${index}`])
    );
    const canonicalPayload = sortedCanonical(payload, featureTokens);
    const hoverRows = popupModule.buildPairwiseMatchHoverRows(payload);
    captured[name] = { payload: canonicalPayload, hoverRows };
    for (const reference of fixture.featureReferences) {
      const exposed = [...payloadFeatureRows(payload), ...payloadMemberRows(payload)]
        .some((row) => row.feature === reference);
      assert.equal(exposed, true, `${name} lost a canonical Feature object reference`);
    }
  }
  // The baselines are JSON: a member's FASTA header fields (`fastaFeature`) leave a field the feature lacks undefined.
  assert.deepStrictEqual(JSON.parse(JSON.stringify(captured)), SYNTHETIC_BASELINES);
});

const expectedHoverSummary = (payload) => ({
  id: payload.id,
  title: payload.title,
  subtitle: payload.subtitle,
  fill: payload.fill,
  rows: popupModule.buildPairwiseMatchHoverRows(payload)
});

test('lightweight hover summaries match accepted compact payload-derived rows', () => {
  const fixtures = buildSyntheticCases();
  for (const [name, fixture] of Object.entries(fixtures)) {
    if (name === 'orthogroup') continue;
    const payload = popupModule.buildMatchPopupPayload(fixture.element, fixture.options);
    assert.deepStrictEqual(
      popupModule.buildPairwiseMatchHoverSummary(fixture.element, fixture.options),
      expectedHoverSummary(payload),
      `${name} hover summary changed`
    );
  }

  const state = makeSyntheticState();
  const orthogroupElement = elementFrom({
    ...commonAttributes,
    'data-match-kind': 'orthogroup',
    'data-orthogroup-id': 'og-main',
    'data-group-scope': 'global'
  });
  const orthogroupOptions = {
    featureLookup: state.featureLookup,
    sourceFeatures: state.sourceFeatures,
    orthogroups: [state.globalGroup]
  };
  const orthogroupPayload = popupModule.buildMatchPopupPayload(
    orthogroupElement,
    orthogroupOptions
  );
  assert.deepStrictEqual(
    popupModule.buildPairwiseMatchHoverSummary(orthogroupElement, orthogroupOptions),
    expectedHoverSummary(orthogroupPayload)
  );
});

const createPairwiseMetrics = () => ({
  contextCreations: 0,
  fastaCalls: 0,
  sequenceReadNames: 0,
  identityCalls: []
});

const throwingCollection = (message) => new Proxy([], {
  get() {
    throw new Error(message);
  }
});

test('hover cannot touch full payload, Feature, source, member, FASTA, or sequence work', () => {
  const metrics = createPairwiseMetrics();
  globalThis.__GBDRAW_PAIRWISE_TEST_METRICS__ = metrics;
  const pairwiseElement = elementFrom({
    ...commonAttributes,
    'data-match-kind': 'pairwise',
    'data-orthogroup-id': 'og-main'
  });
  const pairwiseSummary = instrumentedPopupModule.buildPairwiseMatchHoverSummary(
    pairwiseElement,
    {
      featureLookup: throwingCollection('hover traversed rendered Features'),
      sourceFeatures: throwingCollection('hover traversed source Features'),
      orthogroups: () => {
        throw new Error('ordinary hover traversed orthogroups');
      }
    }
  );
  assert.deepStrictEqual(Object.keys(pairwiseSummary), ['id', 'title', 'subtitle', 'fill', 'rows']);

  const members = throwingCollection('orthogroup hover expanded members');
  const orthogroupElement = elementFrom({
    ...commonAttributes,
    'data-match-kind': 'orthogroup',
    'data-orthogroup-id': 'og-main',
    'data-group-scope': 'global'
  });
  const orthogroupSummary = instrumentedPopupModule.buildPairwiseMatchHoverSummary(
    orthogroupElement,
    {
      orthogroups: [{
        id: 'og-main',
        scope: 'global',
        name: 'Main enzymes',
        member_count: 2,
        members
      }]
    }
  );
  assert.equal(orthogroupSummary.rows.at(-1).value, '2');
  assert.equal(metrics.contextCreations, 0);
  assert.equal(metrics.fastaCalls, 0);
  assert.equal(metrics.sequenceReadNames, 0);
  delete globalThis.__GBDRAW_PAIRWISE_TEST_METRICS__;
});

const countedArray = (values, onTraversal) => new Proxy(values, {
  get(target, property, receiver) {
    if (property === Symbol.iterator) {
      return function iterator() {
        onTraversal();
        return target[Symbol.iterator]();
      };
    }
    return Reflect.get(target, property, receiver);
  }
});

// The popup opens without a base (D5): it names one read per member and per span, and Python reads them when a
// sequence action is pressed (app/sequence-reads.js), so opening it builds no FASTA.
test('one full payload request traverses each input once and builds no FASTA when it opens', () => {
  const state = makeSyntheticState();
  // The catalog rows carry the identity their sequences are read by (and, as every catalog row, no sequence).
  for (const [feature, recordKey, biologicalFeatureId] of [
    [state.renderedQuery, 'query-unit', 'query-bio'], [state.sourceQuery, 'query-unit', 'query-bio'],
    [state.renderedSubject, 'subject-unit', 'subject-bio'], [state.sourceSubject, 'subject-unit', 'subject-bio']
  ]) {
    Object.assign(feature, { recordKey, biologicalFeatureId });
  }
  const traversalCounts = {
    rendered: 0,
    source: 0,
    orthogroup: 0,
    members: 0
  };
  class CountedFeatureMap extends Map {
    values() {
      traversalCounts.rendered += 1;
      return super.values();
    }
  }
  const featureLookup = new CountedFeatureMap(state.featureLookup);
  const sourceFeatures = countedArray(state.sourceFeatures, () => {
    traversalCounts.source += 1;
  });
  const group = {
    ...state.globalGroup,
    members: countedArray([state.queryMember, state.subjectMember], () => {
      traversalCounts.members += 1;
    })
  };
  const orthogroups = countedArray([group], () => {
    traversalCounts.orthogroup += 1;
  });
  const element = elementFrom({
    ...commonAttributes,
    'data-match-kind': 'pairwise',
    'data-orthogroup-id': 'og-main',
    'data-group-scope': 'global'
  });
  const metrics = createPairwiseMetrics();
  globalThis.__GBDRAW_PAIRWISE_TEST_METRICS__ = metrics;
  const payload = instrumentedPopupModule.buildMatchPopupPayload(element, {
    featureLookup,
    sourceFeatures,
    orthogroups,
    matchSequenceRecords: syntheticSequenceRecords()
  });
  delete globalThis.__GBDRAW_PAIRWISE_TEST_METRICS__;

  assert.deepStrictEqual(traversalCounts, {
    rendered: 1,
    source: 1,
    orthogroup: 1,
    members: 1
  });
  assert.equal(metrics.contextCreations, 1);
  const memberRows = payload.group.memberRows;
  const resolvedMemberCount = memberRows.filter((row) => row.feature).length;
  assert.equal(resolvedMemberCount, 2);
  // One read named per resolved member, by its catalog identity; no FASTA text is built.
  assert.equal(metrics.sequenceReadNames, resolvedMemberCount);
  assert.equal(metrics.fastaCalls, 0);
  assert.deepStrictEqual(memberRows.map((row) => row.read), [
    { recordKey: 'query-unit', biologicalFeatureId: 'query-bio', bases: 9 },
    { recordKey: 'subject-unit', biologicalFeatureId: 'subject-bio', bases: 9 }
  ]);
  assert.equal(payload.group.memberSequenceCount, 2);
  assert.deepStrictEqual(payload.group.memberRead, { status: '', error: null, retry: null });
  // The spans name their reads by the drawn unit and the displayed coordinates.
  assert.deepStrictEqual(payload.sequenceBundle.entries.map((entry) => entry.read), [
    { kind: 'record', recordKey: 'query-unit', start: 10, end: 18 },
    { kind: 'record', recordKey: 'subject-unit', start: 30, end: 38 }
  ]);
  const payloadText = JSON.stringify(payload);
  assert.doesNotMatch(payloadText, /"(?:nt|aa|member(?:Nt|Aa))Fasta"|combinedFasta|"fasta"|_sequence"|Sequence":"/);

  const identityCounts = new Map();
  metrics.identityCalls.forEach(({ source, mode }) => {
    if (!identityCounts.has(source)) identityCounts.set(source, new Map());
    const modes = identityCounts.get(source);
    modes.set(mode, (modes.get(mode) || 0) + 1);
  });
  identityCounts.forEach((modes) => {
    modes.forEach((count, mode) => {
      assert.ok(count <= 1, `identity mode ${mode} normalized one object ${count} times`);
    });
  });
});

// OV-1010: a member's null or blank bound is missing, not 0.
test('group member coordinates read a missing bound as missing (OV-1010)', () => {
  for (const missing of [null, '', '  ']) {
    const state = makeSyntheticState();
    state.subjectMember.start = missing;
    const payload = popupModule.buildMatchPopupPayload(elementFrom({
      ...commonAttributes,
      'data-match-kind': 'orthogroup',
      'data-orthogroup-id': 'og-main'
    }), {
      featureLookup: state.featureLookup,
      sourceFeatures: state.sourceFeatures,
      orthogroups: [state.globalGroup]
    });
    assert.deepEqual(payload.group.memberRows.map((row) => row.coordinates), ['10..18 (+)', '38 (-)'], JSON.stringify(missing));
  }
});

test('OV-930: hover summaries of similarity bands read the groups once per groups array', () => {
  let traversals = 0;
  const group = (id, scope, name) => ({ id, scope, name, member_count: 2, members: [] });
  const orthogroups = countedArray([
    group('og-main', 'global', 'Main enzymes'),
    group('og-main', 'record_local', 'Local enzymes'),
    group('og-dup', 'global', 'First copy'),
    group('og-dup', 'global', 'Second copy'),
    group('og-other', 'global', 'Other enzymes')
  ], () => {
    traversals += 1;
  });
  const titleOf = (orthogroupId, groupScope) => instrumentedPopupModule.buildPairwiseMatchHoverSummary(
    elementFrom({
      ...commonAttributes,
      'data-match-kind': 'orthogroup',
      'data-orthogroup-id': orthogroupId,
      ...(groupScope ? { 'data-group-scope': groupScope } : {})
    }),
    { orthogroups }
  ).title;
  for (let index = 0; index < 20; index += 1) {
    assert.match(titleOf('og-main', 'global'), /Main enzymes/);
  }
  assert.match(titleOf('og-main', 'record_local'), /Local enzymes/);
  assert.match(titleOf('og-other', 'global'), /Other enzymes/);
  // A group is named only when exactly one group has that ID (and scope).
  assert.doesNotMatch(titleOf('og-dup', 'global'), /copy/);
  assert.doesNotMatch(titleOf('og-main', ''), /enzymes/);
  assert.equal(traversals, 1);
});

test('OV-911: popups on the same Result read its Features and source Features once', () => {
  const state = makeSyntheticState();
  const traversalCounts = { rendered: 0, source: 0 };
  class CountedFeatureMap extends Map {
    values() {
      traversalCounts.rendered += 1;
      return super.values();
    }
  }
  const options = {
    featureLookup: new CountedFeatureMap(state.featureLookup),
    sourceFeatures: countedArray(state.sourceFeatures, () => {
      traversalCounts.source += 1;
    }),
    orthogroups: [state.globalGroup, state.localGroup]
  };
  const element = elementFrom({
    ...commonAttributes,
    'data-match-kind': 'pairwise',
    'data-orthogroup-id': 'og-main',
    'data-group-scope': 'global'
  });
  const first = instrumentedPopupModule.buildMatchPopupPayload(element, options);
  for (let index = 0; index < 5; index += 1) {
    assert.deepStrictEqual(instrumentedPopupModule.buildMatchPopupPayload(element, options), first);
  }
  assert.deepStrictEqual(traversalCounts, { rendered: 1, source: 1 });
  // A new Result (a new Feature map) is read again.
  const next = instrumentedPopupModule.buildMatchPopupPayload(element, {
    ...options,
    featureLookup: new CountedFeatureMap(state.featureLookup)
  });
  assert.deepStrictEqual(next, first);
  assert.deepStrictEqual(traversalCounts, { rendered: 2, source: 2 });
});

test('source contracts keep heavy work out of hover and migrated inner loops', async () => {
  const svgActionsSource = await readFile(new URL(
    'gbdraw/web/js/app/feature-editor/svg-actions.js',
    repoRoot
  ), 'utf8');
  const activateStart = svgActionsSource.indexOf('const activateMatchHover =');
  const activateEnd = svgActionsSource.indexOf('handlerState.beginPreviewTransformInteraction', activateStart);
  const activateSource = svgActionsSource.slice(activateStart, activateEnd);
  assert.match(activateSource, /buildMatchHoverSummary\(matchEl\)/);
  assert.doesNotMatch(activateSource, /buildMatchPayload|ensureFeatureLookup/);

  const payloadStart = popupSource.indexOf('export const buildMatchPopupPayload =');
  const payloadEnd = popupSource.indexOf('const formatPairwiseMatchHoverRows =', payloadStart);
  const payloadSource = popupSource.slice(payloadStart, payloadEnd);
  assert.equal(
    (payloadSource.match(/createPairwisePayloadContext\s*\(/g) || []).length,
    1,
    'the full payload must create exactly one request context'
  );
  assert.doesNotMatch(popupSource, /memberFastaText\s*\(/);
  assert.doesNotMatch(popupSource, /globalThis|MutationObserver|requestIdleCallback/);

  for (const [startMarker, endMarker] of [
    ['const buildOrthogroupMemberRows =', 'const resolveFeatureSectionProteinIds ='],
    ['const buildBlockOrthogroups =', 'const buildFeatureRows =']
  ]) {
    const start = popupSource.indexOf(startMarker);
    const end = popupSource.indexOf(endMarker, start);
    const innerLoopSource = popupSource.slice(start, end);
    assert.doesNotMatch(
      innerLoopSource,
      /getFeatureLookupValues|sourceFeatures|orthogroups\.filter|featureLookup\.values/
    );
  }
});

const decodeXml = (value) => String(value || '')
  .replaceAll('&quot;', '"')
  .replaceAll('&apos;', "'")
  .replaceAll('&lt;', '<')
  .replaceAll('&gt;', '>')
  .replaceAll('&amp;', '&')
  .replace(/&#(\d+);/g, (_match, decimal) => String.fromCodePoint(Number(decimal)))
  .replace(/&#x([\da-f]+);/gi, (_match, hexadecimal) => String.fromCodePoint(Number.parseInt(hexadecimal, 16)));

const pairwiseElementsFromSvg = (svg) => Array.from(
  svg.matchAll(/<path\b(?=[^>]*data-gbdraw-match-id=)[^>]*>/g),
  (match) => {
    const attributes = {};
    for (const attribute of match[0].matchAll(/([\w:-]+)="([^"]*)"/g)) {
      attributes[attribute[1]] = decodeXml(attribute[2]);
    }
    return { attributes, element: elementFrom(attributes) };
  }
);

const exactFixturePath = new URL(
  'gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw',
  repoRoot
);
const exactSaved = await readSavedSession(new Blob([await readFile(exactFixturePath)]));
const exactSession = savedDrawingView(exactSaved, exactSaved.ui.activeDrawingId);
const exactAdmission = admitFeatureCatalog(
  exactSession.editorState.featureCatalog,
  exactSession.results,
  { mode: exactSession.renderRequest.mode }
);
const exactState = exactAdmission.featureState;
const exactFeatureLookup = new Map(
  exactState.extractedFeatures.map((feature) => [feature.svg_id, feature])
);
// The Gallery Session's Result and the request it was drawn from: its Linear spans read drawn units.
const exactSequenceRecords = matchSequenceRecordsOf(
  exactAdmission.catalog.items.find((item) => item.resultIndex === 0),
  { renderRequest: exactSession.renderRequest, resources: {} }
);
assert.ok(exactSequenceRecords, 'the exact fixture names no records for its spans');
const exactOptions = {
  featureLookup: exactFeatureLookup,
  sourceFeatures: exactState.biologicalFeatures,
  // The app passes the catalog groups, collinear presentation groups included
  // (svg-actions.js groupsForMatch); adding collinearGroups would list them twice.
  orthogroups: exactState.orthogroups,
  matchSequenceRecords: exactSequenceRecords
};
const exactFeatureTokens = new WeakMap([
  ...exactState.extractedFeatures.map((feature) => [feature, `rendered:${feature.svg_id}`]),
  ...exactState.biologicalFeatures.map((feature) => [
    feature,
    `source:${feature.recordKey}:${feature.biologicalFeatureId}`
  ])
]);
const exactElements = pairwiseElementsFromSvg(exactSession.results[0].content);
const exactFirst = exactElements[0];
const exactLater = exactElements.find((entry) => (
  entry.attributes['data-collinearity-block-id']
  && entry.attributes['data-collinearity-block-id']
    !== exactFirst.attributes['data-collinearity-block-id']
));

const EXACT_BASELINES = {
  first: {
    target: {
      matchId: 'comparison1_match1',
      blockId: 'block_0024',
      orthogroupIds: 'og_34;og_35;og_36;og_37;og_38;og_39;og_40;og_41;og_42;og_43;og_44;og_45;og_46;og_47;og_48'
    },
    // The Gallery Session uploads the pinned NZ_CP006932.1 revision (tutorial-data), whose WP_025208444.1 product
    // is "GGDEF domain-containing protein"; the earlier Session's copy read "DHH family phosphoesterase".
    // D5: the payload names the reads of the spans and the members' sequences, and holds no FASTA text.
    digest: 'b6e1f43ee2121aabca28c4cb9db4d979f42be5f4822ae524835bfda33fd7429b',
    hoverRows: [
      { id: 'Kind', label: 'Kind', value: 'Collinear block' },
      { id: 'Identity', label: 'Identity', value: '43.9684414414' },
      { id: 'Query', label: 'Query', value: '61879..78669' },
      { id: 'Subject', label: 'Subject', value: '19161..35565' },
      { id: 'Collinear groups', label: 'Collinear groups', value: '15' },
      { id: 'Block', label: 'Block', value: 'block_0024' }
    ],
    sectionTitles: [],
    blockOrthogroupCount: 15,
    memberRowCounts: Array(15).fill(5),
    featureRowCount: 30,
    exposedFeatureCount: 105,
    sequenceAvailability: [true, true],
    sequenceLengths: [16791,  16405],
    sequenceReads: [
      { kind: 'record', recordKey: 'linear-seq-2ccc54ef-3738-475f-9303-b61e2655f5b9', start: 61879, end: 78669 },
      { kind: 'record', recordKey: 'linear-seq-0a4f5060-b01d-4f8d-b5c8-c938c8083c7c', start: 19161, end: 35565 }
    ]
  },
  laterDistinctBlock: {
    target: {
      matchId: 'comparison1_match2',
      blockId: 'block_0155',
      orthogroupIds: 'og_274;og_275;og_276;og_277;og_278;og_279;og_280;og_281;og_282;og_283;og_284;og_285;og_286;og_287;og_288;og_289;og_290;og_291;og_292;og_293'
    },
    digest: '5460321d024aafb70dbc93db79e734587ab46bacd23ed58c426b17131ed71ba4',
    hoverRows: [
      { id: 'Kind', label: 'Kind', value: 'Collinear block' },
      { id: 'Identity', label: 'Identity', value: '64.1326369725' },
      { id: 'Query', label: 'Query', value: '484749..494272' },
      { id: 'Subject', label: 'Subject', value: '586163..595833' },
      { id: 'Collinear groups', label: 'Collinear groups', value: '20' },
      { id: 'Block', label: 'Block', value: 'block_0155' }
    ],
    sectionTitles: [],
    blockOrthogroupCount: 20,
    memberRowCounts: Array(20).fill(5),
    featureRowCount: 40,
    exposedFeatureCount: 140,
    sequenceAvailability: [true, true],
    sequenceLengths: [9524,  9671],
    sequenceReads: [
      { kind: 'record', recordKey: 'linear-seq-2ccc54ef-3738-475f-9303-b61e2655f5b9', start: 484749, end: 494272 },
      { kind: 'record', recordKey: 'linear-seq-0a4f5060-b01d-4f8d-b5c8-c938c8083c7c', start: 586163, end: 595833 }
    ]
  }
};

test('exact large-fixture payloads retain canonical digests, structure, and Feature references', () => {
  assert.ok(exactFirst, 'the exact fixture has no first pairwise target');
  assert.ok(exactLater, 'the exact fixture has no later distinct-block target');
  const captured = {};
  for (const [name, target] of [['first', exactFirst], ['laterDistinctBlock', exactLater]]) {
    const payload = popupModule.buildMatchPopupPayload(target.element, exactOptions);
    assert.deepStrictEqual(
      popupModule.buildPairwiseMatchHoverSummary(target.element, exactOptions),
      expectedHoverSummary(payload),
      `${name} exact-fixture hover summary changed`
    );
    const canonical = sortedCanonical(payload, exactFeatureTokens);
    const featureRows = payloadFeatureRows(payload);
    const memberRows = payloadMemberRows(payload);
    const exposedFeatures = [...featureRows, ...memberRows]
      .map((row) => row.feature)
      .filter(Boolean);
    exposedFeatures.forEach((feature) => {
      assert.ok(exactFeatureTokens.has(feature), `${name} exposed a non-canonical Feature object`);
    });
    captured[name] = {
      target: {
        matchId: target.attributes['data-gbdraw-match-id'],
        blockId: target.attributes['data-collinearity-block-id'],
        orthogroupIds: target.attributes['data-orthogroup-id']
      },
      digest: sha256(canonical),
      hoverRows: popupModule.buildPairwiseMatchHoverRows(payload),
      sectionTitles: payload.sections.map((entry) => entry.title),
      blockOrthogroupCount: payload.blockOrthogroupCount,
      memberRowCounts: payload.blockOrthogroups.map((group) => group.memberRows.length),
      featureRowCount: featureRows.length,
      exposedFeatureCount: exposedFeatures.length,
      sequenceAvailability: payload.sequenceBundle.entries.map((entry) => entry.available),
      sequenceLengths: payload.sequenceBundle.entries.map((entry) => entry.sequenceLength),
      sequenceReads: payload.sequenceBundle.entries.map((entry) => entry.read)
    };
  }
  assert.deepStrictEqual(captured, EXACT_BASELINES);
});
