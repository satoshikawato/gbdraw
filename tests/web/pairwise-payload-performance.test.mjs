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
  // The popup reads the text of each match span reason from the generated table.
  'gbdraw/web/js/services/drawn-attributes.generated.js',
  'gbdraw/web/js/services/feature-utils.js',
  'gbdraw/web/js/services/feature-sequence-fasta.js',
  // The popup names the reads of member and span sequences (D5); Python reads them when a sequence action is pressed.
  'gbdraw/web/js/services/catalog-sequences.js',
  // catalog-sequences.js keys a feature's reply by the draft owner's identity key.
  'gbdraw/web/js/services/feature-placement.js',
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
    "import { formatFastaEntry, memberFastaFeature } from './feature-sequence-fasta.js';",
    "import { formatFastaEntry as formatFastaEntryBase, memberFastaFeature } from './feature-sequence-fasta.js';"
  ], [
    "import { featureSequenceReads, unreadSequence } from './catalog-sequences.js';",
    "import { featureSequenceReads as featureSequenceReadsBase, unreadSequence } from './catalog-sequences.js';"
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
    { kind: 'record', recordKey: 'query-unit', recordId: 'query-record', start: 10, end: 18 },
    { kind: 'record', recordKey: 'subject-unit', recordId: 'subject-record', start: 30, end: 38 }
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

// LOW-10 (D5 review r1): a Circular comparison band carries no record index; it names its reference by the record ID
// of its hit, which may be an alias (name, accession) or a record without feature rows. The popup sends that ID, and
// Python finds the drawn record by ID, then alias, as the match registry did before D5.
test('a Circular comparison band reads its reference by the record ID its hit names', () => {
  const attributes = {
    ...commonAttributes,
    'data-match-kind': 'homology',
    'data-source-index': '0',
    'data-reference-side': 'query',
    'data-query-record-id': 'NC_000001'
  };
  delete attributes['data-query-record-index'];
  delete attributes['data-subject-record-index'];
  const records = { unitKeys: ['unit-a', 'unit-b'], comparisonSupplied: () => true };
  const payload = popupModule.buildMatchPopupPayload(elementFrom(attributes), { matchSequenceRecords: records });
  assert.deepStrictEqual(payload.sequenceBundle.entries.map((entry) => [entry.read, entry.unavailableReason]), [
    [{ kind: 'reference', recordId: 'NC_000001', recordKey: null, start: 10, end: 18 }, ''],
    [{ kind: 'comparison', sourceIndex: 0, recordId: 'subject-record', start: 30, end: 38 }, '']
  ]);
});

// Every coordinate is Python's to judge after the record, as the popup judged it before D5 (the record, then "whole",
// then "greater than or equal to 1"): a coordinate that is not a whole number goes as null, a whole one as itself.
test('a coordinate below 1 or not whole is read by Python after the record, in the popup\'s old order', () => {
  const spansOf = (qstart) => popupModule.buildMatchPopupPayload(elementFrom({ ...commonAttributes, 'data-qstart': qstart }), {
    matchSequenceRecords: syntheticSequenceRecords()
  }).sequenceBundle.entries[0];
  for (const [qstart, start] of [['0', 0], ['-5', -5], ['1.5', null], ['ten', null], ['1e21', 1e21]]) {
    const entry = spansOf(qstart);
    assert.deepStrictEqual([entry.read, entry.unavailableReason], [{ kind: 'record', recordKey: 'query-unit', recordId: 'query-record', start, end: 18 }, ''], qstart);
    // The popup row shows no length for a span whose coordinates are not whole (index.html hides a 0 length).
    if (start === null) assert.equal(entry.sequenceLength, 0, qstart);
  }
  // A comparison source that supplied no FASTA is said first, whatever its coordinates (base resolved the source
  // before it read the coordinates).
  const homology = popupModule.buildMatchPopupPayload(elementFrom({
    ...commonAttributes, 'data-match-kind': 'homology', 'data-source-index': '4', 'data-reference-side': 'query',
    'data-sstart': '2.5'
  }), { matchSequenceRecords: syntheticSequenceRecords({ comparisonSupplied: false }) }).sequenceBundle.entries[1];
  assert.deepStrictEqual([homology.read, homology.unavailableReason], [null, 'Comparison sequence was not supplied for this BLAST source.']);
});

// D5 review r3: what the popup knows at open, as base c75daa1de showed it (r4fix/probe-base.out,
// probe-base-reference.out): a Linear record index without a drawn unit is an unavailable source; a Circular reference
// whose record index names no unit does not find the record it names, and one that names a unit is searched there
// only; a comparison band without a source index searches every source; the row shows the drawing's own text.
test('the popup finds a span\'s record and shows its row as before D5', () => {
  const entriesOf = (attributes, records = syntheticSequenceRecords()) => {
    // A stub element reads an undefined attribute as missing.
    const element = elementFrom({ ...commonAttributes, ...attributes });
    return popupModule.buildMatchPopupPayload(element, { matchSequenceRecords: records }).sequenceBundle.entries;
  };
  const [missingUnit] = entriesOf({ 'data-query-record-index': '7' });
  assert.deepStrictEqual([missingUnit.read, missingUnit.unavailableReason], [null, 'The sequence source for this match is unavailable.']);
  const homology = { 'data-match-kind': 'homology', 'data-reference-side': 'query', 'data-source-index': '0' };
  const [unitless] = entriesOf({ ...homology, 'data-query-record-index': '7' });
  assert.deepStrictEqual([unitless.read, unitless.unavailableReason], [null, 'No sequence source matched record “query-record”.']);
  const [reference] = entriesOf(homology);
  assert.deepStrictEqual(reference.read, { kind: 'reference', recordId: 'query-record', recordKey: 'query-unit', start: 10, end: 18 });
  const sourceless = { ...homology, 'data-source-index': undefined };
  assert.deepStrictEqual(entriesOf(sourceless)[1].read,
    { kind: 'comparison', sourceIndex: null, recordId: 'subject-record', start: 30, end: 38 });
  const [, unsupplied] = entriesOf(sourceless, syntheticSequenceRecords({ comparisonSupplied: false }));
  assert.deepStrictEqual([unsupplied.read, unsupplied.unavailableReason], [null, 'Comparison sequence was not supplied for this BLAST source.']);
  for (const [qstart, row] of [['abc', 'abc..18'], ['', '..18'], ['1e1', '1e1..18']]) {
    const [entry] = entriesOf({ 'data-qstart': qstart });
    assert.equal(`${entry.span.start}..${entry.span.end}`, row, qstart);
    assert.ok(!/NaN|_0-/.test(entry.filename), entry.filename);
  }
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
    digest: 'ced4e768e946a7c4b144b0e491aaf2bd760513dab3123bd19f8dac6def3e1919',
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
      { kind: 'record', recordKey: 'linear-seq-2ccc54ef-3738-475f-9303-b61e2655f5b9', recordId: 'AP027078.1', start: 61879, end: 78669 },
      { kind: 'record', recordKey: 'linear-seq-0a4f5060-b01d-4f8d-b5c8-c938c8083c7c', recordId: 'AP027131.1', start: 19161, end: 35565 }
    ]
  },
  laterDistinctBlock: {
    target: {
      matchId: 'comparison1_match2',
      blockId: 'block_0155',
      orthogroupIds: 'og_274;og_275;og_276;og_277;og_278;og_279;og_280;og_281;og_282;og_283;og_284;og_285;og_286;og_287;og_288;og_289;og_290;og_291;og_292;og_293'
    },
    digest: '8ac71b136226e4e16bdf84059dfcf98f10114cd2361cc6d59dea55c42ad42436',
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
      { kind: 'record', recordKey: 'linear-seq-2ccc54ef-3738-475f-9303-b61e2655f5b9', recordId: 'AP027078.1', start: 484749, end: 494272 },
      { kind: 'record', recordKey: 'linear-seq-0a4f5060-b01d-4f8d-b5c8-c938c8083c7c', recordId: 'AP027131.1', start: 586163, end: 595833 }
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
