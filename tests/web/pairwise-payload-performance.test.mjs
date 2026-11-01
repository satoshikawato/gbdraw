import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { mkdtemp, mkdir, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import { test } from 'node:test';
import { pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';

const repoRoot = new URL('../..', import.meta.url);
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-pairwise-payload-'));

const modulePaths = [
  'gbdraw/web/js/app/pairwise-match-popup.js',
  'gbdraw/web/js/services/feature-utils.js',
  'gbdraw/web/js/services/feature-sequence-fasta.js',
  'gbdraw/web/js/services/match-sequences.js',
  'gbdraw/web/js/services/record-options.js',
  'gbdraw/web/js/services/record-draw-selection.js',
  'gbdraw/web/js/services/drawing/feature-edits.js',
  'gbdraw/web/js/services/drawing/paint.js',
  'gbdraw/web/js/services/drawing/slice-handle.js',
  'gbdraw/web/js/services/json-clone.js',
  'gbdraw/web/js/utils/reactivity.js',
  'gbdraw/web/js/utils/download-names.js',
  'gbdraw/web/js/services/record-source-coordinates.js',
  'gbdraw/web/js/services/genbank-header.js',
  'gbdraw/web/js/services/conservation-series.js',
  'gbdraw/web/js/utils/conservation-colors.generated.js',
  'gbdraw/web/js/utils/color-utils.js',
  'gbdraw/web/js/utils/named-colors.js',
  'gbdraw/web/js/utils/vue-raw.js',
  'gbdraw/web/js/services/losat-normalization.js',
  'gbdraw/web/js/contracts/request-shape.generated.js',
  'gbdraw/web/js/utils/contract-shape.js',
  'gbdraw/web/js/utils/field-reading.js',
  'gbdraw/web/js/utils/optional-positive-number.js',
  'gbdraw/web/js/utils/error-normalization.js',
  'gbdraw/web/js/services/file-content-cache.js',
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
  const source = relativePath.endsWith('/file-content-cache.js')
    ? 'export const readFileText = (file) => file.text();\n'
    : await readFile(new URL(relativePath, repoRoot), 'utf8');
  await writeFile(target, source, 'utf8');
}
// match-sequences.js reads restored Session files through record-discovery.js, which this test does not use.
await writeFile(join(tempRoot, 'gbdraw/web/js/services/record-discovery.js'), [
  "const unused = () => { throw new Error('record-discovery.js is not part of this test'); };",
  'export const gffFastaReaderSelector = unused;',
  'export const readAllSequenceRecords = unused;',
  'export const readComparisonSequenceRecords = unused;',
  ''
].join('\n'), 'utf8');

const popupPath = join(tempRoot, 'gbdraw/web/js/app/pairwise-match-popup.js');
const instrumentedPopupPath = join(
  tempRoot,
  'gbdraw/web/js/app/pairwise-match-popup-instrumented.js'
);
const popupSource = await readFile(popupPath, 'utf8');
const instrumentedPopupSource = popupSource
  .replace(
    "import { buildFeatureSequenceFastas } from '../services/feature-sequence-fasta.js';",
    "import { buildFeatureSequenceFastas as buildFeatureSequenceFastasBase } from '../services/feature-sequence-fasta.js';"
  )
  .replace(
    "import { buildMatchSequenceBundle } from '../services/match-sequences.js';",
    "import { buildMatchSequenceBundle as buildMatchSequenceBundleBase } from '../services/match-sequences.js';"
  )
  .replace(
    '  featureIdentity,\n  identityMatches,\n  renderedFeatureIdentity,',
    '  featureIdentity as featureIdentityBase,\n  identityMatches,\n  renderedFeatureIdentity as renderedFeatureIdentityBase,'
  )
  .replace(
    'export const PAIRWISE_MATCH_SELECTOR = [',
    `const pairwiseMetrics = () => globalThis.__GBDRAW_PAIRWISE_TEST_METRICS__;
const buildFeatureSequenceFastas = (...args) => {
  if (pairwiseMetrics()) pairwiseMetrics().fastaCalls += 1;
  return buildFeatureSequenceFastasBase(...args);
};
const buildMatchSequenceBundle = (...args) => {
  if (pairwiseMetrics()) pairwiseMetrics().sequenceBundleCalls += 1;
  return buildMatchSequenceBundleBase(...args);
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
  )
  .replace(
    'const createPairwisePayloadContext = ({',
    `const createPairwisePayloadContext = (options = {}) => {
  if (pairwiseMetrics()) pairwiseMetrics().contextCreations += 1;
  return createPairwisePayloadContextImplementation(options);
};

const createPairwisePayloadContextImplementation = ({`
  );
assert.notEqual(instrumentedPopupSource, popupSource, 'pairwise instrumentation was not applied');
await writeFile(instrumentedPopupPath, instrumentedPopupSource, 'utf8');

const popupModule = await import(pathToFileURL(join(
  tempRoot,
  'gbdraw/web/js/app/pairwise-match-popup.js'
)));
const instrumentedPopupModule = await import(pathToFileURL(instrumentedPopupPath));
const { admitFeatureCatalog } = await import(pathToFileURL(join(
  tempRoot,
  'gbdraw/web/js/services/feature-catalog.js'
)));
const { createSequenceSourceRegistry } = await import(pathToFileURL(join(
  tempRoot,
  'gbdraw/web/js/services/match-sequences.js'
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
    product: 'query enzyme',
    nucleotide_sequence: 'ATGAAATAA',
    amino_acid_sequence: 'MK'
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
    product: 'query enzyme',
    nucleotideSequence: 'ATGAAATAA',
    aminoAcidSequence: 'MK'
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

const buildSyntheticCases = () => {
  const state = makeSyntheticState();
  const sources = createSequenceSourceRegistry([
    {
      key: 'linear:record:0',
      recordId: 'query-record',
      aliases: [],
      sequence: 'AAAAAAAAAATGAAATAACCCCC',
      origin: 'linear-record',
      recordIndex: 0
    },
    {
      key: 'linear:record:1',
      recordId: 'subject-record',
      aliases: [],
      sequence: 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCATGAAATAAGGGG',
      origin: 'linear-record',
      recordIndex: 1
    }
  ]);
  const fullOptions = {
    featureLookup: state.featureLookup,
    sourceFeatures: state.sourceFeatures,
    orthogroups: [state.globalGroup, state.localGroup],
    resolveSequenceSource: sources.resolve
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
      options: { resolveSequenceSource: () => null },
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
  assert.deepStrictEqual(captured, SYNTHETIC_BASELINES);
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
  sequenceBundleCalls: 0,
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
  assert.equal(metrics.sequenceBundleCalls, 0);
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

test('one full payload request traverses each input once and builds one FASTA pair per member', () => {
  const state = makeSyntheticState();
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
  const sequenceRegistry = createSequenceSourceRegistry([
    {
      key: 'linear:record:0',
      recordId: 'query-record',
      aliases: [],
      sequence: 'AAAAAAAAAATGAAATAACCCCC',
      origin: 'linear-record',
      recordIndex: 0
    },
    {
      key: 'linear:record:1',
      recordId: 'subject-record',
      aliases: [],
      sequence: 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCATGAAATAAGGGG',
      origin: 'linear-record',
      recordIndex: 1
    }
  ]);
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
    resolveSequenceSource: sequenceRegistry.resolve
  });
  delete globalThis.__GBDRAW_PAIRWISE_TEST_METRICS__;

  assert.deepStrictEqual(traversalCounts, {
    rendered: 1,
    source: 1,
    orthogroup: 1,
    members: 1
  });
  assert.equal(metrics.contextCreations, 1);
  assert.equal(metrics.sequenceBundleCalls, 1);
  const memberRows = payload.group.memberRows;
  const resolvedMemberCount = memberRows.filter((row) => row.feature).length;
  assert.equal(resolvedMemberCount, 2);
  assert.equal(metrics.fastaCalls, resolvedMemberCount);

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
  'gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz',
  repoRoot
);
const exactSession = JSON.parse(gunzipSync(await readFile(exactFixturePath)).toString('utf8'));
const exactAdmission = admitFeatureCatalog(
  exactSession.editorState.featureCatalog,
  exactSession.results,
  { mode: exactSession.renderRequest.mode }
);
const exactState = exactAdmission.featureState;
const exactFeatureLookup = new Map(
  exactState.extractedFeatures.map((feature) => [feature.svg_id, feature])
);
const exactSequenceRegistry = createSequenceSourceRegistry(exactState.sequenceSources);
const exactOptions = {
  featureLookup: exactFeatureLookup,
  sourceFeatures: exactState.biologicalFeatures,
  // The app passes the catalog groups, collinear presentation groups included
  // (svg-actions.js groupsForMatch); adding collinearGroups would list them twice.
  orthogroups: exactState.orthogroups,
  resolveSequenceSource: exactSequenceRegistry.resolve
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
    digest: 'aacd129d83820a2532174227fd21d685a98f3c0c10ecf0a4a5af13df13d061b6',
    hoverRows: [
      { id: 'Kind', label: 'Kind', value: 'Collinearity block' },
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
    combinedSequenceLength: 33898
  },
  laterDistinctBlock: {
    target: {
      matchId: 'comparison1_match2',
      blockId: 'block_0155',
      orthogroupIds: 'og_274;og_275;og_276;og_277;og_278;og_279;og_280;og_281;og_282;og_283;og_284;og_285;og_286;og_287;og_288;og_289;og_290;og_291;og_292;og_293'
    },
    digest: '6cc891d0fc0d80ce1e305b0da66085989c6f646cce74dadbcbc3d0cdbc7f06c0',
    hoverRows: [
      { id: 'Kind', label: 'Kind', value: 'Collinearity block' },
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
    combinedSequenceLength: 19668
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
      combinedSequenceLength: payload.sequenceBundle.combinedFasta.length
    };
  }
  assert.deepStrictEqual(captured, EXACT_BASELINES);
});
