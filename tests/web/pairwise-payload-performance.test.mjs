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
  'gbdraw/web/js/app/feature-utils.js',
  'gbdraw/web/js/app/feature-sequence-fasta.js',
  'gbdraw/web/js/app/match-sequences.js',
  'gbdraw/web/js/app/conservation-series.js',
  'gbdraw/web/js/app/color-utils.js',
  'gbdraw/web/js/app/losat-normalization.js',
  'gbdraw/web/js/services/file-content-cache.js',
  'gbdraw/web/js/services/feature-catalog.js',
  'gbdraw/web/js/services/feature-identity.js',
  'gbdraw/web/js/services/orthogroup-feature-metadata.js',
  'gbdraw/web/js/services/runtime-test-hooks.js'
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

const popupPath = join(tempRoot, 'gbdraw/web/js/app/pairwise-match-popup.js');
const instrumentedPopupPath = join(
  tempRoot,
  'gbdraw/web/js/app/pairwise-match-popup-instrumented.js'
);
const popupSource = await readFile(popupPath, 'utf8');
const instrumentedPopupSource = popupSource
  .replace(
    "import { buildFeatureSequenceFastas } from './feature-sequence-fasta.js';",
    "import { buildFeatureSequenceFastas as buildFeatureSequenceFastasBase } from './feature-sequence-fasta.js';"
  )
  .replace(
    "import { buildMatchSequenceBundle } from './match-sequences.js';",
    "import { buildMatchSequenceBundle as buildMatchSequenceBundleBase } from './match-sequences.js';"
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
  'gbdraw/web/js/app/match-sequences.js'
)));

const elementFrom = (attributes) => ({
  style: { fill: attributes.fill || '' },
  getAttribute: (name) => attributes[name] || ''
});

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
  'data-gbdraw-pairwise-match-id': 'match-synthetic',
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
    const payload = popupModule.buildPairwiseMatchPayload(fixture.element, fixture.options);
    const featureTokens = new WeakMap(
      fixture.featureReferences.map((feature, index) => [feature, `${name}:${index}`])
    );
    const canonicalPayload = sortedCanonical(payload, featureTokens);
    const hoverRows = popupModule.buildPairwiseMatchHoverRows(payload);
    captured[name] = { payload: canonicalPayload, hoverRows };
    for (const reference of fixture.featureReferences) {
      const exposed = payload.sections
        .flatMap((entry) => entry.featureRows || [])
        .some((row) => row.feature === reference)
        || (payload.blockOrthogroups || [])
          .flatMap((group) => group.memberRows || [])
          .some((row) => row.feature === reference)
        || payload.sections
          .flatMap((entry) => entry.memberRows || [])
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
    const payload = popupModule.buildPairwiseMatchPayload(fixture.element, fixture.options);
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
  const orthogroupPayload = popupModule.buildPairwiseMatchPayload(
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
  const payload = instrumentedPopupModule.buildPairwiseMatchPayload(element, {
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
  const memberRows = payload.sections
    .find((entry) => entry.title === 'Similarity group')
    .memberRows;
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
  const payloadEnd = popupSource.indexOf('// Compatibility export retained', payloadStart);
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
  svg.matchAll(/<path\b(?=[^>]*data-gbdraw-pairwise-match-id=)[^>]*>/g),
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
  orthogroups: [...exactState.orthogroups, ...exactState.collinearGroups],
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
    digest: '7f55fa2faa8ee3a3fb798334b9b4ac74954e55ffafc5d15e3b7cf57c43f7093f',
    hoverRows: [
      { label: 'Kind', value: 'collinear' },
      { label: 'Identity', value: '43.9684414414' },
      { label: 'Query', value: '61879..78669' },
      { label: 'Subject', value: '19161..35565' },
      { label: 'Similarity groups', value: '15' },
      { label: 'Block', value: 'block_0024' }
    ],
    sectionTitles: ['Summary', 'Similarity groups covered', 'Collinearity', 'Query', 'Subject'],
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
    digest: '08e079f6f37ce48da62331aed504303e2e6239094c77e4f32f227ea55d02d8c6',
    hoverRows: [
      { label: 'Kind', value: 'collinear' },
      { label: 'Identity', value: '64.1326369725' },
      { label: 'Query', value: '484749..494272' },
      { label: 'Subject', value: '586163..595833' },
      { label: 'Similarity groups', value: '20' },
      { label: 'Block', value: 'block_0155' }
    ],
    sectionTitles: ['Summary', 'Similarity groups covered', 'Collinearity', 'Query', 'Subject'],
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
    const payload = popupModule.buildPairwiseMatchPayload(target.element, exactOptions);
    assert.deepStrictEqual(
      popupModule.buildPairwiseMatchHoverSummary(target.element, exactOptions),
      expectedHoverSummary(payload),
      `${name} exact-fixture hover summary changed`
    );
    const canonical = sortedCanonical(payload, exactFeatureTokens);
    const featureRows = payload.sections.flatMap((entry) => entry.featureRows || []);
    const memberRows = [
      ...payload.sections.flatMap((entry) => entry.memberRows || []),
      ...payload.blockOrthogroups.flatMap((group) => group.memberRows || [])
    ];
    const exposedFeatures = [...featureRows, ...memberRows]
      .map((row) => row.feature)
      .filter(Boolean);
    exposedFeatures.forEach((feature) => {
      assert.ok(exactFeatureTokens.has(feature), `${name} exposed a non-canonical Feature object`);
    });
    captured[name] = {
      target: {
        matchId: target.attributes['data-gbdraw-pairwise-match-id'],
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
