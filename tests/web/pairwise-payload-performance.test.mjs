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
  'gbdraw/web/js/services/record-source-coordinates.js',
  'gbdraw/web/js/services/genbank-header.js',
  'gbdraw/web/js/services/conservation-series.js',
  'gbdraw/web/js/utils/color-utils.js',
  'gbdraw/web/js/utils/named-colors.js',
  'gbdraw/web/js/services/losat-normalization.js',
  'gbdraw/web/js/contracts/request-shape.generated.js',
  'gbdraw/web/js/utils/contract-shape.js',
  'gbdraw/web/js/utils/field-reading.js',
  'gbdraw/web/js/services/file-content-cache.js',
  'gbdraw/web/js/services/feature-catalog.js',
  'gbdraw/web/js/contracts/schema.generated.js',
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

const SYNTHETIC_BASELINES_GZIP = 'H4sIAAAAAAACA+1c227bOBp+FYKLAVrUduxkZpG66AJOephsO+mMk14s1oFBS7TNGZnUiFQ63jbPtff7ZAseRJGSZSuJU7uN0ItGFPmT/P7zT8qfIUtCQlGyhP3PMEbLiKFQ/hmwKCIUo4SI5UnEgj/OQtiHsAWnJIpgH/6td3j0409/hy04S1gayyETFM4w7EM2ay8QobAFA5ZScYn/ErAPD8EMU8wBoeAQJDhgSchhC4aYBwmJBWFU0w8JjyO0PEcLSevPFCdLgOl/lgsMW5CEHv0FXkxwMkCnch7Y7+UtbxAXCPbhP377+Hr4r04PuIRG9Jd3I3f8GxJhqic0xMeq/1j3HyPUmSJkB5yyeGl2NVQbGYlTpnEUmIMnzw7aT0fi14QJTCg4ezUSQxZh2YlOSYhpgEdiwDmZ0QWmAiQYcUbVgDANBDgAlAk8omoJ7cTM0Ot2Or1j8OTZ05EwmxoJ/c/fG08nv+NA2IFH3U7n6Bg8kWu6+Hjyz9enl/lQ0zlHONuhAhQe2pZzUQT5XFiQ3ZX2e91277iA9+Dy7WAwuBwMHNjPxSbYqehMaQ77kH3isP/vzxDVZC9axVgzIucpspwYKkZoMQwQ/RBjCvsiSbGU5Ix15nXOb9iHOW+kgmAk0gRLjRiPzcMQT3FiRmcK1+/Cmxakd0aRijW7s9DFWqg+JOdMrFCoWAupUm4zFragXkPW2yxJtrPI7F/+ZSjCm5bLkircrehtFflcuG+DfM9HvgrMfMmVcJb0xwXUjnch9dWzGtSrbMwpu8YJUpZVqqPAwZySAEVWHSI0wdIkX5AFiZTBBsomg7NXsAWvUZS6Vllyy44IWIydPqcJ49ysDPACObUkQURUEqMb4wQ0uT6cRWyCosxcL5AI5m2+pGKOBQmkOsuWd4TKtzEiySfC5SpYIuZMETrzzbya7I3m58X1TL1NMA1xgsO2eglbkONAOhENSVLE5je1Xot5tmFPvj1kDF9LIwrsc8foOQgVOLlW28/GaPOwkv6K3lqkPbQv0sVCumhJobS1sxBTQcTSIXF83PnJm24QZc4mwnQm5k7f517H123d7iwet3tdr88JEYAHLHF7/XjoLdjOpwZmqpitvKTh1p3WcHm6S8SClGeeTwYWvhesiCNkx1tb6D/wsiRv7S5s5Xy3djNiATKRjOcS1GrPcoOqnnNrcl+zzNPJe7MWhz44AN6ueZXmiGWcm5yidA1rq8z7bPNFwVcguD3zyMjpm22z0FHhU5g9tzve7Cn3KbpYu11fackAytoXSSuwXFFWWn1LMa4ZgmXdjDBnvkSLc8mz+ALtdr616yuIdLYOV6hd3+WItedvc7H2duIJ9hYcpCPe3jzgABRQKIm4eX9rIV9j5VeIuQPKRjHPt7xe0Au4bRB1H/8awu4B57sajZhs4/jPVMrPSUrDCOuscDEhFIc2Zi04dx2/f9GwvXRNxRcVt/GX2iR84SJBNHz5zAtqS9TMKjN6Plcyihr91RRVvGiWnMd3xWkmTMxNfIepSAg2acY1IhGayJ1rDZ8+yK51Tl21Np0QuTTHKjUw62UJwVRk2vkM5kx7rx19/3kL8hhRyT1jQYbMMWxyy1JlescyMGMRDlL1lqZBhJkgIba6eVadEdggjKVJgN8p86LLB32Tz0jbwgVKhHIL0gil1AKch/8ql6iN+zbkYw32ht7Ypzs+6raPtoD/hbVNmgNHNTlQlUDktq6aCz2HC0druHB1k2/k0liGXyQ6OARZO9emWc5ZFZs7S0onmYXJ4/qs5VeTAwDFALmqucx5SgmOShlyG2ZTB9fe1YmGM8nfHJ+vCMtb1QnXqmxLIpknNk2RrSmyNUW2r1Jky1SuKbN9/TKbxb4ptG270FajpOa4my0W1a7u5u0dN2/m7he3VMvZO5u6rQuum/qvSLXeSm/qdDzUDt367YI/n0gfriI9LBApC47y8b60qCEOK93ZCxAEUcoFTrwuH5z404mMokIWeMoiloAFC7EHqR07JlnY5FUNtT4Asjqm6v5QrjHicEWFESDkFxAVDKtKiC45GsxZ4kZp4AD0lG5MUSAksPACLTBgHgTHz39w19tTARGQoSK8akH5IKNNzRYbVPnhkMzwIyMNvBgP2eqJaVBFE60u/f2vbFrBfWylTWu81nNp/wp3OccefeVOqr5S2Z8xCgmdabum0bHW3+q+yftakGpHr40Optc4YjHmYIHkWU0QpSHWpzAzTElgs1uAaGiSKTFHAqAEy4wBIGWWeEfigUiSx+vyiTtJi3viNBL6KWNH9hw72YttkfCPhH8AlT/nJLIWl4jTpskUAFKdim3SGdZKhLwwWdxOV+wwYzdHomhmFYCZuLh227Gv9rjMJu0t22Zwt+3aD5nW52ASry1kCVkBOcvJqwCwZLlBhouk97RgVpzFmWJBaXG2fcPiSppUvTy/K8i44CxQnR2aqXLeq8CjD6p3B645WEf8BSi63xfA9b4vgMvEm5YOdT7YOC5XncIL7ZwfotCx1rM/aKVDzWxSbr6/pQ5t3UZiTmbzOxc99Kbk/yFJF/tYAPG5sa8VEMmE+1VBCuHet1gEQVkasKNSiBbi+xVECjHcntVDtNT7COs28L//Ar1/4GByh4JJOQrZlKc7dtw5O82y5BUlkzXFFm31g0Ko6BZZfL8g36yuzBeT9apCfXYRCoW/owBTMZZRW80LUXltYcvFm+YwtznMbQ5zH+th7h1Kt2Vrl705dayjTiRq13Fz87arU9uCM3LLu71ymXJdqfZKbXrBIjZbfqSOpN/vkJemUXSPG7XZisr+49FepXV6DwmdAf3gHMewRYwSwhkFiSwseeSV/gFCQ/yXM8SXTRvtAU7C0slCRdcNEH4n939rhR7bcu1TFHHHx/geoZaF71ZZ+DygfzgvW8utXmQVSsKB875TdLcPiEWuLw/s8Gp5OEd7bfX2E+KqXMvTOI4IDsGUJUDMCQcn7wcXl0BP0nnYm05lO529+dnYaGVubnntydr3r+pAJU5oMSGzlKX8TZ7w7tDRVX860nwz8v1+M2LM2oaj1VWnqauOVmsep2pJ3Xhc6h+Blo5Max2XFko5fre6x6Q7/+Zj80X4qq8zmg851tYSizZ45UF+80FH80HH3T7o+AYC8X0sde1BTL6/FahqcHYTfH9Tnxvo64Z0GpFAEDp7TcOYESqa0LsJvZvQuwm9m9D7a17FLJnhJvpuou/mc+rmBL45gW9O4Hf0OfV3lOeEaRyRAAn8drvfShfuDm++J+xcDbZrKt8P7pbvB1deArZkxt5D9UXg0iXZ8h3ZbvmObOW114r5V199XXXB734fROYobueryALJLXwV6fJ5J59Gevt5uE8inWmUVZIsCN/u728TNFf2myv7j+nKvqOQzaX9XV7adxnxSK/t3/53DnZ7b/+7+HWE0p52//MIxSVt+H2Em/8Dwch8jVVcAAA=';
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
    digest: 'a1b15f0356b7c7aec35188f9487eb3a9cd33488c698714ae042b757d70090666',
    hoverRows: [
      { label: 'Kind', value: 'collinear' },
      { label: 'Identity', value: '43.9684414414' },
      { label: 'Query', value: '61879..78669' },
      { label: 'Subject', value: '19161..35565' },
      { label: 'Similarity groups', value: '15' },
      { label: 'Block', value: 'block_0024' }
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
    digest: 'a86114d1d7195b48815b02e86ffca8a08e21120c2d29abd2e01bfcb2a15e421a',
    hoverRows: [
      { label: 'Kind', value: 'collinear' },
      { label: 'Identity', value: '64.1326369725' },
      { label: 'Query', value: '484749..494272' },
      { label: 'Subject', value: '586163..595833' },
      { label: 'Similarity groups', value: '20' },
      { label: 'Block', value: 'block_0155' }
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
