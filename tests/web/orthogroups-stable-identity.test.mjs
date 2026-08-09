import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-orthogroup-stable-identity-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await mkdir(join(tempDir, 'app', 'feature-search'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await mkdir(join(tempDir, 'utils'), { recursive: true });

const copyModule = async (sourceRelative, targetRelative) => {
  await writeFile(
    join(tempDir, targetRelative),
    await readFile(join(repoRoot, sourceRelative), 'utf8'),
    'utf8'
  );
};

await copyModule('gbdraw/web/js/app/feature-utils.js', 'app/feature-utils.js');
await copyModule('gbdraw/web/js/app/feature-sequence-fasta.js', 'app/feature-sequence-fasta.js');
await copyModule('gbdraw/web/js/app/losat-normalization.js', 'app/losat-normalization.js');
await copyModule('gbdraw/web/js/app/feature-search/search-core.js', 'app/feature-search/search-core.js');
await copyModule(
  'gbdraw/web/js/services/standalone-interactivity-assets.js',
  'services/standalone-interactivity-assets.js'
);
await copyModule(
  'gbdraw/web/js/services/svg-serialization.js',
  'services/svg-serialization.js'
);
await copyModule(
  'gbdraw/web/js/services/orthogroup-feature-metadata.js',
  'services/orthogroup-feature-metadata.js'
);
await copyModule('gbdraw/web/js/services/feature-identity.js', 'services/feature-identity.js');
await copyModule('gbdraw/web/js/services/text-download.js', 'services/text-download.js');
await copyModule('gbdraw/web/js/utils/clipboard.js', 'utils/clipboard.js');
const standaloneSource = await readFile(
  join(repoRoot, 'gbdraw/web/js/services/standalone-interactivity.js'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'standalone-interactivity.js'),
  `${standaloneSource}\nexport { buildStandaloneBiologicalFeaturePayloads, buildStandaloneOrthogroupPayloads, createStandaloneBiologicalFeatureResolver, getStandaloneFeatureOrthogroupEntry, selectStandaloneCatalogItem, selectStandaloneSequenceSources };\n`,
  'utf8'
);
await writeFile(join(tempDir, 'app', 'feature-editor-svg-actions.js'), `
  export const FEATURE_SELECTOR = '[data-gbdraw-feature-id]';
  export const getFeatureIdentity = (element) => String(
    element?.getAttribute?.('data-gbdraw-feature-id') || ''
  ).replace(/__part\\d+$/, '');
  export const getFeatureElements = (svg, featureId) => (
    svg?.featureElements?.get?.(String(featureId || '')) || []
  );
`, 'utf8');

const orthogroupSource = (await readFile(
  join(repoRoot, 'gbdraw/web/js/app/orthogroups.js'),
  'utf8'
)).replace(
  "from './feature-editor/svg-actions.js';",
  "from './feature-editor-svg-actions.js';"
);
await writeFile(join(tempDir, 'app', 'orthogroups.js'), orthogroupSource, 'utf8');

globalThis.window = {
  Vue: {
    computed: (getter) => ({ get value() { return getter(); } })
  }
};
globalThis.isSecureContext = true;
let copiedText = '';
Object.defineProperty(globalThis, 'navigator', {
  configurable: true,
  value: {
    clipboard: {
      writeText: async (value) => { copiedText = String(value); }
    }
  }
});

const {
  createOrthogroupEditor,
  resolveUniqueOrthogroupMemberForFeature
} = await import(
  pathToFileURL(join(tempDir, 'app', 'orthogroups.js'))
);
const {
  isInternalProteinDisplayId,
  resolveDisplayProteinId
} = await import(pathToFileURL(join(tempDir, 'app', 'feature-utils.js')));
const { runFeatureSearch } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-search', 'search-core.js'))
);
const {
  buildOrthogroupFeatureIndex,
  enrichFeatureWithOrthogroup,
  enrichFeaturesWithOrthogroups
} = await import(
  pathToFileURL(join(tempDir, 'services', 'orthogroup-feature-metadata.js'))
);
const {
  buildStandaloneBiologicalFeaturePayloads,
  buildStandaloneOrthogroupPayloads,
  createStandaloneBiologicalFeatureResolver,
  getStandaloneFeatureOrthogroupEntry,
  selectStandaloneCatalogItem,
  selectStandaloneSequenceSources
} = await import(
  pathToFileURL(join(tempDir, 'services', 'standalone-interactivity.js'))
);
const { STANDALONE_INTERACTIVE_SCRIPT } = await import(
  pathToFileURL(join(tempDir, 'services', 'standalone-interactivity-assets.js'))
);

const embeddedFunctionSource = (name) => {
  const marker = `function ${name}(`;
  const start = STANDALONE_INTERACTIVE_SCRIPT.indexOf(marker);
  assert.notEqual(start, -1, `missing embedded function ${name}`);
  const bodyStart = STANDALONE_INTERACTIVE_SCRIPT.indexOf('{', start);
  let depth = 0;
  for (let index = bodyStart; index < STANDALONE_INTERACTIVE_SCRIPT.length; index += 1) {
    const character = STANDALONE_INTERACTIVE_SCRIPT[index];
    if (character === '{') depth += 1;
    if (character === '}') depth -= 1;
    if (depth === 0) return STANDALONE_INTERACTIVE_SCRIPT.slice(start, index + 1);
  }
  assert.fail(`unterminated embedded function ${name}`);
};

const catalogRuntime = new Function(`
  ${[
    'catalogFeatureKey',
    'catalogCanonicalFeatureIdentity',
    'catalogRenderedFeatureIdentity',
    'catalogStableFeatureIdentity',
    'catalogRecordIndex',
    'catalogStableFeatureId',
    'catalogSequenceSourceStrings',
    'catalogNucleotideSequence',
    'expandCatalogBiologicalFeature',
    'firstCatalogText',
    'firstCatalogValue',
    'catalogQualifier',
    'catalogSourceReferenceAgrees',
    'catalogMatchRoleTextIdentity',
    'catalogMatchEndpointReference',
    'catalogMatchEndpoints',
    'clearCatalogMatchEndpoint',
    'expandCatalogItem',
    'consistentTextIdentity',
    'nonnegativeIntegerIdentity',
    'reverseComplementMatchSequence'
  ].map(embeddedFunctionSource).join('\n')}
  return { expandCatalogItem };
`)();
const embeddedCatalogItem = {
  recordKeys: ['record-a'],
  biologicalFeatures: [{
    recordKey: 'record-a',
    biologicalFeatureId: 'biological-a',
    stableFeatureId: 'stable-a',
    sourceFeatureIndex: 4,
    record_id: 'record-a',
    nucleotide_sequence: 'ACGT'
  }],
  features: [{
    recordKey: 'record-a',
    biologicalFeatureId: 'biological-a',
    svgId: 'rendered-a'
  }],
  orthogroups: [],
  comparisonMatches: [{
    id: 'match-a',
    match_kind: 'pairwise',
    queryRecordKey: 'record-a',
    queryBiologicalFeatureId: 'biological-a',
    query_record_id: 'record-a',
    query_record_index: 0,
    query_feature_index: 4
  }]
};
const expandedCatalogItem = catalogRuntime.expandCatalogItem(embeddedCatalogItem, 'rich');
assert.equal(expandedCatalogItem.features[0].svg_id, 'rendered-a');
assert.equal(expandedCatalogItem.features[0].sourceFeatureIndex, 4);
assert.equal(expandedCatalogItem.features[0].feature_index, 4);
assert.equal(expandedCatalogItem.matches[0]._gbdraw_query_endpoint_resolved, true);
assert.equal(expandedCatalogItem.matches[0].query_feature_svg_id, 'rendered-a');
assert.equal(expandedCatalogItem.matches[0].query_stable_feature_svg_id, 'stable-a');
const duplicateLocationCatalog = {
  ...embeddedCatalogItem,
  biologicalFeatures: [
    embeddedCatalogItem.biologicalFeatures[0],
    {
      ...embeddedCatalogItem.biologicalFeatures[0],
      biologicalFeatureId: 'biological-b',
      sourceFeatureIndex: 7
    }
  ],
  features: [
    embeddedCatalogItem.features[0],
    {
      ...embeddedCatalogItem.features[0],
      biologicalFeatureId: 'biological-b',
      svgId: 'rendered-b'
    }
  ],
  comparisonMatches: [{
    id: 'match-duplicates',
    match_kind: 'pairwise',
    queryFeatureReferences: [{
      recordKey: 'record-a',
      biologicalFeatureId: 'biological-a'
    }, {
      recordKey: 'record-a',
      biologicalFeatureId: 'biological-b'
    }]
  }]
};
const expandedDuplicateLocationCatalog = catalogRuntime.expandCatalogItem(
  duplicateLocationCatalog,
  'rich'
);
assert.equal(
  expandedDuplicateLocationCatalog.matches[0].query_feature_svg_id,
  'rendered-a;rendered-b'
);
assert.equal(
  expandedDuplicateLocationCatalog.matches[0].query_stable_feature_svg_id,
  'stable-a;stable-a'
);
assert.equal(
  expandedDuplicateLocationCatalog.matches[0].query_feature_index,
  '4;7'
);
const missingDuplicateSourceIndex = structuredClone(duplicateLocationCatalog);
delete missingDuplicateSourceIndex.biologicalFeatures[1].sourceFeatureIndex;
assert.throws(
  () => catalogRuntime.expandCatalogItem(missingDuplicateSourceIndex, 'rich'),
  /unique source feature indexes/
);
assert.throws(
  () => catalogRuntime.expandCatalogItem({
    ...embeddedCatalogItem,
    biologicalFeatures: [
      ...embeddedCatalogItem.biologicalFeatures,
      { ...embeddedCatalogItem.biologicalFeatures[0], product: 'wrong duplicate' }
    ]
  }, 'rich'),
  /Duplicate or invalid catalog biological feature identity/
);
assert.throws(
  () => catalogRuntime.expandCatalogItem({
    ...embeddedCatalogItem,
    features: [
      ...embeddedCatalogItem.features,
      { ...embeddedCatalogItem.features[0], svgId: 'rendered-copy' }
    ]
  }, 'rich'),
  /Duplicate or unresolved catalog rendered feature reference/
);
assert.throws(
  () => catalogRuntime.expandCatalogItem({
    ...embeddedCatalogItem,
    orthogroups: [{
      id: 'og-a',
      orthogroupId: 'og-b',
      members: []
    }]
  }, 'rich'),
  /Invalid or duplicate catalog orthogroup identity/
);
assert.throws(
  () => catalogRuntime.expandCatalogItem({
    ...embeddedCatalogItem,
    orthogroups: [
      { id: 'og-duplicate', members: [] },
      { orthogroupId: 'og-duplicate', members: [] }
    ]
  }, 'rich'),
  /Invalid or duplicate catalog orthogroup identity/
);
const expandedCatalogOrthogroup = catalogRuntime.expandCatalogItem({
  ...embeddedCatalogItem,
  orthogroups: [{
    id: 'og-consistent',
    orthogroupId: 'og-consistent',
    orthogroup_id: 'og-consistent',
    members: [{
      recordKey: 'record-a',
      biologicalFeatureId: 'biological-a'
    }]
  }]
}, 'rich');
assert.equal(expandedCatalogOrthogroup.orthogroups[0].id, 'og-consistent');
assert.equal(expandedCatalogOrthogroup.features[0].orthogroup_id, 'og-consistent');
assert.throws(
  () => catalogRuntime.expandCatalogItem({
    ...embeddedCatalogItem,
    biologicalFeatures: [
      ...embeddedCatalogItem.biologicalFeatures,
      {
        recordKey: 'record-a',
        biologicalFeatureId: 'biological-b',
        stableFeatureId: 'stable-b',
        record_id: 'record-a',
        nucleotide_sequence: 'TGCA'
      }
    ],
    features: [
      ...embeddedCatalogItem.features,
      {
        recordKey: 'record-a',
        biologicalFeatureId: 'biological-b',
        svgId: 'rendered-a'
      }
    ]
  }, 'rich'),
  /Duplicate or unresolved catalog rendered feature reference/
);
assert.throws(
  () => catalogRuntime.expandCatalogItem({
    ...embeddedCatalogItem,
    biologicalFeatures: [{
      ...embeddedCatalogItem.biologicalFeatures[0],
      record_key: 'different-record'
    }]
  }, 'rich'),
  /Invalid catalog biological feature identity/
);
const conflictingCatalogEndpoint = catalogRuntime.expandCatalogItem({
  ...embeddedCatalogItem,
  comparisonMatches: [{
    ...embeddedCatalogItem.comparisonMatches[0],
    query_feature_svg_id: 'tampered-rendered-id'
  }]
}, 'rich').matches[0];
assert.equal(conflictingCatalogEndpoint._gbdraw_query_endpoint_resolved, false);
assert.equal(conflictingCatalogEndpoint.query_feature_svg_id, undefined);
assert.equal(conflictingCatalogEndpoint.query_record_index, undefined);
assert.equal(conflictingCatalogEndpoint.query_record_id, undefined);
for (const conflictingEndpointIdentity of [
  { query_record_index: 1 },
  { queryBiologicalFeatureId: 'different-biological-feature' },
  { queryStableFeatureId: 'different-stable-feature' }
]) {
  const unresolvedEndpoint = catalogRuntime.expandCatalogItem({
    ...embeddedCatalogItem,
    comparisonMatches: [{
      ...embeddedCatalogItem.comparisonMatches[0],
      ...conflictingEndpointIdentity
    }]
  }, 'rich').matches[0];
  assert.equal(unresolvedEndpoint._gbdraw_query_endpoint_resolved, false);
  assert.equal(unresolvedEndpoint.query_feature_svg_id, undefined);
  assert.equal(unresolvedEndpoint.query_record_index, undefined);
}
const catalogEndpointRuntime = new Function('featuresById', 'sequenceSources', `
  ${[
    'catalogFeatureKey',
    'catalogCanonicalFeatureIdentity',
    'catalogMatchRoleTextIdentity',
    'consistentTextIdentity',
    'nonnegativeIntegerIdentity',
    'featureRecordIdentity',
    'biologicalFeatureStableSvgId',
    'resolvedCatalogMatchFeature',
    'matchSourceAliases',
    'resolveEmbeddedMatchSource'
  ].map(embeddedFunctionSource).join('\n')}
  return { resolvedCatalogMatchFeature, resolveEmbeddedMatchSource };
`)(
  new Map([['rendered-a', expandedCatalogItem.features[0]]]),
  [{ origin: 'linear-record', recordIndex: 0, recordId: 'record-a', sequence: 'ACGT' }]
);
assert.equal(
  catalogEndpointRuntime.resolvedCatalogMatchFeature(
    expandedCatalogItem.matches[0],
    'query'
  ),
  expandedCatalogItem.features[0]
);
assert.equal(
  catalogEndpointRuntime.resolveEmbeddedMatchSource(
    expandedCatalogItem.matches[0],
    'query'
  ).source.sequence,
  'ACGT'
);
assert.equal(
  catalogEndpointRuntime.resolvedCatalogMatchFeature(
    conflictingCatalogEndpoint,
    'query'
  ),
  null
);
assert.equal(
  catalogEndpointRuntime.resolveEmbeddedMatchSource(
    conflictingCatalogEndpoint,
    'query'
  ).source,
  null
);

assert.equal(isInternalProteinDisplayId(`h_${'c'.repeat(26)}`), true);
assert.equal(isInternalProteinDisplayId(`f_${'d'.repeat(64)}`), true);
assert.equal(
  resolveDisplayProteinId(null, {
    proteinId: `h_${'c'.repeat(26)}`,
    label: 'display fallback'
  }),
  'display fallback'
);
assert.equal(
  resolveDisplayProteinId(null, { proteinId: `h_${'c'.repeat(26)}` }),
  ''
);
const runtimeHandle = `h_${'e'.repeat(26)}`;
const featureAnalysisId = `f_${'f'.repeat(64)}`;
const unsupportedHistoricalTransportId = `record@instance|alias~f_${'a'.repeat(64)}`;
assert.equal(
  resolveDisplayProteinId(
    {
      displayProteinId: runtimeHandle,
      sourceProteinId: featureAnalysisId,
      qualifiers: { protein_id: [unsupportedHistoricalTransportId] },
      locusTag: 'WP_SAFE.1'
    },
    { label: runtimeHandle }
  ),
  'WP_SAFE.1'
);
assert.equal(
  resolveDisplayProteinId(
    {
      displayProteinId: runtimeHandle,
      sourceProteinId: featureAnalysisId,
      qualifiers: { protein_id: [unsupportedHistoricalTransportId] }
    },
    { label: runtimeHandle }
  ),
  ''
);
assert.equal(
  resolveDisplayProteinId({
    qualifiers: { protein_id: [runtimeHandle, 'WP_ARRAY_FALLBACK.1'] }
  }),
  'WP_ARRAY_FALLBACK.1'
);

const ref = (value) => ({ value });
const visibleElement = {
  attrs: new Map([
    ['data-gbdraw-feature-id', 'stable-x_record_1'],
    ['data-gbdraw-stable-feature-id', 'stable-x'],
    ['data-gbdraw-record-index', '0'],
    ['stroke', '#111111'],
    ['stroke-width', '0.5']
  ]),
  getAttribute(name) { return this.attrs.get(name) ?? null; },
  hasAttribute(name) { return this.attrs.has(name); },
  setAttribute(name, value) { this.attrs.set(name, String(value)); },
  removeAttribute(name) { this.attrs.delete(name); }
};
const svg = {
  featureElements: new Map([['stable-x_record_1', [visibleElement]]]),
  querySelectorAll(selector) {
    if (selector === '[data-gbdraw-feature-id]') return [visibleElement];
    if (selector.includes('[data-og-original-stroke]')) {
      return visibleElement.hasAttribute('data-og-original-stroke') ? [visibleElement] : [];
    }
    return [];
  }
};

const group = {
  id: 'og_1',
  members: [
    {
      recordIndex: 0,
      recordId: 'record-a',
      featureSvgId: 'stable-x',
      renderedFeatureSvgId: 'stable-x_record_1',
      proteinId: `h_${'a'.repeat(26)}`
    },
    {
      recordIndex: 1,
      recordId: 'record-b',
      featureSvgId: 'stable-x',
      renderedFeatureSvgId: 'stable-x_record_2',
      proteinId: `h_${'b'.repeat(26)}`
    },
    {
      recordIndex: 2,
      recordId: 'record-c',
      recordKey: 'record-key-c',
      biologicalFeatureId: 'biological-c',
      featureSvgId: 'canonical-stable-id',
      stableFeatureSvgId: 'canonical-stable-id',
      renderedFeatureSvgId: 'reversed-display-stable-id_record_3',
      proteinId: `h_${'c'.repeat(26)}`,
      sourceProteinId: 'CAG34720.1'
    }
  ]
};
const orthogroupIndex = buildOrthogroupFeatureIndex([group]);
const generatedFeatures = enrichFeaturesWithOrthogroups([
  {
    fileIdx: 0,
    svg_id: 'stable-x_record_1',
    stable_svg_id: 'stable-x',
    rendered_svg_id: 'stable-x_record_1'
  },
  {
    fileIdx: 2,
    recordKey: 'record-key-c',
    biologicalFeatureId: 'biological-c',
    svg_id: 'reversed-display-stable-id_record_3',
    stable_svg_id: 'canonical-stable-id',
    rendered_svg_id: 'reversed-display-stable-id_record_3'
  },
  {
    fileIdx: 1,
    svg_id: 'stable-x_record_2',
    stable_svg_id: 'stable-x',
    rendered_svg_id: 'stable-x_record_2'
  }
], orthogroupIndex);
assert.deepEqual(
  generatedFeatures.map((feature) => feature.orthogroupId),
  ['og_1', 'og_1', 'og_1']
);
assert.equal(generatedFeatures[0].proteinId, `h_${'a'.repeat(26)}`);
assert.equal(
  enrichFeatureWithOrthogroup(
    orthogroupIndex,
    { fileIdx: 1, svg_id: 'stable-x_record_2', stable_svg_id: 'stable-x' }
  ).proteinId,
  `h_${'b'.repeat(26)}`
);
const ambiguousRecordlessFeature = { svg_id: 'stable-x', stable_svg_id: 'stable-x' };
assert.equal(
  enrichFeatureWithOrthogroup(orthogroupIndex, ambiguousRecordlessFeature),
  ambiguousRecordlessFeature
);
assert.equal(generatedFeatures[1].orthogroupMemberCount, 3);
assert.equal(generatedFeatures[1].orthogroupRecordCoverage, 3);
assert.equal(
  generatedFeatures[1].orthogroupMember.renderedFeatureSvgId,
  'reversed-display-stable-id_record_3'
);
assert.equal(generatedFeatures[1].sourceProteinId, 'CAG34720.1');
const conflictingIdentityGroup = {
  id: 'og_conflict',
  members: [
    {
      recordIndex: 0,
      stableFeatureSvgId: 'stable-a',
      renderedFeatureSvgId: 'rendered-a'
    },
    {
      recordIndex: 0,
      stableFeatureSvgId: 'stable-b',
      renderedFeatureSvgId: 'rendered-b'
    }
  ]
};
const conflictingIdentityFeature = {
  fileIdx: 0,
  svg_id: 'rendered-a',
  stable_svg_id: 'stable-b'
};
assert.equal(
  enrichFeatureWithOrthogroup(
    buildOrthogroupFeatureIndex([conflictingIdentityGroup]),
    conflictingIdentityFeature
  ),
  conflictingIdentityFeature
);
const roleSeparatedIndex = buildOrthogroupFeatureIndex([
  {
    id: 'og_role_a',
    members: [{
      recordIndex: 0,
      stableFeatureSvgId: 'stable-a',
      featureSvgId: 'stable-a',
      renderedFeatureSvgId: 'collision'
    }]
  },
  {
    id: 'og_role_b',
    members: [{
      recordIndex: 0,
      stableFeatureSvgId: 'collision',
      featureSvgId: 'collision',
      renderedFeatureSvgId: 'rendered-b'
    }]
  }
]);
assert.equal(
  enrichFeatureWithOrthogroup(
    roleSeparatedIndex,
    { fileIdx: 0, stable_svg_id: 'collision' }
  ).orthogroupId,
  'og_role_b'
);
assert.equal(
  enrichFeatureWithOrthogroup(
    roleSeparatedIndex,
    { fileIdx: 0, stable_svg_id: 'stable-a', svg_id: 'collision' }
  ).orthogroupId,
  'og_role_a'
);
const clickedOrthogroupMember = {
  recordIndex: 0,
  featureIndex: 7,
  recordKey: 'record-key-a',
  biologicalFeatureId: 'biological-a',
  stableFeatureSvgId: 'stable-a',
  featureSvgId: 'stable-a',
  renderedFeatureSvgId: 'collision'
};
const wrongFirstOrthogroupMember = {
  ...clickedOrthogroupMember,
  stableFeatureSvgId: 'collision',
  featureSvgId: 'collision',
  renderedFeatureSvgId: 'other'
};
const clickedOrthogroupFeature = {
  fileIdx: 0,
  feature_index: 7,
  record_key: 'record-key-a',
  biological_feature_id: 'biological-a',
  stable_svg_id: 'stable-a',
  rendered_svg_id: 'collision',
  svg_id: 'collision'
};
assert.equal(
  resolveUniqueOrthogroupMemberForFeature(
    clickedOrthogroupFeature,
    [wrongFirstOrthogroupMember, clickedOrthogroupMember]
  ),
  clickedOrthogroupMember
);
assert.equal(
  resolveUniqueOrthogroupMemberForFeature(
    { ...clickedOrthogroupFeature, stableFeatureId: 'collision' },
    [clickedOrthogroupMember]
  ),
  null
);
assert.equal(
  resolveUniqueOrthogroupMemberForFeature(
    { stable_svg_id: 'stable-a', rendered_svg_id: 'collision', svg_id: 'collision' },
    [clickedOrthogroupMember]
  ),
  null
);
assert.equal(
  resolveUniqueOrthogroupMemberForFeature(
    clickedOrthogroupFeature,
    [clickedOrthogroupMember, { ...clickedOrthogroupMember }]
  ),
  null
);
assert.deepEqual(
  enrichFeatureWithOrthogroup(
    buildOrthogroupFeatureIndex([{
      id: 'og_conflicting_member_aliases',
      members: [{
        recordIndex: 0,
        stableFeatureSvgId: 'stable-a',
        featureSvgId: 'collision'
      }]
    }]),
    { fileIdx: 0, stable_svg_id: 'stable-a' }
  ),
  { fileIdx: 0, stable_svg_id: 'stable-a' }
);
const wrongRecordFeature = { fileIdx: 9, stable_svg_id: 'canonical-stable-id' };
assert.equal(
  enrichFeatureWithOrthogroup(orthogroupIndex, wrongRecordFeature),
  wrongRecordFeature
);
for (const invalidIdentity of [
  { fileIdx: 0, stable_svg_id: 'stable-x', rendered_svg_id: 'unknown-rendered' },
  { fileIdx: 0, recordIndex: -1, stable_svg_id: 'stable-x' },
  { fileIdx: '0.9', stable_svg_id: 'stable-x' },
  { fileIdx: 0, stable_svg_id: 'stable-x', stableFeatureId: 'wrong-stable' },
  { stable_svg_id: 'canonical-stable-id' }
]) {
  assert.deepEqual(
    enrichFeatureWithOrthogroup(orthogroupIndex, invalidIdentity),
    invalidIdentity
  );
}
const duplicateIdentityIndex = buildOrthogroupFeatureIndex([{
  id: 'og_duplicate',
  members: [
    { recordIndex: 0, stableFeatureSvgId: 'duplicate-stable' },
    { recordIndex: 0, stableFeatureSvgId: 'duplicate-stable' }
  ]
}]);
const duplicateIdentityFeature = { fileIdx: 0, stable_svg_id: 'duplicate-stable' };
assert.deepEqual(
  enrichFeatureWithOrthogroup(duplicateIdentityIndex, duplicateIdentityFeature),
  duplicateIdentityFeature
);
const canonicalMember = {
  recordIndex: 0,
  stableFeatureSvgId: 'canonical-source',
  recordKey: 'record-key-a',
  biologicalFeatureId: 'biological-a'
};
const canonicalIndex = buildOrthogroupFeatureIndex([{
  id: 'og_canonical',
  members: [canonicalMember]
}]);
assert.equal(
  enrichFeatureWithOrthogroup(canonicalIndex, {
    fileIdx: 0,
    stable_svg_id: 'canonical-source',
    recordKey: 'record-key-a',
    biologicalFeatureId: 'biological-a'
  }).orthogroupId,
  'og_canonical'
);
const conflictingCanonicalFeature = {
  fileIdx: 0,
  stable_svg_id: 'canonical-source',
  recordKey: 'record-key-a',
  biologicalFeatureId: 'missing-biological'
};
assert.deepEqual(
  enrichFeatureWithOrthogroup(canonicalIndex, conflictingCanonicalFeature),
  conflictingCanonicalFeature
);

const aliasedGroupFeature = { fileIdx: 4, stable_svg_id: 'aliased-group-feature' };
const aliasedGroupIndex = buildOrthogroupFeatureIndex([{
  orthogroupId: 'og_aliased',
  orthogroup_id: 'og_aliased',
  members: [{ recordIndex: 4, stableFeatureSvgId: 'aliased-group-feature' }]
}]);
assert.equal(
  enrichFeatureWithOrthogroup(aliasedGroupIndex, aliasedGroupFeature).orthogroupId,
  'og_aliased'
);
for (const invalidGroups of [
  [{
    id: 'og_conflicting_group_id',
    orthogroupId: 'og_other_group_id',
    members: [{ recordIndex: 4, stableFeatureSvgId: 'aliased-group-feature' }]
  }],
  [{
    members: [{ recordIndex: 4, stableFeatureSvgId: 'aliased-group-feature' }]
  }],
  [
    {
      id: 'og_duplicate_group_id',
      members: [{ recordIndex: 4, stableFeatureSvgId: 'duplicate-group-feature-a' }]
    },
    {
      orthogroupId: 'og_duplicate_group_id',
      members: [{ recordIndex: 5, stableFeatureSvgId: 'duplicate-group-feature-b' }]
    }
  ]
]) {
  const invalidGroupIndex = buildOrthogroupFeatureIndex(invalidGroups);
  for (const feature of [
    aliasedGroupFeature,
    { fileIdx: 4, stable_svg_id: 'duplicate-group-feature-a' },
    { fileIdx: 5, stable_svg_id: 'duplicate-group-feature-b' }
  ]) {
    assert.deepEqual(
      enrichFeatureWithOrthogroup(invalidGroupIndex, feature),
      feature
    );
  }
}

const standaloneRecordZero = {
  fileIdx: 0,
  svg_id: 'rendered-zero_record_1',
  stable_svg_id: 'shared-stable'
};
const standaloneRecordOne = {
  fileIdx: 1,
  svg_id: 'rendered-one_record_2',
  stable_svg_id: 'shared-stable'
};
const standaloneCanonicalFeature = {
  recordKey: 'standalone-record-key',
  biologicalFeatureId: 'standalone-biological-feature'
};
const standaloneResolver = createStandaloneBiologicalFeatureResolver([
  standaloneRecordZero,
  standaloneRecordOne,
  standaloneCanonicalFeature,
  { fileIdx: 0, svg_id: 'rendered-conflict_record_1', stable_svg_id: 'conflict-stable' }
]);
assert.equal(
  standaloneResolver({ stableSvgId: 'shared-stable', recordIndex: 0 }),
  standaloneRecordZero
);
assert.equal(
  standaloneResolver({ renderedSvgId: 'rendered-one_record_2', recordIndex: 1 }),
  standaloneRecordOne
);
assert.equal(
  standaloneResolver({
    recordKey: 'standalone-record-key',
    biologicalFeatureId: 'standalone-biological-feature'
  }),
  standaloneCanonicalFeature
);
assert.equal(standaloneResolver({ stableSvgId: 'shared-stable' }), null);
assert.equal(
  standaloneResolver({ stableSvgId: 'shared-stable', recordIndex: 9 }),
  null
);
assert.equal(
  standaloneResolver({
    stableSvgId: 'shared-stable',
    renderedSvgId: 'rendered-conflict_record_1',
    recordIndex: 0
  }),
  null
);
assert.equal(
  standaloneResolver({
    stableSvgId: 'shared-stable',
    renderedSvgId: 'unknown-rendered',
    recordIndex: 0
  }),
  null
);
assert.equal(
  standaloneResolver({ stableSvgId: 'shared-stable', recordIndex: '0.9' }),
  null
);
assert.equal(
  standaloneResolver({ stableSvgId: 'shared-stable', recordIndex: -1 }),
  null
);
assert.equal(
  standaloneResolver({
    stableSvgId: 'shared-stable',
    recordIndex: '9007199254740993'
  }),
  null
);
const duplicateStandaloneResolver = createStandaloneBiologicalFeatureResolver([
  { fileIdx: 0, stable_svg_id: 'duplicate-stable' },
  { fileIdx: 0, stable_svg_id: 'duplicate-stable' }
]);
assert.equal(
  duplicateStandaloneResolver({ stableSvgId: 'duplicate-stable', recordIndex: 0 }),
  null
);
const standaloneRoleA = {
  fileIdx: 0,
  stable_svg_id: 'standalone-stable-a',
  svg_id: 'standalone-collision'
};
const standaloneRoleB = {
  fileIdx: 0,
  stable_svg_id: 'standalone-collision',
  svg_id: 'standalone-rendered-b'
};
const roleSeparatedStandaloneResolver = createStandaloneBiologicalFeatureResolver([
  standaloneRoleA,
  standaloneRoleB
]);
assert.equal(
  roleSeparatedStandaloneResolver({
    stableSvgId: 'standalone-collision',
    recordIndex: 0
  }),
  standaloneRoleB
);
assert.equal(
  roleSeparatedStandaloneResolver({
    renderedSvgId: 'standalone-collision',
    recordIndex: 0
  }),
  standaloneRoleA
);
assert.equal(
  roleSeparatedStandaloneResolver({
    stableSvgId: 'standalone-collision',
    renderedSvgId: 'standalone-collision',
    recordIndex: 0
  }),
  null
);

const standaloneEntry = { orthogroupId: 'og_record_one' };
const standaloneCanonicalEntry = { orthogroupId: 'og_canonical' };
const standaloneIndex = new Map([
  ['1:shared-stable', standaloneEntry],
  ['shared-stable', standaloneEntry],
  ['canonical:standalone-record-key\0standalone-biological-feature', standaloneCanonicalEntry]
]);
const standaloneContext = { featureOrthogroupIndex: standaloneIndex };
assert.equal(
  getStandaloneFeatureOrthogroupEntry(
    { fileIdx: 1, stable_svg_id: 'shared-stable' },
    standaloneContext
  ),
  standaloneEntry
);
assert.equal(
  getStandaloneFeatureOrthogroupEntry(
    {
      recordKey: 'standalone-record-key',
      biologicalFeatureId: 'standalone-biological-feature'
    },
    standaloneContext
  ),
  standaloneCanonicalEntry
);
const catalogSelectionItem = { resultIndex: 0, resultName: 'result-zero' };
assert.equal(
  selectStandaloneCatalogItem({
    featureCatalog: { schema: 3, items: [catalogSelectionItem] },
    catalogResultIndex: '0',
    catalogResultName: 'result-zero'
  }),
  catalogSelectionItem
);
for (const invalidCatalogIndex of [
  '-1',
  '0.5',
  'bogus',
  '9007199254740992'
]) {
  assert.equal(
    selectStandaloneCatalogItem({
      featureCatalog: { schema: 3, items: [catalogSelectionItem] },
      catalogResultIndex: invalidCatalogIndex,
      catalogResultName: 'result-zero'
    }),
    null
  );
}
const linearSequenceSource = {
  origin: 'linear-record',
  recordIndex: 1,
  recordId: 'record-one',
  sequence: 'ACGT'
};
assert.deepEqual(
  selectStandaloneSequenceSources(
    [{ match_kind: 'pairwise', query_record_index: '1', query_record_id: 'record-one' }],
    [linearSequenceSource]
  ),
  [linearSequenceSource]
);
for (const invalidRecordIndex of [
  '-1',
  '1.5',
  'bogus',
  '9007199254740992'
]) {
  assert.deepEqual(
    selectStandaloneSequenceSources(
      [{ match_kind: 'pairwise', query_record_index: invalidRecordIndex }],
      [linearSequenceSource]
    ),
    []
  );
}
assert.deepEqual(
  selectStandaloneSequenceSources(
    [{ match_kind: 'pairwise', query_record_index: '1' }],
    [{ ...linearSequenceSource, recordIndex: '9007199254740992' }]
  ),
  []
);
assert.deepEqual(
  selectStandaloneSequenceSources(
    [{ match_kind: 'pairwise', query_record_index: '1', query_record_id: 'wrong-record' }],
    [linearSequenceSource]
  ),
  []
);
for (const unresolved of [
  { fileIdx: 0, stable_svg_id: 'shared-stable' },
  { stable_svg_id: 'shared-stable' },
  { fileIdx: 0, recordIndex: 1, stable_svg_id: 'shared-stable' },
  { fileIdx: 1, recordIndex: -1, stable_svg_id: 'shared-stable' },
  { fileIdx: 1, stable_svg_id: 'shared-stable', rendered_svg_id: 'unknown' },
  { fileIdx: 1, protein_id: 'CAG34720.1' }
]) {
  assert.equal(
    getStandaloneFeatureOrthogroupEntry(unresolved, standaloneContext),
    null
  );
}
const directGroupMember = {
  recordIndex: 1,
  stableFeatureSvgId: 'shared-stable'
};
const directContext = {
  ...standaloneContext,
  orthogroups: [{ id: 'og_direct', members: [directGroupMember] }]
};
assert.equal(
  getStandaloneFeatureOrthogroupEntry({
    fileIdx: 1,
    stable_svg_id: 'shared-stable',
    orthogroupId: 'og_direct',
    orthogroupMember: { ...directGroupMember }
  }, directContext)?.orthogroupId,
  'og_direct'
);
assert.equal(
  getStandaloneFeatureOrthogroupEntry({
    fileIdx: 1,
    stable_svg_id: 'shared-stable',
    orthogroupId: 'og_stale_direct',
    orthogroupMember: {
      recordIndex: 0,
      stableFeatureSvgId: 'shared-stable'
    }
  }, {
    ...standaloneContext,
    orthogroups: [{
      id: 'og_stale_direct',
      members: [{ recordIndex: 0, stableFeatureSvgId: 'shared-stable' }]
    }]
  }),
  null
);
assert.equal(
  getStandaloneFeatureOrthogroupEntry({
    fileIdx: 1,
    featureIndex: 7,
    stable_svg_id: 'shared-stable',
    orthogroupId: 'og_direct'
  }, directContext),
  null
);
assert.equal(
  getStandaloneFeatureOrthogroupEntry({
    fileIdx: 1,
    stable_svg_id: 'shared-stable',
    orthogroupId: 'og_direct',
    orthogroup_id: 'og_conflict'
  }, directContext),
  null
);

assert.deepEqual(
  buildStandaloneOrthogroupPayloads([], {
    orthogroups: [{ id: 'og_a', orthogroupId: 'og_b', members: [] }]
  }),
  []
);
assert.deepEqual(
  buildStandaloneOrthogroupPayloads([], {
    orthogroups: [
      { id: 'og_duplicate', members: [] },
      { id: 'og_duplicate', members: [] }
    ]
  }),
  []
);
const conflictingStandaloneMemberGroups = buildStandaloneOrthogroupPayloads([], {
  orthogroups: [{
    id: 'og_conflicting_member',
    members: [{
      recordKey: 'record-key',
      biologicalFeatureId: 'bio-id',
      stableFeatureSvgId: 'stable-a',
      featureSvgId: 'stable-b'
    }]
  }]
});
assert.equal(conflictingStandaloneMemberGroups.length, 1);
assert.deepEqual(conflictingStandaloneMemberGroups[0].members, []);

assert.deepEqual(
  buildStandaloneBiologicalFeaturePayloads({
    features: [
      { fileIdx: 0, stable_svg_id: 'duplicate-biological', product: 'first' },
      { fileIdx: 0, stable_svg_id: 'duplicate-biological', product: 'second' }
    ],
    orthogroups: []
  }, []),
  []
);

const fiveMemberGroup = {
  id: 'og_1',
  members: Array.from({ length: 5 }, (_, recordIndex) => {
    const stableFeatureSvgId = recordIndex === 2
      ? 'canonical-stable-id'
      : `stable-${recordIndex}`;
    const renderedFeatureSvgId = recordIndex === 2
      ? 'reversed-display-stable-id_record_3'
      : `${stableFeatureSvgId}_record_${recordIndex + 1}`;
    return {
      recordIndex,
      stableFeatureSvgId,
      featureSvgId: stableFeatureSvgId,
      renderedFeatureSvgId
    };
  })
};
const fiveGeneratedFeatures = enrichFeaturesWithOrthogroups(
  fiveMemberGroup.members.map((member) => ({
    fileIdx: member.recordIndex,
    svg_id: member.renderedFeatureSvgId,
    stable_svg_id: member.stableFeatureSvgId,
    rendered_svg_id: member.renderedFeatureSvgId
  })),
  buildOrthogroupFeatureIndex([fiveMemberGroup])
);
const orthogroupSearchResult = runFeatureSearch({
  features: fiveGeneratedFeatures,
  renderedFeatureIds: new Set(fiveGeneratedFeatures.map((feature) => feature.svg_id)),
  query: 'og_1',
  field: 'orthogroup',
  useRegex: false,
  orthogroups: [fiveMemberGroup]
});
assert.deepEqual(orthogroupSearchResult.matches, [
  'stable-0_record_1',
  'stable-1_record_2',
  'reversed-display-stable-id_record_3',
  'stable-3_record_4',
  'stable-4_record_5'
]);
assert.match(
  orthogroupSearchResult.matchDetails['reversed-display-stable-id_record_3'][0].value,
  /og_1/
);

const searchableFeature = {
  svg_id: 'searchable-feature',
  displayLabel: 'Edited beta subunit',
  recordId: 'NC_000913.3',
  type: 'CDS',
  location_parts: [{ display: 'join(11..40, 80..120)' }],
  qualifiers: {
    product: ['DNA-directed RNA polymerase subunit beta'],
    note: ['core enzyme'],
    translation: ['MPEPTIDE']
  },
  nucleotideSequence: 'ATGGCN'
};
const search = (query, field = 'all', popupMode = 'rich') => runFeatureSearch({
  features: [searchableFeature],
  renderedFeatureIds: new Set(['searchable-feature']),
  query,
  field,
  popupMode
});
for (const [query, field] of [
  ['polymerase', 'all'],
  ['NC_000913.3', 'record-id'],
  ['join(11..40, 80..120)', 'location'],
  ['ATGGNN', 'nucleotide'],
  ['MPEPTIDE', 'amino-acid']
]) {
  assert.deepEqual(search(query, field).matches, ['searchable-feature']);
}
assert.deepEqual(search('core enzyme', 'all', 'simple').matches, []);

const publicProteinOnlyFeature = { protein_id: 'CAG34720.1', fileIdx: 2 };
assert.equal(
  enrichFeatureWithOrthogroup(orthogroupIndex, publicProteinOnlyFeature),
  publicProteinOnlyFeature
);
const staleFeature = {
  svg_id: 'not-a-group-member',
  orthogroupId: 'og_stale',
  orthogroup_member: { feature_svg_id: 'not-a-group-member' }
};
assert.deepEqual(
  enrichFeatureWithOrthogroup(orthogroupIndex, staleFeature),
  { svg_id: 'not-a-group-member' }
);
const state = {
  orthogroups: ref([group]),
  orthogroupNameOverrides: {},
  orthogroupDescriptionOverrides: {},
  selectedOrthogroupId: ref('og_1'),
  orthogroupSearch: ref(''),
  orthogroupSortMode: ref('id'),
  selectedOrthogroupAlignmentFeature: ref(''),
  svgContainer: ref({ querySelector: () => svg }),
  showRightDrawer: ref(false),
  rightDrawerTab: ref('orthogroups'),
  showFeaturePanel: ref(false),
  showLegendPanel: ref(false),
  linearSeqs: [
    { definition: 'record-a' },
    { definition: 'record-b' },
    { definition: 'record-c' }
  ],
  extractedFeatures: ref(generatedFeatures),
  biologicalFeatures: ref([
    {
      fileIdx: 0,
      record_id: 'record-a',
      svg_id: 'stable-x',
      stable_svg_id: 'stable-x',
      locus_tag: 'LOC_A',
      nucleotide_sequence: 'AAAA',
      amino_acid_sequence: 'MK'
    },
    {
      fileIdx: 1,
      record_id: 'record-b',
      svg_id: 'stable-x',
      stable_svg_id: 'stable-x',
      locus_tag: 'LOC_B',
      nucleotide_sequence: 'CCCC',
      amino_acid_sequence: 'MP'
    },
    {
      fileIdx: 2,
      record_id: 'record-c',
      recordKey: 'record-key-c',
      biologicalFeatureId: 'biological-c',
      svg_id: 'canonical-stable-id',
      stable_svg_id: 'canonical-stable-id',
      qualifiers: { protein_id: ['CAG34720.1'] },
      nucleotide_sequence: 'GGGG',
      amino_acid_sequence: 'MQ'
    }
  ])
};

const editor = createOrthogroupEditor({ state, runAnalysis: null });
const aliasedEditorGroup = {
  orthogroupId: 'og_editor_alias',
  orthogroup_id: 'og_editor_alias',
  members: group.members
};
state.orthogroups.value = [aliasedEditorGroup];
assert.equal(editor.getOrthogroupById('og_editor_alias'), aliasedEditorGroup);
assert.equal(editor.resolveOrthogroupName(aliasedEditorGroup), 'og_editor_alias');
assert.equal(editor.orthogroupCount.value, 1);
state.orthogroups.value = [{
  id: 'og_editor_conflict',
  orthogroupId: 'og_editor_other',
  members: group.members
}];
assert.equal(editor.getOrthogroupById('og_editor_conflict'), null);
assert.equal(editor.getOrthogroupById('og_editor_other'), null);
assert.equal(editor.orthogroupCount.value, 0);
assert.deepEqual(editor.filteredOrthogroups.value, []);
state.orthogroups.value = [
  { id: 'og_editor_duplicate', members: group.members },
  { orthogroup_id: 'og_editor_duplicate', members: group.members }
];
assert.equal(editor.getOrthogroupById('og_editor_duplicate'), null);
assert.deepEqual(editor.getEnrichedOrthogroupMembers('og_editor_duplicate'), []);
assert.equal(editor.orthogroupCount.value, 0);
state.orthogroups.value = [group];
editor.openOrthogroupInDrawer(generatedFeatures[1].orthogroupId);
assert.equal(state.selectedOrthogroupId.value, 'og_1');
assert.equal(state.rightDrawerTab.value, 'orthogroups');
assert.equal(state.showRightDrawer.value, true);
const members = editor.getEnrichedOrthogroupMembers(group);
assert.deepEqual(members.map((member) => member.nucleotideSequence), ['AAAA', 'CCCC', 'GGGG']);
assert.deepEqual(members.map((member) => member.displayProteinId), ['LOC_A', 'LOC_B', 'CAG34720.1']);
const wrongRecordMember = editor.getEnrichedOrthogroupMembers({
  id: 'og_wrong_record',
  members: [{ recordIndex: 1, stableFeatureSvgId: 'canonical-stable-id' }]
})[0];
assert.equal(wrongRecordMember.nucleotideSequence, undefined);
assert.equal(wrongRecordMember.aminoAcidSequence, undefined);
const canonicalOnlyMember = editor.getEnrichedOrthogroupMembers({
  id: 'og_canonical_only',
  members: [{ recordKey: 'record-key-c', biologicalFeatureId: 'biological-c' }]
})[0];
assert.equal(canonicalOnlyMember.nucleotideSequence, 'GGGG');
assert.equal(canonicalOnlyMember.aminoAcidSequence, 'MQ');
for (const invalidMember of [
  { stableFeatureSvgId: 'canonical-stable-id' },
  {
    recordIndex: 2,
    stableFeatureSvgId: 'canonical-stable-id',
    renderedFeatureSvgId: 'unknown-rendered'
  },
  { recordIndex: -1, stableFeatureSvgId: 'canonical-stable-id' },
  { recordIndex: '2.5', stableFeatureSvgId: 'canonical-stable-id' },
  {
    recordIndex: 2,
    record_index: 1,
    stableFeatureSvgId: 'canonical-stable-id'
  }
]) {
  const unresolved = editor.getEnrichedOrthogroupMembers({
    id: 'og_invalid_identity',
    members: [invalidMember]
  })[0];
  assert.equal(unresolved.nucleotideSequence, undefined);
  assert.equal(unresolved.aminoAcidSequence, undefined);
}
const originalBiologicalFeatures = state.biologicalFeatures.value;
state.biologicalFeatures.value = [
  ...originalBiologicalFeatures,
  { ...originalBiologicalFeatures[0] }
];
const duplicateSequenceMember = editor.getEnrichedOrthogroupMembers({
  id: 'og_duplicate_sequence_identity',
  members: [{ recordIndex: 0, stableFeatureSvgId: 'stable-x' }]
})[0];
assert.equal(duplicateSequenceMember.nucleotideSequence, undefined);
assert.equal(duplicateSequenceMember.aminoAcidSequence, undefined);
state.biologicalFeatures.value = originalBiologicalFeatures;
assert.equal(editor.getOrthogroupSequenceCount(group, 'nt'), 3);
await editor.copyOrthogroupSequences(group, 'nt');
assert.match(copiedText, /AAAA/);
assert.match(copiedText, /CCCC/);
assert.match(copiedText, /GGGG/);
await editor.copyOrthogroupSequences(group, 'aa');
assert.match(copiedText, />LOC_A\b/);
assert.match(copiedText, />LOC_B\b/);
assert.match(copiedText, />CAG34720\.1\b/);
assert.doesNotMatch(copiedText, /h_[a-z2-7]{26}/);

const membersWithoutStableIds = editor.getEnrichedOrthogroupMembers({
  id: 'og_no_stable_ids',
  members: [
    {
      recordId: 'record-c',
      proteinId: runtimeHandle,
      sourceProteinId: 'WP_VISIBLE.1',
      label: featureAnalysisId
    },
    {
      recordId: 'record-d',
      proteinId: runtimeHandle,
      displayProteinId: featureAnalysisId,
      sourceProteinId: unsupportedHistoricalTransportId,
      label: runtimeHandle
    }
  ]
});
assert.deepEqual(
  membersWithoutStableIds.map((member) => member.displayProteinId),
  ['WP_VISIBLE.1', '']
);

const downloadedFilenames = [];
globalThis.document = {
  createElement: () => ({
    addEventListener() {},
    click() { downloadedFilenames.push(this.download); }
  }),
  body: {
    appendChild() {},
    removeChild() {}
  }
};
Object.defineProperty(globalThis.URL, 'createObjectURL', {
  configurable: true,
  value: () => 'blob:orthogroup-test'
});
Object.defineProperty(globalThis.URL, 'revokeObjectURL', {
  configurable: true,
  value: () => {}
});
editor.downloadOrthogroupMemberSequence(
  {
    proteinId: runtimeHandle,
    sourceProteinId: 'WP_DOWNLOAD.1',
    aminoAcidSequence: 'MK'
  },
  'aa',
  { id: 'og_download' }
);
editor.downloadOrthogroupMemberSequence(
  {
    proteinId: runtimeHandle,
    stableFeatureSvgId: featureAnalysisId,
    aminoAcidSequence: 'MK'
  },
  'aa',
  { id: 'og_download' }
);
assert.deepEqual(downloadedFilenames, [
  'og_download_WP_DOWNLOAD.1_aa.faa',
  'og_download_member_aa.faa'
]);
assert.doesNotMatch(downloadedFilenames.join('\n'), /(?:h_[a-z2-7]{26}|f_[0-9a-f]{64})/);

const indexSource = await readFile(join(repoRoot, 'gbdraw/web/index.html'), 'utf8');
const appSetupSource = await readFile(join(repoRoot, 'gbdraw/web/js/app/app-setup.js'), 'utf8');
assert.match(
  appSetupSource,
  /resolveUniqueOrthogroupMemberForFeature\(cf\?\.feat, members\)/
);
assert.doesNotMatch(appSetupSource, /const currentFeatureIds = new Set/);
assert.doesNotMatch(
  indexSource,
  /member\.displayProteinId\s*\|\|\s*member\.(?:sourceProteinId|locusTag|label)/
);
assert.doesNotMatch(
  indexSource,
  /currentMember\.displayProteinId\s*\|\|\s*clickedOrthogroupDetail\.currentMember\.sourceProteinId/
);
const standaloneAssetSource = await readFile(
  join(repoRoot, 'gbdraw/web/js/services/standalone-interactivity-assets.js'),
  'utf8'
);
const embeddedCatalogRecordIndex = new Function(
  `${embeddedFunctionSource('catalogRecordIndex')}; return catalogRecordIndex;`
)();
assert.equal(embeddedCatalogRecordIndex({ recordKeys: [] }, { record_idx: '2' }), 2);
for (const invalidIndex of ['-1', '1.5', 'bogus', '9007199254740992']) {
  assert.equal(
    embeddedCatalogRecordIndex({ recordKeys: [] }, { record_idx: invalidIndex }),
    null
  );
}
const embeddedCatalogSequenceSourceStrings = new Function(
  `${embeddedFunctionSource('catalogSequenceSourceStrings')}; return catalogSequenceSourceStrings;`
)();
assert.deepEqual(
  embeddedCatalogSequenceSourceStrings({ sequenceSources: [{ sequence: 'ACGT' }] }),
  ['ACGT']
);
assert.throws(
  () => embeddedCatalogSequenceSourceStrings({ sequenceSources: [{ sequence: 'AC GT' }] }),
  /Invalid catalog sequence source/
);
const embeddedEscapeRegExpClassChar = new Function(
  `${embeddedFunctionSource('escapeRegExpClassChar')}; return escapeRegExpClassChar;`
)();
for (const character of ['\\', ']', '[', '^', '-']) {
  assert.equal(embeddedEscapeRegExpClassChar(character), `\\${character}`);
}

const embeddedMatchIdFunctions = new Function(`
  var ambiguousMatchElementIds = new Set();
  ${embeddedFunctionSource('matchAttr')}
  ${embeddedFunctionSource('getElementMatchId')}
  return { getElementMatchId, ambiguousMatchElementIds };
`)();
const matchElement = (attributes) => ({
  getAttribute(name) { return attributes[name] ?? null; }
});
assert.equal(
  embeddedMatchIdFunctions.getElementMatchId(matchElement({
    'data-gbdraw-match-id': 'match-1',
    'data-gbdraw-pairwise-match-id': 'match-1'
  })),
  'match-1'
);
assert.equal(
  embeddedMatchIdFunctions.getElementMatchId(matchElement({
    'data-gbdraw-match-id': 'match-1',
    'data-gbdraw-pairwise-match-id': 'match-2'
  })),
  ''
);
embeddedMatchIdFunctions.ambiguousMatchElementIds.add('match-1');
assert.equal(
  embeddedMatchIdFunctions.getElementMatchId(matchElement({
    'data-gbdraw-match-id': 'match-1',
    'data-gbdraw-pairwise-match-id': 'match-1'
  })),
  ''
);

const embeddedIdentityHelpers = [
  embeddedFunctionSource('consistentTextIdentity'),
  embeddedFunctionSource('nonnegativeIntegerIdentity'),
  embeddedFunctionSource('matchSourceAliases'),
  embeddedFunctionSource('resolveEmbeddedMatchSource')
].join('\n');
const createEmbeddedMatchSourceResolver = (sequenceSources) => new Function(
  'sequenceSources',
  `${embeddedIdentityHelpers}; return resolveEmbeddedMatchSource;`
)(sequenceSources);
const embeddedLinearSource = {
  origin: 'linear-record',
  recordIndex: 0,
  recordId: 'record-a',
  aliases: ['alias-a'],
  sequence: 'ACGT'
};
let embeddedMatchSourceResolver = createEmbeddedMatchSourceResolver([embeddedLinearSource]);
assert.equal(
  embeddedMatchSourceResolver({
    match_kind: 'pairwise',
    query_record_index: 0,
    query_record_id: 'alias-a'
  }, 'query').source,
  embeddedLinearSource
);
assert.equal(
  embeddedMatchSourceResolver({
    match_kind: 'pairwise',
    query_record_index: 0,
    query_record_id: 'record-b'
  }, 'query').source,
  null
);
assert.equal(
  embeddedMatchSourceResolver({
    match_kind: 'pairwise',
    query_record_index: 0,
    query_record_id: 'record-a',
    queryRecordId: 'record-b'
  }, 'query').source,
  null
);
assert.equal(
  embeddedMatchSourceResolver({
    match_kind: 'pairwise',
    query_record_index: '0.5',
    query_record_id: 'record-a'
  }, 'query').source,
  null
);
embeddedMatchSourceResolver = createEmbeddedMatchSourceResolver([
  embeddedLinearSource,
  { ...embeddedLinearSource, recordId: 'record-b', aliases: [] }
]);
assert.equal(
  embeddedMatchSourceResolver({
    match_kind: 'pairwise',
    query_record_index: 0
  }, 'query').source,
  null
);
embeddedMatchSourceResolver = createEmbeddedMatchSourceResolver([
  embeddedLinearSource,
  { ...embeddedLinearSource, recordId: 'record-b', aliases: ['record-a'] }
]);
assert.equal(
  embeddedMatchSourceResolver({
    match_kind: 'pairwise',
    query_record_index: 0,
    query_record_id: 'record-a'
  }, 'query').source,
  null
);

const scopedFeature = {};
const otherScopedFeature = {};
const embeddedGroupMemberResolver = new Function(
  'biologicalFeatureForMember',
  `${embeddedFunctionSource('getFeatureOrthogroupMember')}; return getFeatureOrthogroupMember;`
)((value) => value?.resolvedFeature || null);
assert.equal(
  embeddedGroupMemberResolver(
    { resolvedFeature: scopedFeature },
    { members: [{ resolvedFeature: scopedFeature }, { resolvedFeature: otherScopedFeature }] }
  )?.resolvedFeature,
  scopedFeature
);
assert.equal(
  embeddedGroupMemberResolver(
    { resolvedFeature: scopedFeature },
    { members: [{ resolvedFeature: scopedFeature }, { resolvedFeature: scopedFeature }] }
  ),
  null
);
assert.doesNotMatch(standaloneAssetSource, /\bbiologicalFeaturesById\b/);
assert.match(standaloneAssetSource, /biologicalFeaturesByCanonicalId/);
assert.match(
  standaloneAssetSource,
  /!hasSourceIdentity && !recordKey\.value/
);
assert.match(standaloneAssetSource, /resultIndexIdentity/);
assert.match(standaloneAssetSource, /Match sequence source identity is invalid/);
assert.doesNotMatch(standaloneAssetSource, /recordIndexText \? Number\(recordIndexText\)/);
assert.doesNotMatch(
  standaloneAssetSource,
  /if \(hasBiologicalFeatureCatalog\) return null;/
);
assert.match(standaloneAssetSource, /rendered\.length === 1/);
assert.match(standaloneAssetSource, /var ambiguousMatchIds = new Set\(\)/);
assert.match(standaloneAssetSource, /matchesById\.delete\(id\)/);
assert.match(standaloneAssetSource, /var ambiguousMatchElementIds = new Set\(\)/);
assert.match(
  standaloneAssetSource,
  /consistentTextIdentity\(feature, \['orthogroup_id', 'orthogroupId'\]\)/
);
assert.match(
  standaloneAssetSource,
  /getFeatureOrthogroupMember\(renderedFeature, group\)/
);
assert.doesNotMatch(
  standaloneAssetSource,
  /return explicit\.value \|\| memberFeatureSvgId\(memberOrRow\)/
);
assert.match(
  standaloneAssetSource,
  /function featureForMember\(memberOrRow\) \{\s*return biologicalFeatureForMember\(memberOrRow\);/
);
assert.match(
  STANDALONE_INTERACTIVE_SCRIPT,
  /if \(!\/\^\\d\+\$\/\.test\(text\)/
);
assert.match(STANDALONE_INTERACTIVE_SCRIPT, /\|\| \/\\s\/\.test\(source\.sequence\)/);

editor.highlightOrthogroupById('og_1');
assert.equal(visibleElement.getAttribute('stroke'), '#2563eb');
assert.equal(visibleElement.getAttribute('stroke-width'), '2.4');
assert.equal(svg.featureElements.has('stable-x_record_2'), false);

state.orthogroups.value = [{
  id: 'og_hidden_only',
  members: [
    { recordIndex: 1, recordId: 'record-b', featureSvgId: 'stable-x', proteinId: 'hidden' }
  ]
}];
editor.highlightOrthogroupById('og_hidden_only');
assert.equal(visibleElement.getAttribute('stroke'), '#111111');
assert.equal(visibleElement.getAttribute('stroke-width'), '0.5');

console.log('orthogroup stable identity tests passed');
