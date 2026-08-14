import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-session-feature-metadata-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });

const copyModule = async (sourceRelative, targetRelative) => {
  await writeFile(
    join(tempDir, targetRelative),
    await readFile(join(repoRoot, sourceRelative), 'utf8'),
    'utf8'
  );
};

await copyModule('gbdraw/web/js/app/session-feature-metadata.js', 'app/session-feature-metadata.js');
await copyModule('gbdraw/web/js/app/feature-metadata-extraction.js', 'app/feature-metadata-extraction.js');
await copyModule('gbdraw/web/js/app/losat-normalization.js', 'app/losat-normalization.js');
await copyModule('gbdraw/web/js/services/diagram-generation.js', 'services/diagram-generation.js');
await copyModule('gbdraw/web/js/services/diagram-worker-protocol.js', 'services/diagram-worker-protocol.js');
await copyModule('gbdraw/web/js/services/error-normalization.js', 'services/error-normalization.js');
await copyModule('gbdraw/web/js/services/file-content-cache.js', 'services/file-content-cache.js');
await copyModule('gbdraw/web/js/services/depth-file-codec.js', 'services/depth-file-codec.js');
await copyModule('gbdraw/web/js/services/json-clone.js', 'services/json-clone.js');
await copyModule('gbdraw/web/js/services/feature-identity.js', 'services/feature-identity.js');
await copyModule(
  'gbdraw/web/js/services/orthogroup-feature-metadata.js',
  'services/orthogroup-feature-metadata.js'
);
await copyModule('gbdraw/web/js/services/pyodide-assets.js', 'services/pyodide-assets.js');
await copyModule('gbdraw/web/js/services/runtime-capabilities.js', 'services/runtime-capabilities.js');
await copyModule('gbdraw/web/js/services/runtime-test-hooks.js', 'services/runtime-test-hooks.js');
await copyModule(
  'gbdraw/web/js/services/session-resource-backing.js',
  'services/session-resource-backing.js'
);
await copyModule('gbdraw/web/js/services/session-feature-metadata.js', 'services/session-feature-metadata.js');
await copyModule('gbdraw/web/js/services/svg-result-ingestion.js', 'services/svg-result-ingestion.js');
await copyModule('gbdraw/web/js/services/svg-result-normalization.js', 'services/svg-result-normalization.js');
await copyModule('gbdraw/web/js/services/svg-sanitization.js', 'services/svg-sanitization.js');
await copyModule('gbdraw/web/js/services/svg-serialization.js', 'services/svg-serialization.js');
await copyModule('gbdraw/web/js/config.js', 'config.js');

const {
  alignRecoveredFeatureIdsToRenderedSvg,
  buildSessionFeatureRecoveryPlan,
  classifyFeatureMetadataState,
  migrateFeatureOverrideState,
  normalizeRecordIndex
} = await import(pathToFileURL(join(tempDir, 'app', 'session-feature-metadata.js')));
const { collectRenderedFeatureIdentitiesFromSvgRoot } = await import(
  pathToFileURL(join(tempDir, 'services', 'session-feature-metadata.js'))
);
const { ingestSvgResult } = await import(
  pathToFileURL(join(tempDir, 'services', 'svg-result-ingestion.js'))
);
const { extractFeatureMetadataForPreview } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-metadata-extraction.js'))
);
const { normalizeGenerationResponse } = await import(pathToFileURL(join(tempDir, 'services', 'diagram-generation.js')));

{
  const gffFile = { name: 'record.gff3' };
  const fastaFile = { name: 'record.fna' };
  let extractionRequest = null;
  const metadata = await extractFeatureMetadataForPreview({
    mode: 'circular',
    cInputType: 'gff',
    lInputType: 'gb',
    circularFile: gffFile,
    circularFastaFile: fastaFile,
    readFeatureExtractionDataImpl: async (request) => {
      extractionRequest = request;
      return {
        features: [{ id: 'feature-a', svg_id: 'stable-a', type: 'CDS' }],
        record_ids: ['record-a'],
        selector_safety_scope: []
      };
    }
  });
  assert.equal(extractionRequest.path, '/input.gff');
  assert.equal(extractionRequest.file, gffFile);
  assert.equal(extractionRequest.format, 'gff');
  assert.equal(extractionRequest.fastaPath, '/input.fasta');
  assert.equal(extractionRequest.fastaFile, fastaFile);
  assert.equal(extractionRequest.mode, 'circular');
  assert.equal(metadata.extractedFeatures.length, 1);
}

{
  const requests = [];
  const metadata = await extractFeatureMetadataForPreview({
    mode: 'linear',
    cInputType: 'gb',
    lInputType: 'gff',
    linearSeqs: [
      { gff: { name: 'a.gff3' }, fasta: { name: 'a.fna' } },
      { gff: { name: 'b.gff3' }, fasta: { name: 'b.fna' } }
    ],
    readFeatureExtractionDataImpl: async (request) => {
      requests.push(request);
      const index = requests.length - 1;
      return {
        features: [{ id: `feature-${index}`, svg_id: `stable-${index}`, record_id: `record-${index}` }],
        biological_features: [
          { id: `feature-${index}`, svg_id: `stable-${index}`, record_id: `record-${index}` },
          { id: `hidden-${index}`, svg_id: `hidden-stable-${index}`, record_id: `record-${index}` }
        ],
        record_ids: [`record-${index}`],
        selector_safety_scope: []
      };
    }
  });
  assert.deepEqual(requests.map((request) => request.path), ['/seq_0.gff', '/seq_1.gff']);
  assert.deepEqual(requests.map((request) => request.fastaPath), ['/seq_0.fasta', '/seq_1.fasta']);
  assert.deepEqual(
    metadata.extractedFeatures.map((feature) => feature.svg_id),
    ['stable-0', 'stable-1']
  );
  assert.deepEqual(
    metadata.biologicalFeatures.map((feature) => feature.svg_id),
    ['stable-0', 'hidden-stable-0', 'stable-1', 'hidden-stable-1']
  );
  assert.equal(requests.every((request) => request.includeBiologicalFeatures === true), true);
}

{
  const legacyResults = [{ name: 'out.svg', content: '<svg />' }];
  assert.deepEqual(normalizeGenerationResponse(legacyResults), { results: legacyResults, metadata: {} });
  const metadata = { trackSlotGeometry: { schema: 1, mode: 'linear', records: [] } };
  assert.deepEqual(
    normalizeGenerationResponse({ results: legacyResults, metadata }),
    { results: legacyResults, metadata }
  );
  const errorPayload = { error: { type: 'OutputError', message: 'No output files generated.' } };
  assert.deepEqual(normalizeGenerationResponse(errorPayload), { results: errorPayload, metadata: {} });
}

const parseAttributes = (source) => {
  const attrs = {};
  for (const match of source.matchAll(/([:\w-]+)\s*=\s*(["'])(.*?)\2/g)) {
    attrs[match[1]] = match[3];
  }
  return attrs;
};

class FakeElement {
  constructor(tagName, attrs) {
    this.tagName = tagName;
    this.attrs = attrs;
  }

  getAttribute(name) {
    return Object.hasOwn(this.attrs, name) ? this.attrs[name] : null;
  }

  setAttribute(name, value) {
    this.attrs[name] = String(value);
  }

  removeAttribute(name) {
    delete this.attrs[name];
  }
}

class FakeDocument {
  constructor(elements) {
    this.elements = elements;
  }

  querySelector(selector) {
    return selector === 'parsererror' ? null : null;
  }

  querySelectorAll(selector) {
    if (selector === '[id]') {
      return this.elements.filter((element) => element.getAttribute('id'));
    }
    if (!String(selector).includes('data-gbdraw-feature-id') && !String(selector).includes('[id^="f"]')) {
      return [];
    }
    return this.elements.filter((element) => {
      if (element.getAttribute('data-gbdraw-feature-id')) return true;
      const id = element.getAttribute('id') || '';
      return ['path', 'polygon', 'rect'].includes(element.tagName) && id.startsWith('f');
    });
  }
}

globalThis.DOMParser = class {
  parseFromString(source) {
    const elements = [];
    for (const match of String(source || '').matchAll(/<(path|polygon|rect|g)\b([^>]*)>/g)) {
      elements.push(new FakeElement(match[1], parseAttributes(match[2])));
    }
    return new FakeDocument(elements);
  }
};

const parsedFeatureRoot = (source) => {
  const parsed = new globalThis.DOMParser().parseFromString(source, 'image/svg+xml');
  return {
    localName: 'svg',
    content: source,
    cloneNode() { return this; },
    getAttribute() { return null; },
    setAttribute() {},
    removeAttribute() {},
    querySelectorAll: (selector) => parsed.querySelectorAll(selector)
  };
};

const collectRenderedFeatureIdentitiesFromSvg = (source) =>
  collectRenderedFeatureIdentitiesFromSvgRoot(parsedFeatureRoot(source));

class IngressDomParser {
  parseFromString(source) {
    return {
      documentElement: parsedFeatureRoot(source),
      querySelector: () => null
    };
  }
}

globalThis.XMLSerializer = class {
  serializeToString(root) { return root.content; }
};

const committedResults = (content) => [ingestSvgResult(
  { name: 'preview.svg', content },
  {
    sanitizer: { sanitize: (value) => value },
    parser: IngressDomParser
  }
)];

const svgWithFeature = ({
  renderedId = 'rendered-a',
  stableId = 'stable-a',
  recordIndex = 0,
  recordId = 'record-a',
  elementId = renderedId
} = {}) => `
  <svg>
    <path
      id="${elementId}"
      data-gbdraw-feature-id="${renderedId}"
      data-gbdraw-stable-feature-id="${stableId}"
      data-gbdraw-record-index="${recordIndex}"
      data-gbdraw-record-id="${recordId}"
    />
  </svg>`;

{
  const identities = collectRenderedFeatureIdentitiesFromSvg(svgWithFeature({
    renderedId: 'rendered-a__part2',
    stableId: 'stable-a',
    recordIndex: 3,
    recordId: 'record-A',
    elementId: 'dom-feature-a'
  }));
  const identity = identities.byRenderedId.get('rendered-a');
  assert.deepEqual(identity, {
    renderedId: 'rendered-a',
    stableId: 'stable-a',
    recordIndex: 3,
    recordId: 'record-A',
    elementId: 'dom-feature-a'
  });
  assert.equal(identities.renderedIds.has('rendered-a'), true);
  assert.equal(identities.byStableId.get('stable-a')[0], identity);
  assert.equal(identities.totalRenderedCount, 1);
}

{
  const ambiguousSvg = `
    <svg>
      <path data-gbdraw-feature-id="rendered-r" data-gbdraw-stable-feature-id="stable-a" data-gbdraw-record-index="0" data-gbdraw-record-id="record-a" />
      <path data-gbdraw-feature-id="rendered-r" data-gbdraw-stable-feature-id="stable-b" data-gbdraw-record-index="1" data-gbdraw-record-id="record-b" />
    </svg>
  `;
  const identities = collectRenderedFeatureIdentitiesFromSvg(ambiguousSvg);
  assert.equal(identities.renderedIds.has('rendered-r'), true);
  assert.equal(identities.ambiguousRenderedIds.has('rendered-r'), true);
  assert.equal(identities.byRenderedId.has('rendered-r'), false);
  assert.equal(identities.byStableId.has('stable-a'), false);
  assert.equal(identities.byStableId.has('stable-b'), false);

  const recovered = {
    id: 'feature-a',
    svg_id: 'rendered-r',
    stable_svg_id: 'stable-a',
    record_idx: 0,
    record_id: 'record-a'
  };
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [recovered],
    renderedIdentities: identities
  });
  assert.equal(aligned.exactCount, 0);
  assert.equal(aligned.alignedCount, 0);
  assert.equal(aligned.unresolvedCount, 1);
  assert.equal(aligned.features[0], recovered);
  assert.equal(aligned.features[0].rendered_svg_id, undefined);
  assert.equal(aligned.svgIdMap['rendered-r'], undefined);

  const state = classifyFeatureMetadataState({
    results: committedResults(ambiguousSvg),
    selectedResultIndex: 0,
    extractedFeatures: [recovered]
  });
  assert.equal(state.state, 'stale');
  assert.equal(state.renderedCount, 1);
  assert.equal(state.totalRenderedCount, 1);
  assert.equal(state.exactMatchingCount, 0);
  assert.equal(state.aliasMatchingCount, 0);
  assert.equal(state.missingExactCount, 1);

  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'circular',
      cInputType: 'gff',
      selectedResultIndex: 0,
      results: committedResults(ambiguousSvg),
      featureState: { extractedFeatures: [recovered], biologicalFeatures: [] },
      editorState: {}
    },
    featureVisibilityTsv: ''
  });
  assert.equal(plan.status, 'unrecoverable');
  assert.equal(plan.validation.state, 'stale');
}

{
  const identities = collectRenderedFeatureIdentitiesFromSvg(`
    <svg>
      <path id="part-1" data-gbdraw-feature-id="rendered-r__part1" data-gbdraw-stable-feature-id="stable-a" data-gbdraw-record-index="0" data-gbdraw-record-id="record-a" />
      <path id="part-2" data-gbdraw-feature-id="rendered-r__part2" data-gbdraw-stable-feature-id="stable-a" data-gbdraw-record-index="0" data-gbdraw-record-id="record-a" />
    </svg>
  `);
  assert.equal(identities.ambiguousRenderedIds.has('rendered-r'), false);
  assert.equal(identities.byRenderedId.get('rendered-r')?.stableId, 'stable-a');
}

{
  const state = classifyFeatureMetadataState({
    results: committedResults(svgWithFeature({ renderedId: 'rendered-a' })),
    selectedResultIndex: 0,
    extractedFeatures: [{ id: 'feature-a', svg_id: 'rendered-a' }]
  });
  assert.equal(state.state, 'ready');
  assert.equal(state.exactMatchingCount, 1);
  assert.equal(state.missingExactCount, 0);
}

{
  const state = classifyFeatureMetadataState({
    results: committedResults(svgWithFeature({ renderedId: 'rendered-a', stableId: 'stable-a' })),
    selectedResultIndex: 0,
    extractedFeatures: [{ id: 'feature-a', svg_id: 'stable-a', stable_svg_id: 'stable-a' }]
  });
  assert.equal(state.state, 'alignable');
  assert.equal(state.exactMatchingCount, 0);
  assert.equal(state.aliasMatchingCount, 1);
  assert.equal(state.missingExactCount, 1);
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(
    svgWithFeature({ renderedId: 'rendered-a', stableId: 'stable-a' })
  );
  const feature = { id: 'feature-a', svg_id: 'stable-a', stable_svg_id: 'stable-a' };
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [feature],
    renderedIdentities
  });
  assert.equal(aligned.changedCount, 0);
  assert.equal(aligned.alignedCount, 0);
  assert.equal(aligned.unresolvedCount, 1);
  assert.equal(aligned.features[0], feature);
  assert.equal(aligned.svgIdMap['stable-a'], undefined);
}

assert.equal(normalizeRecordIndex('4'), 4);
for (const invalidRecordIndex of [
  -1,
  '4.5',
  'bogus',
  '9007199254740992',
  true
]) {
  assert.equal(normalizeRecordIndex(invalidRecordIndex), null);
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(`
    <svg>
      <path id="r1" data-gbdraw-feature-id="stable-x_record_1" data-gbdraw-stable-feature-id="stable-x" data-gbdraw-record-index="0" />
      <path id="r2" data-gbdraw-feature-id="stable-x_record_2" data-gbdraw-stable-feature-id="stable-x" data-gbdraw-record-index="1" />
    </svg>
  `);
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [
      { id: 'feature-a', svg_id: 'stable-x', stable_svg_id: 'stable-x', fileIdx: 0 },
      { id: 'feature-b', svg_id: 'stable-x', stable_svg_id: 'stable-x', fileIdx: 1 }
    ],
    renderedIdentities,
    writeRenderedSvgIdToSvgId: true
  });
  assert.equal(aligned.ambiguousCount, 0);
  assert.equal(aligned.unresolvedCount, 0);
  assert.equal(aligned.features[0].svg_id, 'stable-x_record_1');
  assert.equal(aligned.features[1].svg_id, 'stable-x_record_2');
  assert.equal(aligned.features[0].stable_svg_id, 'stable-x');
  assert.equal(aligned.features[1].stable_svg_id, 'stable-x');
  assert.equal(aligned.features[0].rendered_svg_id, 'stable-x_record_1');
  assert.equal(aligned.features[1].rendered_svg_id, 'stable-x_record_2');
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(`
    <svg>
      <path data-gbdraw-feature-id="stable-y_record_1" data-gbdraw-stable-feature-id="stable-y" />
      <path data-gbdraw-feature-id="stable-y_record_2" data-gbdraw-stable-feature-id="stable-y" />
    </svg>
  `);
  const feature = { id: 'feature-a', svg_id: 'stable-y', stable_svg_id: 'stable-y' };
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [feature],
    renderedIdentities
  });
  assert.equal(aligned.ambiguousCount, 1);
  assert.equal(aligned.changedCount, 0);
  assert.equal(aligned.features[0], feature);
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(
    svgWithFeature({ renderedId: 'rendered-z', stableId: 'stable-z' })
  );
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [{ id: 'feature-a', svg_id: 'missing-a', stable_svg_id: 'missing-a' }],
    renderedIdentities
  });
  assert.equal(aligned.unresolvedCount, 1);
  assert.equal(aligned.changedCount, 0);
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(
    svgWithFeature({ renderedId: 'rendered-a', stableId: 'stable-a' })
  );
  const feature = { id: 'feature-a', svg_id: 'rendered-a', stable_svg_id: 'stable-a' };
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [feature],
    renderedIdentities
  });
  assert.equal(aligned.exactCount, 1);
  assert.notEqual(aligned.features[0], feature);
  assert.equal(aligned.features[0].svg_id, 'rendered-a');
  assert.equal(aligned.features[0].rendered_svg_id, 'rendered-a');
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(
    svgWithFeature({ renderedId: 'rendered-a', stableId: 'stable-a', recordIndex: 0 })
  );
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [
      {
        id: 'wrong-stable',
        svg_id: 'rendered-a',
        stable_svg_id: 'stable-wrong',
        record_idx: 0
      },
      {
        id: 'wrong-record',
        svg_id: 'rendered-a',
        stable_svg_id: 'stable-a',
        record_idx: 1
      },
      {
        id: 'conflicting-stable-aliases',
        svg_id: 'rendered-a',
        stable_svg_id: 'stable-a',
        stableFeatureId: 'stable-wrong',
        record_idx: 0
      },
      {
        id: 'conflicting-record-aliases',
        svg_id: 'rendered-a',
        stable_svg_id: 'stable-a',
        record_idx: 0,
        recordIndex: 1
      },
      {
        id: 'conflicting-rendered-aliases',
        svg_id: 'stable-a',
        stable_svg_id: 'stable-a',
        rendered_svg_id: 'rendered-a',
        renderedFeatureSvgId: 'rendered-wrong',
        record_idx: 0
      }
    ],
    renderedIdentities
  });
  assert.equal(aligned.exactCount, 0);
  assert.equal(aligned.alignedCount, 0);
  assert.equal(aligned.unresolvedCount + aligned.ambiguousCount, 5);
  assert.equal(aligned.features.every((feature) => feature.rendered_svg_id === undefined), true);
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(`
    <svg>
      <path data-gbdraw-feature-id="rendered-a" data-gbdraw-stable-feature-id="stable-a" />
      <path data-gbdraw-feature-id="rendered-b" data-gbdraw-stable-feature-id="stable-b" data-gbdraw-record-index="1" />
    </svg>
  `);
  const exact = { id: 'feature-a', svg_id: 'rendered-a', stable_svg_id: 'stable-a' };
  const changed = {
    id: 'feature-b',
    svg_id: 'stable-b',
    stable_svg_id: 'stable-b',
    fileIdx: 1
  };
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [exact, changed],
    renderedIdentities
  });
  assert.notEqual(aligned.features[0], exact);
  assert.equal(aligned.features[0].rendered_svg_id, 'rendered-a');
  assert.notEqual(aligned.features[1], changed);
  assert.equal(aligned.features[1].svg_id, 'stable-b');
  assert.equal(aligned.features[1].rendered_svg_id, 'rendered-b');
}

{
  const migrated = migrateFeatureOverrideState({
    featureState: {
      extractedFeatures: [{ id: 'new-id', svg_id: 'new-svg' }],
      featureColorOverrides: { 'old-id': { color: '#111111' } },
      featureVisibilityOverrides: { 'old-svg': 'off' },
      labelTextFeatureOverrides: { 'old-svg': 'renamed' },
      labelTextFeatureOverrideSources: { 'old-id': 'manual' },
      labelVisibilityOverrides: { 'old-svg': 'hidden' }
    },
    editorState: {
      featureStrokes: {
        overrides: { 'old-svg': { strokeColor: '#222222', strokeWidth: 2 } }
      }
    },
    migration: {
      featureIdMap: { 'old-id': 'new-id' },
      svgIdMap: { 'old-svg': 'new-svg' }
    }
  });
  assert.deepEqual(migrated.featureState.featureColorOverrides, { 'new-id': { color: '#111111' } });
  assert.deepEqual(migrated.featureState.featureVisibilityOverrides, { 'new-svg': 'off' });
  assert.deepEqual(migrated.featureState.labelTextFeatureOverrides, { 'new-svg': 'renamed' });
  assert.deepEqual(migrated.featureState.labelTextFeatureOverrideSources, { 'new-id': 'manual' });
  assert.deepEqual(migrated.featureState.labelVisibilityOverrides, { 'new-svg': 'hidden' });
  assert.deepEqual(migrated.editorState.featureStrokes.overrides, {
    'new-svg': { strokeColor: '#222222', strokeWidth: 2 }
  });
}

{
  const visible = {
    id: 'visible',
    svg_id: 'rendered-a',
    stable_svg_id: 'stable-a',
    stable_feature_id: 'stable-a',
    rendered_svg_id: 'rendered-a',
    renderedSvgId: 'rendered-a',
    nucleotide_sequence: 'AAAA'
  };
  const hidden = {
    id: 'hidden',
    svg_id: 'stable-hidden',
    stable_svg_id: 'stable-hidden',
    stable_feature_id: 'stable-hidden',
    nucleotide_sequence: 'CCCC'
  };
  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'linear',
      lInputType: 'gb',
      selectedResultIndex: 0,
      results: committedResults(svgWithFeature({ renderedId: 'rendered-a', stableId: 'stable-a' })),
      featureState: {
        extractedFeatures: [visible, hidden],
        biologicalFeatures: [visible, hidden]
      },
      editorState: {}
    },
    featureVisibilityTsv: ''
  });
  assert.equal(plan.status, 'aligned');
  assert.deepEqual(plan.recoveredFeatureState.extractedFeatures.map((feature) => feature.id), ['visible']);
  assert.deepEqual(plan.recoveredFeatureState.biologicalFeatures.map((feature) => feature.id), ['visible', 'hidden']);
}

{
  const visible = {
    id: 'visible',
    svg_id: 'stable-a',
    stable_svg_id: 'stable-a',
    rendered_feature_svg_id: 'rendered-a'
  };
  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'linear',
      lInputType: 'gb',
      selectedResultIndex: 0,
      results: committedResults(svgWithFeature({ renderedId: 'rendered-a', stableId: 'stable-a' })),
      featureState: {
        extractedFeatures: [visible],
        biologicalFeatures: []
      },
      editorState: {}
    },
    featureVisibilityTsv: ''
  });
  assert.equal(plan.status, 'aligned');
  assert.equal(plan.recoveredFeatureState.extractedFeatures[0].svg_id, 'rendered-a');
}

{
  let extractionCalls = 0;
  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'linear',
      lInputType: 'gb',
      selectedResultIndex: 0,
      results: [],
      linearSeqs: [{ gb: { name: 'unused.gb' } }],
      featureState: {
        extractedFeatures: [],
        biologicalFeatures: []
      },
      editorState: {}
    },
    featureVisibilityTsv: '',
    extractFeatureMetadataForPreviewImpl: async () => {
      extractionCalls += 1;
      throw new Error('Feature extraction must not run without a preview result.');
    }
  });
  assert.equal(plan.status, 'ready');
  assert.equal(plan.reason, 'not-needed');
  assert.equal(extractionCalls, 0);
}

console.log('session feature metadata tests passed');
