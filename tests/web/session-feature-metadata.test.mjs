import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { cp, mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-session-feature-metadata-'));
// The slice editors and svg-serialization resolve a drawing's ID by its mode.
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'drawing-view.js'), join(tempDir, 'services', 'drawing-view.js'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'services', 'drawing'), { recursive: true });
await mkdir(join(tempDir, 'contracts'), { recursive: true });
await mkdir(join(tempDir, 'utils'), { recursive: true });

const copyModule = async (sourceRelative, targetRelative) => {
  await writeFile(
    join(tempDir, targetRelative),
    await readFile(join(repoRoot, sourceRelative), 'utf8'),
    'utf8'
  );
};

await copyModule('gbdraw/web/js/services/session-feature-recovery.js', 'services/session-feature-recovery.js');
await copyModule('gbdraw/web/js/services/feature-metadata-extraction.js', 'services/feature-metadata-extraction.js');
// feature-metadata-extraction.js reads the Linear Region fields through record-options.js.
for (const name of ['services/record-options.js', 'services/record-draw-selection.js', 'utils/download-names.js',
  'services/drawing/feature-edits.js', 'services/drawing/paint.js', 'services/drawing/slice-handle.js', 'utils/reactivity.js']) {
  await copyModule(`gbdraw/web/js/${name}`, name);
}
await copyModule('gbdraw/web/js/services/losat-normalization.js', 'services/losat-normalization.js');
// losat-normalization.js reads the comparison settings from the generated request description.
await copyModule('gbdraw/web/js/contracts/request-shape.generated.js', 'contracts/request-shape.generated.js');
await copyModule('gbdraw/web/js/utils/contract-shape.js', 'utils/contract-shape.js');
await copyModule('gbdraw/web/js/utils/field-reading.js', 'utils/field-reading.js');
await copyModule('gbdraw/web/js/services/feature-dom.js', 'services/feature-dom.js');
await copyModule('gbdraw/web/js/services/legend-svg.js', 'services/legend-svg.js');
await copyModule('gbdraw/web/js/services/result-paint-bases.js', 'services/result-paint-bases.js');
await copyModule('gbdraw/web/js/services/svg-transform.js', 'services/svg-transform.js');
await copyModule('gbdraw/web/js/services/diagram-generation.js', 'services/diagram-generation.js');
await copyModule('gbdraw/web/js/services/label-fit.js', 'services/label-fit.js');
await copyModule('gbdraw/web/js/services/font-store.js', 'services/font-store.js');
await copyModule('gbdraw/web/js/services/font-families.js', 'services/font-families.js');
await copyModule('gbdraw/web/js/utils/font-family-rules.generated.js', 'utils/font-family-rules.generated.js');
await copyModule('gbdraw/web/js/contracts/schema.generated.js', 'contracts/schema.generated.js');
await copyModule('gbdraw/web/js/services/bounded-json-transport.js', 'services/bounded-json-transport.js');
await copyModule(
  'gbdraw/web/js/services/current-worker-result-source.js',
  'services/current-worker-result-source.js'
);
await copyModule(
  'gbdraw/web/js/services/diagram-resource-staging.js',
  'services/diagram-resource-staging.js'
);
await copyModule(
  'gbdraw/web/js/services/canonical-resource-references.js',
  'services/canonical-resource-references.js'
);
await copyModule(
  'gbdraw/web/js/services/resource-payload-owner.js',
  'services/resource-payload-owner.js'
);
await copyModule('gbdraw/web/js/services/diagram-worker-protocol.js', 'services/diagram-worker-protocol.js');
await copyModule('gbdraw/web/js/utils/error-normalization.js', 'utils/error-normalization.js');
await copyModule('gbdraw/web/js/services/byte-utils.js', 'services/byte-utils.js');
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
await copyModule('gbdraw/web/js/utils/own-entry.js', 'utils/own-entry.js');
await copyModule('gbdraw/web/js/utils/vue-raw.js', 'utils/vue-raw.js');
await copyModule('gbdraw/web/js/services/drawn-groups.js', 'services/drawn-groups.js');
await copyModule('gbdraw/web/js/services/svg-result-normalization.js', 'services/svg-result-normalization.js');
await copyModule('gbdraw/web/js/services/svg-sanitization.js', 'services/svg-sanitization.js');
await copyModule('gbdraw/web/js/services/svg-serialization.js', 'services/svg-serialization.js');
await copyModule('gbdraw/web/js/config.js', 'config.js');
// The Legend naming rule and the key mirror of the color rules, which the Legend reader and the admission import.
for (const name of ['services/legend-row-names.js', 'services/specific-color-rules.js', 'services/file-imports.js', 'services/card-row-marks.js',
  'services/rule-matchers.js', 'utils/color-utils.js', 'utils/i18n.js', 'utils/named-colors.js', 'utils/tsv-cell.js']) {
  await copyModule(`gbdraw/web/js/${name}`, name);
}

const {
  alignRecoveredFeatureIdsToRenderedSvg,
  buildSessionFeatureRecoveryPlan,
  classifyFeatureMetadataState,
  extractSessionSourceFeatures,
  migrateFeatureOverrideState,
  normalizeRecordIndex
} = await import(pathToFileURL(join(tempDir, 'services', 'session-feature-recovery.js')));
const { collectRenderedFeatureIdentitiesFromSvgRoot } = await import(
  pathToFileURL(join(tempDir, 'services', 'session-feature-metadata.js'))
);
const { admitLegacyImportedResults, createLegacyImportResultSource } = await import(
  pathToFileURL(join(tempDir, 'services', 'svg-result-ingestion.js'))
);
const { extractFeatureMetadataForPreview } = await import(
  pathToFileURL(join(tempDir, 'services', 'feature-metadata-extraction.js'))
);
const {
  DiagramGenerationCanceledError,
  deserializeWorkerError,
  normalizeGenerationResponse
} = await import(pathToFileURL(join(tempDir, 'services', 'diagram-generation.js')));

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
  assert.deepEqual(normalizeGenerationResponse(legacyResults), {
    results: legacyResults,
    metadata: {},
    artifactIdentity: null
  });
  const metadata = { trackSlotGeometry: { schema: 1, mode: 'linear', records: [] }, annotationWarnings: [{ code: 'feature_selector_unmatched', setId: 's', annotationId: 'a', recordId: 'r', recordIndex: 0, missingCount: 1, message: 'Skipped annotation: 1 feature selector(s) unmatched.', resultIndex: 0, resultName: 'out.svg' }] };
  // The feature catalog travels beside the metadata (step 3a).
  const featureCatalog = { schema: 4, items: [] };
  assert.deepEqual(
    normalizeGenerationResponse({ results: legacyResults, metadata, featureCatalog }),
    { results: legacyResults, metadata: { ...metadata, featureCatalog }, artifactIdentity: null }
  );
  // A catalog inside the metadata is the removed shape.
  assert.throws(
    () => normalizeGenerationResponse({ results: legacyResults, metadata: { ...metadata, featureCatalog }, featureCatalog }),
    (error) => error.code === 'RESULT_INVALID' && error.stage === 'result-admission'
  );
  assert.throws(
    () => normalizeGenerationResponse({
      results: legacyResults, metadata, featureCatalog, featureCatalogRevision: { revision: 'e:1', sha256: 'x', bytes: 1 }
    }),
    (error) => error.code === 'RESULT_INVALID' && error.stage === 'result-admission'
  );
  const errorPayload = { error: { type: 'OutputError', message: 'No output files generated.' } };
  assert.deepEqual(normalizeGenerationResponse(errorPayload), {
    results: errorPayload,
    metadata: {},
    artifactIdentity: null
  });
  // Step 3b S: a Result is its document pieces; the children texts arrive in
  // one UTF-8 buffer and are decoded once each.
  const { SEGMENTS } = await import(pathToFileURL(join(tempDir, 'services', 'current-worker-result-source.js')));
  const children = new TextEncoder().encode('<text>µ</text>');
  const revision = createHash('sha256').update(children).digest('hex');
  const decoded = normalizeGenerationResponse({
    results: [{
      name: 'binary.svg',
      document: { head: '<svg>', tail: '</svg>', pieces: [
        { kind: 'raw', text: '<defs />' },
        { kind: 'content', start: '<g>', revision, bytes: children.byteLength, offset: 0, length: children.byteLength }
      ] },
      segments: children
    }],
    metadata: new TextEncoder().encode(JSON.stringify(metadata)),
    featureCatalog: new TextEncoder().encode(JSON.stringify(featureCatalog))
  });
  assert.deepEqual(decoded, {
    results: [{
      name: 'binary.svg',
      document: { head: '<svg>', tail: '</svg>', pieces: [
        { kind: 'raw', text: '<defs />' },
        { kind: 'content', start: '<g>', revision, bytes: children.byteLength, children: '<text>µ</text>' }
      ] },
      [SEGMENTS]: new Map()
    }],
    metadata: { ...metadata, featureCatalog },
    artifactIdentity: null
  });
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

const committedResults = (content) => admitLegacyImportedResults(
  createLegacyImportResultSource([{ name: 'preview.svg', content }]),
  {
    sanitizer: { sanitize: (value) => value },
    parser: IngressDomParser
  }
);

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
  // Rendered-ID Feature visibility and label edits keep their saved keys here:
  // services/feature-edit-migration.js moves them onto source identities.
  assert.deepEqual(migrated.featureState.featureVisibilityOverrides, { 'old-svg': 'off' });
  assert.deepEqual(migrated.featureState.labelTextFeatureOverrides, { 'old-svg': 'renamed' });
  assert.deepEqual(migrated.featureState.labelTextFeatureOverrideSources, { 'old-id': 'manual' });
  assert.deepEqual(migrated.featureState.labelVisibilityOverrides, { 'old-svg': 'hidden' });
  assert.equal(migrated.skippedOverrideCount, 0);
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

// OV-39: a Session without a feature catalog reads its sources again. A source
// that Python reports it cannot read falls back to the saved metadata (null);
// a runtime failure (Worker start, transport, cancellation) fails the Load
// instead of dropping the Session's feature edits as unmatched.
{
  const snapshot = { mode: 'linear', lInputType: 'gb', linearSeqs: [{ gb: { name: 'a.gb' } }] };
  const reading = (read) => extractSessionSourceFeatures({ snapshot, readFeatureExtractionDataImpl: read });
  const failing = (error) => reading(async () => { throw error; });
  const read = await reading(async () => ({ features: [{ id: 'a', svg_id: 'stable-a', record_id: 'A' }], record_ids: ['A'] }));
  assert.deepEqual(read.extractedFeatures.map((feature) => feature.svg_id), ['stable-a']);
  assert.equal(await extractSessionSourceFeatures({ snapshot: { ...snapshot, linearSeqs: [{}] } }), null);
  assert.equal(await reading(async () => ({ error: 'unreadable' })), null);
  for (const [code, stage] of [['RECORD_SELECTION', 'helper'], ['REGION_INVALID', 'helper'], ['INPUT_UNREADABLE', 'helper'],
    ['NO_RECORDS', 'helper'], ['VALIDATION_UNCLASSIFIED', 'helper'], ['RENDER_FAILED', 'render']]) {
    assert.equal(await failing(deserializeWorkerError({ code, operation: 'feature-extraction', stage })), null, code);
  }
  for (const error of [
    deserializeWorkerError({ code: 'WORKER_INIT', stage: 'initialization' }, { operation: 'feature-extraction' }),
    deserializeWorkerError({ code: 'RESULT_INVALID', operation: 'feature-extraction', stage: 'result-admission' }),
    deserializeWorkerError({ code: 'UNKNOWN', operation: 'feature-extraction', stage: 'unknown' }),
    new DiagramGenerationCanceledError()
  ]) {
    await assert.rejects(failing(error), (thrown) => thrown === error, error.code || error.name);
  }
  // OV-315: a falsy rejection is an unknown failure too (OV-70). It fails the
  // Load with the caught value, not with a TypeError from the catch block.
  for (const value of [undefined, null, 0, '']) {
    await assert.rejects(failing(value), (thrown) => thrown === value, `falsy rejection ${String(value)}`);
  }
}

// OV-751: a Session below 40 whose features are read again aligns them with
// the saved SVG. A Circular SVG names a feature without its record (the record
// index and ID sit on the record group), and a Linear feature read again names
// its input (`fileIdx`) and its record within that input (`record_idx`) and
// carries the extraction's predicted rendered ID, which is the source hash
// without the `_record_<n>` suffix the SVG draws.
{
  const circularPlan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'circular',
      cInputType: 'gb',
      files: { c_gb: { name: 'c.gb' } },
      selectedResultIndex: 0,
      results: committedResults('<svg><g data-gbdraw-record-id="R1" data-gbdraw-record-index="0"><path id="hc" data-gbdraw-feature-id="hc" /></g></svg>'),
      featureState: {},
      editorState: {}
    },
    featureVisibilityTsv: '',
    readFeatureExtractionDataImpl: async () => {
      const feature = { id: 'f0', svg_id: 'hc', stable_svg_id: 'hc', rendered_feature_svg_id: 'hc', record_idx: 0, record_id: 'R1', type: 'CDS' };
      return { features: [feature], biological_features: [feature], record_ids: ['R1'] };
    }
  });
  assert.equal(circularPlan.status, 'recovered', circularPlan.reason);
  assert.deepEqual(circularPlan.recoveredFeatureState.extractedFeatures.map((feature) => feature.rendered_svg_id), ['hc']);
}

{
  // Input 0 holds two records and input 1 one: the SVG numbers them 0, 1, 2.
  const svg = `<svg>
    <path id="a_record_1" data-gbdraw-feature-id="a_record_1" data-gbdraw-stable-feature-id="a" data-gbdraw-record-index="0" data-gbdraw-record-id="RA" />
    <path id="a2_record_2" data-gbdraw-feature-id="a2_record_2" data-gbdraw-stable-feature-id="a2" data-gbdraw-record-index="1" data-gbdraw-record-id="RA2" />
    <path id="b_record_3" data-gbdraw-feature-id="b_record_3" data-gbdraw-stable-feature-id="b" data-gbdraw-record-index="2" data-gbdraw-record-id="RB" />
  </svg>`;
  const extracted = {
    '/seq_0.gb': {
      features: [
        { id: 'f0', svg_id: 'a', stable_svg_id: 'a', rendered_feature_svg_id: 'a', record_idx: 0, record_id: 'RA', type: 'CDS' },
        { id: 'f1', svg_id: 'a2', stable_svg_id: 'a2', rendered_feature_svg_id: 'a2', record_idx: 1, record_id: 'RA2', type: 'CDS' }
      ],
      record_ids: ['RA', 'RA2']
    },
    '/seq_1.gb': {
      features: [{ id: 'f0', svg_id: 'b', stable_svg_id: 'b', rendered_feature_svg_id: 'b', record_idx: 0, record_id: 'RB', type: 'CDS' }],
      record_ids: ['RB']
    }
  };
  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'linear',
      lInputType: 'gb',
      linearSeqs: [{ gb: { name: 'a.gb' } }, { gb: { name: 'b.gb' } }],
      selectedResultIndex: 0,
      results: committedResults(svg),
      featureState: {},
      editorState: {}
    },
    featureVisibilityTsv: '',
    readFeatureExtractionDataImpl: async ({ path }) => extracted[path]
  });
  assert.equal(plan.status, 'recovered', plan.reason);
  const features = plan.recoveredFeatureState.extractedFeatures;
  assert.deepEqual(features.map((feature) => feature.rendered_svg_id), ['a_record_1', 'a2_record_2', 'b_record_3']);
  assert.deepEqual(features.map((feature) => feature.rendered_feature_svg_id), ['a_record_1', 'a2_record_2', 'b_record_3']);
  // The features keep the extraction's reading of their records, which the
  // rendered-ID edit migration reads (services/feature-edit-migration.js).
  assert.deepEqual(features.map((feature) => [feature.fileIdx, feature.record_idx]), [[0, 0], [0, 1], [1, 0]]);
}

{
  // A cropped, reverse-complemented or rotated record draws its features under
  // their hash in the drawn record (`drawn_selector.hash`), not their source
  // hash; a biological feature takes that hash from the drawn feature.
  const source = { svg_id: 'src', stable_svg_id: 'src', record_idx: 0, record_id: 'R1', type: 'CDS' };
  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'circular',
      cInputType: 'gb',
      files: { c_gb: { name: 'c.gb' } },
      selectedResultIndex: 0,
      results: committedResults('<svg><path id="drw" data-gbdraw-feature-id="drw" /></svg>'),
      featureState: {},
      editorState: {}
    },
    featureVisibilityTsv: '',
    readFeatureExtractionDataImpl: async () => ({
      features: [{ ...source, id: 'f0', rendered_feature_svg_id: 'src', drawn_selector: { hash: 'drw' } }],
      biological_features: [{ ...source, id: 'f0' }],
      record_ids: ['R1']
    })
  });
  assert.equal(plan.status, 'recovered', plan.reason);
  const [drawn] = plan.recoveredFeatureState.extractedFeatures;
  assert.equal(drawn.rendered_svg_id, 'drw');
  assert.equal(drawn.stable_svg_id, 'src');
  assert.deepEqual(plan.recoveredFeatureState.biologicalFeatures.map((feature) => feature.rendered_svg_id), ['drw']);
}

{
  // Saved metadata whose IDs all match the SVG does not lose features to an
  // alignment that cannot place them: Load keeps it as saved.
  const saved = [
    { id: 'f0', svg_id: 'a_record_1', stable_svg_id: 'a', fileIdx: 0, record_idx: 0, record_id: 'RA' },
    { id: 'f1', svg_id: 'b_record_2', stable_svg_id: 'b', fileIdx: 1, record_idx: 0, record_id: 'RB' }
  ];
  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'linear',
      lInputType: 'gb',
      selectedResultIndex: 0,
      results: committedResults(`<svg>
        <path id="a_record_1" data-gbdraw-feature-id="a_record_1" data-gbdraw-stable-feature-id="a" data-gbdraw-record-index="0" data-gbdraw-record-id="RA" />
        <path id="b_record_2" data-gbdraw-feature-id="b_record_2" data-gbdraw-stable-feature-id="b" data-gbdraw-record-index="1" data-gbdraw-record-id="RB" />
      </svg>`),
      featureState: { extractedFeatures: saved },
      editorState: {}
    },
    featureVisibilityTsv: ''
  });
  assert.notEqual(plan.status, 'aligned');
  assert.equal(plan.recoveredFeatureState, undefined);
}

// OV-751: two features that differ only by their record share one
// hash. When the drawn element carries no record (an older Linear SVG with
// bare IDs, or a Circular figure), both agree with it; neither is aligned and
// the Load lists its feature edits, as it did before OV-751.
{
  const feature = (product, recordIdx, extra = {}) => ({ id: 'f0', svg_id: 'h', stable_svg_id: 'h', rendered_feature_svg_id: 'h',
    record_idx: recordIdx, record_id: 'contig_1', type: 'CDS', product, ...extra });
  const linearPlan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'linear', lInputType: 'gb', linearSeqs: [{ gb: { name: '0.gb' } }, { gb: { name: '1.gb' } }], selectedResultIndex: 0,
      results: committedResults('<svg><path id="h" data-gbdraw-feature-id="h" /><path id="h" data-gbdraw-feature-id="h" /></svg>'),
      featureState: {}, editorState: {}
    },
    featureVisibilityTsv: '',
    readFeatureExtractionDataImpl: async ({ path }) => ({
      '/seq_0.gb': { features: [feature('A', 0)], record_ids: ['contig_1'] },
      '/seq_1.gb': { features: [feature('B', 0)], record_ids: ['contig_1'] }
    })[path]
  });
  assert.equal(linearPlan.status, 'failed', linearPlan.reason);
  assert.equal(linearPlan.recoveredFeatureState, undefined);
  const circularPlan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'circular', cInputType: 'gb', files: { c_gb: { name: 'c.gb' } }, selectedResultIndex: 0,
      results: committedResults('<svg><g data-gbdraw-record-index="0"><path id="h" data-gbdraw-feature-id="h" /></g>'
        + '<g data-gbdraw-record-index="1"><path id="h" data-gbdraw-feature-id="h" /></g></svg>'),
      featureState: {}, editorState: {}
    },
    featureVisibilityTsv: '',
    readFeatureExtractionDataImpl: async () => {
      const features = [feature('A', 0), feature('B', 1, { id: 'f1' })];
      return { features, biological_features: features, record_ids: ['contig_1', 'contig_1'] };
    }
  });
  assert.equal(circularPlan.status, 'failed', circularPlan.reason);
  assert.equal(circularPlan.recoveredFeatureState, undefined);
}

// OV-751: two features of one record that share a hash (a 0.13.0
// Linear Session such as hash-color-rule-linear-rc.v30, whose SVG draws both
// under one ID) both agree with that element. They differ by no record, so
// both stay aligned, as before OV-751.
{
  const twin = (id, product) => ({ id, svg_id: 'h', stable_svg_id: 'h', rendered_feature_svg_id: 'h',
    fileIdx: 0, record_idx: 0, record_id: 'TESTA', type: 'CDS', product });
  const twins = [twin('f1', 'TWIN_A'), twin('f2', 'dupA')];
  const drawn = '<path id="h_record_1" data-gbdraw-feature-id="h_record_1" data-gbdraw-stable-feature-id="h" '
    + 'data-gbdraw-record-index="0" data-gbdraw-record-id="TESTA" />';
  const plan = await buildSessionFeatureRecoveryPlan({
    snapshot: {
      mode: 'linear', lInputType: 'gb', linearSeqs: [{ gb: { name: 'a.gb' } }], selectedResultIndex: 0,
      results: committedResults(`<svg>${drawn}${drawn}</svg>`),
      featureState: {}, editorState: {}
    },
    featureVisibilityTsv: '',
    readFeatureExtractionDataImpl: async () => ({ features: twins, record_ids: ['TESTA'] })
  });
  assert.equal(plan.status, 'recovered', plan.reason);
  assert.deepEqual(plan.recoveredFeatureState.extractedFeatures.map((feature) => feature.rendered_svg_id), ['h_record_1', 'h_record_1']);
}

console.log('session feature metadata tests passed');
