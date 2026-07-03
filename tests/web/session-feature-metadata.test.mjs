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
await copyModule('gbdraw/web/js/services/diagram-generation.js', 'services/diagram-generation.js');
await copyModule('gbdraw/web/js/services/pyodide-assets.js', 'services/pyodide-assets.js');
await copyModule('gbdraw/web/js/config.js', 'config.js');

const {
  alignRecoveredFeatureIdsToRenderedSvg,
  classifyFeatureMetadataState,
  collectRenderedFeatureIdentitiesFromSvg,
  migrateFeatureOverrideState
} = await import(pathToFileURL(join(tempDir, 'app', 'session-feature-metadata.js')));
const { normalizeGenerationResponse } = await import(pathToFileURL(join(tempDir, 'services', 'diagram-generation.js')));

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
}

class FakeDocument {
  constructor(elements) {
    this.elements = elements;
  }

  querySelector(selector) {
    return selector === 'parsererror' ? null : null;
  }

  querySelectorAll() {
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
  const state = classifyFeatureMetadataState({
    results: [{ content: svgWithFeature({ renderedId: 'rendered-a' }) }],
    selectedResultIndex: 0,
    extractedFeatures: [{ id: 'feature-a', svg_id: 'rendered-a' }]
  });
  assert.equal(state.state, 'ready');
  assert.equal(state.exactMatchingCount, 1);
  assert.equal(state.missingExactCount, 0);
}

{
  const state = classifyFeatureMetadataState({
    results: [{ content: svgWithFeature({ renderedId: 'rendered-a', stableId: 'stable-a' }) }],
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
  assert.equal(aligned.changedCount, 1);
  assert.equal(aligned.alignedCount, 1);
  assert.equal(aligned.features[0].svg_id, 'rendered-a');
  assert.equal(aligned.features[0].stable_svg_id, 'stable-a');
  assert.notEqual(aligned.features[0], feature);
  assert.equal(aligned.svgIdMap['stable-a'], 'rendered-a');
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
    renderedIdentities
  });
  assert.equal(aligned.ambiguousCount, 0);
  assert.equal(aligned.unresolvedCount, 0);
  assert.equal(aligned.features[0].svg_id, 'stable-x_record_1');
  assert.equal(aligned.features[1].svg_id, 'stable-x_record_2');
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
  assert.equal(aligned.features[0], feature);
}

{
  const renderedIdentities = collectRenderedFeatureIdentitiesFromSvg(`
    <svg>
      <path data-gbdraw-feature-id="rendered-a" data-gbdraw-stable-feature-id="stable-a" />
      <path data-gbdraw-feature-id="rendered-b" data-gbdraw-stable-feature-id="stable-b" />
    </svg>
  `);
  const exact = { id: 'feature-a', svg_id: 'rendered-a', stable_svg_id: 'stable-a' };
  const changed = { id: 'feature-b', svg_id: 'stable-b', stable_svg_id: 'stable-b' };
  const aligned = alignRecoveredFeatureIdsToRenderedSvg({
    features: [exact, changed],
    renderedIdentities
  });
  assert.equal(aligned.features[0], exact);
  assert.notEqual(aligned.features[1], changed);
  assert.equal(aligned.features[1].svg_id, 'rendered-b');
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

console.log('session feature metadata tests passed');
