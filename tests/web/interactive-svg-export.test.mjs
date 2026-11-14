// The Web Interactive SVG export (D5): the capture selects the shown Result's stored catalog item, the export asks
// Python for that item's sequences (`readEmbedJson`, part by part, from the encoder the CLI export uses) before it
// builds the file, and the file embeds the item with what Python returned under the export schema
// (INTERACTIVE_EMBED_SCHEMA), not the stored catalog schema. A failed read writes no file; a read that another Result outlived writes nothing and says nothing.
import assert from 'node:assert/strict';
import test from 'node:test';

import { FEATURE_CATALOG_SCHEMA, INTERACTIVE_EMBED_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';

// The few DOM calls the export makes on an SVG clone.
class FakeElement {
  constructor(name) {
    this.localName = name;
    this.attributes = new Map();
    this.childNodes = [];
    this.parentNode = null;
    this.textContent = '';
    this.style = { setProperty: () => {}, removeProperty: () => {} };
  }

  get firstChild() { return this.childNodes[0] ?? null; }

  setAttribute(name, value) { this.attributes.set(name, String(value)); }

  getAttribute(name) { return this.attributes.has(name) ? this.attributes.get(name) : null; }

  removeAttribute(name) { this.attributes.delete(name); }

  appendChild(child) { return this.insertBefore(child, null); }

  insertBefore(child, before) {
    child.parentNode = this;
    const index = before ? this.childNodes.indexOf(before) : -1;
    if (index < 0) this.childNodes.push(child); else this.childNodes.splice(index, 0, child);
    return child;
  }

  removeChild(child) { this.childNodes.splice(this.childNodes.indexOf(child), 1); }

  querySelector(selector) {
    return selector === 'defs' ? this.childNodes.find((child) => child.localName === 'defs') ?? null : null;
  }

  querySelectorAll() { return []; }
}
globalThis.document = { createElementNS: (_namespace, name) => new FakeElement(name) };
globalThis.CSS = { escape: (value) => value };
const rawTargets = new WeakMap();
globalThis.window = {
  Vue: {
    computed: (getter) => ({ get value() { return getter(); } }),
    reactive: (value) => value,
    ref: (value) => ({ value }),
    toRaw: (value) => rawTargets.get(value) ?? value
  }
};

const { captureSvgExport } = await import('../../gbdraw/web/js/services/svg-serialization.js');
const { embeddedCatalogItemJson, enrichSvgWithStandaloneInteractivity } = await import('../../gbdraw/web/js/services/standalone-interactivity.js');
const { createResultExport } = await import('../../gbdraw/web/js/app/result-actions.js');

const storedItem = {
  resultIndex: 0,
  resultName: 'diagram.svg',
  recordKeys: ['record-a'],
  features: [{ svgId: 'f0001', recordKey: 'record-a', biologicalFeatureId: 'feature-a', drawnSelector: null }],
  biologicalFeatures: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', stableFeatureId: 'feature-a', type: 'CDS',
    qualifiers: { gene: ['g'], translation: ['MK'] } }],
  orthogroups: [],
  annotations: [],
  comparisonMatches: [{ id: 'match-a', match_kind: 'pairwise', query_record_index: 0, query_record_id: 'acc-a' }]
};
// What Python's `embed_catalog_item` returns for it: the rows' sequences in the compact export encoding and the
// whole-record sources the matches read; the page gets the rows and the sources as text (`embeddedSequences`).
const embedItem = {
  ...storedItem,
  biologicalFeatures: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', stableFeatureId: 'feature-a', type: 'CDS',
    qualifiers: { gene: ['g'] }, sequenceSourceIndex: 0, amino_acid_sequence: 'MK', translationFromAminoAcidSequence: true }],
  sequenceSources: [{ key: 'linear:record:0', origin: 'linear-record', recordIndex: 0, recordId: 'acc-a', sequence: 'ATGAAA' }]
};

const exportState = (items) => ({
  selectedResultIndex: { value: 0 },
  results: { value: [{ name: 'diagram.svg' }] },
  svgContainer: { value: null },
  svgContent: { value: '<svg />' },
  downloadDpi: { value: 300 },
  downloadWhiteBackground: { value: false },
  richFeaturePopup: { value: true },
  featureCatalog: { value: { schema: FEATURE_CATALOG_SCHEMA, items } },
  editableLabels: { value: [] }
});
const drawing = {
  featureOverrides: { a: { recordKey: 'record-a', biologicalFeatureId: 'feature-a', labelText: 'renamed' } },
  labelTextBulkOverrides: {}, orthogroupNameOverrides: {}, orthogroupDescriptionOverrides: {}
};

test('the capture holds the stored catalog item of the shown Result, found by its index and name', () => {
  const { interactivity } = captureSvgExport(exportState([storedItem]), drawing, { interactive: true });
  assert.equal(interactivity.catalogItem, storedItem);
  assert.deepEqual(interactivity.labelTextFeatureOverrides, { f0001: 'renamed' });
  for (const other of [{ ...storedItem, resultName: 'other.svg' }, { ...storedItem, resultIndex: 1 }]) {
    const captured = captureSvgExport(exportState([other]), drawing, { interactive: true }).interactivity;
    assert.equal(captured.catalogItem, null);
    assert.deepEqual(captured.labelTextFeatureOverrides, {});
  }
  assert.equal(captureSvgExport(exportState([storedItem]), drawing).interactivity, undefined);
  // The page holds the catalog behind a Vue proxy; the Worker gets the raw item (a proxy cannot be posted).
  const state = exportState([]);
  const raw = { schema: FEATURE_CATALOG_SCHEMA, items: [storedItem] };
  rawTargets.set(state.featureCatalog.value, raw);
  assert.equal(captureSvgExport(state, drawing, { interactive: true }).interactivity.catalogItem, storedItem);
});

const embeddedSequences = {
  biologicalFeatures: JSON.stringify(embedItem.biologicalFeatures), sequenceSources: JSON.stringify(embedItem.sequenceSources)
};

test('the file embeds the item Python encoded, under the export schema, byte for byte', () => {
  assert.equal(INTERACTIVE_EMBED_SCHEMA, 5);
  assert.notEqual(INTERACTIVE_EMBED_SCHEMA, FEATURE_CATALOG_SCHEMA);
  const svg = new FakeElement('svg');
  assert.equal(enrichSvgWithStandaloneInteractivity(svg, { catalogItem: storedItem, embeddedSequences }), true);
  const metadata = svg.childNodes.find((child) => child.localName === 'metadata');
  assert.equal(metadata.getAttribute('data-schema'), String(INTERACTIVE_EMBED_SCHEMA));
  assert.equal(metadata.getAttribute('data-result-index'), '0');
  assert.equal(metadata.getAttribute('data-result-name'), 'diagram.svg');
  assert.equal(metadata.textContent, JSON.stringify({ schema: INTERACTIVE_EMBED_SCHEMA, items: [embedItem] }));
  // Without Python's sequences the export writes nothing: the stored item holds no sequence.
  for (const options of [{ catalogItem: storedItem }, { embeddedSequences }]) {
    assert.equal(enrichSvgWithStandaloneInteractivity(new FakeElement('svg'), options), false);
  }
});

test('an item whose matches read no source embeds no sequenceSources key, as Python writes it', () => {
  const rows = embedItem.biologicalFeatures.map(({ sequenceSourceIndex, ...row }) => ({ ...row, nucleotide_sequence: 'ATG' }));
  const sequences = { biologicalFeatures: JSON.stringify(rows), sequenceSources: null };
  assert.equal(embeddedCatalogItemJson(storedItem, sequences), JSON.stringify({ ...storedItem, biologicalFeatures: rows }));
});

/**
 * The export action over stubbed owners: `readEmbedJson` answers as given, and the export service records the
 * snapshot it would write.
 * @param {(item: object) => Promise<object | null>} readEmbedJson
 */
const exportAction = (readEmbedJson) => {
  const written = [];
  const asked = [];
  const errorLog = { value: null };
  const actions = createResultExport({
    state: { loadedFonts: { value: [] } },
    errorLog,
    showNotice: () => {},
    captureSvgExport: (state, shown, options) => captureSvgExport(exportState([storedItem]), drawing, options),
    loadExportService: async () => ({
      downloadInteractiveSVG: async (snapshot) => {
        written.push(snapshot);
        return { filename: 'diagram.interactive.svg', fontWarnings: [] };
      }
    }),
    runDiagramHelperOperation: async () => ({ result: null }),
    DIAGRAM_HELPER_OPERATIONS: { EMBED_FONTS: 'embedFonts', PDF_FONTS: 'pdfFonts', CHECK_EXPORT_LIMIT: 'checkExportLimit' },
    displayedArtifactDrawing: () => drawing,
    embedInExports: { value: false },
    shownInputSummary: { value: null },
    shownLinearWrapLeaf: () => null,
    shownCommittedRequest: () => null,
    awaitLinearRecordDiscovery: async () => {},
    committedLinearRecordsOf: () => [],
    readEmbedJson: async (item) => {
      asked.push(item);
      return readEmbedJson(item);
    }
  });
  return { actions, written, asked, errorLog };
};

test('Interactive SVG export asks Python for the shown item and writes what it returns', async () => {
  const { actions, written, asked, errorLog } = exportAction(async () => embeddedSequences);
  const outcome = await actions.downloadInteractiveSVG();
  assert.equal(outcome.filename, 'diagram.interactive.svg');
  assert.deepEqual(asked, [storedItem]);
  assert.equal(written.length, 1);
  assert.equal(written[0].interactivity.embeddedSequences, embeddedSequences);
  assert.equal(written[0].interactivity.catalogItem, storedItem);
  assert.equal(errorLog.value, null);
});

test('a read that another Result outlived writes nothing and says nothing', async () => {
  const { actions, written, errorLog } = exportAction(async () => null);
  assert.deepEqual(await actions.downloadInteractiveSVG(), { status: 'stale' });
  assert.equal(written.length, 0);
  assert.equal(errorLog.value, null);
  assert.equal(actions.canRetryInteractiveSvgExport.value, false);
});

test('a failed read writes no file and shows why, with Retry', async () => {
  const failure = Object.assign(new Error('RESOURCE_INVALID/'), {
    code: 'RESOURCE_INVALID', stage: 'helper', operation: 'readEmbedSequences', context: {}
  });
  const { actions, written, errorLog } = exportAction(async () => { throw failure; });
  const outcome = await actions.downloadInteractiveSVG();
  assert.equal(outcome.status, 'error');
  assert.equal(written.length, 0);
  assert.equal(errorLog.value, outcome.error);
  assert.equal(errorLog.value.code, 'RESOURCE_INVALID');
  assert.equal(errorLog.value.operation, 'export-svg');
  assert.equal(actions.canRetryInteractiveSvgExport.value, true);
});
