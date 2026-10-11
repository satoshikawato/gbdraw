// The legend-entry and record-definition popups of the preview: what a click
// reaches, what each popup shows, and the setting a record edit writes.
import assert from 'node:assert/strict';
import test from 'node:test';

import {
  buildLegendEntryView,
  buildRecordDefinitionView,
  createPreviewItemPopup,
  previewItemTarget,
  readDrawnDefinitionLines
} from '../../gbdraw/web/js/app/preview-item-popup.js';

class FakeElement {
  constructor(attributes = {}, { parent = null, text = '', children = [] } = {}) {
    this.attributes = attributes;
    this.parent = parent;
    this.textContent = text;
    this.children = children;
  }

  getAttribute(name) { return this.attributes[name] ?? null; }

  hasAttribute(name) { return name in this.attributes; }

  matches(selector) {
    return selector.split(',').some((part) => {
      const query = part.trim();
      if (query === 'g[data-legend-key]') return this.hasAttribute('data-legend-key');
      const role = query.match(/^\[data-gbdraw-role="([^"]+)"\]$/)?.[1];
      if (role) return this.getAttribute('data-gbdraw-role') === role;
      const attribute = query.match(/^\[([\w-]+)\]$/)?.[1];
      return Boolean(attribute && this.hasAttribute(attribute));
    });
  }

  closest(selector) {
    for (let node = this; node; node = node.parent) if (node.matches(selector)) return node;
    return null;
  }

  querySelectorAll(selector) {
    const found = [];
    const visit = (node) => node.children.forEach((child) => {
      if ((selector === 'text' && child.attributes.tag === 'text') || (selector !== 'text' && child.matches(selector))) found.push(child);
      visit(child);
    });
    visit(this);
    return found;
  }
}

const withElement = (run) => {
  const previous = globalThis.Element;
  globalThis.Element = FakeElement;
  try { return run(); } finally { globalThis.Element = previous; }
};

test('a click on a Legend row or a record definition names its popup item', () => withElement(() => {
  const svg = { contains: () => true };
  const legend = new FakeElement({ 'data-legend-key': 'tRNA' });
  const swatch = new FakeElement({}, { parent: legend });
  assert.deepEqual(previewItemTarget(swatch, svg)?.item, { kind: 'legend-entry', key: 'tRNA' });

  const definition = new FakeElement({
    'data-gbdraw-role': 'record-definition-row', 'data-gbdraw-record-index': '1', 'data-gbdraw-record-id': 'BGC0000709'
  });
  const line = new FakeElement({}, { parent: definition });
  assert.deepEqual(previewItemTarget(line, svg)?.item, { kind: 'record-definition', recordIndex: 1, recordId: 'BGC0000709' });

  // The comparison legend (similarity bar) opens no preview popup.
  const comparison = new FakeElement({ 'data-gbdraw-role': 'comparison-legend' });
  const comparisonRow = new FakeElement({ 'data-legend-key': 'Similarity' }, { parent: comparison });
  assert.equal(previewItemTarget(comparisonRow, svg), null);
  assert.equal(previewItemTarget(new FakeElement({}), svg), null);
  assert.equal(previewItemTarget(line, { contains: () => false }), null);
}));

test('the legend-entry view names the rule behind the color and the feature count', () => {
  const entry = { caption: 'Core biosynthetic genes', color: '#d62728' };
  const ruled = buildLegendEntryView({
    entry, index: 2, color: '#d62728', featureCount: 16, featureTypes: new Set(['CDS']),
    rules: [{ feat: 'CDS', qual: 'gene_kind', val: 'biosynthetic$' }]
  });
  assert.equal(ruled.subtitle, '16 features');
  assert.equal(ruled.colorSource, 'From this color rule:');
  assert.deepEqual(ruled.ruleLines, ['CDS where gene_kind matches “biosynthetic$”']);
  assert.equal(ruled.dialogLabel, 'Legend entry: Core biosynthetic genes');

  const byType = buildLegendEntryView({
    entry: { caption: 'tRNA' }, index: 0, color: '#ff0000', featureCount: null, featureTypes: new Set(['tRNA']), rules: []
  });
  assert.equal(byType.subtitle, 'Counting features…');
  assert.equal(byType.colorSource, 'Default color of tRNA features');
  const single = buildLegendEntryView({
    entry: { caption: 'x' }, index: 0, color: '#000000', featureCount: 1, featureTypes: new Set(), rules: []
  });
  assert.equal(single.subtitle, '1 feature');
  assert.equal(single.colorSource, 'Set on individual features');
});

test('the record-definition view lists the facts of the displayed Result', () => {
  const view = buildRecordDefinitionView({
    scope: 'linear', recordIndex: 0, recordId: 'BGC0000708',
    drawn: { name: 'Streptomyces sp.', subtitle: 'neocarzinostatin', accession: 'BGC0000708', length: '28,450 bp' },
    record: { recordKey: 'uid-1', presentation: { reverseComplement: true }, display: { startCoordinate: 1200 }, region: { start: 10, end: 2000 } },
    draft: null,
    alignment: { recordKey: 'uid-1', status: 'aligned', anchorLabel: 'neoU' },
    alignmentShown: true
  });
  assert.deepEqual(view.facts, [
    { label: 'Accession', value: 'BGC0000708' },
    { label: 'Organism', value: 'Streptomyces sp.' },
    { label: 'Length', value: '28,450 bp' },
    { label: 'Reversed', value: 'Yes' },
    { label: 'Display start', value: '1,200' },
    { label: 'Region', value: '10..2,000' },
    { label: 'Alignment', value: 'Aligned on neoU' }
  ]);
  assert.equal(view.subtitle, 'neocarzinostatin · BGC0000708');
  const circular = buildRecordDefinitionView({
    scope: 'circular', recordIndex: 0, recordId: 'NC_012920.1', drawn: {}, record: null, draft: null, alignment: null, alignmentShown: false
  });
  assert.equal(circular.facts.find((fact) => fact.label === 'Region').value, 'Whole record');
  assert.equal(circular.facts.some((fact) => fact.label === 'Alignment'), false);
});

test('a Circular definition without line kinds reads its name, accession, and length lines', () => {
  const texts = ['Homo sapiens', 'mitochondrion', 'NC_012920.1', '16,569 bp', '44.36% GC']
    .map((line) => new FakeElement({ tag: 'text' }, { text: line }));
  const group = new FakeElement({}, { children: texts });
  assert.deepEqual(readDrawnDefinitionLines(group, 'NC_012920.1'), {
    length: '16,569 bp', name: 'Homo sapiens', subtitle: 'mitochondrion', accession: 'NC_012920.1'
  });
  const linear = new FakeElement({}, {
    children: [
      new FakeElement({ 'data-definition-line-kind': 'name' }, { text: 'Streptomyces' }),
      new FakeElement({ 'data-definition-line-kind': 'length' }, { text: '28 kbp' })
    ]
  });
  assert.deepEqual(readDrawnDefinitionLines(linear, 'x'), { name: 'Streptomyces', length: '28 kbp' });
});

const createPopup = ({ mode = 'linear', busy = null } = {}) => {
  const steps = [];
  const form = { species: 'Homo sapiens', strain: '', circular_reverse: false };
  const sequence = { uid: 'uid-1', definition: 'Old name', record_subtitle: '', region_reverse: false };
  const entries = [{ caption: 'tRNA', originalCaption: 'tRNA', color: '#ff0000' }];
  const state = {
    clickedPairwiseMatch: { value: null },
    svgContainer: { value: null },
    mode: { value: mode },
    linearSeqs: [sequence],
    extractedFeatures: { value: [{ type: 'tRNA' }] },
    activeDrawing: () => ({ form, legendEntries: { value: entries } })
  };
  const popup = createPreviewItemPopup({
    state,
    computed: (getter) => ({ get value() { return getter(); } }),
    readLegendEntryColor: (entry) => entry.color,
    readLegendRowRules: () => [],
    countLegendCaptionFeatures: () => 22,
    getCommittedRequest: () => ({ records: [{ recordKey: 'uid-1', presentation: { label: 'Old name' } }] }),
    readAlignmentInspector: () => null,
    readCircularReverseAvailable: () => true,
    runUndoable: (label, apply) => { steps.push(label); apply(); return true; },
    sessionBusy: () => busy,
    renameLegendRow: (index, caption) => {
      steps.push('Rename legend item');
      Object.assign(entries[index], { caption, originalCaption: caption });
      return true;
    }
  });
  return { popup, state, steps, form, sequence };
};

test('the popup view follows the open item and counts once the matches are ready', () => {
  const { popup, state } = createPopup();
  assert.equal(popup.view.value, null);
  state.clickedPairwiseMatch.value = { previewItem: { kind: 'legend-entry', key: 'tRNA' }, matchesReady: false };
  assert.equal(popup.view.value.subtitle, 'Counting features…');
  state.clickedPairwiseMatch.value.matchesReady = true;
  assert.equal(popup.view.value.subtitle, '22 features');
  state.clickedPairwiseMatch.value = { previewItem: { kind: 'legend-entry', key: 'gone' } };
  assert.equal(popup.view.value, null);
});

test('a record edit writes the sidebar setting as one History step', () => {
  const linear = createPopup();
  linear.state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 0, recordId: 'R1' } };
  const view = linear.popup.view.value;
  assert.equal(view.draft.name, 'Old name');
  linear.popup.setRecordField(view, 'name', 'New name');
  linear.popup.setRecordField(view, 'reverse', true);
  linear.popup.setRecordField(view, 'reverse', true);
  assert.equal(linear.sequence.definition, 'New name');
  assert.equal(linear.sequence.region_reverse, true);
  assert.deepEqual(linear.steps, ['Edit record definition', 'Change record orientation']);

  const circular = createPopup({ mode: 'circular' });
  circular.state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 0, recordId: 'R1' } };
  circular.popup.setRecordField(circular.popup.view.value, 'subtitle', 'mitochondrion');
  assert.equal(circular.form.strain, 'mitochondrion');

  const busy = createPopup({ busy: 'Session loading' });
  busy.state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 0, recordId: 'R1' } };
  assert.equal(busy.popup.setRecordField(busy.popup.view.value, 'name', 'X'), null);
  assert.deepEqual(busy.steps, []);
});

test('the legend popup follows its row when a rename changes the row key', async () => {
  const { popup, state, steps } = createPopup();
  state.clickedPairwiseMatch.value = { previewItem: { kind: 'legend-entry', key: 'tRNA' }, matchesReady: true };
  await popup.renameLegendRow(0, 'Transfer RNA');
  assert.deepEqual(steps, ['Rename legend item']);
  assert.deepEqual(state.clickedPairwiseMatch.value.previewItem, { kind: 'legend-entry', key: 'Transfer RNA' });
  assert.equal(popup.view.value.title, 'Transfer RNA');
});
