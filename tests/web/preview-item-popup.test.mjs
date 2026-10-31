// The legend-entry and record-definition popups of the preview: what a click
// reaches, what each popup shows, and the setting a record edit writes.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { createPreviewItemPopup, defaultColorScopeLine, previewItemTarget } from '../../gbdraw/web/js/app/preview-item-popup.js';
import { legendEntryKey } from '../../gbdraw/web/js/services/drawing/paint.js';
import { trackLegendRows } from '../../gbdraw/web/js/services/track-legend-rows.js';
import { SOURCE_LOCALE, setUiLocale } from '../../gbdraw/web/js/utils/i18n.js';

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

const createPopup = ({
  mode = 'linear', busy = null, circularBlock = 'not-inspected', recordLabel = '',
  entries = [{ caption: 'tRNA', originalCaption: 'tRNA', color: '#ff0000' }],
  renameLegendRow = null, legendRenameSettled = async () => {}, readShownLegendName = undefined,
  requestRecords = [{ recordKey: 'uid-1', presentation: { label: 'Old name' } }],
  recordRows = [{ scope: 'linear', sourceUid: 'uid-1', recordId: 'R1', canonicalRecordKey: 'uid-1' }],
  rules = [], featureTypes = ['tRNA'], featureCount = 22, definitionGroup = null, alignment = null, trackKinds = {},
  countLegendCaptionFeatures = () => featureCount, defaultColorTypes = {}, userDefaultColors = {},
  readLegendTrackKind = (entry) => trackKinds[entry.originalCaption || entry.caption] ?? null,
  palettes = {}, appliedPalette = 'default'
} = {}) => {
  const steps = [];
  const form = {
    species: 'Homo sapiens', strain: '', circular_reverse: false, circular_record_label: recordLabel, circular_record_subtitle: ''
  };
  const sequence = { uid: 'uid-1', definition: 'Old name', record_subtitle: '', region_reverse: false };
  const state = {
    clickedPairwiseMatch: { value: null },
    svgContainer: { value: definitionGroup && { querySelector: () => definitionGroup } },
    mode: { value: mode },
    linearSeqs: [sequence],
    extractedFeatures: { value: featureTypes.map((type) => ({ type })) },
    paletteDefinitions: { value: palettes },
    appliedPaletteName: { value: appliedPalette },
    activeDrawing: () => ({ form, legendEntries: { value: entries } })
  };
  const popup = createPreviewItemPopup({
    state,
    computed: (getter) => ({ get value() { return getter(); } }),
    readLegendEntryColor: (entry) => entry.color,
    readShownLegendName,
    readLegendRowRules: () => rules,
    countLegendCaptionFeatures,
    readLegendTrackKind,
    readLegendRowDefaultColorType: (entry) => defaultColorTypes[entry.originalCaption || entry.caption] ?? null,
    readUserDefaultColor: (type) => userDefaultColors[type] ?? null,
    getCommittedRequest: () => ({ records: requestRecords }),
    readCommittedRecordRows: () => recordRows,
    readAlignmentInspector: () => alignment && { records: [alignment] },
    readCircularOwnRecordBlock: () => circularBlock,
    runUndoable: (label, apply) => { steps.push(label); apply(); return true; },
    sessionBusy: () => busy,
    renameLegendRow: renameLegendRow || ((index, caption) => {
      steps.push('Rename legend item');
      Object.assign(entries[index], { caption, originalCaption: caption });
      return true;
    }),
    legendRenameSettled
  });
  return { popup, state, steps, form, sequence, entries };
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
  assert.deepEqual(linear.steps, ['Edit record definition', 'Change Reverse complement']);

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

// A rename onto another row's name asks Merge, Keep both, or Cancel; the
// rename returns once the dialog is open. The popup keeps its row until the
// choice is answered, then follows it.
test('the legend popup follows its row through the name-conflict choice', async () => {
  const conflict = (choice) => {
    let answer;
    // The target row comes first, so a lookup by the typed name finds the target.
    const entries = [
      { caption: 'Core biosynthetic genes', originalCaption: 'Core biosynthetic genes' },
      { caption: 'Regulatory genes', originalCaption: 'Regulatory genes' }
    ];
    const settled = new Promise((resolve) => { answer = resolve; });
    const created = createPopup({ entries, renameLegendRow: async () => true, legendRenameSettled: () => settled });
    created.state.clickedPairwiseMatch.value = { previewItem: { kind: 'legend-entry', key: 'Regulatory genes' }, matchesReady: true };
    const renamed = created.popup.renameLegendRow(1, 'Core biosynthetic genes');
    const choose = async () => {
      assert.equal(created.popup.view.value.title, 'Regulatory genes', 'the popup waits for the choice');
      // PD-OI-061 scenario revision 5: Keep both keeps the row's key and writes the typed name.
      if (choice === 'keep both') Object.assign(entries[1], { caption: 'Core biosynthetic genes' });
      if (choice === 'merge') entries.splice(1, 1);
      answer();
      await renamed;
      return created.state.clickedPairwiseMatch.value.previewItem.key;
    };
    return choose;
  };
  assert.equal(await (conflict('cancel'))(), 'Regulatory genes');
  assert.equal(await (conflict('keep both'))(), 'Regulatory genes', 'Keep both stays on the renamed row, by its key');
  assert.equal(await (conflict('merge'))(), 'Core biosynthetic genes');
});

// A Circular record with its own Record label draws that label, not Species,
// so the popup edits the field the Result draws.
test('a Circular record popup edits the Record label when the record has its own', () => {
  const own = createPopup({ mode: 'circular', circularBlock: '', recordLabel: 'Label X' });
  own.state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 0, recordId: 'R1' } };
  const view = own.popup.view.value;
  assert.equal(view.draft.nameLabel, 'Record label');
  assert.equal(view.draft.name, 'Label X');
  own.popup.setRecordField(view, 'name', 'Label Y');
  own.popup.setRecordField(view, 'subtitle', 'plasmid');
  assert.equal(own.form.circular_record_label, 'Label Y');
  assert.equal(own.form.circular_record_subtitle, 'plasmid');
  assert.equal(own.form.species, 'Homo sapiens');

  // OV-432: a Record label drawn while its fields are unavailable overrides
  // Species, so the popup offers no name field.
  const shared = createPopup({ mode: 'circular', recordLabel: 'Label X' });
  shared.state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 0, recordId: 'R1' } };
  assert.equal(shared.popup.view.value.draft, null);
  shared.form.circular_record_label = '';
  assert.equal(shared.popup.view.value.draft.nameLabel, 'Species');
  assert.equal(shared.popup.view.value.draft.reverseAvailable, false);
});

// The popup says why the Record label cannot be edited, in the sentence the
// sidebar shows for the same cause where it has one.
test('a Circular record popup without its own fields names the cause', () => {
  const sidebar = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const hint = (circularBlock) => {
    const { popup, state } = createPopup({ mode: 'circular', circularBlock, recordLabel: 'Label X' });
    state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 0, recordId: 'R1' } };
    assert.equal(popup.view.value.draft, null);
    return popup.view.value.editHint;
  };
  assert.equal(hint('not-inspected'), 'Inspect source records in the sidebar to edit this record\'s label.');
  assert.equal(hint('multi-record-canvas'), 'Multi-Record Canvas uses a grid. Turn it off to select one record.');
  assert.equal(hint('no-single-record'), 'Select one record to edit crop, orientation, and title lines.');
  for (const cause of ['multi-record-canvas', 'no-single-record']) {
    assert.ok(sidebar.includes(`>{{ t('${hint(cause)}') }}</p>`), `the sidebar shows the ${cause} sentence`);
  }
});

// A Linear file that draws several records: the drawn record id names the
// request record and its File card, not the drawn index.
test('a Linear record popup finds its request record and File card by record id', () => {
  const requestRecords = [
    { recordKey: 'uid-1', cardinality: 'all', presentation: { label: 'Old name' } },
    { recordKey: 'uid-2', cardinality: 'exactly_one', presentation: { label: 'Other', reverseComplement: true } }
  ];
  const recordRows = [
    { scope: 'linear', sourceUid: 'uid-1', recordId: 'R1', canonicalRecordKey: 'uid-1' },
    { scope: 'linear', sourceUid: 'uid-1', recordId: 'R2', canonicalRecordKey: 'uid-1' },
    { scope: 'linear', sourceUid: 'uid-2', recordId: 'R3', canonicalRecordKey: 'uid-2' }
  ];
  const { popup, state } = createPopup({ requestRecords, recordRows });
  // The second drawn record (index 1) is R2 of the first file, not request record 1.
  state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 1, recordId: 'R2' } };
  const view = popup.view.value;
  assert.equal(view.draft.sequenceUid, 'uid-1');
  assert.equal(view.facts.find((fact) => fact.label === 'Reverse complement').value, 'No');
  state.clickedPairwiseMatch.value = { previewItem: { kind: 'record-definition', recordIndex: 1, recordId: 'unknown' } };
  assert.equal(popup.view.value.draft, null);
});

// A loaded Session has read no record list yet: the drawn index names the
// record while every request record draws exactly one.
test('a Linear record popup without record rows uses the drawn index of one-record files only', () => {
  const item = { previewItem: { kind: 'record-definition', recordIndex: 0, recordId: 'R1' } };
  const one = createPopup({
    requestRecords: [{ recordKey: 'uid-1', cardinality: 'exactly_one', presentation: { label: 'Old name' } }], recordRows: []
  });
  one.state.clickedPairwiseMatch.value = item;
  assert.equal(one.popup.view.value.draft.sequenceUid, 'uid-1');
  const whole = createPopup({ requestRecords: [{ recordKey: 'uid-1', cardinality: 'all' }], recordRows: [] });
  whole.state.clickedPairwiseMatch.value = item;
  assert.equal(whole.popup.view.value.draft, null);
});

/** Opens the popup of `item` and returns its view. */
const openView = (created, previewItem, matchesReady = true) => {
  created.state.clickedPairwiseMatch.value = { previewItem, matchesReady };
  return created.popup.view.value;
};

test('the legend-entry view names the rule behind the color and the feature count', () => {
  const ruled = openView(createPopup({
    entries: [{ caption: 'Core biosynthetic genes', color: '#d62728' }], featureTypes: ['CDS'], featureCount: 16,
    rules: [{ feat: 'CDS', qual: 'gene_kind', val: 'biosynthetic$' }]
  }), { kind: 'legend-entry', key: 'Core biosynthetic genes' });
  assert.equal(ruled.subtitle, '16 features');
  assert.equal(ruled.colorSource, 'From this color rule:');
  assert.deepEqual(ruled.ruleLines, ['CDS where gene_kind matches “biosynthetic$”']);
  assert.equal(ruled.dialogLabel, 'Legend item: Core biosynthetic genes');

  assert.equal(ruled.name, 'Core biosynthetic genes');

  // The popup shows a row by the name the Legend draws, never by its key.
  const numbered = openView(createPopup({
    entries: [{ caption: 'Shared [#112233]', originalCaption: 'Shared [#112233]', color: '#112233' }],
    readShownLegendName: (caption) => (caption === 'Shared [#112233]' ? 'Shared (1)' : caption)
  }), { kind: 'legend-entry', key: 'Shared [#112233]' });
  assert.deepEqual([numbered.name, numbered.title, numbered.dialogLabel], ['Shared (1)', 'Shared (1)', 'Legend item: Shared (1)']);
  assert.equal(JSON.stringify(numbered).includes('[#112233]') && !numbered.identity.includes('[#112233]'), false, 'the key is its identity only');

  const byType = openView(createPopup(), { kind: 'legend-entry', key: 'tRNA' }, false);
  assert.equal(byType.subtitle, 'Counting features…');
  assert.equal(byType.colorSource, 'Default color of tRNA features');
  // OV-1281 / OIC-029: once the matches are ready, the view carries the one line of what Color sets; the owner shapes it.
  assert.deepEqual(defaultColorScopeLine('tRNA', null), { kind: 'sets', type: 'tRNA', color: null });
  assert.deepEqual(defaultColorScopeLine('tRNA', '#112233'), { kind: 'replaces', type: 'tRNA', color: '#112233' });
  assert.equal(defaultColorScopeLine(null, '#112233'), null);
  const drawn = openView(createPopup({ defaultColorTypes: { tRNA: 'tRNA' } }), { kind: 'legend-entry', key: 'tRNA' });
  // OV-1336: once the line is shown it is the only text; the caption does not repeat it.
  assert.deepEqual([drawn.colorSource, drawn.defaultColorLine], [null, { kind: 'sets', type: 'tRNA', color: null }]);
  // N-2: a type that already has a user default color says it replaces it.
  assert.deepEqual(openView(createPopup({ defaultColorTypes: { tRNA: 'tRNA' }, userDefaultColors: { tRNA: '#112233' } }),
    { kind: 'legend-entry', key: 'tRNA' }).defaultColorLine, { kind: 'replaces', type: 'tRNA', color: '#112233' });
  const early = openView(createPopup({ defaultColorTypes: { tRNA: 'tRNA' } }), { kind: 'legend-entry', key: 'tRNA' }, false);
  assert.deepEqual([early.colorSource, early.defaultColorLine], ['Default color of tRNA features', null],
    'before the matches are ready the popup does not say what Color sets');
  const single = openView(createPopup({ entries: [{ caption: 'x', color: '#000000' }], featureTypes: [], featureCount: 1 }),
    { kind: 'legend-entry', key: 'x' });
  assert.equal(single.subtitle, '1 feature');
  assert.equal(single.colorSource, 'Color set on individual features');
});

// OV-1466: a rule for any feature is one sentence, so a locale can write the
// words for "any feature" where and how its sentence needs them.
test('a color rule for any feature is one whole sentence in every locale (OV-1466)', async () => {
  const lines = () => openView(createPopup({
    entries: [{ caption: 'Ruled', color: '#d62728' }],
    rules: [{ feat: '', qual: '', val: '' }, { feat: '', qual: 'product', val: 'kinase' }]
  }), { kind: 'legend-entry', key: 'Ruled' }).ruleLines;
  assert.deepEqual(lines(), ['Any feature color rule', 'Any feature where product matches “kinase”']);
  try {
    await setUiLocale('ja', async () => ({ default: Object.freeze({
      'Any feature color rule': 'A', 'Any feature where {qualifier} matches “{pattern}”': 'B {qualifier} {pattern}',
      '{feature} color rule': '{feature} C', '{feature} where {qualifier} matches “{pattern}”': '{feature} D'
    }) }));
    assert.deepEqual(lines(), ['A', 'B product kinase']);
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

test('a legend-entry view of color none shows None; its picker starts at the palette color (OV-1386)', () => {
  const none = (matchesReady) => openView(createPopup({
    entries: [{ caption: 'CDS', originalCaption: 'CDS', color: 'none' }], featureTypes: ['CDS'],
    defaultColorTypes: { CDS: 'CDS' }, userDefaultColors: { CDS: 'none' }, palettes: { default: { CDS: '#54BCF8' } }
  }), { kind: 'legend-entry', key: 'CDS' }, matchesReady);
  assert.deepEqual(none(true).swatch, { none: true, fill: 'transparent', picker: '#54bcf8', label: 'None' });
  assert.deepEqual(none(false).swatch, { none: true, fill: 'transparent', picker: '#cccccc', label: 'None' },
    'before the matches are ready the type is unknown: the Palette pane\'s last fallback');
  assert.deepEqual(openView(createPopup(), { kind: 'legend-entry', key: 'tRNA' }).swatch,
    { none: false, fill: '#ff0000', picker: '#ff0000', label: '#ff0000' });
});

test('the legend-entry view of a track row names its plot, not features (OV-1272)', () => {
  let counted = 0;
  const created = createPopup({
    entries: [
      { caption: 'GC content', originalCaption: 'GC content', color: '#a0a0a0' },
      { caption: 'GC skew (+)', originalCaption: 'GC skew (+)', color: '#6dded3' },
      { caption: 'Depth', originalCaption: 'Depth', color: '#333333' }
    ],
    featureTypes: ['CDS'], featureCount: 0,
    // A rule whose caption spells a track row does not draw it.
    rules: [{ feat: 'CDS', qual: 'product', val: 'x' }],
    trackKinds: { 'GC content': 'content', 'GC skew (+)': 'skew', Depth: 'depth' },
    countLegendCaptionFeatures: () => { counted += 1; return 0; }
  });
  const views = ['GC content', 'GC skew (+)', 'Depth'].map((key) => openView(created, { kind: 'legend-entry', key }));
  assert.deepEqual(views.map((view) => view.subtitle), ['Content plot', 'Skew plot', 'Depth plot']);
  views.forEach((view) => {
    assert.equal(view.colorSource, 'Changes this legend swatch only. The plot keeps the color set in Settings.');
    assert.deepEqual(view.ruleLines, []);
  });
  assert.deepEqual([views[0].title, views[0].name, views[0].dialogLabel], ['GC content', 'GC content', 'Legend item: GC content']);
  // Not counted as features while the matches prepare either.
  const preparing = openView(created, { kind: 'legend-entry', key: 'GC content' }, false);
  assert.equal(preparing.subtitle, 'Content plot');
  assert.equal(counted, 0, 'a track row counts no features');
});

// Review N6: the row kind comes from Python's keys (trackLegendRows), never from the name. A feature row named
// "GC content" holds that key, and the content track takes "GC content (2)" (D-15-65). The port is wired as
// app-setup.js `readLegendTrackKind` wires it: displayedTrackLegendRows(...).get(legendEntryKey(entry))?.kind.
test('a feature row named like a track row keeps the feature view (OV-1272)', () => {
  const pythonRows = new Map([['GC content', { color: '#ff0000' }], ['GC content (2)', { color: '#a0a0a0' }]]);
  const rows = trackLegendRows(
    [{ slotId: null, kind: 'content', label: 'GC content', positiveColor: null, negativeColor: null }], pythonRows);
  const created = createPopup({
    entries: [
      { caption: 'GC content', originalCaption: 'GC content', color: '#ff0000' },
      { caption: 'GC content (2)', originalCaption: 'GC content (2)', color: '#a0a0a0' }
    ],
    featureTypes: ['CDS'], featureCount: 3,
    readLegendTrackKind: (entry) => rows.get(legendEntryKey(entry))?.kind ?? null
  });
  const feature = openView(created, { kind: 'legend-entry', key: 'GC content' });
  const track = openView(created, { kind: 'legend-entry', key: 'GC content (2)' });
  assert.equal(feature.subtitle, '3 features');
  assert.notEqual(feature.colorSource, track.colorSource);
  assert.equal(track.subtitle, 'Content plot');
  assert.equal(track.colorSource, 'Changes this legend swatch only. The plot keeps the color set in Settings.');
});

test('the record-definition view lists the facts of the displayed Result', () => {
  const definitionGroup = new FakeElement({}, {
    children: [
      new FakeElement({ 'data-definition-line-kind': 'name' }, { text: 'Streptomyces sp.' }),
      new FakeElement({ 'data-definition-line-kind': 'subtitle' }, { text: 'neocarzinostatin' }),
      new FakeElement({ 'data-definition-line-kind': 'accession' }, { text: 'BGC0000708' }),
      new FakeElement({ 'data-definition-line-kind': 'length' }, { text: '28,450 bp' })
    ]
  });
  const view = openView(createPopup({
    definitionGroup,
    requestRecords: [{
      recordKey: 'uid-1', presentation: { reverseComplement: true }, display: { startCoordinate: 1200 }, region: { start: 10, end: 2000 }
    }],
    recordRows: [{ scope: 'linear', sourceUid: 'uid-1', recordId: 'BGC0000708', canonicalRecordKey: 'uid-1' }],
    alignment: { recordKey: 'uid-1', status: 'aligned', anchorLabel: 'neoU' }
  }), { kind: 'record-definition', recordIndex: 0, recordId: 'BGC0000708' });
  assert.deepEqual(view.facts, [
    { label: 'Accession', value: 'BGC0000708' },
    { label: 'Organism', value: 'Streptomyces sp.' },
    { label: 'Length', value: '28,450 bp' },
    { label: 'Reverse complement', value: 'Yes' },
    { label: 'Display start', value: '1,200' },
    { label: 'Region', value: '10..2,000' },
    { label: 'Alignment', value: 'Aligned on neoU' }
  ]);
  assert.equal(view.subtitle, 'neocarzinostatin · BGC0000708');
  const circular = openView(createPopup({ mode: 'circular', requestRecords: [] }),
    { kind: 'record-definition', recordIndex: 0, recordId: 'NC_012920.1' });
  assert.equal(circular.facts.find((fact) => fact.label === 'Region').value, 'Whole record');
  assert.equal(circular.facts.some((fact) => fact.label === 'Alignment'), false);
});

// OV-2328: the reason a record was not moved is its label as the UI shows it
// elsewhere, inserted unchanged; lower-casing a translated label damages
// German nouns and Turkish I.
test('a record the alignment did not move shows its reason label unchanged', () => {
  const alignmentOf = (alignment) => openView(createPopup({
    requestRecords: [{ recordKey: 'uid-1' }],
    recordRows: [{ scope: 'linear', sourceUid: 'uid-1', recordId: 'BGC0000708', canonicalRecordKey: 'uid-1' }],
    alignment: { recordKey: 'uid-1', ...alignment }
  }), { kind: 'record-definition', recordIndex: 0, recordId: 'BGC0000708' })
    .facts.find((fact) => fact.label === 'Alignment').value;
  assert.equal(alignmentOf({ status: 'skipped', rationaleLabel: 'Skipped by user' }), 'Not moved (Skipped by user)');
  assert.equal(alignmentOf({ status: 'skipped', rationaleLabel: '' }), 'Not moved (skipped)');
});

// OV-502: only a Linear definition marks which line is the record's name. A
// Circular definition's name lines may be an organism, a strain, a replicon or
// a Record label, wrapped or not: the popup shows them as its heading, as the
// Result draws them, and lists none of them as a fact.
test('the record popup names a definition line as Organism only where the Result marks it', () => {
  const circular = (/** @type {string[]} */ lines, label = '') => openView(createPopup({
    mode: 'circular',
    requestRecords: [{ recordKey: 'uid-1', presentation: label ? { label } : {} }],
    definitionGroup: new FakeElement({}, {
      children: lines.map((line) => new FakeElement({ tag: 'text' }, { text: line }))
    })
  }), { kind: 'record-definition', recordIndex: 0, recordId: 'NC_012920.1' });
  const factLabels = (/** @type {{ facts: { label: string }[] }} */ view) => view.facts.map((fact) => fact.label);
  const circularFacts = ['Accession', 'Length', 'Reverse complement', 'Display start', 'Region'];

  const labelled = circular(['Human mitochondrion', 'NC_012920.1', '16,569 bp', '44.36% GC'], '<i>Human</i> mitochondrion');
  assert.deepEqual(factLabels(labelled), circularFacts);
  assert.equal(labelled.title, 'Human mitochondrion');
  assert.equal(labelled.subtitle, 'NC_012920.1');

  const wrapped = circular(['Candidatus Pelagibacter', 'ubique HTCC1062', 'NC_012920.1', '16,569 bp'],
    'Candidatus Pelagibacter ubique HTCC1062');
  assert.deepEqual(factLabels(wrapped), circularFacts);
  assert.equal(wrapped.title, 'Candidatus Pelagibacter');
  assert.equal(wrapped.subtitle, 'ubique HTCC1062 · NC_012920.1');

  const replicon = circular(['Chromosome 1', 'NC_012920.1', '16,569 bp', '44.36% GC']);
  assert.deepEqual(factLabels(replicon), circularFacts);
  assert.equal(replicon.title, 'Chromosome 1');

  // The record ID alone (a Multi-Record Canvas record without a replicon), also
  // when the request holds a Record label the definition does not draw.
  const idOnly = circular(['NC_012920.1', '16,569 bp', '44.36% GC'], 'Not drawn');
  assert.deepEqual(factLabels(idOnly), circularFacts);
  assert.equal(idOnly.title, 'NC_012920.1');
  assert.equal(idOnly.subtitle, '');

  const linear = openView(createPopup({
    requestRecords: [{ recordKey: 'uid-1', presentation: { label: 'Typed name' } }],
    recordRows: [{ scope: 'linear', sourceUid: 'uid-1', recordId: 'R1', canonicalRecordKey: 'uid-1' }],
    definitionGroup: new FakeElement({}, {
      children: [
        new FakeElement({ 'data-definition-line-kind': 'name' }, { text: 'Typed name' }),
        new FakeElement({ 'data-definition-line-kind': 'accession' }, { text: 'R1' })
      ]
    })
  }), { kind: 'record-definition', recordIndex: 0, recordId: 'R1' });
  assert.deepEqual(linear.facts[1], { label: 'Organism', value: 'Typed name' });
  const linearIdOnly = openView(createPopup({
    recordRows: [{ scope: 'linear', sourceUid: 'uid-1', recordId: 'R1', canonicalRecordKey: 'uid-1' }],
    definitionGroup: new FakeElement({}, {
      children: [new FakeElement({ 'data-definition-line-kind': 'accession' }, { text: 'R1' })]
    })
  }), { kind: 'record-definition', recordIndex: 0, recordId: 'R1' });
  assert.equal(linearIdOnly.facts.some((fact) => fact.label === 'Organism'), false);
  assert.equal(linearIdOnly.title, 'R1');
});

test('a Circular definition without line kinds reads its name, accession, and length lines', () => {
  const definitionGroup = new FakeElement({}, {
    children: ['Homo sapiens', 'mitochondrion', 'NC_012920.1', '16,569 bp', '44.36% GC']
      .map((line) => new FakeElement({ tag: 'text' }, { text: line }))
  });
  const view = openView(createPopup({ mode: 'circular', requestRecords: [], definitionGroup }),
    { kind: 'record-definition', recordIndex: 0, recordId: 'NC_012920.1' });
  assert.equal(view.title, 'Homo sapiens');
  assert.equal(view.subtitle, 'mitochondrion · NC_012920.1');
  assert.equal(view.facts.find((fact) => fact.label === 'Length').value, '16,569 bp');
});
