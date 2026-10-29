// PD-OI-103 rev 2 (U2 CONTROL GROUPS), PD-OI-059 rev 3, D-150, D-151: the G2
// units go live. A committed value of a G2 row reaches the Result through the
// automatic rerender: its writes are the committed request with only the
// unit's leaves replaced by what a Generate of the same draft sends (R3), and
// every other pending draft value, of another unit too, stays out (N-16).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};

const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const config = await import('../../gbdraw/web/js/services/config.js');
const request = await import('../../gbdraw/web/js/services/session-request.js');
const owner = await import('../../gbdraw/web/js/services/setting-application.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { createAnnotationSet } = await import('../../gbdraw/web/js/services/annotation-state.js');
const { createDefaultLinearTrackSlots } = await import('../../gbdraw/web/js/services/linear-track-slot-model.js');
const { RECORD_POPUP_FIELDS } = await import('../../gbdraw/web/js/app/preview-item-popup.js');

const record = (id) => `LOCUS       ${id.padEnd(16)}          40 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Live setting fixture ${id}.
ACCESSION   ${id}
VERSION     ${id}
FEATURES             Location/Qualifiers
     CDS             1..30
                     /product="test protein"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
const genbank = (ids) => {
  const text = ids.map(record).join('');
  return { name: 'live.gb', type: 'text/plain', size: new TextEncoder().encode(text).byteLength, lastModified: 0, data: btoa(text) };
};
const seq = (uid, id) => ({
  uid, gb: genbank([id]), depth: null, losat_gencode: 1, region_record_id: '', region_start: null, region_end: null,
  region_reverse: false, definition: '', record_subtitle: '', file_definition: `File ${id}`, file_subtitle: ''
});
const inputs = { linearSeqs: [seq('record-1', 'LIVE1'), seq('record-2', 'LIVE2')] };

// The committed Session a Generate of `mode` builds from that mode's drawing now.
const generate = (mode) => {
  const drawing = modeDrawing(state, mode);
  config.installShownDrawing(mode);
  modeDrawing(state, 'circular').inputType.value = 'gb';
  modeDrawing(state, 'linear').inputType.value = 'gb';
  const filesData = mode === 'circular'
    ? { c_gb: genbank(['LIVE1']), linearSeqs: [] }
    : { linearSeqs: inputs.linearSeqs, linearComparisons: [] };
  const comparisonPlanSnapshot = mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
        losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp.mode
      })
    : null;
  const built = request.buildCanonicalRenderRequest({ state, drawing, filesData, comparisonPlanSnapshot });
  return { renderRequest: built.renderRequest, resources: built.resources, webFiles: built.webFiles || {} };
};
const project = (committed, mode, requestPaths) => request.projectLiveSettingWrites({
  committed, drawing: modeDrawing(state, mode), requestPaths, inputs
});
// The committed Session after an apply: the request writes in its request, the
// Web file writes in its Web file record (F-2).
const applied = (committed, writes) => {
  const renderRequest = structuredClone(committed.renderRequest);
  request.applyLiveSettingWrites(renderRequest, writes);
  return { ...committed, renderRequest, webFiles: request.applyLiveWebFileWrites(committed.webFiles, writes) };
};
const at = (root, path) => path.reduce((node, key) => node?.[key], root);
const put = (root, path, value) => {
  const parent = at(root, path.slice(0, -1));
  if (value === undefined) delete parent[path.at(-1)];
  else parent[path.at(-1)] = structuredClone(value);
};
const restore = [];
const set = (object, key, value) => {
  const had = Object.hasOwn(object, key);
  const old = object[key];
  restore.push(() => { if (had) object[key] = old; else delete object[key]; });
  object[key] = value;
};
const undoAll = () => { while (restore.length) restore.pop()(); };

const DEF = 'diagramOptions.configOverrides.objects.definition.font_size';
const STYLES = 'diagramOptions.configOverrides.objects.definition.linear.line_styles';
const SPECIES = 'diagramOptions.species';
const WINDOW = 'diagramOptions.window';
const ANNOTATIONS = 'diagramOptions.annotations.sets';
const HEIGHT = 'diagramOptions.configOverrides.canvas.linear.comparison_height';
const RADIUS = 'diagramOptions.tracks.centerReservedRadius';
const SLOTS = 'diagramOptions.tracks.linearTrackSlots';
const RECORDS = 'records';
const linearOverride = (key) => [['renderRequest', 'diagramOptions', 'configOverrides', `${key}.short`],
  ['renderRequest', 'diagramOptions', 'configOverrides', `${key}.long`]];

// One annotation set with a region item, and the custom stack with a GC row,
// so that every G2 unit has a committed value to replace.
const fixtures = (mode) => {
  const drawing = modeDrawing(state, mode);
  // The annotation sets are the slice's array: the fixture replaces its items.
  const sets = drawing.annotationSets;
  const before = sets.splice(0, sets.length);
  restore.push(() => sets.splice(0, sets.length, ...before));
  sets.push(...[createAnnotationSet({
    id: 'islands', annotations: [
      { id: 'r1', target: { kind: 'coordinateSpan', start: 2, end: 20 }, mark: 'band' },
      { id: 'r2', target: { kind: 'coordinateSpan', start: 22, end: 30 }, mark: 'band', legendLabel: 'Island' }
    ]
  })]);
  if (mode === 'linear') {
    set(drawing.adv, 'linear_track_slots_enabled', true);
    set(drawing.adv, 'linear_track_slots', createDefaultLinearTrackSlots({ showGc: true }));
  }
  set(drawing.adv, 'linear_definition_line_styles', {});
  // A single Circular record draws its record presentation.
  if (mode === 'circular') set(drawing.form, 'multi_record_canvas', false);
};

// Each unit: its request paths, the draft change, and the committed-Session
// leaves a Generate of the change writes (paths from the Session root).
const UNITS = {
  'definition-text': {
    circular: { paths: [DEF, STYLES], change: (d) => set(d.adv, 'def_font_size', 23),
      leaves: () => [['renderRequest', 'diagramOptions', 'configOverrides', 'objects.definition.circular.font_size']] },
    linear: { paths: [DEF, STYLES], change: (d) => {
      set(d.adv, 'def_font_size', 23);
      set(d.adv, 'linear_definition_line_styles', { name: { font_size: 15, font_weight: 'bold', fill: '#123456' } });
    }, leaves: () => [...linearOverride('objects.definition.linear.font_size'),
      ...['font_size', 'font_weight', 'fill'].map((field) => ['renderRequest', 'diagramOptions', 'configOverrides',
        `objects.definition.linear.line_styles.name.${field}`])] }
  },
  species: {
    circular: { paths: [SPECIES], change: (d) => set(d.form, 'species', 'Live species'),
      leaves: () => [['renderRequest', 'diagramOptions', 'species']] }
  },
  'dinucleotide-window': Object.fromEntries(['circular', 'linear'].map((mode) => [mode, {
    paths: [WINDOW], change: (d) => set(d.adv, 'window_size', 7), leaves: () => [['renderRequest', 'diagramOptions', 'window']]
  }])),
  'annotation-style': Object.fromEntries(['circular', 'linear'].map((mode) => [mode, {
    paths: [ANNOTATIONS], change: (d) => {
      const item = d.annotationSets[0].annotations[0];
      set(item, 'style', { fill: '#ff0000', stroke: '#00ff00' });
    }, leaves: () => [['renderRequest', 'diagramOptions', 'annotations', 'sets', 0, 'annotations', 0, 'style']],
    // The same item's Hatch and label apply on Generate.
    samePending: (d) => {
      const item = d.annotationSets[0].annotations[0];
      set(item, 'style', { ...item.style, hatch: { angle: 45, spacing: 4, color: '#000000', width: 1, cross: false } });
      set(item, 'label', 'Pending label');
    }
  }])),
  'comparison-height': {
    linear: { paths: [HEIGHT], change: (d) => set(d.adv, 'comparison_height', 31), leaves: () => [] }
  },
  'track-geometry': {
    circular: { paths: [RADIUS, SLOTS], change: (d) => set(d.adv, 'center_reserved_radius', 120),
      leaves: () => [['renderRequest', 'diagramOptions', 'tracks', 'centerReservedRadius']] },
    linear: { paths: [RADIUS, SLOTS], change: (d) => {
      set(d.adv.linear_track_slots[1], 'height', '44');
    }, leaves: () => [['renderRequest', 'diagramOptions', 'tracks', 'linearTrackSlots', 1, 'height']],
    // The same row's spacing and renderer apply on Generate.
    samePending: (d) => {
      set(d.adv.linear_track_slots[1], 'spacing', '9');
      set(d.adv.linear_track_slots[1], 'renderer', 'dinucleotide_skew');
    } }
  },
  'record-names': {
    circular: { paths: [RECORDS], change: (d) => {
      set(d.form, 'circular_record_label', 'Live label');
      set(d.form, 'circular_record_subtitle', 'Live subtitle');
    }, leaves: () => [['renderRequest', 'records', 0, 'presentation', 'label'], ['renderRequest', 'records', 0, 'presentation', 'subtitle']] },
    linear: { paths: [RECORDS], change: () => set(inputs.linearSeqs[1], 'definition', 'Live definition'),
      leaves: () => [['renderRequest', 'records', 1, 'presentation', 'label'], ['webFiles', 'linearRecordMetadata', 1, 'recordDefinition']] }
  }
};
// Draft values that apply on Generate (G3 and generate rows) and stay pending.
const pendingGenerateOnly = (mode) => {
  const drawing = modeDrawing(state, mode);
  set(drawing.adv, 'label_font_size', 29);
  set(drawing.adv, 'step_size', 3);
  if (mode === 'circular') {
    set(drawing.form, 'strain', 'Pending strain');
    set(drawing.form, 'circular_reverse', !drawing.form.circular_reverse);
  } else {
    set(inputs.linearSeqs[1], 'record_subtitle', 'Pending subtitle');
    set(inputs.linearSeqs[1], 'region_reverse', true);
    set(inputs.linearSeqs[1], 'losat_gencode', 11);
  }
};

test('the G2 rows are live, and each request path has its projector', () => {
  const g2 = ['definition-text', 'species', 'dinucleotide-window', 'annotation-style', 'comparison-height',
    'track-geometry', 'record-names'];
  for (const id of g2) {
    const row = owner.SETTING_APPLICATION.find((candidate) => candidate.id === id);
    assert.equal(row?.application, 'live', id);
    for (const path of Object.values(row.requestPaths)) assert.ok(request.LIVE_SETTING_REQUEST_PATHS.includes(path), path);
  }
  const paths = (id) => owner.SETTING_APPLICATION.find((row) => row.id === id).requestPaths;
  assert.deepEqual(paths('record-names'), {
    'files.linearSeqs': RECORDS, 'files.linearSeqs.*.definition': RECORDS,
    'config.form.circular_record_label': RECORDS, 'config.form.circular_record_subtitle': RECORDS
  });
  // A History step of an annotation item's fill, a Linear slot's height and a
  // record's definition requests the unit (`*` matches an index).
  const step = (...path) => owner.liveSettingRequestPaths([{ path }]);
  assert.deepEqual(step('config', 'annotationSets', 0, 'annotations', 2, 'style', 'fill'), [ANNOTATIONS]);
  assert.deepEqual(step('config', 'adv', 'linear_track_slots', 3, 'height'), [SLOTS]);
  assert.deepEqual(step('files', 'linearSeqs', 1, 'definition'), [RECORDS]);
  assert.deepEqual(step('config', 'form', 'strain'), []);
  assert.deepEqual(step('config', 'adv', 'step_size'), []);
  // Strain and the GC step stay on Generate.
  assert.equal(owner.settingApplicationNote('species-strain').application, 'generate');
  assert.equal(owner.settingApplicationNote('dinucleotide').application, 'generate');
});

for (const [unit, modes] of Object.entries(UNITS)) {
  for (const [mode, spec] of Object.entries(modes)) {
    test(`${unit} (${mode}): the live writes equal a Generate of the draft, and other pending values stay out`, () => {
      try {
        fixtures(mode);
        if (unit === 'comparison-height') {
          // A Python-written request with pairs carries no match style (F-10).
          const drawing = modeDrawing(state, 'linear');
          set(drawing.adv, 'comparison_height', '');
        }
        const committed = generate(mode);
        if (unit === 'comparison-height') {
          committed.renderRequest.comparisons = [{ kind: 'pairwiseComparison', pairs: [{ query: 0, subject: 1 }] }];
        }
        spec.change(modeDrawing(state, mode));
        const draft = generate(mode);
        // Only this unit's change pending: the live Session equals the Generate.
        const live = applied(committed, project(committed, mode, request.LIVE_SETTING_REQUEST_PATHS));
        if (unit === 'comparison-height') {
          assert.equal(live.renderRequest.diagramOptions.configOverrides['canvas.linear.comparison_height'], 31);
        } else {
          assert.deepEqual(live.renderRequest, draft.renderRequest);
          assert.deepEqual(live.webFiles, draft.webFiles);
          for (const leaf of spec.leaves()) assert.notDeepEqual(at(committed, leaf), at(draft, leaf), leaf.join('.'));
        }
        // Every other unit's change and the Generate-only values pending too:
        // this unit's projection writes only its own leaves.
        for (const [other, otherModes] of Object.entries(UNITS)) {
          if (other !== unit && otherModes[mode]) otherModes[mode].change(modeDrawing(state, mode));
        }
        pendingGenerateOnly(mode);
        spec.samePending?.(modeDrawing(state, mode));
        const busy = applied(committed, project(committed, mode, spec.paths));
        const expected = structuredClone(committed);
        if (unit === 'comparison-height') {
          expected.renderRequest.diagramOptions.configOverrides['canvas.linear.comparison_height'] = 31;
        }
        for (const leaf of spec.leaves()) put(expected, leaf, at(draft, leaf));
        assert.deepEqual(busy.renderRequest, expected.renderRequest);
        assert.deepEqual(busy.webFiles, expected.webFiles);
      } finally {
        undoAll();
      }
    });
  }
}

test('a value Generate refuses fails the apply (INPUT_INVALID)', () => {
  try {
    fixtures('linear');
    const committed = generate('linear');
    const refused = (path) => assert.throws(() => project(committed, 'linear', [path]), (error) => error.code === 'INPUT_INVALID', path);
    set(modeDrawing(state, 'linear').adv.linear_track_slots[1], 'height', '-4');
    refused(SLOTS);
    set(modeDrawing(state, 'linear').adv, 'window_size', 'abc');
    refused(WINDOW);
    committed.renderRequest.comparisons = [{ kind: 'pairwiseComparison', pairs: [] }];
    set(modeDrawing(state, 'linear').adv, 'comparison_height', -1);
    refused(HEIGHT);
  } finally {
    undoAll();
  }
});

test('a features or Depth row, an unknown slot, a new item, and an item with a Legend label keep their committed values', () => {
  try {
    fixtures('linear');
    const drawing = modeDrawing(state, 'linear');
    const committed = generate('linear');
    set(drawing.adv.linear_track_slots[0], 'height', '50');
    drawing.adv.linear_track_slots.push({ id: 'new_row', renderer: 'dinucleotide_skew', height: '12', enabled: true });
    restore.push(() => drawing.adv.linear_track_slots.pop());
    assert.deepEqual(project(committed, 'linear', [SLOTS]), []);
    // O-8: a Legend label draws the item's fill in a Legend row, so that item waits.
    set(drawing.annotationSets[0].annotations[1], 'style', { fill: '#ff0000' });
    drawing.annotationSets[0].annotations.push({ id: 'r3', target: { kind: 'coordinateSpan', start: 1, end: 3 }, style: { fill: '#000000' } });
    restore.push(() => drawing.annotationSets[0].annotations.pop());
    assert.deepEqual(project(committed, 'linear', [ANNOTATIONS]), []);
    // A set Default fill is not live: an item that follows it keeps its drawn colors.
    set(drawing.annotationSets[0], 'defaultStyle', { ...drawing.annotationSets[0].defaultStyle, fill: '#abcdef' });
    assert.deepEqual(project(committed, 'linear', [ANNOTATIONS]), []);
  } finally {
    undoAll();
  }
});

test('a null committed item style takes the set default with the new colors', () => {
  try {
    fixtures('circular');
    const committed = generate('circular');
    const item = committed.renderRequest.diagramOptions.annotations.sets[0].annotations[0];
    assert.equal(item.style, null);
    set(modeDrawing(state, 'circular').annotationSets[0].annotations[0], 'style', { fill: '#ff0000' });
    const draft = generate('circular');
    const live = applied(committed, project(committed, 'circular', [ANNOTATIONS]));
    assert.deepEqual(live.renderRequest.diagramOptions.annotations, draft.renderRequest.diagramOptions.annotations);
  } finally {
    undoAll();
  }
});

test('Match height applies to a committed comparison only, LOSAT-only intent and Python-written pairs included', () => {
  try {
    const drawing = modeDrawing(state, 'linear');
    const committed = generate('linear');
    assert.equal(Object.hasOwn(committed.renderRequest.diagramOptions, 'pairwiseMatchStyle'), false);
    set(drawing.adv, 'comparison_height', 31);
    // No comparison intent: Match height waits for Generate.
    assert.equal(request.liveSettingApplies(committed.renderRequest, HEIGHT, drawing, inputs), false);
    assert.deepEqual(project(committed, 'linear', [HEIGHT]), []);
    // The Web's comparison intent (pairs or LOSAT) sends the match style.
    const web = structuredClone(committed);
    web.renderRequest.diagramOptions.pairwiseMatchStyle = 'ribbon';
    assert.equal(request.liveSettingApplies(web.renderRequest, HEIGHT, drawing, inputs), true);
    // A Python-written request omits the default style but carries its pairs.
    const python = structuredClone(committed);
    python.renderRequest.comparisons = [{ kind: 'pairwiseComparison', pairs: [] }];
    assert.equal(request.liveSettingApplies(python.renderRequest, HEIGHT, drawing, inputs), true);
    assert.equal(project(python, 'linear', [HEIGHT]).length, 1);
    // A CLI request without -b or a protein mode writes a disabled protein
    // pipeline that draws no comparison: Match height waits.
    const cli = structuredClone(committed);
    cli.renderRequest.comparisons = [{ kind: 'generatedProteinComparison', mode: 'none', pairs: [] }];
    assert.equal(request.liveSettingApplies(cli.renderRequest, HEIGHT, drawing, inputs), false);
    assert.deepEqual(project(cli, 'linear', [HEIGHT]), []);
  } finally {
    undoAll();
  }
});

test('Circular record names apply to a single record only; a batch or grid of one waits', () => {
  const drawing = modeDrawing(state, 'circular');
  for (const grouping of ['batch', 'grid']) {
    assert.equal(request.liveSettingApplies({ mode: 'circular', grouping, records: [] }, RECORDS, drawing, inputs), false, grouping);
  }
  assert.equal(request.liveSettingApplies({ mode: 'circular', grouping: 'single', records: [] }, RECORDS, drawing, inputs), true);
});

test('Linear record names follow the draft records in order; a removed or moved record waits, and split members each take the label', () => {
  try {
    const committed = generate('linear');
    const applies = (seqs) => request.liveSettingApplies(committed.renderRequest, RECORDS, modeDrawing(state, 'linear'),
      { linearSeqs: seqs }, committed.resources);
    const [first, second] = inputs.linearSeqs;
    assert.equal(applies([first, second]), true);
    assert.equal(applies([first, second, seq('record-3', 'LIVE3')]), true, 'an added record leaves the drawn ones in order');
    assert.equal(applies([second, first]), false, 'moved');
    assert.equal(applies([first]), false, 'removed');
    // A replaced File: the committed resource of the record does not hold it,
    // whatever its definition says (D-193 D).
    assert.equal(applies([first, { ...second, gb: genbank(['LIVE2']) }]), false, 'replaced by a File of the same definition');
    assert.equal(applies([first, { ...second, file_definition: 'Another file' }]), true, 'the same File');
    assert.equal(request.liveSettingApplies(committed.renderRequest, RECORDS, modeDrawing(state, 'linear'), inputs), false,
      'no committed resources name the File');
    // F-18: Generate splits a card into `<uid>:<k>` members, each with its own
    // metadata entry; each takes the card's label.
    // The resources keep their payload identity (a clone would not).
    const split = { ...structuredClone({ renderRequest: committed.renderRequest, webFiles: committed.webFiles }), resources: committed.resources };
    const card = split.renderRequest.records[1];
    split.renderRequest.records.splice(1, 1, { ...card, recordKey: 'record-2:1' }, { ...card, recordKey: 'record-2:2' });
    const entry = split.webFiles.linearRecordMetadata[1];
    split.webFiles.linearRecordMetadata.splice(1, 1, { ...entry, recordKey: 'record-2:1' }, { ...entry, recordKey: 'record-2:2' });
    set(second, 'definition', 'Split label');
    const live = applied(split, project(split, 'linear', [RECORDS]));
    assert.deepEqual(live.renderRequest.records.map((r) => r.presentation.label), ['File LIVE1', 'Split label', 'Split label']);
    assert.deepEqual(live.webFiles.linearRecordMetadata.map((e) => e.recordDefinition), ['', 'Split label', 'Split label']);
    assert.deepEqual(live.webFiles.linearRecordMetadata.map((e) => e.recordKey), ['record-1', 'record-2:1', 'record-2:2']);
  } finally {
    undoAll();
  }
});

test('a Linear card\'s record selection applies on Generate: its record names wait until then', () => {
  try {
    // A two-record File whose card draws its first record by the definition
    // inferred from that record (no File or record Definition).
    const card = { ...seq('record-2', 'LIVE2'), gb: genbank(['LIVE2', 'LIVE3']), file_definition: '',
      region_record_id: '#1', inferred_definition: 'Live setting fixture LIVE2.' };
    set(inputs.linearSeqs, 1, card);
    const committed = generate('linear');
    const drawing = modeDrawing(state, 'linear');
    const note = () => owner.settingApplicationNote('record-names', (path) => request.liveSettingNoteApplies(
      committed.renderRequest, path, drawing, inputs, committed.resources)).text;
    assert.equal(note(), 'Live edit');
    // The selector picks the second record; the inferred definition follows it
    // (setLinearRecordSelector). Generate draws the second record with it.
    set(card, 'region_record_id', '#2');
    set(card, 'inferred_definition', 'Live setting fixture LIVE3.');
    const draft = generate('linear');
    assert.notDeepEqual(draft.renderRequest.records[1].selector, committed.renderRequest.records[1].selector);
    const live = applied(committed, project(committed, 'linear', request.LIVE_SETTING_REQUEST_PATHS));
    assert.deepEqual(live.renderRequest.records, committed.renderRequest.records, 'the drawn record keeps its name');
    assert.deepEqual(live.webFiles, committed.webFiles);
    assert.equal(note(), 'Applies on Generate');
    // The same record by a crop of it: the card still selects it.
    set(card, 'region_record_id', '#1');
    set(card, 'region_start', 2);
    set(card, 'region_end', 20);
    set(card, 'definition', 'Cropped name');
    assert.equal(note(), 'Live edit');
    assert.equal(applied(committed, project(committed, 'linear', [RECORDS])).renderRequest.records[1].presentation.label,
      'Cropped name');
    // A split member (F-18) was drawn from a card that selects no record.
    const split = { ...structuredClone({ renderRequest: committed.renderRequest, webFiles: committed.webFiles }), resources: committed.resources };
    split.renderRequest.records.splice(1, 1, { ...split.renderRequest.records[1], recordKey: 'record-2:1',
      selector: { kind: 'recordIndex', index: 0 } });
    split.webFiles.linearRecordMetadata.splice(1, 1, { ...split.webFiles.linearRecordMetadata[1], recordKey: 'record-2:1' });
    const splitApplies = () => request.liveSettingApplies(split.renderRequest, RECORDS, drawing, inputs, split.resources);
    assert.equal(splitApplies(), false, 'the card now selects one record');
    set(card, 'region_record_id', '');
    set(card, 'region_start', null);
    set(card, 'region_end', null);
    assert.equal(splitApplies(), true);
  } finally {
    undoAll();
  }
});

test('a live write creates an optional object the committed request holds as null', () => {
  const target = { diagramOptions: { configOverrides: null } };
  request.applyLiveSettingWrites(target, [[['diagramOptions', 'configOverrides', 'objects.definition.font_size'], 23]]);
  assert.deepEqual(target.diagramOptions, { configOverrides: { 'objects.definition.font_size': 23 } });
  const removal = { diagramOptions: { configOverrides: null } };
  request.applyLiveSettingWrites(removal, [[['diagramOptions', 'configOverrides', 'canvas.linear.comparison_height'], undefined]]);
  assert.deepEqual(removal.diagramOptions, { configOverrides: null }, 'a removal has nothing to remove');
});

test('a live write never enters an array by index (F-1)', () => {
  const target = { records: [{ presentation: {} }], diagramOptions: { window: 5 } };
  assert.throws(() => request.applyLiveSettingWrites(target, [[['records', '0', 'presentation', 'label'], 'x']]));
  assert.throws(() => request.applyLiveSettingWrites(target, [[['diagramOptions', 'window', 'x'], 1]]), undefined, 'a scalar is no parent');
  request.applyLiveSettingWrites(target, [[['diagramOptions', 'tracks', 'centerReservedRadius'], 3]]);
  assert.deepEqual(target.diagramOptions, { window: 5, tracks: { centerReservedRadius: 3 } });
});

test('a row with paths of the other mode shows Live edit in each mode', () => {
  const circular = { mode: 'circular', grouping: 'single', records: [] };
  const linear = { mode: 'linear', grouping: 'single', records: [], diagramOptions: {}, comparisons: [] };
  const note = (id, committed) => owner.settingApplicationNote(id, (path) => request.liveSettingNoteApplies(
    committed, path, modeDrawing(state, committed.mode), inputs));
  for (const id of ['definition-text', 'track-geometry', 'dinucleotide-window']) {
    assert.equal(note(id, circular).application, 'live', `${id} circular`);
  }
  assert.equal(note('definition-text', circular).text, 'Live edit');
  assert.equal(note('track-geometry', linear).application, 'live');
  assert.equal(note('definition-text', linear).text, 'Live edit');
  assert.equal(note('comparison-height', linear).application, 'generate', 'no committed comparison');
  // The other mode's path is still not applied.
  assert.equal(request.liveSettingApplies(circular, STYLES, modeDrawing(state, 'circular'), inputs), false);
});

test('an annotation item\'s and a Linear row\'s note read the condition their apply reads (D-193 A2)', () => {
  try {
    fixtures('linear');
    const drawing = modeDrawing(state, 'linear');
    const items = drawing.annotationSets[0].annotations;
    const committed = generate('linear');
    // r4 is new since the Generate.
    items.push({ id: 'r4', target: { kind: 'coordinateSpan', start: 4, end: 6 } });
    restore.push(() => items.pop());
    const note = (id, item, shown = committed) => owner.settingApplicationNote(id, (path) => request.liveSettingNoteApplies(
      shown.renderRequest, path, modeDrawing(state, shown.renderRequest.mode), inputs, shown.resources, item)).text;
    // Without an item the rows say Live edit; each item says what its apply does.
    assert.equal(note('annotation-style', null), 'Live edit');
    assert.equal(note('track-geometry', null), 'Live edit');
    for (const [itemId, text] of Object.entries({ r1: 'Live edit', r2: 'Applies on Generate', r4: 'Applies on Generate' })) {
      assert.equal(note('annotation-style', { setId: 'islands', itemId }), text, itemId);
      // The note says Live edit exactly when the apply writes the item's new colors.
      const item = items.find((candidate) => candidate.id === itemId);
      const style = item.style;
      item.style = { fill: '#ff0000', stroke: '#00ff00' };
      assert.equal(project(committed, 'linear', [ANNOTATIONS]).length, text === 'Live edit' ? 1 : 0, itemId);
      item.style = style;
    }
    assert.equal(note('annotation-style', { setId: 'another', itemId: 'r1' }), 'Applies on Generate', 'a set the Result lacks');
    // A set Legend label draws every item of the set in a Legend row (O-8).
    const labelled = structuredClone(committed);
    labelled.renderRequest.diagramOptions.annotations.sets[0].legendLabel = 'Islands';
    assert.equal(note('annotation-style', { setId: 'islands', itemId: 'r1' }, labelled), 'Applies on Generate');
    // A Linear row: the GC Height is live; the features row and a new or renamed row wait.
    for (const [slotId, text] of Object.entries({ gc_content: 'Live edit', features: 'Applies on Generate', new_row: 'Applies on Generate' })) {
      assert.equal(note('track-geometry', { slotId }), text, slotId);
    }
    const stringSlots = structuredClone(committed);
    stringSlots.renderRequest.diagramOptions.tracks.linearTrackSlots = ['gc_content:dinucleotide_content'];
    assert.equal(note('track-geometry', { slotId: 'gc_content' }, stringSlots), 'Applies on Generate', 'a spec string row');
    // A Circular Result draws no Linear row and waits for none (the row rule).
    assert.equal(note('track-geometry', { slotId: 'features' }, { renderRequest: { mode: 'circular', grouping: 'single', records: [] } }), 'Live edit');
  } finally {
    undoAll();
  }
});

test('the record popup fields name their setting row and its draft path', () => {
  const rows = new Map(owner.SETTING_APPLICATION.map((row) => [row.id, row]));
  const drafts = { linear: 'files.linearSeqs.*', circularRecord: 'config.form', circularDrawing: 'config.form' };
  for (const [target, fields] of Object.entries(RECORD_POPUP_FIELDS)) {
    for (const [field, { key, setting }] of Object.entries(fields)) {
      const row = rows.get(setting);
      assert.ok(row, `${target}.${field}: ${setting}`);
      if (row.application === 'live') {
        assert.ok(Object.hasOwn(row.requestPaths, `${drafts[target]}.${key}`), `${target}.${field}: ${key} in ${setting}`);
      } else {
        assert.equal(row.application, 'generate', `${target}.${field}`);
      }
    }
  }
  assert.deepEqual(Object.fromEntries(Object.entries(RECORD_POPUP_FIELDS).map(([target, fields]) => [target,
    Object.fromEntries(Object.entries(fields).map(([field, { setting }]) => [field, setting]))])), {
    linear: { name: 'record-names', subtitle: 'inputs', reverse: 'inputs' },
    circularRecord: { name: 'record-names', subtitle: 'record-names', reverse: 'inputs' },
    circularDrawing: { name: 'species', subtitle: 'species-strain', reverse: 'inputs' }
  });
});
