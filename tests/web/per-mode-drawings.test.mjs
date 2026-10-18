// Per-mode drawings (plan OV-80 §8 "Per-mode values", §5 Depth series count):
// a value set in one mode's drawing stays out of the other drawing, and each
// mode's request carries its own drawing's values.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readdirSync } from 'node:fs';
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

const { state, createDefaultAdv, createDefaultForm } = await import('../../gbdraw/web/js/state.js');
const { installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { MODE_SCOPED_SETTINGS } = await import('../../gbdraw/web/js/mode-scoped-settings.generated.js');
const MODE_SCOPED_ROWS = MODE_SCOPED_SETTINGS.rows;
const { DRAWING_SLICES, SLICE_READERS } = await import('./helpers/drawing-state.mjs');
const { createAnnotationSet } = await import('../../gbdraw/web/js/services/annotation-state.js');

const MODES = ['circular', 'linear'];
const genbankText = `LOCUS       PERMODE                    40 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Per-mode drawing fixture.
ACCESSION   PERMODE
VERSION     PERMODE
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
     CDS             1..30
                     /product="test protein"
ORIGIN
        1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//
`;
const file = (name, text) => ({
  name, type: 'text/plain', size: new TextEncoder().encode(text).byteLength, lastModified: 0, data: btoa(text)
});
const genbank = () => file('per-mode.gb', genbankText);
const depthFile = (name) => file(name, 'position\tdepth\n1\t10\n2\t12\n');

const resetDrawing = (mode) => {
  const drawing = modeDrawing(state, mode);
  Object.assign(drawing.form, createDefaultForm());
  Object.assign(drawing.adv, createDefaultAdv(mode));
  drawing.linearRecordLayoutEnabled.value = false;
  drawing.linearComparisonPlan.mode = 'none';
  drawing.linearComparisonPlan.edges.splice(0);
};
// Generate builds the request of the shown mode from that mode's drawing.
const requestOf = (mode, { linearDepth = null } = {}) => {
  const drawing = modeDrawing(state, mode);
  installShownDrawing(mode);
  state.cInputType.value = 'gb';
  state.lInputType.value = 'gb';
  state.circularRecordList.value = [];
  const filesData = mode === 'circular'
    ? { c_gb: genbank(), linearSeqs: [] }
    : {
        linearSeqs: [{
          uid: 'record-1', gb: genbank(), depth: linearDepth, losat_gencode: 1, region_record_id: '',
          region_start: null, region_end: null, region_reverse: false
        }],
        linearComparisons: []
      };
  const comparisonPlanSnapshot = mode === 'linear'
    ? resolveLinearComparisonPlan({
        plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
        losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp.mode
      })
    : null;
  return buildCanonicalRenderRequest({ state, drawing, filesData, comparisonPlanSnapshot });
};
const SAMPLE_PATHS = {
  circular: ['labels.font_size.long', 'objects.axis.circular.stroke_width.long', 'canvas.strandedness',
    'objects.gc_content.mode', 'objects.depth.min_depth'],
  linear: ['labels.font_size.linear.long', 'objects.axis.linear.stroke_width.long', 'canvas.strandedness',
    'objects.gc_content.mode', 'objects.depth.min_depth']
};
const samples = (mode, canonical) => SAMPLE_PATHS[mode]
  .map((path) => canonical.renderRequest.diagramOptions.configOverrides[path]);

test('no object of one drawing is shared with the other drawing', () => {
  // A shared nested object (a track slot list, a Legend map) would carry an
  // edit across modes; every registry row then holds a value per mode.
  const seen = new Map();
  const walk = (value, path, mode) => {
    if (!value || typeof value !== 'object') return;
    const raw = Object.hasOwn(value, 'value') && Object.keys(value).length === 1 ? value.value : value;
    if (!raw || typeof raw !== 'object' || Object.isFrozen(raw)) return;
    const other = seen.get(raw);
    assert.ok(!other || other.mode === mode, `${path} (${mode}) is ${other?.path} (${other?.mode})`);
    if (other) return;
    seen.set(raw, { path, mode });
    for (const [key, child] of Object.entries(raw)) walk(child, `${path}.${key}`, mode);
  };
  for (const mode of MODES) {
    for (const [key, member] of Object.entries(modeDrawing(state, mode))) walk(member, key, mode);
  }
  assert.ok(MODE_SCOPED_ROWS.length > 100, 'the registry names the per-mode values');
});

test('every setting shared by both modes before PR-1 is set per drawing', () => {
  // The registry's form and advanced settings that both modes shared before
  // PR-1 (every row but a one-mode `own` row): each drawing now holds its own.
  const rows = MODE_SCOPED_ROWS.filter((row) => row.migrate !== 'own'
    && (row.domain === 'config.form' || row.domain === 'config.adv'));
  assert.ok(rows.length > 50, 'the registry names the shared settings');
  const sentinel = (value) => (typeof value === 'number' ? value + 7
    : typeof value === 'boolean' ? !value
      : typeof value === 'string' ? `${value}-circular`
        : Array.isArray(value) ? [...value, 'circular'] : value === null ? 'circular' : { circular: true });
  const checked = rows.flatMap((row) => {
    const container = row.domain === 'config.form' ? 'form' : 'adv';
    const circular = modeDrawing(state, 'circular')[container];
    const linear = modeDrawing(state, 'linear')[container];
    if (!Object.hasOwn(circular, row.path)) return [];
    const before = JSON.stringify(linear[row.path]);
    const saved = circular[row.path];
    circular[row.path] = sentinel(saved);
    const kept = JSON.stringify(linear[row.path]) === before;
    circular[row.path] = saved;
    assert.ok(kept, `${row.domain}.${row.path} reached the Linear drawing`);
    return [row.path];
  });
  assert.ok(checked.length > 50, `checked ${checked.length} shared settings`);
});

test('each mode request carries its own drawing values', () => {
  MODES.forEach(resetDrawing);
  const linearDefaults = samples('linear', requestOf('linear'));
  const circular = modeDrawing(state, 'circular');
  circular.adv.label_font_size = 13;
  circular.adv.axis_stroke_width = 3;
  circular.form.separate_strands = false;
  circular.adv.gc_content_mode = 'percent';
  circular.adv.depth_min = 4;
  assert.deepEqual(samples('circular', requestOf('circular')), [13, 3, false, 'percent', 4]);
  assert.deepEqual(samples('linear', requestOf('linear')), linearDefaults, 'the Linear request keeps the Linear values');
  const linear = modeDrawing(state, 'linear');
  linear.adv.label_font_size = 7;
  linear.adv.depth_min = 2;
  linear.adv.gc_content_mode = 'deviation';
  const fromLinear = samples('linear', requestOf('linear'));
  assert.deepEqual([fromLinear[0], fromLinear[3], fromLinear[4]], [7, 'deviation', 2]);
  assert.deepEqual(samples('circular', requestOf('circular')), [13, 3, false, 'percent', 4]);
  installShownDrawing('circular');
  MODES.forEach(resetDrawing);
});

test('two Circular Depth series and one Linear Depth file both build their requests (Depth series count)', () => {
  MODES.forEach(resetDrawing);
  // Two series configured while Circular had two Depth files.
  modeDrawing(state, 'circular').form.show_depth = true;
  modeDrawing(state, 'circular').adv.depth_tracks.splice(0, modeDrawing(state, 'circular').adv.depth_tracks.length,
    { label: 'Illumina', color: '#112233' }, { label: 'Nanopore', color: '#445566' });
  modeDrawing(state, 'linear').form.show_depth = true;
  installShownDrawing('linear');
  const linear = requestOf('linear', { linearDepth: [depthFile('a.depth.tsv')] });
  assert.equal(linear.renderRequest.diagramOptions.depthTracks?.length, 1);
  assert.equal(modeDrawing(state, 'circular').adv.depth_tracks.length, 2, 'the Linear request leaves the Circular series');
  installShownDrawing('circular');
  MODES.forEach(resetDrawing);
});

test('records turned OFF in one mode stay ON in the other (record selection)', () => {
  MODES.forEach(resetDrawing);
  const circular = modeDrawing(state, 'circular');
  const linear = modeDrawing(state, 'linear');
  installShownDrawing('circular');
  state.cInputType.value = 'gb';
  state.circularRecordList.value = [{ selector: '#1', record_id: 'a' }, { selector: '#2', record_id: 'b' }];
  circular.form.multi_record_canvas = true;
  const circularKeys = () => buildCanonicalRenderRequest({
    state, drawing: circular, filesData: { c_gb: genbank(), linearSeqs: [] }
  }).renderRequest.records.map((record) => record.recordKey);
  const savedCards = [...state.linearSeqs];
  state.linearSeqs.splice(0, Infinity, { uid: 'card-1' }, { uid: 'card-2' });
  linear.recordsOff.splice(0, Infinity, 'card-2');
  assert.deepEqual(circularKeys(), ['record-1', 'record-2'], 'a Linear OFF card leaves the Circular request alone');
  circular.recordsOff.splice(0, Infinity, '#2');
  assert.deepEqual(circularKeys(), ['record-1']);
  assert.deepEqual([...linear.recordsOff], ['card-2']);
  assert.deepEqual(linear.drawnLinearSeqs.value.map((card) => card.uid), ['card-1']);
  circular.recordsOff.splice(0);
  linear.recordsOff.splice(0);
  state.linearSeqs.splice(0, Infinity, ...savedCards);
  state.circularRecordList.value = [];
  MODES.forEach(resetDrawing);
});

// M6 isolation, node half (plan-r2 §4): every write a slice offers, its
// `install` and `reset` included, changes its own drawing and leaves the sha256
// of the other drawing's captures as it was. The table names each writer that
// a slice module exports, so a new writer fails until it is listed here; an
// export that only reads or builds a value, and each slice's `capture`, is
// listed in `SLICE_READERS` (shared with the R2 writer guard) and, given a
// handle, changes neither drawing. Each writer runs after its slice holds the
// seed, whose rows the writers address by identity (their own view's rows).
const SLICE_DIRECTORY = new URL('../../gbdraw/web/js/services/drawing/', import.meta.url);
const NOT_WRITERS = new Set(SLICE_READERS.filter((name) => !name.includes('.')));
const legendRow = (caption, originalCaption, color) => ({ caption, originalCaption, color, xPos: 0, yPos: 0, featureIds: [] });
const SEEDS = {
  canvas: () => ({ canvasPadding: { top: 1, right: 2, bottom: 3, left: 4 } }),
  visibilityRules: () => ({ featureVisibilityManualRules: [
    { id: 'rule-1', featureType: 'CDS', qualifier: 'locus_tag', value: '^A$', action: 'off' }
  ] }),
  colorRules: () => ({ manualSpecificRules: [{ feat: 'CDS', qual: 'product', val: 'kinase', color: '#123456', cap: 'Kinase' }] }),
  palette: () => ({
    selectedPalette: 'mint', currentColors: { CDS: '#123456' }, pendingPaletteName: 'queued', pendingPaletteColors: { CDS: '#abcdef' }
  }),
  labelTables: () => ({
    labelTextBulkOverrides: { 'hypothetical protein': 'HP' },
    canonicalLabelOverrideRows: [{ recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: '^kinase$', labelText: 'Kinase' }],
    manualPriorityRules: [{ feat: 'CDS', order: 'gene,product' }],
    filterMode: 'Whitelist',
    manualBlacklist: 'putative',
    manualWhitelist: [{ feat: 'CDS', qual: 'product', key: 'kinase' }]
  }),
  annotationSets: () => ({ annotationSets: [{ id: 'set-1', legendLabel: 'Set', annotations: [
    { id: 'span', target: { kind: 'coordinateSpan', start: 10, end: 20 }, label: 'Span' },
    { id: 'gene', target: { kind: 'feature', selectors: [{ key: 'gene', value: 'dnaA' }] }, label: 'Gene' }
  ] }] }),
  groupNames: () => ({
    orthogroupNameOverrides: { g1: 'Name' }, orthogroupDescriptionOverrides: { g1: 'Description' },
    orthogroupDormantOverrides: { '["member"]': { name: 'Dormant', description: 'Dormant description' } }
  }),
  featureEdits: () => ({
    featureOverrides: { '["record-1","CDS_1"]': {
      recordKey: 'record-1', biologicalFeatureId: 'CDS_1', featureVisibility: 'off', labelVisibility: null, labelText: 'A', labelSourceText: null
    } },
    featurePlacementOverrides: { '["record-1","CDS_2"]': {
      recordKey: 'record-1', biologicalFeatureId: 'CDS_2', placement: { kind: 'lane', side: 'outward', level: 1 }
    } }
  }),
  paint: () => ({
    legendEntries: [legendRow('Renamed CDS', 'CDS', '#abcdef'), legendRow('tRNA', 'tRNA', '#fedcba')],
    deletedLegendEntries: [legendRow('rRNA', 'rRNA', '#111111')],
    dormantLegendEntries: [{ caption: 'Coverage', originalCaption: 'depth', color: '#222222', featureIds: [] }],
    addedLegendCaptions: ['Kinase'],
    legendColorOverrides: { CDS: '#abcdef' },
    legendStrokeOverrides: { CDS: { strokeColor: '#333333', strokeWidth: 1.5 } },
    featureColorOverrides: { 'record-1:CDS_1': { color: '#444444', caption: 'Kinase' } },
    featureStrokeOverrides: { 'record-1:CDS_1': { strokeColor: '#555555', strokeWidth: 2 } }
  }),
  settings: (drawing) => {
    const settings = DRAWING_SLICES.settings.capture(drawing.slices.settings);
    return { ...settings, form: { ...settings.form, plot_title: 'Isolation title' }, losatProgram: 'tblastx' };
  }
};
const stroke = { strokeColor: '#666666', strokeWidth: 3 };
const firstSet = (drawing) => drawing.annotationSets[0];
const row = (drawing, index = 0) => firstSet(drawing).annotations[index];
// Each writer of each slice, as its callers call it on the drawing's handle.
const WRITES = {
  canvas: { setCanvasPaddingSide: (h) => modules.canvas.setCanvasPaddingSide(h, 'top', 9) },
  visibilityRules: { replaceVisibilityRules: (h) => modules.visibilityRules.replaceVisibilityRules(h, []) },
  colorRules: {
    replaceColorRules: (h) => modules.colorRules.replaceColorRules(h, [{ feat: 'tRNA', qual: 'product', val: 'x', color: '#000001', cap: '' }])
  },
  palette: {
    selectPalette: (h) => modules.palette.selectPalette(h, 'other', { CDS: '#000001' }),
    replaceDefaultColors: (h) => modules.palette.replaceDefaultColors(h, { tRNA: '#000002' }),
    setDefaultColors: (h) => modules.palette.setDefaultColors(h, { rRNA: '#000003' }),
    queuePalette: (h) => modules.palette.queuePalette(h, 'next', { CDS: '#000004' }),
    restoreArtifactQueuedPalette: (h) => modules.palette.restoreArtifactQueuedPalette(h, 'artifact', { CDS: '#000005' }),
    clearQueuedPalette: (h) => modules.palette.clearQueuedPalette(h)
  },
  labelTables: {
    replaceBulkLabelEdits: (h) => modules.labelTables.replaceBulkLabelEdits(h, { kinase: 'K' }),
    clearLabelTextEdits: (h) => modules.labelTables.clearLabelTextEdits(h),
    setLabelFilterMode: (h) => modules.labelTables.setLabelFilterMode(h, 'Blacklist'),
    setBlacklistText: (h) => modules.labelTables.setBlacklistText(h, 'unknown'),
    appendBlacklistWords: (h) => modules.labelTables.appendBlacklistWords(h, ['hypothetical']),
    addWhitelistRules: (h) => modules.labelTables.addWhitelistRules(h, [{ feat: 'tRNA', qual: 'product', key: 'trn' }]),
    removeWhitelistRule: (h) => modules.labelTables.removeWhitelistRule(h, 0),
    setWhitelistRuleField: (h, drawing) => modules.labelTables.setWhitelistRuleField(h, drawing.manualWhitelist[0], 'key', 'ligase'),
    setPriorityRules: (h) => modules.labelTables.setPriorityRules(h, [{ feat: 'tRNA', order: 'product' }]),
    removePriorityRule: (h) => modules.labelTables.removePriorityRule(h, 0)
  },
  annotationSets: {
    addAnnotationSet: (h) => modules.annotationSets.addAnnotationSet(h, createAnnotationSet({ id: 'set-2' })),
    removeAnnotationSet: (h, drawing) => modules.annotationSets.removeAnnotationSet(h, firstSet(drawing)),
    renameAnnotationSet: (h, drawing) => modules.annotationSets.renameAnnotationSet(h, firstSet(drawing), 'set-renamed'),
    setAnnotationSetLegendLabel: (h, drawing) => modules.annotationSets.setAnnotationSetLegendLabel(h, firstSet(drawing), 'Label'),
    addAnnotations: (h, drawing) => modules.annotationSets.addAnnotations(h, firstSet(drawing), [
      { id: 'added', target: { kind: 'coordinateSpan', start: 1, end: 2 }, label: 'Added' }
    ]),
    removeAnnotation: (h, drawing) => modules.annotationSets.removeAnnotation(h, firstSet(drawing), row(drawing)),
    removeAnnotations: (h, drawing) => modules.annotationSets.removeAnnotations(h, [row(drawing, 1)]),
    setAnnotationFields: (h, drawing) => modules.annotationSets.setAnnotationFields(h, firstSet(drawing), row(drawing), { label: 'Edited' }),
    setAnnotationStyleField: (h, drawing) => modules.annotationSets.setAnnotationStyleField(
      h, firstSet(drawing), row(drawing), 'color', '#000006', { color: '#000000' }
    ),
    replaceAnnotationTarget: (h, drawing) => modules.annotationSets.replaceAnnotationTarget(
      h, firstSet(drawing), row(drawing), { kind: 'coordinateSpan', start: 3, end: 4 }
    ),
    setAnnotationTargetBound: (h, drawing) => modules.annotationSets.setAnnotationTargetBound(h, firstSet(drawing), row(drawing), 'start', 5),
    setAnnotationSelectorField: (h, drawing) => modules.annotationSets.setAnnotationSelectorField(
      h, firstSet(drawing), row(drawing, 1), 'value', 'dnaN'
    ),
    bindAnnotationRecords: (h, drawing) => modules.annotationSets.bindAnnotationRecords(h, [
      { annotation: row(drawing), record: { kind: 'recordId', value: 'record-1' }, metadata: { bound: true } }
    ])
  },
  groupNames: {
    setGroupName: (h) => modules.groupNames.setGroupName(h, 'g2', 'Second'),
    setGroupDescription: (h) => modules.groupNames.setGroupDescription(h, 'g2', 'Second description'),
    clearDormantGroupNames: (h) => modules.groupNames.clearDormantGroupNames(h),
    clearGroupRename: (h) => modules.groupNames.clearGroupRename(h, 'g1')
  },
  featureEdits: {
    setFeatureOverrideRow: (h) => modules.featureEdits.setFeatureOverrideRow(h, '["record-1","CDS_1"]', null),
    clearFeatureOverrides: (h) => modules.featureEdits.clearFeatureOverrides(h),
    setFeaturePlacementRow: (h) => modules.featureEdits.setFeaturePlacementRow(h, '["record-1","CDS_2"]', null)
  },
  paint: {
    replaceLegendRows: (h, drawing) => modules.paint.replaceLegendRows(h, [drawing.legendEntries.value[1]]),
    deleteLegendRow: (h) => modules.paint.deleteLegendRow(h, 0),
    restoreLegendRows: (h, drawing) => modules.paint.restoreLegendRows(h, new Set([0]),
      [...drawing.legendEntries.value, drawing.deletedLegendEntries.value[0]]),
    replaceDormantLegendRows: (h) => modules.paint.replaceDormantLegendRows(h, []),
    replaceAddedLegendCaptions: (h) => modules.paint.replaceAddedLegendCaptions(h, ['Ligase']),
    renameListedLegendRow: (h, drawing) => modules.paint.renameListedLegendRow(h, drawing.legendEntries.value[1],
      { caption: 'Transfer RNA', originalCaption: 'tRNA' }),
    recolorLegendRows: (h) => modules.paint.recolorLegendRows(h, ['tRNA'], '#000007'),
    rekeyLegendPaint: (h) => modules.paint.rekeyLegendPaint(h,
      { deletedEntries: [], dormantEntries: [], colorOverrides: {}, strokeOverrides: {}, addedCaptions: [] }),
    retireLegendPaint: (h) => modules.paint.retireLegendPaint(h, new Set(['CDS'])),
    resetPaintEdits: (h) => modules.paint.resetPaintEdits(h),
    clearSessionPaint: (h) => modules.paint.clearSessionPaint(h),
    setLegendColor: (h) => modules.paint.setLegendColor(h, 'tRNA', '#000008'),
    setLegendStroke: (h) => modules.paint.setLegendStroke(h, 'tRNA', stroke),
    setFeatureFill: (h) => modules.paint.setFeatureFill(h, 'record-1:CDS_2', '#000009'),
    replaceFeatureFills: (h) => modules.paint.replaceFeatureFills(h, {}),
    setFeatureStroke: (h) => modules.paint.setFeatureStroke(h, 'record-1:CDS_2', stroke),
    clearStrokes: (h) => modules.paint.clearStrokes(h)
  },
  settings: {}
};
// The slice modules (a module that exports a slice), by the slice's name.
const sliceNameOf = new Map(Object.entries(DRAWING_SLICES).map(([name, slice]) => [slice, name]));
const modules = Object.fromEntries((await Promise.all(readdirSync(SLICE_DIRECTORY).filter((file) => file.endsWith('.js'))
  .map((file) => import(new URL(file, SLICE_DIRECTORY).href))))
  .flatMap((module) => Object.values(module).filter((value) => sliceNameOf.has(value))
    .map((slice) => [sliceNameOf.get(slice), module])));
const hashOf = (drawing) => createHash('sha256').update(JSON.stringify(Object.entries(DRAWING_SLICES)
  .map(([name, slice]) => slice.capture(drawing.slices[name])))).digest('hex');

test('every slice writer, install, and reset changes its own drawing only (M6 isolation, sha256 of the captures)', () => {
  assert.deepEqual(Object.keys(modules).sort(), Object.keys(DRAWING_SLICES).sort(), 'a module exports each slice');
  const listed = new Set();
  const captures = [];
  for (const [name, module] of Object.entries(modules)) {
    captures.push(`${Object.keys(module).find((key) => module[key] === DRAWING_SLICES[name])}.capture`);
    const functions = Object.entries(module).filter(([, value]) => typeof value === 'function').map(([key]) => key);
    assert.deepEqual(functions.filter((key) => !NOT_WRITERS.has(key)).sort(), Object.keys(WRITES[name]).sort(),
      `${name}: the table lists every writer the module exports`);
    functions.filter((key) => NOT_WRITERS.has(key)).forEach((key) => listed.add(key));
  }
  assert.deepEqual([...NOT_WRITERS].filter((key) => !listed.has(key)), [], 'each listed non-writer is a slice export');
  assert.deepEqual(SLICE_READERS.filter((key) => key.includes('.')).sort(), captures.sort(),
    'each slice\'s capture is listed under the name its module exports the slice by');
  const saved = Object.fromEntries(MODES.map((mode) => [mode, Object.fromEntries(Object.entries(DRAWING_SLICES)
    .map(([name, slice]) => [name, slice.capture(modeDrawing(state, mode).slices[name])]))]));
  for (const mode of MODES) {
    const drawing = modeDrawing(state, mode);
    const other = modeDrawing(state, MODES.find((candidate) => candidate !== mode));
    for (const [name, slice] of Object.entries(DRAWING_SLICES)) {
      const handle = drawing.slices[name];
      const seed = SEEDS[name](drawing);
      const writes = [
        ['install', (h) => slice.install(h, seed), (h) => slice.reset(h)],
        ['reset', (h) => slice.reset(h)],
        ['capture', (h) => slice.capture(h), undefined, 'reads'],
        ...Object.entries(WRITES[name]).map(([writer, write]) => [writer, (h) => write(h, drawing)]),
        ...[...listed].filter((key) => Object.hasOwn(modules[name], key))
          // A builder of a value may reject a handle; it still writes nothing.
          .map((key) => [key, (h) => { try { modules[name][key](h); } catch { /* not a handle reader */ } }, null, 'reads'])
      ];
      for (const [writer, write, prepare = (h) => slice.install(h, seed), reads] of writes) {
        if (prepare) prepare(handle);
        const own = hashOf(drawing);
        const others = hashOf(other);
        write(handle);
        if (reads) assert.equal(hashOf(drawing), own, `${mode} ${name}.${writer} only reads`);
        else assert.notEqual(hashOf(drawing), own, `${mode} ${name}.${writer} changes its drawing`);
        assert.equal(hashOf(other), others, `${mode} ${name}.${writer} leaves the other drawing's captures`);
      }
    }
  }
  for (const mode of MODES) {
    Object.entries(DRAWING_SLICES).forEach(([name, slice]) => slice.install(modeDrawing(state, mode).slices[name], saved[mode][name]));
  }
});
