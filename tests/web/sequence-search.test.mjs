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
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { createAnnotationSet } = await import('../../gbdraw/web/js/services/annotation-state.js');
const { createDefaultLinearTrackSlots } = await import('../../gbdraw/web/js/services/linear-track-slot-model.js');
const { RECORD_POPUP_FIELDS } = await import('../../gbdraw/web/js/app/preview-item-popup.js');
const { DEFINITION_LINE_STYLE_KINDS } = await import('../../gbdraw/web/js/services/definition-line-style-state.js');

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
// A drawing's Linear card holds the `InputFileRef` of its File (`state.inputFiles`);
// Generate and the live projectors read its bytes through `fileOf` (OV-2205).
const fileRef = (file) => state.inputFiles.registerFile(file);
const fileOf = state.inputFiles.fileOf;
const seq = (uid, id) => ({
  uid, gb: fileRef(genbank([id])), depth: null, losat_gencode: 1, region_record_id: '', region_start: null, region_end: null,
  region_reverse: false, definition: '', record_subtitle: '', file_definition: `File ${id}`, file_subtitle: ''
});
const inputs = { linearSeqs: [seq('record-1', 'LIVE1'), seq('record-2', 'LIVE2')], fileOf };

// The committed Session a Generate of `mode` builds from that mode's drawing now.
const generate = (mode) => {
  const drawing = modeDrawing(state, mode);
  config.installShownDrawing(mode);
  modeDrawing(state, 'circular').inputType.value = 'gb';
  modeDrawing(state, 'linear').inputType.value = 'gb';
  const filesData = mode === 'circular'
    ? { c_gb: genbank(['LIVE1']), linearSeqs: [] }
    : { linearSeqs: inputs.linearSeqs.map((card) => ({ ...card, gb: fileOf(card.gb) })), linearComparisons: [] };
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


const search = await import('../../gbdraw/web/js/services/sequence-search.js');
const { resultRecordBindings } = await import('../../gbdraw/web/js/app/annotations/record-selector.js');
const ANNOTATIONS = 'diagramOptions.annotations.sets';

// The set a search adds, as the Python helper writes it (`encode_annotation_set`) and the panel reads it.
const searchSet = (id) => createAnnotationSet({
  id,
  defaultStyle: { fill: '#c8607e', stroke: '#c8607e', fillOpacity: 0.85 },
  annotations: [{
    id: 'q1#1', label: '', mark: 'highlight', lane: null, legendLabel: 'q1', style: null, strand: '-',
    target: { kind: 'coordinateSpan', record: { kind: 'recordId', value: 'LIVE1' }, start: 5, end: 20, coordinateSpace: 'source' },
    metadata: { query: 'q1', program: 'blastn', identity: 100, alignment_length: 16, evalue: 1e-5, bitscore: 30 }
  }]
});
// `recordSelector` names a bound row's record in the committed request (the shown Result's, by binding key).
const projectSearch = (committed, mode, ids, recordSelector = () => null) => request.projectLiveSettingWrites({
  committed, drawing: modeDrawing(state, mode), requestPaths: [ANNOTATIONS], inputs, searchSets: { ids, recordSelector }
});

for (const mode of ['circular', 'linear']) {
  test(`${mode}: a search's set is drawn at once as Generate sends it, and the other pending sets wait`, () => {
    const drawing = modeDrawing(state, mode);
    const before = [...drawing.annotationSets];
    try {
      const committed = generate(mode);
      assert.equal(committed.renderRequest.diagramOptions.annotations, undefined);
      // A set added by hand applies on Generate; the search's set by its id applies at once.
      drawing.annotationSets.push(createAnnotationSet({ id: 'by-hand' }));
      drawing.annotationSets.push(searchSet('Search 1'));
      assert.deepEqual(projectSearch(committed, mode, []), []);
      const live = applied(committed, projectSearch(committed, mode, ['Search 1']));
      drawing.annotationSets.splice(drawing.annotationSets.findIndex((set) => set.id === 'by-hand'), 1);
      const draft = generate(mode);
      assert.deepEqual(live.renderRequest.diagramOptions.annotations, draft.renderRequest.diagramOptions.annotations);
      // Each hit keeps the strand Python wrote, which places it on its strand's lane (D-49, D-50).
      assert.equal(live.renderRequest.diagramOptions.annotations.sets[0].annotations[0].strand, '-');
      // Its Undo removes it at once: the request draws no annotation set again.
      drawing.annotationSets.pop();
      const undone = applied(live, projectSearch(live, mode, ['Search 1']));
      assert.equal(undone.renderRequest.diagramOptions.annotations, undefined);
      assert.deepEqual(undone.renderRequest, committed.renderRequest);
    } finally {
      drawing.annotationSets.splice(0, drawing.annotationSets.length, ...before);
    }
  });
}

test('a History step names the sets a search added: its own step, Undo and Redo; not a set added or renamed by hand', () => {
  const set = JSON.parse(JSON.stringify(searchSet('Search 1')));
  const searched = new Map([['Search 1', JSON.stringify(set)]]);
  const path = ['config', 'annotationSets'];
  assert.deepEqual(search.searchSetIdsOfStep([{ path, before: [], after: [set] }], searched), ['Search 1']);
  assert.deepEqual(search.searchSetIdsOfStep([{ path, before: [set], after: [] }], searched), []);
  assert.deepEqual(search.searchSetIdsOfStep([{ path, before: [], after: [{ ...set, legendLabel: 'x' }] }], searched), []);
  assert.deepEqual(search.searchSetIdsOfStep([{ path, before: [], after: [{ ...set, id: 'Search 2' }] }], searched), []);
  assert.deepEqual(search.searchSetIdsOfStep([{ path: ['config', 'form', 'species'], before: 'a', after: 'b' }], searched), []);
  assert.equal(search.nextSearchSetId(['Search 1', 'Search 3']), 'Search 2');
});

test('the search reads each drawn record once over the bases it draws, and names its record as the Result draws it', () => {
  const axis = (drawnLength, recordBase = 1, step = 1) => ({
    kind: 'linear', rows: [], source: { anchorBase: recordBase, recordBase, step, sourceLength: 100, drawnLength, isCircular: false }
  });
  const records = [
    { resultIndex: 0, recordIndex: 0, recordKey: 'a', recordId: 'X', axis: axis(100) },
    { resultIndex: 1, recordIndex: 0, recordKey: 'a', recordId: 'X', axis: axis(100) },
    { resultIndex: 0, recordIndex: 1, recordKey: 'b', recordId: 'X', axis: axis(30, 41) }
  ];
  assert.deepEqual(search.searchSpans(records, 600), [{ recordKey: 'a', start: 1, end: 100 }, { recordKey: 'b', start: 41, end: 70 }]);
  const read = { fasta: '', spans: [
    { recordKey: 'a', recordId: 'X', start: 1, end: 100, length: 100, wrapsOrigin: false, originWindow: false, recordLength: 100, header: 'X:1-100', sequence: '' },
    { recordKey: 'b', recordId: 'X', start: 41, end: 70, length: 30, wrapsOrigin: false, originWindow: false, recordLength: 100, header: 'X:41-70', sequence: '' }
  ] };
  // Two drawn records share the ID X: each is named by its `#<n>` among the Result's records (`resultRecordBindings`).
  const bindings = resultRecordBindings(records, []);
  assert.deepEqual(search.searchSubjects(read, bindings), [
    { subjectId: 'X:1-100', record: { kind: 'recordIndex', index: 0 }, recordKey: 'a', sourceStart: 1, originWindow: false, recordLength: 100 },
    { subjectId: 'X:41-70', record: { kind: 'recordIndex', index: 1 }, recordKey: 'b', sourceStart: 41, originWindow: false, recordLength: 100 }
  ]);
  assert.throws(() => search.searchSubjects(read, resultRecordBindings([], [])), (error) => error.code === 'RESULT_INVALID',
    'a span of a record the Result does not draw does not match the read');
});

// OV-2571, D-52: a record that a span may cross the origin of (circular, drawn whole), in either mode, is also read
// over its origin window, which Python places (origin_window_span); both subjects name the record by its catalog
// record, and the window is marked as one.
test('the search also reads the origin window of each record drawn whole across its origin', () => {
  const axis = (kind, drawnLength, isCircular) => ({
    kind, rows: [], source: { anchorBase: 1, recordBase: 1, step: 1, sourceLength: 100, drawnLength, isCircular }
  });
  const records = [
    { resultIndex: 0, recordIndex: 0, recordKey: 'c', recordId: 'C', axis: axis('circular', 100, true) },
    { resultIndex: 0, recordIndex: 1, recordKey: 'l', recordId: 'L', axis: axis('linear', 100, true) },
    { resultIndex: 0, recordIndex: 2, recordKey: 'n', recordId: 'N', axis: axis('circular', 100, false) }
  ];
  assert.deepEqual(search.searchSpans(records, 30), [
    { recordKey: 'c', start: 1, end: 100 }, { recordKey: 'l', start: 1, end: 100 }, { recordKey: 'n', start: 1, end: 100 },
    { recordKey: 'c', originWindow: 30 }, { recordKey: 'l', originWindow: 30 }
  ]);
  const read = { fasta: '', spans: [
    { recordKey: 'c', recordId: 'C', start: 1, end: 100, length: 100, wrapsOrigin: false, originWindow: false, recordLength: 100, header: 'C:1-100', sequence: '' },
    { recordKey: 'c', recordId: 'C', start: 71, end: 30, length: 60, wrapsOrigin: true, originWindow: true, recordLength: 100, header: 'C:71-30', sequence: '' }
  ] };
  assert.deepEqual(search.searchSubjects(read, resultRecordBindings(records, [])), [
    { subjectId: 'C:1-100', record: { kind: 'recordId', value: 'C' }, recordKey: 'c', sourceStart: 1, originWindow: false, recordLength: 100 },
    { subjectId: 'C:71-30', record: { kind: 'recordId', value: 'C' }, recordKey: 'c', sourceStart: 71, originWindow: true, recordLength: 100 }
  ]);
});

test('the helper replies are checked against their typed shapes', () => {
  assert.deepEqual(search.readSearchQueriesReply({ queries: [{ name: 'q', length: 4, fasta: '>q\nACGT\n' }] }),
    [{ name: 'q', length: 4, fasta: '>q\nACGT\n' }]);
  assert.throws(() => search.readSearchQueriesReply({ queries: [{ name: 'q' }] }), (error) => error.code === 'RESULT_INVALID');
  assert.throws(() => search.readSequenceSearchReply({ annotationSet: null, hits: [{}] }), (error) => error.code === 'RESULT_INVALID');
  assert.deepEqual(search.panelQueries([{ name: 'q', length: 4, fasta: '' }], [{ name: 'q', length: 4, fasta: '', color: '#000000' }])[0].color, '#000000');
});

// PD-OI-117 (Q-B10 A): with a custom track stack, or the Circular track sizes (an explicit stack too), the search's
// set is sent with its marker and the stack unchanged; Python places it, so the live request equals Generate's.
const STACKS = [
  ['circular', 'custom stack', (drawing) => set(drawing.adv, 'circular_track_slots_enabled', true)],
  ['circular', 'track sizes', (drawing) => set(drawing.adv, 'feature_width_circular', 30)],
  ['linear', 'custom stack', (drawing) => set(drawing.adv, 'linear_track_slots_enabled', true)]
];
for (const [mode, label, explicit] of STACKS) {
  test(`${mode}, ${label}: a search's set is drawn at once with the stack unchanged`, () => {
    const drawing = modeDrawing(state, mode);
    const before = [...drawing.annotationSets];
    try {
      explicit(drawing);
      const committed = generate(mode);
      const tracks = structuredClone(committed.renderRequest.diagramOptions.tracks);
      assert.ok(tracks.circularTrackSlots || tracks.linearTrackSlots, 'an explicit stack');
      drawing.annotationSets.push(createAnnotationSet({ ...searchSet('Search 1'), keepHighlightOverlay: true }));
      const live = applied(committed, projectSearch(committed, mode, ['Search 1']));
      const draft = generate(mode);
      assert.deepEqual(live.renderRequest, draft.renderRequest);
      assert.deepEqual(draft.renderRequest.diagramOptions.tracks, tracks);
      assert.equal(draft.renderRequest.diagramOptions.annotations.sets[0].keepHighlightOverlay, true);
    } finally {
      drawing.annotationSets.splice(0, drawing.annotationSets.length, ...before);
      undoAll();
    }
  });
}

// Receipt web.sequence-search: a search's set keeps Region Annotations' timing, where colours apply at once; its
// rows' Legend labels (one row per query, or "Query position") do not hold them back for Generate, as they do for
// another set's row (O-8). Every edit of a search row's colours, the gradient's two colours included, draws live as
// Generate sends it.
for (const mode of ['circular', 'linear']) {
  test(`${mode}: the colours of a search's rows apply at once though they have a Legend row`, () => {
    const drawing = modeDrawing(state, mode);
    const before = [...drawing.annotationSets];
    try {
      const gradient = { startColor: '#f6c9d8', endColor: '#5b1a4a', direction: 'reverse', queryStart: 0.25, queryEnd: 0.5 };
      const searched = createAnnotationSet({ ...searchSet('Search 1'), keepHighlightOverlay: true });
      const labelled = createAnnotationSet({ ...searchSet('by-hand'), annotations: searchSet('by-hand').annotations.map((item) => ({ ...item, strand: undefined })) });
      searched.annotations[0].style = createAnnotationSet({ defaultStyle: { fill: '#5b1a4a', stroke: '#5b1a4a', gradient } }).defaultStyle;
      drawing.annotationSets.push(searched, labelled);
      const committed = generate(mode);
      const items = { search: drawing.annotationSets.at(-2).annotations[0], other: drawing.annotationSets.at(-1).annotations[0] };
      const style = (item, fill, gradientColors) => ({ ...(item.style ?? createAnnotationSet({}).defaultStyle), fill, stroke: fill,
        ...(gradientColors ? { gradient: { ...gradient, ...gradientColors } } : {}) });
      const note = (setId, itemId) => request.liveSettingNoteApplies(
        committed.renderRequest, ANNOTATIONS, drawing, inputs, committed.resources, { setId, itemId });
      assert.equal(note('Search 1', 'q1#1'), true);
      assert.equal(note('by-hand', 'q1#1'), false, 'a Legend-labelled row of another set waits (O-8)');
      set(items.search, 'style', style(items.search, '#123456', { startColor: '#ffffff', endColor: '#000000' }));
      set(items.other, 'style', style(items.other, '#654321'));
      const live = applied(committed, project(committed, mode, [ANNOTATIONS]));
      const [liveSearch, liveOther] = live.renderRequest.diagramOptions.annotations.sets;
      const draftSets = generate(mode).renderRequest.diagramOptions.annotations.sets;
      assert.deepEqual(liveSearch, draftSets[0]);
      assert.equal(liveSearch.annotations[0].style.gradient.startColor, '#ffffff');
      assert.equal(liveSearch.annotations[0].style.gradient.endColor, '#000000');
      assert.deepEqual(liveOther, committed.renderRequest.diagramOptions.annotations.sets[1], 'the other set waits');
    } finally {
      drawing.annotationSets.splice(0, drawing.annotationSets.length, ...before);
      undoAll();
    }
  });
}

// Receipt web.sequence-search: a row click outlines that hit; a chip outlines its query's hits; All outlines none.
// The outline names hits by the annotation id Python writes on the drawn hit.
test('the outline names the clicked row\'s hit, else the chosen query\'s hits, by annotation id', () => {
  const hit = (annotationId, query) => ({ annotationId, query });
  const hits = [hit('a#1', 'a'), hit('b#1', 'b'), hit('a#2', 'a')];
  assert.deepEqual(search.outlinedHitIds(hits, '', ''), []);
  assert.deepEqual(search.outlinedHitIds(hits, 'a', ''), ['a#1', 'a#2']);
  assert.deepEqual(search.outlinedHitIds(hits, 'a', 'a#2'), ['a#2']);
  assert.deepEqual(search.outlinedHitIds(hits, '', 'b#1'), ['b#1']);
  assert.deepEqual(search.outlinedHitIds(hits, '', 'gone#1'), [], 'a row no longer listed outlines nothing');
});

// PD-OI-115 Use as search query (PD-OI-117): the span's bases, read by Python and named by their header, join the
// query text, which Python reads again; the folded panel opens. The same span is not added twice, and a removed
// query stays removed when another span is added. The folded line names the program once.
test('Use as search query adds a span as a query once; a removed query stays removed', async () => {
  const { createSequenceSearch } = await import('../../gbdraw/web/js/app/sequence-search.js');
  /** The query reader's reply for a FASTA text, as Python's `read_search_queries` reads it (one query per header). */
  const readFasta = (text) => text.split('>').filter(Boolean).map((entry) => {
    const [header, ...lines] = entry.split('\n');
    return { name: header.trim(), length: lines.join('').length, fasta: `>${entry}` };
  });
  const reads = [];
  const panel = createSequenceSearch({
    state: { activeDrawing: () => ({ id: 'circular', mode: 'circular', annotationSets: [] }) },
    readAxisRecords: () => [],
    readRecordSpans: async () => { throw new Error('no search here'); },
    runLosatJobs: async () => ({ status: 'canceled' }),
    runHelper: async (_operation, { files }) => {
      const text = new TextDecoder().decode(files[0].bytes);
      reads.push(text);
      return { result: { queries: readFasta(text) } };
    },
    helperOperations: { readQueries: 'readSearchQueries', buildHitSet: 'searchHitAnnotationSet' },
    addSearchSet: () => {},
    shownRoot: () => null,
    closeRegionPopup: () => {},
    revealSection: async () => {}
  });
  await panel.readQueries('>q1\nACGT\n');
  panel.setFolded(true);
  await panel.addQuery({ name: 'NC_1:5-8', fasta: '>NC_1:5-8\nTTGA\n' });
  assert.equal(panel.panel.folded, false);
  assert.deepEqual(panel.panel.queries.map((query) => query.name), ['q1', 'NC_1:5-8']);
  assert.equal(reads.at(-1), '>q1\nACGT\n>NC_1:5-8\nTTGA\n');
  await panel.addQuery({ name: 'NC_1:5-8', fasta: '>NC_1:5-8\nTTGA\n' });
  assert.equal(reads.length, 2, 'the same span is not read or added twice');
  panel.removeQuery('q1');
  await panel.addQuery({ name: 'NC_1:9-12', fasta: '>NC_1:9-12\nCCCC\n' });
  assert.deepEqual(panel.panel.queries.map((query) => query.name), ['NC_1:5-8', 'NC_1:9-12']);
  assert.equal(panel.programText.value, 'blastn, task megablast');
  panel.panel.task = 'blastn';
  assert.equal(panel.programText.value, 'blastn');
});

// OV-2579: a click on Search blurs the query box, whose change event starts the read of the pasted text before the
// click lands. While that read runs, Search stays on (the search waits for the read), so the click searches.
test('Search stays on while a pasted text is read, so a click right after the paste searches', async () => {
  const { createSequenceSearch } = await import('../../gbdraw/web/js/app/sequence-search.js');
  /** @type {() => void} */
  let finish = () => {};
  const read = new Promise((resolve) => { finish = resolve; });
  const panel = createSequenceSearch({
    state: { activeDrawing: () => ({ id: 'circular', mode: 'circular', annotationSets: [] }) },
    readAxisRecords: () => [{ recordKey: 'record-1' }],
    readRecordSpans: async () => { throw new Error('no search here'); },
    runLosatJobs: async () => ({ status: 'canceled' }),
    runHelper: async () => {
      await read;
      return { result: { queries: [{ name: 'q', length: 4, fasta: '>q\nACGT\n' }] } };
    },
    helperOperations: { readQueries: 'readSearchQueries', buildHitSet: 'searchHitAnnotationSet' },
    addSearchSet: () => {},
    shownRoot: () => null,
    closeRegionPopup: () => {},
    revealSection: async () => {}
  });
  const reading = panel.readQueries('>q\nACGT\n');
  assert.equal(panel.searchDisabledReason.value, '', 'on while the text is read');
  finish();
  await reading;
  assert.equal(panel.searchDisabledReason.value, '');
  assert.deepEqual(panel.panel.queries.map(({ name }) => name), ['q']);
});

// A search through the panel with stand-in ports: a Result of two records that share the ID X (cards `card-a` and
// `card-b`, whose record settings owners bind them by `card-a:1` and `card-b:1`), one query, and the reply Python gives
// for one hit on `card-b`. `ports` replaces any port; given as a function, it reads the panel (`ref.panel`).
const { ANNOTATION_RECORD_BINDING_KEY, annotationOptionsPayload } = await import('../../gbdraw/web/js/services/annotation-state.js');
/** @param {string} uid @param {string[]} bindingKeys */
const owner = (uid, bindingKeys) => ({ key: uid, requestKeys: [uid], ownsExpansions: true, bindingKeys, displaySource: uid, displaySelector: null });
const runPanelSearch = async (ports = {}) => {
  const { createSequenceSearch } = await import('../../gbdraw/web/js/app/sequence-search.js');
  const drawing = { id: 'linear', mode: 'linear', annotationSets: [] };
  const added = [];
  const sent = [];
  const axis = { kind: 'linear', rows: [], source: { anchorBase: 1, recordBase: 1, step: 1, sourceLength: 100, drawnLength: 100, isCircular: false } };
  const span = { recordKey: 'card-b', recordId: 'X', start: 1, end: 100, length: 100, wrapsOrigin: false, originWindow: false, recordLength: 100, header: 'X:1-100', sequence: 'ACGT' };
  const reply = (input) => ({
    annotationSet: {
      id: input.setId, defaultStyle: { fill: '#c8607e', stroke: '#c8607e', fillOpacity: 0.85 }, keepHighlightOverlay: true,
      annotations: [{
        id: 'q#1', label: '', mark: 'highlight', lane: null, legendLabel: 'q', style: null, strand: '+',
        target: { kind: 'coordinateSpan', record: input.subjects[0].record, start: 5, end: 20, coordinateSpace: 'source' },
        metadata: { query: 'q', program: 'blastn', identity: '100', alignment_length: '16', evalue: '1e-05', bitscore: '30' }
      }]
    },
    hits: [{ annotationId: 'q#1', query: 'q', recordKey: 'card-b', start: 5, end: 20, strand: '+', identity: 100, alignmentLength: 16, evalue: 1e-5, bitscore: 30, qstart: 1, qend: 16 }],
    rowsRead: 1, rowsKept: 1, thresholds: [], unmatchedRows: 0, unmatchedSubjectIds: []
  });
  const ref = {};
  const axisRecords = [
    { resultIndex: 0, recordIndex: 0, recordKey: 'card-a', recordId: 'X', axis },
    { resultIndex: 0, recordIndex: 1, recordKey: 'card-b', recordId: 'X', axis }
  ];
  const panel = createSequenceSearch({
    state: { activeDrawing: () => drawing },
    readAxisRecords: () => axisRecords,
    readRecordSpans: async () => ({ fasta: '>X:1-100\nACGT\n', spans: [span] }),
    committedSets: () => [],
    readRecordBindings: () => resultRecordBindings(axisRecords, [owner('card-a', ['card-a:1']), owner('card-b', ['card-b:1'])]),
    runLosatJobs: async () => ({ status: 'ok', results: [{ text: 'q\tX:1-100\t100\t16\t0\t0\t1\t16\t5\t20\t1e-5\t30\n' }] }),
    runHelper: async (operation, payload) => {
      if (operation === 'readSearchQueries') return { result: { queries: [{ name: 'q', length: 16, fasta: '>q\nACGT\n' }] } };
      sent.push(payload.search);
      return { result: reply(payload.search) };
    },
    helperOperations: { readQueries: 'readSearchQueries', buildHitSet: 'searchHitAnnotationSet' },
    addSearchSet: (_drawingId, set) => { added.push(set); drawing.annotationSets.push(set); },
    shownRoot: () => null,
    closeRegionPopup: () => {},
    revealSection: async () => {},
    ...(typeof ports === 'function' ? ports(ref) : ports)
  });
  ref.panel = panel;
  await panel.readQueries('>q\nACGT\n');
  await panel.search();
  return { panel, drawing, added, sent };
};

// Review N1 and r2 B1 (#1667): a search's rows bind their record as Add to annotation set binds a row, by the key of
// its record settings owner (so a row of a record turned OFF waits with it at Generate), and name it as the shown
// Result draws it (`#2` of two records named X), never through the draft's record catalog.
test('a search binds each hit by its record settings owner and names its record as the Result draws it', async () => {
  const { panel, added, sent } = await runPanelSearch();
  assert.equal(panel.panel.searchError, null);
  const second = { kind: 'recordIndex', index: 1 };
  assert.deepEqual(sent[0].subjects.map(({ record }) => record), [second]);
  const [row] = added[0].annotations;
  assert.equal(row.metadata[ANNOTATION_RECORD_BINDING_KEY], 'card-b:1');
  assert.deepEqual(row.target.record, second);
  assert.deepEqual(annotationOptionsPayload(added, [], ['card-b:1']).sets[0].annotations, []);
  assert.equal(annotationOptionsPayload(added, [], []).sets[0].annotations.length, 1);
});

// Review r2 B1: the binding owner of a Result's records. A record binds by its record settings owner's key (a card
// that draws every record of its file binds its k-th record by its k-th key; an OFF card's owner binds too; a record
// no owner holds binds by none) and is named among the Result's records as the catalog names the draft's: by its
// record ID when no other drawn record has it, else `#<n>`.
test('a Result record binds by its settings owner and is named among the Result records', () => {
  const owners = [owner('u1', ['f1::a']), owner('u2', ['f2::a']), owner('u3', ['f3::a', 'f3::b'])];
  const bindings = resultRecordBindings([
    { resultIndex: 0, recordIndex: 0, recordKey: 'u1', recordId: 'contig1' },
    { resultIndex: 0, recordIndex: 1, recordKey: 'u2', recordId: 'contig1' },
    { resultIndex: 0, recordIndex: 2, recordKey: 'u3:1', recordId: 'NC_1' },
    { resultIndex: 1, recordIndex: 2, recordKey: 'u3:1', recordId: 'NC_1' },
    { resultIndex: 0, recordIndex: 3, recordKey: 'u3:2', recordId: 'NC_2' },
    { resultIndex: 0, recordIndex: 4, recordKey: 'gone', recordId: 'NC_9' }
  ], owners);
  assert.deepEqual(['u1', 'u2', 'u3:1', 'u3:2', 'gone'].map(bindings.bindingKeyOf), ['f1::a', 'f2::a', 'f3::a', 'f3::b', '']);
  assert.deepEqual(['u1', 'u2', 'u3:2', 'other'].map(bindings.selectorOf), [
    { kind: 'recordIndex', index: 0 }, { kind: 'recordIndex', index: 1 }, { kind: 'recordId', value: 'NC_2' }, null
  ]);
  assert.deepEqual(['f2::a', 'f3::a', '', 'f9::a'].map(bindings.selectorOfBinding), [
    { kind: 'recordIndex', index: 1 }, { kind: 'recordId', value: 'NC_1' }, null, null
  ]);
});

// Review r2 B1: the live draw names a search row's record among the committed request's records by its binding key,
// not by the selector the draft holds (it follows the draft's records, which may differ from the Result's: a record
// OFF, the Record choice changed, a Redo after a Generate). A row of a record that request does not draw is not
// drawn, as Generate draws no row of an OFF record; the live request then equals Generate's.
test('linear: the live draw names the record of a search row by its binding key among the committed records', () => {
  const drawing = modeDrawing(state, 'linear');
  const before = [...drawing.annotationSets];
  try {
    const committed = generate('linear');
    const set = searchSet('Search 1');
    const [row] = set.annotations;
    // The draft names LIVE1 (the draft's own selector); the row's key is LIVE2's in the committed request.
    set.annotations = [
      { ...row, metadata: { ...row.metadata, [ANNOTATION_RECORD_BINDING_KEY]: 'card-2' } },
      { ...row, id: 'q1#2', metadata: { ...row.metadata, [ANNOTATION_RECORD_BINDING_KEY]: 'card-off' } }
    ];
    drawing.annotationSets.push(set);
    const recordSelector = (key) => (key === 'card-2' ? { kind: 'recordId', value: 'LIVE2' } : null);
    const live = applied(committed, projectSearch(committed, 'linear', ['Search 1'], recordSelector));
    const [drawn] = live.renderRequest.diagramOptions.annotations.sets[0].annotations;
    assert.deepEqual(live.renderRequest.diagramOptions.annotations.sets[0].annotations.map(({ id }) => id), ['q1#1']);
    assert.deepEqual(drawn.target.record, { kind: 'recordId', value: 'LIVE2' });
    // Generate of the draft once the record select resolved the key (`reconcileAnnotationRecordBindings`) and with the
    // OFF row left out.
    set.annotations = [{ ...set.annotations[0], target: { ...row.target, record: { kind: 'recordId', value: 'LIVE2' } } }];
    assert.deepEqual(live.renderRequest.diagramOptions.annotations, generate('linear').renderRequest.diagramOptions.annotations);
  } finally {
    drawing.annotationSets.splice(0, drawing.annotationSets.length, ...before);
  }
});

// Review B2 (#1667): a set id the draft freed (its set renamed or deleted, which waits for Generate) is not reused
// while the committed Session still draws it, so the first search's hits stay drawn until Generate and the second
// search's Undo removes only its own hits.
for (const mode of ['circular', 'linear']) {
  test(`${mode}: a new search set id is unique against the draft and the committed Session`, async () => {
    const { added } = await runPanelSearch({ committedSets: () => [{ id: 'Search 1' }] });
    assert.equal(added[0].id, 'Search 2');
    const drawing = modeDrawing(state, mode);
    const before = [...drawing.annotationSets];
    try {
      drawing.annotationSets.push(searchSet('Search 1'));
      const live = applied(generate(mode), projectSearch(generate(mode), mode, ['Search 1']));
      // Renamed in Region Annotations: waits for Generate. The next search adds "Search 2".
      drawing.annotationSets[drawing.annotationSets.length - 1].id = 'COX1 hits';
      drawing.annotationSets.push(searchSet('Search 2'));
      const second = applied(live, projectSearch(live, mode, ['Search 2']));
      assert.deepEqual(second.renderRequest.diagramOptions.annotations.sets.map(({ id }) => id), ['Search 1', 'Search 2']);
      drawing.annotationSets.pop();
      const undone = applied(second, projectSearch(second, mode, ['Search 2']));
      assert.deepEqual(undone.renderRequest.diagramOptions.annotations, live.renderRequest.diagramOptions.annotations);
    } finally {
      drawing.annotationSets.splice(0, drawing.annotationSets.length, ...before);
    }
  });
}

// Receipt: a search that fails adds nothing and keeps the drawing; the panel says why.
test('a search whose helper fails adds nothing and shows the reason', async () => {
  const { panel, added } = await runPanelSearch({
    runHelper: async (operation) => {
      if (operation === 'readSearchQueries') return { result: { queries: [{ name: 'q', length: 16, fasta: '>q\nACGT\n' }] } };
      throw new Error('The search subject ID names more than one searched sequence.');
    }
  });
  assert.deepEqual(added, []);
  assert.equal(panel.panel.outcome, null);
  assert.notEqual(panel.panel.searchError, null);
});

// Review B1 (#1667), receipt: the panel's queries and settings are not saved; after a Load the panel is empty, and a
// History step of the loaded document is no search step.
test('a Session Load empties the panel: queries, settings, hit table, chips and outline', async () => {
  const { panel, added } = await runPanelSearch();
  panel.panel.task = 'blastn';
  panel.panel.colorMode = 'gradient';
  panel.chooseChip('q');
  // Outside Vue the computed values are refs.
  panel.selectHit(panel.hitRows.value[0]);
  assert.notEqual(panel.outline.value, null);
  const step = [{ path: ['config', 'annotationSets'], before: [], after: [JSON.parse(JSON.stringify(added[0]))] }];
  assert.deepEqual(panel.searchSetIdsOf(step), ['Search 1']);
  panel.reset();
  const { queryText, queries, outcome, folded, chip, selectedHit, task, colorMode, searchError, queryError } = panel.panel;
  assert.deepEqual({ queryText, queries, outcome, folded, chip, selectedHit, task, colorMode, searchError, queryError }, {
    queryText: '', queries: [], outcome: null, folded: false, chip: '', selectedHit: '', task: 'megablast', colorMode: 'query',
    searchError: null, queryError: null
  });
  assert.deepEqual(panel.hitRows.value, []);
  assert.equal(panel.outline.value, null);
  assert.deepEqual(panel.searchSetIdsOf(step), []);
});

// Review r2 N-r2-1, receipt "after a Load the panel is empty": a Session Load counts from its start (`loadStarts`,
// before it installs the loaded document); a search it overtook adds nothing, and its failure or busy reason is not
// written into the panel the Load empties (`reset`, at the Load's end).
test('a search a Session Load overtakes adds nothing and leaves the emptied panel empty', async () => {
  const replied = await runPanelSearch((ref) => ({
    runHelper: async (operation) => {
      if (operation === 'readSearchQueries') return { result: { queries: [{ name: 'q', length: 16, fasta: '>q\nACGT\n' }] } };
      ref.panel.loadStarts();
      return { result: { annotationSet: null, hits: [], rowsRead: 0, rowsKept: 0, thresholds: [], unmatchedRows: 0, unmatchedSubjectIds: [] } };
    }
  }));
  assert.deepEqual(replied.added, []);
  assert.equal(replied.panel.panel.outcome, null, 'the reply of the previous document is not shown');
  const failed = await runPanelSearch((ref) => ({
    runHelper: async (operation) => {
      if (operation === 'readSearchQueries') return { result: { queries: [{ name: 'q', length: 16, fasta: '>q\nACGT\n' }] } };
      ref.panel.loadStarts();
      ref.panel.reset();
      throw new Error('The search subject ID names more than one searched sequence.');
    }
  }));
  assert.equal(failed.panel.panel.searchError, null);
  const busy = await runPanelSearch((ref) => ({
    runLosatJobs: async () => {
      ref.panel.loadStarts();
      ref.panel.reset();
      return { status: 'busy', reason: 'A Session load is in progress.' };
    }
  }));
  assert.equal(busy.panel.panel.busy, '');
  assert.equal(busy.panel.panel.searchError, null);
});

// Review r2 N-r2-2 (N5): the hit set reply is checked in full, its threshold tallies and unmatched subject IDs too.
test('a hit set reply with a malformed threshold tally or unmatched subject ID is refused', () => {
  const good = {
    annotationSet: null, hits: [], rowsRead: 2, rowsKept: 0, unmatchedRows: 0, unmatchedSubjectIds: ['X:1-9'],
    thresholds: [{ name: 'E-value', threshold: 0.01, rowsLeftOut: 2, bestValue: 0.5 }]
  };
  assert.deepEqual(search.readSequenceSearchReply(good), good);
  for (const bad of [
    { thresholds: [{ name: 'E-value' }] }, { thresholds: [{ ...good.thresholds[0], rowsLeftOut: '2' }] }, { thresholds: undefined },
    { unmatchedSubjectIds: [1] }, { unmatchedSubjectIds: undefined }
  ]) {
    assert.throws(() => search.readSequenceSearchReply({ ...good, ...bad }), (error) => error.code === 'RESULT_INVALID', JSON.stringify(bad));
  }
});
