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
const projectSearch = (committed, mode, searchSetIds) => request.projectLiveSettingWrites({
  committed, drawing: modeDrawing(state, mode), requestPaths: [ANNOTATIONS], inputs, searchSetIds
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
  assert.equal(search.nextSearchSetId([{ id: 'Search 1' }, { id: 'Search 3' }]), 'Search 2');
});

test('the search reads each drawn record once over the bases it draws, and names its record by ID or place', () => {
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
    { recordKey: 'a', recordId: 'X', start: 1, end: 100, length: 100, wrapsOrigin: false, recordLength: 100, header: 'X:1-100', sequence: '' },
    { recordKey: 'b', recordId: 'X', start: 41, end: 70, length: 30, wrapsOrigin: false, recordLength: 100, header: 'X:41-70', sequence: '' }
  ] };
  assert.deepEqual(search.searchSubjects(read, records), [
    { subjectId: 'X:1-100', record: { kind: 'recordIndex', index: 0 }, recordKey: 'a', sourceStart: 1, wrapsOrigin: false, recordLength: 100 },
    { subjectId: 'X:41-70', record: { kind: 'recordIndex', index: 1 }, recordKey: 'b', sourceStart: 41, wrapsOrigin: false, recordLength: 100 }
  ]);
});

// OV-2571, D-52: a record that a span may cross the origin of (circular, drawn whole), in either mode, is also read
// over its origin window, which Python places (origin_window_span); its two subjects still name the record by ID.
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
    { recordKey: 'c', recordId: 'C', start: 1, end: 100, length: 100, wrapsOrigin: false, recordLength: 100, header: 'C:1-100', sequence: '' },
    { recordKey: 'c', recordId: 'C', start: 71, end: 30, length: 60, wrapsOrigin: true, recordLength: 100, header: 'C:71-30', sequence: '' }
  ] };
  assert.deepEqual(search.searchSubjects(read, records), [
    { subjectId: 'C:1-100', record: { kind: 'recordId', value: 'C' }, recordKey: 'c', sourceStart: 1, wrapsOrigin: false, recordLength: 100 },
    { subjectId: 'C:71-30', record: { kind: 'recordId', value: 'C' }, recordKey: 'c', sourceStart: 71, wrapsOrigin: true, recordLength: 100 }
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
