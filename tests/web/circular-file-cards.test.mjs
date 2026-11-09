// The Circular input is a list of File cards, one request source per file, as
// the CLI gives it (circular-multi-file D-01, D-08). A CLI Session of several
// `--gbk` files loads as one card per file, keyed by the CLI's record key, the
// Web builds the CLI's request for the same files, and Save writes the cards,
// not one concatenated File (OV-2838).
import { acceptNodeFeatureRecoveryRow, installSessionImportWorker } from './helpers/session-import-node.mjs';
import { pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtemp, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { gunzipSync } from 'node:zlib';
import { test } from 'node:test';
import { createHash } from 'node:crypto';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';

globalThis.window = {
  Vue: {
    ref: value => ({ value }), reactive: value => value,
    computed: getter => ({ get value() { return getter(); } }), nextTick: async () => {}
  },
  DOMPurify: { sanitize: value => value }
};
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({ addEventListener: () => {}, click: () => {}, remove: () => {}, parentNode: null })
};
installFakeSvgDom();
globalThis.File = class File extends Blob {
  constructor(parts, name, options = {}) {
    super(parts, options);
    this.name = name;
    this.lastModified = options.lastModified ?? 0;
  }
};
globalThis.alert = () => {};

installSessionImportWorker();
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { getCommittedCanonicalRenderRequest, serializeActiveRenderFiles, setUnmanagedConfigOverrideValidator } = sessionConfig;
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { readFileBytes } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { readCircularSourceRecords } = await import('../../gbdraw/web/js/services/record-discovery.js');
const { circularSourceInputs } = await import('../../gbdraw/web/js/services/input-types.js');
const { circularDiscoveryForInput } = await import('../../gbdraw/web/js/services/record-discovery.js');
setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));

const root = process.cwd();
const twoRecords = path.join(root, 'tests/fixtures/web_batch_two_records.gb');
const mito = path.join(root, 'tests/fixtures/sessions/cli-web-mito.gb');
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const load = (bytes) => importSession({ target: {
  files: [new File([bytes], 'current.json', { type: 'application/json' })], value: 'selected'
} }, { transformLegacyResultSvg, confirmLoadRepairs: acceptNodeFeatureRecoveryRow });

// What a request record says about which record it draws.
const shapeOf = ({ recordKey, cardinality, selector }) => ({ recordKey, cardinality, selector: selector ?? null });
// The records of the request the next Generate builds (run-analysis.js): it
// reads each card's records first (the one read, `readCircularSourceRecords`).
const nextGenerateRecords = async () => {
  const drawing = state.activeDrawing();
  const reads = await readCircularSourceRecords(circularSourceInputs(drawing.inputType.value, drawing.circularSources),
    (ref) => state.inputFiles.fileOf(ref));
  assert.deepEqual(reads.map(read => read.error), reads.map(() => null));
  drawing.circularRecordList.value = reads.flatMap(read => read.records);
  const filesData = await serializeActiveRenderFiles(state, drawing);
  return buildCanonicalRenderRequest({ state, drawing, filesData }).renderRequest.records;
};

const cliSession = async (directory, args) => {
  const output = path.join(directory, 'cli');
  const file = `${output}.gbdraw-session.json.gz`;
  execFileSync('python', ['-m', 'gbdraw.cli', 'circular', ...args, '-o', output, '--session_output', file], {
    cwd: directory, env: { ...process.env, PYTHONPATH: root }, stdio: 'pipe', timeout: 1_800_000
  });
  return gunzipSync(await readFile(file));
};

// Failing-first test 1 and 4 of DESIGN section 7; OV-2838.
await test('a CLI Session of two --gbk files loads as two File cards and the Web request is the CLI request', async () => {
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-circular-cards-'));
  try {
    const bytes = await cliSession(directory, ['--gbk', twoRecords, mito, '--multi_record_canvas']);
    const session = JSON.parse(bytes);
    const view = savedDrawingView(session, session.ui.activeDrawingId);
    const cliRecords = view.renderRequest.records;
    // The CLI commits the records it resolved from one ALL input per file:
    // `<record key>:<k>` for a file of several records, the key alone for one.
    assert.deepEqual(cliRecords.map(record => record.recordKey), ['record-1:1', 'record-1:2', 'record-2']);
    assert.equal(cliRecords[0].source.resourceId, cliRecords[1].source.resourceId);
    assert.notEqual(cliRecords[0].source.resourceId, cliRecords[2].source.resourceId);
    // The CLI writes one card per file, no composite.
    const bindings = view.webFiles.bindings;
    assert.equal(Object.hasOwn(bindings, 'c_gb'), false, 'no c_gb binding');
    assert.deepEqual(bindings.circularSources.map(row => row.uid), ['record-1', 'record-2']);
    assert.equal(JSON.stringify(session).includes('"composite"'), false);

    const result = await load(bytes);
    assert.equal(result.status, 'ok', result.error?.stack);
    const drawing = modeDrawing(state, 'circular');
    assert.deepEqual(drawing.circularSources.map(row => row.uid), ['record-1', 'record-2']);
    const cardBytes = await Promise.all(drawing.circularSources.map(async row => hash(await readFileBytes(state.inputFiles.fileOf(row.gb)))));
    assert.deepEqual(cardBytes, [hash(await readFile(twoRecords)), hash(await readFile(mito))]);

    // A Generate from the cards builds the CLI's records: the same keys,
    // cardinalities and selectors, each file its own source (OV-2838: the
    // composite gave one concatenated source whose records were keyed record-1..3).
    const records = await nextGenerateRecords();
    assert.deepEqual(records.map(shapeOf), cliRecords.map(shapeOf));
    assert.equal(records[0].source.resourceId, records[1].source.resourceId);
    assert.notEqual(records[0].source.resourceId, records[2].source.resourceId);
    assert.deepEqual(getCommittedCanonicalRenderRequest().records.map(shapeOf), cliRecords.map(shapeOf));
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
});

// OV-2844 (D-13): right after Load, Generate reads the records of every card
// before it builds its request, so a CLI Session drawn as separate figures
// (batch) gives the CLI's records, never one ALL record per card, which Python
// refuses in a batch.
await test('Generate right after loading a CLI Session in batch grouping writes the CLI records', async () => {
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-circular-cards-batch-'));
  try {
    const bytes = await cliSession(directory, ['--gbk', twoRecords, mito]);
    const session = JSON.parse(bytes);
    const cliRequest = savedDrawingView(session, session.ui.activeDrawingId).renderRequest;
    assert.equal(cliRequest.grouping, 'batch');
    assert.equal((await load(bytes)).status, 'ok');
    // Load seeds no records from the committed request.
    assert.deepEqual(modeDrawing(state, 'circular').circularRecordList.value, []);
    // OV-2881: the cards' read waits for Generate (as dev's one input did), so
    // no automatic read after Load writes the Record Order into the draft.
    const loaded = modeDrawing(state, 'circular');
    assert.equal(loaded.circularRecordDiscovery.status, 'deferred');
    assert.equal(circularDiscoveryForInput(loaded).current, true);
    const records = await nextGenerateRecords();
    assert.deepEqual(records.map(shapeOf), cliRequest.records.map(shapeOf));
    assert.equal(records.some(record => record.cardinality === 'all'), false);
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
});

// PD-OI-044 scenario revision 5: a records-table crop of a Circular grid or batch is saved as its file and range;
// the Web loads the range by the record's card and draws it in the next Generate, as the CLI drew it.
for (const grouping of ['grid', 'batch']) {
  await test(`a CLI ${grouping} Session of records-table crops loads their ranges and Generate writes the CLI records`, async () => {
    const directory = await mkdtemp(path.join(tmpdir(), `gbdraw-circular-ranges-${grouping}-`));
    try {
      const { writeFile } = await import('node:fs/promises');
      const table = path.join(directory, 'records.tsv');
      await writeFile(table, 'gbk\trecord_id\tregion\n'
        + `${twoRecords}\tTESTA\t101-500:rc\n${twoRecords}\tTESTB\t\n${mito}\t\t1001-3000\n`);
      const bytes = await cliSession(directory, ['--records_table', table, ...(grouping === 'grid' ? ['--multi_record_canvas'] : [])]);
      const session = JSON.parse(bytes);
      const cliRequest = savedDrawingView(session, session.ui.activeDrawingId).renderRequest;
      assert.equal(cliRequest.grouping, grouping);
      const windowsOf = ({ windows }) => windows.map(({ start, end, reverseComplement }) => ({ start, end, reverseComplement }));
      assert.deepEqual(cliRequest.records.map(windowsOf), [
        [{ start: 101, end: 500, reverseComplement: true }], [], [{ start: 1001, end: 3000, reverseComplement: false }]
      ]);
      assert.equal((await load(bytes)).status, 'ok');
      const records = await nextGenerateRecords();
      // A table row is keyed by its row (`record-<n>`); the Web keys each record by its card (`cardUids`).
      const selectorOf = ({ cardinality, selector }) => ({ cardinality, selector: selector ?? null });
      assert.deepEqual(records.map(selectorOf), cliRequest.records.map(selectorOf));
      assert.deepEqual(records.map(windowsOf), cliRequest.records.map(windowsOf));
    } finally {
      await rm(directory, { recursive: true, force: true });
    }
  });
}

// The Web side of the cards (DESIGN section 7, failing-first tests 2, 3, 6, 7,
// 8 and 9): a record's identity is its card's uid with its place in that
// card's file, so edits of other cards move none of its rows.
const {
  drawnRecordIdentities, recordKeyMoves, resliceRecordRows
} = await import('../../gbdraw/web/js/services/record-draw-selection.js');
const { readConservationSubject } = await import('../../gbdraw/web/js/services/circular-conservation-keys.js');
const { createCircularSourceEdits } = await import('../../gbdraw/web/js/app/circular-record-controls.js');
const { createDrawingState } = await import('../../gbdraw/web/js/state.js');

const entries = (...ids) => ids.map((recordId, sourceIndex) => ({ sourceIndex, selector: `#${sourceIndex + 1}`, recordId }));
const circularSources = (cards) => ({ mode: 'circular', cards: new Map(cards) });
const allOf = (...uids) => uids.map((recordKey) => ({ recordKey, cardinality: 'all', windows: [] }));
const CARD_A = ['card-a', entries('chr', 'pA')];
const CARD_B = ['card-b', entries('pB')];
const CARD_C = ['card-c', entries('pC')];
// The moves that take a row away from a record's own key (the keys an older Session wrote may still bring rows
// to it, the last test).
const movesOf = (moves, keys) => [...moves].filter(([from, to]) => from !== to && keys.includes(from));

await test('adding, removing or moving a card moves no row of the other cards', () => {
  const before = drawnRecordIdentities(allOf('card-a', 'card-b'), circularSources([CARD_A, CARD_B]));
  assert.deepEqual(before.map(({ requestKey }) => requestKey), ['card-a:1', 'card-a:2', 'card-b']);
  const keys = ['card-a:1', 'card-a:2', 'card-b'];
  for (const cards of [[CARD_B], [CARD_B, CARD_A], [CARD_A, CARD_B, CARD_C], [CARD_C, CARD_A, CARD_B]]) {
    const after = drawnRecordIdentities(allOf(...cards.map(([uid]) => uid)), circularSources(cards));
    assert.deepEqual(movesOf(recordKeyMoves(before, after), keys), [], cards.map(([uid]) => uid).join());
  }
});

await test('turning a record of a card OFF moves no row', () => {
  const sources = circularSources([CARD_A, CARD_B]);
  const before = drawnRecordIdentities(allOf('card-a', 'card-b'), sources);
  // A card with an OFF record draws its ON records one by one, by identity.
  const after = drawnRecordIdentities([
    { recordKey: 'card-a:1', cardinality: 'exactly_one', selector: { kind: 'recordIndex', index: 0 }, windows: [] },
    { recordKey: 'card-b', cardinality: 'all', windows: [] }
  ], sources);
  assert.deepEqual(after.map(({ identity }) => identity), ['card-a:1', 'card-b']);
  assert.deepEqual(movesOf(recordKeyMoves(before, after), ['card-a:1', 'card-a:2', 'card-b']), []);
});

await test('a record-major Depth matrix follows its records by identity', () => {
  const previous = ['card-a:1', 'card-a:2', 'card-b'];
  const rows = [['a1'], ['a2'], ['b']];
  const empty = () => [];
  assert.deepEqual(resliceRecordRows(rows, previous, ['card-b', 'card-a:1', 'card-a:2'], empty), [['b'], ['a1'], ['a2']]);
  assert.deepEqual(resliceRecordRows(rows, previous, ['card-b'], empty), [['b']]);
  assert.deepEqual(resliceRecordRows(rows, previous, [...previous, 'card-c'], empty), [...rows, []]);
  // A matrix shorter than the list: the records past its end have no row yet.
  assert.deepEqual(resliceRecordRows([['a1']], previous, ['card-b', 'card-a:1'], empty), [[], ['a1']]);
});

await test('the ring subject is every card whole, in card order, and refuses an ID two cards hold', async () => {
  const fasta = (name, text) => new File([text], name, { type: 'text/plain' });
  const chromosome = fasta('chr.fasta', '>chr\nACGTACGTAC\n');
  const plasmid = fasta('p.fasta', '>p1\nGGGCCC\n>p2\nTTTAAA\n');
  const subject = await readConservationSubject([chromosome, plasmid], 'fasta');
  const one = await readConservationSubject([chromosome], 'fasta');
  const other = await readConservationSubject([plasmid], 'fasta');
  assert.equal(subject.fasta, one.fasta + other.fasta);
  // The hash the CLI's ring of the same files has: the SHA-256 of that FASTA.
  assert.equal(subject.hash, hash(subject.fasta));
  await assert.rejects(readConservationSubject([plasmid, fasta('again.fasta', '>p2\nCCCC\n')], 'fasta'),
    (error) => error.code === 'COMPARISON_INPUT' && error.context.reason === 'UNIQUE_IDS');
});

await test('one pick of several files is one History step that adds a card per file', async () => {
  const drawing = createDrawingState('circular-cards-history', 'circular');
  const steps = [];
  const edits = createCircularSourceEdits({
    state,
    CIRCULAR_RECORD_FORM_DEFAULTS: { circular_record_selector: '', circular_region_start: null, circular_region_end: null,
      circular_reverse: false, circular_record_label: '', circular_record_subtitle: '' },
    changeRecordSources: (_drawing, mutate) => mutate(),
    pendingDiscoveryOf: () => ({ circularUpload: null }),
    discoverDrawingRecords: () => {},
    shownAction: (_mode, action) => action(drawing),
    runUndoable: async (label, fn) => { steps.push(label); return fn(); }
  });
  const before = drawing.circularSources.length;
  await edits.addCircularSourceFiles(['a.gb', 'b.gb', 'c.gb'].map((name) => new File(['LOCUS'], name)));
  assert.deepEqual(steps, ['Add File']);
  assert.deepEqual(drawing.circularSources.slice(before).map((card) => card.gb?.name), ['a.gb', 'b.gb', 'c.gb']);
  // GFF3 + FASTA: the files of one pick pair by stem (D-07), still one step.
  drawing.inputType.value = 'gff';
  const gffBefore = drawing.circularSources.length;
  await edits.addCircularSourceFiles(['a.fasta', 'a.gff3', 'b.gff3'].map((name) => new File(['x'], name)));
  assert.deepEqual(steps, ['Add File', 'Add File']);
  assert.deepEqual(drawing.circularSources.slice(gffBefore).map((card) => [card.gb, card.gff?.name ?? null, card.fasta?.name ?? null]),
    [[null, 'a.gff3', 'a.fasta'], [null, 'b.gff3', null]]);
});

// #1809 review 1 (M2, D-159): a Session whose Linear card holds genome.fna as its
// FASTA and annotation.txt as its GFF3, and an empty Circular drawing in GFF3 +
// FASTA. The empty zone offers both files, and each one chosen goes to the slot
// of the kind the Session binds it as, whatever its name says (by name,
// annotation.txt would be a FASTA).
await test('a Session file chosen in the empty Add files zone goes to the slot of its kind, not by its name', async () => {
  const { bindLinearSeqFile, sessionFilesOfKind } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
  const { createCircularFileCards } = await import('../../gbdraw/web/js/app/circular-file-cards.js');
  const linear = createDrawingState('circular-cards-session-kind-linear', 'linear');
  const row = linear.linearSeqs[0].uid;
  bindLinearSeqFile(state.inputFiles, linear.slices.inputs, row, 'fasta', new File(['>c\nACGT\n'], 'genome.fna'));
  bindLinearSeqFile(state.inputFiles, linear.slices.inputs, row, 'gff', new File(['##gff-version 3\n'], 'annotation.txt'));
  const [fasta, gff] = [linear.linearSeqs[0].fasta, linear.linearSeqs[0].gff];
  const drawing = createDrawingState('circular-cards-session-kind', 'circular');
  drawing.inputType.value = 'gff';
  // The Session of the two drawings: what the zone offers and the writer reads the kinds of.
  const session = { inputFiles: state.inputFiles, drawingList: () => [linear, drawing] };
  const edits = createCircularSourceEdits({
    state: session,
    CIRCULAR_RECORD_FORM_DEFAULTS: { circular_record_selector: '', circular_region_start: null, circular_region_end: null,
      circular_reverse: false, circular_record_label: '', circular_record_subtitle: '' },
    changeRecordSources: (_drawing, mutate) => mutate(),
    pendingDiscoveryOf: () => ({ circularUpload: null }),
    discoverDrawingRecords: () => {},
    shownAction: (_mode, action) => action(drawing),
    runUndoable: async (_label, fn) => fn()
  });
  const cards = createCircularFileCards({
    shownDrawingOfMode: () => drawing,
    circularRecordDiscoveryState: { value: { status: 'idle', error: '', failedSources: [] } },
    addCircularSourceFiles: edits.addCircularSourceFiles,
    removeCircularSourceCard: async () => {}
  });
  // The empty zone offers the Session's files of both kinds the input type reads.
  assert.deepEqual(cards.addZone.value.fileKinds, ['gff', 'fasta']);
  assert.deepEqual(sessionFilesOfKind(session, cards.addZone.value.fileKinds).map((ref) => ref.name), ['annotation.txt', 'genome.fna']);
  await cards.addFiles([gff]);
  await cards.addFiles([fasta]);
  assert.deepEqual(drawing.circularSources.map((card) => [card.gff?.name ?? null, card.fasta?.name ?? null]),
    [['annotation.txt', null], [null, 'genome.fna']]);
  assert.equal(drawing.circularSources[0].gff.fileId, gff.fileId, 'the chosen ref itself is bound; no file is stored twice');
  // Once a card shows, the zone offers none: the cards' own uploaders do (D-156).
  assert.deepEqual(cards.addZone.value.fileKinds, []);
});

await test('rows an older Session left under record-<n> move to their record on the first Generate', () => {
  // Session 50 and older named a Circular record by its place among all the drawing's records.
  // The cards as Load installed them (OV-2867).
  const next = drawnRecordIdentities(allOf('card-a', 'card-b'), { ...circularSources([CARD_A, CARD_B]), legacyOrder: ['card-a', 'card-b'] });
  const moves = recordKeyMoves([], next);
  assert.equal(moves.get('record-1'), 'card-a:1');
  assert.equal(moves.get('record-2'), 'card-a:2');
  assert.equal(moves.get('record-3'), 'card-b');
  assert.deepEqual([...moves].filter(([from]) => from.startsWith('circular-pB-')).map(([, to]) => to), ['card-b']);
});

// OV-2843 (D-12): a Session 50 that drew a composite of a.gb (two records) and
// b.gb as one concatenated source keyed its records record-1..3 over all of
// them. Its cards are `source-<n>`, a key no Session 50 has, so a row under
// record-2 (a.gb's second record) moves to that record on the first Generate,
// not to b.gb's record, as it would if the cards were record-1 and record-2.
await test('a row a composite Session keyed record-2 lands on its record after the first Generate', () => {
  const cards = [['source-1', entries('chr', 'pA')], ['source-2', entries('pB')]];
  const moves = recordKeyMoves([], drawnRecordIdentities(allOf('source-1', 'source-2'), { ...circularSources(cards), legacyOrder: ['source-1', 'source-2'] }));
  assert.equal(moves.get('record-2'), 'source-1:2');
  assert.equal(moves.get('record-3'), 'source-2');
  const colliding = [['record-1', entries('chr', 'pA')], ['record-2', entries('pB')]];
  const collided = recordKeyMoves([], drawnRecordIdentities(allOf('record-1', 'record-2'), circularSources(colliding)));
  assert.notEqual(collided.get('record-2'), 'record-1:2');
});

// OV-2880: the Circular annotation source key names each card by its uid and
// its files' fingerprints (name, size), never by a per-page file ID, so a
// Session saved and loaded again (new file IDs) keeps its annotations' records;
// one card has dev's key, so Sessions dev wrote still match.
const { createDrawingRecordCatalog } = await import('../../gbdraw/web/js/app/drawing-record-catalog.js');
const { annotationSourceKey } = await import('../../gbdraw/web/js/app/annotations/record-catalog.js');
const { addCircularSource, bindCircularSourceFile, createCircularSource } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const cardDrawing = (id, cards) => {
  const drawing = createDrawingState(id, 'circular');
  cards.forEach(([uid, name, text]) => {
    addCircularSource(drawing.slices.inputs, createCircularSource({ uid }));
    bindCircularSourceFile(state.inputFiles, drawing.slices.inputs, uid, 'gb', new File([text], name));
  });
  return drawing;
};
const { circularAnnotationSourceKey } = createDrawingRecordCatalog({
  state, readInputRecords: async () => [], linearRecordCatalogStateFor: () => ({ status: 'idle', records: [], error: '' }),
  linearRecordErrorModelFor: () => '', refreshLinearRecords: async () => ({ status: 'ok' })
});

await test('the Circular annotation source key survives new file IDs and is dev\'s key for one card', () => {
  const one = [['record-1', 'a.gb', 'LOCUS a']];
  const saved = cardDrawing('annotation-key-saved', one);
  const loaded = cardDrawing('annotation-key-loaded', one);
  assert.notEqual(saved.circularSources[0].gb.fileId, loaded.circularSources[0].gb.fileId);
  assert.equal(circularAnnotationSourceKey(loaded), circularAnnotationSourceKey(saved));
  assert.equal(circularAnnotationSourceKey(saved), annotationSourceKey({
    scope: 'circular', inputType: 'gb', primaryFile: saved.circularSources[0].gb, pairedFile: null
  }));
  const two = [['source-1', 'a.gb', 'LOCUS a'], ['source-2', 'b.gb', 'LOCUS bb']];
  assert.equal(circularAnnotationSourceKey(cardDrawing('annotation-key-two-a', two)),
    circularAnnotationSourceKey(cardDrawing('annotation-key-two-b', two)));
  // Another file in a card, or the cards in another order, is another source.
  const keyOf = (cards, id) => circularAnnotationSourceKey(cardDrawing(id, cards));
  assert.notEqual(keyOf([two[0], ['source-2', 'c.gb', 'LOCUS bb']], 'annotation-key-other'), keyOf(two, 'annotation-key-two-c'));
  assert.notEqual(keyOf([two[1], two[0]], 'annotation-key-moved'), keyOf(two, 'annotation-key-two-d'));
});

// OV-2867: Load keeps the cards an older Session moved its Circular input to
// (not saved), so its flat record positions count over them in that order.
await test('Load of a Session 50 or older keeps the cards its record positions count over', async () => {
  const { cardsTakeLegacyRecordPositions } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
  const bytes = gunzipSync(await readFile(path.join(root, 'tests/fixtures/sessions/composite-circular-three-files.v44-schema8.gbdraw-session.json.gz')));
  assert.equal((await load(bytes)).status, 'ok');
  const drawing = modeDrawing(state, 'circular');
  const uids = drawing.circularSources.map((card) => card.uid);
  // The Session 44 drew a composite of three files as one concatenated source (D-12: `source-<n>`).
  assert.deepEqual(uids, ['source-1', 'source-2', 'source-3']);
  assert.equal(cardsTakeLegacyRecordPositions(drawing.circularSources), true);
  assert.deepEqual(drawing.files.c_legacy_record_cards?.map(({ uid, primaryFile }) => [uid, primaryFile]),
    drawing.circularSources.map(({ uid, gb }) => [uid, gb]));
});

// OV-2869: a Session 46-50 Web composite keeps the `#N` its move cannot place
// (drawings-vectors.json). A Save before any read writes them through as a
// Session 51; its Load keeps the cards as loaded too, so the first read
// places them over those cards.
await test('a Session 51 that holds a kept flat #N keeps its cards as loaded', async () => {
  const vectors = JSON.parse(await readFile(path.join(root, 'tests/fixtures/sessions/drawings-vectors.json'), 'utf8'));
  const vector = vectors.upgrades.find(({ name }) => name.includes('as a Web composite with #N selectors (OV-2869)'));
  const input = JSON.parse(gunzipSync(await readFile(path.join(root, vector.input))).toString('utf8'));
  for (const { path: [...parents], value } of vector.edits) {
    const key = parents.pop();
    parents.reduce((holder, part) => holder[part], input)[key] = value;
  }
  assert.equal((await load(Buffer.from(JSON.stringify(input)))).status, 'ok');
  const drawing = modeDrawing(state, 'circular');
  assert.deepEqual([...drawing.recordsOff], ['#2']);
  assert.ok(drawing.files.c_legacy_record_cards);
  const result = await sessionConfig.exportSession('flat-selectors', { shownDrawing: state.activeDrawing, askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result));
  const saved = JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
  assert.deepEqual(saved.drawings.find(({ id }) => id === 'circular').config.recordsOff, ['#2']);
  assert.equal((await load(Buffer.from(JSON.stringify(saved)))).status, 'ok');
  const reloaded = modeDrawing(state, 'circular');
  assert.equal(reloaded.form.circular_record_selector, '#3');
  assert.deepEqual(reloaded.files.c_legacy_record_cards?.map(({ uid }) => uid), reloaded.circularSources.map(({ uid }) => uid));
});
