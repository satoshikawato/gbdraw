// OV-3592 through the shipped chain: Load (config.js importSession) gives each comparison-plan edge the name its
// request stored for the pair's table (`loadedSourceName`), and the draft request a Save or a Generate writes
// (`projectDraftRequest`: `serializeActiveRenderFiles`, then `buildCanonicalRenderRequest`) names each table so, also
// after Undo and Redo. Two comparisons share one table: the user named one ("hits.tsv"); gbdraw made the other
// (`sourceName: null`), which no file name may replace.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { sessionLoadEvent, readSavedSession } from './helpers/session-documents.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { body: { appendChild: () => {} }, createElement: () => ({ click: () => {}, addEventListener: () => {}, parentNode: null }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
const { createArtifactSlotStore } = await import('../../gbdraw/web/js/services/artifact-slot.js');
const { createHistoryRouter } = await import('../../gbdraw/web/js/services/history.js');
const { createHistorySnapshotService } = await import('../../gbdraw/web/js/services/history-snapshot.js');
const { setComparisonEdgeFile } = await import('../../gbdraw/web/js/services/linear-comparisons.js');

/**
 * A Web-saved Session 51 of three Linear records whose two pairs bind one table (as two empty tables, or one table
 * a Python Save stored once, do): built from drawings-linear-pair-losat.v51 (two records and one LOSAT table).
 */
const sharedTableSession = async () => {
  const session = await readSavedSession(new Blob([
    readFileSync(new URL('../fixtures/sessions/drawings-linear-pair-losat.v51.gbdraw', import.meta.url))
  ]));
  const [drawing] = session.drawings;
  session.drawings = [drawing];
  session.ui.activeDrawingId = drawing.id;
  const { records, comparisons: [table] } = drawing.renderRequest;
  records.push({ ...structuredClone(records[1]), recordKey: 'linear-seq-third', presentation: { ...records[1].presentation, gridRow: 3 } });
  const { bindings } = drawing.webFiles;
  bindings.linearSeqs.push({ ...structuredClone(bindings.linearSeqs[1]), uid: 'linear-seq-third' });
  drawing.webFiles.linearRecordMetadata.push({ recordKey: 'linear-seq-third', losatGencode: 1 });
  drawing.config.linearRecordLayout.rows.push({ uid: 'linear-seq-third', row: 3 });
  const uids = bindings.linearSeqs.map((sequence) => sequence.uid);
  drawing.renderRequest.comparisons = [
    { ...table, sourceName: 'hits.tsv' },
    { ...table, queryRecordIndex: 1, subjectRecordIndex: 2, sourceName: null }
  ];
  // The plan's own edge IDs: a Web file binding names a pair's table by its edge, the request by the pair's records.
  drawing.config.linearComparisonPlan = {
    mode: 'selected', defaultSource: 'losat', edges: ['edge-ab', 'edge-bc'].map((id, index) => ({
      id, queryUid: uids[index], subjectUid: uids[index + 1], included: true, fileActive: true,
      losatFilenameActive: false, source: 'upload', losatFilename: ''
    }))
  };
  bindings.linearComparisons = [['edge-ab', 'hits.tsv'], ['edge-bc', 'comparison-losat-1.tsv']].map(([id, name]) => ({
    id, file: { resourceId: table.resourceId, name, type: 'text/tab-separated-values', lastModified: 0 }
  }));
  // A drawing without a Result: the feature catalog and run metadata describe the two records the fixture drew.
  drawing.results = [];
  delete drawing.runMetadata;
  delete drawing.editorState.featureCatalog;
  return session;
};

/** Each nucleotide comparison of the shown drawing's draft request, by its query record: the name it stores. */
const draftNames = async () => {
  const parser = globalThis.DOMParser;
  globalThis.DOMParser = class { parseFromString() { return { documentElement: { localName: 'none' } }; } };
  try {
    const { renderRequest } = (await sessionConfig.projectDraftRequest(state, state.activeDrawing(), {
      linearRecordCatalog: null, recordDisplayRows: [], omittedAnnotationRecordKeys: [],
      readOfDrawing: createArtifactSlotStore().readOfDrawing
    })).envelope;
    return renderRequest.comparisons
      .map((comparison) => [comparison.queryRecordIndex, Object.hasOwn(comparison, 'sourceName') ? comparison.sourceName : 'absent'])
      .sort(([left], [right]) => left - right);
  } finally {
    globalThis.DOMParser = parser;
  }
};

const STORED = [[0, 'hits.tsv'], [1, null]];

test('a Load, then Undo and Redo, keep the name the request stored for each comparison of a shared table (OV-3592)', async () => {
  const loaded = await importSession(await sessionLoadEvent(await sharedTableSession()), { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const drawing = state.activeDrawing();
  assert.deepEqual(drawing.linearComparisonPlan.edges.map((edge) => [edge.id, edge.loadedSourceName]), [
    ['edge-ab', { sourceName: 'hits.tsv' }], ['edge-bc', { sourceName: null }]
  ]);
  assert.deepEqual(await draftNames(), STORED);

  // The service is given the store's registerFile alone: a restore gives each plan edge its file back through it.
  /** @type {unknown[]} */
  const registered = [];
  /** @type {typeof state.inputFiles.registerFile} */
  const registerFile = (file) => { registered.push(file); return state.inputFiles.registerFile(file); };
  const snapshots = createHistorySnapshotService({
    state, registerFile, currentResults: () => state.results.value,
    buildConfigData: sessionConfig.buildConfigData, applyConfigData: sessionConfig.applyCapturedConfigData
  });
  const history = createHistoryRouter({
    drawingIds: state.drawingList().map(({ id }) => id),
    activeDrawingId: () => state.activeDrawingId.value,
    managerOptions: (drawingId) => snapshots.drawingHistory(drawingId),
    sessionHistory: snapshots.sessionHistory()
  });
  await history.initializeIntentBaseline();
  // A settings edit: its Undo installs the plan again and gives each edge back the file it held (retained files).
  await history.runUndoable('Window', () => { state.activeDrawing().adv.window_size = 1234; });
  await history.undo();
  assert.deepEqual(await draftNames(), STORED, 'Undo of a settings edit');
  await history.redo();
  assert.deepEqual(await draftNames(), STORED, 'Redo of a settings edit');
  // A file the user chooses for the second pair is named by its name; Undo restores the loaded file with the name
  // the request stored for it (the step's files), and Redo the user's file.
  const chosen = state.inputFiles.registerFile(new File(['q\ts\n'], 'mine.tsv', { type: 'text/tab-separated-values' }));
  await history.runUndoable('Table', () => {
    setComparisonEdgeFile(state.activeDrawing().linearComparisonPlan.edges[1], chosen);
  });
  assert.deepEqual(await draftNames(), [[0, 'hits.tsv'], [1, 'mine.tsv']]);
  registered.length = 0;
  await history.undo();
  assert.deepEqual(await draftNames(), STORED, 'Undo of a file choice');
  /** @type {(file: unknown) => unknown} */
  const fileIdOf = (file) => (file && typeof file === 'object' && 'fileId' in file ? file.fileId : null);
  const edgeFileIds = state.activeDrawing().linearComparisonPlan.edges.map(({ file }) => fileIdOf(file));
  assert.ok(edgeFileIds.every((id) => id !== null && registered.map(fileIdOf).includes(id)),
    'the Undo gives each plan edge its file through the port the service was given');
  await history.redo();
  assert.deepEqual(await draftNames(), [[0, 'hits.tsv'], [1, 'mine.tsv']], 'Redo of a file choice');
});
