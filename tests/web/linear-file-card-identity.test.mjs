// Two Linear File cards whose files have the same bytes stay two File cards after Save and Load (OV-3703). Save
// stores one resource for equal bytes, so each Linear row binding names its File card (`fileUid`); Load gives the
// rows of one card one file and each card its own source. What reads the File grouping then acts on one card:
// the card list, Move, Replace, Remove and the LOSAT job plan; and the request names one resource per card, as before
// Save (OV-3710). A Session without `fileUid` keeps the rule of one source per resource.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import {
  gallerySessionDocument, readSavedSession, sessionLoadEvent, sessionLoadEventOfBytes
} from './helpers/session-documents.mjs';
import { validateWebConfigOverrides, pythonLegacySessionUpgrader } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const { createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { groupLinearSourceRecords, planLinearSourceRemoval, planLosatSourceJobs } = await import(
  '../../gbdraw/web/js/services/linear-sources.js'
);
const { createLinearSourceFileReplacement } = await import('../../gbdraw/web/js/app/linear-source-edits.js');
const { planLinearSourceRowMove } = await import('../../gbdraw/web/js/services/linear-record-layout.js');
const { webFileBindingFaults } = await import('../../gbdraw/web/js/services/session-resource-backing.js');
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
const { readOfDrawing } = (await import('../../gbdraw/web/js/services/artifact-slot.js')).createArtifactSlotStore();

const TWO_RECORDS = readFileSync(new URL('../fixtures/web_batch_two_records.gb', import.meta.url), 'utf8');
const LAMBDA = await gallerySessionDocument('lambda_basic_linear');
const linear = () => modeDrawing(state, 'linear');
const load = async (document) => {
  const loaded = await importSession(await sessionLoadEvent(document), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.ok(['ok', 'legacy'].includes(loaded.status), JSON.stringify(loaded.error));
};
const save = async () => {
  const result = await exportSession('cards', { askChoice: async () => 'download' });
  assert.equal(result?.status, 'saved', JSON.stringify(result?.error ?? result));
  return readSavedSession(result.blob);
};
/** The Linear drawing's File cards, each as the uids of its rows. */
const cards = () => groupLinearSourceRecords(linear().linearSeqs).map(({ records }) => records.map(({ sequence }) => sequence.uid));
/** The Linear drawing's row bindings as Save wrote them. */
const linearBindings = (document) => document.drawings.find((drawing) => drawing.mode === 'linear').webFiles.bindings.linearSeqs;

/**
 * The Linear drawing with card A (a.gb, its two records) and card B (b.gb, one record), the two files chosen
 * separately with the same bytes; Save then stores one resource for them.
 */
const twoCardsOfEqualBytes = async () => {
  await load(structuredClone(LAMBDA));
  const fileA = state.inputFiles.registerFile(new File([TWO_RECORDS], 'a.gb', { type: 'text/plain', lastModified: 1000 }));
  const fileB = state.inputFiles.registerFile(new File([TWO_RECORDS], 'b.gb', { type: 'text/plain', lastModified: 1000 }));
  linear().linearSeqs.splice(0, Infinity,
    createLinearSeq({ uid: 'card-a', gb: fileA, region_record_id: 'TESTA' }),
    createLinearSeq({ uid: 'card-a-2', gb: fileA, region_record_id: 'TESTB' }),
    createLinearSeq({ uid: 'card-b', gb: fileB, region_record_id: 'TESTA' }));
  assert.deepEqual(cards(), [['card-a', 'card-a-2'], ['card-b']]);
};

test('two File cards of equal bytes stay two after Save and Load', async () => {
  await twoCardsOfEqualBytes();
  const document = await save();
  const rows = linearBindings(document);
  assert.deepEqual(rows.map(({ fileUid }) => fileUid), ['card-a', 'card-a', 'card-b'], 'each row names its File card');
  assert.equal(new Set(rows.map(({ gb }) => gb.resourceId)).size, 1, 'equal bytes are stored once');
  await load(document);
  assert.deepEqual(cards(), [['card-a', 'card-a-2'], ['card-b']]);
  const [a, a2, b] = linear().linearSeqs;
  assert.equal(a.gb, a2.gb, 'the rows of one card share its file');
  assert.notEqual(a.gb.sourceId, b.gb.sourceId);
  assert.deepEqual([a.gb.name, b.gb.name], ['a.gb', 'b.gb']);
  // A second Save of the loaded drawing writes the same cards.
  assert.deepEqual(linearBindings(await save()).map(({ fileUid }) => fileUid), ['card-a', 'card-a', 'card-b']);
});

test('after Load, Move, Replace and Remove act on one card, and the LOSAT plan searches each card as its own file', async () => {
  await twoCardsOfEqualBytes();
  await load(await save());
  // Move (`planLinearSourceRowMove`): card B moves above card A, with its own rows.
  const move = planLinearSourceRowMove({
    sourceGroups: groupLinearSourceRecords(linear().linearSeqs), entries: linear().linearRecordRows, sourceIndex: 1, direction: -1
  });
  assert.equal(move.allowed, true, JSON.stringify(move));
  // Remove (`planLinearSourceRemoval`): card B alone.
  const removal = planLinearSourceRemoval({ sequences: linear().linearSeqs, sourceUid: 'card-b', intent: 'delete' });
  assert.equal(removal.allowed, true);
  assert.deepEqual(removal.allowed && removal.removedUids, ['card-b']);
  // The LOSAT job of a card A record against card B's record searches card B's file, not card A's other record
  // (the database per File).
  const { jobs } = planLosatSourceJobs({
    sequences: linear().linearSeqs, specs: [{ queryIndex: 0, subjectIndex: 2 }], buildArgs: () => []
  });
  assert.deepEqual(jobs.map(({ scope, queryIndexes, subjectIndexes }) => ({ scope, queryIndexes, subjectIndexes })),
    [{ scope: 'between-sources', queryIndexes: [0, 1], subjectIndexes: [2] }]);
  // Replace of card B's file keeps card A.
  const { replaceLinearSourceFile } = createLinearSourceFileReplacement({
    state, shownAction: (_mode, action) => action(linear()),
    pendingDiscoveryOf: () => ({ expansions: new Map(), uploads: new Set() }),
    applyLinearSeqMutation: (drawing, items) => { drawing.linearSeqs.splice(0, Infinity, ...items); },
    askLinearFileRemoval: () => { throw new Error('no question asked'); },
    sourceReplacementRetiresPresentation: () => false
  });
  replaceLinearSourceFile(linear(), 'card-b', 'gb', new File(['LOCUS other\n//\n'], 'c.gb'));
  assert.deepEqual(cards(), [['card-a', 'card-a-2'], ['card-b']]);
  assert.deepEqual(linear().linearSeqs.map(({ gb }) => gb.name), ['a.gb', 'a.gb', 'c.gb']);
});

test('the request after Load reads each File card from a resource of its own, as before Save', async () => {
  // The request builder that Generate and a draft Save share (`buildCanonicalRenderRequest` through
  // `projectDraftRequest`): Python reads each resource as one source file, so one resource for both cards would
  // make the two files one genome to its LOSAT planner (OV-3710).
  const sources = async () => {
    const { envelope, filesByResource } = await sessionConfig.projectDraftRequest(state, linear(), {
      linearRecordCatalog: null, recordDisplayRows: [], omittedAnnotationRecordKeys: [], readOfDrawing
    });
    const ids = envelope.renderRequest.records.map(({ source }) => source.resourceId);
    return { ids, names: [...new Set(ids)].map((id) => filesByResource.get(id)?.name) };
  };
  await twoCardsOfEqualBytes();
  const before = await sources();
  assert.equal(new Set(before.ids).size, 2, JSON.stringify(before));
  await load(await save());
  const after = await sources();
  assert.equal(after.ids[0], after.ids[1], 'the records of card A read one resource');
  assert.notEqual(after.ids[0], after.ids[2], JSON.stringify(after));
  assert.deepEqual(after.names, ['a.gb', 'b.gb']);
});

test('the positive fixture, a Web Save of two cards of equal bytes, loads as two File cards', async () => {
  const fixture = readFileSync(new URL('../fixtures/sessions/web-linear-equal-file-cards.v51.gbdraw', import.meta.url));
  const loaded = await importSession(await sessionLoadEventOfBytes(fixture), { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.ok(['ok', 'legacy'].includes(loaded.status), JSON.stringify(loaded.error));
  assert.equal(cards().length, 2);
  assert.deepEqual(linear().linearSeqs.map(({ gb }) => gb.name), ['a.gb', 'b.gb']);
});

test('a Session whose rows name no File card keeps one source per resource', async () => {
  await twoCardsOfEqualBytes();
  const document = await save();
  linearBindings(document).forEach((row) => { delete row.fileUid; });
  await load(document);
  assert.deepEqual(cards(), [['card-a', 'card-a-2', 'card-b']]);
});

test('a File card name that is no text is dropped from its row', async () => {
  await twoCardsOfEqualBytes();
  const document = await save();
  const drawing = document.drawings.find(({ mode }) => mode === 'linear');
  const bindings = drawing.webFiles.bindings;
  for (const fileUid of ['', ' card-a', 7, null]) {
    bindings.linearSeqs[0].fileUid = fileUid;
    assert.deepEqual(webFileBindingFaults(drawing.webFiles, document.resources, 'linearSeqs'),
      [{ at: ['linearSeqs', 0, 'fileUid'], field: true }], JSON.stringify(fileUid));
  }
  bindings.linearSeqs[0].fileUid = 'card-a';
  assert.deepEqual(webFileBindingFaults(drawing.webFiles, document.resources, 'linearSeqs'), []);
});
