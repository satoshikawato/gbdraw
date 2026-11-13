// #1777 review 7: the draft's similarity rings through a rerender of the committed request, a Generate after a Load
// that lost a ring's genome, a ring label with a comma, and a ring edit beside another drawing's LOSAT cache rows.
// - OV-3194: a rerender of the committed request (a Rotate, an Align, a placement: R1 (c)) replaces only the committed
//   Session; the draft's ring inputs, series, label list and stack rows stay as each kind of ring edit left them.
// - OV-3195: a Generate after a replay Load that lost ring 0's genome gives each ring the genetic code, label and row
//   it searched it with, and the next Generate and Add Seq read the same rings.
// - OV-3199: a Generate that searched no ring, after a replay Load that lost every ring's genome, leaves the rings.
// - OV-3193: a ring label with a comma names no later ring; OV-3200: nor does an older Session's label list in a replay.
// - OV-3196: a ring edit clears only its drawing's LOSAT cache-info rows.
// Generate runs as production builds it: the rings its ring keys read (`searchedRingSources`), the request's rings
// (`losatRequestRings`) and the request builder (`buildCanonicalRenderRequest`), then the composition root's adoption.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { isDeepStrictEqual } from 'node:util';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { readSavedSession, sessionLoadEvent } from './helpers/session-documents.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
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

const { state } = await import('../../gbdraw/web/js/state.js');
await (await import('./helpers/drawing-state.mjs')).addStartLinearDrawing();
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const {
  createCircularConservationSources, createCircularConservationSeriesPanel, createCircularConservationComparisonFiles
} = await import('../../gbdraw/web/js/app/circular-conservation.js');
const { createCircularTrackSlotEditor } = await import('../../gbdraw/web/js/app/circular-track-slots.js');
const { ringGenomeFileLists } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const {
  circularRingSources, isManagedConservationSlot, orderedConservationSources
} = await import('../../gbdraw/web/js/services/conservation-series.js');
const ringKeys = await import('../../gbdraw/web/js/services/circular-conservation-keys.js');
const { readFileText } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { ringIdentityDefects } = await import('./helpers/ring-identity.mjs');
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
// Each drawing's LOSAT result rows, as the History snapshot service keeps them (D-221): an edit reaches one drawing's.
/** @type {Map<string, readonly object[]>} */
const losatRowsByDrawing = new Map();
/** @type {import('../../gbdraw/web/js/services/losat-cache.js').EditLosatRows} */
const editLosatRows = (drawingId, edit) => { losatRowsByDrawing.set(drawingId, edit(losatRowsByDrawing.get(drawingId) ?? [])); };

const circular = () => modeDrawing(state, 'circular');
sessionConfig.installShownDrawing('circular');
const shownAction = (mode, action) => action(modeDrawing(state, mode));
const askingStackPorts = () => ({
  askChoice: async () => null, runUndoable: (_label, fn) => fn(), busy: () => false, rows: () => [], mode: 'circular'
});
const editor = createCircularTrackSlotEditor({ state });
const sources = createCircularConservationSources({
  state, editLosatRows, askingStackPorts,
  syncCircularConservationSlots: (drawing) => editor.syncCircularConservationSlots(drawing),
  rekeyCircularConservationRows: (drawing, keyChanges) => editor.rekeyCircularConservationRows(drawing, keyChanges),
  circularConservationRowsLeaving: (drawing, entries, keyChanges) => editor.circularConservationRowsLeaving(drawing, entries, keyChanges)
});
const panel = createCircularConservationSeriesPanel({
  state, shownAction, shownDrawingOfMode: (mode) => modeDrawing(state, mode), editLosatRows,
  circularConservationFastaInput: { value: null },
  reconcileCircularConservation: sources.reconcileCircularConservation,
  refreshCircularConservationSlots: sources.refreshCircularConservationSlots,
  askRingSourceRemoval: sources.askRingSourceRemoval
});
const files = createCircularConservationComparisonFiles({
  state, shownAction, editLosatRows, runUndoable: async () => undefined,
  setRingGenomeLists: sources.setRingGenomeLists,
  setCircularConservationSourceFiles: sources.setCircularConservationSourceFiles,
  reconcileCircularConservation: sources.reconcileCircularConservation,
  refreshCircularConservationSlots: sources.refreshCircularConservationSlots,
  askRingSourceRemoval: sources.askRingSourceRemoval
});

/** A text file as the uploader hands it: `arrayBuffer` lets Generate's serializer read it. */
const textFile = (name, text) => ({
  name, type: 'text/plain', size: text.length, lastModified: 0, data: btoa(text),
  arrayBuffer: async () => new TextEncoder().encode(text).buffer
});
const genome = (same, bases) => (same ? textFile('sequence.fasta', `>r\n${bases}\n`) : textFile(`g_${bases}.fasta`, `>r\n${bases}\n`));
const REFERENCE_GENBANK = 'LOCUS       REFERENCE                   8 bp    DNA     linear   UNK 01-JAN-1980\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaccggtt\n//\n';
const addSeq = (file) => files.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [file], value: 'x' } }));
const ringRows = (drawing) => drawing.adv.circular_track_slots.filter((slot) => isManagedConservationSlot(slot));
const bases = async (file) => (await readFileText(state.inputFiles.fileOf(file))).split('\n')[1];

/**
 * A LOSAT Generate of the draft's rings as production runs it: the rings its ring keys read, their tables, the request's
 * rings and the request builder, then the composition root's adoption of the request with Generate's own rings.
 * @param {import('../../gbdraw/web/js/state.js').DrawingState} drawing
 */
const generate = async (drawing) => {
  const searched = ringKeys.searchedRingSources(drawing);
  const before = JSON.stringify(searched.map((ring) => ring.sourceKey));
  const filesData = await sessionConfig.serializeActiveRenderFiles(state, drawing);
  const tables = searched.map((ring, index) => ({
    name: `${ring.fileName.replace(/\.[^.]+$/, '')}.circular_conservation.losatn.tsv`,
    text: `q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t2${index}\n`
  }));
  const canonical = buildCanonicalRenderRequest({
    state, drawing, filesData,
    resolvedCircularConservation: ringKeys.losatRequestRings(searched, tables, ringGenomeFileLists(filesData.c_conservation_fastas))
  });
  // Like production, a Generate with no ring genome searches no ring and adopts no ring result.
  const ringAdoption = searched.length === 0 ? null
    : { searched, builtFromDraft: before === JSON.stringify(ringKeys.searchedRingSources(drawing).map((ring) => ring.sourceKey)) };
  sources.adoptCircularRingResult(drawing, (installRings) => sessionConfig.adoptCanonicalRenderArtifacts(canonical, {
    adoptOwnedRequest: true, drawing, conservationPairs: ringKeys.NO_CONSERVATION_PAIRS, installRings
  }), ringAdoption);
  assert.equal(drawing.files.c_conservation_blasts_source, 'losat-cache', 'the Generate adopted a LOSAT-cache result');
  return canonical;
};
/**
 * A rerender of the committed request (a Rotate): the candidate of the committed request, adopted as
 * `runCommittedCanonicalCandidateInternal` adopts it, through the composition root's adoption (no ring adoption).
 * @param {import('../../gbdraw/web/js/state.js').DrawingState} drawing
 * @param {Parameters<typeof sessionConfig.adoptCanonicalRenderArtifacts>[0]} canonical
 */
const rerender = (drawing, canonical) => sources.adoptCircularRingResult(drawing, (installRings) => sessionConfig.adoptCanonicalRenderArtifacts(
  canonical, { adoptOwnedRequest: true, drawing, conservationPairs: 'committed', installRings }
), null);

/** Everything of the draft's rings a rerender of the committed request must leave as it is. */
const draftRings = (drawing) => JSON.parse(JSON.stringify({
  series: drawing.circularConservation.series,
  labels: drawing.circularConservation.labels,
  fastas: ringGenomeFileLists(drawing.files.c_conservation_fastas).map((list) => list.map((file) => file.fileId)),
  blasts: drawing.files.c_conservation_blasts.map((file) => file.fileId),
  blastsSource: drawing.files.c_conservation_blasts_source,
  rows: ringRows(drawing).map((slot) => [slot.id, slot.width, slot.params.series_key, slot.params.label])
}));

const placeReference = (drawing) => {
  drawing.circularSources.splice(0, drawing.circularSources.length,
    { uid: 'record-1', gb: state.inputFiles.registerFile(textFile('reference.gb', REFERENCE_GENBANK)), gff: null, fasta: null });
  drawing.circularRecordList.value = [{ uid: 'record-1', k: 1, identity: 'record-1', selector: '#1', record_id: 'REFERENCE',
    record_length: 8, detectedTopology: 'linear' }];
  drawing.circularConservation.enabled = true;
  drawing.adv.circular_track_slots_enabled = true;
};
const clearRings = (drawing) => {
  while (drawing.circularConservation.series.length > 0) files.removeCircularConservationSource(0);
  sources.setRingGenomeLists(drawing, [], { setsFile: false, series: 'removed' });
  drawing.circularConservation.labels = '';
  drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length,
    ...drawing.adv.circular_track_slots.filter((slot) => !isManagedConservationSlot(slot)));
};
/** Rings GA, GB, GC of genomes AAAAA, CCCCC, GGGGG with codes 11, 4, 2 and rows 10, 11, 12 px wide. */
const threeRings = async (same) => {
  const drawing = circular();
  clearRings(drawing);
  placeReference(drawing);
  for (const ring of ['AAAAA', 'CCCCC', 'GGGGG']) await addSeq(genome(same, ring));
  drawing.circularConservation.series.forEach((entry, index) => Object.assign(entry, {
    label: `G${'ABC'[index]}`, color: ['#110000', '#001100', '#000011'][index], losat_gencode: [11, 4, 2][index]
  }));
  sources.refreshCircularConservationSlots(drawing);
  ringRows(drawing).forEach((slot, index) => { slot.width = `${10 + index}px`; });
  return drawing;
};

/** Each kind of ring edit after a Generate, as the review's q4 probe made them. */
const EDITS = {
  'Remove series 1': () => files.removeCircularConservationSource(1),
  'Remove series 0': () => files.removeCircularConservationSource(0),
  'Move series 0 down': () => panel.moveCircularConservationSeries(0, 1),
  'a label and colour': (drawing) => {
    Object.assign(drawing.circularConservation.series[0], { label: 'EDITED', color: '#abcdef' });
    panel.refreshCircularConservationSlots();
  },
  'a genetic code': (drawing) => { drawing.circularConservation.series[2].losat_gencode = 25; },
  'Add Seq': () => addSeq(textFile('new.fasta', '>r\nTTTTT\n')),
  'a ring row\'s Add file': () => files.addRingRowFiles(1, /** @type {any} */ ({ target: { files: [textFile('plasmid.fasta', '>p\nACACA\n')], value: 'x' } })),
  'a stack row\'s width': (drawing) => { ringRows(drawing)[2].width = '30px'; }
};

test('a rerender of the committed request leaves the draft\'s rings as each ring edit left them (#1777 review 7, OV-3194)', async () => {
  /** Each edit after which the rerender changed the draft's rings or broke their identity. */
  const changed = [];
  for (const same of [false, true]) {
    for (const [name, edit] of Object.entries(EDITS)) {
      const drawing = await threeRings(same);
      const canonical = await generate(drawing);
      await edit(drawing);
      const edited = draftRings(drawing);
      rerender(drawing, canonical);
      const where = `${name} (${same ? 'sequence.fasta x3' : 'distinct files'})`;
      if (!isDeepStrictEqual(draftRings(drawing), edited)) changed.push(`${where}: the draft's rings changed`);
      if (ringIdentityDefects(drawing).length > 0) changed.push(`${where}: ${ringIdentityDefects(drawing).join('; ')}`);
      clearRings(drawing);
    }
  }
  assert.deepEqual(changed, [], 'the draft\'s rings after a rerender of the committed request');
});

test('a Generate after a replay Load that lost ring 0\'s genome gives each ring its own code, label and row (#1777 review 7, OV-3195)', async () => {
  let drawing = await threeRings(false);
  await generate(drawing);
  const rows = Object.fromEntries(ringRows(drawing).map((slot) => [slot.params.label, [slot.id, slot.width]]));
  const result = await exportSession('ring loss', { askChoice: async () => 'download' });
  const saved = await readSavedSession(result.blob);
  const savedDrawing = saved.drawings.find((candidate) => candidate.mode === 'circular');
  savedDrawing.webFiles.bindings.c_conservation_fastas[0].forEach((file, index) => { file.resourceId = `missing-0-${index}`; });
  const loaded = await importSession(await sessionLoadEvent(saved),
    { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok');
  drawing = circular();
  sources.reconcileLoadedCircularConservation(drawing, true);
  /** The rings as drawn and as the next Generate searches them: label, genome bases, genetic code, stack row. */
  const rings = async () => {
    const drawn = await Promise.all(orderedConservationSources(circularRingSources(drawing), drawing.circularConservation)
      .map(async (ring) => {
        const row = ringRows(drawing).find((slot) => slot.params.series_key === ring.sourceKey);
        return [ring.label, ring.losat_gencode, row ? [row.id, row.width] : null];
      }));
    const searched = await Promise.all(ringKeys.searchedRingSources(drawing)
      .map(async (ring) => [ring.label, (await Promise.all(ring.files.map(bases))).join('+'), ring.losat_gencode]));
    return { drawn, searched };
  };
  assert.deepEqual((await rings()).searched, [['GB', 'CCCCC', 4], ['GC', 'GGGGG', 2]], 'the Load lost ring 0\'s genome');
  await generate(drawing);
  const afterGenerate = {
    drawn: [['GB', 4, rows.GB], ['GC', 2, rows.GC]],
    searched: [['GB', 'CCCCC', 4], ['GC', 'GGGGG', 2]]
  };
  assert.deepEqual(await rings(), afterGenerate, 'the Generate after the Load');
  assert.equal(ringGenomeFileLists(drawing.files.c_conservation_fastas).length, drawing.files.c_conservation_blasts.length,
    'one genome list per table');
  assert.deepEqual(ringIdentityDefects(drawing), []);
  await generate(drawing);
  assert.deepEqual(await rings(), afterGenerate, 'the next Generate');
  await addSeq(textFile('g_TTTTT.fasta', '>r\nTTTTT\n'));
  const added = await rings();
  assert.deepEqual(added.drawn.slice(0, 2), afterGenerate.drawn, 'Add Seq keeps the rings before it');
  assert.deepEqual(added.searched, [...afterGenerate.searched, ['g_TTTTT', 'TTTTT', 1]], 'then Add Seq');
  assert.deepEqual(ringIdentityDefects(drawing), []);
  clearRings(drawing);
});

test('a Generate after a replay Load that lost every ring\'s genome leaves each ring\'s code, label and row (#1777 review 8, OV-3199)', async () => {
  let drawing = await threeRings(false);
  await generate(drawing);
  const result = await exportSession('every ring lost', { askChoice: async () => 'download' });
  const saved = await readSavedSession(result.blob);
  saved.drawings.find((candidate) => candidate.mode === 'circular').webFiles.bindings.c_conservation_fastas
    .forEach((list, place) => list.forEach((file, index) => { file.resourceId = `missing-${place}-${index}`; }));
  const loaded = await importSession(await sessionLoadEvent(saved),
    { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok');
  drawing = circular();
  sources.reconcileLoadedCircularConservation(drawing, true);
  assert.deepEqual(ringKeys.searchedRingSources(drawing), [], 'the Load lost every ring\'s genome');
  const view = () => ({
    series: drawing.circularConservation.series.map((entry) => [entry.label, entry.losat_gencode]),
    rows: ringRows(drawing).map((slot) => [slot.id, slot.width, slot.params.label])
  });
  const loadedRings = view();
  assert.deepEqual(loadedRings.series, [['GA', 11], ['GB', 4], ['GC', 2]]);
  assert.deepEqual(loadedRings.rows.map(([, width]) => width), ['10px', '11px', '12px']);
  await generate(drawing);
  assert.deepEqual(view(), loadedRings, 'a Generate that searched no ring');
  assert.deepEqual(ringIdentityDefects(drawing), []);
  clearRings(drawing);
});

test('a ring label with a comma names no later ring (#1777 review 7, OV-3193)', async () => {
  const drawing = circular();
  clearRings(drawing);
  placeReference(drawing);
  for (const ring of ['AAAAA', 'CCCCC']) await addSeq(genome(false, ring));
  drawing.circularConservation.series[0].label = 'E. coli, K-12';
  await addSeq(genome(false, 'GGGGG'));
  await addSeq(genome(false, 'TTTTT'));
  assert.deepEqual(drawing.circularConservation.series.map((entry) => entry.label), ['E. coli, K-12', 'g_CCCCC', 'g_GGGGG', 'g_TTTTT']);
  files.removeCircularConservationSource(1);
  await addSeq(genome(false, 'ACACA'));
  assert.deepEqual(drawing.circularConservation.series.map((entry) => entry.label), ['E. coli, K-12', 'g_GGGGG', 'g_TTTTT', 'g_ACACA']);
  clearRings(drawing);
});

test('a replay Load folds an older Session\'s ring label list into the series, so it names no later ring (#1777 review 8, OV-3200)', async () => {
  let drawing = await threeRings(false);
  drawing.circularConservation.series[0].label = 'E. coli, K-12';
  await generate(drawing);
  const result = await exportSession('older list', { askChoice: async () => 'download' });
  const saved = await readSavedSession(result.blob);
  const savedRings = saved.drawings.find((candidate) => candidate.mode === 'circular').config.circularConservation;
  // dev and main write the series labels' comma join at a LOSAT Generate's adoption.
  savedRings.labels = savedRings.series.map((entry) => entry.label).join(',');
  const loaded = await importSession(await sessionLoadEvent(saved),
    { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok');
  drawing = circular();
  sources.reconcileLoadedCircularConservation(drawing, true);
  assert.equal(drawing.files.c_conservation_blasts_source, 'losat-cache', 'the Load replays the LOSAT tables');
  assert.equal(drawing.circularConservation.labels, '', 'the Load read the list into the series');
  await addSeq(textFile('g_TTTTT.fasta', '>r\nTTTTT\n'));
  assert.deepEqual(drawing.circularConservation.series.map((entry) => entry.label), ['E. coli, K-12', 'GB', 'GC', 'g_TTTTT']);
  assert.deepEqual(ringIdentityDefects(drawing), []);
  clearRings(drawing);
});

test('a ring edit clears only its drawing\'s LOSAT cache-info rows (#1777 review 7, OV-3196)', async () => {
  const drawing = circular();
  const linear = modeDrawing(state, 'linear');
  clearRings(drawing);
  placeReference(drawing);
  const ringRow = { key: 'ring-key', filename: 'g_AAAAA.circular_conservation.losatn.tsv', display: true };
  const linearRow = { key: 'pair-key', filename: 'pair.tsv', display: true, edgeKey: 'edge-1' };
  const linearRows = [linearRow];
  for (const [name, edit] of /** @type {const} */ ([
    ['Add Seq', () => addSeq(genome(false, 'AAAAA'))],
    ['a ring row\'s Add file', () => files.addRingRowFiles(0, /** @type {any} */ ({ target: { files: [textFile('plasmid.fasta', '>p\nACACA\n')], value: 'x' } }))],
    ['Remove series', () => files.removeCircularConservationSource(0)]
  ])) {
    losatRowsByDrawing.set(drawing.id, [ringRow]);
    losatRowsByDrawing.set(linear.id, linearRows);
    await edit();
    assert.deepEqual(losatRowsByDrawing.get(drawing.id), [], `${name}: the ring drawing's rows go`);
    assert.equal(losatRowsByDrawing.get(linear.id), linearRows, `${name}: the Linear drawing's rows stay`);
  }
  losatRowsByDrawing.clear();
  clearRings(drawing);
});
