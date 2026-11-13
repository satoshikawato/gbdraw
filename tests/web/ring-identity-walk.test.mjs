// Guard: the rings' identity invariant through every ring transition the Web offers (#1777 round 8, review 7). A
// deterministic walk (fixed seeds, bounded steps) over Add Seq, Remove series, a ring row's Add file and file
// remove, Move series, a ring's label (some with commas) and colour edit, a ring's LOSAT genetic code edit, a stack
// row's width edit, a Generate as production builds it (the rings its ring keys read, `losatRequestRings`, the request
// builder, the adoption with the rings it searched), a rerender of the committed request after any of these steps,
// committed or not (R1 (c): the draft stays as it is), the replacement of a replayed table, the source radio's round
// trip, a Session 51 Save -> Load (also one that lost a ring's genome file, in a LOSAT-cache replay too), and the track
// stack off and on, over rings of NCBI `sequence.fasta` files of one name, size and date mixed with rings of distinct
// files. After every step it checks the invariant (`ringIdentityDefects`: each series entry names one source by key,
// file name and place; a replay's tables sit beside their genome lists; each drawn ring's one stack row is its own row
// of the walk's model) and that each ring keeps its label, its colour, its genetic code (also as the next Generate
// reads the rings) and its stack row's id and width. Every expected value comes from the walk's own model. Among its
// steps it must draw a Remove of one of two copies and an Add Seq while a ring result is replayed and its rows hold
// edits (OV-3306, OV-3307).
// Undo and Redo are not walked: a History restore installs the series, the stack rows and the file bindings it
// captured together (R11), after a step this walk checks, and calls no reconcile or re-key.
import assert from 'node:assert/strict';
import { test } from 'node:test';
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
const {
  bindDrawingInputFileList, bindDrawingInputFileTree, ringGenomeFileLists, setConservationBlastsSource
} = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const {
  circularRingSources, isManagedConservationSlot, orderedConservationSources, ringDrawnPositions, ringGenomeFilesInRingOrder
} = await import('../../gbdraw/web/js/services/conservation-series.js');
const { NO_CONSERVATION_PAIRS, losatRequestRings, searchedRingSources } = await import('../../gbdraw/web/js/services/circular-conservation-keys.js');
const { readFileText } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { ringIdentityDefects } = await import('./helpers/ring-identity.mjs');
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { exportSession, importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const circular = () => modeDrawing(state, 'circular');
sessionConfig.installShownDrawing('circular');
const shownAction = (mode, action) => action(modeDrawing(state, mode));
const askingStackPorts = () => ({
  askChoice: async () => null, runUndoable: (_label, fn) => fn(), busy: () => false, rows: () => [], mode: 'circular'
});
const editor = createCircularTrackSlotEditor({ state });
const sources = createCircularConservationSources({
  state, editLosatRows: () => {}, askingStackPorts,
  syncCircularConservationSlots: (drawing) => editor.syncCircularConservationSlots(drawing),
  rekeyCircularConservationRows: (drawing, keyChanges) => editor.rekeyCircularConservationRows(drawing, keyChanges),
  circularConservationRowsLeaving: (drawing, entries, keyChanges) => editor.circularConservationRowsLeaving(drawing, entries, keyChanges)
});
const panel = createCircularConservationSeriesPanel({
  state, shownAction, shownDrawingOfMode: (mode) => modeDrawing(state, mode), editLosatRows: () => {},
  circularConservationFastaInput: { value: null },
  reconcileCircularConservation: sources.reconcileCircularConservation,
  refreshCircularConservationSlots: sources.refreshCircularConservationSlots,
  askRingSourceRemoval: sources.askRingSourceRemoval
});
const files = createCircularConservationComparisonFiles({
  state, shownAction, editLosatRows: () => {}, runUndoable: async () => undefined,
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
/** Each walk file's own bases, so a file names the ring that holds it however Generate rebinds it. */
let fileCount = 0;
const nextBases = () => {
  fileCount += 1;
  return Array.from({ length: 5 }, (_, digit) => 'ACGT'[Math.floor(fileCount / 4 ** digit) % 4]).join('');
};
/** NCBI's name, one size and one date (`sequence.fasta`, 9 B, 0), or a file of its own name and size. */
const genomeFile = (same) => {
  const bases = nextBases();
  return same ? textFile('sequence.fasta', `>r\n${bases}\n`) : textFile(`genome_${bases}.fasta`, `>genome ${bases}\n${bases}${bases}\n`);
};
const basesOfText = (text) => text.split('\n')[1].slice(0, 5);
const REFERENCE_GENBANK = 'LOCUS       REFERENCE                   8 bp    DNA     linear   UNK 01-JAN-1980\nFEATURES             Location/Qualifiers\nORIGIN\n        1 aaccggtt\n//\n';

/** A seeded generator (mulberry32), so every run walks the same steps. */
const seeded = (seed) => () => {
  seed = (seed + 0x6d2b79f5) | 0;
  let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
};

/**
 * The walk's model: each ring by its own id, with its files' bases, its label, own colour (or null) and scheme index
 * (PD-OI-110) and LOSAT genetic code, its
 * stack row's id and width once the stack has made it, and whether a Load lost its genome while its replayed table
 * stays (`lost`: drawn from its table, not searched by the next Generate).
 * @typedef {{ id: number, bases: Set<string>, label: string, color: string | null, schemeIndex: number, code: number,
 *   row: { id: string, width: unknown } | null, lost: boolean }} ModelRing
 */

const ringRows = (drawing) => drawing.adv.circular_track_slots.filter((slot) => isManagedConservationSlot(slot));

/** Everything of the draft's rings a rerender of the committed request leaves as it is (R1 (c)). */
const draftRings = (drawing) => JSON.parse(JSON.stringify({
  series: drawing.circularConservation.series,
  labels: drawing.circularConservation.labels,
  enabled: drawing.circularConservation.enabled,
  fastas: ringGenomeFileLists(drawing.files.c_conservation_fastas).map((list) => list.map((file) => file.fileId)),
  blasts: drawing.files.c_conservation_blasts.map((file) => file.fileId),
  blastsSource: drawing.files.c_conservation_blasts_source,
  rows: ringRows(drawing).map((slot) => [slot.id, slot.width, slot.params.series_key])
}));

/** The drawing's rings as drawn, in order: each with its genome's bases, label, colour, genetic code and stack row. */
const drawnRings = async (drawing) => {
  const genomes = ringGenomeFilesInRingOrder(drawing);
  const drawn = orderedConservationSources(circularRingSources(drawing), drawing.circularConservation);
  return Promise.all(drawn.map(async (source, index) => ({
    bases: await Promise.all((genomes[index] || []).map(async (file) => basesOfText(await readFileText(state.inputFiles.fileOf(file))))),
    label: source.label,
    color: source.ownColor,
    schemeIndex: source.schemeIndex,
    code: source.losat_gencode,
    row: ringRows(drawing).find((slot) => String(slot.params.series_key) === source.sourceKey) || null
  })));
};

const runWalk = async (seed, steps) => {
  const random = seeded(seed);
  const pick = (count) => Math.floor(random() * count);
  // A Session Load installs new drawings: the walk reads the shown Circular drawing again after it.
  let drawing = circular();
  /** @type {ModelRing[]} */
  let model = [];
  let ringCount = 0;
  /** @type {string[]} */
  const trail = [];
  /** @type {Map<string, number>} */
  const kinds = new Map();
  /**
   * The committed request (the last Generate's, or the one a Load installed); a rerender of it (a Rotate) replaces only
   * the committed Session, whatever the draft's rings did since.
   * @type {Parameters<typeof sessionConfig.adoptCanonicalRenderArtifacts>[0] | null}
   */
  let committed = null;
  drawing.circularSources.splice(0, drawing.circularSources.length,
    { uid: 'record-1', gb: state.inputFiles.registerFile(textFile('reference.gb', REFERENCE_GENBANK)), gff: null, fasta: null });
  drawing.circularRecordList.value = [{ uid: 'record-1', k: 1, identity: 'record-1', selector: '#1', record_id: 'REFERENCE',
    record_length: 8, detectedTopology: 'linear' }];
  drawing.circularConservation.enabled = true;
  drawing.adv.circular_track_slots_enabled = true;

  /** The drawn rings in the model's order, each found by the bases of its genome files. */
  const check = async (step) => {
    const where = `seed ${seed}, step ${step}: ${trail.join(' > ')}`;
    assert.deepEqual(ringIdentityDefects(drawing, { rows: model.map((ring) => (ring.row ? ring.row.id : null)) }), [],
      `the rings' identity invariant (${where})`);
    const drawn = await drawnRings(drawing);
    assert.equal(drawn.length, model.length, `one drawn ring per ring (${where})`);
    // PD-OI-110, fix round 13: no two rings share a scheme index, so no two draw one scheme colour.
    const schemeIndexes = drawn.map((ring) => ring.schemeIndex);
    assert.equal(new Set(schemeIndexes).size, schemeIndexes.length, `no two rings share a scheme index ${JSON.stringify(schemeIndexes)} (${where})`);
    drawn.forEach((ring, index) => {
      const expected = model[index];
      assert.ok(expected.lost
        ? ring.bases.length === 0
        : ring.bases.length > 0 && ring.bases.every((bases) => expected.bases.has(bases)) && ring.bases.length === expected.bases.size,
      `ring ${index} draws ring ${expected.id}'s files (${where})`);
      assert.deepEqual([ring.label, ring.color, ring.schemeIndex, ring.code], [expected.label, expected.color, expected.schemeIndex, expected.code],
        `ring ${expected.id}'s label, colour and genetic code (${where})`);
      if (drawing.adv.circular_track_slots_enabled !== true) return;
      assert.ok(ring.row, `ring ${expected.id} has its stack row (${where})`);
      if (expected.row) {
        assert.deepEqual([ring.row.id, ring.row.width], [expected.row.id, expected.row.width], `ring ${expected.id}'s stack row (${where})`);
      } else {
        expected.row = { id: String(ring.row.id), width: ring.row.width };
      }
    });
    // The rings the next Generate searches, in order, and the genetic code it searches each with, as its ring keys
    // read them (`circularConservationKeys`).
    const searchedRings = await Promise.all(searchedRingSources(drawing).map(async (source) => ({
      bases: await Promise.all(source.files.map(async (file) => basesOfText(await readFileText(state.inputFiles.fileOf(file))))),
      code: source.losat_gencode
    })));
    const searchable = model.filter((ring) => !ring.lost);
    assert.equal(searchedRings.length, searchable.length, `the next Generate searches every ring with a genome (${where})`);
    searchedRings.forEach((ring, index) => {
      const expected = searchable[index];
      assert.ok(ring.bases.length === expected.bases.size && ring.bases.every((bases) => expected.bases.has(bases)),
        `the next Generate searches ring ${expected.id}'s files at ${index} (${where})`);
      assert.equal(ring.code, expected.code, `the next Generate searches ring ${expected.id} with its genetic code (${where})`);
    });
  };
  const replayed = () => drawing.files.c_conservation_blasts_source === 'losat-cache' && drawing.files.c_conservation_blasts.length > 0;
  /** Counts a kind of step the walk must draw. */
  const drew = (kind) => kinds.set(kind, (kinds.get(kind) || 0) + 1);
  /** A replay's ring rows with an edit to keep: the rows a ring transition there re-keys (OV-3306, OV-3307). */
  const replayRows = () => replayed() && drawing.adv.circular_track_slots_enabled === true && model.some((ring) => ring.row);
  const newRingOf = (index, bases) => {
    const entry = drawing.circularConservation.series[index];
    ringCount += 1;
    return { id: ringCount, bases: new Set([bases]), label: entry.label, color: entry.color, schemeIndex: entry.schemeIndex, code: entry.losat_gencode, row: null, lost: false };
  };
  /** A change of the genome files leaves a LOSAT-cache replay: a ring whose genome a Load lost goes with its table. */
  const leaveReplay = () => { model = model.filter((ring) => !ring.lost); };

  /** @type {Record<string, { weight: number, can: () => boolean, run: () => Promise<void> | void }>} */
  const STEPS = {
    'add seq': {
      weight: 4, can: () => model.length < 6,
      run: async () => {
        const file = genomeFile(random() < 0.6);
        const bases = basesOfText(atob(file.data));
        // OV-3307: each kept ring's key moves from its table to its genome files; its row keeps its id and width.
        if (replayRows()) drew('add seq to a replay with ring rows');
        await files.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [file], value: 'x' } }));
        const index = drawing.circularConservation.series.length - 1;
        assert.equal(drawing.circularConservation.series[index].label, file.name.replace(/\.[^.]+$/, ''),
          `a new ring takes its file's name (seed ${seed}: ${trail.join(' > ')})`);
        leaveReplay();
        model.push(newRingOf(index, bases));
      }
    },
    'remove series': {
      weight: 2, can: () => model.length > 0,
      run: () => {
        const index = pick(model.length);
        // OV-3306: copies of one file share a table name; the kept copy's key moves to the freed copy number, and its
        // row stays its own.
        const names = drawing.circularConservation.series.map(({ fileName }) => fileName);
        if (replayRows() && names.filter((name) => name === names[index]).length > 1) drew('remove a copy from a replay with ring rows');
        files.removeCircularConservationSource(index);
        model = model.filter((_, ring) => ring !== index);
      }
    },
    'row add file': {
      weight: 2, can: () => model.length > 0,
      run: () => {
        const index = pick(model.length);
        const file = random() < 0.6 ? textFile('plasmid.fasta', `>p\n${nextBases()}\n`) : genomeFile(false);
        files.addRingRowFiles(index, /** @type {any} */ ({ target: { files: [file], value: 'x' } }));
        const ring = model[index];
        // A ring whose genome a Load lost takes the added file as its genome.
        if (ring.lost) Object.assign(ring, { bases: new Set(), lost: false });
        ring.bases.add(basesOfText(atob(file.data)));
        leaveReplay();
      }
    },
    'row remove file': {
      weight: 2, can: () => model.some((ring) => !ring.lost && ring.bases.size > 1),
      run: async () => {
        const candidates = model.map((ring, index) => (!ring.lost && ring.bases.size > 1 ? index : -1)).filter((index) => index >= 0);
        const index = candidates[pick(candidates.length)];
        const genome = ringGenomeFilesInRingOrder(drawing)[index];
        const fileIndex = pick(genome.length);
        const bases = basesOfText(await readFileText(state.inputFiles.fileOf(genome[fileIndex])));
        files.removeRingRowFile(index, fileIndex);
        model[index].bases.delete(bases);
        leaveReplay();
      }
    },
    'move series': {
      weight: 2, can: () => model.length > 1,
      run: () => {
        const index = pick(model.length - 1);
        panel.moveCircularConservationSeries(index, 1);
        model.splice(index, 2, model[index + 1], model[index]);
      }
    },
    'edit label and colour': {
      weight: 1, can: () => model.length > 0,
      run: () => {
        const index = pick(model.length);
        const entry = drawing.circularConservation.series[index];
        // Half of the labels hold a comma, which no stored list may split into another ring's name (OV-3193).
        entry.label = random() < 0.5 ? `ring ${model[index].id}, strain K-${model[index].id}` : `ring ${model[index].id}`;
        entry.color = `#00${String(model[index].id).padStart(2, '0')}aa`;
        panel.refreshCircularConservationSlots();
        Object.assign(model[index], { label: entry.label, color: entry.color });
      }
    },
    'set a genetic code': {
      // The ring's LOSAT genetic code field (`circularConservation.series[i].losat_gencode`). A rerender of the
      // committed rings keeps it with its ring (OV-3178); the request does not carry it.
      weight: 2, can: () => model.length > 0,
      run: () => {
        const index = pick(model.length);
        const code = 2 + pick(24);
        drawing.circularConservation.series[index].losat_gencode = code;
        model[index].code = code;
      }
    },
    'widen a stack row': {
      weight: 2, can: () => model.some((ring) => ring.row) && drawing.adv.circular_track_slots_enabled === true,
      run: () => {
        const withRows = model.filter((candidate) => candidate.row);
        const ring = withRows[pick(withRows.length)];
        const row = ringRows(drawing).find((slot) => slot.id === ring.row.id);
        row.width = `${10 + ring.id + pick(5)}px`;
        ring.row.width = row.width;
      }
    },
    'generate': {
      // Also when a Load lost every ring's genome: such a Generate searches no ring (OV-3199).
      weight: 3, can: () => model.length > 0,
      run: async () => {
        // As production runs it: the rings its ring keys read, their tables, the request's rings, the request builder,
        // and the adoption with the rings it searched.
        const searched = searchedRingSources(drawing);
        const filesData = await sessionConfig.serializeActiveRenderFiles(state, drawing);
        const tables = searched.map((ring, index) => ({
          name: `${ring.fileName.replace(/\.[^.]+$/, '')}.circular_conservation.losatn.tsv`,
          text: `q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t2${index}\n`
        }));
        const canonical = buildCanonicalRenderRequest({
          state, drawing, filesData,
          resolvedCircularConservation: losatRequestRings(searched, tables, ringGenomeFileLists(filesData.c_conservation_fastas))
        });
        // Production searches only when a ring has genome files (run-analysis.js), and hands the adoption the rings it
        // searched; a Generate that searched none hands none, and the draft's rings stay as they are (OV-3199).
        const ringAdoption = searched.length > 0 ? { searched, builtFromDraft: true } : null;
        const before = draftRings(drawing);
        sources.adoptCircularRingResult(drawing, (installRings) => sessionConfig.adoptCanonicalRenderArtifacts(canonical, {
          adoptOwnedRequest: true, drawing, conservationPairs: NO_CONSERVATION_PAIRS, installRings
        }), ringAdoption);
        if (ringAdoption) {
          assert.ok(replayed(), 'the Generate adopted its LOSAT-cache result');
          // The rings it did not search (a genome a Load lost) go: the Result draws the rings it searched.
          model = model.filter((ring) => !ring.lost);
        } else {
          assert.deepEqual(draftRings(drawing), before, `the draft's rings after a Generate that searched none (seed ${seed}: ${trail.join(' > ')})`);
          kinds.set('a Generate that searched no ring', (kinds.get('a Generate that searched no ring') || 0) + 1);
        }
        committed = canonical;
      }
    },
    'rerender the committed request': {
      // A rerender of the committed request (a Rotate, an Align, a placement) after any step, committed or not: it
      // replaces only the committed Session and leaves the draft's rings as they are (R1 (c), OV-3194).
      weight: 3, can: () => committed !== null,
      run: () => {
        const before = draftRings(drawing);
        sources.adoptCircularRingResult(drawing, (installRings) => sessionConfig.adoptCanonicalRenderArtifacts(
          /** @type {NonNullable<typeof committed>} */ (committed), { adoptOwnedRequest: true, drawing, conservationPairs: 'committed', installRings }
        ), null);
        assert.deepEqual(draftRings(drawing), before, `the draft's rings after a rerender (seed ${seed}: ${trail.join(' > ')})`);
      }
    },
    'replace a replayed table': {
      weight: 1, can: () => replayed(),
      run: () => {
        // The replacement of an input file that the drawing uses (`createInputReplacement`): the table at one place.
        const tables = [...drawing.files.c_conservation_blasts];
        const place = pick(tables.length);
        tables[place] = textFile(tables[place].name, `q\ts\t99\t8\t0\t0\t1\t8\t1\t8\t1e-9\t3${place}\n`);
        sources.replaceCircularConservationBlasts(drawing, tables);
      }
    },
    'source radio round trip': {
      weight: 1, can: () => replayed(),
      run: () => {
        panel.changeCircularConservationSource('upload');
        panel.changeCircularConservationSource('losat');
      }
    },
    'save and load': {
      weight: 2, can: () => true,
      run: async () => {
        // Save Session, then Load of what it wrote, then the Load's reconcile (`afterImport`). One of three Loads of
        // LOSAT rings could not read one ring's genome files: outside a replay that ring goes; in a LOSAT-cache replay
        // its table stays, so it is drawn, and the next Generate does not search it.
        const result = await exportSession('ring walk', { askChoice: async () => 'download' });
        const saved = await readSavedSession(result.blob);
        const withGenome = model.map((ring, index) => (ring.lost ? -1 : index)).filter((index) => index >= 0);
        // In a replay, half the Loads that lose a genome loses every ring's genome, so the next Generate
        // searches no ring (OV-3199).
        const roll = withGenome.length > 0 ? random() : 1;
        const lostRings = roll < 0.17 && replayed() ? withGenome : roll < 0.34 ? [withGenome[pick(withGenome.length)]] : [];
        const savedDrawing = saved.drawings.find((candidate) => candidate.mode === 'circular');
        for (const lost of lostRings) {
          const place = drawing.circularConservation.series[lost].sourceIndex;
          savedDrawing.webFiles.bindings.c_conservation_fastas[place].forEach((file, index) => { file.resourceId = `missing-${place}-${index}`; });
          if (replayed()) model[lost].lost = true;
        }
        if (lostRings.length > 0) {
          if (replayed()) {
            kinds.set('a replay Load that lost a genome', (kinds.get('a replay Load that lost a genome') || 0) + 1);
          } else {
            model = model.filter((_, ring) => !lostRings.includes(ring));
          }
          trail[trail.length - 1] += ` (lost ring ${lostRings.join(', ')})`;
          kinds.set('a Load that lost a genome', (kinds.get('a Load that lost a genome') || 0) + 1);
        }
        const before = drawing;
        const loaded = await importSession(await sessionLoadEvent(saved),
          { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
        assert.equal(loaded.status, 'ok', `the Load (seed ${seed}: ${trail.join(' > ')})`);
        drawing = circular();
        assert.notEqual(drawing, before, 'the Load installed new drawings');
        sources.reconcileLoadedCircularConservation(drawing, true);
        committed = sessionConfig.getCommittedCanonicalSession();
      }
    },
    'stack off and on': {
      weight: 2, can: () => true,
      run: () => {
        editor.setCircularTrackSlotsEnabled(!(drawing.adv.circular_track_slots_enabled === true));
      }
    }
  };
  // PD-OI-110 (fix round 13): the ring colour steps draw from their own generator after a step, so the walk above
  // draws the same steps as without them.
  const colourRandom = seeded(seed + 110);
  const colourPick = (count) => Math.floor(colourRandom() * count);
  /** @type {Record<string, { can: () => boolean, run: () => void }>} */
  const COLOUR_STEPS = {
    're-assign colours': {
      // 「今の並びで色を振り直す」 under 10 colors: every ring takes its drawn position as its scheme index, one with
      // its own colour too (it keeps the colour).
      can: () => model.length > 1,
      run: () => {
        panel.setCircularConservationColorScheme('ten');
        panel.reassignCircularConservationColors();
        const slots = drawing.adv.circular_track_slots_enabled === true ? drawing.adv.circular_track_slots : null;
        const positions = ringDrawnPositions(model.length, slots);
        model.forEach((ring, index) => { ring.schemeIndex = positions[index]; });
      }
    },
    'use the scheme colour': {
      can: () => model.some((ring) => ring.color !== null),
      run: () => {
        const own = model.map((ring, index) => (ring.color !== null ? index : -1)).filter((index) => index >= 0);
        const index = own[colourPick(own.length)];
        panel.useCircularConservationSchemeColor(index);
        model[index].color = null;
      }
    },
  };
  const names = Object.keys(STEPS);
  for (let step = 1; step <= steps; step += 1) {
    const open = names.filter((name) => STEPS[name].can());
    const total = open.reduce((sum, name) => sum + STEPS[name].weight, 0);
    let roll = random() * total;
    const name = open.find((candidate) => (roll -= STEPS[candidate].weight) < 0) || open[open.length - 1];
    trail.push(name);
    kinds.set(name, (kinds.get(name) || 0) + 1);
    await STEPS[name].run();
    await check(step);
    const colourOpen = Object.keys(COLOUR_STEPS).filter((candidate) => COLOUR_STEPS[candidate].can());
    if (colourOpen.length === 0 || colourRandom() >= 0.2) continue;
    const colourStep = colourOpen[colourPick(colourOpen.length)];
    trail.push(colourStep);
    kinds.set(colourStep, (kinds.get(colourStep) || 0) + 1);
    COLOUR_STEPS[colourStep].run();
    await check(step);
  }
  // The walk ends with the stack on, so the rows of every ring are checked once more.
  if (drawing.adv.circular_track_slots_enabled !== true) {
    editor.setCircularTrackSlotsEnabled(true);
    await check(steps + 1);
  }
  while (drawing.circularConservation.series.length > 0) files.removeCircularConservationSource(0);
  bindDrawingInputFileTree(state.inputFiles, drawing.slices.inputs, 'c_conservation_fastas', []);
  drawing.adv.circular_track_slots.splice(0, drawing.adv.circular_track_slots.length,
    ...drawing.adv.circular_track_slots.filter((slot) => !isManagedConservationSlot(slot)));
  return kinds;
};

test('each ring keeps one identity, its label, colour, genetic code and stack row through a walk of ring transitions', async () => {
  /** @type {Map<string, number>} */
  const drawnKinds = new Map();
  for (const seed of [1777, 3165, 3166, 3167, 8, 20261110]) {
    for (const [kind, count] of await runWalk(seed, 80)) drawnKinds.set(kind, (drawnKinds.get(kind) || 0) + count);
  }
  console.log(`ring walk steps: ${JSON.stringify(Object.fromEntries(drawnKinds))}`);
  for (const kind of ['add seq', 'remove series', 'row add file', 'row remove file', 'move series', 'edit label and colour',
    're-assign colours', 'use the scheme colour', 'set a genetic code', 'widen a stack row', 'generate', 'rerender the committed request', 'replace a replayed table', 'source radio round trip',
    'save and load', 'a Load that lost a genome', 'a replay Load that lost a genome', 'a Generate that searched no ring',
    'add seq to a replay with ring rows', 'remove a copy from a replay with ring rows', 'stack off and on']) {
    assert.ok((drawnKinds.get(kind) || 0) >= 3, `the walk draws ${kind} (${drawnKinds.get(kind) || 0} times)`);
  }
});
