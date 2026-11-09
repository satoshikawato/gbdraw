// Guard: the rings' identity invariant through every ring transition the Web offers (#1777 round 8). A
// deterministic walk (fixed seeds, bounded steps) over Add Seq, Remove series, a ring row's Add file and file
// remove, Move series, a ring's label and colour edit, a ring's LOSAT genetic code edit, a stack row's width edit, Generate's adoption of a LOSAT result
// (the real request and adoption path), the replacement of a replayed table, the source radio's round trip, a Session 51 Save -> Load (also one that lost
// a ring's genome file), and the track stack off and on, over rings of NCBI `sequence.fasta` files of one name, size
// and date mixed with rings of distinct files. After every step it checks the invariant (`ringIdentityDefects`: each
// series entry names one source by key, file name and place; each managed stack row carries the key the one matcher
// gives its entry) and that each ring keeps its label, its colour, its genetic code (also as the next Generate reads
// the rings, round 9) and its stack row's id and width.
// Undo and Redo are not walked: a History restore installs the series, the stack rows and the file bindings it
// captured together (R11), after a step this walk checks, and calls no reconcile or re-key.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';
import { modeDrawing } from './helpers/drawing-state.mjs';
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
  circularRingSources, isManagedConservationSlot, orderedConservationSources, ringGenomeFilesInRingOrder
} = await import('../../gbdraw/web/js/services/conservation-series.js');
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
  state, writeLosatCacheInfo: () => {}, askingStackPorts,
  syncCircularConservationSlots: (drawing) => editor.syncCircularConservationSlots(drawing),
  rekeyCircularConservationRows: (drawing, keyChanges) => editor.rekeyCircularConservationRows(drawing, keyChanges),
  circularConservationRowsLeaving: (drawing, entries, keyChanges) => editor.circularConservationRowsLeaving(drawing, entries, keyChanges),
  keepCircularConservationRowsThrough: (drawing, replaceSources) => editor.keepCircularConservationRowsThrough(drawing, replaceSources)
});
const panel = createCircularConservationSeriesPanel({
  state, shownAction, shownDrawingOfMode: (mode) => modeDrawing(state, mode), writeLosatCacheInfo: () => {},
  circularConservationFastaInput: { value: null },
  reconcileCircularConservation: sources.reconcileCircularConservation,
  refreshCircularConservationSlots: sources.refreshCircularConservationSlots,
  askRingSourceRemoval: sources.askRingSourceRemoval
});
const files = createCircularConservationComparisonFiles({
  state, shownAction, writeLosatCacheInfo: () => {}, runUndoable: async () => undefined,
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
 * The walk's model: each ring by its own id, with its files' bases, its label, colour and LOSAT genetic code, and its
 * stack row's id and width once the stack has made it.
 * @typedef {{ id: number, bases: Set<string>, label: string, color: string, code: number,
 *   row: { id: string, width: unknown } | null }} ModelRing
 */

const ringRows = (drawing) => drawing.adv.circular_track_slots.filter((slot) => isManagedConservationSlot(slot));

/** The drawing's rings as drawn, in order: each with its genome's bases, label, colour, genetic code and stack row. */
const drawnRings = async (drawing) => {
  const genomes = ringGenomeFilesInRingOrder(drawing);
  const drawn = orderedConservationSources(circularRingSources(drawing), drawing.circularConservation);
  return Promise.all(drawn.map(async (source, index) => ({
    bases: await Promise.all((genomes[index] || []).map(async (file) => basesOfText(await readFileText(state.inputFiles.fileOf(file))))),
    label: source.label,
    color: source.color,
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
   * The last Generate's request and its rings, while the draft's rings differ from them only by a move, a label or
   * colour edit or a row width; a rerender of the committed Session adopts them again.
   * @type {{ canonical: Parameters<typeof sessionConfig.adoptCanonicalRenderArtifacts>[0],
   *   rings: { id: number, label: string, color: string }[] } | null}
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
    assert.deepEqual(ringIdentityDefects(drawing), [], `the rings' identity invariant (${where})`);
    const drawn = await drawnRings(drawing);
    assert.equal(drawn.length, model.length, `one drawn ring per ring (${where})`);
    drawn.forEach((ring, index) => {
      const expected = model[index];
      assert.ok(ring.bases.length > 0 && ring.bases.every((bases) => expected.bases.has(bases)) && ring.bases.length === expected.bases.size,
        `ring ${index} draws ring ${expected.id}'s files (${where})`);
      assert.deepEqual([ring.label, ring.color, ring.code], [expected.label, expected.color, expected.code],
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
    const searched = orderedConservationSources(ringGenomeFileLists(drawing.files.c_conservation_fastas), drawing.circularConservation);
    const searchedRings = await Promise.all(searched.map(async (source) => ({
      bases: await Promise.all(source.files.map(async (file) => basesOfText(await readFileText(state.inputFiles.fileOf(file))))),
      code: source.losat_gencode
    })));
    assert.equal(searchedRings.length, model.length, `the next Generate searches every ring (${where})`);
    searchedRings.forEach((ring, index) => {
      const expected = model[index];
      assert.ok(ring.bases.length === expected.bases.size && ring.bases.every((bases) => expected.bases.has(bases)),
        `the next Generate searches ring ${expected.id}'s files at ${index} (${where})`);
      assert.equal(ring.code, expected.code, `the next Generate searches ring ${expected.id} with its genetic code (${where})`);
    });
  };
  const replayed = () => drawing.files.c_conservation_blasts_source === 'losat-cache' && drawing.files.c_conservation_blasts.length > 0;
  const newRingOf = (index, bases) => {
    const entry = drawing.circularConservation.series[index];
    ringCount += 1;
    return { id: ringCount, bases: new Set([bases]), label: entry.label, color: entry.color, code: entry.losat_gencode, row: null };
  };

  /** @type {Record<string, { weight: number, can: () => boolean, run: () => Promise<void> | void }>} */
  const STEPS = {
    'add seq': {
      weight: 4, can: () => model.length < 6,
      run: async () => {
        committed = null;
        const file = genomeFile(random() < 0.6);
        const bases = basesOfText(atob(file.data));
        await files.addCircularConservationComparisonFile(/** @type {any} */ ({ target: { files: [file], value: 'x' } }));
        const index = drawing.circularConservation.series.length - 1;
        assert.equal(drawing.circularConservation.series[index].label, file.name.replace(/\.[^.]+$/, ''),
          `a new ring takes its file's name (seed ${seed}: ${trail.join(' > ')})`);
        model.push(newRingOf(index, bases));
      }
    },
    'remove series': {
      weight: 2, can: () => model.length > 0,
      run: () => {
        committed = null;
        const index = pick(model.length);
        files.removeCircularConservationSource(index);
        model = model.filter((_, ring) => ring !== index);
      }
    },
    'row add file': {
      weight: 2, can: () => model.length > 0,
      run: () => {
        committed = null;
        const index = pick(model.length);
        const file = random() < 0.6 ? textFile('plasmid.fasta', `>p\n${nextBases()}\n`) : genomeFile(false);
        files.addRingRowFiles(index, /** @type {any} */ ({ target: { files: [file], value: 'x' } }));
        model[index].bases.add(basesOfText(atob(file.data)));
      }
    },
    'row remove file': {
      weight: 2, can: () => model.some((ring) => ring.bases.size > 1),
      run: async () => {
        committed = null;
        const candidates = model.map((ring, index) => (ring.bases.size > 1 ? index : -1)).filter((index) => index >= 0);
        const index = candidates[pick(candidates.length)];
        const genome = ringGenomeFilesInRingOrder(drawing)[index];
        const fileIndex = pick(genome.length);
        const bases = basesOfText(await readFileText(state.inputFiles.fileOf(genome[fileIndex])));
        files.removeRingRowFile(index, fileIndex);
        model[index].bases.delete(bases);
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
        entry.label = `ring ${model[index].id}`;
        entry.color = `#00${String(model[index].id).padStart(2, '0')}aa`;
        panel.refreshShownCircularConservationSlots();
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
      weight: 3, can: () => model.length > 0,
      run: async () => {
        const filesData = await sessionConfig.serializeActiveRenderFiles(state, drawing);
        const genomes = ringGenomeFileLists(filesData.c_conservation_fastas).filter((list) => list.length > 0);
        const rings = orderedConservationSources(circularRingSources(drawing), drawing.circularConservation);
        const canonical = buildCanonicalRenderRequest({
          state, drawing, filesData,
          resolvedCircularConservation: rings.map((ring, index) => ({
            name: `${ring.fileName.replace(/\.[^.]+$/, '')}.circular_conservation.losatn.tsv`,
            text: `q\ts\t100\t8\t0\t0\t1\t8\t1\t8\t1e-9\t2${index}\n`,
            sourceIndex: ring.sourceIndex, label: ring.label, color: ring.color, fastas: genomes[ring.sourceIndex] || []
          }))
        });
        sources.adoptCircularRingResult(drawing, () => sessionConfig.adoptCanonicalRenderArtifacts(canonical, {
          adoptOwnedRequest: true, drawing, conservationPairs: 'committed'
        }), true);
        assert.ok(replayed(), 'the Generate adopted its LOSAT-cache result');
        committed = { canonical, rings: model.map(({ id, label, color }) => ({ id, label, color })) };
      }
    },
    'rerender the committed rings': {
      // A rerender of the committed Session (a Rotate): the adoption of the last Generate's request when the draft no
      // longer is that request; each ring comes back in its committed order, label and colour, in its own row, with the
      // genetic code it has in the draft, which the request does not carry (OV-3178). Most
      // of these steps move two rings first, as the review's case did (a move, then a Rotate).
      weight: 3, can: () => committed !== null,
      run: () => {
        const { canonical, rings } = /** @type {NonNullable<typeof committed>} */ (committed);
        if (model.length > 1 && random() < 0.7) {
          const index = pick(model.length - 1);
          panel.moveCircularConservationSeries(index, 1);
          trail[trail.length - 1] += ` (after moving ring ${index} down)`;
        }
        sources.adoptCircularRingResult(drawing, () => sessionConfig.adoptCanonicalRenderArtifacts(canonical, {
          adoptOwnedRequest: true, drawing, conservationPairs: 'committed'
        }), false);
        model = rings.map(({ id, label, color }) => Object.assign(/** @type {ModelRing} */ (model.find((ring) => ring.id === id)), { label, color }));
      }
    },
    'replace a replayed table': {
      weight: 1, can: () => replayed(),
      run: () => {
        committed = null;
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
        committed = null;
        panel.changeCircularConservationSource('upload');
        panel.changeCircularConservationSource('losat');
      }
    },
    'save and load': {
      weight: 2, can: () => true,
      run: async () => {
        committed = null;
        // Save Session, then Load of what it wrote, then the Load's reconcile (`afterImport`). One of three Loads of
        // LOSAT rings could not read one ring's genome files: that ring goes.
        const result = await exportSession('ring walk', { askChoice: async () => 'download' });
        const saved = JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
        const lost = !replayed() && model.length > 0 && random() < 0.34 ? pick(model.length) : -1;
        if (lost >= 0) {
          const place = drawing.circularConservation.series[lost].sourceIndex;
          const savedDrawing = saved.drawings.find((candidate) => candidate.mode === 'circular');
          savedDrawing.webFiles.bindings.c_conservation_fastas[place].forEach((file, index) => { file.resourceId = `missing-${place}-${index}`; });
          model = model.filter((_, ring) => ring !== lost);
          trail[trail.length - 1] += ` (lost ring ${lost})`;
          kinds.set('a Load that lost a genome', (kinds.get('a Load that lost a genome') || 0) + 1);
        }
        const before = drawing;
        const loaded = await importSession({ target: { files: [new File([JSON.stringify(saved)], 'walk.gbdraw.json', { type: 'application/json' })], value: 'x' } },
          { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
        assert.equal(loaded.status, 'ok', `the Load (seed ${seed}: ${trail.join(' > ')})`);
        drawing = circular();
        assert.notEqual(drawing, before, 'the Load installed new drawings');
        sources.reconcileLoadedCircularConservation(drawing, true);
      }
    },
    'stack off and on': {
      weight: 2, can: () => true,
      run: () => {
        editor.setCircularTrackSlotsEnabled(!(drawing.adv.circular_track_slots_enabled === true));
      }
    }
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
    'set a genetic code', 'widen a stack row', 'generate', 'rerender the committed rings', 'replace a replayed table', 'source radio round trip', 'save and load',
    'a Load that lost a genome', 'stack off and on']) {
    assert.ok((drawnKinds.get(kind) || 0) >= 3, `the walk draws ${kind} (${drawnKinds.get(kind) || 0} times)`);
  }
});
