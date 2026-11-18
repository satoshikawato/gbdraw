// Every Session fixture (tests/fixtures/sessions/, Sessions 27 to 51 in JSON
// files and Session file containers) and every Gallery Session goes through
// Load, Save and Load in the Web engine, on the real app (`createAppSetup`):
// - G3: the first Load lists no repair (a repair there is a producer bug);
// - G2b: the document the Save writes reads again with no document repair;
// - the second Load lists no repair, and each Load leaves the same state: the
//   drawing tabs, the shown drawing, and every slice of every drawing, read
//   by the one slice reader (`drawingSliceTexts`, the reader the random walk
//   and the Gallery round-trip audit hash) with the runtime identities a Load
//   gives anew renamed, and the tables of each drawing's saved comparison by
//   their bytes.
// Each state is read after the record read a user makes after a Load (Reload
// Records in each drawing, Product Contract PD-OI-044: a Load of a preview
// reads no records), so a Load that loses a value that read rewrites fails.
// Each Circular drawing's read must end ready (or idle, with no input); the
// fixtures whose read needs the diagram Worker, which node lacks, are listed
// (`RECORDS_NEED_DIAGRAM_WORKER`): their read fails in both Loads alike, so
// they compare the state before the record read.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { createHash } from 'node:crypto';
import { readFileSync, readdirSync } from 'node:fs';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { readSavedSession } from './helpers/session-documents.mjs';
import { drawingSliceTexts } from './helpers/drawing-slice-hashes.mjs';
import { circularDiscoveryForInput } from '../../gbdraw/web/js/services/record-discovery.js';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }), querySelector: () => null, documentElement: { setAttribute: () => {} } };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
// A Load prefetches the Worker's boot files unless the browser saves data; a node test fetches nothing.
Object.defineProperty(globalThis.navigator, 'connection', { value: { saveData: true }, configurable: true });
installSessionImportWorker();

const config = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup, transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { readSessionDocument } = await import('../../gbdraw/web/js/services/session-authority.js');
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
const { setStoredProteinResultFrameConverter } = await import('../../gbdraw/web/js/services/main-session-comparison-frame.js');
const { callerSession } = await import('./helpers/session-caller.mjs');
const { DRAWING_SLICES } = await import('./helpers/drawing-state.mjs');
const { storedResourceBytes } = await import('../../gbdraw/web/js/services/session-resource-backing.js');

const app = createAppSetup();
// The app runs these ports in its diagram Worker; node runs them in Python.
// A Session 27-44 moves to a Session 49 in Python, as the composition root's
// Worker port does it (DESIGN 4.3).
setLegacySessionUpgrader(pythonLegacySessionUpgrader());
config.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
// The Load service with the app's rollback port: the app's own Load waits for
// the preview of a shown Result, which a node test does not mount. The test
// keeps the artifact slots of the drawings a Load does not show and gives
// them to the Save, as the app's slot store does (`installLoadedArtifactSlots`,
// `readKeptArtifacts`).
const { exportSession, importSession } = callerSession(config, state, app.sessionImportRollbackState);
// The app's asynchronous work after an action (a record read, a switch) runs in a few turns.
const settle = async () => { for (let i = 0; i < 20; i += 1) await new Promise((resolve) => setImmediate(resolve)); };

/** @typedef {{ status: string, error?: unknown, repairs: { code: string, path: string }[] }} LoadOutcome */
/** @typedef {{ drawingId: string }} KeptSlot An artifact slot of a drawing a Load did not show. */

/**
 * A Load of `bytes` as a file named `name`, answering the repair dialog with
 * Load, and the artifact slots it kept.
 * @param {Uint8Array | Blob} bytes @param {string} name
 * @returns {Promise<{ loaded: LoadOutcome, kept: KeptSlot[] }>}
 */
const load = async (bytes, name) => {
  /** @type {KeptSlot[]} */
  let kept = [];
  const loaded = await importSession({ target: { files: [new File([bytes], name)], value: '' } }, {
    transformLegacyResultSvg, confirmLoadRepairs: async () => true,
    installLoadedArtifactSlots: (/** @type {{ opening: KeptSlot | null, stashed: KeptSlot[] }} */ { opening, stashed }) => {
      kept = /** @type {KeptSlot[]} */ ([opening, ...stashed].filter(Boolean));
    }
  });
  return { loaded, kept };
};

// Reload Records in each drawing, shown by its drawing ID, as a user does
// after a Load; the drawing shown before is shown again. Returns each Circular
// drawing whose read did not end ready (`circularDiscoveryForInput`; idle is a
// drawing with no input), with its status and error code.
/** @returns {Promise<string[]>} */
const readRecords = async () => {
  /** @type {string[]} */
  const unread = [];
  const shown = state.activeDrawingId.value;
  /** @param {string} id */
  const show = async (id) => {
    if (state.activeDrawingId.value === id) return;
    const switched = app.selectDrawing(id);
    assert.equal(switched?.status, 'ok', `show drawing ${id}: ${JSON.stringify(switched)}`);
    await settle();
  };
  for (const drawing of state.drawingList()) {
    await show(drawing.id);
    const outcome = drawing.mode === 'circular' ? await app.inspectCircularSourceRecords() : await app.refreshLinearRecordSelectors();
    assert.notEqual(outcome?.status, 'busy', `Reload Records in drawing ${drawing.id}: ${JSON.stringify(outcome)}`);
    await settle();
    if (drawing.mode !== 'circular') continue;
    const { status, error } = circularDiscoveryForInput(state.drawingOf(drawing.id));
    if (status !== 'ready' && status !== 'idle') unread.push(`${drawing.id} ${status} ${typeof error === 'object' && error ? error.code : error}`);
  }
  await show(shown);
  return unread;
};

// The mark a Load sets on the Circular File cards an older Session's flat
// record selectors count over (`files.c_legacy_record_cards`, OV-2867): the
// record read after the Load places those selectors by record identity
// (`placeFlatRecordSelectors`), so after it no setting depends on the mark,
// and a Save writes none. The state is compared as a Save writes it.
//
// What Inherit copies (`importedComparisonIntent.source`) holds the saved
// tables as resource descriptors, which say how the file stored the bytes: a
// Session 27-50 holds them as base64 texts in its JSON, a Session file as its
// entries, which Load reads from there (`copyContainerInputs`), and a Save
// writes both as entries (`storedResourceBytes`). Each is compared as Inherit
// puts it into a request: its fields but the stored text, and the sha256 of
// the bytes a Save stores, read from either form (`comparisonSourceTables`).
const LOAD_ONLY = Object.freeze({
  inputs: [{ path: 'files.c_legacy_record_cards' }],
  settings: [{ path: 'importedComparisonIntent.source.resources' }]
});
/**
 * The tables of each drawing's saved comparison, by resource ID: each
 * descriptor's fields but its stored text, and the sha256 of its stored bytes.
 * @returns {Promise<Record<string, Record<string, unknown>>>}
 */
const comparisonSourceTables = async () => Object.fromEntries(await Promise.all(state.drawingList().map(async ({ id }) => {
  const resources = state.drawingOf(id).importedComparisonIntent.source?.resources ?? {};
  return [id, Object.fromEntries(await Promise.all(Object.entries(resources).map(async ([resourceId, descriptor]) => {
    const { data: _data, encoding: _encoding, ...fields } = /** @type {Record<string, unknown>} */ (descriptor);
    const bytes = new Uint8Array(await new Response(await storedResourceBytes(/** @type {never} */ (descriptor)).open()).arrayBuffer());
    return [resourceId, { ...fields, sha256: createHash('sha256').update(bytes).digest('hex') }];
  })))];
})));
/**
 * The state a Load leaves, after the record read: the drawing tabs, the shown
 * drawing, every slice of every drawing (runtime identities renamed), and the
 * tables of each drawing's saved comparison.
 * @returns {Promise<{ active: string, tabs: { id: string, name: string }[], drawings: Record<string, Record<string, unknown>>,
 *   comparisonTables: Record<string, Record<string, unknown>> }>}
 */
const loadedState = async () => {
  const { active, renamed } = drawingSliceTexts(state, DRAWING_SLICES, { omit: LOAD_ONLY });
  return {
    active,
    tabs: state.drawingTabs.value.map(({ id, name }) => ({ id, name })),
    drawings: Object.fromEntries(Object.entries(renamed).map(([id, slices]) => [id,
      Object.fromEntries(Object.entries(slices).map(([name, text]) => [name, JSON.parse(text)]))])),
    comparisonTables: await comparisonSourceTables()
  };
};
/**
 * The leaves where two states differ, each `path: before -> after`.
 * @param {unknown} before @param {unknown} after @param {string} [at]
 * @returns {string[]}
 */
const differences = (before, after, at = '') => {
  if (JSON.stringify(before) === JSON.stringify(after)) return [];
  if (before && after && typeof before === 'object' && typeof after === 'object') {
    const keys = [...new Set([...Object.keys(before), ...Object.keys(after)])];
    return keys.flatMap((key) => differences(/** @type {Record<string, unknown>} */ (before)[key],
      /** @type {Record<string, unknown>} */ (after)[key], `${at}.${key}`));
  }
  const shown = (/** @type {unknown} */ value) => (value === undefined ? '(none)' : JSON.stringify(value).slice(0, 200));
  return [`${at.slice(1)}: ${shown(before)} -> ${shown(after)}`];
};

// Two fixtures are lossy on purpose. Node supplies no diagram Worker: a
// Session older than 40 reads its features again through it, so the fixtures
// below list that recovery here (and only because the Worker is missing: each
// asks for it). The browser specs load them (OV-751): gui-audit-20260930-inputs
// and live-generate-parity assert no dialog for seven; session-loading-feedback
// answers the one table row of whitelist-tab-keyword.v39 and loads the two
// mode-switch v30 fixtures without a dialog: their recovery reads the inputs of
// the Result's mode (OV-677). Five fixtures needed a Worker service to load at
// all until Python converted a Session 27-44 (OV-1958, D-15-140); they load in
// node now, the cropped or reversed ones asking the diagram Worker for their
// features. The browser specs feature-identity-overrides and
// gui-audit-20260930-comparisons load them.
const LOSSY_FIXTURES = {
  'feature-edits-circular-copies.v44.gbdraw-session.json.gz': ['feature-visibility-narrowed'],
  'whitelist-tab-keyword.v39.gbdraw-session.json.gz': ['table-row-joined']
};
const FEATURES_NEED_DIAGRAM_WORKER = new Set([
  'BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz', 'BGC0000708-BGC0000713.v29.gbdraw-session.json.gz',
  'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz', 'cli-linear-protein.v30.gbdraw-session.json.gz', 'conservation-fasta.v39.gbdraw-session.json.gz',
  'lambda_basic_linear.v31-schema1.gbdraw-session.json.gz', 'precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz', 'scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz',
  'scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz', 'whitelist-tab-keyword.v39.gbdraw-session.json.gz',
  'mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz', 'mode-switch-linear-result-circular-shown.v30.gbdraw-session.json.gz',
  'legacy-font-size-zero-circular-cli.v30.gbdraw-session.json.gz', 'legacy-font-size-zero-linear-cli.v30.gbdraw-session.json.gz',
  'feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz', 'hash-color-rule-linear-rc.v30.gbdraw-session.json.gz',
  'cli-web-mito-circular.v28-cli.gbdraw-session.json.gz',
  // gui-audit-20260930-comparisons loads it in the browser (OV-1181).
  'BGC0000708-BGC0000709-rc.v39.gbdraw-session.json.gz'
]);

// The slice reader names a runtime identity by one name across a drawing's
// slices: a Load that gives a drawing new cards and keeps rows naming other
// cards reads differently from a Load whose rows name its own cards, though
// each slice alone reads the same (OV-3714's channel: the request places a
// card by the row that names it).
test('a drawing whose record rows name cards it does not hold reads differently from one whose rows name its cards', () => {
  const asCaptured = { capture: (/** @type {unknown} */ handle) => handle };
  /** @param {string[]} cards @param {string[]} rowCards */
  const renamed = (cards, rowCards) => drawingSliceTexts({
    activeDrawing: () => ({ id: 'linear' }),
    drawingList: () => [{ id: 'linear', slices: {
      settings: { linearRecordRows: rowCards.map((uid) => ({ uid, row: 1 })) }, inputs: { linearSeqs: cards.map((uid) => ({ uid })) }
    } }]
  }, { settings: asCaptured, inputs: asCaptured }).renamed.linear;
  const first = renamed(['linear-seq-a1', 'linear-seq-b1'], ['linear-seq-a1', 'linear-seq-b1']);
  assert.deepEqual(renamed(['linear-seq-a2', 'linear-seq-b2'], ['linear-seq-a2', 'linear-seq-b2']), first);
  assert.notDeepEqual(renamed(['linear-seq-a2', 'linear-seq-b2'], ['linear-seq-c2', 'linear-seq-d2']), first);
});

// The fixtures whose record read needs the diagram Worker: a Circular FASTA
// card's records are listed by Python (`listSequenceRecords`), so in node the
// read ends with WORKER_INIT in both Loads. The browser specs read them
// (circular-multi-file-input).
const RECORDS_NEED_DIAGRAM_WORKER = new Set(['cli-circular-fasta.v51.gbdraw']);

// The Sessions: every Session file (`.gbdraw`) and every JSON file (`.json`,
// `.json.gz`) that holds a Session envelope (`format` 'gbdraw-session' and an
// integer `version`, `SessionEnvelope` of contracts/session.generated.js; a
// provenance file names the format but holds no version), at any depth of the
// folders. tests/fixtures/performance/ is left out for CI time: its one
// Session, an 11 MB Session 40 of the Gallery's hepatoplasmataceae_collinear
// (whose current Session file the corpus holds), takes 23 s and 2.4 GB of
// memory for its round trip alone in node.
/** @param {string} dir */
const sessionsIn = async (dir) => {
  const root = new URL(dir, import.meta.url);
  /** @type {{ name: string, bytes: Buffer, session: Record<string, any> }[]} */
  const found = [];
  for (const path of readdirSync(root, { recursive: true, encoding: 'utf8' }).sort()) {
    if (!/\.(gbdraw|json|json\.gz)$/.test(path)) continue;
    const bytes = readFileSync(new URL(path, root));
    const session = await readSavedSession(new Blob([bytes]));
    if (path.endsWith('.gbdraw') || (session?.format === 'gbdraw-session' && Number.isInteger(session.version))) found.push({ name: path, bytes, session });
  }
  return found;
};
// A Session 27-39 with a Result is saved only after a Generate (PD-OI-079).
/** @param {Record<string, any>} session */
const savedOnlyAfterGenerate = (session) => session.version < 40 && Array.isArray(session.results) && session.results.length > 0;

test('every Session fixture and Gallery Session loads without a repair (G3), its Save reads again without a document repair (G2b), and a second Load gives the same drawings', async () => {
  const corpus = [...await sessionsIn('../fixtures/sessions/'), ...await sessionsIn('../../gbdraw/web/gallery/sessions/')];
  assert.ok(corpus.length > 100, String(corpus.length));
  const versions = new Set();
  const ImportWorker = globalThis.Worker;
  let askedDiagramWorker = false;
  globalThis.Worker = class extends ImportWorker {
    /** @param {string | URL} url */
    constructor(url) { askedDiagramWorker ||= !String(url).endsWith('/session-import-worker.js'); super(url); }
  };
  // The protein result conversion answers "nothing changed"; it is asked for
  // exactly the Sessions 46-47 that stored protein results of a reversed
  // record (Python's conversion of a Session 27-44 moves its own, OV-1969).
  const askedProteinConversion = new Set();
  let current = '';
  setStoredProteinResultFrameConverter(async () => { askedProteinConversion.add(current); return { result: { resources: {} } }; });
  // What each fixture did that the check forbids, all reported at the end.
  /** @type {{ fixture: string, step: string, found: unknown }[]} */
  const findings = [];
  /** @param {string} fixture @param {string} step @param {unknown} found */
  const note = (fixture, step, found) => { findings.push({ fixture, step, found }); };
  /** @param {LoadOutcome} loaded */
  const rows = (loaded) => loaded.repairs.map(({ code, path }) => `${code} ${path}`);
  // The record read of a Load: every Circular drawing ends ready, or, for a
  // fixture whose read needs the diagram Worker, asks it and fails to start it.
  // Each Linear record row then names a card of its drawing (a row naming
  // another card is a stale row that a Save would write).
  /** @param {string} name @param {string} step */
  const readRecordsOf = async (name, step) => {
    askedDiagramWorker = false;
    const unread = await readRecords();
    for (const { id, linearSeqs, linearRecordRows } of state.drawingList()) {
      const cards = new Set(linearSeqs.map(({ uid }) => uid));
      const stale = linearRecordRows.filter(({ uid }) => !cards.has(uid));
      if (stale.length) note(name, step, `${id}: record rows name no card of the drawing: ${JSON.stringify(stale)}`);
    }
    const needsWorker = RECORDS_NEED_DIAGRAM_WORKER.has(name);
    if (askedDiagramWorker !== needsWorker) note(name, step, `the record read asked the diagram Worker: ${askedDiagramWorker}`);
    if (needsWorker ? !unread.length || unread.some((entry) => !entry.endsWith(' error WORKER_INIT')) : unread.length) {
      note(name, step, `record read not ready: ${unread.join('; ') || '(none failed)'}`);
    }
  };
  try {
    for (const { name, bytes, session } of corpus) {
      versions.add(session.version);
      askedDiagramWorker = false;
      current = name;
      const first = await load(bytes, name);
      if (first.loaded.status !== 'ok') { note(name, 'Load', first.loaded.error ?? first.loaded.status); continue; }
      const nodeOnly = FEATURES_NEED_DIAGRAM_WORKER.has(name) ? ['feature-recovery-failed'] : [];
      const expected = [...(LOSSY_FIXTURES[name] ?? []), ...nodeOnly];
      if (JSON.stringify(first.loaded.repairs.map(({ code }) => code)) !== JSON.stringify(expected)) note(name, 'Load (G3)', rows(first.loaded));
      if (askedDiagramWorker !== nodeOnly.length > 0) note(name, 'Load', `asked the diagram Worker: ${askedDiagramWorker}`);
      await readRecordsOf(name, 'Reload Records after the Load');
      const before = await loadedState();
      const written = await exportSession(`round trip ${name}`, { askChoice: async () => 'download', readKeptArtifacts: () => first.kept })
        .catch((error) => ({ status: 'error', error }));
      if (savedOnlyAfterGenerate(session)) {
        // Node cannot Generate: such a Session gets the checks of the first Load.
        if (written.error?.code !== 'SESSION_SAVE_REQUIRES_GENERATE') note(name, 'Save (PD-OI-079)', written.error?.code ?? written.status);
        continue;
      }
      if (written.status !== 'saved') { note(name, 'Save', written.error?.code ?? written); continue; }
      const documentRepairs = readSessionDocument(await readSavedSession(written.blob)).repairs;
      if (documentRepairs.length) note(name, 'Save (G2b)', documentRepairs);
      current = `${name} (saved)`;
      askedDiagramWorker = false;
      const second = await load(written.blob, `${name}.gbdraw`);
      if (second.loaded.status !== 'ok') { note(name, 'Load of the saved file', second.loaded.error ?? second.loaded.status); continue; }
      if (second.loaded.repairs.length) note(name, 'Load of the saved file', rows(second.loaded));
      if (askedDiagramWorker) note(name, 'Load of the saved file', 'asked the diagram Worker');
      await readRecordsOf(name, 'Reload Records after the Load of the saved file');
      const found = differences(before, await loadedState());
      if (found.length) note(name, 'Load, Save and Load', found);
    }
  } finally {
    globalThis.Worker = ImportWorker;
    setStoredProteinResultFrameConverter(null);
  }
  assert.deepEqual([...askedProteinConversion].sort(), [
    'BGC0000708-BGC0000713.v47.gbdraw-session.json.gz', 'majanivirus_collinear-rc.v47.gbdraw-session.json.gz'
  ]);
  // Sessions 27 and 28 included (the oldest versions Load reads).
  assert.ok(versions.has(27) && versions.has(28), [...versions].sort().join());
  assert.deepEqual(findings, [], 'each fixture loads, saves and loads again to the same drawings');
});
