// Every Session fixture (tests/fixtures/sessions/, Sessions 27 to 51 in JSON
// files and Session file containers) and every Gallery Session goes through
// Load, Save and Load in the Web engine, on the real app (`createAppSetup`):
// - G3: the first Load lists no repair (a repair there is a producer bug);
// - G2b: the document the Save writes reads again with no document repair;
// - the second Load lists no repair, and each Load leaves the same state: the
//   drawing tabs, the shown drawing, and every slice of every drawing, read
//   by the one slice reader (`drawingSliceTexts`, the reader the random walk
//   and the Gallery round-trip audit hash) with the runtime identities a Load
//   gives anew renamed.
// Each state is read after the record read a user makes after a Load (Reload
// Records in each drawing, Product Contract PD-OI-044: a Load of a preview
// reads no records), so a Load that loses a value that read rewrites fails.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { readFileSync, readdirSync } from 'node:fs';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { readSavedSession } from './helpers/session-documents.mjs';
import { drawingSliceTexts } from './helpers/drawing-slice-hashes.mjs';

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
// after a Load; the drawing shown before is shown again.
const readRecords = async () => {
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
  }
  await show(shown);
};

// The mark a Load sets on the Circular File cards an older Session's flat
// record selectors count over (`files.c_legacy_record_cards`, OV-2867): the
// record read after the Load places those selectors by record identity
// (`placeFlatRecordSelectors`), so after it no setting depends on the mark,
// and a Save writes none. The state is compared as a Save writes it.
const LOAD_ONLY = Object.freeze({ inputs: [{ path: 'files.c_legacy_record_cards' }] });
/**
 * The state a Load leaves, after the record read: the drawing tabs, the shown
 * drawing, and every slice of every drawing (runtime identities renamed).
 * @returns {{ active: string, tabs: { id: string, name: string }[], drawings: Record<string, Record<string, unknown>> }}
 */
const loadedState = () => {
  const { active, renamed } = drawingSliceTexts(state, DRAWING_SLICES, { omit: LOAD_ONLY });
  return {
    active,
    tabs: state.drawingTabs.value.map(({ id, name }) => ({ id, name })),
    drawings: Object.fromEntries(Object.entries(renamed).map(([id, slices]) => [id,
      Object.fromEntries(Object.entries(slices).map(([name, text]) => [name, JSON.parse(text)]))]))
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

test('every Session fixture and Gallery Session loads without a repair (G3), its Save reads again without a document repair (G2b), and a second Load gives the same drawings', async () => {
  /** @param {string} dir @param {(name: string) => boolean} keep */
  const sessions = (dir, keep) => readdirSync(new URL(dir, import.meta.url)).filter(keep).sort()
    .map((name) => /** @type {const} */ ([name, new URL(`${dir}${name}`, import.meta.url)]));
  const corpus = [
    ...sessions('../fixtures/sessions/', (name) => /\.(json(\.gz)?|gbdraw)$/.test(name) && !/provenance|expected|evidence|vectors/.test(name)),
    ...sessions('../../gbdraw/web/gallery/sessions/', (name) => /\.(json(\.gz)?|gbdraw)$/.test(name))
  ];
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
  /** @type {string[]} */
  const notSaved = [];
  // What each fixture did that the check forbids, all reported at the end.
  /** @type {{ fixture: string, step: string, found: unknown }[]} */
  const findings = [];
  /** @param {string} fixture @param {string} step @param {unknown} found */
  const note = (fixture, step, found) => { findings.push({ fixture, step, found }); };
  /** @param {LoadOutcome} loaded */
  const rows = (loaded) => loaded.repairs.map(({ code, path }) => `${code} ${path}`);
  try {
    for (const [name, url] of corpus) {
      const bytes = readFileSync(url);
      const session = await readSavedSession(new Blob([bytes]));
      versions.add(session.version);
      askedDiagramWorker = false;
      current = name;
      const first = await load(bytes, name);
      if (first.loaded.status !== 'ok') { note(name, 'Load', first.loaded.error ?? first.loaded.status); continue; }
      const nodeOnly = FEATURES_NEED_DIAGRAM_WORKER.has(name) ? ['feature-recovery-failed'] : [];
      const expected = [...(LOSSY_FIXTURES[name] ?? []), ...nodeOnly];
      if (JSON.stringify(first.loaded.repairs.map(({ code }) => code)) !== JSON.stringify(expected)) note(name, 'Load (G3)', rows(first.loaded));
      if (askedDiagramWorker !== nodeOnly.length > 0) note(name, 'Load', `asked the diagram Worker: ${askedDiagramWorker}`);
      await readRecords();
      const before = loadedState();
      // A Session 27-39 with a Result is saved only after a Generate (PD-OI-079).
      const written = await exportSession(`round trip ${name}`, { askChoice: async () => 'download', readKeptArtifacts: () => first.kept })
        .catch((error) => (error?.code === 'SESSION_SAVE_REQUIRES_GENERATE' && session.version < 40 ? null : { status: 'error', error }));
      if (!written) { notSaved.push(name); continue; }
      if (written.status !== 'saved') { note(name, 'Save', written.error?.code ?? written); continue; }
      const documentRepairs = readSessionDocument(await readSavedSession(written.blob)).repairs;
      if (documentRepairs.length) note(name, 'Save (G2b)', documentRepairs);
      current = `${name} (saved)`;
      const second = await load(written.blob, `${name}.gbdraw`);
      if (second.loaded.status !== 'ok') { note(name, 'Load of the saved file', second.loaded.error ?? second.loaded.status); continue; }
      if (second.loaded.repairs.length) note(name, 'Load of the saved file', rows(second.loaded));
      await readRecords();
      const found = differences(before, loadedState());
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
  assert.ok(notSaved.length < corpus.length / 2, notSaved.join());
  assert.deepEqual(findings, [], 'each fixture loads, saves and loads again to the same drawings');
});
