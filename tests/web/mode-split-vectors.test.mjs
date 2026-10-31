// Python's split of a Session 27-44 draft into Session 46 mode slices,
// against the vectors it shares with the Web Load
// (tests/fixtures/sessions/mode-split-vectors.json, tests/test_mode_split_vectors.py).
//
// A flat-draft case runs Python's `split_draft_into_modes` on its input and checks its
// `expect` and `expectAbsent` pointers and, where present, `expectedModes` (the
// slices the Python split must equal). A whole-Session case loads the Session
// as the app does (`importSession`): Python converts it once (D-15-129), and
// the case checks its `/modes/` `expect` pointers on Python's drawings, that no
// flat draft outlives the conversion, and that Load installs Python's drawings
// as a Load of that Session 49 does; tests/test_mode_split_vectors.py checks
// Python's split against the whole case. `loadDefaults` lists the values Load
// fills for the pointers the split leaves absent, checked on Load's drawings.
import { acceptNodeFeatureRecoveryRow, installSessionImportWorker } from './helpers/session-import-node.mjs';
import { pythonDraftSplit, pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};
installFakeSvgDom();
globalThis.alert = () => {};
installSessionImportWorker();

const { sessionLegendEditsByRowKey } = await import('../../gbdraw/web/js/services/session-drawings-upgrade.js');
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const {
  buildModeSliceData,
  setUnmanagedConfigOverrideValidator
} = sessionConfig;
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
// The composition root's transform of an older Session's Results (R13 port).
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
// A Session 27-44 moves to a Session 49 in Python, as the composition root's
// Worker port does it (DESIGN 4.3).
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
// Each conversion Python returned to Load: a fixture case's oracle.
/** @type {Record<string, any>[]} */
const conversions = [];
const pythonUpgrade = pythonLegacySessionUpgrader();
setLegacySessionUpgrader(async (request) => {
  const reply = await pythonUpgrade(request);
  conversions.push(reply.result.document);
  return reply;
});
const { normalizePaletteDefinitions } = await import('../../gbdraw/web/js/utils/color-utils.js');

const REPO = new URL('../../', import.meta.url);
// The Worker's typed override check, run by Python (as session-cli-compatibility does).
setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const VECTORS = JSON.parse(readFileSync(new URL('tests/fixtures/sessions/mode-split-vectors.json', REPO), 'utf8'));
const palettes = JSON.parse(readFileSync(new URL('gbdraw/web/gallery/palettes/palettes.json', REPO), 'utf8'));
state.paletteDefinitions.value = normalizePaletteDefinitions(palettes.palettes || palettes);
const MODES = ['circular', 'linear'];
const ABSENT = Symbol('absent');

/** The value at an RFC 6901 JSON pointer, or ABSENT. */
const resolve = (document, pointer) => {
  let current = document;
  for (const raw of pointer.split('/').slice(1)) {
    const token = raw.replace(/~1/g, '/').replace(/~0/g, '~');
    if (current && typeof current === 'object' && Object.hasOwn(current, token)) current = current[token];
    else return ABSENT;
  }
  return current;
};
const remove = (document, pointer) => {
  const tokens = pointer.split('/').slice(1).map((raw) => raw.replace(/~1/g, '/').replace(/~0/g, '~'));
  const last = tokens.pop();
  const parent = tokens.reduce((current, token) => current[token], document);
  if (Array.isArray(parent)) parent.splice(Number(last), 1);
  else delete parent[last];
};
const sortKeys = (value) => (Array.isArray(value) ? value.map(sortKeys)
  : value && typeof value === 'object'
    ? Object.fromEntries(Object.keys(value).sort().map((key) => [key, sortKeys(value[key])]))
    : value);
const canonical = (value) => JSON.stringify(sortKeys(value));

const readFixture = (path) => {
  const bytes = readFileSync(new URL(path, REPO));
  return JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));
};
const load = async (session) => {
  const outcome = await importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, { transformLegacyResultSvg, confirmLoadRepairs: acceptNodeFeatureRecoveryRow });
  assert.equal(outcome?.status, 'ok', JSON.stringify(outcome?.error));
};
// The drawings as Save writes them (complete slices).
// A one-mode Load made no drawing of the other mode (D-15-120); the offer tab makes it at its defaults,
// which Load used to fill.
const addPairDrawings = await (await import('./helpers/drawing-state.mjs')).modeDrawingsAdder();
const savedModes = () => {
  addPairDrawings();
  return Object.fromEntries(MODES.map((mode) => [mode, JSON.parse(JSON.stringify(
    buildModeSliceData(modeDrawing(state, mode), mode, {
      selectedFeatureRecordIdx: mode === state.mode.value ? Number(state.selectedFeatureRecordIdx.value) || 0 : 0
    })
  ))]));
};

// A Session 46 that carries a flat-draft case's split and nothing else: Load
// fills each slice from the mode's defaults.
const settingsOnlyBase = readFixture('tests/fixtures/sessions/settings-only.v42.json.gz');
const settingsOnlySession = (split) => {
  const { config: _config, features: _features, ...session } = structuredClone(settingsOnlyBase);
  const editorState = { ...session.editorState, legend: {} };
  delete editorState.featureStrokes;
  const ui = Object.fromEntries(Object.entries(session.ui).filter(([key]) => ![
    'layoutPreferences', 'canvasPadding', 'pendingPaletteName', 'pendingPaletteColors', 'linearTypographyLinked'
  ].includes(key)));
  const modes = structuredClone(split.modes);
  // Session validation requires the comparison plan beside a Linear record
  // layout; the split of a real Session writes both (the synthetic cases do not).
  const linear = modes.linear?.config;
  if (linear && Object.hasOwn(linear, 'linearRecordLayout') && !Object.hasOwn(linear, 'linearComparisonPlan')) {
    linear.linearComparisonPlan = { mode: 'none', defaultSource: 'losat', edges: [] };
  }
  // A settings-only Session holds the `form` and `adv` of its shown mode, as
  // Save writes them (OV-591); a synthetic case may hold neither.
  const shownMode = split.ui?.mode || 'circular';
  modes[shownMode] = { ...modes[shownMode], config: { form: {}, adv: {}, ...modes[shownMode]?.config } };
  return {
    ...session, version: 46, editorState,
    ui: { ...ui, ...split.ui, mode: shownMode },
    modes
  };
};

/** Python's split of a flat-draft case, seen as a (partial) Session 46. */
const splitOf = (vector) => pythonDraftSplit(vector.input, vector.context);

/** A drawing as Save writes it (a complete slice). */
const savedSlice = (drawing) => JSON.parse(JSON.stringify(buildModeSliceData(drawing, drawing.mode, {
  selectedFeatureRecordIdx: drawing.id === state.activeDrawing().id ? Number(state.selectedFeatureRecordIdx.value) || 0 : 0
})));

/**
 * A whole-Session case loaded as the app loads it: Python converts the
 * Session 27-44 once (D-15-129; tests/test_mode_split_vectors.py pins that
 * split against this case), and Load installs Python's drawings. Returns
 * Python's conversion and Load's drawings by Python's drawing IDs, as Save
 * writes them (complete slices), keyed by mode.
 */
const loadedOf = async (vector) => {
  const session = vector.fixture ? readFixture(vector.fixture) : structuredClone(vector.input);
  for (const pointer of vector.omit || []) remove(session, pointer);
  conversions.length = 0;
  await load(session);
  assert.equal(conversions.length, 1, 'Load converts a Session 27-44 once, in Python');
  const [converted] = conversions;
  const modes = Object.fromEntries(converted.drawings.map(({ id, mode }) => {
    const drawing = state.drawingOf(id);
    assert.equal(drawing.mode, mode, id);
    return [mode, savedSlice(drawing)];
  }));
  // Python's drawings as the slices of a Session 46 (`modes`).
  const convertedModes = Object.fromEntries(converted.drawings.map(({ mode, config, features, editorState, ui }) => (
    [mode, { config, features, editorState, ui }])));
  return { converted, convertedModes, modes };
};

/**
 * The mismatches of `expect` pointers read on `view` (all, or those `keep`
 * keeps). `expected` defaults to the case's `expect`; a whole-Session case
 * passes `loadExpectOf(vector)`.
 */
const expectMismatches = (vector, view, keep = () => true, expected = vector.expect || {}) => Object.entries(expected)
  .filter(([pointer]) => keep(pointer)).flatMap(([pointer, expected]) => {
    const actual = resolve(view, pointer);
    return actual !== ABSENT && canonical(actual) === canonical(expected) ? [] : [`${pointer}: expected ${
      canonical(expected).slice(0, 200)}, got ${actual === ABSENT ? '<absent>' : canonical(actual).slice(0, 200)}`];
  });
/** The `expectAbsent` pointers that `view` holds (all, or those `keep` keeps). */
const absentMismatches = (vector, view, keep = () => true) => (vector.expectAbsent || []).filter(keep)
  .filter((pointer) => resolve(view, pointer) !== ABSENT)
  .map((pointer) => `${pointer}: expected absent, got ${canonical(resolve(view, pointer)).slice(0, 200)}`);
/** The `loadDefaults` that Load's drawings (`loaded`, as Save writes them) do not hold. */
const defaultMismatches = (vector, loaded) => (vector.loadDefaults || []).flatMap(({ pointer, value, projected }) => {
  const actual = resolve(loaded, pointer);
  if (actual === ABSENT) return [`${pointer}: Load left it absent`];
  // The committed mode's request projection fills it first.
  if (projected) return [];
  return canonical(actual) === canonical(value) ? [] : [`${pointer}: expected ${canonical(value)}, got ${canonical(actual)}`];
});
const isModePointer = (pointer) => pointer.startsWith('/modes/');
/**
 * The `expect` values of a whole-Session case as Web Load sends them to
 * Python: Load's first pass carries a field each drawing reads differently
 * from the stored draft as the read value (D-15-148 (A)), so `loadExpect`
 * holds the read value for such a carried field (e.g. "20" read as 20). The
 * pytest vector converts the stored draft and keeps `expect`.
 */
const loadExpectOf = (vector) => ({ ...(vector.expect || {}), ...(vector.loadExpect || {}) });

// `loadExpect` of F1: Load reads the retired spelling "either" as "strand" (OV-2301). Only that
// spelling moves: a stored value is never replaced by today's default.
test('Load keeps a stored Collinear setting that differs from today\'s default', async () => {
  const session = readFixture('tests/fixtures/sessions/two-mode-project.v44.gbdraw-session.json.gz');
  session.config.losat.blastp.collinearMergeOrientation = 'order';
  const { convertedModes, modes } = await loadedOf({ input: session });
  for (const { config } of [convertedModes.linear, modes.linear]) {
    assert.equal(config.losat.blastp.collinearMergeOrientation, 'order');
    assert.deepEqual(config.losat.blastp.hitLimitsByMode.collinear, { candidateLimit: 5, orthogroupMemberMaxHits: 5 });
  }
});

for (const vector of VECTORS.cases) {
  const wholeSession = Boolean(vector.fixture) || vector.input?.format === 'gbdraw-session';
  test(vector.name, async () => {
    if (wholeSession) {
      // Load's drawings against Python's conversion of the same Session.
      const { converted, convertedModes, modes } = await loadedOf(vector);
      const mismatches = [
        // The values the split moves into each mode, in Python's drawing of
        // that mode. A one-mode Session composes to one drawing (DESIGN 3.2):
        // the split's slice of the other mode is pinned in Python only.
        ...expectMismatches(vector, { modes: convertedModes },
          (pointer) => isModePointer(pointer) && Object.hasOwn(convertedModes, pointer.split('/')[2]),
          loadExpectOf(vector)),
        // The flat draft does not outlive the conversion.
        ...absentMismatches(vector, converted, (pointer) => !isModePointer(pointer)),
        // Load fills what the split left absent.
        ...defaultMismatches(vector, { modes })
      ];
      // Load installs Python's drawings as it installs them from a Session 49 file.
      conversions.length = 0;
      await load(converted);
      assert.equal(conversions.length, 0, 'a Session 49 is not converted');
      const direct = Object.fromEntries(converted.drawings.map(({ id, mode }) => [mode, savedSlice(state.drawingOf(id))]));
      for (const mode of Object.keys(modes)) {
        if (canonical(modes[mode]) !== canonical(direct[mode])) mismatches.push(`/modes/${mode}: differs from a Load of Python's Session 49`);
      }
      assert.deepEqual(mismatches, []);
      return;
    }
    const result = await splitOf(vector);
    const mismatches = [...expectMismatches(vector, result), ...absentMismatches(vector, result)];
    if (vector.expectedModes) {
      for (const mode of MODES) {
        if (canonical(result.modes?.[mode]) !== canonical(vector.expectedModes[mode])) {
          mismatches.push(`/modes/${mode}: differs from expectedModes`);
        }
      }
    }
    assert.deepEqual(mismatches, []);

    // Session 46 Load fills what the split left absent.
    if (!(vector.loadDefaults || []).length) return;
    await load(settingsOnlySession(result));
    assert.deepEqual(defaultMismatches(vector, { modes: savedModes() }), []);
  });
}

// The positive Session 46 fixture (first-parent main 9839de8e): its renamed
// rows' edits load under their row keys (tests/fixtures/sessions/README.md).
test('a Session 46 written by main loads its renamed rows\' edits under their row keys', async () => {
  const fixture = readFixture('tests/fixtures/sessions/legend-renamed-rows.v46.gbdraw-session.json.gz');
  assert.equal(fixture.version, 46);
  // The moves Load reports before it commits, as Python reports them (tests/test_session_io.py).
  assert.deepEqual(sessionLegendEditsByRowKey(fixture, fixture.version).outcomes, [
    ['colorOverrides', 'GC skew (+)', 'GC content'], ['colorOverrides', 'W', 'GC skew (-)'], ['strokeOverrides', 'Z', 'GC skew (+)']
  ].map(([store, fromKey, toKey]) => ({ store, fromKey, toKey, outcome: 'moved', slice: 'circular' })));
  await load(fixture);
  const drawing = modeDrawing(state, 'circular');
  assert.deepEqual({ ...drawing.legendColorOverrides }, { 'GC content': '#ff0000', 'GC skew (-)': '#0000ff' });
  assert.deepEqual(Object.keys(drawing.legendStrokeOverrides), ['GC skew (+)']);
  assert.equal(drawing.legendStrokeOverrides['GC skew (+)'].strokeColor, '#00aa00');
  assert.deepEqual(drawing.legendEntries.value.map(({ caption, originalCaption }) => [caption, originalCaption]).slice(3),
    [['GC skew (+)', 'GC content'], ['Z', 'GC skew (+)']]);
  assert.deepEqual(drawing.deletedLegendEntries.value.map(({ caption, originalCaption }) => [caption, originalCaption]),
    [['W', 'GC skew (-)']]);
});
