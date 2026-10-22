// The JavaScript split of a Session 27-44 draft into Session 46 mode slices,
// against the vectors it shares with the Python twin
// (tests/fixtures/sessions/mode-split-vectors.json, tests/test_mode_split_vectors.py).
//
// A flat-draft case runs `splitDraftIntoModes` on its input. A fixture case
// loads the Session as the app does (`importSession`) and reads the split Load
// made (`session-draft-split`). Each case checks its `expect` and `expectAbsent`
// pointers and, where present, `expectedModes` (the slices the Python split
// must equal). The split writes only migrated values; `loadDefaults` lists the
// values Session 46 Load fills for the pointers the split leaves absent, and
// this test checks them on the drawings Load installs.
import { acceptNodeFeatureRecoveryRow, installSessionImportWorker } from './helpers/session-import-node.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

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

/** @type {Record<string, any>[]} */
const splits = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onSessionLifecycleEvent: (event) => {
    if (event.name === 'session-draft-split') splits.push(event.split);
  }
};

const { sessionLegendEditsByRowKey, splitDraftIntoModes } = await import('../../gbdraw/web/js/services/mode-scoped-migration.js');
const {
  buildModeSliceData, importSession, setUnmanagedConfigOverrideValidator
} = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
// The composition root's transform of an older Session's Results (R13 port).
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
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
const savedModes = () => Object.fromEntries(MODES.map((mode) => [mode, JSON.parse(JSON.stringify(
  buildModeSliceData(state.drawings[mode], mode, {
    selectedFeatureRecordIdx: mode === state.mode.value ? Number(state.selectedFeatureRecordIdx.value) || 0 : 0
  })
))]));

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

/** The split of one case, seen as a (partial) Session 46. */
const splitOf = async (vector) => {
  if (!vector.fixture && vector.input?.format !== 'gbdraw-session') {
    return splitDraftIntoModes(structuredClone(vector.input), {
      committedMode: vector.context.committedMode,
      modeProfiles: vector.context.modeProfiles,
      depthSources: vector.context.depthSources ?? null,
      paletteColors: vector.context.paletteColors ?? null
    });
  }
  const session = vector.fixture ? readFixture(vector.fixture) : structuredClone(vector.input);
  for (const pointer of vector.omit || []) remove(session, pointer);
  splits.length = 0;
  await load(session);
  assert.equal(splits.length, 1, 'Load splits a Session 27-44 draft once');
  return { ...splits[0], version: 46 };
};

for (const vector of VECTORS.cases) {
  test(vector.name, async () => {
    const result = await splitOf(vector);
    const mismatches = [];
    for (const [pointer, expected] of Object.entries(vector.expect || {})) {
      const actual = resolve(result, pointer);
      if (actual === ABSENT || canonical(actual) !== canonical(expected)) {
        mismatches.push(`${pointer}: expected ${canonical(expected).slice(0, 200)}, got ${
          actual === ABSENT ? '<absent>' : canonical(actual).slice(0, 200)}`);
      }
    }
    for (const pointer of vector.expectAbsent || []) {
      const actual = resolve(result, pointer);
      if (actual !== ABSENT) mismatches.push(`${pointer}: expected absent, got ${canonical(actual).slice(0, 200)}`);
    }
    if (vector.expectedModes) {
      for (const mode of MODES) {
        if (canonical(result.modes?.[mode]) !== canonical(vector.expectedModes[mode])) {
          mismatches.push(`/modes/${mode}: differs from expectedModes`);
        }
      }
    }
    assert.deepEqual(mismatches, []);

    // Session 46 Load fills what the split left absent.
    const defaults = vector.loadDefaults || [];
    if (!defaults.length) return;
    if (!vector.fixture) await load(settingsOnlySession(result));
    const loaded = { modes: savedModes() };
    const missing = defaults.flatMap(({ pointer, value, projected }) => {
      const actual = resolve(loaded, pointer);
      if (actual === ABSENT) return [`${pointer}: Load left it absent`];
      // The committed mode's request projection fills it first.
      if (projected) return [];
      return canonical(actual) === canonical(value) ? [] : [`${pointer}: expected ${canonical(value)}, got ${canonical(actual)}`];
    });
    assert.deepEqual(missing, []);
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
  const drawing = state.drawings.circular;
  assert.deepEqual({ ...drawing.legendColorOverrides }, { 'GC content': '#ff0000', 'GC skew (-)': '#0000ff' });
  assert.deepEqual(Object.keys(drawing.legendStrokeOverrides), ['GC skew (+)']);
  assert.equal(drawing.legendStrokeOverrides['GC skew (+)'].strokeColor, '#00aa00');
  assert.deepEqual(drawing.legendEntries.value.map(({ caption, originalCaption }) => [caption, originalCaption]).slice(3),
    [['GC skew (+)', 'GC content'], ['Z', 'GC skew (+)']]);
  assert.deepEqual(drawing.deletedLegendEntries.value.map(({ caption, originalCaption }) => [caption, originalCaption]),
    [['W', 'GC skew (-)']]);
});
