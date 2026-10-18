// A mode switch changes no drawing (plan OV-80 §4.5, §8 "Mode switch"; PR-0b
// open items). Each mode has its own drawing, so the switch sets `mode`, swaps
// the artifact slots, and writes no setting or edit:
//   1. no watcher's source changes on a switch except the ones that watch the
//      mode itself (each drawing holds its own palette colors and padding
//      here, so a watcher that followed the shown drawing's would see its
//      source change), no watcher imports an auxiliary table's file (M3.12:
//      the file input's transition does), watchers.js writes no editor
//      member through a drawing (the slice write scan), and no watcher's
//      source or callback hands a slice handle to a call other than a slice
//      reader (the R2 slice-call scan); both scans are in
//      helpers/drawing-writes.mjs and see neither an alias held under another
//      name nor a write by a helper a callback calls;
//   2. the transition writes no drawing member and does not commit a Result
//      edit;
//   3. the Depth panels read the drawing's series and never write them (OV-109).
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import { test } from 'node:test';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};

const { state } = await import('../../gbdraw/web/js/state.js');
const { installShownDrawing } = await import('../../gbdraw/web/js/services/config.js');
const { createAuxiliaryFileImports, setupWatchers } = await import('../../gbdraw/web/js/app/watchers.js');
const { memberWriteSites, sliceHandleCalls } = await import('./helpers/drawing-writes.mjs');
const { DRAWING_SLICES, SLICE_READERS } = await import('./helpers/drawing-state.mjs');
// The calls in `source` that hand a slice handle to a writer.
const sliceWriterCalls = (source) => sliceHandleCalls(source, Object.keys(DRAWING_SLICES))
  .filter(({ callee }) => !SLICE_READERS.includes(callee)).map(({ callee }) => callee);
const APP_SETUP = readFileSync(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8');
const WATCHERS = readFileSync(new URL('../../gbdraw/web/js/app/watchers.js', import.meta.url), 'utf8');

// Any port of the watcher owner: a function that does nothing.
const port = () => new Proxy({}, { get: () => () => {} });
const snapshot = (value) => {
  if (value && typeof value === 'object') {
    try {
      return { identity: value, json: JSON.stringify(value) };
    } catch {
      return { identity: value, json: null };
    }
  }
  return { identity: value, json: null };
};
const sameSource = (before, after) => before.identity === after.identity && before.json === after.json;

/** @type {{ source: Function, callback: Function, line: string }[]} */
const watchers = [];
/** @type {string[]} */
const closedDialogs = [];
const { resetModeTransientUi, importAuxiliaryTable } = setupWatchers({
  state,
  rulePreparation: port(),
  trackAuxiliaryFileImport: createAuxiliaryFileImports((target) => target).track,
  ref: (value) => ({ value }),
  computed: (getter) => ({ get value() { return getter(); } }),
  watch: (source, callback) => {
    watchers.push({ source: typeof source === 'function' ? source : () => source, callback,
      line: String(source).replace(/\s+/g, ' ').slice(0, 600) });
    return () => {};
  },
  nextTick: async () => {},
  onMounted: () => {},
  legendActions: port(),
  featureActions: port(),
  legendLayout: port(),
  resultsManager: port(),
  runLabelReflow: () => {},
  refreshCircularRecordOrder: () => {},
  refreshLinearRecordSelectors: () => {},
  resetPreviewViewport: () => {},
  resetRightDrawer: () => {},
  closeLabelTextScopeDialog: () => { closedDialogs.push('label text scope'); },
  closeHiddenLabelTextDialog: () => { closedDialogs.push('Label Not Shown'); },
  clearLabelBuildNotices: () => {},
  previewRuntime: port()
});

test('a mode switch changes the source of no drawing watcher', () => {
  assert.ok(watchers.length > 0, 'the watchers were registered');
  assert.deepEqual(memberWriteSites(WATCHERS), [], 'the watchers write no editor member of a drawing');
  assert.deepEqual(watchers.flatMap(({ source, callback }) => sliceWriterCalls(`${source}\n${callback}`)), [],
    'the watchers call no slice writer');

  // Distinct values in each drawing, so a watcher that followed the shown
  // drawing would see its source change.
  modeDrawing(state, 'circular').currentColors.value = { CDS: '#111111' };
  modeDrawing(state, 'linear').currentColors.value = { CDS: '#222222' };
  modeDrawing(state, 'circular').canvasPadding.top = 5;
  modeDrawing(state, 'linear').canvasPadding.top = 9;

  for (const [from, to] of [['circular', 'linear'], ['linear', 'circular']]) {
    installShownDrawing(from);
    const before = watchers.map(({ source }) => snapshot(source()));
    installShownDrawing(to);
    const changed = watchers.filter((watcher, index) => !sameSource(before[index], snapshot(watcher.source())))
      .map(({ line }) => line);
    // Only watchers of the mode itself (and of values derived from it) see a switch.
    const unexpected = changed.filter((line) => !/\bmode\.value\b|state\.mode\b|mode\b\.value/.test(line));
    assert.deepEqual(unexpected, [], `${from} -> ${to}: these watchers would fire on the switch`);
  }
  installShownDrawing('circular');
});

test('no watcher stores the Multi-Record Canvas defaults, and no flag suppresses them (A2 W1)', () => {
  assert.deepEqual(watchers.filter(({ source }) => /multi_record_canvas/.test(String(source))).map(({ line }) => line), []);
  assert.ok(!('suppressCircularMultiRecordDefaults' in state), 'the suppress flag is gone');
});

test('the watcher scans flag a slice writer call and leave slice reads alone', () => {
  const [palette] = Object.keys(DRAWING_SLICES).filter((name) => name === 'palette');
  assert.deepEqual(sliceWriterCalls([
    `watch(() => state.paletteInstantPreviewEnabled.value, () => { clearQueuedPalette(state.activeDrawing().slices.${palette}); });`,
    `watch(() => state.mode.value, () => { const queued = queuedPaletteData(state.activeDrawing().slices.${palette}); });`,
    `watch(() => state.mode.value, () => { if (state.activeDrawing().slices.${palette}) return; });`
  ].join('\n')), ['clearQueuedPalette']);
});

test('an auxiliary table imports through its file input, into the drawing shown when it was chosen', async () => {
  assert.deepEqual(watchers.filter(({ line }) => /\bfiles(?:\[|\.(?:d_color|t_color|qualifier_priority|whitelist|blacklist)\b)/.test(line))
    .map(({ line }) => line), [], 'no watcher imports a table file');
  installShownDrawing('linear');
  const before = modeDrawing(state, 'circular').manualBlacklist.value;
  const imported = importAuxiliaryTable('blacklist', new File(['M3_12_WORD\n'], 'blacklist.txt', { type: 'text/plain' }));
  installShownDrawing('circular');
  assert.equal(await imported, true);
  assert.equal(state.files.blacklist?.name, 'blacklist.txt');
  assert.match(modeDrawing(state, 'linear').manualBlacklist.value, /M3_12_WORD/);
  assert.equal(modeDrawing(state, 'circular').manualBlacklist.value, before);
});

// OV-462: the label owner closes its dialogs, so a Label Not Shown dialog
// that is still open is answered (as Cancel) and not only hidden.
test('the mode switch closes the label dialogs through the label owner\'s ports', () => {
  const dialog = { ...state.hiddenLabelTextDialog };
  resetModeTransientUi();
  assert.deepEqual(closedDialogs, ['label text scope', 'Label Not Shown']);
  assert.deepEqual({ ...state.hiddenLabelTextDialog }, dialog, 'the watcher owner does not write the dialog');
});

// The body of a `const <name> = (...) => { ... };` in app-setup.js.
const functionBody = (name) => {
  const start = APP_SETUP.indexOf(`const ${name} = `);
  assert.ok(start >= 0, `${name} exists`);
  const open = APP_SETUP.indexOf('{', APP_SETUP.indexOf('=>', start));
  let depth = 0;
  for (let index = open; index < APP_SETUP.length; index += 1) {
    if (APP_SETUP[index] === '{') depth += 1;
    if (APP_SETUP[index] === '}' && --depth === 0) return APP_SETUP.slice(open, index + 1);
  }
  throw new Error(`${name} has no end`);
};
// A drawing write: an assignment or an array mutator through a drawing (a
// `drawing` parameter, `state.drawings.<id>` or `[...]`, `activeDrawing()`, or a
// local alias of one of these), or a whole-object writer.
const DRAWING_RECEIVER = String.raw`\bdrawing|\bstate\.drawings(?:\.\w+|\[[^\]]+\])|\bdrawings\.\w+|\bactiveDrawing\(\)`;
const DRAWING_ALIAS = /\b(?:const|let)\s+(\w+)\s*=\s*(?:state\.drawings(?:\.\w+|\[[^\]]+\])|(?:state\.)?activeDrawing\(\))/g;
const WHOLE_OBJECT_WRITE = /\b(replaceReactiveObject|replaceReactiveArray|clearReactiveObject|apply(?:Captured)?ConfigData)\(/;
/** @param {string} text */
const writesDrawing = (text) => {
  const receiver = `(?:${[DRAWING_RECEIVER, ...[...text.matchAll(DRAWING_ALIAS)].map(([, name]) => String.raw`\b${name}`)].join('|')})`;
  return new RegExp(String.raw`(?:${receiver})\.[\w.]+(\.value)?\s*(=(?!=)|\+=|-=)|(?:${receiver})\.[\w.]+\.(splice|push|pop|shift|unshift)\(`)
    .test(text) || WHOLE_OBJECT_WRITE.test(text);
};

test('the drawing-write scan sees each receiver form (A2 review-I2, Low 2)', () => {
  for (const planted of [
    'const drawing = state.drawings.circular; drawing.form.show_depth = true',
    'state.drawings.circular.adv.depth_tracks.splice(0)',
    'state.activeDrawing().form.show_depth = true',
    'const d = state.drawings.circular; d.form.show_depth = true',
    'state.drawings.circular.circularConservation.enabled = false',
    "state.drawings.linear.form.show_depth = false"
  ]) assert.equal(writesDrawing(`watch(() => x, () => { ${planted}; })`), true, planted);
  for (const read of ['const d = state.drawings.circular; if (d.form.show_depth === true) f(d)', 'drawing.form.show_depth !== x'])
    assert.equal(writesDrawing(read), false, read);
});

test('the drawing switch writes no drawing member and commits no Result edit', () => {
  const body = functionBody('selectDrawing');
  assert.equal(writesDrawing(body), false);
  assert.doesNotMatch(body, /commitActiveResultEdit|recordHistory|history\.(record|commit)/);
  // The removed steps stay removed: the Show Depth repair and the profile transition.
  assert.ok(!/modeProfileStateManager|\.transition\(drawing\.adv/.test(APP_SETUP),
    'app-setup.js runs no profile transition on a switch');
});

// The text of every `watch(...)` call in a module (balanced parentheses).
const watchCalls = (source) => [...source.matchAll(/\bwatch\(/g)].map(({ index }) => {
  let depth = 0;
  for (let end = index + 'watch'.length; end < source.length; end += 1) {
    if (source[end] === '(') depth += 1;
    if (source[end] === ')' && --depth === 0) return source.slice(index, end + 1);
  }
  throw new Error(`the watch call at ${index} has no end`);
});
const APP_DIR = new URL('../../gbdraw/web/js/app/', import.meta.url);
const appModules = readdirSync(APP_DIR, { recursive: true }).filter((file) => String(file).endsWith('.js'))
  .map((file) => [String(file), readFileSync(new URL(String(file), APP_DIR), 'utf8')]);
// The settings fits and drafts the former settings watchers wrote: each
// control, Reset, Load, Generate, and setup write them in their transitions
// (A2, D-15-46). A drawing write by a helper a callback calls is not seen,
// except through these names.
const SETTINGS_WRITES = /\b(ensureDepthTrackConfigCount|reconcileCircularConservation|syncCircularConservationSeries|refreshCircularConservationSlots|syncCircularConservationSlots|recordDisplayDrafts|layoutPreferences)\b/;

test('no watcher of an app module writes a drawing (A2: settings watchers)', () => {
  const calls = appModules.flatMap(([file, source]) => watchCalls(source).map((text) => ({ file, text })));
  assert.ok(calls.some(({ file }) => file === 'app-setup.js'), 'the watch calls were found');
  assert.deepEqual(calls.filter(({ text }) => SETTINGS_WRITES.test(text) || writesDrawing(text))
    .map(({ file, text }) => `${file}: ${text.slice(0, 160)}`), []);
});

test('the Depth panels read the drawing series and never write them (OV-109)', () => {
  for (const name of ['rowsForDepthTrackCount', 'circularDepthTrackRows', 'linearDepthTrackRows', 'depthTrackRows',
    'readDepthTrackConfig', 'getDepthTrackLabel', 'getDepthTrackColor']) {
    const body = functionBody(name);
    assert.doesNotMatch(body, /ensureDepthTrack|depth_tracks\s*=(?!=)|depth_tracks\.(splice|push)\(/, `${name} writes the series`);
  }
});
