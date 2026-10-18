// A mode switch changes no drawing (plan OV-80 §4.5, §8 "Mode switch"; PR-0b
// open items). Each mode has its own drawing, so the switch sets `mode`, swaps
// the artifact slots, and writes no setting or edit:
//   1. no watcher's source changes on a switch except the ones that watch the
//      mode itself (each drawing holds its own palette colors and padding
//      here, so a watcher that followed the shown drawing's would see its
//      source change), the watchers write no editor member (the slice write
//      scan, helpers/drawing-writes.mjs), and no watcher imports an auxiliary
//      table's file (M3.12: the file input's transition does);
//   2. the transition writes no drawing member and does not commit a Result
//      edit;
//   3. the Depth panels read the drawing's series and never write them (OV-109).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

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
const { memberWriteSites } = await import('./helpers/drawing-writes.mjs');
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

  // Distinct values in each drawing, so a watcher that followed the shown
  // drawing would see its source change.
  state.drawings.circular.currentColors.value = { CDS: '#111111' };
  state.drawings.linear.currentColors.value = { CDS: '#222222' };
  state.drawings.circular.canvasPadding.top = 5;
  state.drawings.linear.canvasPadding.top = 9;

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

test('an auxiliary table imports through its file input, into the drawing shown when it was chosen', async () => {
  assert.deepEqual(watchers.filter(({ line }) => /\bfiles(?:\[|\.(?:d_color|t_color|qualifier_priority|whitelist|blacklist)\b)/.test(line))
    .map(({ line }) => line), [], 'no watcher imports a table file');
  installShownDrawing('linear');
  const before = state.drawings.circular.manualBlacklist.value;
  const imported = importAuxiliaryTable('blacklist', new File(['M3_12_WORD\n'], 'blacklist.txt', { type: 'text/plain' }));
  installShownDrawing('circular');
  assert.equal(await imported, true);
  assert.equal(state.files.blacklist?.name, 'blacklist.txt');
  assert.match(state.drawings.linear.manualBlacklist.value, /M3_12_WORD/);
  assert.equal(state.drawings.circular.manualBlacklist.value, before);
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
const DRAWING_WRITE = /\bdrawing\.[\w.]+(\.value)?\s*(=(?!=)|\+=|-=)|\bdrawing\.[\w.]+\.(splice|push|pop|shift|unshift)\(|\b(replaceReactiveObject|replaceReactiveArray|clearReactiveObject|apply(?:Captured)?ConfigData)\(/;

test('the mode transition writes no drawing member and commits no Result edit', () => {
  const body = functionBody('transitionDiagramMode');
  assert.doesNotMatch(body, DRAWING_WRITE);
  assert.doesNotMatch(body, /commitActiveResultEdit|recordHistory|history\.(record|commit)/);
  // The removed steps stay removed: the Show Depth repair and the profile transition.
  assert.ok(!/modeProfileStateManager|\.transition\(drawing\.adv/.test(APP_SETUP),
    'app-setup.js runs no profile transition on a switch');
});

test('the Depth panels read the drawing series and never write them (OV-109)', () => {
  for (const name of ['rowsForDepthTrackCount', 'circularDepthTrackRows', 'linearDepthTrackRows', 'depthTrackRows',
    'readDepthTrackConfig', 'getDepthTrackLabel', 'getDepthTrackColor']) {
    const body = functionBody(name);
    assert.doesNotMatch(body, /ensureDepthTrack|depth_tracks\s*=(?!=)|depth_tracks\.(splice|push)\(/, `${name} writes the series`);
  }
});
