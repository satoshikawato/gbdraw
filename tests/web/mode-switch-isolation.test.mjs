// A mode switch changes no drawing (plan OV-80 §4.5, §8 "Mode switch"; PR-0b
// open items). Each mode has its own drawing, so the switch sets `mode`, swaps
// the artifact slots, and writes no setting or edit:
//   1. no watcher's source changes on a switch except the ones that watch the
//      mode itself, so the drawing watchers (palette colors, Canvas padding,
//      bulk label text) never take a switch for an edit;
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
const { setupWatchers } = await import('../../gbdraw/web/js/app/watchers.js');
const APP_SETUP = readFileSync(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8');

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

test('a mode switch changes the source of no drawing watcher', () => {
  /** @type {{ source: Function, callback: Function, line: string }[]} */
  const watchers = [];
  setupWatchers({
    state,
    rulePreparation: port(),
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
    closeLabelTextScopeDialog: () => {},
    clearLabelBuildNotices: () => {},
    previewRuntime: port()
  });
  assert.ok(watchers.length > 10, 'the watchers were registered');
  const drawingWatchers = watchers.filter(({ line }) => /drawing\./.test(line));
  assert.ok(drawingWatchers.length >= 6, 'each drawing has its own palette, padding, and label text watchers');

  // Distinct values in each drawing, so a watcher that followed the shown
  // drawing would see its source change.
  state.drawings.circular.currentColors.value = { CDS: '#111111' };
  state.drawings.linear.currentColors.value = { CDS: '#222222' };
  state.drawings.circular.canvasPadding.top = 5;
  state.drawings.linear.canvasPadding.top = 9;
  state.drawings.circular.labelTextBulkOverrides.hypothetical = 'HP';
  state.drawings.linear.canonicalLabelOverrideRows.value = [{ kept: true }];

  for (const [from, to] of [['circular', 'linear'], ['linear', 'circular']]) {
    state.mode.value = from;
    const before = watchers.map(({ source }) => snapshot(source()));
    state.mode.value = to;
    const changed = watchers.filter((watcher, index) => !sameSource(before[index], snapshot(watcher.source())))
      .map(({ line }) => line);
    // Only watchers of the mode itself (and of values derived from it) see a switch.
    const unexpected = changed.filter((line) => !/\bmode\.value\b|state\.mode\b|mode\b\.value/.test(line));
    assert.deepEqual(unexpected, [], `${from} -> ${to}: these watchers would fire on the switch`);
  }
  // The bulk label text watcher of the Linear drawing did not run, so the
  // Linear drawing keeps its saved label table.
  assert.deepEqual(state.drawings.linear.canonicalLabelOverrideRows.value, [{ kept: true }]);
  state.mode.value = 'circular';
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
const DRAWING_WRITE = /\bdrawing\.[\w.]+(\.value)?\s*(=(?!=)|\+=|-=)|\bdrawing\.[\w.]+\.(splice|push|pop|shift|unshift)\(|\b(replaceReactiveObject|replaceReactiveArray|clearReactiveObject|applyConfigData)\(/;

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
