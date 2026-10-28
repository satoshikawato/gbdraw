// A Session 49 Load reads each drawing's parts once, at the path the file
// holds them (`drawings[i]`): the opened drawing's editor state and LOSAT
// caches as every other drawing's (OV-1791, OV-1797).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { createHistoryRouter } from '../../gbdraw/web/js/services/history.js';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({ addEventListener: () => {}, click: () => {}, parentNode: null })
};
globalThis.alert = () => {};
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const v49Session = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}.v49.gbdraw-session.json.gz`, import.meta.url))));
const load = async (session) => {
  const loaded = await importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, { confirmLoadRepairs: async () => true });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  return loaded;
};
const rowsAt = (loaded, part) => loaded.repairs.filter((repair) => repair.path.includes(part))
  .map((repair) => `${repair.code} ${repair.path}`).sort();

// OV-2131: a Load whose drawings the page does not hold (a second Circular
// drawing, on a page that starts with one Circular drawing) installs them with
// their tabs and gives each a History at the store commit, so a reader of the
// shown drawing's History before the Load's baseline finds it.
test('a Session 49 Load of two Circular drawings and a Linear one into one Circular drawing installs every drawing, its tab and its History', async () => {
  assert.deepEqual(state.drawingTabs.value.map((tab) => tab.id), ['circular'], 'the page starts with one Circular drawing');
  const router = createHistoryRouter({
    drawingIds: state.drawingTabs.value.map((tab) => tab.id),
    activeDrawingId: () => state.activeDrawingId.value,
    maxBytes: Number.MAX_SAFE_INTEGER,
    sessionHistory: { buildIntent: () => ({}), applyIntent: () => {} },
    estimateFileBytes: () => 0,
    retainOnlyFiles: () => {},
    managerOptions: () => ({ buildIntent: () => ({}), applyIntent: () => {}, buildCheckpoint: () => ({}), applyCheckpoint: () => {} })
  });
  const routed = callerSession(sessionConfig, state, await sessionRollbackPort(state, router));
  const session = v49Session('drawings-circular-linear');
  const [circular, linear] = session.drawings;
  session.drawings = [circular, { ...structuredClone(circular), id: 'circular-2', name: 'Plasmid map' }, linear];
  session.ui.activeDrawingId = 'circular-2';
  /** @type {string[]} */
  const historyAtAfterImport = [];
  const loaded = await routed.importSession({
    target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
  }, {
    confirmLoadRepairs: async () => true,
    afterImport: async () => {
      router.canUndo();
      for (const tab of state.drawingTabs.value) historyAtAfterImport.push(router.managerOf(tab.id) ? tab.id : '');
    }
  });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  const stored = session.drawings.map(({ id, name }) => ({ id, name }));
  assert.deepEqual(state.drawingTabs.value.map(({ id, name }) => ({ id, name })), stored);
  assert.deepEqual(state.drawingList().map((drawing) => drawing.id), stored.map(({ id }) => id), 'every tab names a drawing');
  assert.equal(state.activeDrawing().id, 'circular-2', 'the stored shown drawing is shown');
  assert.deepEqual(historyAtAfterImport, stored.map(({ id }) => id), 'each loaded drawing has a History before the baseline');
});

for (const opened of ['circular', 'linear']) {
  test(`a Session 49 Load lists the opened drawing's editor-state fault once, at its drawings[i] path (${opened} opened)`, async () => {
    const session = v49Session('drawings-circular-linear');
    session.ui.activeDrawingId = opened;
    const index = session.drawings.findIndex((entry) => entry.id === opened);
    session.drawings[index].editorState.originalSvgStroke = { color: 'not-a-color' };
    const loaded = await load(session);
    assert.deepEqual(rowsAt(loaded, 'originalSvgStroke'),
      [`setting-replaced drawings[${index}].editorState.originalSvgStroke.color`]);
  });

  test(`a Session 49 Load lists the opened drawing's dropped LOSAT entry once, at its drawings[i] path (${opened} opened)`, async () => {
    const session = v49Session('drawings-circular-linear');
    session.ui.activeDrawingId = opened;
    const index = session.drawings.findIndex((entry) => entry.id === opened);
    session.drawings[index].losatCache.entries = [{ no: 'entry' }];
    const loaded = await load(session);
    assert.deepEqual(rowsAt(loaded, 'losatCache'), [`cache-entry-dropped drawings[${index}].losatCache.entries[0]`]);
  });
}

// N1: the drawings Load does not open read their LOSAT caches through the
// one reader too, each at its own path.
test('a Session 49 Load reads every drawing\'s LOSAT caches, each at its drawings[i] path', async () => {
  const session = v49Session('drawings-circular-linear');
  session.drawings.forEach((entry) => {
    entry.losatCache.entries = [{ no: 'entry' }];
    entry.losatDerivedCache = [];
  });
  const loaded = await load(session);
  assert.deepEqual([...rowsAt(loaded, 'losatCache'), ...rowsAt(loaded, 'losatDerivedCache')].sort(), session.drawings
    .flatMap((_, index) => [`cache-entry-dropped drawings[${index}].losatCache.entries[0]`, `cache-entry-dropped drawings[${index}].losatDerivedCache`])
    .sort());
});

// N2: each loaded drawing installs the inputs paired with it; no drawing
// falls back to the shown drawing's inputs when its own are missing.
test('the Load install of each drawing\'s inputs takes them paired with the drawing, with no fallback', () => {
  const source = readFileSync(new URL('../../gbdraw/web/js/services/config.js', import.meta.url), 'utf8');
  const start = source.indexOf('const restoreLiveFileState = ');
  assert.ok(start >= 0);
  const body = source.slice(start, source.indexOf('\n};\n', start));
  assert.deepEqual(body.match(/\.get\([^)]*\)\s*\?\?[^;\n]*/g) ?? [], []);
});

// OV-1792: a legacy configuration Load returns before `afterImport`, which
// only a Session Load reaches, so it holds no legacy arm.
test('the Load afterImport hook holds no legacy-configuration arm', () => {
  const source = readFileSync(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8');
  const start = source.indexOf('afterImport: async');
  assert.ok(start >= 0);
  assert.equal(source.slice(start, source.indexOf('\n  });\n', start)).includes("'legacy'"), false);
});
