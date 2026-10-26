// A Session 49 Load keeps every drawing's request metadata (`webFiles` beside
// its bindings) from its own entry, never from the drawing the Load opened, so
// a re-Save writes it back (OV-1872).
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { File } from 'node:buffer';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
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
globalThis.document = { body: { appendChild: () => {} }, createElement: () => ({ click: () => {}, addEventListener: () => {}, parentNode: null }) };
installFakeSvgDom();
globalThis.File = File;
globalThis.alert = () => {};
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const v49Session = () => JSON.parse(gunzipSync(readFileSync(new URL(
  '../fixtures/sessions/drawings-v49/drawings-circular-linear.v49.gbdraw-session.json.gz', import.meta.url
))));
// The Load hands its caller the artifact slot of each drawing it does not
// show (app-setup keeps them; Save writes each drawing's request metadata from
// its slot's committed request).
const loadKeepingSlots = async (session) => {
  const stashed = [];
  const loaded = await importSession({
    target: { files: [new File([JSON.stringify(session)], 'session.json', { type: 'application/json' })], value: 'selected' }
  }, { confirmLoadRepairs: async () => true, installLoadedArtifactSlots: (slots) => { stashed.push(...slots.stashed); } });
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  return stashed;
};
/** `webFiles` without its bindings (Save writes the bindings from the drawing's inputs). */
const requestMetadata = (webFiles) => {
  const { bindings: _bindings, resourceOriginalNames: _names, ...metadata } = webFiles ?? {};
  return metadata;
};

for (const opened of ['linear', 'circular']) {
  test(`a Session 49 Load keeps each drawing's own request metadata in the drawing it does not open (${opened} opened)`, async () => {
    const session = v49Session();
    session.ui.activeDrawingId = opened;
    const stashed = await loadKeepingSlots(session);
    assert.equal(state.activeDrawing().id, opened);
    const other = session.drawings.find((entry) => entry.id !== opened);
    const slot = stashed.find((entry) => entry.drawingId === other.id);
    assert.ok(slot, JSON.stringify(stashed.map((entry) => entry.drawingId)));
    assert.deepEqual(requestMetadata(slot.runtimeState.canonical.committedCanonicalSession.webFiles), requestMetadata(other.webFiles));
  });
}
