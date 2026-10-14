// Session admission vectors shared with Python: what the Web's whole Load
// (`importSession`) decides. tests/fixtures/session_admission_vectors.json is
// also read by tests/test_session_admission_vectors.py (`load_session`, and the
// rules S1-S3 between the two columns). Each vector edits a copy of its base
// Session; the Web's verdict is `loaded`, `repaired` with the code@path of each
// row of the Load dialog (`webRepairs`), or `refused`.
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
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

const { importSession } = await import('../../gbdraw/web/js/services/config.js');
// The composition root's transform of an older Session's Results (R13 port).
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');

const REPO = new URL('../../', import.meta.url);
const FIXTURE = JSON.parse(readFileSync(new URL('tests/fixtures/session_admission_vectors.json', REPO), 'utf8'));

const base = (name) => {
  const entry = FIXTURE.bases[name];
  if (entry.document) return structuredClone(entry.document);
  const bytes = readFileSync(new URL(entry.fixture, REPO));
  return JSON.parse((bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8'));
};
const edited = (vector) => {
  const document = base(vector.base);
  for (const edit of vector.edits) {
    const tokens = edit.path.split('/').slice(1).map((raw) => raw.replace(/~1/g, '/').replace(/~0/g, '~'));
    const last = tokens.pop();
    const parent = tokens.reduce((current, token) => current[token], document);
    if (edit.op === 'set') parent[last] = structuredClone(edit.value); // a later edit may write inside it
    else {
      assert.equal(edit.op, 'remove', vector.id);
      delete parent[last];
    }
  }
  return document;
};

// The diagram Worker that recovers the Features of a Session 27-39 does not run
// in node; its row is the environment's, as in the G3 corpus test.
const NODE_ONLY_REPAIRS = new Set(['feature-recovery-failed']);

for (const vector of FIXTURE.vectors) {
  test(`${vector.id}: the Web's Load: ${vector.web} (${vector.finding})`, async () => {
    const outcome = await importSession({
      target: { files: [new Blob([JSON.stringify(edited(vector))], { type: 'application/json' })], value: 'selected' }
    }, { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    const repairs = (outcome?.repairs ?? []).filter(({ code }) => !NODE_ONLY_REPAIRS.has(code)).map(({ code, path }) => `${code}@${path}`);
    const verdict = outcome?.status !== 'ok' ? 'refused' : repairs.length > 0 ? 'repaired' : 'loaded';
    assert.equal(verdict, vector.web, JSON.stringify(outcome?.error ?? repairs));
    assert.deepEqual(repairs, vector.webRepairs ?? []);
  });
}
