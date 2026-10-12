// Session admission vectors shared with Python: what the Web's whole Load
// (`importSession`) decides. tests/fixtures/session_admission_vectors.json is
// also read by tests/test_session_admission_vectors.py (`load_session`). Each
// vector edits a copy of its base Session and states `accept` or `refuse` per
// side; a vector whose sides disagree names its open finding in `difference`.
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
    if (edit.op === 'set') parent[last] = edit.value;
    else {
      assert.equal(edit.op, 'remove', vector.id);
      delete parent[last];
    }
  }
  return document;
};

for (const vector of FIXTURE.vectors) {
  test(`${vector.id}: the Web's Load ${vector.web}s it (${vector.finding})`, async () => {
    const outcome = await importSession({
      target: { files: [new Blob([JSON.stringify(edited(vector))], { type: 'application/json' })], value: 'selected' }
    }, { transformLegacyResultSvg });
    assert.equal(outcome?.status, vector.web === 'accept' ? 'ok' : 'error', JSON.stringify(outcome?.error));
  });
}
