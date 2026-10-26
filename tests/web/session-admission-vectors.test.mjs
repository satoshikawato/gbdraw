// Session admission vectors shared with Python: what the Web's whole Load
// (`importSession`) decides. tests/fixtures/session_admission_vectors.json is
// also read by tests/test_session_admission_vectors.py (`load_session`, and the
// rules S1-S3 between the two columns). Each vector edits a copy of its base
// Session; the Web's verdict is `loaded`, `repaired` with the code@path of each
// repair, listed in the Load dialog or not (`webRepairs`), or `refused`.
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { register } from 'node:module';
import test from 'node:test';
import { gunzipSync } from 'node:zlib';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';

// R4 has one owner (E-OD22): an R4 vector passes only when `restoresNothing`
// decided the refusal, not another producer of the same diagnostic.
register('./helpers/restores-nothing-spy.mjs', import.meta.url);

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

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
// A saved unmanaged override is checked as the browser's Load checks it (PD-OI-009).
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
// A Session 27-44 moves to a Session 49 in Python, as the composition root's
// Worker port does it (DESIGN 4.3).
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
setLegacySessionUpgrader(pythonLegacySessionUpgrader());
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
// in node; its row (at `features`) is the environment's, as in the G3 corpus test.
const nodeOnly = ({ code, path }) => code === 'feature-recovery-failed' && path === 'features';
// The diagnostic reason of each refusal class of the Load repair receipt (S2).
const CLASS_REASONS = { R2: 'SESSION_FORMAT', R3: 'SESSION_TOO_OLD', R4: 'FIELDS' };
const diagnosticText = (error) => (error?.details ?? []).map(({ text }) => text).join('\n');

for (const vector of FIXTURE.vectors) {
  test(`${vector.id}: the Web's Load: ${vector.web} (${vector.finding})`, async () => {
    const verdicts = [];
    globalThis.restoresNothingVerdicts = verdicts;
    const outcome = await importSession({
      target: { files: [new Blob([JSON.stringify(edited(vector))], { type: 'application/json' })], value: 'selected' }
    }, { transformLegacyResultSvg, confirmLoadRepairs: async () => true });
    const repairs = (outcome?.repairs ?? []).filter((repair) => !nodeOnly(repair)).map(({ code, path }) => `${code}@${path}`);
    const verdict = outcome?.status !== 'ok' ? 'refused' : repairs.length > 0 ? 'repaired' : 'loaded';
    assert.equal(verdict, vector.web, JSON.stringify(outcome?.error ?? repairs));
    assert.deepEqual(repairs, vector.webRepairs ?? []);
    if (verdict === 'refused') {
      assert.equal(outcome.error.code, 'INPUT_INVALID', vector.id);
      assert.match(diagnosticText(outcome.error), new RegExp(`reason: ${CLASS_REASONS[vector.webClass]}$`, 'm'), vector.id);
    }
    // R4's verdict is the last one Load asked for; no other class asks it to true.
    assert.equal(verdicts.at(-1) === true, vector.webClass === 'R4', `${vector.id}: restoresNothing verdicts ${JSON.stringify(verdicts)}`);
  });
}
