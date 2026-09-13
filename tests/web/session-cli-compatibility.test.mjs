import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtemp, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { gunzipSync } from 'node:zlib';
import { test } from 'node:test';
import { createHash } from 'node:crypto';
import { installFakeSvgDom } from './fake-svg-dom.mjs';

globalThis.window = {
  Vue: {
    ref: value => ({ value }), reactive: value => value,
    computed: getter => ({ get value() { return getter(); } }), nextTick: async () => {}
  },
  DOMPurify: { sanitize: value => value }
};
globalThis.document = {};
installFakeSvgDom();
globalThis.File = class File extends Blob {
  constructor(parts, name, options = {}) {
    super(parts, options);
    this.name = name;
    this.lastModified = options.lastModified ?? 0;
  }
};
globalThis.alert = () => {};

const { importSession, getCommittedCanonicalRenderRequest, setUnmanagedConfigOverrideValidator } = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { getSessionResourceSource, readFileBytes } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const root = process.cwd();
// Exercise the Worker's actual typed helper without starting a browser runtime.
setUnmanagedConfigOverrideValidator(payload => ({ result: JSON.parse(execFileSync('python', ['-c', `
import json, sys
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
p = json.load(sys.stdin)
print(validate_web_config_overrides_json(p['mode'], json.dumps(p['config']),
    json.dumps(p['configOverrides']), json.dumps(p['managedPaths']), p['requireUnmanagedOnly']))
`], { input: JSON.stringify(payload), encoding: 'utf8', cwd: root })) }));
const mito = path.join(root, 'tests/fixtures/sessions/cli-web-mito.gb');
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const lambda = path.join(root, 'tests/test_inputs/NC_001416.gb');
const cases = [
  ['single', 'circular', ['--gbk', mito, '--labels', 'out'], [mito]],
  ['composite', 'circular', ['--gbk', mito, lambda, '--multi_record_canvas'], [mito, lambda]],
  ['linear', 'linear', ['--gbk', mito, lambda], [mito, lambda]],
  ['gff', 'circular', ['--gff', path.join(root, 'tests/test_inputs/NC_013668.gff3'),
    '--fasta', path.join(root, 'tests/test_inputs/NC_013668.fasta')], []]
];
const load = bytes => importSession({ target: {
  files: [new File([bytes], 'current.json', { type: 'application/json' })], value: 'selected'
} });

for (const [label, mode, args, sourcePaths] of cases) {
  await test(`current CLI ${label} and CLI replay initialize Web from the canonical request`, async () => {
    const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-cli-web-'));
    try {
      let input = args;
      for (const phase of ['fresh', 'replay']) {
        const output = path.join(directory, phase);
        const file = `${output}.gbdraw-session.json.gz`;
        execFileSync('python', ['-m', 'gbdraw.cli', mode, ...input, '-o', output, '--session_output', file], {
          cwd: directory, env: { ...process.env, PYTHONPATH: root }, stdio: 'pipe', timeout: 1_800_000
        });
        const bytes = gunzipSync(await readFile(file));
        const session = JSON.parse(bytes);
        assert.equal(session.version, 42);
        assert.equal(session.renderRequest.schema, 7);
        assert.equal(Object.hasOwn(session, 'config'), false);
        assert.equal(session.webFiles.bindings.schema, 2);
        const result = await load(bytes);
        assert.equal(result.status, 'ok', result.error?.stack);
        assert.equal(state.mode.value, mode);
        assert.equal((mode === 'circular' ? state.cInputType : state.lInputType).value, label === 'gff' ? 'gff' : 'gb');
        assert.equal(state.importedComparisonIntent.disposition, 'EDITABLE');
        assert.deepEqual(getCommittedCanonicalRenderRequest(), session.renderRequest);
        assert.ok(state.results.value.length > 0);
        if (label === 'single') assert.equal(state.form.labels_mode, 'out');
        const projectedFiles = mode === 'linear' ? state.linearSeqs.map(seq => seq.gb)
          : label === 'gff' ? [] : [state.files.c_gb];
        const sources = projectedFiles.filter(Boolean).flatMap(file => {
          const source = getSessionResourceSource(file);
          return source.descriptors || [source];
        });
        assert.equal(sources.length, sourcePaths.length);
        for (const [index, source] of sources.entries()) {
          assert.equal(hash(Buffer.from(source.descriptor.data, 'base64')), hash(await readFile(sourcePaths[index])));
          assert.ok(session.resources[source.resourceId]);
        }
        if (label === 'gff') {
          assert.equal(hash(await readFileBytes(state.files.c_gff)), hash(await readFile(args[1])));
          assert.equal(hash(await readFileBytes(state.files.c_fasta)), hash(await readFile(args[3])));
        }
        // The original malformed adjunct is still invalid, even beside a valid request.
        const before = getCommittedCanonicalRenderRequest();
        const errorLog = console.error;
        console.error = () => {};
        try {
          const rejected = await load(JSON.stringify({ ...session, config: { adv: {} } }));
          assert.equal(rejected.status, 'error');
          assert.match(rejected.error.message, /missing its active form or advanced settings/);
          assert.deepEqual(getCommittedCanonicalRenderRequest(), before);
        } finally {
          console.error = errorLog;
        }
        input = ['--session', file];
      }
    } finally {
      await rm(directory, { recursive: true, force: true });
    }
  });
}
