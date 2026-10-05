import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtemp, readFile, rm, writeFile } from 'node:fs/promises';
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

installSessionImportWorker();

const {
  SESSION_VERSION, importSession, getCommittedCanonicalRenderRequest, getCommittedCanonicalSession,
  serializeActiveRenderFiles, setUnmanagedConfigOverrideValidator
} = await import('../../gbdraw/web/js/services/config.js');
const { CANONICAL_REQUEST_SCHEMA, buildCanonicalRenderRequest } = await import('../../gbdraw/web/js/services/session-request.js');
const { inheritCommittedComparisonIntent } = await import('../../gbdraw/web/js/services/imported-comparison-intent.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/app/linear-comparisons.js');
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
// A non-default --legend checks that the committed legend survives load (SE-07).
const cases = [
  ['single', 'circular', ['--gbk', mito, '--labels', 'out', '--legend', 'upper_left'], [mito], 'upper_left'],
  ['composite', 'circular', ['--gbk', mito, lambda, '--multi_record_canvas', '--legend', 'upper_right'],
    [mito, lambda], 'upper_right'],
  ['linear', 'linear', ['--gbk', mito, lambda, '--legend', 'left'], [mito, lambda], 'left'],
  ['gff', 'circular', ['--gff', path.join(root, 'tests/test_inputs/NC_013668.gff3'),
    '--fasta', path.join(root, 'tests/test_inputs/NC_013668.fasta')], [], 'right']
];
// The committed request owns Linear record identity; a CLI binding uid is only
// an initial value (SE-06).
const fileRecordKeys = request => [...new Set(request.records.map(
  record => record.recordKey.replace(/:[1-9]\d*$/, '')
))];
const load = bytes => importSession({ target: {
  files: [new File([bytes], 'current.json', { type: 'application/json' })], value: 'selected'
} });
const DEFAULT_PLAN = { mode: 'none', defaultSource: 'losat', edges: [] };
// B15: a read-only CLI comparison (-b) is not a Web draft. The replacement draft
// is the Web default (No comparison), so Replace with current controls
// (app-setup.js: a valid plan with comparison intent) waits until the user sets
// up a comparison, and no LOSAT run starts from the draft.
const assertNoReplacementDraft = () => {
  assert.deepEqual(state.linearComparisonPlan, DEFAULT_PLAN);
  const draft = resolveLinearComparisonPlan({
    plan: state.linearComparisonPlan, sequences: state.linearSeqs, layout: [],
    losatProgram: state.losatProgram.value, blastpMode: state.losat.blastp.mode
  });
  assert.equal(draft.valid && draft.hasComparisonIntent, false);
  assert.equal(draft.hasLosatIntent, false);
};

for (const [label, mode, args, sourcePaths, legend] of cases) {
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
        assert.equal(session.version, SESSION_VERSION);
        assert.equal(session.renderRequest.schema, CANONICAL_REQUEST_SCHEMA);
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
        assert.equal(session.renderRequest.diagramOptions.output.legend, legend);
        assert.equal(state.form.legend, legend);
        if (mode === 'linear') {
          assert.deepEqual(state.linearSeqs.map(seq => seq.uid), fileRecordKeys(session.renderRequest));
          // B13: without -b the CLI commits only a disabled protein pipeline
          // (mode `none`, no pairs). The Web draft has no comparison, so the
          // request Generate builds has no comparison and starts no LOSAT run.
          assert.ok(session.renderRequest.comparisons.every(item => (
            item.kind === 'generatedProteinComparison' && item.mode === 'none' && item.pairs.length === 0
          )));
          assert.deepEqual(state.linearComparisonPlan, { mode: 'none', defaultSource: 'losat', edges: [] });
          assert.notEqual(state.losatProgram.value, 'blastp');
          const filesData = await serializeActiveRenderFiles('linear', state);
          const comparisonPlanSnapshot = resolveLinearComparisonPlan({
            plan: state.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
            losatProgram: state.losatProgram.value, blastpMode: state.losat.blastp.mode
          });
          assert.equal(comparisonPlanSnapshot.hasLosatIntent, false);
          const candidate = buildCanonicalRenderRequest({ state, filesData, comparisonPlanSnapshot });
          assert.deepEqual(candidate.renderRequest.comparisons, []);
        }
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
          assert.equal(rejected.error.code, 'INPUT_INVALID');
          assert.deepEqual(rejected.error.context, {field:'config',reason:'FIELDS'});
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

// A v42 CLI Linear BLAST sidecar written by first-parent main (provenance in
// tests/fixtures/sessions/se06-main-linear-blast-cli.provenance.json) keeps its
// comparison read-only and gives each Linear file its committed record identity.
await test('a main v42 CLI Linear BLAST sidecar keeps a read-only comparison with committed record keys', async () => {
  const bytes = gunzipSync(await readFile(path.join(
    root, 'tests/fixtures/sessions/se06-main-linear-blast-cli.v42.gbdraw-session.json.gz'
  )));
  const session = JSON.parse(bytes);
  assert.equal(session.version, 42);
  assert.equal(session.renderRequest.schema, 7);
  assert.deepEqual(session.webFiles.bindings.linearSeqs.map(seq => seq.uid), ['cli-seq-1', 'cli-seq-2']);
  const result = await load(bytes);
  assert.equal(result.status, 'ok', result.error?.stack);
  assert.equal(state.importedComparisonIntent.disposition, 'PRESERVED_READ_ONLY');
  assert.deepEqual(state.linearSeqs.map(seq => seq.uid), ['record-1', 'record-2']);
  assert.deepEqual(fileRecordKeys(getCommittedCanonicalRenderRequest()), ['record-1', 'record-2']);
  assertNoReplacementDraft();
});

// B3: one multi-record GenBank file with -b gives records `record-1:1` and
// `record-1:2`. Each record becomes its own Linear row keyed by that recordKey
// with the request's record selector, so Inherit binds the comparison to the
// named records.
await test('a multi-record CLI Linear BLAST Session inherits its comparison onto the named records', async () => {
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-cli-web-'));
  try {
    const blast = path.join(directory, 'R2c_R3c.tsv');
    await writeFile(blast, 'R2c\tR3c\t100.000\t1000\t0\t0\t2001\t3000\t1\t1000\t0.0\t1847\n');
    const file = path.join(directory, 'multi.gbdraw-session.json.gz');
    execFileSync('python', ['-m', 'gbdraw.cli', 'linear',
      '--gbk', path.join(root, 'tests/fixtures/web_comparison_shared_block.gb'), '-b', blast,
      '-o', path.join(directory, 'multi'), '--session_output', file], {
      cwd: directory, env: { ...process.env, PYTHONPATH: root }, stdio: 'pipe', timeout: 1_800_000
    });
    const bytes = gunzipSync(await readFile(file));
    const session = JSON.parse(bytes);
    assert.deepEqual(session.renderRequest.records.map(record => record.recordKey), ['record-1:1', 'record-1:2']);
    assert.equal(session.webFiles.bindings.linearSeqs.length, 1);
    const result = await load(bytes);
    assert.equal(result.status, 'ok', result.error?.stack);
    assert.equal(state.importedComparisonIntent.disposition, 'PRESERVED_READ_ONLY');
    assertNoReplacementDraft();
    // The candidate Generate builds after Inherit (run-analysis.js: empty plan, committed comparison).
    const filesData = await serializeActiveRenderFiles('linear', state);
    const candidate = buildCanonicalRenderRequest({
      state,
      filesData: { ...filesData, linearCanonicalComparisons: [] },
      comparisonPlanSnapshot: resolveLinearComparisonPlan({
        plan: { mode: 'none', defaultSource: 'losat', edges: [] },
        sequences: filesData.linearSeqs, layout: [], losatProgram: 'blastn', blastpMode: 'orthogroup'
      })
    });
    inheritCommittedComparisonIntent({ candidate, committed: getCommittedCanonicalSession() });
    const records = candidate.renderRequest.records;
    const [comparison] = candidate.renderRequest.comparisons.filter(item => item.kind === 'nucleotideBlast');
    assert.deepEqual([records[comparison.queryRecordIndex].recordKey, records[comparison.subjectRecordIndex].recordKey],
      ['record-1:1', 'record-1:2']);
    assert.deepEqual(records.map(record => [record.cardinality, record.selector]),
      session.renderRequest.records.map(record => ['exactly_one', record.selector]));
    assert.deepEqual(state.linearSeqs.map(seq => [seq.uid, seq.region_record_id]), [['record-1:1', 'R2c'], ['record-1:2', 'R3c']]);
    const [first, second] = state.linearSeqs.map(seq => getSessionResourceSource(seq.gb));
    assert.equal(first.resourceId, second.resourceId);
    assert.equal(hash(Buffer.from(first.descriptor.data, 'base64')),
      hash(await readFile(path.join(root, 'tests/fixtures/web_comparison_shared_block.gb'))));
    // A CLI Session that stored a copy of each drawn record has no selectors;
    // each row takes the `#n` of its expanded recordKey.
    const copies = { ...session, renderRequest: { ...session.renderRequest,
      records: session.renderRequest.records.map(record => ({ ...record, selector: null })) } };
    assert.equal((await load(JSON.stringify(copies))).status, 'ok');
    assert.deepEqual(state.linearSeqs.map(seq => [seq.uid, seq.region_record_id]), [['record-1:1', '#1'], ['record-1:2', '#2']]);
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
});

// B15: an editable CLI protein pipeline is the adjacent LOSATP comparison the
// CLI drew. The projection states that plan (the plan normalizer no longer
// falls back to adjacent), so Generate rebuilds the same pipeline.
await test('a CLI Linear protein Session keeps the adjacent LOSATP plan it drew', async () => {
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-cli-web-'));
  try {
    const file = path.join(directory, 'protein.gbdraw-session.json.gz');
    execFileSync('python', ['-m', 'gbdraw.cli', 'linear', '--gbk', mito, lambda,
      '--losat', 'losatp', '--losatp_mode', 'similarity_groups', '-o', path.join(directory, 'protein'), '--session_output', file], {
      cwd: directory, env: { ...process.env, PYTHONPATH: root }, stdio: 'pipe', timeout: 1_800_000
    });
    const bytes = gunzipSync(await readFile(file));
    const session = JSON.parse(bytes);
    assert.deepEqual(session.renderRequest.comparisons.map(item => [item.kind, item.mode]),
      [['generatedProteinComparison', 'orthogroup']]);
    const result = await load(bytes);
    assert.equal(result.status, 'ok', result.error?.stack);
    assert.equal(state.importedComparisonIntent.disposition, 'EDITABLE');
    assert.deepEqual(state.linearComparisonPlan, { ...DEFAULT_PLAN, mode: 'adjacent' });
    assert.equal(state.losatProgram.value, 'blastp');
    const filesData = await serializeActiveRenderFiles('linear', state);
    const comparisonPlanSnapshot = resolveLinearComparisonPlan({
      plan: state.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
      losatProgram: state.losatProgram.value, blastpMode: state.losat.blastp.mode
    });
    const candidate = buildCanonicalRenderRequest({ state, filesData, comparisonPlanSnapshot });
    assert.deepEqual(candidate.renderRequest.comparisons.map(item => [item.kind, item.mode]),
      [['generatedProteinComparison', 'orthogroup']]);
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
});

// B15: a 0.13.0 CLI protein sidecar (version 30, no renderRequest; provenance in
// tests/fixtures/sessions/cli-linear-protein.v30.provenance.json) is a CLI-only
// draft. The legacy migrator states the adjacent LOSATP comparison its CLI drew,
// so Generate rebuilds it; without --protein_blastp_mode the same sidecar has
// no Web comparison draft (No comparison).
await test('a 0.13.0 CLI Linear protein sidecar keeps the adjacent LOSATP plan it drew', async () => {
  const bytes = gunzipSync(await readFile(path.join(
    root, 'tests/fixtures/sessions/cli-linear-protein.v30.gbdraw-session.json.gz'
  )));
  const session = JSON.parse(bytes);
  assert.equal(session.version, 30);
  assert.equal(Object.hasOwn(session, 'renderRequest'), false);
  const loadDraft = async (document) => {
    const result = await load(JSON.stringify(document));
    assert.equal(result.status, 'ok', result.error?.stack);
    const filesData = await serializeActiveRenderFiles('linear', state);
    const comparisonPlanSnapshot = resolveLinearComparisonPlan({
      plan: state.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
      losatProgram: state.losatProgram.value, blastpMode: state.losat.blastp.mode
    });
    const candidate = buildCanonicalRenderRequest({ state, filesData, comparisonPlanSnapshot });
    return {
      plan: structuredClone(state.linearComparisonPlan),
      program: state.losatProgram.value,
      comparisons: candidate.renderRequest.comparisons.map(item => [item.kind, item.mode])
    };
  };
  assert.deepEqual(await loadDraft(session), {
    plan: { ...DEFAULT_PLAN, mode: 'adjacent' }, program: 'blastp',
    comparisons: [['generatedProteinComparison', 'orthogroup']]
  });
  // The 0.13.0 writer without a protein mode stores LOSAT program blastn.
  const withoutProtein = structuredClone(session);
  withoutProtein.config.cliOptions.rawArgs = ['--gbk', 'P1.gb', 'P2.gb', '-o', 'legacy-protein'];
  withoutProtein.cliInvocation.args = withoutProtein.config.cliOptions.rawArgs;
  withoutProtein.config.losatProgram = 'blastn';
  withoutProtein.config.adv.losatProgram = 'blastn';
  const none = await loadDraft(withoutProtein);
  assert.deepEqual([none.plan, none.comparisons], [DEFAULT_PLAN, []]);
});

// OV-38: the Web writers of Sessions 27–33 (release tags 0.12.0–0.13.0 and
// main before Session 39) saved every schema-4 Circular track slot row with
// `spacing: null`. With Custom Track Slots off the null is lossless, so the
// 0.13.0 Gallery Session loads and the field is not kept; with the slots on, or
// with a non-null value, Load fails and names the field and the row.
await test('a 0.13.0 Gallery Session with Custom Track Slots off loads without the null slot spacing', async () => {
  const bytes = gunzipSync(await readFile(path.join(
    root, 'tests/fixtures/sessions/BGC0000708-BGC0000713.v30.gbdraw-session.json.gz'
  )));
  const session = JSON.parse(bytes);
  assert.equal(session.version, 30);
  assert.equal(session.config.adv.circular_track_slots_enabled, false);
  assert.ok(session.config.adv.circular_track_slots.every(slot => slot.spacing === null));
  const result = await load(bytes);
  assert.equal(result.status, 'ok', JSON.stringify(result.error));
  assert.equal(state.mode.value, 'linear');
  assert.equal(state.adv.circular_track_slots.length, session.config.adv.circular_track_slots.length);
  assert.ok(state.adv.circular_track_slots.every(slot => !Object.hasOwn(slot, 'spacing')));
  for (const [label, edit, row] of [
    ['Custom Track Slots on', adv => { adv.circular_track_slots_enabled = true; }, 1],
    ['a non-null spacing', adv => { adv.circular_track_slots[1].spacing = '4px'; }, 2]
  ]) {
    const document = structuredClone(session);
    edit(document.config.adv);
    const failed = await load(JSON.stringify(document));
    assert.equal(failed.status, 'error', label);
    assert.equal(failed.error.code, 'TRACK_INVALID', label);
    assert.deepEqual(failed.error.context, { field: 'spacing', reason: 'OBSOLETE_TRACK_FIELD', slotIndex: row - 1 }, label);
    assert.match(failed.error.summary, new RegExp(`^The track settings are invalid\\. Track row ${row}\\. Field: spacing\\. `), label);
  }
});
