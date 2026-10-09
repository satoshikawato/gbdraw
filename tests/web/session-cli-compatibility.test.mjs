import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
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
// Save Session downloads its file through a link.
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({ addEventListener: () => {}, click: () => {}, remove: () => {}, parentNode: null })
};
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
  SESSION_VERSION, exportSession, importSession, getCommittedCanonicalRenderRequest, getCommittedCanonicalSession,
  serializeActiveRenderFiles, setUnmanagedConfigOverrideValidator
} = await import('../../gbdraw/web/js/services/config.js');
const {
  CANONICAL_REQUEST_SCHEMA, buildCanonicalRenderRequest, committedFeatureVisibilityMatches, projectCommittedEditorIntent
} = await import('../../gbdraw/web/js/services/session-request.js');
const { inheritCommittedComparisonIntent } = await import('../../gbdraw/web/js/services/imported-comparison-intent.js');
const { resolveLinearComparisonPlan } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { state } = await import('../../gbdraw/web/js/state.js');
// The composition root's transform of an older Session's Results (R13 port).
const { transformLegacyResultSvg } = await import('../../gbdraw/web/js/app/app-setup.js');
const { getSessionResourceSource, readFileBytes } = await import('../../gbdraw/web/js/services/file-content-cache.js');
const { getResourcePayloadOwner } = await import('../../gbdraw/web/js/services/resource-payload-owner.js');
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
} }, { transformLegacyResultSvg });
const DEFAULT_PLAN = { mode: 'none', defaultSource: 'losat', edges: [] };
// B15: a read-only CLI comparison (-b) is not a Web draft. The replacement draft
// is the Web default (No comparison), so Replace with current controls
// (app-setup.js: a valid plan with comparison intent) waits until the user sets
// up a comparison, and no LOSAT run starts from the draft.
const assertNoReplacementDraft = () => {
  assert.deepEqual(state.activeDrawing().linearComparisonPlan, DEFAULT_PLAN);
  const draft = resolveLinearComparisonPlan({
    plan: state.activeDrawing().linearComparisonPlan, sequences: state.linearSeqs, layout: [],
    losatProgram: state.activeDrawing().losatProgram.value, blastpMode: state.activeDrawing().losat.blastp.mode
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
        assert.equal(state.activeDrawing().importedComparisonIntent.disposition, 'EDITABLE');
        assert.deepEqual(getCommittedCanonicalRenderRequest(), session.renderRequest);
        assert.ok(state.results.value.length > 0);
        if (label === 'single') assert.equal(state.activeDrawing().form.labels_mode, 'out');
        assert.equal(session.renderRequest.diagramOptions.output.legend, legend);
        assert.equal(state.activeDrawing().form.legend, legend);
        if (mode === 'linear') {
          assert.deepEqual(state.linearSeqs.map(seq => seq.uid), fileRecordKeys(session.renderRequest));
          // B13: without -b the CLI commits only a disabled protein pipeline
          // (mode `none`, no pairs). The Web draft has no comparison, so the
          // request Generate builds has no comparison and starts no LOSAT run.
          assert.ok(session.renderRequest.comparisons.every(item => (
            item.kind === 'generatedProteinComparison' && item.mode === 'none' && item.pairs.length === 0
          )));
          assert.deepEqual(state.activeDrawing().linearComparisonPlan, { mode: 'none', defaultSource: 'losat', edges: [] });
          assert.notEqual(state.activeDrawing().losatProgram.value, 'blastp');
          const filesData = await serializeActiveRenderFiles('linear', state, state.activeDrawing());
          const comparisonPlanSnapshot = resolveLinearComparisonPlan({
            plan: state.activeDrawing().linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
            losatProgram: state.activeDrawing().losatProgram.value, blastpMode: state.activeDrawing().losat.blastp.mode
          });
          assert.equal(comparisonPlanSnapshot.hasLosatIntent, false);
          const candidate = buildCanonicalRenderRequest({ state, drawing: state.activeDrawing(), filesData, comparisonPlanSnapshot });
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
        // A draft outside `modes` is a Session-format error, even beside a
        // valid request (Session 46 keeps the Web draft in `modes`).
        const before = getCommittedCanonicalRenderRequest();
        const errorLog = console.error;
        console.error = () => {};
        try {
          const rejected = await load(JSON.stringify({ ...session, config: { adv: {} } }));
          assert.equal(rejected.status, 'error');
          assert.equal(rejected.error.code, 'INPUT_INVALID');
          assert.deepEqual(rejected.error.context, {field:'schema',reason:'FIELDS'});
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
  assert.equal(state.activeDrawing().importedComparisonIntent.disposition, 'PRESERVED_READ_ONLY');
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
    assert.equal(state.activeDrawing().importedComparisonIntent.disposition, 'PRESERVED_READ_ONLY');
    assertNoReplacementDraft();
    // The candidate Generate builds after Inherit (run-analysis.js: empty plan, committed comparison).
    const filesData = await serializeActiveRenderFiles('linear', state, state.activeDrawing());
    const candidate = buildCanonicalRenderRequest({
      state,
      drawing: state.activeDrawing(),
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
    assert.equal(state.activeDrawing().importedComparisonIntent.disposition, 'EDITABLE');
    assert.deepEqual(state.activeDrawing().linearComparisonPlan, { ...DEFAULT_PLAN, mode: 'adjacent' });
    assert.equal(state.activeDrawing().losatProgram.value, 'blastp');
    const filesData = await serializeActiveRenderFiles('linear', state, state.activeDrawing());
    const comparisonPlanSnapshot = resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
      losatProgram: state.activeDrawing().losatProgram.value, blastpMode: state.activeDrawing().losat.blastp.mode
    });
    const candidate = buildCanonicalRenderRequest({ state, drawing: state.activeDrawing(), filesData, comparisonPlanSnapshot });
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
    const filesData = await serializeActiveRenderFiles('linear', state, state.activeDrawing());
    const comparisonPlanSnapshot = resolveLinearComparisonPlan({
      plan: state.activeDrawing().linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
      losatProgram: state.activeDrawing().losatProgram.value, blastpMode: state.activeDrawing().losat.blastp.mode
    });
    const candidate = buildCanonicalRenderRequest({ state, drawing: state.activeDrawing(), filesData, comparisonPlanSnapshot });
    return {
      plan: structuredClone(state.activeDrawing().linearComparisonPlan),
      program: state.activeDrawing().losatProgram.value,
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
  assert.equal(state.activeDrawing().adv.circular_track_slots.length, session.config.adv.circular_track_slots.length);
  assert.ok(state.activeDrawing().adv.circular_track_slots.every(slot => !Object.hasOwn(slot, 'spacing')));
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

// D-04: the main b05a6bb8 Web writer (Session 33, request schema 2; provenance in
// tests/fixtures/sessions/scale-interval.provenance.json) saved a Scale Interval of 0
// as the flat override `scale_interval` and drew the automatic interval. Load shows
// the stored 0, and the next request sends it literally under the current path (R7);
// the Python request decoder reads it as automatic (tests/test_scale_interval_domain.py).
await test('a main v33 Web Session with Scale Interval 0 loads into the current request path', async () => {
  const bytes = gunzipSync(await readFile(path.join(
    root, 'tests/fixtures/sessions/scale-interval-zero-circular-web.v33.gbdraw-session.json.gz'
  )));
  const session = JSON.parse(bytes);
  assert.equal(session.version, 33);
  assert.equal(session.renderRequest.diagramOptions.configOverrides.scale_interval, 0);
  const result = await load(bytes);
  assert.equal(result.status, 'ok', result.error?.stack);
  assert.equal(state.mode.value, 'circular');
  assert.equal(state.activeDrawing().adv.scale_interval, 0);
  const filesData = await serializeActiveRenderFiles('circular', state, state.activeDrawing());
  const candidate = buildCanonicalRenderRequest({
    state, drawing: state.activeDrawing(), filesData, comparisonPlanSnapshot: null
  });
  const overrides = candidate.renderRequest.diagramOptions.configOverrides ?? {};
  assert.equal(overrides['objects.scale.interval'], 0);
  assert.equal(Object.hasOwn(overrides, 'scale_interval'), false);
});

// D-04 (review 2 item 1): Web Load checks the stored configuration and keeps the
// stored request as the committed request, which Worker helpers decode directly.
// Sessions that store a Scale Interval of 0 or less (main fe6861f0 CLI and typed API,
// main b05a6bb8 Web; scale-interval.provenance.json) load, and Load Feature Edits TSV
// reads against their committed request.
await test('a Session with a Scale Interval of 0 or less loads and reads a Feature Edits TSV', async () => {
  for (const [name, stored] of [
    ['scale-interval-zero-circular-cli.v44', (request) => request.diagramOptions.config.objects.scale.interval],
    ['scale-interval-negative-linear-api.v44', (request) => request.diagramOptions.configOverrides['objects.scale.interval']],
    // Load rebuilds a schema 2 request from the loaded fields (current path, literal 0).
    ['scale-interval-zero-circular-web.v33', (request) => request.diagramOptions.configOverrides['objects.scale.interval']]
  ]) {
    const fixture = path.join(root, `tests/fixtures/sessions/${name}.gbdraw-session.json.gz`);
    const result = await load(gunzipSync(await readFile(fixture)));
    assert.equal(result.status, 'ok', `${name}: ${JSON.stringify(result.error?.context ?? result.error)}`);
    const committed = getCommittedCanonicalRenderRequest();
    assert.ok(stored(committed) <= 0, name);
    const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-scale-interval-'));
    try {
      await writeFile(path.join(directory, 'request.json'), JSON.stringify(committed));
      await writeFile(path.join(directory, 'table.tsv'),
        'record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text\n');
      const output = execFileSync('python', ['-c', `
import gzip, json, sys
from pathlib import Path
from gbdraw.api import load_session_document, materialize_session
from gbdraw.web_support.feature_override_table import read_feature_override_table_json
fixture, directory = sys.argv[1], Path(sys.argv[2])
document = load_session_document(json.loads(gzip.decompress(Path(fixture).read_bytes())))
with materialize_session(document, output_directory=directory / 'out') as materialized:
    print(read_feature_override_table_json(
        str(directory / 'table.tsv'), (directory / 'request.json').read_text(),
        json.dumps({key: str(value) for key, value in materialized.resource_paths.items()}),
        str(directory / 'out')))
`, fixture, directory], { encoding: 'utf8', cwd: root });
      assert.deepEqual(JSON.parse(output), { rows: [], unmatchedRows: [] }, name);
    } finally {
      await rm(directory, { recursive: true, force: true });
    }
  }
});

// D-01, OV-221: a Session the CLI writes draws the CLI's figure when the Web app
// loads it. Load gives the drawing every table the request holds, so the next
// Generate, a label reflow, and a Web re-save keep them; the Sessions a CLI
// replay and Python save again (no Web draft either) load the same way, and so
// do a case's Sessions 40 and 41 from older CLI writers (`legacySessions`,
// OV-269). The CLI and Python cells of the matrix are
// tests/test_cli_session_cross_surface.py.
const crossSurface = (...args) => JSON.parse(execFileSync('python', [
  path.join(root, 'tests/web/helpers/cli-session-cross-surface.py'), ...args
], { cwd: root, encoding: 'utf8', timeout: 1_800_000, maxBuffer: 16 * 1024 * 1024 }));
// The drawing state each table option of the CLI fills at Load.
// A loaded table holds exactly the rows of its fixture file, which has no header.
const fixtureRows = file => readFileSync(path.join(root, file), 'utf8').split(/\r?\n/).filter(line => line.trim()).length;
const LOADED_TABLES = {
  '-t': (drawing, file) => drawing.manualSpecificRules.length === fixtureRows(file),
  '-d': drawing => drawing.currentColors.value.CDS === '#123456',
  '--feature_visibility_table': (drawing, file) => drawing.featureVisibilityRules.value.length === fixtureRows(file),
  '--label_table': (drawing, file) => drawing.canonicalLabelOverrideRows.value.length === fixtureRows(file),
  '--label_whitelist': (drawing, file) => drawing.filterMode.value === 'Whitelist'
    && drawing.manualWhitelist.length === fixtureRows(file),
  '--label_blacklist': drawing => drawing.filterMode.value === 'Blacklist' && drawing.manualBlacklist.value !== '',
  '--qualifier_priority': (drawing, file) => drawing.manualPriorityRules.length === fixtureRows(file),
  '--feature_override_table': drawing => Object.keys(drawing.featureOverrides).length > 0
};
// A Session of the request and resources the Web renders.
const renderedSession = async ({ renderRequest, resources }) => JSON.stringify({
  format: 'gbdraw-session', version: SESSION_VERSION, createdAt: '2026-10-09T00:00:00.000Z', renderRequest,
  resources: Object.fromEntries(await Promise.all(Object.entries(resources).map(async ([id, descriptor]) => {
    if (typeof descriptor.data === 'string') return [id, descriptor];
    const bytes = await readFileBytes(getResourcePayloadOwner(descriptor));
    return [id, { kind: descriptor.kind || 'web-file', name: descriptor.name, type: descriptor.type || 'application/octet-stream',
      size: bytes.byteLength, lastModified: 0, encoding: 'base64', data: Buffer.from(bytes).toString('base64') }];
  }))),
  results: [], editorState: { featureCatalog: null }
});
// The request the next Generate builds (run-analysis.js), inheriting a
// read-only CLI comparison (-b) as Inherit saved comparison does.
const nextGenerate = async () => {
  const drawing = state.activeDrawing();
  const filesData = await serializeActiveRenderFiles(state.mode.value, state, drawing);
  const inherit = drawing.importedComparisonIntent.disposition === 'PRESERVED_READ_ONLY';
  const candidate = buildCanonicalRenderRequest({ state, drawing,
    filesData: inherit ? { ...filesData, linearCanonicalComparisons: [] } : filesData,
    comparisonPlanSnapshot: state.mode.value === 'linear' ? resolveLinearComparisonPlan({
      plan: drawing.linearComparisonPlan, sequences: filesData.linearSeqs, layout: [],
      losatProgram: drawing.losatProgram.value, blastpMode: drawing.losat.blastp.mode
    }) : null });
  if (inherit) inheritCommittedComparisonIntent({ candidate, committed: getCommittedCanonicalSession() });
  return candidate;
};

// Generate reuses committed LOSATP results only under the visibility rules
// they were drawn with; the CLI stores its rules as `featureVisibilityTable`.
await test('a CLI Session keeps its comparisons under its own visibility rules', () => {
  const rules = '*\tCDS\tproduct\tmajor capsid\toff\n';
  const committed = {
    renderRequest: { diagramOptions: { featureVisibilityTable: { resourceId: 'rules', representation: 'canonicalTsv' } } },
    resources: { rules: { kind: 'canonical-tsv', encoding: 'base64',
      data: Buffer.from(`record_id\tfeature_type\tqualifier\tvalue\taction\n${rules}`).toString('base64') } }
  };
  assert.equal(committedFeatureVisibilityMatches(committed, rules), true);
  assert.equal(committedFeatureVisibilityMatches(committed, ''), false);
});

// OV-273: Load parses and normalizes only a Result without composition
// metadata (the Session 40 writers of main 8228ffab and 7aad9e3e); a current
// writer's Results take the EMPTY admission plan and no application parse,
// and `currentLegacyNormalizationCount` counts the normalized Results.
await test('only a Result without composition metadata takes the legacy composition at Load', async () => {
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-cli-composition-'));
  const admission = async (bytes) => {
    const events = [], metrics = [];
    globalThis.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: event => events.push(event), onStructuralMetric: metric => metrics.push(metric)
    };
    try {
      assert.equal((await load(bytes)).status, 'ok');
    } finally {
      delete globalThis.__GBDRAW_TEST_HOOKS__;
    }
    const phase = events.find(event => event.name === 'svg.admission-started' && event.mutationKind).phase;
    const total = name => metrics.filter(metric => metric.name === name && metric.phase === phase)
      .reduce((sum, metric) => sum + metric.value, 0);
    return {
      plan: events.filter(event => event.name === 'svg.admission-started' && event.phase === phase).map(event => event.mutationKind),
      parses: total('applicationSvgParseCount'), normalized: total('currentLegacyNormalizationCount')
    };
  };
  try {
    const current = path.join(directory, 'current.gbdraw-session.json');
    execFileSync('python', ['-m', 'gbdraw.cli', 'linear', '--gbk', lambda, '-o', path.join(directory, 'current'),
      '-f', 'svg', '--session_output', current], { cwd: directory, env: { ...process.env, PYTHONPATH: root }, stdio: 'pipe' });
    assert.deepEqual(await admission(await readFile(current)), { plan: ['EMPTY'], parses: 0, normalized: 0 });
    const early = gunzipSync(await readFile('tests/fixtures/sessions/cli-linear-tables.v40-early.gbdraw-session.json.gz'));
    assert.equal(early.toString().includes('data-gbdraw-composition-schema='), false);
    // The fake DOM serializes the parsed text unchanged; the Playwright
    // OV-273 case checks the composition the transform writes.
    assert.deepEqual(await admission(early), { plan: ['MUTATING'], parses: 1, normalized: 1 });
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
});

await test('every load and re-save of a CLI Session in Web draws the CLI figure', async () => {
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-cli-cross-'));
  try {
    const checks = [];
    for (const entry of crossSurface('prepare', directory)) {
      const first = checks.length;
      const sources = { 'CLI Session': entry.cli_session, 'CLI re-save': entry.cli_resave, 'Python re-save': entry.python_resave };
      for (const [source, file] of Object.entries(sources)) {
        const result = await load(await readFile(file));
        assert.equal(result.status, 'ok', `${entry.id}, ${source}: ${result.error?.stack}`);
        const drawing = state.activeDrawing();
        for (const [option, loaded] of Object.entries(LOADED_TABLES)) {
          const index = entry.args.indexOf(option);
          if (index >= 0) assert.ok(loaded(drawing, entry.args[index + 1]), `${entry.id}, ${source}: ${option}`);
        }
        const generated = path.join(directory, entry.id, `${source} Generate.json`);
        await writeFile(generated, await renderedSession(await nextGenerate()));
        checks.push({ label: `${entry.id}: ${source}, Generate`, session: generated, via: 'python' });
        if (source !== 'CLI Session') continue;
        // A label reflow renders the committed request with the drawing's tables.
        const reflow = path.join(directory, entry.id, 'reflow.json');
        await writeFile(reflow, await renderedSession(projectCommittedEditorIntent({
          committed: getCommittedCanonicalSession(), state, drawing })));
        checks.push({ label: `${entry.id}: label reflow`, session: reflow, via: 'python' });
        // Save before Generate keeps the CLI request and the tables in the draft.
        const saved = await exportSession(entry.id);
        assert.equal(saved.status, 'saved', saved.error?.message);
        const webSave = path.join(directory, entry.id, 'web-save.json');
        await writeFile(webSave, gunzipSync(Buffer.from(await saved.blob.arrayBuffer())));
        checks.push({ label: `${entry.id}: Web save, CLI`, session: webSave, via: 'cli' },
          { label: `${entry.id}: Web save, Python`, session: webSave, via: 'python' });
        assert.equal((await load(await readFile(webSave))).status, 'ok');
        const resaved = path.join(directory, entry.id, 'web-save Generate.json');
        await writeFile(resaved, await renderedSession(await nextGenerate()));
        checks.push({ label: `${entry.id}: Web save, Web Generate`, session: resaved, via: 'python' });
      }
      checks.slice(first).forEach(item => Object.assign(item, { mode: entry.mode, expected: entry.svg }));
    }
    const checksFile = path.join(directory, 'checks.json');
    await writeFile(checksFile, JSON.stringify(checks));
    const failed = crossSurface('check', checksFile).filter(item => !item.equal);
    assert.deepEqual(failed, []);
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
});
