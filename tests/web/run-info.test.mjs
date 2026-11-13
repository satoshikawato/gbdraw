import assert from 'node:assert/strict';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import { File } from 'node:buffer';
import { spawnSync } from 'node:child_process';
import { cp, readFile, stat, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { delimiter, join } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';
import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { evaluatePython } from './helpers/python-evaluator.mjs';
import { readSavedSession } from './helpers/session-documents.mjs';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-run-info-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js'), join(tempDir, 'js'), { recursive: true });
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}', 'utf8');
const modulePath = join(tempDir, 'js', 'app', 'run-info.js');

const {
  buildRunInfo,
  LOSAT_DATABASE_SCOPE_NOTE,
  LOSAT_EVALUE_SCALE_NOTE,
  LOSAT_RING_EVALUE_SCALE_NOTE,
  buildSourceRecipe,
  quoteShellArg,
  reproducibilityLabel,
  shownRunInfoNoteParts,
  summarizeLosatRuntimes
} = await import(pathToFileURL(modulePath));
const { adoptCurrentSessionResources, createSessionResourceFileView, fileBackedDescriptor } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'session-resource-backing.js'))
);
const { holdReadBytes, readFileContentFacts, readFileText } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'file-content-cache.js'))
);
const { setResourcePayloadOwner } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'resource-payload-owner.js'))
);
const { readCanonicalResourceFileBytes, readCanonicalResourceRecordCount } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'session-request.js'))
);
const { createDiagramResourceTransport } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'diagram-resource-staging.js'))
);

const base64 = (text) => Buffer.from(text).toString('base64');
const resource = (kind, name, text = 'fixture') => ({
  kind,
  name,
  type: 'text/plain',
  size: Buffer.byteLength(text),
  lastModified: 0,
  encoding: 'base64',
  data: base64(text)
});
const presentation = () => ({
  label: null,
  subtitle: null,
  reverseComplement: false,
  gridRow: null,
  gridColumn: null
});
const baseOptions = (mode) => ({
  configOverrides: {
    'canvas.show_gc': true,
    'canvas.show_skew': true,
    'canvas.show_depth': false,
    'canvas.strandedness': false,
    'canvas.resolve_overlaps': false,
    'objects.scale.show': true,
    [mode === 'circular' ? 'labels.circular.scope' : 'labels.linear.scope']: 'none'
  },
  tracks: {
    circularTrackSlots: null,
    circularTrackAxisIndex: null,
    linearTrackSlots: null,
    linearTrackAxisIndex: null,
    centerReservedRadius: null
  },
  output: { legend: 'right', plotTitlePosition: mode === 'circular' ? 'none' : 'bottom' },
  colors: {
    colorTable: null,
    colorTableFile: null,
    defaultColors: null,
    defaultColorsPalette: 'default',
    defaultColorsFile: null
  },
  selectedFeaturesSet: ['CDS', 'rRNA'],
  featureShapes: { repeat_region: 'underlay' },
  dinucleotide: 'GC',
  window: null,
  step: null,
  depthWindow: null,
  depthStep: null,
  plotTitle: null,
  plotTitleFontSize: null
});
// A committed request is of the current schema: Load reads a stored one to it.
const canonical = ({ mode, records, resources, comparisons = [], webFiles = {} }) => ({
  renderRequest: {
    schema: CANONICAL_REQUEST_SCHEMA,
    mode,
    grouping: 'single',
    records: records.map((record) => ({ display: { isCircular: null, startCoordinate: null }, ...record })),
    diagramOptions: { featurePlacements: [], featureOverrides: [], ...baseOptions(mode) },
    layout: {},
    comparisons,
    output: {
      prefix: mode === 'circular' ? 'input' : 'comparison',
      formats: ['svg'],
      overwrite: false,
      interactiveMetadataPolicy: 'auto'
    }
  },
  resources,
  webFiles
});
const assertCliParserAccepts = async (recipe) => {
  const visibleArgs = recipe.args.map((token) => recipe.fileMetadata.get(token)?.name || token);
  const parsed = await evaluatePython([
    'import contextlib, io',
    'from gbdraw.circular import _get_args as circular_args',
    'from gbdraw.linear import _get_args as linear_args',
    'stderr = io.StringIO()',
    'try:',
    '    with contextlib.redirect_stderr(stderr):',
    '        (linear_args if payload["mode"] == "linear" else circular_args)(payload["args"])',
    '    status = 0',
    'except SystemExit as error:',
    '    status = error.code',
    'result = {"status": status, "stderr": stderr.getvalue()}'
  ].join('\n'), { mode: recipe.mode, args: visibleArgs });
  assert.equal(parsed.status, 0, parsed.stderr);
};

const materializeAndRunRecipe = async (session, recipe, info, fixture) => {
  const runDir = await mkdtemp(join(tmpdir(), 'gbdraw-source-recipe-'));
  const generatedFiles = Array.isArray(info.sourceRecipe.generatedFiles)
    ? info.sourceRecipe.generatedFiles
    : recipe.generatedFiles;
  const namesBySlot = new Map(
    info.sourceRecipe.invocation.fileBindings.map(({ slot, name }) => [slot, name])
  );
  for (const [path, metadata] of recipe.fileMetadata) {
    if (metadata.kind === 'generated') {
      const generated = generatedFiles.find((file) => file.path === path);
      assert.ok(generated, `${fixture}: missing generated file ${path}`);
      await writeFile(join(runDir, generated.name), generated.data);
      continue;
    }
    const name = namesBySlot.get(metadata.slot) || metadata.name;
    const descriptor = session.resources[metadata.resourceId];
    assert.ok(descriptor, `${fixture}: missing resource ${metadata.resourceId}`);
    const data = descriptor.encoding === 'base64'
      ? Buffer.from(descriptor.data, 'base64')
      : descriptor.data;
    await writeFile(join(runDir, name), data);
  }
  const args = info.sourceRecipe.commandArgs.slice(2);
  const executed = spawnSync(
    process.env.PYTHON || 'python',
    ['-m', 'gbdraw.cli', recipe.mode, ...args],
    {
      cwd: runDir,
      encoding: 'utf8',
      timeout: 120_000,
      env: {
        ...process.env,
        PYTHONPATH: [repoRoot, process.env.PYTHONPATH].filter(Boolean).join(delimiter)
      }
    }
  );
  assert.equal(executed.status, 0, `${fixture}: ${executed.stderr || executed.stdout}`);
  const outputIndex = args.indexOf('-o');
  assert.notEqual(outputIndex, -1, `${fixture}: emitted command lacks -o`);
  const outputPath = join(runDir, `${args[outputIndex + 1]}.svg`);
  const output = await stat(outputPath);
  assert.ok(output.size > 0, `${fixture}: CLI emitted an empty SVG`);
  return outputPath;
};

const sourceFeatureIds = (genbank) => evaluatePython([
  'import io',
  'from Bio import SeqIO',
  'from gbdraw.features.source import capture_source_features',
  'record = SeqIO.read(io.StringIO(payload), "genbank")',
  'result = [entry.biological_feature_id for entry in capture_source_features(record).entries]'
].join('\n'), genbank);

assert.equal(quoteShellArg('simple.tsv'), 'simple.tsv');
assert.equal(quoteShellArg('two words.tsv'), "'two words.tsv'");
assert.equal(quoteShellArg("Bob's genome.gbk"), "'Bob'\\''s genome.gbk'");
assert.equal(quoteShellArg(''), "''");

let circularGallerySession = null;
let linearGallerySession = null;
for (const { fixture, mode } of [
  { fixture: 'HmmtDNA_ATskew.gbdraw', mode: 'circular' },
  { fixture: 'lambda_basic_linear.gbdraw', mode: 'linear' }
]) {
  const saved = await readSavedSession(new Blob([await readFile(
    join(repoRoot, 'gbdraw', 'web', 'gallery', 'sessions', fixture)
  )]));
  const session = savedDrawingView(saved, saved.ui.activeDrawingId);
  if (mode === 'circular') circularGallerySession = session;
  else linearGallerySession = session;
  const sourceRecipe = await buildSourceRecipe(session);
  assert.equal(sourceRecipe.available, true, `${fixture}: ${sourceRecipe.unavailableReason}`);
  assert.equal(sourceRecipe.mode, mode);
  assert.ok(sourceRecipe.args.includes('--gbk'), `${fixture}: expected direct GenBank input`);
  assert.ok(sourceRecipe.args.includes('-o'), `${fixture}: expected output prefix`);
  assert.ok(sourceRecipe.args.includes('-f'), `${fixture}: expected output formats`);
  await assertCliParserAccepts(sourceRecipe);
  const info = buildRunInfo({ mode: sourceRecipe.mode, sourceRecipe });
  await materializeAndRunRecipe(session, sourceRecipe, info, fixture);
}

// A FASTA given without GFF3 is `--fasta` alone in both modes (PD-OI-068 rev 2); a records table names it in its
// `fasta` column. The recipe runs: the CLI draws the same input.
test('source recipe of a FASTA alone gives --fasta in both modes and a fasta records-table column', async () => {
  const text = `>alpha\n${'ACGT'.repeat(30)}\n`;
  for (const mode of ['circular', 'linear']) {
    for (const labelled of [false, true]) {
      const session = canonical({
        mode,
        records: [{
          recordKey: 'alpha',
          cardinality: 'all',
          source: { kind: 'fasta', fastaResourceId: 'source' },
          selector: null,
          windows: [],
          // A Circular record label needs a records table.
          presentation: { ...presentation(), label: labelled ? 'Alpha' : null }
        }],
        resources: { source: resource('fasta', 'alpha.fasta', text) },
        webFiles: { resourceOriginalNames: { source: 'alpha.fasta' } }
      });
      const recipe = await buildSourceRecipe(session);
      assert.equal(recipe.available, true, `${mode}: ${recipe.unavailableReason}`);
      const usesTable = mode === 'circular' && labelled;
      assert.equal(recipe.args.includes('--records_table'), usesTable, mode);
      assert.equal(recipe.args.includes('--gbk') || recipe.args.includes('--gff'), false, mode);
      if (usesTable) {
        const [header, row] = recipe.generatedFiles[0].data.split('\n');
        assert.deepEqual(header.split('\t').slice(0, 2), ['fasta', 'record_label']);
        assert.deepEqual(row.split('\t').slice(0, 2), ['alpha.fasta', 'Alpha']);
      } else {
        const fasta = recipe.args.indexOf('--fasta');
        assert.notEqual(fasta, -1, mode);
        assert.equal(recipe.fileMetadata.get(recipe.args[fasta + 1])?.name, 'alpha.fasta');
      }
      await assertCliParserAccepts(recipe);
      const info = buildRunInfo({ mode: recipe.mode, sourceRecipe: recipe });
      await materializeAndRunRecipe(session, recipe, info, `fasta ${mode}${labelled ? ' labelled' : ''}`);
    }
  }
});

test('source recipe reuses materialized source text and preserves record selection', async () => {
  for (const format of ['genbank', 'fasta']) {
    for (const recordCount of [1, 2]) {
      const text = Array.from({ length: recordCount }, (_, index) => (
        format === 'genbank' ? `LOCUS       record${index}\n//\n` : `>record${index}\nACGT\n`
      )).join('');
      const session = canonical({
        mode: 'circular',
        records: [{
          recordKey: 'selected',
          cardinality: 'exactly_one',
          source: format === 'genbank'
            ? { kind: 'genbank', resourceId: 'source' }
            : { kind: 'gffFasta', gffResourceId: 'gff', fastaResourceId: 'source' },
          selector: { kind: 'recordIndex', index: 0 },
          windows: [],
          presentation: presentation()
        }],
        resources: {
          source: resource(format, `source.${format === 'genbank' ? 'gbk' : 'fasta'}`, text),
          gff: resource('gff', 'source.gff', '##gff-version 3\n'),
          unused: { ...resource('web-file', 'unused.txt'), data: '%%%' }
        },
        webFiles: { resourceOriginalNames: { source: `source.${format}`, gff: 'source.gff' } }
      });
      const expected = await buildSourceRecipe(session);
      assert.equal(expected.available, true, expected.unavailableReason);
      assert.equal(expected.args.includes('--records_table'), recordCount > 1);
      const readRecipe = () => buildSourceRecipe({
        ...session,
        readResourceRecordCount: (resourceId, kind) => (
          readCanonicalResourceRecordCount(session.resources, resourceId, kind)
        )
      });
      assert.deepEqual(await readRecipe(), expected);

      const table = adoptCurrentSessionResources(session.resources);
      const view = createSessionResourceFileView(table, 'source');
      setResourcePayloadOwner(session.resources.source, view);
      assert.equal(await readFileText(view), text);
      const nativeAtob = globalThis.atob;
      let decodes = 0;
      globalThis.atob = (...args) => {
        if (args[0] === session.resources.source.data) decodes += 1;
        return nativeAtob(...args);
      };
      try {
        assert.deepEqual(await readRecipe(), expected);
        assert.deepEqual(await readRecipe(), expected);
        assert.equal(decodes, 0, `${format}: recipe must reuse the materialized source`);
      } finally {
        globalThis.atob = nativeAtob;
      }
    }
  }
});

test('source recipe counts survive transfer and follow same-name source replacements', async () => {
  for (const input of ['native', 'lazy']) {
    for (const format of ['genbank', 'fasta']) {
      const transport = createDiagramResourceTransport();
      let previousToken;
      for (const recordCount of [1, 2, 1]) {
        const text = Array.from({ length: recordCount }, (_, index) => (
          format === 'genbank' ? `LOCUS       record${index}\n//\n` : `>record${index}\nACGT\n`
        )).join('').padEnd(128, ' ');
        const descriptor = resource(format, 'same-source.txt', text);
        const session = canonical({
          mode: 'circular',
          records: [{
            recordKey: 'selected',
            cardinality: 'exactly_one',
            source: format === 'genbank'
              ? { kind: 'genbank', resourceId: 'source' }
              : { kind: 'gffFasta', gffResourceId: 'gff', fastaResourceId: 'source' },
            selector: { kind: 'recordIndex', index: recordCount - 1 },
            windows: [],
            presentation: presentation()
          }],
          resources: {
            source: descriptor,
            gff: resource('gff', 'same-source.gff', '##gff-version 3\n'),
            unused: { ...resource('web-file', 'unused.txt'), data: '%%%' }
          },
          webFiles: { resourceOriginalNames: { source: 'same-source.txt', gff: 'same-source.gff' } }
        });
        const expected = await buildSourceRecipe(session);
        assert.equal(expected.available, true, expected.unavailableReason);
        assert.equal(expected.args.includes('--records_table'), recordCount === 2);
        if (recordCount === 2) {
          assert.match(expected.generatedFiles[0].data, /#2/);
        }
        let reads = 0;
        let decodes = 0;
        const file = input === 'native'
          ? new class extends File {
            async arrayBuffer() {
              reads += 1;
              return super.arrayBuffer();
            }
          }([text], descriptor.name, { lastModified: 7 })
          : createSessionResourceFileView(adoptCurrentSessionResources(session.resources), 'source');
        const readRecipe = () => buildSourceRecipe({
          ...session,
          readResourceRecordCount: (resourceId, kind) => (
            readCanonicalResourceRecordCount(session.resources, resourceId, kind)
          )
        });
        setResourcePayloadOwner(descriptor, file);
        // Count decodes through the resource owner, independent of the decoder API.
        globalThis.__GBDRAW_TEST_HOOKS__ = {
          onStructuralMetric: ({ name, resourceId }) => {
            if (name === 'base64DecodeCount' && resourceId === 'source') decodes += 1;
          }
        };
        try {
          assert.equal(await readFileText(file), text);
          if (input === 'native') {
            // The discovery read ends once its facts are hashed; a Generate that started
            // before that would share its bytes instead of reading the File again.
            await readFileContentFacts(file);
            await new Promise((resolve) => setImmediate(resolve));
          }
          // Generate holds its reads: the record count and the transfer share one read of the File.
          const releaseReads = holdReadBytes();
          assert.deepEqual(await readRecipe(), expected);
          const first = await transport.prepare({ request: session.renderRequest, resources: session.resources });
          const sourceTransfer = first.stagedResources.find(({ resourceId }) => resourceId === 'source');
          assert.ok(sourceTransfer);
          assert.notEqual(sourceTransfer.cacheToken, previousToken);
          previousToken = sourceTransfer.cacheToken;
          assert.equal(await sourceTransfer.blob.text(), text);
          first.commit();
          releaseReads();
          for (let run = 0; run < 2; run += 1) {
            // Canonical projection creates a fresh descriptor for the same payload owner.
            session.resources.source = setResourcePayloadOwner({ ...descriptor }, file);
            assert.deepEqual(await readRecipe(), expected);
            const warm = await transport.prepare({ request: session.renderRequest, resources: session.resources });
            assert.equal(warm.stagedResources.length, 0);
            assert.equal(warm.resourceManifest.find(({ resourceId }) => resourceId === 'source').cacheToken, previousToken);
            warm.commit();
          }
          // A File is read once to discover it and once more by the first Generate; a Worker cache hit reads nothing.
          assert.equal(reads, input === 'native' ? 2 : 0);
          assert.equal(decodes, input === 'lazy' ? 1 : 0);
          // Discovery may refill released content; recipe construction must not do so.
          assert.equal(await readFileText(file), text);
          assert.deepEqual(await readRecipe(), expected);
          assert.equal(reads, input === 'native' ? 3 : 0);
          assert.equal(decodes, input === 'lazy' ? 2 : 0);
        } finally {
          delete globalThis.__GBDRAW_TEST_HOOKS__;
        }
      }
    }
  }
  const locus = (id) => `LOCUS       ${id}                    10 bp    DNA     linear   UNK 01-JAN-2000\n//\n`;
  const resources = { source: resource('genbank', 'mutable.gb', locus('one')) };
  assert.equal(await readCanonicalResourceRecordCount(resources, 'source', 'genbank'), 1);
  Object.assign(resources.source, resource('genbank', 'mutable.gb', `${locus('one')}${locus('two')}`));
  assert.equal(await readCanonicalResourceRecordCount(resources, 'source', 'genbank'), 2);
  // Biopython starts a record only at "LOCUS" plus seven spaces, so this malformed line is not one.
  Object.assign(resources.source, resource('genbank', 'mutable.gb', `${locus('one')}LOCUS AB1 100 bp DNA\n//\n`));
  assert.equal(await readCanonicalResourceRecordCount(resources, 'source', 'genbank'), 1);
});

test('source recipe preserves selectedFeaturesSet empty, invalid, and non-empty semantics', async () => {
  for (const value of ['absent', null]) {
    const defaulted = structuredClone(circularGallerySession);
    if (value === 'absent') delete defaulted.renderRequest.diagramOptions.selectedFeaturesSet;
    else defaulted.renderRequest.diagramOptions.selectedFeaturesSet = value;
    const defaultedRecipe = await buildSourceRecipe(defaulted);
    assert.equal(defaultedRecipe.available, true, defaultedRecipe.unavailableReason);
    assert.equal(defaultedRecipe.args.includes('-k'), false);
  }

  const empty = structuredClone(circularGallerySession);
  empty.renderRequest.diagramOptions.selectedFeaturesSet = [];
  const emptyRecipe = await buildSourceRecipe(empty);
  assert.equal(emptyRecipe.available, true, emptyRecipe.unavailableReason);
  assert.equal(emptyRecipe.args.includes('-k'), false);
  await assertCliParserAccepts(emptyRecipe);
  const emptyInfo = buildRunInfo({ mode: 'circular', sourceRecipe: emptyRecipe });
  await materializeAndRunRecipe(empty, emptyRecipe, emptyInfo, 'empty-selected-features');

  const invalid = structuredClone(circularGallerySession);
  invalid.renderRequest.diagramOptions.selectedFeaturesSet = [''];
  const invalidRecipe = await buildSourceRecipe(invalid);
  assert.equal(invalidRecipe.available, false);
  assert.match(invalidRecipe.unavailableReason, /selected feature|non-empty/i);

  const nonEmpty = structuredClone(circularGallerySession);
  nonEmpty.renderRequest.diagramOptions.selectedFeaturesSet = ['CDS', 'tRNA'];
  const nonEmptyRecipe = await buildSourceRecipe(nonEmpty);
  assert.equal(nonEmptyRecipe.available, true, nonEmptyRecipe.unavailableReason);
  const featureOptionIndex = nonEmptyRecipe.args.indexOf('-k');
  assert.notEqual(featureOptionIndex, -1);
  assert.equal(nonEmptyRecipe.args[featureOptionIndex + 1], 'CDS,tRNA');
  await assertCliParserAccepts(nonEmptyRecipe);
});

test('source recipe preserves canonical string track slots', async () => {
  const stringSlot = structuredClone(circularGallerySession);
  stringSlot.renderRequest.diagramOptions.tracks.circularTrackSlots = ['ticks:ticks'];
  stringSlot.renderRequest.diagramOptions.tracks.circularTrackAxisIndex = 0;
  const stringRecipe = await buildSourceRecipe(stringSlot);
  assert.equal(stringRecipe.available, true, stringRecipe.unavailableReason);
  const slotOptionIndex = stringRecipe.args.indexOf('--circular_track_slot');
  assert.notEqual(slotOptionIndex, -1);
  assert.equal(stringRecipe.args[slotOptionIndex + 1], 'ticks:ticks');
  assert.equal(stringRecipe.args.includes('gc_skew:dinucleotide_skew'), false);
  await assertCliParserAccepts(stringRecipe);
  const stringInfo = buildRunInfo({ mode: 'circular', sourceRecipe: stringRecipe });
  await materializeAndRunRecipe(stringSlot, stringRecipe, stringInfo, 'ticks-string-slot');
});

test('source recipe serializes validated object slots without dropping child fields', async () => {
  const objectSlot = structuredClone(circularGallerySession);
  objectSlot.renderRequest.diagramOptions.tracks.circularTrackSlots = [{
    kind: 'circularTrackSlot',
    id: 'ticks',
    renderer: 'ticks',
    enabled: true,
    side: null,
    radius: null,
    width: null,
    z: 0,
    params: { tick_label_layout: 'label_out_tick_in' },
    innerGapPx: null,
    outerGapPx: null
  }];
  objectSlot.renderRequest.diagramOptions.tracks.circularTrackAxisIndex = 0;
  const objectRecipe = await buildSourceRecipe(objectSlot);
  assert.equal(objectRecipe.available, true, objectRecipe.unavailableReason);
  const optionIndex = objectRecipe.args.indexOf('--circular_track_slot');
  assert.match(
    objectRecipe.args[optionIndex + 1],
    /^ticks:ticks@tick_label_layout=label_out_tick_in$/
  );
  await assertCliParserAccepts(objectRecipe);

  const unsupportedChild = structuredClone(objectSlot);
  unsupportedChild.renderRequest.diagramOptions.tracks.circularTrackSlots[0].params.nt = 'GC';
  const unsupportedRecipe = await buildSourceRecipe(unsupportedChild);
  assert.equal(unsupportedRecipe.available, false);
  assert.match(unsupportedRecipe.unavailableReason, /params\.nt|losslessly/i);
});

// The slot editor's Extend to writes to_slot; the recipe carries it as the CLI reads it.
test('source recipe carries an annotation slot\'s to_slot in both modes', async () => {
  for (const mode of ['circular', 'linear']) {
    const session = structuredClone(mode === 'circular' ? circularGallerySession : linearGallerySession);
    const row = (id, renderer, side, z, params) => (mode === 'circular'
      ? { kind: 'circularTrackSlot', id, renderer, enabled: true, side, radius: null, width: null, z, params, innerGapPx: null, outerGapPx: null }
      : { kind: 'linearTrackSlot', id, renderer, enabled: true, side, height: null, spacing: null, z, params });
    const tracks = session.renderRequest.diagramOptions.tracks;
    tracks[`${mode}TrackSlots`] = [
      row('features', 'features', 'overlay', 0, mode === 'circular' ? { lane_direction: 'split' } : {}),
      row('gc', 'dinucleotide_content', mode === 'circular' ? 'inside' : 'below', 0, { nt: 'GC' }),
      row('regions', 'annotations', 'overlay', -1,
        { set_id: 'regions', anchor_slot: 'features', to_slot: 'gc', layer: 'underlay' })
    ];
    tracks[`${mode}TrackAxisIndex`] = 0;
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, `${mode}: ${recipe.unavailableReason}`);
    const tokens = recipe.args.filter((_, index) => recipe.args[index - 1] === `--${mode}_track_slot`);
    const regions = tokens.find((token) => token.startsWith('regions:'));
    const parsed = await evaluatePython([
      `from gbdraw.tracks.${mode} import parse_${mode}_track_slot`,
      `result = dict(parse_${mode}_track_slot(payload).params)`
    ].join('\n'), regions);
    assert.equal(parsed.to_slot, 'gc', `${mode}: ${regions}`);
    assert.equal(parsed.anchor_slot, 'features');
  }
});

test('source recipe distinguishes absent, null, and explicit empty track slots', async () => {
  const absent = structuredClone(circularGallerySession);
  delete absent.renderRequest.diagramOptions.tracks.circularTrackSlots;
  delete absent.renderRequest.diagramOptions.tracks.circularTrackAxisIndex;
  assert.equal((await buildSourceRecipe(absent)).available, true);

  const nullSlots = structuredClone(circularGallerySession);
  nullSlots.renderRequest.diagramOptions.tracks.circularTrackSlots = null;
  nullSlots.renderRequest.diagramOptions.tracks.circularTrackAxisIndex = null;
  assert.equal((await buildSourceRecipe(nullSlots)).available, true);

  const emptySlots = structuredClone(circularGallerySession);
  emptySlots.renderRequest.diagramOptions.tracks.circularTrackSlots = [];
  emptySlots.renderRequest.diagramOptions.tracks.circularTrackAxisIndex = null;
  const emptyRecipe = await buildSourceRecipe(emptySlots);
  const emptyInfo = buildRunInfo({
    mode: 'circular',
    sourceRecipe: emptyRecipe,
    exactReplayArgs: ['--session', '/exact-session'],
    fileMetadata: {
      '/exact-session': {
        name: 'exact-session.gbdraw-session.json',
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      }
    }
  });
  assert.equal(emptyRecipe.available, false);
  assert.match(emptyRecipe.unavailableReason, /empty.*track|track.*empty|no tracks/i);
  assert.equal(emptyInfo.sourceRecipe.available, false);
  assert.equal(emptyInfo.exactReplay.available, true);

  for (const slot of ['future:future', 'ticks:ticks@nt=GC']) {
    const invalid = structuredClone(circularGallerySession);
    invalid.renderRequest.diagramOptions.tracks.circularTrackSlots = [slot];
    assert.equal((await buildSourceRecipe(invalid)).available, false, slot);
  }
});

test('allocated helper names propagate into comparisons.tsv and the executable recipe', async () => {
  const genomeTemplate = linearGallerySession.resources['record-1-genbank'];
  const collision = canonical({
    mode: 'linear',
    records: ['alpha', 'middle', 'omega'].map((recordKey, index) => ({
      recordKey,
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: `genome-${index + 1}` },
      selector: null,
      windows: [],
      presentation: presentation()
    })),
    resources: {
      'genome-1': { ...genomeTemplate, name: 'genome-1.gbk' },
      'genome-2': { ...genomeTemplate, name: 'genome-2.gbk' },
      'genome-3': { ...genomeTemplate, name: 'genome-3.gbk' },
      'comparison-nucleotide-1': resource(
        'nucleotide-blast',
        'comparison-nucleotide-1.tsv',
        'NC_001416\tNC_001416\t100\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n'
      )
    },
    comparisons: [{
      kind: 'nucleotideBlast',
      resourceId: 'comparison-nucleotide-1',
      queryRecordIndex: 1,
      subjectRecordIndex: 2
    }],
    webFiles: {
      resourceOriginalNames: {
        'genome-1': 'alpha.gbk',
        'genome-2': 'middle.gbk',
        'genome-3': 'omega.gbk'
      }
    }
  });
  collision.generatedFileNameHints = new Map([
    ['generatedFiles.losat_blasts[0]', 'alpha.gbk']
  ]);
  const recipe = await buildSourceRecipe(collision);
  const info = buildRunInfo({ mode: 'linear', sourceRecipe: recipe });

  assert.equal(recipe.available, true, recipe.unavailableReason);
  await assertCliParserAccepts(recipe);
  const generatedFiles = recipe.generatedFiles;
  assert.ok(Array.isArray(generatedFiles));
  const blast = generatedFiles.find(
    ({ slot }) => slot === 'generatedFiles.source_recipe.comparison_blasts[0]'
  );
  const comparisonsTable = generatedFiles.find(
    ({ slot }) => slot === 'generatedFiles.source_recipe.comparisons_table'
  );
  assert.equal(blast?.name, 'alpha-2.gbk');
  assert.equal(comparisonsTable?.name, 'comparisons.tsv');
  // An uploaded BLAST table Generate described by its File holds no text (OV-2612):
  // the recipe reads its bytes through the committed resources' reader.
  const blastText = Buffer.from(collision.resources['comparison-nucleotide-1'].data, 'base64').toString();
  const { encoding: _encoding, data: _data, ...blastFields } = collision.resources['comparison-nucleotide-1'];
  // The reader gives the File's cached bytes, which Generate then transfers to the Worker.
  const cachedBlastBytes = new TextEncoder().encode(blastText);
  const fileBackedResources = { ...collision.resources,
    'comparison-nucleotide-1': fileBackedDescriptor(blastFields, async () => cachedBlastBytes) };
  const fileBacked = await buildSourceRecipe({ ...collision, resources: fileBackedResources,
    readResourceFileBytes: (id) => readCanonicalResourceFileBytes(fileBackedResources, id) });
  structuredClone(cachedBlastBytes.buffer, { transfer: [cachedBlastBytes.buffer] });
  assert.equal(fileBacked.available, true, fileBacked.unavailableReason);
  const fileBackedBlast = fileBacked.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.comparison_blasts[0]');
  assert.equal(new TextDecoder().decode(fileBackedBlast.data), blastText);
  assert.equal(Object.hasOwn(fileBackedBlast, 'resourceId'), false);
  assert.equal((await buildSourceRecipe({ ...collision, resources: fileBackedResources })).available, false,
    'without a reader of its bytes the recipe is unavailable');
  const [header, row] = comparisonsTable.data.trim().split('\n').map((line) => line.split('\t'));
  assert.equal(row[header.indexOf('blast')], 'alpha-2.gbk');
  assert.doesNotMatch(comparisonsTable.data, /^alpha\.gbk\t/m);
  assert.match(info.sourceRecipe.command, /--comparisons_table comparisons\.tsv/);
  assert.deepEqual(
    info.sourceRecipe.helperFiles.map(({ name }) => name).sort(),
    ['alpha-2.gbk', 'comparisons.tsv']
  );
  assert.deepEqual(
    generatedFiles.map(({ name }) => name).sort(),
    info.sourceRecipe.helperFiles.map(({ name }) => name).sort()
  );
  assert.doesNotMatch(
    `${info.sourceRecipe.command}\n${comparisonsTable.data}\n${generatedFiles.map(({ name }) => name).join('\n')}`,
    /comparison-nucleotide-1|run-info-resource/
  );
  await materializeAndRunRecipe(collision, recipe, info, 'collision-dependent-helper');
});

const circularCanonical = canonical({
  mode: 'circular',
  records: [{
    recordKey: 'record-1',
    cardinality: 'exactly_one',
    source: { kind: 'genbank', resourceId: 'record-1-genbank' },
    selector: null,
    windows: [],
    presentation: presentation()
  }],
  resources: {
    'record-1-genbank': resource(
      'genbank', 'record-1-genbank-input_file.gbk', 'LOCUS       input\n//\n'
    )
  },
  webFiles: {
    circularOutputPrefixExplicit: true,
    resourceOriginalNames: { 'record-1-genbank': 'input file.gbk' },
    circularInputOriginalName: 'input file.gbk'
  }
});

test('circular source recipe names comparison genomes with --conservation_sequence', async () => {
  const session = structuredClone(circularCanonical);
  session.renderRequest.diagramOptions.conservationBlastFiles = [
    { resourceId: 'conservation-blast-files-1', representation: 'file' }
  ];
  session.renderRequest.diagramOptions.conservationFastaFiles = [
    [{ resourceId: 'conservation-fasta-files-1', representation: 'file' }]
  ];
  session.resources['conservation-blast-files-1'] = resource('conservation-blast-file', 'ring.tsv', '');
  session.resources['conservation-fasta-files-1'] = resource('conservation-fasta-file', 'ring.gb', 'LOCUS       ring\n//\n');
  const sourceRecipe = await buildSourceRecipe(session);
  assert.equal(sourceRecipe.available, true, sourceRecipe.unavailableReason);
  assert.ok(sourceRecipe.args.includes('--conservation_sequence'));
  assert.ok(!sourceRecipe.args.includes('--conservation_fasta'));
  await assertCliParserAccepts(sourceRecipe);
});

test('circular source recipe gives a ring of several genome files one --conservation_genome', async () => {
  const session = structuredClone(circularCanonical);
  const ref = (resourceId) => ({ resourceId, representation: 'file' });
  session.renderRequest.diagramOptions.conservationBlastFiles = [ref('conservation-blast-files-1'), ref('conservation-blast-files-2')];
  session.renderRequest.diagramOptions.conservationFastaFiles = [
    [ref('conservation-fasta-files-1-1'), ref('conservation-fasta-files-1-2')],
    []
  ];
  session.resources['conservation-blast-files-1'] = resource('conservation-blast-file', 'ring.tsv', '');
  session.resources['conservation-blast-files-2'] = resource('conservation-blast-file', 'ring2.tsv', '');
  session.resources['conservation-fasta-files-1-1'] = resource('conservation-fasta-file', 'chromosome.gb', 'LOCUS       c\n//\n');
  session.resources['conservation-fasta-files-1-2'] = resource('conservation-fasta-file', 'plasmid.gb', 'LOCUS       p\n//\n');
  const sourceRecipe = await buildSourceRecipe(session);
  assert.equal(sourceRecipe.available, true, sourceRecipe.unavailableReason);
  const visible = sourceRecipe.args.map((token) => sourceRecipe.fileMetadata.get(token)?.name || token);
  const first = visible.indexOf('--conservation_genome');
  const [flag, chromosome, plasmid, secondFlag, noGenome] = visible.slice(first, first + 5);
  assert.deepEqual([flag, secondFlag, noGenome], ['--conservation_genome', '--conservation_genome', '']);
  assert.ok(chromosome && plasmid && chromosome !== plasmid);
  assert.ok(!visible.includes('--conservation_sequence'));
  await assertCliParserAccepts(sourceRecipe);
});

{
  const sourceRecipe = await buildSourceRecipe(circularCanonical);
  assert.equal(sourceRecipe.available, true);
  await assertCliParserAccepts(sourceRecipe);
  const info = buildRunInfo({
    mode: 'circular',
    sourceRecipe,
    exactReplayArgs: ['--session', '/canonical-render-session.gbdraw-session.json'],
    fileMetadata: {
      '/canonical-render-session.gbdraw-session.json': {
        name: "Bob's replay.gbdraw-session.json",
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      }
    },
    elapsedMs: 321,
    resultCount: 1
  });

  assert.match(info.sourceRecipe.command, /--gbk 'input file\.gbk'/);
  assert.equal(info.command, info.sourceRecipe.command);
  assert.equal(
    info.exactReplay.command,
    "gbdraw circular --session 'Bob'\\''s replay.gbdraw-session.json'"
  );
  assert.notEqual(info.sourceRecipe.command, info.exactReplay.command);
  assert.equal(info.elapsedMs, 321);
  assert.equal(info.resultCount, 1);
  assert.equal(info.reproducibility.level, 'exact-uploaded-files');
}

{
  const apostropheName = structuredClone(circularCanonical);
  apostropheName.webFiles.resourceOriginalNames['record-1-genbank'] = "Bob's genome.gbk";
  const info = buildRunInfo({
    mode: 'circular',
    sourceRecipe: await buildSourceRecipe(apostropheName)
  });
  assert.match(info.command, /--gbk 'Bob'\\''s genome\.gbk'/);
}

{
  const multiRecord = canonical({
    mode: 'circular',
    records: ['alpha', 'beta'].map((recordId, index) => ({
      recordKey: `record-${index + 1}`,
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'multi-genbank' },
      selector: { kind: 'recordId', value: recordId },
      windows: [],
      presentation: {
        ...presentation(),
        label: recordId.toUpperCase()
      }
    })),
    resources: {
      'multi-genbank': resource(
        'genbank',
        'multi-genbank-two-records.gbk',
        'LOCUS       alpha\n//\nLOCUS       beta\n//\n'
      )
    },
    webFiles: {
      resourceOriginalNames: { 'multi-genbank': 'two-records.gbk' }
    }
  });
  multiRecord.renderRequest.layout = {
    multiRecordSizeMode: 'equal',
    multiRecordPositions: ['alpha@1', 'beta@2']
  };
  const sourceRecipe = await buildSourceRecipe(multiRecord);

  assert.equal(sourceRecipe.available, true, sourceRecipe.unavailableReason);
  await assertCliParserAccepts(sourceRecipe);
  assert.ok(sourceRecipe.args.includes('--records_table'));
  assert.ok(sourceRecipe.args.includes('--multi_record_canvas'));
  assert.equal(sourceRecipe.args.includes('--multi_record_position'), false);
  const recordsTable = sourceRecipe.generatedFiles.find((file) => file.name === 'records.tsv');
  assert.ok(recordsTable);
  assert.match(recordsTable.data, /two-records\.gbk\tALPHA\t\talpha\t\t0\t1\t1\t/);
  assert.match(recordsTable.data, /two-records\.gbk\tBETA\t\tbeta\t\t0\t2\t2\t/);

  // One circle: the recipe asks for one circle with the stored settings, not for the grid (D-45 A).
  multiRecord.renderRequest.layout = {
    multiRecordSizeMode: 'equal', recordLayout: 'one_circle', oneCircleGapDegrees: 2.5, oneCircleMinLength: 0,
    oneCircleRecordBand: false
  };
  const oneCircleRecipe = await buildSourceRecipe(multiRecord);
  assert.equal(oneCircleRecipe.available, true, oneCircleRecipe.unavailableReason);
  await assertCliParserAccepts(oneCircleRecipe);
  assert.equal(oneCircleRecipe.args.includes('--multi_record_canvas'), false);
  const oneCircleArgs = oneCircleRecipe.args.slice(oneCircleRecipe.args.indexOf('--one_circle'));
  assert.deepEqual(oneCircleArgs.slice(0, 6), ['--one_circle', '--one_circle_gap_degrees', '2.5',
    '--one_circle_min_length', '0', '--one_circle_no_record_band']);
}

const linearCanonical = canonical({
  mode: 'linear',
  records: [
    {
      recordKey: 'alpha',
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'record-1-genbank' },
      selector: null,
      windows: [],
      presentation: presentation()
    },
    {
      recordKey: 'beta',
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'record-2-genbank' },
      selector: null,
      windows: [],
      presentation: presentation()
    }
  ],
  resources: {
    'record-1-genbank': resource(
      'genbank', 'record-1-genbank-alpha.gbk', 'LOCUS       alpha\n//\n'
    ),
    'record-2-genbank': resource(
      'genbank', 'record-2-genbank-beta.gbk', 'LOCUS       beta\n//\n'
    ),
    'comparison-nucleotide-1': resource('nucleotide-blast', 'comparison-nucleotide-1-alpha_vs_beta.tsv')
  },
  comparisons: [{
    kind: 'nucleotideBlast',
    resourceId: 'comparison-nucleotide-1',
    queryRecordIndex: 0,
    subjectRecordIndex: 1
  }],
  webFiles: {
    resourceOriginalNames: {
      'record-1-genbank': 'alpha.gbk',
      'record-2-genbank': 'beta.gbk',
      'comparison-nucleotide-1': 'alpha_vs_beta.tsv'
    }
  }
});

{
  const positioned = structuredClone(linearCanonical);
  positioned.renderRequest.records.forEach((record) => {
    record.display = { isCircular: null, startCoordinate: null };
  });
  positioned.renderRequest.layout = {
    recordGapPx: 24,
    multiRecordPositions: ['#1@1', '#2@2'],
    recordTranslations: [],
    similarityAlignment: null
  };
  const directPositions = await buildSourceRecipe(positioned);
  assert.equal(directPositions.available, true, directPositions.unavailableReason);
  await assertCliParserAccepts(directPositions);
  assert.deepEqual(directPositions.args.filter((arg) => /^#[12]@/.test(arg)), ['#1@1', '#2@2']);
  positioned.renderRequest.records[0].presentation.gridColumn = 1;
  const tablePositions = await buildSourceRecipe(positioned);
  assert.equal(tablePositions.available, true, tablePositions.unavailableReason);
  assert.match(tablePositions.generatedFiles.find((file) => file.name === 'records.tsv').data, /alpha\.gbk.*\t1\t1\t\t\n/);
  positioned.renderRequest.layout.recordTranslations = [
    { recordKey: 'alpha', x: 1, y: 0 }, { recordKey: 'beta', x: 0, y: 0 }
  ];
  const shifted = await buildSourceRecipe(positioned);
  assert.equal(shifted.available, false);
  assert.match(shifted.unavailableReason, /record position offsets/);
  const partialRows = structuredClone(linearCanonical);
  partialRows.renderRequest.records[0].presentation.gridRow = 2;
  const partialRecipe = await buildSourceRecipe(partialRows);
  assert.equal(partialRecipe.available, true, partialRecipe.unavailableReason);
  assert.ok(partialRecipe.args.includes('#1@2'));
  assert.equal(partialRecipe.args.includes('#2@2'), false);
}

{
  const direct = await buildSourceRecipe(linearCanonical);
  assert.equal(direct.available, true);
  await assertCliParserAccepts(direct);
  const info = buildRunInfo({
    mode: 'linear',
    sourceRecipe: direct,
    exactReplayArgs: ['--session', '/linear-canonical-session'],
    fileMetadata: {
      '/linear-canonical-session': {
        name: 'comparison.gbdraw-session.json',
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      }
    }
  });
  assert.match(info.command, /--gbk alpha\.gbk beta\.gbk/);
  assert.match(info.command, /-b alpha_vs_beta\.tsv/);
  assert.equal(
    info.exactReplay.command,
    'gbdraw linear --session comparison.gbdraw-session.json'
  );

  const restored = structuredClone(linearCanonical);
  restored.webFiles.bindings = {
    schema: 1,
    linearSeqs: [
      { gb: { resourceId: 'record-1-genbank', name: 'alpha.gbk' } },
      { gb: { resourceId: 'record-2-genbank', name: 'beta.gbk' } }
    ],
    linearComparisons: [{
      id: 'edge-1',
      file: { resourceId: 'comparison-nucleotide-1', name: 'alpha_vs_beta.tsv' }
    }]
  };
  delete restored.webFiles.resourceOriginalNames;
  const roundTrip = await buildSourceRecipe(restored);
  assert.equal(roundTrip.available, true);
  assert.deepEqual(roundTrip.args, direct.args);
}

{
  const generatedHelper = structuredClone(linearCanonical);
  delete generatedHelper.webFiles.resourceOriginalNames['comparison-nucleotide-1'];
  generatedHelper.resources['comparison-nucleotide-1'].name =
    'comparison-nucleotide-1-browser-comparison.tsv';
  generatedHelper.generatedFileNameHints = new Map([
    ['generatedFiles.losat_blasts[0]', 'browser-comparison.tsv']
  ]);
  const sourceRecipe = await buildSourceRecipe(generatedHelper);
  const info = buildRunInfo({ mode: 'linear', sourceRecipe });

  assert.equal(sourceRecipe.available, true);
  assert.match(info.command, /-b browser-comparison\.tsv/);
  assert.doesNotMatch(info.command, /comparison-nucleotide-1-browser/);
  assert.equal(sourceRecipe.generatedFiles.length, 1);
  assert.equal(sourceRecipe.generatedFiles[0].name, 'browser-comparison.tsv');
  assert.ok(info.sourceRecipe.helperFiles.every(({ name }) => !name.includes('comparison-nucleotide-1')));
  assert.equal(info.reproducibility.level, 'session-recommended');
  assert.match(info.reproducibility.notes.join('\n'), /Download reproducibility files/);

  const invalidHint = structuredClone(generatedHelper);
  invalidHint.generatedFileNameHints = new Map([
    ['generatedFiles.losat_blasts[0]', '..']
  ]);
  const unavailable = await buildSourceRecipe(invalidHint);
  assert.equal(unavailable.available, false);
  assert.match(unavailable.unavailableReason, /trustworthy visible filename/);
}

{
  const duplicateUploads = structuredClone(linearCanonical);
  duplicateUploads.webFiles.resourceOriginalNames['record-1-genbank'] = 'Genome.gbk';
  duplicateUploads.webFiles.resourceOriginalNames['record-2-genbank'] = 'genome.GBK';
  const sourceRecipe = await buildSourceRecipe(duplicateUploads);

  assert.equal(sourceRecipe.available, false);
  assert.match(sourceRecipe.unavailableReason, /multiple required files.*genome\.GBK/i);
}

{
  const sourceRecipe = {
    available: true,
    mode: 'linear',
    args: ['--gbk', '/uploaded', '-b', '/generated-a', '/generated-b'],
    fileMetadata: new Map([
      ['/uploaded', { name: 'Records.tsv', slot: 'resources.uploaded', kind: 'uploaded' }],
      ['/generated-a', { name: 'records.tsv', slot: 'generatedFiles.first', kind: 'generated' }],
      ['/generated-b', { name: 'records.tsv', slot: 'generatedFiles.second', kind: 'generated' }]
    ]),
    generatedFiles: [
      { path: '/generated-a', name: 'records.tsv', slot: 'generatedFiles.first', data: 'a' },
      { path: '/generated-b', name: 'records.tsv', slot: 'generatedFiles.second', data: 'b' }
    ]
  };
  const info = buildRunInfo({
    mode: 'linear',
    sourceRecipe,
    exactReplayArgs: ['--session', '/exact'],
    fileMetadata: new Map([
      ['/exact', {
        name: 'records-2.tsv',
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      }]
    ])
  });

  assert.deepEqual(
    info.sourceRecipe.commandArgs.slice(-3),
    ['-b', 'records-2.tsv', 'records-3.tsv']
  );
  assert.deepEqual(
    info.sourceRecipe.helperFiles.map(({ name }) => name),
    ['records-2.tsv', 'records-3.tsv']
  );
  assert.deepEqual(info.exactReplay.helperFiles.map(({ name }) => name), ['records-2-2.tsv']);
  assert.match(info.exactReplay.command, /--session records-2-2\.tsv$/);
  const allNames = info.helperFiles.map(({ name }) => name);
  assert.equal(new Set(allNames.map((name) => name.toLowerCase())).size, allNames.length);
}

{
  const cardinalityAll = structuredClone(linearCanonical);
  cardinalityAll.renderRequest.records = [cardinalityAll.renderRequest.records[0]];
  cardinalityAll.renderRequest.records[0].cardinality = 'all';
  cardinalityAll.renderRequest.comparisons = [];
  delete cardinalityAll.resources['record-2-genbank'];
  delete cardinalityAll.resources['comparison-nucleotide-1'];
  delete cardinalityAll.webFiles.resourceOriginalNames['record-2-genbank'];
  delete cardinalityAll.webFiles.resourceOriginalNames['comparison-nucleotide-1'];
  const sourceRecipe = await buildSourceRecipe(cardinalityAll);

  assert.equal(sourceRecipe.available, true, sourceRecipe.unavailableReason);
  assert.ok(sourceRecipe.semanticCoverage.consumedPaths.includes('records[0].cardinality'));
}

{
  const gridColumn = structuredClone(linearCanonical);
  gridColumn.renderRequest.records.forEach((record, index) => {
    record.presentation.gridRow = 1;
    record.presentation.gridColumn = index + 1;
  });
  const sourceRecipe = await buildSourceRecipe(gridColumn);

  assert.equal(sourceRecipe.available, true, sourceRecipe.unavailableReason);
  const table = sourceRecipe.generatedFiles.find((file) => file.name === 'records.tsv');
  assert.ok(table);
  // Columns `topology` and `display_start` follow, empty.
  assert.match(table.data, /alpha\.gbk.*\t1\t1\t\t\n/);
  assert.match(table.data, /beta\.gbk.*\t1\t2\t\t\n/);
}

{
  const unboundGeneratedProtein = structuredClone(linearCanonical);
  unboundGeneratedProtein.renderRequest.comparisons = [{
    kind: 'generatedProteinComparison',
    mode: 'pairwise',
    pairs: [{ queryRecordIndex: 0, subjectRecordIndex: 1 }],
    settings: {}
  }];
  const sourceRecipe = await buildSourceRecipe(unboundGeneratedProtein);

  assert.equal(sourceRecipe.available, false);
  assert.match(sourceRecipe.unavailableReason, /comparison.*lossless|binding/i);
}

{
  const unknownSemanticField = structuredClone(circularCanonical);
  unknownSemanticField.renderRequest.records[0].presentation.futurePlacement = 'spiral';
  const sourceRecipe = await buildSourceRecipe(unknownSemanticField);

  assert.equal(sourceRecipe.available, false);
  assert.match(sourceRecipe.unavailableReason, /presentation\.futurePlacement/);

  const displayMetadata = structuredClone(circularCanonical);
  displayMetadata.webFiles.displayOnlyNote = 'not part of render semantics';
  assert.equal((await buildSourceRecipe(displayMetadata)).available, true);
}

{
  const batch = structuredClone(circularCanonical);
  batch.renderRequest.grouping = 'batch';
  batch.renderRequest.output = [
    { ...batch.renderRequest.output, prefix: 'first' },
    { ...batch.renderRequest.output, prefix: 'second' }
  ];
  const sourceRecipe = await buildSourceRecipe(batch);

  assert.equal(sourceRecipe.available, false);
  assert.match(sourceRecipe.unavailableReason, /batch|output/i);
}

test('a one-record Circular batch, as the CLI writes it, has the Source recipe of its one output (OV-3191)', async () => {
  const batch = structuredClone(circularCanonical);
  batch.renderRequest.grouping = 'batch';
  batch.renderRequest.output = [{ ...batch.renderRequest.output, prefix: 'yeast_c' }];
  const single = structuredClone(circularCanonical);
  single.renderRequest.output.prefix = 'yeast_c';
  const recipe = await buildSourceRecipe(batch);

  assert.equal(recipe.available, true, recipe.unavailableReason);
  assert.deepEqual(recipe.args, (await buildSourceRecipe(single)).args);
  assert.ok(recipe.semanticCoverage.consumedPaths.includes('output[0].prefix'));
  await assertCliParserAccepts(recipe);

  const unknownOutputField = structuredClone(batch);
  unknownOutputField.renderRequest.output[0].futureName = 'x';
  assert.match((await buildSourceRecipe(unknownOutputField)).unavailableReason, /output\[0\]\.futureName/);
});

{
  const missingProvenance = structuredClone(circularCanonical);
  missingProvenance.format = 'gbdraw-session';
  missingProvenance.version = 32;
  missingProvenance.renderRequest.schema = 2;
  missingProvenance.webFiles = {};
  const sourceRecipe = await buildSourceRecipe(missingProvenance);
  const info = buildRunInfo({
    mode: 'circular',
    sourceRecipe,
    exactReplayArgs: ['--session', '/canonical-render-session.gbdraw-session.json'],
    fileMetadata: {
      '/canonical-render-session.gbdraw-session.json': {
        name: 'legacy-replay.gbdraw-session.json',
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      }
    }
  });

  assert.equal(sourceRecipe.available, false);
  assert.match(sourceRecipe.unavailableReason, /source-file provenance/i);
  assert.equal(info.sourceRecipe.command, '');
  assert.match(info.sourceRecipe.unavailableReason, /source-file provenance/i);
  assert.equal(
    info.exactReplay.command,
    'gbdraw circular --session legacy-replay.gbdraw-session.json'
  );
  assert.equal(reproducibilityLabel(info.reproducibility.level), 'Source recipe unavailable');
}

{
  const info = buildRunInfo({
    mode: 'circular',
    args: ['-o', 'my diagram', '-d', '/combined_d.tsv', '--track_type', 'middle', '--gbk', '/input.gb'],
    fileMetadata: new Map([
      ['/input.gb', { name: 'input file.gbk', slot: 'files.c_gb', kind: 'uploaded' }],
      ['/combined_d.tsv', { name: 'combined_d.tsv', slot: 'generatedFiles.combined_d', kind: 'generated' }]
    ]),
    elapsedMs: 1234.5,
    resultCount: 1,
    startedAtIso: '2026-06-23T00:00:00.000Z'
  });

  assert.equal(
    info.command,
    "gbdraw circular -o 'my diagram' -d combined_d.tsv --track_type middle --gbk 'input file.gbk'"
  );
  assert.equal(info.invocation.args.includes('-f'), false);
  assert.deepEqual(info.helperFiles, [
    { path: '/combined_d.tsv', name: 'combined_d.tsv', slot: 'generatedFiles.combined_d' }
  ]);
  assert.equal(info.reproducibility.level, 'requires-helper-files');
  assert.equal(reproducibilityLabel(info.reproducibility.level), 'Source recipe needs helper files');
  assert.match(info.reproducibility.notes.join('\n'), /combined_d\.tsv/);
  assert.match(info.reproducibility.notes.join('\n'), /a \.gbdraw file/);
  assert.equal(info.invocation.sessionExportable, false);
}

{
  const info = buildRunInfo({
    mode: 'linear',
    sourceRecipe: {
      available: false,
      args: [],
      fileMetadata: new Map(),
      generatedFiles: [],
      unavailableReason: 'Source recipe unavailable: source-file provenance is missing.',
      unavailableReasonParts: [{ text: 'Source recipe unavailable: source-file provenance is missing.', verbatim: true }]
    },
    exactReplayArgs: ['--session', '/canonical-render-session.gbdraw-session.json'],
    fileMetadata: {
      '/canonical-render-session.gbdraw-session.json': {
        name: 'comparison.gbdraw-session.json',
        slot: 'generatedFiles.canonical_render_session',
        kind: 'generated'
      }
    },
    elapsedMs: 500,
    resultCount: 1
  });

  assert.equal(
    info.exactReplay.command,
    'gbdraw linear --session comparison.gbdraw-session.json'
  );
  assert.equal(info.reproducibility.level, 'source-unavailable');
  assert.equal(
    reproducibilityLabel(info.reproducibility.level),
    'Source recipe unavailable'
  );
  assert.match(info.reproducibility.notes.join('\n'), /exact replay/i);
}

{
  const info = buildRunInfo({
    mode: 'linear',
    args: ['--gbk', '/seq_0.gb', '/seq_1.gb', '-b', '/blast_0.txt'],
    fileMetadata: {
      '/seq_0.gb': { name: 'alpha.gb', slot: 'files.linearSeqs[0].gb', kind: 'uploaded' },
      '/seq_1.gb': { name: 'beta.gb', slot: 'files.linearSeqs[1].gb', kind: 'uploaded' },
      '/blast_0.txt': { name: 'alpha_vs_beta.tsv', slot: 'files.linearSeqs[0].blast', kind: 'uploaded' }
    },
    elapsedMs: 900,
    resultCount: 2
  });

  assert.equal(info.command, 'gbdraw linear --gbk alpha.gb beta.gb -b alpha_vs_beta.tsv');
  assert.deepEqual(info.invocation.fileBindings, [
    { argIndex: 1, slot: 'files.linearSeqs[0].gb', name: 'alpha.gb' },
    { argIndex: 2, slot: 'files.linearSeqs[1].gb', name: 'beta.gb' },
    { argIndex: 4, slot: 'files.linearSeqs[0].blast', name: 'alpha_vs_beta.tsv' }
  ]);
  assert.equal(info.reproducibility.level, 'exact-uploaded-files');
  assert.equal(reproducibilityLabel(info.reproducibility.level), 'Source recipe ready');
  assert.equal(info.invocation.sessionExportable, true);
}

{
  const info = buildRunInfo({
    mode: 'linear',
    args: ['--gbk', '/seq_0.gb', '/seq_1.gb', '-b', '/blast_0.txt'],
    fileMetadata: {
      '/seq_0.gb': { name: 'alpha.gb', slot: 'files.linearSeqs[0].gb', kind: 'uploaded' },
      '/seq_1.gb': { name: 'beta.gb', slot: 'files.linearSeqs[1].gb', kind: 'uploaded' },
      '/blast_0.txt': { name: 'alpha_vs_beta.tsv', slot: 'generatedFiles.losat_blasts[0]', kind: 'generated' }
    },
    elapsedMs: 900,
    resultCount: 1
  });

  assert.equal(info.reproducibility.level, 'session-recommended');
  assert.equal(reproducibilityLabel(info.reproducibility.level), 'Source recipe needs helper files');
  assert.match(info.reproducibility.notes.join('\n'), /Download reproducibility files/);
  assert.match(info.reproducibility.notes.join('\n'), /save the Session \(a \.gbdraw file\) and load it again/);
  assert.equal(info.sessionCommand, '');
  assert.equal(info.invocation.sessionExportable, false);
}


// TR-07: a slot legend label that the CLI slot grammar cannot carry must make
// the recipe unavailable instead of producing a command the CLI rejects.
test('slot legend labels with a comma or an inline comment marker stay lossless (TR-07)', async () => {
  const outcomes = [];
  for (const label of ['GC skew (1 kb, AT-rich)', 'a #b']) {
    const session = structuredClone(circularGallerySession);
    session.renderRequest.diagramOptions.tracks.circularTrackSlots = [{
      kind: 'circularTrackSlot', id: 'ticks', renderer: 'ticks', enabled: true, side: null,
      radius: null, width: null, z: 0,
      params: { tick_label_layout: 'label_out_tick_in', legend_label: label },
      innerGapPx: null, outerGapPx: null
    }];
    session.renderRequest.diagramOptions.tracks.circularTrackAxisIndex = 0;
    const recipe = await buildSourceRecipe(session);
    outcomes.push(recipe.available
      ? { label, token: recipe.args[recipe.args.indexOf('--circular_track_slot') + 1] }
      : { label, unavailable: true });
  }
  const tokens = outcomes.filter(({ token }) => token).map(({ token }) => token);
  const labels = await evaluatePython([
    'from gbdraw.tracks.circular import parse_circular_track_slot',
    'def label(token):',
    '    try:',
    '        return parse_circular_track_slot(token).params.get("legend_label")',
    '    except ValueError as error:',
    '        return f"rejected: {error}"',
    'result = [label(token) for token in payload]'
  ].join('\n'), tokens);
  outcomes.filter(({ token }) => token).forEach((outcome, index) => { outcome.parsed = labels[index]; });
  for (const outcome of outcomes) {
    assert.ok(outcome.unavailable || outcome.parsed === outcome.label, JSON.stringify(outcome));
  }
  assert.deepEqual(outcomes.map(({ unavailable }) => Boolean(unavailable)), [true, true]);
  // A label the CLI grammar can carry still yields a lossless token.
  const plain = structuredClone(circularGallerySession);
  plain.renderRequest.diagramOptions.tracks.circularTrackSlots = [{
    kind: 'circularTrackSlot', id: 'ticks', renderer: 'ticks', enabled: true, side: null,
    radius: null, width: null, z: 0,
    params: { tick_label_layout: 'label_out_tick_in', legend_label: 'GC skew 1 kb#AT=rich' },
    innerGapPx: null, outerGapPx: null
  }];
  plain.renderRequest.diagramOptions.tracks.circularTrackAxisIndex = 0;
  const plainRecipe = await buildSourceRecipe(plain);
  assert.equal(plainRecipe.available, true, plainRecipe.unavailableReason);
});

// GE-03 (Linear): with the scale font set and no ruler-label font, the CLI makes
// the ruler labels follow the scale font while the Web keeps the configured
// ruler-label defaults, so the recipe is unavailable with a reason.
// OV-490: Python reads a Linear slot that names no dinucleotide with the request's.
test('a Linear slot without a dinucleotide takes the request\'s in the source recipe (OV-490)', async () => {
  const session = structuredClone(linearGallerySession);
  session.renderRequest.diagramOptions.dinucleotide = 'AT';
  const slot = (id, renderer, side) => ({ kind: 'linearTrackSlot', id, renderer, enabled: true, side, height: null, spacing: null, z: 0, params: {} });
  session.renderRequest.diagramOptions.tracks.linearTrackSlots = [slot('features', 'features', 'overlay'), slot('skew', 'dinucleotide_skew', 'below')];
  session.renderRequest.diagramOptions.tracks.linearTrackAxisIndex = 0;
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const tokens = recipe.args.filter((_, index) => recipe.args[index - 1] === '--linear_track_slot');
  assert.match(tokens.find((token) => token.startsWith('skew:')) || '', /[@,]nt=AT(,|$)/, tokens.join(' | '));
});

// OV-3157: Python reads a Linear slot that names no side from its place beside the Axis index. The recipe writes
// that side, so its command draws the request's stack instead of failing on a side that conflicts with the index.
test('a slot without a side takes the side or lane of its place in the source recipe (OV-3157, OV-3186)', async () => {
  const session = structuredClone(linearGallerySession);
  const slot = (id, renderer) => ({ kind: 'linearTrackSlot', id, renderer, enabled: true, side: null, height: null, spacing: null, z: 0, params: {} });
  session.renderRequest.diagramOptions.tracks.linearTrackSlots = [
    slot('gc', 'dinucleotide_content'), slot('features', 'features'), slot('sk', 'dinucleotide_skew')
  ];
  session.renderRequest.diagramOptions.tracks.linearTrackAxisIndex = 2;
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const tokens = recipe.args.filter((_, index) => recipe.args[index - 1] === '--linear_track_slot');
  assert.deepEqual(tokens.map((token) => [token.split(':')[0], (token.match(/[@,]side=(\w+)/) || [])[1]]),
    [['gc', 'above'], ['features', 'above'], ['sk', 'below']], tokens.join(' | '));
  await assertCliParserAccepts(recipe);
  // OV-3186: a Circular features slot without a side or a lane takes the lane of its place, under every track type
  // (Middle wrote `split`, a different figure; Tuckin wrote `inside`, which the CLI refuses beside this index).
  for (const trackType of ['middle', 'tuckin', 'spreadout']) {
    const circular = structuredClone(circularGallerySession);
    const ring = (id, renderer) => ({ kind: 'circularTrackSlot', id, renderer, enabled: true, side: null, radius: null, width: null,
      z: 0, params: {}, innerGapPx: null, outerGapPx: null });
    circular.renderRequest.diagramOptions.configOverrides['canvas.circular.track_type'] = trackType;
    circular.renderRequest.diagramOptions.tracks.circularTrackSlots = [
      ring('gc', 'dinucleotide_content'), ring('features', 'features'), ring('sk', 'dinucleotide_skew')
    ];
    circular.renderRequest.diagramOptions.tracks.circularTrackAxisIndex = 2;
    const circularRecipe = await buildSourceRecipe(circular);
    assert.equal(circularRecipe.available, true, circularRecipe.unavailableReason);
    const features = circularRecipe.args.find((arg, index) => circularRecipe.args[index - 1] === '--circular_track_slot' && arg.startsWith('features:'));
    assert.equal((features.match(/lane_direction=(\w+)/) || [])[1], 'outside', `${trackType}: ${features}`);
    await assertCliParserAccepts(circularRecipe);
  }
});

// PD-OI-078 (rule L): a request that carries the track narrowing settings (a Python
// `config_overrides` Session) writes them as the CLI flags that set them.
test('the Circular track narrowing settings have their CLI flags in the Source recipe', async () => {
  for (const [narrowing, minPx] of [[false, 44], [true, 30]]) {
    const session = structuredClone(circularGallerySession);
    Object.assign(session.renderRequest.diagramOptions.configOverrides, {
      'canvas.circular.track_narrowing': narrowing,
      'canvas.circular.track_narrowing_min_px': minPx
    });
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    assert.equal(recipe.args[recipe.args.indexOf('--track_narrowing_min_px') + 1], String(minPx), JSON.stringify(recipe.args));
    assert.equal(recipe.args.includes('--no_track_narrowing'), !narrowing);
    await assertCliParserAccepts(recipe);
  }
});

// A request with italic gene names writes `--italic_gene_labels`; one without the override writes nothing.
test('Italic gene names has its CLI flag in the Source recipe of either mode', async () => {
  for (const [session, mode] of [[circularGallerySession, 'circular'], [linearGallerySession, 'linear']]) {
    for (const italic of [true, undefined]) {
      const copy = structuredClone(session);
      const overrides = copy.renderRequest.diagramOptions.configOverrides;
      delete overrides[`labels.${mode}.italic_gene_labels`];
      if (italic !== undefined) overrides[`labels.${mode}.italic_gene_labels`] = italic;
      const recipe = await buildSourceRecipe(copy);
      assert.equal(recipe.available, true, recipe.unavailableReason);
      assert.equal(recipe.args.includes('--italic_gene_labels'), italic === true, JSON.stringify(recipe.args));
      await assertCliParserAccepts(recipe);
    }
  }
});

// Wrap Rows writes `--wrap_length` with the request's value, Auto as `auto`; no wrap length, null and "none" (off)
// write nothing. Any other stored text (a Python API Session's "5kb", " 20000 ") is not written as typed: the recipe
// is unavailable with the reason, since the CLI parser may not read it back (review 3 F3, R12).
test('a Linear wrap length has its CLI flag in the Source recipe', async () => {
  const withLength = (length) => {
    const copy = structuredClone(linearGallerySession);
    const overrides = copy.renderRequest.diagramOptions.configOverrides;
    delete overrides['canvas.linear.wrap_length'];
    if (length !== undefined) overrides['canvas.linear.wrap_length'] = length;
    return copy;
  };
  for (const [length, flag] of [['auto', 'auto'], [20000, '20000'], [undefined, undefined], [null, undefined], ['none', undefined]]) {
    const recipe = await buildSourceRecipe(withLength(length));
    assert.equal(recipe.available, true, recipe.unavailableReason);
    const at = recipe.args.indexOf('--wrap_length');
    assert.equal(at === -1 ? undefined : recipe.args[at + 1], flag, JSON.stringify(recipe.args));
    await assertCliParserAccepts(recipe);
  }
  for (const length of ['5kb', ' 20000 ', 'Auto', 1.5]) {
    const recipe = await buildSourceRecipe(withLength(length));
    assert.equal(recipe.available, false, JSON.stringify(length));
    assert.match(recipe.unavailableReason, /wrap length/, JSON.stringify(length));
  }
});

// A figure size writes `--figure_width` and `--figure_height` with the request's value in either mode; Auto writes
// nothing (issue #1427).
test('a figure size has its CLI flags in the Source recipe of either mode', async () => {
  for (const session of [circularGallerySession, linearGallerySession]) {
    for (const [width, height] of [['180mm', '5in'], ['1200px', undefined], [undefined, undefined]]) {
      const copy = structuredClone(session);
      const overrides = copy.renderRequest.diagramOptions.configOverrides;
      if (width !== undefined) overrides['canvas.figure_width'] = width;
      if (height !== undefined) overrides['canvas.figure_height'] = height;
      const recipe = await buildSourceRecipe(copy);
      assert.equal(recipe.available, true, recipe.unavailableReason);
      for (const [flag, value] of [['--figure_width', width], ['--figure_height', height]]) {
        const at = recipe.args.indexOf(flag);
        assert.equal(at === -1 ? undefined : recipe.args[at + 1], value, JSON.stringify(recipe.args));
      }
      await assertCliParserAccepts(recipe);
    }
  }
});

// A fixed identity colour scale writes `--identity_scale MIN MAX` from the mode's leaf; no scale writes nothing.
test('a fixed identity scale has its CLI flag in the Source recipe of either mode', async () => {
  const leaves = { circular: 'objects.conservation.identity_scale', linear: 'objects.blast_match.identity_scale' };
  for (const [session, mode] of [[circularGallerySession, 'circular'], [linearGallerySession, 'linear']]) {
    for (const scale of [[40, 95.5], undefined]) {
      const copy = structuredClone(session);
      const overrides = copy.renderRequest.diagramOptions.configOverrides;
      delete overrides[leaves[mode]];
      if (scale !== undefined) overrides[leaves[mode]] = scale;
      const recipe = await buildSourceRecipe(copy);
      assert.equal(recipe.available, true, recipe.unavailableReason);
      const at = recipe.args.indexOf('--identity_scale');
      assert.deepEqual(at < 0 ? null : recipe.args.slice(at + 1, at + 3), scale ? ['40', '95.5'] : null, JSON.stringify(recipe.args));
      await assertCliParserAccepts(recipe);
    }
  }
});

// PD-OI-109: a request that hides the Circular definition writes `--hide_definition`;
// one that shows it (no override) writes nothing.
test('a hidden Circular definition has its CLI flag in the Source recipe', async () => {
  for (const show of [false, undefined]) {
    const session = structuredClone(circularGallerySession);
    const overrides = session.renderRequest.diagramOptions.configOverrides;
    delete overrides['objects.definition.circular.show'];
    if (show !== undefined) overrides['objects.definition.circular.show'] = show;
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    assert.equal(recipe.args.includes('--hide_definition'), show === false, JSON.stringify(recipe.args));
    await assertCliParserAccepts(recipe);
  }
});

// PD-OI-110: a ring without a color of its own ('') is written `auto`, which the CLI reads back as '';
// with no ring colored the flag is left out. The ring color scheme override has its CLI flag.
test('ring colors and the ring color scheme read back from the Source recipe', async () => {
  for (const [colors, scheme] of [[['', '#123456'], 'ten'], [['', ''], undefined]]) {
    const session = structuredClone(circularCanonical);
    const options = session.renderRequest.diagramOptions;
    options.conservationBlastFiles = colors.map((_, index) => ({ resourceId: `ring-${index}`, representation: 'file' }));
    options.conservationColors = colors;
    colors.forEach((_, index) => { session.resources[`ring-${index}`] = resource('conservation-blast-file', `ring${index}.tsv`, ''); });
    if (scheme) options.configOverrides = { ...(options.configOverrides || {}), 'objects.conservation.color_scheme': scheme };
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    await assertCliParserAccepts(recipe);
    const visibleArgs = recipe.args.map((token) => recipe.fileMetadata.get(token)?.name || token);
    const readBack = await evaluatePython([
      'from gbdraw.circular import _get_args, conservation_colors_from_cli',
      'args = _get_args(payload["args"])',
      'result = {"colors": conservation_colors_from_cli(args.conservation_colors), "scheme": args.conservation_color_scheme}'
    ].join('\n'), { args: visibleArgs });
    assert.deepEqual(readBack, { colors: colors.some(Boolean) ? colors : null, scheme: scheme ?? null }, JSON.stringify(recipe.args));
  }
});

test('a Linear scale font without a ruler-label font is not a lossless recipe (GE-03)', async () => {
  const session = structuredClone(linearGallerySession);
  Object.assign(session.renderRequest.diagramOptions.configOverrides, {
    'objects.scale.style': 'ruler',
    'objects.scale.font_size.short': 10,
    'objects.scale.font_size.long': 10
  });
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, false, JSON.stringify(recipe.args || []));
  assert.match(recipe.unavailableReason, /ruler/i);
  // Linked fonts send both sizes, so the CLI receives both flags.
  Object.assign(session.renderRequest.diagramOptions.configOverrides, {
    'objects.scale.ruler_label_font_size.short': 10,
    'objects.scale.ruler_label_font_size.long': 10
  });
  const linked = await buildSourceRecipe(session);
  assert.equal(linked.available, true, linked.unavailableReason);
  for (const flag of ['--scale_font_size', '--ruler_label_font_size']) {
    assert.equal(linked.args[linked.args.indexOf(flag) + 1], '10', JSON.stringify(linked.args));
  }
});

// Design Q4 (PR-Q4-1): identical features share an original-coordinate hash, so
// the placement row must name the complete biological feature ID with hash=.
// CI-02: a Circular record order changed with Up/Down is the order of
// `layout.multiRecordPositions`; the records table carries it as `column`, and
// the CLI draws the records in that order.
test('a reordered Circular multi-record canvas keeps its order in the records table recipe', async () => {
  const text = await readFile(join(repoRoot, 'tests', 'fixtures', 'web_batch_two_records.gb'), 'utf8');
  const session = canonical({
    mode: 'circular',
    records: ['TESTA', 'TESTB'].map((recordId) => ({
      recordKey: recordId,
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'source' },
      selector: { kind: 'recordId', value: recordId },
      windows: [],
      presentation: presentation()
    })),
    resources: { source: resource('genbank', 'source-two.gbk', text) },
    webFiles: { resourceOriginalNames: { source: 'two.gbk' } }
  });
  session.renderRequest.layout = { multiRecordSizeMode: 'equal', multiRecordPositions: ['#2@1', '#1@1'] };
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find((file) => file.name === 'records.tsv').data;
  assert.match(table, /two\.gbk\t\t\tTESTA\t\t0\t1\t1\t2\t\t\n/);
  assert.match(table, /two\.gbk\t\t\tTESTB\t\t0\t2\t1\t1\t\t\n/);
  const svg = await readFile(await materializeAndRunRecipe(
    session, recipe, buildRunInfo({ mode: recipe.mode, sourceRecipe: recipe }), 'reordered records'
  ), 'utf8');
  const drawn = [...new Set([...svg.matchAll(/data-gbdraw-record-id="([^"]+)"/g)].map((match) => match[1]))];
  assert.deepEqual(drawn, ['TESTB', 'TESTA']);
});

// OV-3313: the records table reads a Record Order token's row as the renderer reads it (`recordOrderRows`): `+1`,
// `01` and ` 1 ` are row 1; `1.0`, `0x1` and `1e0`, which the CLI refuses, give no records table.
test('the records table recipe reads a Record Order row as the renderer reads it', async () => {
  const text = await readFile(join(repoRoot, 'tests', 'fixtures', 'web_batch_two_records.gb'), 'utf8');
  const recipeOf = async (row) => {
    const session = canonical({
      mode: 'circular',
      records: ['TESTA', 'TESTB'].map((recordId) => ({
        recordKey: recordId,
        cardinality: 'exactly_one',
        source: { kind: 'genbank', resourceId: 'source' },
        selector: { kind: 'recordId', value: recordId },
        windows: [],
        presentation: presentation()
      })),
      resources: { source: resource('genbank', 'source-two.gbk', text) },
      webFiles: { resourceOriginalNames: { source: 'two.gbk' } }
    });
    session.renderRequest.layout = { multiRecordSizeMode: 'equal', multiRecordPositions: [`#2@${row}`, '#1@1'] };
    return buildSourceRecipe(session);
  };
  const table = (recipe) => recipe.generatedFiles.find((file) => file.name === 'records.tsv').data;
  const expected = table(await recipeOf('1'));
  // The row is read from the token as stored (review 1 L2): Python's `strip()` takes U+001F and U+0085, which `trim()`
  // leaves, and leaves U+FEFF, which `trim()` takes.
  for (const row of ['+1', '01', ' 1 ', '1\u001f', '1\u0085']) {
    const recipe = await recipeOf(row);
    assert.equal(recipe.available, true, `${row}: ${recipe.unavailableReason}`);
    assert.equal(table(recipe), expected, row);
  }
  for (const row of ['1.0', '0x1', '1e0', '1\ufeff']) assert.equal((await recipeOf(row)).available, false, row);
});

test('placement recipe names one of two identical features and the CLI replays it', async () => {
  const genbank = [
    'LOCUS       dup                      120 bp    DNA     linear   UNK 01-JAN-1980',
    'DEFINITION  Identical CDS fixture.',
    'FEATURES             Location/Qualifiers',
    '     CDS             21..80',
    '                     /locus_tag="first"',
    '     CDS             21..80',
    '                     /locus_tag="second"',
    'ORIGIN',
    '        1 atgatgatga tgatgatgat gatgatgatg atgatgatga tgatgatgat gatgatgatg',
    '       61 atgatgatga tgatgatgat gatgatgatg atgatgatga tgatgatgat gatgatgatg',
    '//',
    ''
  ].join('\n');
  const [first, second] = await sourceFeatureIds(genbank);
  assert.equal(first.replace(/~0$/, '~1'), second);
  const session = canonical({
    mode: 'circular',
    records: [{
      recordKey: 'dup',
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'source' },
      selector: null,
      windows: [],
      presentation: presentation(),
      display: { isCircular: null, startCoordinate: null }
    }],
    resources: { source: resource('genbank', 'source-dup.gbk', genbank) },
    webFiles: { resourceOriginalNames: { source: 'dup.gbk' } }
  });
  session.renderRequest.diagramOptions.featurePlacements = [
    { recordKey: 'dup', biologicalFeatureId: second, placement: { kind: 'main' } }
  ];
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find(
    ({ slot }) => slot === 'generatedFiles.source_recipe.feature_placements'
  );
  assert.match(table.data, new RegExp(`\\n#1\\thash=${second}\\tmain\\t`));
  await materializeAndRunRecipe(
    session, recipe, buildRunInfo({ mode: recipe.mode, sourceRecipe: recipe }), 'identical features'
  );
});

// Design Q4 (PR-Q4-3, Owner Q2 = A): per-feature edits by source identity travel
// as --feature_override_table rows, so the Source recipe draws the Web Result.
const IDENTITY_GENBANK = [
  'LOCUS       dup                      120 bp    DNA     linear   UNK 01-JAN-1980',
  'DEFINITION  Identical CDS fixture.',
  'FEATURES             Location/Qualifiers',
  '     CDS             21..80',
  '                     /locus_tag="first"',
  '     CDS             21..80',
  '                     /locus_tag="second"',
  '     CDS             91..117',
  '                     /locus_tag="third"',
  'ORIGIN',
  '        1 atgatgatga tgatgatgat gatgatgatg atgatgatga tgatgatgat gatgatgatg',
  '       61 atgatgatga tgatgatgat gatgatgatg atgatgatga tgatgatgat gatgatgatg',
  '//',
  ''
].join('\n');
const identityEditSession = (mode, featureOverrides) => {
  const session = canonical({
    mode,
    records: [{
      recordKey: 'dup',
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'source' },
      selector: null,
      windows: [],
      presentation: presentation(),
      display: { isCircular: null, startCoordinate: null }
    }],
    resources: { source: resource('genbank', 'source-dup.gbk', IDENTITY_GENBANK) },
    webFiles: { resourceOriginalNames: { source: 'dup.gbk' } }
  });
  session.renderRequest.diagramOptions.selectedFeaturesSet = ['CDS'];
  Object.assign(session.renderRequest.diagramOptions, { featurePlacements: [], featureOverrides });
  return session;
};
const edit = (biologicalFeatureId, values) => ({
  recordKey: 'dup', biologicalFeatureId, featureVisibility: null, labelVisibility: null, labelText: null,
  ...values
});

// Asserts that the Web renderer on the session's request draws `cliSvg`, and
// that with `diagramOptions[field]` set to `plainValue` it draws something else.
// The CLI keys its records `record-<n>`, so the record key attribute differs.
const webRenderMatches = async (session, cliSvg, field, plainValue) => {
  const path = join(tempDir, `web-oracle-${session.renderRequest.mode}-${field}.json`);
  await writeFile(path, JSON.stringify({ session, field, plainValue }), 'utf8');
  const compared = spawnSync(process.env.PYTHON || 'python', ['-c', [
    'import copy, json, sys, tempfile',
    'from pathlib import Path',
    'from gbdraw.web_support.request_render import render_embedded_canonical_web_request',
    'from tests.utils.svg_compare import compare_svgs',
    'from tests.utils.web_render_result import web_result_svg',
    'data = json.loads(Path(sys.argv[1]).read_text())',
    'session = data["session"]',
    'plain = copy.deepcopy(session["renderRequest"])',
    'plain["diagramOptions"][data["field"]] = data["plainValue"]',
    'outcome = {}',
    'with tempfile.TemporaryDirectory() as directory:',
    '    for name, request in (("edited", session["renderRequest"]), ("plain", plain)):',
    '        result = render_embedded_canonical_web_request(',
    '            request, resources=session["resources"], workspace=Path(directory) / name)',
    '        web = Path(directory) / f"{name}.svg"',
    '        web.write_text(web_result_svg(result["results"][0]), encoding="utf-8")',
    '        compared = compare_svgs(web, Path(sys.argv[2]), ignored_attributes={"data-record-key"})',
    '        outcome[name] = bool(compared.equal)',
    '        if name == "edited": differences = compared.differences',
    'print(json.dumps({"outcome": outcome, "differences": differences}))'
  ].join('\n'), path, cliSvg], {
    cwd: repoRoot,
    encoding: 'utf8',
    env: { ...process.env, PYTHONPATH: [repoRoot, process.env.PYTHONPATH].filter(Boolean).join(delimiter) }
  });
  assert.equal(compared.status, 0, compared.stderr);
  const { outcome, differences } = JSON.parse(compared.stdout.trim().split('\n').at(-1));
  assert.deepEqual(outcome, { edited: true, plain: false }, `${field}: ${differences.join('\n')}`);
};

// OV-3009: a comparison ring row of a custom stack carries the Web's bookkeeping (managed, series_key, label, color,
// fileName), which Python takes and does not read (gbdraw/tracks/circular.py CIRCULAR_SLOT_UNREAD_PARAM_KEYS). The
// Source recipe leaves it out, and the CLI draws the Web Result.
test('a custom stack with a managed comparison ring row gives a Source recipe that the CLI replays', async () => {
  const testa = `${(await readFile(join(repoRoot, 'tests', 'fixtures', 'web_batch_two_records.gb'), 'utf8'))
    .split(/^\/\/\s*$/m)[0].trim()}\n//\n`;
  const blast = ['TESTB\tTESTA\t100.0\t600\t0\t0\t301\t900\t301\t900\t1e-80\t900',
    'TESTB\tTESTA\t85.0\t400\t20\t0\t2001\t2400\t2001\t2400\t1e-60\t500'].join('\n') + '\n';
  const session = canonical({
    mode: 'circular',
    records: [{
      recordKey: 'testa', cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'source' },
      selector: null, windows: [], presentation: presentation()
    }],
    resources: { source: resource('genbank', 'source-testa.gbk', testa), ring: resource('conservation-blast-file', 'ring-b.tsv', blast) },
    webFiles: { resourceOriginalNames: { source: 'testa.gbk', ring: 'ring-b.tsv' } }
  });
  const options = session.renderRequest.diagramOptions;
  options.conservationBlastFiles = [{ resourceId: 'ring', representation: 'file' }];
  options.conservationLabels = ['Ring B'];
  options.conservationColors = [''];
  const row = (id, renderer, side, params) => ({
    kind: 'circularTrackSlot', id, renderer, enabled: true, side, radius: null, width: null, z: 0, params,
    innerGapPx: null, outerGapPx: null
  });
  const plainTracks = structuredClone(options.tracks);
  options.tracks.circularTrackSlots = [
    row('features', 'features', 'overlay', { lane_direction: 'split' }),
    row('conservation_1', 'sequence_conservation', 'inside', {
      managed: 'circular_conservation', series_key: 'ring-b.tsv:0', source_index: '0', track_index: '1',
      label: 'Ring B', color: '', fileName: 'ring-b.tsv'
    })
  ];
  options.tracks.circularTrackAxisIndex = 1;
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const tokens = recipe.args.filter((_, index) => recipe.args[index - 1] === '--circular_track_slot');
  const ring = tokens.find((token) => token.startsWith('conservation_1:'));
  assert.match(ring, /^conservation_1:sequence_conservation@/, JSON.stringify(tokens));
  assert.doesNotMatch(ring, /managed|series_key|fileName|label=|color=/, ring);
  await assertCliParserAccepts(recipe);
  const cliSvg = await materializeAndRunRecipe(
    session, recipe, buildRunInfo({ mode: recipe.mode, sourceRecipe: recipe }), 'managed ring row'
  );
  await webRenderMatches(session, cliSvg, 'tracks', plainTracks);
});

test('feature edit recipe writes --feature_override_table and the CLI draws the Web Result', async () => {
  const [first, second, third] = await sourceFeatureIds(IDENTITY_GENBANK);
  assert.equal(first.replace(/~0$/, '~1'), second);
  for (const mode of ['linear', 'circular']) {
    const session = identityEditSession(mode, [
      edit(first, { featureVisibility: 'off' }),
      edit(second, { labelVisibility: 'on', labelText: ' Second "copy" ' }),
      edit(third, { featureVisibility: 'exclude_matching', labelVisibility: 'off', labelText: 'Third' })
    ]);
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    const tableIndex = recipe.args.indexOf('--feature_override_table');
    assert.notEqual(tableIndex, -1, JSON.stringify(recipe.args));
    const table = recipe.generatedFiles.find(
      ({ slot }) => slot === 'generatedFiles.source_recipe.feature_overrides'
    );
    assert.equal(table.path, recipe.args[tableIndex + 1]);
    assert.equal(table.data, [
      'record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text',
      `#1\thash=${first}\toff\t\t`,
      `#1\thash=${second}\t\ton\t" Second ""copy"" "`,
      `#1\thash=${third}\texclude_matching\toff\tThird`,
      ''
    ].join('\n'));
    const cliSvg = await materializeAndRunRecipe(
      session, recipe, buildRunInfo({ mode, sourceRecipe: recipe }), `${mode} feature edits`
    );
    // The oracle is the Web renderer on the same request and resources; without
    // the edits it must differ, so the comparison sees edits that reach the drawing.
    await webRenderMatches(session, cliSvg, 'featureOverrides', []);
  }
});

// PD-OI-069 revision 2: a "This feature only" fill is a `hash` row naming the
// feature instance `<feature id>@<record key>`. The recipe names the record by
// its CLI record key `record-<n>`, so the CLI colors the same one copy.
const twoCardColorSession = (rows) => {
  const card = (recordKey) => ({
    recordKey,
    cardinality: 'exactly_one',
    source: { kind: 'genbank', resourceId: 'source' },
    selector: null,
    windows: [],
    presentation: presentation(),
    display: { isCircular: null, startCoordinate: null }
  });
  const session = canonical({
    mode: 'linear',
    records: [card('linear-seq-a'), card('linear-seq-b')],
    resources: {
      source: resource('genbank', 'source-dup.gbk', IDENTITY_GENBANK),
      colors: resource('colors-color-table-file', 'specific-colors.tsv', rows.map((row) => `${row.join('\t')}\n`).join(''))
    },
    webFiles: { resourceOriginalNames: { source: 'dup.gbk' } }
  });
  session.renderRequest.diagramOptions.selectedFeaturesSet = ['CDS'];
  session.renderRequest.diagramOptions.colors.colorTableFile = { resourceId: 'colors', representation: 'file' };
  return session;
};

test('a feature instance color row names its record by the CLI record key and the CLI draws the Web Result', async () => {
  const [, , third] = await sourceFeatureIds(IDENTITY_GENBANK);
  const session = twoCardColorSession([
    ['CDS', 'hash', `^${third}@linear-seq-b$`, '#ff0000', 'Third of B'],
    ['CDS', 'product', 'kept', '#00ff00', 'Kept']
  ]);
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.specific_colors');
  assert.ok(table, JSON.stringify(recipe.generatedFiles.map(({ slot }) => slot)));
  assert.equal(recipe.args[recipe.args.indexOf('-t') + 1], table.path);
  assert.equal(recipe.args.filter((token) => token === '-t').length, 1);
  assert.equal(table.data, [
    `CDS\thash\t^${third}@record-2$\t#ff0000\tThird of B`,
    'CDS\tproduct\tkept\t#00ff00\tKept',
    ''
  ].join('\n'));
  const cliSvg = await materializeAndRunRecipe(
    session, recipe, buildRunInfo({ mode: 'linear', sourceRecipe: recipe }), 'feature instance color'
  );
  const svg = await readFile(cliSvg, 'utf8');
  // One copy only: the CLI fills the third feature of the second record and not its twin in the first.
  const filled = new Set([...svg.matchAll(/<path [^>]*data-gbdraw-feature-id="([^"]+)"[^>]*fill="#ff0000"/gi)]
    .map((match) => match[1]));
  assert.deepEqual([...filled], [`${third}_record_2`], 'the CLI draws the instance color on one copy');
  const colors = session.renderRequest.diagramOptions.colors;
  await webRenderMatches(session, cliSvg, 'colors', { ...colors, colorTableFile: null });
});

// A card that draws every record of its file keys them `<card key>:<k>`, as the
// CLI keys its record input `record-<n>:<k>` (gbdraw/api/record_planning.py).
// The CLI run pins the composed key against Python's matcher.
test('a feature instance color row of an expanded record names it by the CLI record key and its record number', async () => {
  const [, , third] = await sourceFeatureIds(IDENTITY_GENBANK);
  const session = twoCardColorSession([['CDS', 'hash', `^${third}@linear-seq-a:2$`, '#ff0000', 'Third of 2']]);
  session.renderRequest.records = [{ ...session.renderRequest.records[0], cardinality: 'all' }];
  Object.assign(session.resources.source, resource('genbank', 'source-dup.gbk', `${IDENTITY_GENBANK}${IDENTITY_GENBANK}`));
  const recipe = await buildSourceRecipe(session);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.specific_colors');
  assert.equal(table?.data, `CDS\thash\t^${third}@record-1:2$\t#ff0000\tThird of 2\n`);
  const cliSvg = await materializeAndRunRecipe(
    session, recipe, buildRunInfo({ mode: 'linear', sourceRecipe: recipe }), 'expanded feature instance color'
  );
  const svg = await readFile(cliSvg, 'utf8');
  const filled = [...svg.matchAll(/<path [^>]*data-gbdraw-feature-id="([^"]+)"[^>]*fill="#ff0000"/gi)].map((match) => match[1]);
  assert.equal(filled.length, 1, `the CLI draws the instance color on one copy: ${filled}`);
  const colors = session.renderRequest.diagramOptions.colors;
  await webRenderMatches(session, cliSvg, 'colors', { ...colors, colorTableFile: null });
});

test('a feature instance color row of a table whose bytes a File holds is read through readResourceText', async () => {
  const [first] = await sourceFeatureIds(IDENTITY_GENBANK);
  const text = `CDS\thash\t^${first}@linear-seq-a$\t#ff0000\tFirst of A\n`;
  const session = twoCardColorSession([]);
  // The descriptor holds no bytes, as for a file the user chose (its File holds them).
  delete session.resources.colors.data;
  const recipe = await buildSourceRecipe({
    ...session, readResourceText: (resourceId) => (resourceId === 'colors' ? text : '')
  });
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.specific_colors');
  assert.equal(table?.data, `CDS\thash\t^${first}@record-1$\t#ff0000\tFirst of A\n`);
});

test('a table without feature instance rows is passed as it is', async () => {
  const [first] = await sourceFeatureIds(IDENTITY_GENBANK);
  const recipe = await buildSourceRecipe(twoCardColorSession([['CDS', 'hash', first, '#ff0000', 'Shared']]));
  assert.equal(recipe.available, true, recipe.unavailableReason);
  assert.equal(recipe.generatedFiles.some(({ slot }) => slot === 'generatedFiles.source_recipe.specific_colors'), false);
  assert.notEqual(recipe.args.indexOf('-t'), -1);
});

// A row whose record this render does not draw (a removed card, a changed
// record key) colors nothing in it, so the recipe's table leaves it out.
test('a feature instance color row naming a record this render does not draw stays out of the recipe table', async () => {
  const [first] = await sourceFeatureIds(IDENTITY_GENBANK);
  const gone = ['CDS', 'hash', `^${first}@linear-seq-gone$`, '#ff0000', 'Gone'];
  let recipe = await buildSourceRecipe(twoCardColorSession([gone, ['CDS', 'product', 'kept', '#00ff00', 'Kept']]));
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.specific_colors');
  assert.equal(table?.data, 'CDS\tproduct\tkept\t#00ff00\tKept\n');
  // With no row left, the recipe names no table.
  recipe = await buildSourceRecipe(twoCardColorSession([gone]));
  assert.equal(recipe.available, true, recipe.unavailableReason);
  assert.equal(recipe.args.includes('-t'), false, JSON.stringify(recipe.args));
});

// Python matches an instance row whole, written bare or as an exact pattern.
test('a bare feature instance color row names its record by the CLI record key too', async () => {
  const [, second] = await sourceFeatureIds(IDENTITY_GENBANK);
  const recipe = await buildSourceRecipe(twoCardColorSession([['CDS', 'hash', `${second}@LINEAR-seq-b`, '#ff0000', 'B']]));
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.specific_colors');
  assert.equal(table?.data, `CDS\thash\t^${second}@record-2$\t#ff0000\tB\n`);
});

test('a hash pattern that names feature instances by the Web record keys makes the Source recipe unavailable', async () => {
  const [first] = await sourceFeatureIds(IDENTITY_GENBANK);
  const recipe = await buildSourceRecipe(twoCardColorSession([['CDS', 'hash', `^${first}@linear-seq-.$`, '#ff0000', 'Both']]));
  assert.equal(recipe.available, false);
  assert.match(recipe.unavailableReason, /names feature instances by a pattern/);
});

test('feature edit rows that the table cannot read back make the Source recipe unavailable', async () => {
  const [first] = await sourceFeatureIds(IDENTITY_GENBANK);
  for (const row of [
    edit(first, { labelText: 'tab\tinside' }),
    edit(first, { featureVisibility: 'hidden' }),
    edit(first, {}),
    { ...edit(first, { labelVisibility: 'on' }), recordKey: 'other' }
  ]) {
    const recipe = await buildSourceRecipe(identityEditSession('linear', [row]));
    assert.equal(recipe.available, false, JSON.stringify(row));
    assert.match(recipe.unavailableReason, /feature edit/i, JSON.stringify(row));
  }
});

test('a feature edit or placement whose feature the source lacks makes the Source recipe unavailable', async () => {
  const [first, second] = await sourceFeatureIds(IDENTITY_GENBANK);
  const session = identityEditSession('linear', [edit(first, { featureVisibility: 'off' })]);
  session.renderRequest.diagramOptions.featurePlacements = [
    { recordKey: 'dup', biologicalFeatureId: second, placement: { kind: 'main' } }
  ];
  for (const [kinds, expected] of [
    [['feature_visibility'], false], [['placement'], false], [null, true]
  ]) {
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    const featureIdentityNotices = kinds ? [{
      recordKey: 'dup', biologicalFeatureId: 'fdeadbeef', status: 'unresolved', kinds, resultIndex: 0
    }, {
      recordKey: 'dup', biologicalFeatureId: first, status: 'crop_excluded', kinds: ['label_text'], resultIndex: 0
    }] : [{
      recordKey: 'dup', biologicalFeatureId: first, status: 'crop_excluded', kinds: ['label_text'], resultIndex: 0
    }];
    const info = buildRunInfo({ mode: 'linear', sourceRecipe: recipe, featureIdentityNotices });
    assert.equal(info.sourceRecipe.available, expected, JSON.stringify(kinds));
    if (!expected) {
      assert.match(info.sourceRecipe.unavailableReason, /source does not have/);
      assert.equal(info.sourceRecipe.command, '');
    }
  }
});

// PD-OI-018 revision 4 and D-40, OV-3143: Run Info states the LOSAT E-value scale.
test('Run Info states the LOSAT E-value scale only when LOSAT comparisons are present', () => {
  const recipe = { mode: 'linear', available: true, args: ['--gbk', '/a.gb'], fileMetadata: [] };
  const withLosat = buildRunInfo({ mode: 'linear', sourceRecipe: recipe, losatComparisons: true });
  assert(withLosat.reproducibility.notes.includes(LOSAT_EVALUE_SCALE_NOTE));
  assert.match(LOSAT_EVALUE_SCALE_NOTE, /LOSATN and TLOSATX E-values use one scale for each pair of genomes/);
  assert.match(LOSAT_EVALUE_SCALE_NOTE, /LOSATP E-values use the proteins of the subject genome as the database size/);
  assert.doesNotMatch(LOSAT_EVALUE_SCALE_NOTE, /subject source file/);
  // The CLI and Python API use the same scale (design D7).
  assert.doesNotMatch(LOSAT_EVALUE_SCALE_NOTE, /CLI/);
  const without = buildRunInfo({ mode: 'linear', sourceRecipe: recipe });
  assert(!without.reproducibility.notes.includes(LOSAT_EVALUE_SCALE_NOTE));
  // A Result made before stored the database note, which stays true of it and is still shown.
  assert.deepEqual(shownRunInfoNoteParts([{ text: LOSAT_DATABASE_SCOPE_NOTE }]), [{ text: LOSAT_DATABASE_SCOPE_NOTE }]);
});

// PD-OI-018 revision 6 (OV-3319): a Circular Result with LOSAT rings states the rings' scale, the reference genome
// against each ring's genome, in the words of the Linear note.
test('Run Info states the similarity ring E-value scale of a Circular Result with LOSAT rings', () => {
  const recipe = { mode: 'circular', available: true, args: ['--gbk', '/a.gb'], fileMetadata: [] };
  const withRings = buildRunInfo({ mode: 'circular', sourceRecipe: recipe, losatComparisons: true });
  assert(withRings.reproducibility.notes.includes(LOSAT_RING_EVALUE_SCALE_NOTE));
  assert(!withRings.reproducibility.notes.includes(LOSAT_EVALUE_SCALE_NOTE));
  assert.match(LOSAT_RING_EVALUE_SCALE_NOTE, /the reference genome is the files of the drawn records and the ring's genome is the files of that ring/);
  assert.match(LOSAT_RING_EVALUE_SCALE_NOTE, /with every record of those files counted, drawn or not/);
  assert.doesNotMatch(LOSAT_RING_EVALUE_SCALE_NOTE, /CLI/);
  assert(!buildRunInfo({ mode: 'circular', sourceRecipe: recipe }).reproducibility.notes.includes(LOSAT_RING_EVALUE_SCALE_NOTE));
  assert.deepEqual(shownRunInfoNoteParts([{ text: LOSAT_RING_EVALUE_SCALE_NOTE }]), [{ text: LOSAT_RING_EVALUE_SCALE_NOTE }]);
});

test('Run Info lists the search runtime of Web, CLI and unrecorded raw entries (D9/D10)', () => {
  const cacheMap = new Map([
    ['web', { program: 'blastn', runtime: { kind: 'losat', source: 'wasm', version: null, program: 'blastn' } }],
    ['web-copy', { program: 'blastn', runtime: { kind: 'losat', source: 'wasm', version: null, program: 'blastn' } }],
    ['cli', {
      program: 'tblastx',
      runtime: {
        kind: 'losat', version: '0.1.0', source: 'bundled',
        path: 'gbdraw/bin/linux-x86_64/losat', program: 'tblastx', cli: 'v1'
      }
    }],
    ['ncbi', { program: 'blastp', runtime: { kind: 'ncbi-blast', version: '2.16.0+', source: 'path', path: '/usr/bin/blastp', program: 'blastp' } }],
    ['old', { program: 'blastn' }],
    ['hidden', { program: 'blastp', runtime: { kind: 'losat', source: 'wasm', version: null, program: 'blastp' } }]
  ]);
  const info = ['web', 'web-copy', 'cli', 'ncbi', 'old'].map((key) => ({ key, display: true }))
    .concat([{ key: 'hidden', display: false }, { key: 'missing', display: true }]);
  assert.deepEqual(summarizeLosatRuntimes(info, cacheMap).map(({ text }) => text), [
    'blastn: LOSAT, version not recorded (wasm)',
    'tblastx: LOSAT 0.1.0 (bundled, gbdraw/bin/linux-x86_64/losat)',
    'blastp: NCBI BLAST+ 2.16.0+ (path, /usr/bin/blastp)',
    'blastn: runtime not recorded'
  ]);
  const runInfo = buildRunInfo({
    mode: 'linear',
    sourceRecipe: { available: true, args: ['--gbk', 'a.gb'], generatedFiles: [] },
    losatRuntimes: summarizeLosatRuntimes(info.slice(2, 3), cacheMap)
  });
  assert.deepEqual(runInfo.losatRuntimes, [{
    program: 'tblastx', kind: 'losat', version: '0.1.0', source: 'bundled',
    path: 'gbdraw/bin/linux-x86_64/losat',
    text: 'tblastx: LOSAT 0.1.0 (bundled, gbdraw/bin/linux-x86_64/losat)',
    textParts: [{ text: 'tblastx: LOSAT 0.1.0 (bundled, gbdraw/bin/linux-x86_64/losat)', verbatim: true }]
  }]);
  assert.deepEqual(buildRunInfo({ mode: 'circular' }).losatRuntimes, []);
});

// OV-1678: the Run info keeps its texts in English (the same in every UI
// language) beside the parts they are made of, and the Run info tab shows the
// parts in the UI language.
test('OV-1678: in Simplified Chinese the Run info shows its notes, files description, unavailable reason and runtimes translated', async () => {
  const { setUiLocale, SOURCE_LOCALE } = await import(pathToFileURL(join(tempDir, 'js', 'utils', 'i18n.js')));
  const { errorText } = await import(pathToFileURL(join(tempDir, 'js', 'utils', 'error-text.js')));
  const { default: zh } = await import(pathToFileURL(join(tempDir, 'js', 'utils', 'i18n-catalogs', 'zh-Hans.js')));
  const cacheMap = new Map([['old', { program: 'blastn' }]]);
  const build = () => buildRunInfo({
    mode: 'linear',
    // Two uploaded inputs with one visible name: the reason names it.
    sourceRecipe: { mode: 'linear', available: true, args: ['--gbk', '/a/x.gb', '/b/x.gb'], generatedFiles: [] },
    fileMetadata: { '/a/x.gb': { name: 'x.gb' }, '/b/x.gb': { name: 'x.gb' } },
    exactReplayArgs: ['--session', '/canonical-render-session.gbdraw-session.json'],
    startedAtIso: '2026-10-23T00:00:00.000Z',
    losatComparisons: true,
    losatRuntimes: summarizeLosatRuntimes([{ key: 'old', display: true }], cacheMap)
  });
  const english = build();
  const reason = 'Source recipe unavailable: multiple uploaded inputs use the visible filename "x.gb".';
  assert.equal(english.sourceRecipe.unavailableReason, reason);
  assert.deepEqual(english.reproducibility.notes, [
    reason, 'Exact replay remains available from the committed canonical session.', LOSAT_EVALUE_SCALE_NOTE
  ]);
  assert.equal(english.downloadBundle.description, 'The Source recipe is unavailable. No Exact replay file is available.');
  assert.equal(english.losatRuntimes[0].text, 'blastn: runtime not recorded');
  // In English the parts show the stored text.
  assert.equal(errorText(english.sourceRecipe.unavailableReasonParts), reason);
  assert.deepEqual(english.reproducibility.noteParts.map((part) => errorText([part])), english.reproducibility.notes);
  await setUiLocale('zh-Hans');
  try {
    // The Run info stays English in every UI language.
    assert.deepEqual(build(), english);
    const shownReason = zh['Source recipe unavailable: {reason}'].replace('{reason}',
      zh['multiple uploaded inputs use the visible filename "{name}".'].replace('{name}', 'x.gb'));
    assert.equal(errorText(english.sourceRecipe.unavailableReasonParts), shownReason);
    assert.deepEqual(english.reproducibility.noteParts.map((part) => errorText([part])), [
      shownReason, zh['Exact replay remains available from the committed canonical session.'],
      // A new note reads in English until the translation batch before a release.
      zh[LOSAT_EVALUE_SCALE_NOTE] ?? LOSAT_EVALUE_SCALE_NOTE
    ]);
    assert.equal(errorText(english.downloadBundle.descriptionParts),
      `${zh['The Source recipe is unavailable.']}${zh['No Exact replay file is available.']}`);
    assert.equal(errorText(english.losatRuntimes[0].textParts), zh['{program}: runtime not recorded'].replace('{program}', 'blastn'));
  } finally {
    await setUiLocale(SOURCE_LOCALE);
  }
});

// ui-copy-diet (2026-11-02): the Run info keeps every note; the tab shows only the notes a reader acts on, and
// says what a Source recipe helper needs in one line.
test('the Run info tab shows the notes a reader acts on, from the stored notes', async () => {
  const { errorText } = await import(pathToFileURL(join(tempDir, 'js', 'utils', 'error-text.js')));
  const helper = (count) => buildRunInfo({
    mode: 'linear',
    args: ['--gbk', '/seq_0.gb', '/seq_1.gb', ...Array.from({ length: count }, (_, index) => ['-b', `/blast_${index}.txt`]).flat()],
    fileMetadata: {
      '/seq_0.gb': { name: 'alpha.gb', slot: 'files.linearSeqs[0].gb', kind: 'uploaded' },
      '/seq_1.gb': { name: 'beta.gb', slot: 'files.linearSeqs[1].gb', kind: 'uploaded' },
      ...Object.fromEntries(Array.from({ length: count }, (_, index) => [`/blast_${index}.txt`,
        { name: `pair_${index}.tsv`, slot: `generatedFiles.losat_blasts[${index}]`, kind: 'generated' }]))
    },
    losatComparisons: true
  }).reproducibility;
  for (const count of [1, 2]) {
    const stored = helper(count);
    assert.equal(stored.notes.length, 3, 'the stored notes keep the helper and Session restore notes');
    assert.deepEqual(shownRunInfoNoteParts(stored.noteParts).map((part) => errorText([part])), [
      count === 1
        ? 'Needs pair_0.tsv from "Download reproducibility files"'
        : 'Needs pair_0.tsv, pair_1.tsv from "Download reproducibility files"',
      LOSAT_EVALUE_SCALE_NOTE
    ]);
  }
  const plain = buildRunInfo({ mode: 'circular', args: ['--gbk', '/seq.gb'], fileMetadata: { '/seq.gb': { name: 'seq.gb', kind: 'uploaded' } } });
  assert.deepEqual(plain.reproducibility.notes, ['The Source recipe can be rerun with the uploaded file names shown here.']);
  assert.deepEqual(shownRunInfoNoteParts(plain.reproducibility.noteParts), []);
  const unread = { text: 'Some browser virtual paths could not be mapped to visible file names, so this command is only an approximation.' };
  assert.deepEqual(shownRunInfoNoteParts([{ text: 'Source recipe unavailable: {reason}', params: { reason: { text: 'x' } } }, unread]), [unread]);
});

test('the download description counts its helper files', () => {
  const exactFile = (name) => ({ name, slot: 'generatedFiles.canonical_render_session', kind: 'generated' });
  const description = (paths) => buildRunInfo({ mode: 'circular', exactReplayArgs: ['--session', ...paths],
    fileMetadata: Object.fromEntries(paths.map((path) => [path, exactFile(`${path.slice(1)}.json`)])) }).downloadBundle.descriptionParts.at(-1);
  assert.equal(description(['/one']).text, 'Exact replay file: {files}.');
  assert.equal(description(['/one', '/two']).text, 'Exact replay files: {files}.');
});

// R12, issue #1427: the Legend row edits and the canvas padding of a request are written as `--legend_rows` (a
// generated table) and `--canvas_padding`, and read back through the CLI's own readers to the same leaves; a hidden
// row with an edited order of the others needs OV-2514. A leaf the table cannot carry makes the recipe unavailable.
const legendLeavesReadBack = async (recipe) => {
  const runDir = await mkdtemp(join(tmpdir(), 'gbdraw-legend-rows-'));
  const args = [];
  for (const token of recipe.args) {
    const generated = recipe.generatedFiles.find((file) => file.path === token);
    if (generated) {
      const file = join(runDir, generated.name);
      await writeFile(file, generated.data, 'utf8');
      args.push(file);
    } else args.push(recipe.fileMetadata.get(token)?.name || token);
  }
  return evaluatePython([
    'from gbdraw.circular import _get_args',
    'from gbdraw.cli_utils.common import legend_row_overrides, figure_size_overrides',
    'args = _get_args(payload["args"])',
    'leaves = {**legend_row_overrides(args), "canvas.padding": figure_size_overrides(args)["canvas.padding"]}',
    'result = {path: value for path, value in leaves.items() if value is not None}'
  ].join('\n'), { args });
};

test('Legend row edits and canvas padding read back from the Source recipe (R12, OV-2514)', async () => {
  const leaves = {
    'objects.legends.row_names': [['rRNA', 'Ribosomal RNA'], ['CDS', 'Genes "all" (5\')']],
    'objects.legends.hidden_rows': ['tRNA'],
    'objects.legends.row_order': ['rRNA', 'CDS'],
    'canvas.padding': [10, 20.5, 0, 40]
  };
  for (const set of [leaves, { 'objects.legends.hidden_rows': ['tRNA', 'CDS'] }, { 'canvas.padding': [0, 8, 0, 0] }]) {
    const session = structuredClone(circularCanonical);
    session.renderRequest.diagramOptions.configOverrides = { ...(session.renderRequest.diagramOptions.configOverrides || {}), ...set };
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    assert.equal(recipe.args.includes('--legend_rows'), Object.keys(set).some((path) => path.startsWith('objects.legends.')));
    assert.deepEqual(await legendLeavesReadBack(recipe), set, JSON.stringify(recipe.generatedFiles.map((file) => file.data)));
  }
  const none = await buildSourceRecipe(circularCanonical);
  assert.equal(none.args.includes('--legend_rows') || none.args.includes('--canvas_padding'), false);
  for (const set of [
    // A cell the table cannot carry: a tab or a line break, a key with outer spaces, a blank name.
    { 'objects.legends.row_names': [['CDS', 'Genes\tall']] },
    { 'objects.legends.row_names': [['CDS', 'Genes\nall']] },
    { 'objects.legends.hidden_rows': [' tRNA'] },
    { 'objects.legends.row_names': [['CDS', ' ']] },
    // One key on two lines, and a shown row the edited order does not place.
    { 'objects.legends.row_names': [['CDS', 'A'], ['CDS', 'B']] },
    { 'objects.legends.row_order': ['CDS', 'CDS'] },
    { 'objects.legends.row_order': ['rRNA'], 'objects.legends.row_names': [['CDS', 'Genes']] },
    // Names or hidden rows in an order the table's lines cannot give back.
    { 'objects.legends.row_order': ['rRNA', 'CDS'], 'objects.legends.row_names': [['CDS', 'A'], ['rRNA', 'B']] },
    { 'canvas.padding': [1, 2, 3] },
    { 'canvas.padding': [1, 2, -3, 4] }
  ]) {
    const session = structuredClone(circularCanonical);
    session.renderRequest.diagramOptions.configOverrides = { ...(session.renderRequest.diagramOptions.configOverrides || {}), ...set };
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, false, JSON.stringify(set));
  }
});

// PD-OI-122: a window across the origin is written `START-END:origin` (before `:rc`), on the command line and in the
// records table, so the CLI draws the window the Web drew; a records table writes every window of a record, its rows
// sharing a record_group value. Python's region parser reads each written region back as the request's window.
test('the Source recipe writes a window across the origin with :origin and every window of a record', async () => {
  const genbank = await readFile(join(repoRoot, 'tests', 'fixtures', 'web_comparison_shared_block.gb'), 'utf8');
  const windows = [
    { start: 2501, end: 800, reverseComplement: true, column: null },
    { start: 1001, end: 1600, reverseComplement: false, column: null }
  ];
  const parsedBack = (regions) => evaluatePython([
    'from gbdraw.io.regions import parse_region_spec',
    'result = [[spec.start, spec.end, spec.reverse_complement] for spec in map(parse_region_spec, payload)]'
  ].join('\n'), regions);
  const expected = windows.map(({ start, end, reverseComplement }) => [start, end, reverseComplement]);
  // One window across the origin goes on the command line; a record drawn in several windows needs a records table.
  for (const [isCircular, drawn] of [[null, windows.slice(0, 1)], [true, windows.slice(0, 1)], [null, windows]]) {
    const session = canonical({
      mode: 'linear',
      records: [{
        recordKey: 'r3c',
        cardinality: 'exactly_one',
        source: { kind: 'genbank', resourceId: 'source' },
        selector: { kind: 'recordId', value: 'R3c' },
        windows: drawn,
        presentation: presentation(),
        // A topology the Web sets is written in a records table.
        display: { isCircular, startCoordinate: null }
      }],
      resources: { source: resource('genbank', 'shared.gb', genbank) },
      webFiles: { resourceOriginalNames: { source: 'shared.gb' } }
    });
    const recipe = await buildSourceRecipe(session);
    assert.equal(recipe.available, true, recipe.unavailableReason);
    const table = recipe.generatedFiles.find((file) => file.name === 'records.tsv');
    let regions;
    if (!table) {
      assert.equal(isCircular, null);
      assert.equal(drawn.length, 1);
      regions = recipe.args.flatMap((token, index) => (token === '--region' ? [recipe.args[index + 1]] : []));
      assert.deepEqual(regions, ['R3c:2501-800:origin:rc']);
    } else {
      const [header, ...lines] = table.data.trim().split('\n').map((line) => line.split('\t'));
      const column = (name) => lines.map((line) => line[header.indexOf(name)]);
      regions = column('region');
      assert.deepEqual(regions, ['2501-800:origin:rc', '1001-1600'].slice(0, drawn.length));
      assert.deepEqual(column('record_group'), drawn.length > 1 ? ['1', '1'] : [undefined]);
      assert.deepEqual(column('topology'), drawn.map(() => (isCircular ? 'circular' : '')));
    }
    assert.deepEqual(await parsedBack(regions), expected.slice(0, drawn.length));
    await assertCliParserAccepts(recipe);
    const info = buildRunInfo({ mode: recipe.mode, sourceRecipe: recipe });
    await materializeAndRunRecipe(session, recipe, info, `origin window ${drawn.length} ${isCircular}`);
  }
});

// Review S4 of PR #1734: a request's comparison indexes count drawn units (Python's unit plan, which Generate passes
// as `drawnUnits`), while the CLI's `#n` counts records. A records table writes record A's two windows as one record,
// so a comparison of B with C is `#2→#3`; a comparison of a record drawn in several windows has no comparison file the
// CLI takes (it refuses one, K-2), so the recipe is unavailable and says why.
test('the Source recipe names compared records by record, not by drawn unit, and refuses a windowed record\'s file', async () => {
  const genbank = await readFile(join(repoRoot, 'tests', 'fixtures', 'web_comparison_shared_block.gb'), 'utf8');
  const record = (recordKey, id, windows = []) => ({
    recordKey,
    cardinality: 'exactly_one',
    source: { kind: 'genbank', resourceId: `source-${recordKey}` },
    selector: { kind: 'recordId', value: id },
    windows,
    presentation: presentation()
  });
  const twoWindows = [
    { start: 1, end: 1000, reverseComplement: false, column: null },
    { start: 2001, end: 3000, reverseComplement: false, column: null }
  ];
  const drawnUnits = { recordKeys: ['a', 'a', 'b', 'c'], windowOf: ['a', 'a', null, null] };
  const session = (comparisons) => canonical({
    mode: 'linear',
    records: [record('a', 'R2c', twoWindows), record('b', 'R3c'), record('c', 'R2c')],
    resources: {
      'source-a': resource('genbank', 'a.gb', genbank),
      'source-b': resource('genbank', 'b.gb', genbank),
      'source-c': resource('genbank', 'c.gb', genbank),
      blast: resource('nucleotide-blast', 'blast.tsv', 'R3c\tR2c\t100\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n')
    },
    comparisons: comparisons.map(([queryRecordIndex, subjectRecordIndex]) => ({
      kind: 'nucleotideBlast', resourceId: 'blast', queryRecordIndex, subjectRecordIndex
    })),
    webFiles: { resourceOriginalNames: { 'source-a': 'a.gb', 'source-b': 'b.gb', 'source-c': 'c.gb', blast: 'blast.tsv' } }
  });

  // A compared with B: the unit pairs (0, 2) and (1, 2), which are A-C and B-C read as records.
  const windowed = await buildSourceRecipe({ ...session([[0, 2], [1, 2]]), drawnUnits });
  assert.equal(windowed.available, false);
  assert.match(windowed.unavailableReason, /several windows/);

  // B compared with C: unit pair (2, 3) is records #2 and #3.
  const plain = session([[2, 3]]);
  const recipe = await buildSourceRecipe({ ...plain, drawnUnits });
  assert.equal(recipe.available, true, recipe.unavailableReason);
  // Run Info names the generated files and writes the comparisons table with those names, in the recipe's list.
  const info = buildRunInfo({ mode: recipe.mode, sourceRecipe: recipe });
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.comparisons_table');
  const [header, ...rows] = table.data.trim().split('\n').map((line) => line.split('\t'));
  assert.deepEqual(rows.map((row) => [row[header.indexOf('query')], row[header.indexOf('subject')]]), [['#2', '#3']]);
  await assertCliParserAccepts(recipe);
  await materializeAndRunRecipe(plain, recipe, info, 'units by record');

  // Without the unit plan a request that draws a record in several windows cannot be read as records.
  assert.equal((await buildSourceRecipe(plain)).available, false);
});

// OV-3072: a Linear File card that draws a file of several records whole, placed in a row beside a record drawn in two
// windows, is written as a records table row per record (#1, #2) at the columns Python numbers (the Worker helper
// `planRecordPlaces`, here the real `plan_record_places_json`), and the CLI draws the Web Result.
/**
 * A Linear request with a placed File card `m` of two records at row 1, column 1, a record `w` in two windows
 * (columns 2 and 3) in row 1, and a record `c` at row 2, column 1; the Source recipe's place plan is the real Python's.
 * @param {{ colorRows?: string[][] | null, comparisons?: [number, number][] }} [extra]
 */
const placedFileCardSession = async ({ colorRows = null, comparisons = [] } = {}) => {
  const text = await readFile(join(repoRoot, 'tests', 'fixtures', 'web_batch_two_records.gb'), 'utf8');
  const [testa, testb] = text.split(/^\/\/\s*$/m).slice(0, 2).map((part) => `${part.trim()}\n//\n`);
  const at = (gridRow, gridColumn) => ({ ...presentation(), gridRow, gridColumn });
  const resources = {
    both: resource('genbank', 'both.gbk', text), one: resource('genbank', 'one.gbk', testa),
    other: resource('genbank', 'other.gbk', testb)
  };
  if (colorRows) {
    resources.colors = resource('colors-color-table-file', 'specific-colors.tsv', colorRows.map((row) => `${row.join('\t')}\n`).join(''));
  }
  if (comparisons.length) {
    resources.blast = resource('nucleotide-blast', 'blast.tsv', 'TESTA\tTESTB\t100\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n');
  }
  const session = canonical({
    mode: 'linear',
    records: [
      { recordKey: 'm', cardinality: 'all', source: { kind: 'genbank', resourceId: 'both' }, selector: null,
        windows: [], presentation: at(1, 1) },
      { recordKey: 'w', cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'one' }, selector: null,
        windows: [
          { start: 1, end: 1000, reverseComplement: false, column: 2 },
          { start: 2001, end: 3000, reverseComplement: false, column: 3 }
        ], presentation: at(1, null) },
      { recordKey: 'c', cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'other' }, selector: null,
        windows: [], presentation: at(2, 1) }
    ],
    resources,
    comparisons: comparisons.map(([query, subject]) => ({
      kind: 'nucleotideBlast', resourceId: 'blast', queryRecordIndex: query, subjectRecordIndex: subject
    })),
    webFiles: { resourceOriginalNames: { both: 'both.gbk', one: 'one.gbk', other: 'other.gbk', blast: 'blast.tsv' } }
  });
  session.renderRequest.layout = { recordGapPx: 24, multiRecordPositions: null, recordTranslations: [], similarityAlignment: null };
  if (colorRows) {
    session.renderRequest.diagramOptions.selectedFeaturesSet = ['CDS'];
    session.renderRequest.diagramOptions.colors.colorTableFile = { resourceId: 'colors', representation: 'file' };
  }
  const counts = { both: 2, one: 1, other: 1 };
  /** @type {unknown[]} */
  const asked = [];
  const options = {
    ...session,
    readResourceRecordCount: (id) => counts[id],
    planRecordPlaces: async (payload) => {
      asked.push(payload);
      return { result: await evaluatePython([
        'import json',
        'from gbdraw.web_support.linear_units import plan_record_places_json',
        'result = json.loads(plan_record_places_json(json.dumps(payload["records"])))'
      ].join('\n'), payload) };
    }
  };
  return { session, testb, asked, options };
};

/**
 * The Web's render of `session` and its comparison with the CLI's `cliSvg` (record keys aside).
 * @returns {{ equal: boolean, differences: string[], webSvg: string }}
 */
const compareWithWebRender = async (session, cliSvg, tag) => {
  const path = join(tempDir, `web-oracle-${tag}.json`);
  const webPath = join(tempDir, `web-oracle-${tag}.svg`);
  await writeFile(path, JSON.stringify(session), 'utf8');
  const compared = spawnSync(process.env.PYTHON || 'python', ['-c', [
    'import json, sys, tempfile',
    'from pathlib import Path',
    'from gbdraw.web_support.request_render import render_embedded_canonical_web_request',
    'from tests.utils.svg_compare import compare_svgs',
    'from tests.utils.web_render_result import web_result_svg',
    'session = json.loads(Path(sys.argv[1]).read_text())',
    'with tempfile.TemporaryDirectory() as directory:',
    '    result = render_embedded_canonical_web_request(',
    '        session["renderRequest"], resources=session["resources"], workspace=Path(directory) / "web")',
    '    web = Path(sys.argv[3])',
    '    web.write_text(web_result_svg(result["results"][0]), encoding="utf-8")',
    '    compared = compare_svgs(web, Path(sys.argv[2]), ignored_attributes={"data-record-key"})',
    'print(json.dumps({"equal": bool(compared.equal), "differences": compared.differences}))'
  ].join('\n'), path, cliSvg, webPath], {
    cwd: repoRoot,
    encoding: 'utf8',
    env: { ...process.env, PYTHONPATH: [repoRoot, process.env.PYTHONPATH].filter(Boolean).join(delimiter) }
  });
  assert.equal(compared.status, 0, compared.stderr);
  return { ...JSON.parse(compared.stdout.trim().split('\n').at(-1)), webSvg: await readFile(webPath, 'utf8') };
};

test('a placed File card of several records gives a records table row per record that the CLI draws as the Web does', async () => {
  const { session, asked, options } = await placedFileCardSession();
  const recipe = await buildSourceRecipe(options);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  assert.deepEqual(asked, [{ records: [
    { gridRow: 1, gridColumn: 1, windowColumns: [], drawnRecords: 2 },
    { gridRow: 1, gridColumn: null, windowColumns: [2, 3], drawnRecords: 1 },
    { gridRow: 2, gridColumn: 1, windowColumns: [], drawnRecords: 1 }
  ] }]);
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.records_table');
  const rows = String(table.data).trim().split('\n').map((line) => line.split('\t'));
  const header = rows[0];
  const cells = (row) => Object.fromEntries(header.map((column, index) => [column, row[index]]));
  assert.deepEqual(rows.slice(1).map(cells).map(({ gbk, record_id: id, region, order, row, column }) => (
    [gbk, id, region, order, row, column])), [
    ['both.gbk', '#1', '', '1', '1', '1'],
    ['both.gbk', '#2', '', '2', '1', '2'],
    ['one.gbk', '', '1-1000', '3', '1', '3'],
    ['one.gbk', '', '2001-3000', '3', '1', '4'],
    ['other.gbk', '', '', '4', '2', '1']
  ]);
  await assertCliParserAccepts(recipe);
  const cliSvg = await materializeAndRunRecipe(
    session, recipe, buildRunInfo({ mode: 'linear', sourceRecipe: recipe }), 'placed file card of several records'
  );
  const { equal, differences } = await compareWithWebRender(session, cliSvg, 'placed-file-card');
  assert.equal(equal, true, differences.join('\n'));
});

test('a colour rule naming one record of a placed File card names that record\'s row in the Source recipe', async () => {
  const text = await readFile(join(repoRoot, 'tests', 'fixtures', 'web_batch_two_records.gb'), 'utf8');
  const [, testb] = text.split(/^\/\/\s*$/m).slice(0, 2).map((part) => `${part.trim()}\n//\n`);
  const [, second] = await sourceFeatureIds(testb);
  const { session, options } = await placedFileCardSession({ colorRows: [
    ['CDS', 'hash', `^${second}@m:2$`, '#ff0000', 'Second of m2'],
    ['CDS', 'hash', `^${second}@c$`, '#0000ff', 'Second of c']
  ] });
  const recipe = await buildSourceRecipe(options);
  assert.equal(recipe.available, true, recipe.unavailableReason);
  const table = recipe.generatedFiles.find(({ slot }) => slot === 'generatedFiles.source_recipe.specific_colors');
  // `m:2` is the File card's second record, the records table's second row; `c` is its fourth input (the windows of
  // `w` share the third).
  assert.equal(String(table.data), [
    `CDS\thash\t^${second}@record-2$\t#ff0000\tSecond of m2`,
    `CDS\thash\t^${second}@record-4$\t#0000ff\tSecond of c`,
    ''
  ].join('\n'));
  const cliSvg = await materializeAndRunRecipe(
    session, recipe, buildRunInfo({ mode: 'linear', sourceRecipe: recipe }), 'colour rule of a placed file card'
  );
  const { equal, differences, webSvg } = await compareWithWebRender(session, cliSvg, 'placed-file-card-colours');
  assert.equal(equal, true, differences.join('\n'));
  const filled = (svg, color) => [...svg.matchAll(new RegExp(`<path [^>]*data-gbdraw-feature-id="([^"]+)"[^>]*fill="${color}"`, 'gi'))]
    .map((match) => match[1]);
  const cli = await readFile(cliSvg, 'utf8');
  assert.equal(filled(cli, '#ff0000').length, 1);
  assert.equal(filled(cli, '#0000ff').length, 1);
  assert.deepEqual([filled(cli, '#ff0000'), filled(cli, '#0000ff')], [filled(webSvg, '#ff0000'), filled(webSvg, '#0000ff')]);
});

test('a placed File card of several records with a comparison has no Source recipe yet (OV-3272)', async () => {
  // The Web's unit plan counts the File card as one unit (m, the two windows of w, c), while Python reads a
  // comparison's indexes over the drawn records (m:1, m:2, the windows of w, c): without the refusal, this m -> c
  // comparison would be written as `#1 -> #3`, and the CLI's third record is w.
  const { options } = await placedFileCardSession({ comparisons: [[0, 3]] });
  const drawnUnits = { recordKeys: ['m', 'w', 'w', 'c'], windowOf: [null, 'w', 'w', null] };
  const recipe = await buildSourceRecipe({ ...options, drawnUnits });
  assert.equal(recipe.available, false);
  assert.match(recipe.unavailableReason, /a comparison cannot be bound losslessly to its source records/);
});

// OV-3118: a request with feature density tracks writes each track as --feature_density_track words in Python's
// spelling (`FEATURE_DENSITY_TRACK_WORDS`, gbdraw/analysis/feature_density.py `feature_density_track_words`), and
// the slots that name them by id; Python reads the words back as the request's tracks. A value the words cannot
// carry (blank around it) makes the recipe unavailable.
test('feature density tracks have their CLI words in the Source recipe of either mode (OV-3118)', async () => {
  const tracks = [
    { id: 'genes', query: { featureType: 'gene', qualifier: null, value: null, minLength: null, maxLength: null },
      metric: 'count', style: 'heatmap', window: null, step: null, color: null, caption: null },
    { id: 'ltr', query: { featureType: 'repeat_region', qualifier: 'rpt_type', value: 'long_terminal=repeat', minLength: 100, maxLength: 7000 },
      metric: 'coverage', style: 'histogram', window: 50000, step: 10000, color: '#d9944a', caption: 'LTR, % of 50 kb' }
  ];
  for (const [session, mode] of [[circularGallerySession, 'circular'], [linearGallerySession, 'linear']]) {
    const request = structuredClone(session);
    request.renderRequest.diagramOptions.featureDensityTracks = structuredClone(tracks);
    const recipe = await buildSourceRecipe(request);
    assert.equal(recipe.available, true, `${mode}: ${recipe.unavailableReason}`);
    const groups = recipe.args.flatMap((arg, index) => {
      if (arg !== '--feature_density_track') return [];
      const end = recipe.args.findIndex((next, at) => at > index && next.startsWith('--'));
      return [recipe.args.slice(index + 1, end < 0 ? undefined : end)];
    });
    assert.deepEqual(groups[0], ['id=genes', 'type=gene', 'metric=count', 'style=heatmap'], mode);
    const read = await evaluatePython([
      'from gbdraw.analysis.feature_density import feature_density_track_words, parse_feature_density_words',
      'result = [list(feature_density_track_words(parse_feature_density_words(words, where="recipe"))) for words in payload]'
    ].join('\n'), groups);
    assert.deepEqual(read, groups, `${mode}: the words read back as themselves`);
    assert.equal(groups[1].includes('value=long_terminal=repeat'), true, JSON.stringify(groups[1]));
    assert.equal(groups[1].includes('caption=LTR, % of 50 kb'), true, JSON.stringify(groups[1]));
    await assertCliParserAccepts(recipe);

    // A custom stack's rows name their track (E4 recipe probe: the rows lost track=, and the CLI refused the recipe).
    const stacked = structuredClone(request);
    const row = (id, renderer, params = {}) => (mode === 'circular'
      ? { kind: 'circularTrackSlot', id, renderer, enabled: true, side: null, radius: null, width: null, z: 0, params, innerGapPx: null, outerGapPx: null }
      : { kind: 'linearTrackSlot', id, renderer, enabled: true, side: renderer === 'features' ? 'overlay' : 'below', height: null, spacing: null, z: 0, params });
    stacked.renderRequest.diagramOptions.tracks[`${mode}TrackSlots`] = [
      row('features', 'features'), row('density_genes', 'feature_density', { track: 'genes' }), row('ltr_cover', 'feature_density', { track: 'ltr' })
    ];
    stacked.renderRequest.diagramOptions.tracks[`${mode}TrackAxisIndex`] = 0;
    const stackedRecipe = await buildSourceRecipe(stacked);
    assert.equal(stackedRecipe.available, true, `${mode}: ${stackedRecipe.unavailableReason}`);
    const rows = stackedRecipe.args.filter((_, index) => stackedRecipe.args[index - 1] === `--${mode}_track_slot`);
    assert.deepEqual(rows.filter((token) => token.includes(':feature_density')).map((token) => token.match(/[@,]track=([^,]*)/)?.[1]),
      ['genes', 'ltr'], `${mode}: ${rows.join(' | ')}`);

    const blank = structuredClone(request);
    blank.renderRequest.diagramOptions.featureDensityTracks[1].caption = ' padded ';
    const unavailable = await buildSourceRecipe(blank);
    assert.equal(unavailable.available, false, mode);
    assert.match(unavailable.unavailableReason, /feature density track/);
  }
});
