import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { cp, readFile, stat, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { delimiter, join } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-run-info-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js'), join(tempDir, 'js'), { recursive: true });
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}', 'utf8');
const modulePath = join(tempDir, 'js', 'app', 'run-info.js');

const {
  buildRunInfo,
  buildSourceRecipe,
  isCliInvocationSessionExportable,
  quoteShellArg,
  reproducibilityLabel
} = await import(pathToFileURL(modulePath));
const { adoptCurrentSessionResources, createSessionResourceFileView } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'session-resource-backing.js'))
);
const { readFileText } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'file-content-cache.js'))
);
const { setResourcePayloadOwner } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'resource-payload-owner.js'))
);
const { readCanonicalResourceText } = await import(
  pathToFileURL(join(tempDir, 'js', 'services', 'session-request.js'))
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
const canonical = ({ mode, records, resources, comparisons = [], webFiles = {} }) => ({
  renderRequest: {
    schema: 6,
    mode,
    grouping: 'single',
    records,
    diagramOptions: baseOptions(mode),
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
const assertCliParserAccepts = (recipe) => {
  const visibleArgs = recipe.args.map((token) => recipe.fileMetadata.get(token)?.name || token);
  const script = [
    'import json, sys',
    'from gbdraw.circular import _get_args as circular_args',
    'from gbdraw.linear import _get_args as linear_args',
    '(linear_args if sys.argv[1] == "linear" else circular_args)(json.loads(sys.argv[2]))'
  ].join('; ');
  const parsed = spawnSync(
    process.env.PYTHON || 'python',
    ['-c', script, recipe.mode, JSON.stringify(visibleArgs)],
    { cwd: repoRoot, encoding: 'utf8' }
  );
  assert.equal(parsed.status, 0, parsed.stderr || parsed.stdout);
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
  const output = await stat(join(runDir, `${args[outputIndex + 1]}.svg`));
  assert.ok(output.size > 0, `${fixture}: CLI emitted an empty SVG`);
};

assert.equal(quoteShellArg('simple.tsv'), 'simple.tsv');
assert.equal(quoteShellArg('two words.tsv'), "'two words.tsv'");
assert.equal(quoteShellArg("Bob's genome.gbk"), "'Bob'\\''s genome.gbk'");
assert.equal(quoteShellArg(''), "''");

let circularGallerySession = null;
let linearGallerySession = null;
for (const { fixture, mode } of [
  { fixture: 'HmmtDNA_ATskew.gbdraw-session.json', mode: 'circular' },
  { fixture: 'lambda_basic_linear.gbdraw-session.json', mode: 'linear' }
]) {
  const session = JSON.parse(await readFile(
    join(repoRoot, 'gbdraw', 'web', 'gallery', 'sessions', fixture),
    'utf8'
  ));
  if (mode === 'circular') circularGallerySession = session;
  else linearGallerySession = session;
  const sourceRecipe = await buildSourceRecipe(session);
  assert.equal(sourceRecipe.available, true, `${fixture}: ${sourceRecipe.unavailableReason}`);
  assert.equal(sourceRecipe.mode, mode);
  assert.ok(sourceRecipe.args.includes('--gbk'), `${fixture}: expected direct GenBank input`);
  assert.ok(sourceRecipe.args.includes('-o'), `${fixture}: expected output prefix`);
  assert.ok(sourceRecipe.args.includes('-f'), `${fixture}: expected output formats`);
  assertCliParserAccepts(sourceRecipe);
  const info = buildRunInfo({ mode: sourceRecipe.mode, sourceRecipe });
  await materializeAndRunRecipe(session, sourceRecipe, info, fixture);
}

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
          region: null,
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
        readResourceText: (resourceId) => readCanonicalResourceText(session.resources, resourceId)
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
  assertCliParserAccepts(emptyRecipe);
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
  assertCliParserAccepts(nonEmptyRecipe);
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
  assertCliParserAccepts(stringRecipe);
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
  assertCliParserAccepts(objectRecipe);

  const unsupportedChild = structuredClone(objectSlot);
  unsupportedChild.renderRequest.diagramOptions.tracks.circularTrackSlots[0].params.nt = 'GC';
  const unsupportedRecipe = await buildSourceRecipe(unsupportedChild);
  assert.equal(unsupportedRecipe.available, false);
  assert.match(unsupportedRecipe.unavailableReason, /params\.nt|losslessly/i);
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
      region: null,
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
  assertCliParserAccepts(recipe);
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
    region: null,
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

{
  const sourceRecipe = await buildSourceRecipe(circularCanonical);
  assert.equal(sourceRecipe.available, true);
  assertCliParserAccepts(sourceRecipe);
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
      region: null,
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
  assertCliParserAccepts(sourceRecipe);
  assert.ok(sourceRecipe.args.includes('--records_table'));
  assert.ok(sourceRecipe.args.includes('--multi_record_canvas'));
  assert.equal(sourceRecipe.args.includes('--multi_record_position'), false);
  const recordsTable = sourceRecipe.generatedFiles.find((file) => file.name === 'records.tsv');
  assert.ok(recordsTable);
  assert.match(recordsTable.data, /two-records\.gbk\tALPHA\t\talpha\t\t0\t1\t1\t/);
  assert.match(recordsTable.data, /two-records\.gbk\tBETA\t\tbeta\t\t0\t2\t2\t/);
}

const linearCanonical = canonical({
  mode: 'linear',
  records: [
    {
      recordKey: 'alpha',
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'record-1-genbank' },
      selector: null,
      region: null,
      presentation: presentation()
    },
    {
      recordKey: 'beta',
      cardinality: 'exactly_one',
      source: { kind: 'genbank', resourceId: 'record-2-genbank' },
      selector: null,
      region: null,
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
  const direct = await buildSourceRecipe(linearCanonical);
  assert.equal(direct.available, true);
  assertCliParserAccepts(direct);
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
  assert.match(table.data, /alpha\.gbk.*\t1\t1\n/);
  assert.match(table.data, /beta\.gbk.*\t1\t2\n/);
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
  assert.match(info.reproducibility.notes.join('\n'), /\.gbdraw-session\.json/);
  assert.equal(isCliInvocationSessionExportable(info.invocation), false);
}

{
  const info = buildRunInfo({
    mode: 'linear',
    sourceRecipe: {
      available: false,
      args: [],
      fileMetadata: new Map(),
      generatedFiles: [],
      unavailableReason: 'Source recipe unavailable: source-file provenance is missing.'
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
  assert.equal(isCliInvocationSessionExportable(info.invocation), true);
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
  assert.match(info.reproducibility.notes.join('\n'), /session restore/);
  assert.equal(info.sessionCommand, '');
  assert.equal(isCliInvocationSessionExportable(info.invocation), false);
}
