import assert from 'node:assert/strict';
import test from 'node:test';
import { cp, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { gunzipSync } from 'node:zlib';

const repoRoot = process.cwd();
const sourceRoot = join(repoRoot, 'gbdraw', 'web', 'js');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-vibrio-coverage-'));
await cp(sourceRoot, join(tempRoot, 'js'), { recursive: true });
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');

const { analyzeCatalogSequenceSourceCoverage } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'match-sequences.js'))
);
const { projectCanonicalSessionRequest } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'session-request.js'))
);
const { adoptCurrentSessionResources } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'session-resource-backing.js'))
);
const { groupLinearSourceRecords } = await import(
  pathToFileURL(join(tempRoot, 'js', 'app', 'linear-sources.js'))
);
const { planLinearSourceRowMove } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'linear-record-layout.js'))
);

const fixturePath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'vibrio-harveyi-group-collinear.gbdraw-session.json.gz'
);
const fixture = JSON.parse(gunzipSync(await readFile(fixturePath)).toString('utf8'));

test('real Vibrio catalog covers every sequence consumer with a valid sparse catalog', () => {
  const coverage = analyzeCatalogSequenceSourceCoverage({
    mode: fixture.renderRequest.mode,
    catalogFeatureState: fixture.editorState.featureCatalog,
    renderRequest: fixture.renderRequest
  });

  assert.equal(coverage.complete, true);
  assert.deepEqual(coverage.missingConsumers, []);
  assert.deepEqual(coverage.invalidCatalogSources, []);
  assert.deepEqual(
    coverage.resolvedConsumers.map(({ expectedSource }) => expectedSource.recordIndex),
    [0, 1, 2, 3]
  );
  assert.deepEqual(coverage.displayedRecordsWithoutConsumers, []);
});

test('real Vibrio session still projects its embedded records and presentation settings', () => {
  const projectedSession = projectCanonicalSessionRequest({
    renderRequest: fixture.renderRequest,
    resources: fixture.resources,
    webFiles: fixture.webFiles,
    legacyFiles: fixture.files,
    storedConfig: fixture.config,
    fileBindings: fixture.cliInvocation?.fileBindings,
    linearTrackSlotSchemaVersion: Number(fixture.version) <= 32 ? 1 : 2,
    sessionResourceTable: adoptCurrentSessionResources(fixture.resources)
  });

  // Two multi-record GenBank Files, one per row: File order can move and the
  // Record Layout is not custom.
  const { linearSeqs } = projectedSession.files;
  assert.deepEqual(
    linearSeqs.map((sequence) => [sequence.gb.name, sequence.region_record_id]),
    [
      ['GCF_000196095.1_ASM19609v1_genomic.gbff', 'NC_004603.1'],
      ['GCF_000196095.1_ASM19609v1_genomic.gbff', 'NC_004605.1'],
      ['GCF_000354175.2_ASM35417v2_genomic.gbff', 'NC_022349.1'],
      ['GCF_000354175.2_ASM35417v2_genomic.gbff', 'NC_022359.1']
    ]
  );
  const sourceGroups = groupLinearSourceRecords(linearSeqs);
  assert.deepEqual(sourceGroups.map((group) => group.records.length), [2, 2]);
  const move = planLinearSourceRowMove({
    sourceGroups,
    entries: projectedSession.config.linearRecordLayout.rows,
    sourceIndex: 0,
    direction: 1
  });
  assert.equal(move.allowed, true, move.reason);
  assert.equal(projectedSession.config.adv.block_stroke_width, 0);
  assert.equal(projectedSession.config.adv.line_stroke_width, 1);
  assert.equal(projectedSession.config.adv.axis_stroke_width, 2);
  assert.equal(projectedSession.config.adv.def_font_size, 16);
  assert.equal(projectedSession.config.filterMode, 'None');
});

test('real current Vibrio Gallery session retains its collinearity unit gap', () => {
  const comparison = fixture.renderRequest.comparisons.find(
    (candidate) => candidate.kind === 'generatedProteinComparison'
  );

  assert.equal(
    comparison.settings.collinearityParams.parameters.maxUnitGap,
    2
  );
});
