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
  pathToFileURL(join(tempRoot, 'js', 'app', 'match-sequences.js'))
);
const { projectCanonicalSessionRequest } = await import(
  pathToFileURL(join(tempRoot, 'js', 'services', 'session-request.js'))
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
    linearTrackSlotSchemaVersion: Number(fixture.version) <= 32 ? 1 : 2
  });

  assert.equal(projectedSession.files.linearSeqs.length, 4);
  assert.equal(
    projectedSession.files.linearSeqs[0].gb.name,
    'NC_004603.1__GCF_000196095.1_ASM19609v1_genomic.gbff'
  );
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
