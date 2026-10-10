import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { execFileSync } from 'node:child_process';
import { mkdtemp, readFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { gunzipSync } from 'node:zlib';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const { resolveLinearComparisonPlan } = await import(
  '../../gbdraw/web/js/services/linear-comparisons.js'
);

const {
  applyDerivedCachePublicationPolicy,
  createGallerySessionPublication
} = await import('../../gbdraw/web/js/services/gallery-session-publication.js');
const {
  CANONICAL_REQUEST_SCHEMA,
  assertCanonicalRenderRequestsEquivalent,
  buildCanonicalRenderRequest,
  buildCanonicalRequestState,
  compareCanonicalRenderRequests,
  promoteCanonicalRenderRequestToCurrent,
  projectCanonicalSessionRequest
} = await import('../../gbdraw/web/js/services/session-request.js');
const { promoteGallerySessionToCurrent } = await import(
  '../../gbdraw/web/js/services/gallery-session-migration.js'
);
const { FEATURE_CATALOG_SCHEMA } = await import('../../gbdraw/web/js/services/feature-catalog.js');
// Session 46 keeps each diagram mode's draft in `modes` (PD-OI-086).
const { MODE_SCOPED_SESSION_VERSION: CURRENT_SESSION_VERSION } = await import(
  '../../gbdraw/web/js/services/session-authority.js'
);

const {
  admitGallerySession,
  finalizeGallerySessionPublication,
  prepareGallerySessionForPublication,
  validateGalleryPublicationReadiness
} = createGallerySessionPublication({
  promoteSession: promoteGallerySessionToCurrent,
  assertRequestsEquivalent: assertCanonicalRenderRequestsEquivalent,
  buildRequest: buildCanonicalRenderRequest,
  buildRequestState: buildCanonicalRequestState,
  promoteRequest: promoteCanonicalRenderRequestToCurrent,
  projectRequest: projectCanonicalSessionRequest,
  resolveComparisonPlan: resolveLinearComparisonPlan
});

const sessionRoot = 'gbdraw/web/gallery/sessions';
const examples = JSON.parse(
  await readFile('gbdraw/web/gallery/examples.json', 'utf8')
);
const sessionNames = examples.map((example) => String(example.session).split('/').pop());

assert.equal(sessionNames.length, 10);
assert.equal(new Set(sessionNames).size, 10);

const loadSessionFile = async (path) => {
  const bytes = await readFile(path);
  const decoded = bytes[0] === 0x1f && bytes[1] === 0x8b
    ? gunzipSync(bytes)
    : bytes;
  return JSON.parse(decoded.toString('utf8'));
};
const loadSession = (name) => loadSessionFile(`${sessionRoot}/${name}`);
// The released Session 44 Gallery files (first-parent main fe6861f0), kept to
// exercise the promotion of released Gallery Sessions to the current writer.
const loadReleasedSession = (name) => loadSessionFile(`tests/fixtures/sessions/${name}`);

for (const name of sessionNames) {
  const source = await loadSession(name);
  const committedBefore = JSON.stringify(source.renderRequest);
  const result = await prepareGallerySessionForPublication(source);
  assert.equal(result.session.version, CURRENT_SESSION_VERSION, name);
  assert.equal(result.session.editorState.featureCatalog.schema, FEATURE_CATALOG_SCHEMA, name);
  // Publication writes the slice of the Session's mode only; the other mode
  // takes its defaults (plan 4.3). The draft has no home outside `modes`.
  const mode = result.session.renderRequest.mode;
  assert.deepEqual(Object.keys(result.session.modes), [mode], name);
  for (const field of ['config', 'features']) assert.equal(Object.hasOwn(result.session, field), false, name);
  const slice = result.session.modes[mode];
  assert.equal(typeof slice.ui.layoutPreferences, 'object', name);
  // A Session written by the CLI carries no Web draft rows.
  assert.deepEqual(slice.features?.featureOverrides ?? {}, {}, name);
  for (const field of ['featureVisibilityOverrides', 'labelVisibilityOverrides', 'labelTextFeatureOverrides',
    'labelTextFeatureOverrideSources']) assert.equal(Object.hasOwn(slice.features || {}, field), false, name);
  assert.equal(result.session.renderRequest.schema, CANONICAL_REQUEST_SCHEMA, name);
  const plan = result.session.renderRequest.layout?.similarityAlignment;
  if (plan) {
    assert.equal(plan.schema, 2, name);
    assert.equal(Object.hasOwn(plan, 'mode'), false, name);
    assert.deepEqual(new Set(plan.records.map(({ recordKey }) => recordKey)),
      new Set(result.session.renderRequest.records.map(({ recordKey }) => recordKey)), name);
    for (const decision of plan.records) {
      assert.deepEqual(Object.keys(decision).sort(),
        ['anchor', 'rationale', 'recordKey', 'status'], name);
    }
  }
  assert.equal(
    Object.hasOwn(result.session.orthogroupState || {}, 'selectedOrthogroupAlignmentFeature'),
    false,
    name
  );
  assert.equal(result.equivalence.equivalent, true, name);
  assert.equal(JSON.stringify(source.renderRequest), committedBefore, name);
  assert.equal(result.session.cliInvocation, source.cliInvocation, name);
  assert.deepEqual(
    result.session.renderRequest.output.interactiveMetadataPolicy,
    source.renderRequest.output.interactiveMetadataPolicy,
    name
  );
  const readiness = await validateGalleryPublicationReadiness(result.session);
  assert.equal(readiness.equivalence.equivalent, true, name);
  // The slice holds no field of the other mode alone (registry `own`).
  const { form: publishedForm, adv: publishedAdv } = slice.config;
  if (mode === 'circular') {
    for (const field of ['linear_show_replicon', 'linear_accession_visibility', 'linear_length_visibility']) {
      assert.equal(Object.hasOwn(publishedAdv, field), false, `${name} ${field}`);
    }
  } else {
    assert.equal(Object.hasOwn(publishedForm, 'multi_record_canvas'), false, name);
  }
  if (name === 'tobacco-chloroplast.gbdraw-session.json') {
    assert.equal(slice.config.rules.length, 71, name);
    assert.deepEqual(slice.config.qualifierPriorityRules, [
      { feat: 'CDS', order: 'gene,old_locus_tag' }
    ], name);
  }
}

const lambda = await loadReleasedSession('lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
const admittedLambda = admitGallerySession(lambda);
assert.equal(admittedLambda.version, CURRENT_SESSION_VERSION);
assert.equal(admittedLambda.editorState.featureCatalog.schema, FEATURE_CATALOG_SCHEMA);
assert.deepEqual(admittedLambda.results, lambda.results);
assert.equal(lambda.version, 44);
assert.equal(lambda.editorState.featureCatalog.schema, 4);
const releasedCurrent = await loadReleasedSession('HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz');
assert.equal(releasedCurrent.version, 44);
assert.equal(releasedCurrent.renderRequest.schema, 8);
assert.equal(releasedCurrent.editorState.featureCatalog.schema, 4);
const admittedReleasedCurrent = admitGallerySession(releasedCurrent);
// Admission promotes the released schema-8 request to the current writer.
assert.equal(admittedReleasedCurrent.renderRequest.schema, CANONICAL_REQUEST_SCHEMA);
assert.equal(admittedReleasedCurrent.editorState.featureCatalog.schema, FEATURE_CATALOG_SCHEMA);
assert.deepEqual(admittedReleasedCurrent.results, releasedCurrent.results);
const alteredProvenance = structuredClone(lambda);
alteredProvenance.cliInvocation = {
  ...alteredProvenance.cliInvocation,
  args: ['--untrusted-publication-option', '--record_label', 'not authority']
};
const [lambdaPrepared, alteredPrepared] = await Promise.all([
  prepareGallerySessionForPublication(lambda),
  prepareGallerySessionForPublication(alteredProvenance)
]);
assert.deepEqual(
  lambdaPrepared.session.renderRequest.records.map(({ cardinality }) => cardinality),
  ['exactly_one'],
  'publication must preserve materialized record cardinality'
);
assert.equal(
  lambdaPrepared.equivalence.actual.digest,
  alteredPrepared.equivalence.actual.digest
);
assert.deepEqual(alteredPrepared.session.cliInvocation, alteredProvenance.cliInvocation);

// Historical full configs omit the new zero-tolerance default; CLI replay writes it.
const fullConfigRequest = structuredClone(lambdaPrepared.session.renderRequest);
fullConfigRequest.diagramOptions.config = { canvas: {} };
delete fullConfigRequest.diagramOptions.configOverrides;
for (const tolerance of [0, 1]) {
  const replayed = structuredClone(fullConfigRequest);
  replayed.diagramOptions.config.canvas.feature_overlap_tolerance_bp = tolerance;
  const comparison = await compareCanonicalRenderRequests({
    expectedRequest: fullConfigRequest, expectedResources: lambdaPrepared.session.resources,
    actualRequest: replayed, actualResources: lambdaPrepared.session.resources
  });
  assert.equal(comparison.equivalent, tolerance === 0);
  if (tolerance) assert.equal(comparison.differences[0].path,
    '$.diagramOptions.config.canvas.feature_overlap_tolerance_bp');
}

const lambdaWithUnusedComparisonDefaults = structuredClone(lambdaPrepared.session.renderRequest);
Object.assign(lambdaWithUnusedComparisonDefaults.diagramOptions, {
  evalue: 0.01,
  bitscore: 50,
  identity: 0,
  alignmentLength: 0
});
assert.equal((await compareCanonicalRenderRequests({
  expectedRequest: lambdaPrepared.session.renderRequest,
  expectedResources: lambdaPrepared.session.resources,
  actualRequest: lambdaWithUnusedComparisonDefaults,
  actualResources: lambdaPrepared.session.resources
})).equivalent, true);

// Without -b or a protein mode the CLI writes its protein settings as a
// disabled pipeline with no pairs; it draws no comparison (OV-224).
assert.deepEqual(lambdaPrepared.session.renderRequest.comparisons, []);
const lambdaWithCliProteinSettings = structuredClone(lambdaPrepared.session.renderRequest);
lambdaWithCliProteinSettings.comparisons = [{"kind": "generatedProteinComparison", "mode": "none", "pairs": [], "settings": {"collinearityParams": {"kind": "lossless", "parameters": {"minAnchors": 1, "maxUnitGap": 0, "maxDiagonalDrift": 0, "maxConflicts": 1, "mergeOrientation": "either"}}, "collinearityUnitMode": "auto", "collinearityAnchorMode": "rbh", "collinearitySearchScope": "adjacent", "collinearityColorMode": "orientation", "losatpBin": "losat", "ncbiBlastpBin": null, "losatpThreads": null, "proteinBlastpMaxHits": 5, "proteinBlastpCandidateLimit": null, "orthogroupMembershipMode": "anchor_core_v1", "orthogroupMemberMaxHits": null, "collinearInferOrthogroups": true, "collinearMaxParalogLinksPerOrthogroup": 2}}];
assert.equal((await compareCanonicalRenderRequests({
  expectedRequest: lambdaWithCliProteinSettings,
  expectedResources: lambdaPrepared.session.resources,
  actualRequest: lambdaPrepared.session.renderRequest,
  actualResources: lambdaPrepared.session.resources
})).equivalent, true);

const lambdaWithDefaultColorFileAlias = structuredClone(lambdaPrepared.session.renderRequest);
lambdaWithDefaultColorFileAlias.diagramOptions.colors.defaultColorsFile =
  lambdaWithDefaultColorFileAlias.diagramOptions.colors.defaultColors;
lambdaWithDefaultColorFileAlias.diagramOptions.colors.defaultColors = null;
const lambdaWithDefaultColorFileAliasResources = structuredClone(lambdaPrepared.session.resources);
lambdaWithDefaultColorFileAliasResources[
  lambdaWithDefaultColorFileAlias.diagramOptions.colors.defaultColorsFile.resourceId
].kind = 'default-colors-file';
assert.equal((await compareCanonicalRenderRequests({
  expectedRequest: lambdaPrepared.session.renderRequest,
  expectedResources: lambdaPrepared.session.resources,
  actualRequest: lambdaWithDefaultColorFileAlias,
  actualResources: lambdaWithDefaultColorFileAliasResources
})).equivalent, true);

// OV-398: a published CLI request states the default definition-line weight
// `normal`; the Web draft holds it as null after Load and omits it.
const nameWeight = 'objects.definition.linear.line_styles.name.font_weight';
const defaultWeight = structuredClone(lambdaPrepared.session.renderRequest);
delete defaultWeight.diagramOptions.configOverrides[nameWeight];
for (const [weight, equivalent] of [['normal', true], ['bold', false]]) {
  const explicitWeight = structuredClone(defaultWeight);
  explicitWeight.diagramOptions.configOverrides[nameWeight] = weight;
  assert.equal((await compareCanonicalRenderRequests({
    expectedRequest: explicitWeight,
    expectedResources: lambdaPrepared.session.resources,
    actualRequest: defaultWeight,
    actualResources: lambdaPrepared.session.resources
  })).equivalent, equivalent, weight);
}

// A Web draft after Load states every default feature rendering; a CLI
// request states only the ones it changes (OV-398).
for (const [shape, equivalent] of [['arrow', true], ['rectangle', false]]) {
  const statedShapes = structuredClone(lambdaPrepared.session.renderRequest);
  statedShapes.diagramOptions.featureShapes = { ...statedShapes.diagramOptions.featureShapes, CDS: shape };
  const omittedShapes = structuredClone(lambdaPrepared.session.renderRequest);
  delete omittedShapes.diagramOptions.featureShapes.CDS;
  assert.equal((await compareCanonicalRenderRequests({
    expectedRequest: omittedShapes,
    expectedResources: lambdaPrepared.session.resources,
    actualRequest: statedShapes,
    actualResources: lambdaPrepared.session.resources
  })).equivalent, equivalent, shape);
}

const ungeneratedDraft = structuredClone(lambda);
ungeneratedDraft.config.adv.arrow_shaft_width_ratio = 0.5;
await assert.rejects(
  prepareGallerySessionForPublication(ungeneratedDraft),
  /shaft_width_ratio/
);

for (const field of ['cli_circular_track_order', 'cli_circular_track_slots']) {
  const invalid = structuredClone(lambda);
  invalid.config.adv[field] = [];
  assert.throws(
    () => admitGallerySession(invalid),
    new RegExp(`config\\.adv.*${field}`)
  );
}

for (const version of [27, 30, 34, 38, 43]) {
  assert.throws(
    () => admitGallerySession({ ...lambda, version }),
    /supports current version 46 or historical versions 31-33\/39-44/
  );
}

const changedRequest = structuredClone(lambdaPrepared.session.renderRequest);
changedRequest.diagramOptions.configOverrides['objects.features.arrow_geometry.shaft_width_ratio'] = 0.5;
const comparison = await compareCanonicalRenderRequests({
  expectedRequest: lambdaPrepared.session.renderRequest,
  expectedResources: lambdaPrepared.session.resources,
  actualRequest: changedRequest,
  actualResources: lambdaPrepared.session.resources
});
assert.equal(comparison.equivalent, false);
assert.ok(comparison.differences.some(
  (difference) => difference.path.endsWith('.objects.features.arrow_geometry.shaft_width_ratio')
));

const changedMetadataPolicy = structuredClone(lambdaPrepared.session.renderRequest);
changedMetadataPolicy.output.interactiveMetadataPolicy = 'omit';
const metadataPolicyComparison = await compareCanonicalRenderRequests({
  expectedRequest: lambdaPrepared.session.renderRequest,
  expectedResources: lambdaPrepared.session.resources,
  actualRequest: changedMetadataPolicy,
  actualResources: lambdaPrepared.session.resources
});
assert.equal(metadataPolicyComparison.equivalent, false);
assert.ok(metadataPolicyComparison.differences.some(
  (difference) => difference.path.endsWith('.interactiveMetadataPolicy')
));

const comparisonTableRequest = {
  schema: 5,
  mode: 'linear',
  grouping: 'single',
  records: [],
  diagramOptions: {},
  comparisons: [{
    kind: 'precomputedProteinComparison',
    resourceId: 'comparison-table',
    queryRecordIndex: 0,
    subjectRecordIndex: 1
  }],
  output: { formats: ['svg'], interactiveMetadataPolicy: 'auto' }
};
const comparisonResource = (row) => ({
  'comparison-table': {
    kind: 'canonical-tsv',
    encoding: 'base64',
    data: Buffer.from(`query\tsubject\tevalue\n${row}\n`, 'utf8').toString('base64')
  }
});
const replayNumberSpelling = await compareCanonicalRenderRequests({
  expectedRequest: comparisonTableRequest,
  expectedResources: comparisonResource('q\ts\t1.1200000000000001e-132'),
  actualRequest: comparisonTableRequest,
  actualResources: comparisonResource('q\ts\t1.12e-132'),
  normalizeReplayGeneratedResources: true
});
assert.equal(replayNumberSpelling.equivalent, true);
const changedComparisonValue = await compareCanonicalRenderRequests({
  expectedRequest: comparisonTableRequest,
  expectedResources: comparisonResource('q\ts\t1.12e-132'),
  actualRequest: comparisonTableRequest,
  actualResources: comparisonResource('q\ts\t1.13e-132'),
  normalizeReplayGeneratedResources: true
});
assert.equal(changedComparisonValue.equivalent, false);
const changedComparisonIdentity = await compareCanonicalRenderRequests({
  expectedRequest: comparisonTableRequest,
  expectedResources: comparisonResource('q\ts\t1.12e-132'),
  actualRequest: comparisonTableRequest,
  actualResources: comparisonResource('q2\ts\t1.12e-132'),
  normalizeReplayGeneratedResources: true
});
assert.equal(changedComparisonIdentity.equivalent, false);

const finalized = await finalizeGallerySessionPublication({
  prepared: lambdaPrepared.session,
  replayed: structuredClone(lambdaPrepared.session)
});
assert.equal(finalized.equivalence.equivalent, true);
assert.deepEqual(finalized.session.cliInvocation, lambda.cliInvocation);

const replayWithUnreferencedResource = structuredClone(lambdaPrepared.session);
replayWithUnreferencedResource.resources['unused-replay-resource'] = {
  kind: 'canonical-tsv',
  encoding: 'base64',
  data: Buffer.from('unused\n', 'utf8').toString('base64')
};
const finalizedWithoutUnusedResource = await finalizeGallerySessionPublication({
  prepared: lambdaPrepared.session,
  replayed: replayWithUnreferencedResource
});
assert.equal(
  finalizedWithoutUnusedResource.session.resources['unused-replay-resource'],
  undefined
);

const replayWithArtifactResource = structuredClone(lambdaPrepared.session);
replayWithArtifactResource.resources['replay-artifact-resource'] = {
  kind: 'runtime-artifact',
  name: 'artifact.txt',
  type: 'text/plain',
  size: 9,
  lastModified: 0,
  encoding: 'base64',
  data: Buffer.from('artifact\n', 'utf8').toString('base64')
};
replayWithArtifactResource.runMetadata = {
  ...(replayWithArtifactResource.runMetadata || {}),
  retainedArtifact: { resourceId: 'replay-artifact-resource' }
};
const finalizedWithArtifactResource = await finalizeGallerySessionPublication({
  prepared: lambdaPrepared.session,
  replayed: replayWithArtifactResource
});
assert.deepEqual(
  finalizedWithArtifactResource.session.resources['replay-artifact-resource'],
  replayWithArtifactResource.resources['replay-artifact-resource']
);

const regenerableCache = {
  renderRequest: { comparisons: [{ kind: 'generatedProteinComparison' }] },
  proteinIdentityManifest: { schema: 2 },
  losatCache: {
    entries: [{
      schema: 4,
      kind: 'raw-losat',
      program: 'blastp',
      idEncoding: 'runtime-handle-v1',
      queryProteinSetHash: 'query-hash',
      subjectProteinSetHash: 'subject-hash'
    }]
  },
  losatDerivedCache: { entries: [{ payload: 'x'.repeat(100) }] }
};
const compactPublication = applyDerivedCachePublicationPolicy(
  regenerableCache,
  { limitBytes: 32 }
);
assert.deepEqual(compactPublication.losatDerivedCache.entries, []);
assert.equal(regenerableCache.losatDerivedCache.entries.length, 1);
const unprovenCache = structuredClone(regenerableCache);
unprovenCache.renderRequest.comparisons = [];
assert.equal(
  applyDerivedCachePublicationPolicy(unprovenCache, { limitBytes: 32 }),
  unprovenCache
);

const collisionReplay = structuredClone(lambdaPrepared.session);
const collisionId = Object.keys(collisionReplay.resources)[0];
collisionReplay.resources[collisionId].data = 'QQ==';
assert.rejects(
  finalizeGallerySessionPublication({
    prepared: lambdaPrepared.session,
    replayed: collisionReplay
  }),
  new RegExp(`resources\\.${collisionId}`)
);

// OV-217: a Web-authored Session 46 slice may hold a record display draft with
// no `scope`; the committed mode's slice must still hold it after publication.
{
  const scopeless = await loadSession('vibrio-harveyi-group-collinear.gbdraw-session.json.gz');
  const mode = scopeless.renderRequest.mode;
  const source = scopeless.modes[mode].config.recordDisplayDrafts;
  assert.ok(source.length > 0);
  assert.ok(source.every((draft) => !Object.hasOwn(draft, 'scope')));
  const republished = await prepareGallerySessionForPublication(scopeless);
  const kept = republished.session.modes[mode].config.recordDisplayDrafts;
  assert.deepEqual(kept.map(({ scope, ...draft }) => draft), source);
  assert.ok(kept.every((draft) => draft.scope === undefined || draft.scope === mode));
  assert.equal(kept.length, source.length);
}

// A Session the CLI writes publishes as the Gallery entry of its command: a
// Circular batch request with one output per record (OV-267), the resolved
// default colors and the editor tables the CLI stores beside the files it read
// (OV-266), and a read-only `-b` comparison, which the rebuild inherits as
// Generate's Inherit does (D-04, OV-268). The published slice holds the CLI's
// visibility and label rows, which Web Load reads from it (OV-299). Circular
// uses the output format of the Gallery commands (interactive_svg).
{
  const { cases } = JSON.parse(await readFile('tests/fixtures/cli_session_cross_surface/cases.json', 'utf8'));
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-gallery-cli-'));
  try {
    for (const id of ['circular_records_tables', 'linear_tables', 'linear_blast']) {
      const { mode, args } = cases.find((entry) => entry.id === id);
      const file = path.join(directory, `${id}.gbdraw-session.json`);
      execFileSync('python', ['-m', 'gbdraw.cli', mode, ...args, '-o', path.join(directory, id),
        '-f', mode === 'circular' ? 'interactive_svg' : 'svg',
        '--session_output', file], { env: { ...process.env, PYTHONPATH: process.cwd() }, stdio: 'pipe', timeout: 1_800_000 });
      const source = JSON.parse(await readFile(file, 'utf8'));
      const { session, equivalence } = await prepareGallerySessionForPublication(source);
      assert.deepEqual(equivalence.differences, [], id);
      assert.equal((await validateGalleryPublicationReadiness(session)).equivalence.equivalent, true, id);
      // Load reads a stored draft's per-feature edits from the draft, so the
      // published slice holds the CLI's visibility and label rows (OV-299).
      const rows = async (option) => (args.includes(option) ? (await readFile(args[args.indexOf(option) + 1], 'utf8'))
        .split(/\r?\n/).filter((line) => line.trim()).length : 0);
      const { features } = session.modes[mode];
      assert.equal(features.featureVisibilityManualRules.length, await rows('--feature_visibility_table'), id);
      assert.equal(features.labelOverrideRows.length, await rows('--label_table'), id);
      if (mode === 'circular') {
        assert.deepEqual([source.renderRequest.grouping, Array.isArray(session.renderRequest.output)], ['batch', true], id);
        assert.deepEqual(session.renderRequest.output.map((output) => output.interactiveMetadataPolicy),
          source.renderRequest.output.map((output) => output.interactiveMetadataPolicy), id);
      }
      if (id === 'linear_blast') {
        assert.deepEqual(session.renderRequest.comparisons, source.renderRequest.comparisons, id);
        assert.deepEqual(session.modes.linear.config.linearComparisonPlan, { mode: 'none', defaultSource: 'losat', edges: [] }, id);
      }
    }
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
}
