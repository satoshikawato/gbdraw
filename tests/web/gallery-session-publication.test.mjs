import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { readFile } from 'node:fs/promises';
import { gunzipSync } from 'node:zlib';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const { resolveLinearComparisonPlan } = await import(
  '../../gbdraw/web/js/app/linear-comparisons.js'
);

const {
  applyDerivedCachePublicationPolicy,
  createGallerySessionPublication
} = await import('../../gbdraw/web/js/services/gallery-session-publication.js');
const {
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

const loadSession = async (name) => {
  const bytes = await readFile(`${sessionRoot}/${name}`);
  const decoded = bytes[0] === 0x1f && bytes[1] === 0x8b
    ? gunzipSync(bytes)
    : bytes;
  return JSON.parse(decoded.toString('utf8'));
};

for (const name of sessionNames) {
  const source = await loadSession(name);
  const committedBefore = JSON.stringify(source.renderRequest);
  const result = await prepareGallerySessionForPublication(source);
  assert.equal(result.session.version, 40, name);
  assert.equal(result.session.renderRequest.schema, 6, name);
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
  if (name === 'tobacco-chloroplast.gbdraw-session.json') {
    assert.equal(result.session.config.rules.length, 71, name);
    assert.deepEqual(result.session.config.qualifierPriorityRules, [
      { feat: 'CDS', order: 'gene,old_locus_tag' }
    ], name);
  }
}

const lambda = await loadSession('lambda_basic_linear.gbdraw-session.json');
assert.equal(admitGallerySession(lambda), lambda);
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
  'publication must preserve schema-5 materialized record cardinality'
);
assert.equal(
  lambdaPrepared.equivalence.actual.digest,
  alteredPrepared.equivalence.actual.digest
);
assert.deepEqual(alteredPrepared.session.cliInvocation, alteredProvenance.cliInvocation);

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

for (const version of [27, 30, 34, 38, 41]) {
  assert.throws(
    () => admitGallerySession({ ...lambda, version }),
    /supports current version 40 or historical versions 31-33\/39/
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
