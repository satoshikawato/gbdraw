import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { execFileSync } from 'node:child_process';
import { mkdtemp, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { gunzipSync } from 'node:zlib';
import path from 'node:path';
import { readSavedSession } from './helpers/session-documents.mjs';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const { pythonDraftSplit } = await import('./helpers/python-evaluator.mjs');
const { pythonLinearUnitPlan } = await import('./helpers/linear-unit-plan.mjs');
// Inherit carries a saved comparison by Python's plan of the drawn units (`planLinearUnits`).
const planLinearUnits = async (payload) => ({ result: pythonLinearUnitPlan(payload) });
// Publication's draft split is Python's (the bridge runs `split_draft_into_modes`); each call's context is kept.
const splitContexts = [];
const splitDraft = async (draft, context) => {
  splitContexts.push(context);
  return pythonDraftSplit(draft, context);
};

const { resolveLinearComparisonPlan } = await import(
  '../../gbdraw/web/js/services/linear-comparisons.js'
);

const {
  applyDerivedCachePublicationPolicy,
  createGallerySessionPublication,
  galleryPublicationView
} = await import('../../gbdraw/web/js/services/gallery-session-publication.js');
const {
  CANONICAL_REQUEST_SCHEMA,
  assertCanonicalRenderRequestsEquivalent,
  buildCanonicalRenderRequest,
  buildCanonicalRequestState,
  compareCanonicalRenderRequests,
  readPublishedRenderRequest,
  projectCanonicalSessionRequest
} = await import('../../gbdraw/web/js/services/session-request.js');
const { FEATURE_CATALOG_SCHEMA } = await import('../../gbdraw/web/js/services/feature-catalog.js');
const { buildSimilarityAlignmentResetReceipt, validateSimilarityAlignmentResetReceipt } = await import(
  '../../gbdraw/web/js/services/session-active-config-contract.js'
);
const { SESSION_VERSION: CURRENT_SESSION_VERSION } = await import(
  '../../gbdraw/web/js/contracts/schema.generated.js'
);

const {
  admitGallerySession,
  finalizeGallerySessionPublication,
  prepareGallerySessionForPublication,
  validateGalleryPublicationReadiness
} = createGallerySessionPublication({
  assertRequestsEquivalent: assertCanonicalRenderRequestsEquivalent,
  buildRequest: buildCanonicalRenderRequest,
  buildRequestState: buildCanonicalRequestState,
  readRequest: readPublishedRenderRequest,
  projectRequest: projectCanonicalSessionRequest,
  resolveComparisonPlan: resolveLinearComparisonPlan,
  planLinearUnits,
  splitDraft
});

const sessionRoot = 'gbdraw/web/gallery/sessions';
const paletteDefinitions = JSON.parse(await readFile('gbdraw/web/gallery/palettes/palettes.json', 'utf8')).palettes;
// The rows of a published Session's Default colors table.
const projectedRows = (session) => {
  const colors = session.renderRequest.diagramOptions.colors;
  const ref = colors.defaultColorsFile || colors.defaultColors;
  if (!ref) return {};
  const text = Buffer.from(session.resources[ref.resourceId].data, 'base64').toString('utf8');
  return Object.fromEntries(text.split('\n').filter((line) => line.includes('\t')).map((line) => line.split('\t')));
};
const examples = JSON.parse(
  await readFile('gbdraw/web/gallery/examples.json', 'utf8')
);
const sessionNames = examples.map((example) => String(example.session).split('/').pop());

assert.equal(sessionNames.length, 12);
assert.equal(new Set(sessionNames).size, 12);

// A Session file as the publication bridge reads it (`tools/publish_gallery_session.mjs`): a Session 51's one
// drawn drawing as its one-drawing view, an older Session as it is.
const loadSessionFile = async (path) => {
  return galleryPublicationView(await readSavedSession(new Blob([await readFile(path)]))).view;
};
const loadSession = (name) => loadSessionFile(`${sessionRoot}/${name}`);
// The released Session 44 Gallery files (first-parent main fe6861f0), kept to
// exercise the promotion of released Gallery Sessions to the current writer.
// A released Session stores the request of its time (schema 9 for Sessions 46-50); publication reads it through
// Load's request reader, which promotes it, and publishes the current schema only.
const loadReleasedSession = (name) => loadSessionFile(`tests/fixtures/sessions/${name}`);

for (const name of sessionNames) {
  const source = await loadSession(name);
  const committedBefore = JSON.stringify(source.renderRequest);
  splitContexts.length = 0;
  const result = await prepareGallerySessionForPublication(source);
  // The draft is split once, by the injected owner, in the committed mode.
  assert.deepEqual(splitContexts.map(({ committedMode }) => committedMode), [result.session.renderRequest.mode], name);
  assert.equal(result.session.version, CURRENT_SESSION_VERSION, name);
  assert.equal(result.session.editorState.featureCatalog.schema, FEATURE_CATALOG_SCHEMA, name);
  // Publication writes the slice of the Session's mode only; the other mode
  // takes its defaults (plan 4.3). The draft has no home outside `modes`.
  const mode = result.session.renderRequest.mode;
  assert.deepEqual(Object.keys(result.session.modes), [mode], name);
  for (const field of ['config', 'features']) assert.equal(Object.hasOwn(result.session, field), false, name);
  const slice = result.session.modes[mode];
  assert.equal(typeof slice.ui.layoutPreferences, 'object', name);
  // A slice's colors are complete, as Load reads them: the Default colors rows over the draft palette's colors, not
  // over the request's rows alone (OV-3257: a Gallery table holds the user's rows only).
  const sliceColors = slice.config.colors || {};
  if (Object.keys(sliceColors).length) {
    const palette = paletteDefinitions[slice.config.palette || 'default'];
    assert.deepEqual(Object.keys(palette).filter((key) => !Object.hasOwn(sliceColors, key)), [], name);
    for (const key of ['rRNA', 'tRNA', 'skew_high', 'pairwise_match_max']) {
      if (!Object.hasOwn(projectedRows(result.session), key)) assert.equal(sliceColors[key].toLowerCase(), palette[key].toLowerCase(), `${name} ${key}`);
    }
  }
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
  if (name === 'tobacco-chloroplast.gbdraw') {
    assert.equal(slice.config.rules.length, 19, name);
    assert.deepEqual(slice.config.qualifierPriorityRules, [
      { feat: 'CDS', order: 'gene,old_locus_tag' }
    ], name);
  }
}

// A Session 46 is admitted as the current Session: its Legend edits move from
// the shown caption to the row key (tests/fixtures/sessions/legend-row-key-vectors.json).
{
  const source = structuredClone(await loadSession(sessionNames[0]));
  const mode = source.renderRequest.mode;
  const slice = source.modes[mode];
  slice.editorState = { ...slice.editorState, legend: {
    ...slice.editorState?.legend,
    entries: [{ caption: 'Z', originalCaption: 'GC content', color: '#a0a0a0', featureIds: [] }],
    colorOverrides: { Z: '#ff0000' }
  } };
  const admitted = await admitGallerySession({ ...source, version: 46 });
  assert.equal(admitted.version, CURRENT_SESSION_VERSION);
  assert.deepEqual(admitted.modes[mode].editorState.legend.colorOverrides, { 'GC content': '#ff0000' });
}

// OV-1145: a Session 47 is admitted as current unless it stored protein
// results of a reverse-complemented record, which only Python converts.
{
  const source = structuredClone(await loadSession(sessionNames[0]));
  assert.equal((await admitGallerySession({ ...source, version: 47 })).version, CURRENT_SESSION_VERSION);
  const stored = await loadReleasedSession('BGC0000708-BGC0000713.v47.gbdraw-session.json.gz');
  await assert.rejects(admitGallerySession(stored), /cannot convert the protein results that Session 47 stored/);
}

// A Session 46-50 stores a request schema 9, whose cropped Linear card kept its record selector in the region (the
// v0.14.5 Web), and an alignment reset receipt bound to that request as stored. Admission promotes the request and
// rebinds the receipt as Load does (`readStoredSimilarityAlignmentResetReceipt`), so the published Session's receipt
// binds the published request.
{
  const source = structuredClone(await loadSession('BGC0000708-BGC0000713.gbdraw'));
  const window = { start: 1, end: 5000, reverseComplement: false, column: null };
  const selector = { kind: 'recordId', value: 'record-1-id' };
  const promoted = structuredClone(source.renderRequest);
  promoted.records[0] = { ...promoted.records[0], selector, windows: [window] };
  // The binding the v0.14.5 Web wrote: the record's own selector, null, and the region's span.
  const asStored = structuredClone(promoted);
  asStored.records[0].selector = null;
  const canonical = (renderRequest) => ({ renderRequest, resources: source.resources });
  const receipt = await buildSimilarityAlignmentResetReceipt({ before: canonical(asStored), after: canonical(asStored) });
  const { windows: _windows, ...unwindowed } = promoted.records[0];
  const storedRequest = { ...structuredClone(source.renderRequest), schema: 9, records: promoted.records.map((record, index) => {
    const { windows: _recordWindows, ...rest } = record;
    return index === 0
      ? { ...unwindowed, selector: null, region: { selector, start: window.start, end: window.end, reverseComplement: false } }
      : { ...rest, region: null };
  }) };
  const stored = { ...source, version: 50, renderRequest: storedRequest,
    editorState: { ...source.editorState, alignmentResetReceipt: receipt } };
  const admitted = await admitGallerySession(stored);
  assert.equal(admitted.renderRequest.schema, CANONICAL_REQUEST_SCHEMA);
  assert.deepEqual(admitted.renderRequest.records[0].selector, selector);
  const rebound = admitted.editorState.alignmentResetReceipt;
  assert.notEqual(rebound.binding, receipt.binding);
  assert.equal(await validateSimilarityAlignmentResetReceipt(rebound, canonical(admitted.renderRequest)), rebound);
  // A receipt that binds neither request is refused, as publication repairs nothing.
  const staleReceipt = { ...receipt, binding: '0'.repeat(64) };
  await assert.rejects(admitGallerySession({ ...stored, editorState: { ...stored.editorState, alignmentResetReceipt: staleReceipt } }),
    /cannot read the Similarity alignment reset receipt/);
}

// Gallery publication admits Sessions 46-51 only; an older Session is refused with the way to convert it.
for (const name of ['lambda_basic_linear.v31-schema1.gbdraw-session.json.gz', 'lambda_basic_linear.v44-schema8.gbdraw-session.json.gz']) {
  const released = await loadReleasedSession(name);
  await assert.rejects(admitGallerySession(released),
    /Gallery publication supports Sessions 46-51; Load a Session 27-44 in the app and save it \(Python converts it at Load\); received /);
}

const lambda = await loadSession('lambda_basic_linear.gbdraw');
const admittedLambda = await admitGallerySession(lambda);
assert.equal(admittedLambda.version, CURRENT_SESSION_VERSION);
assert.equal(admittedLambda.editorState.featureCatalog.schema, FEATURE_CATALOG_SCHEMA);
assert.deepEqual(admittedLambda.results, lambda.results);
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

// The same Default colors rows as a file or as a canonical table are two requests: a canonical table is complete (a
// missing type draws the fallback grey), a file applies over the palette (BGC holds a file table).
const bgcPrepared = await prepareGallerySessionForPublication(await loadSession('BGC0000708-BGC0000713.gbdraw'));
const bgcWithCanonicalDefaultColors = structuredClone(bgcPrepared.session.renderRequest);
bgcWithCanonicalDefaultColors.diagramOptions.colors.defaultColors = {
  ...bgcWithCanonicalDefaultColors.diagramOptions.colors.defaultColorsFile, representation: 'canonicalTsv'
};
bgcWithCanonicalDefaultColors.diagramOptions.colors.defaultColorsFile = null;
const bgcWithCanonicalDefaultColorsResources = structuredClone(bgcPrepared.session.resources);
bgcWithCanonicalDefaultColorsResources[
  bgcWithCanonicalDefaultColors.diagramOptions.colors.defaultColors.resourceId
].kind = 'canonical-tsv';
assert.equal((await compareCanonicalRenderRequests({
  expectedRequest: bgcPrepared.session.renderRequest,
  expectedResources: bgcPrepared.session.resources,
  actualRequest: bgcWithCanonicalDefaultColors,
  actualResources: bgcWithCanonicalDefaultColorsResources
})).equivalent, false);

// Overrides stored before the axis label style existed drew numbered ticks:
// they equal a rebuilt request that states numbered, and differ from range.
const unstatedAxisLabels = structuredClone(lambdaPrepared.session.renderRequest);
for (const leaf of ['objects.depth.axis_labels', 'objects.gc_content.axis_labels']) {
  delete unstatedAxisLabels.diagramOptions.configOverrides[leaf];
}
for (const [style, equivalent] of [['numbered', true], ['range', false]]) {
  const stated = structuredClone(unstatedAxisLabels);
  stated.diagramOptions.configOverrides['objects.depth.axis_labels'] = style;
  stated.diagramOptions.configOverrides['objects.gc_content.axis_labels'] = style;
  assert.equal((await compareCanonicalRenderRequests({
    expectedRequest: unstatedAxisLabels, expectedResources: lambdaPrepared.session.resources,
    actualRequest: stated, actualResources: lambdaPrepared.session.resources
  })).equivalent, equivalent, style);
}

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
// With one changed rendering, the CLI states only that one and the Web draft
// states it beside every default (review 3).
const { DEFAULT_FEATURE_RENDERINGS } = await import('../../gbdraw/web/js/utils/feature-rendering.js');
for (const [webShape, equivalent] of [['arrow', true], ['rectangle', false]]) {
  const cliShapes = structuredClone(lambdaPrepared.session.renderRequest);
  cliShapes.diagramOptions.featureShapes = { misc_feature: 'arrow' };
  const webShapes = structuredClone(lambdaPrepared.session.renderRequest);
  webShapes.diagramOptions.featureShapes = { ...DEFAULT_FEATURE_RENDERINGS, misc_feature: webShape };
  assert.equal((await compareCanonicalRenderRequests({
    expectedRequest: cliShapes,
    expectedResources: lambdaPrepared.session.resources,
    actualRequest: webShapes,
    actualResources: lambdaPrepared.session.resources
  })).equivalent, equivalent, `misc_feature ${webShape}`);
}

const ungeneratedDraft = structuredClone(lambda);
ungeneratedDraft.modes.linear.config.adv.arrow_shaft_width_ratio = 0.5;
await assert.rejects(
  prepareGallerySessionForPublication(ungeneratedDraft),
  /shaft_width_ratio/
);

for (const field of ['cli_circular_track_order', 'cli_circular_track_slots']) {
  const invalid = structuredClone(lambda);
  invalid.modes.linear.config.adv[field] = [];
  await assert.rejects(
    admitGallerySession(invalid),
    { code: 'INPUT_INVALID', context: { field: 'schema', reason: 'FIELDS' } }
  );
}

for (const version of [27, 30, 34, 38, 43]) {
  await assert.rejects(
    admitGallerySession({ ...lambda, version }),
    /Gallery publication supports Sessions 46-51;/
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
  schema: 9,
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
const changedComparisonValue = await compareCanonicalRenderRequests({
  expectedRequest: comparisonTableRequest,
  expectedResources: comparisonResource('q\ts\t1.12e-132'),
  actualRequest: comparisonTableRequest,
  actualResources: comparisonResource('q\ts\t1.13e-132')
});
assert.equal(changedComparisonValue.equivalent, false);
const changedComparisonIdentity = await compareCanonicalRenderRequests({
  expectedRequest: comparisonTableRequest,
  expectedResources: comparisonResource('q\ts\t1.12e-132'),
  actualRequest: comparisonTableRequest,
  actualResources: comparisonResource('q2\ts\t1.12e-132')
});
assert.equal(changedComparisonIdentity.equivalent, false);

const finalized = await finalizeGallerySessionPublication({
  prepared: lambdaPrepared.session,
  replayed: structuredClone(lambdaPrepared.session)
});
assert.equal(finalized.equivalence.equivalent, true);
assert.deepEqual(finalized.session.cliInvocation, lambda.cliInvocation);

// OV-1144: the published legacy protein artifacts are the replay's, also when the replay writes none.
{
  const evidence = (/** @type {string} */ key) => ({
    proteinDerivedEvidence: { schema: 1, entries: [{ schema: 1, kind: 'derived-losatp-payload', key, payload: {} }] }
  });
  const prepared = { ...structuredClone(lambdaPrepared.session), legacyArtifacts: evidence('prepared') };
  const withoutLegacy = await finalizeGallerySessionPublication({
    prepared, replayed: structuredClone(lambdaPrepared.session)
  });
  assert.equal(Object.hasOwn(withoutLegacy.session, 'legacyArtifacts'), false);
  const withLegacy = await finalizeGallerySessionPublication({
    prepared, replayed: { ...structuredClone(lambdaPrepared.session), legacyArtifacts: evidence('replayed') }
  });
  assert.deepEqual(withLegacy.session.legacyArtifacts, evidence('replayed'));
}

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
  const scopeless = await loadSession('vibrio-harveyi-group-collinear.gbdraw');
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

// A Session the Web app saves keeps each imported table file as its input
// binding, while its request reads the table that Generate wrote from the
// draft: the Default colors rows that differ from the palette, or the rules
// left after an edit. Publication rebuilds from the request's tables, not
// from the imported files (OV-367).
{
  const source = await loadSession('BGC0000708-BGC0000713.gbdraw');
  const tsvResource = (kind, name, text) => ({ kind, name, type: 'text/tab-separated-values',
    size: Buffer.byteLength(text), lastModified: 0, encoding: 'base64', data: Buffer.from(text, 'utf8').toString('base64') });
  const binding = (resourceId, name) => ({ resourceId, name, type: 'text/tab-separated-values', lastModified: 0 });
  const tableText = (session, ref) => Buffer.from(session.resources[ref.resourceId].data, 'base64').toString('utf8');
  // The Linear File card rows a Save writes beside the table bindings: one per request record, bound to its
  // GenBank resource (OV-1819).
  const linearSeqs = source.renderRequest.records.map((record) => ({
    uid: record.recordKey,
    gb: { resourceId: record.source.resourceId, name: `${record.recordKey}.gbk`, type: 'text/plain', lastModified: 0 },
    gff: null, fasta: null, depth: null, losat_gencode: 1,
    definition: '', record_subtitle: record.presentation.subtitle || '', file_definition: '', file_subtitle: '',
    inferred_definition: record.presentation.label,
    region_record_id: '', region_start: null, region_end: null,
    region_reverse: Boolean(record.presentation.reverseComplement)
  }));
  const generatedDefaultColors = 'CDS\t#d3d3d3\nrRNA\t#71ee7d\ntRNA\t#e8b441\ntmRNA\t#ded44e\n'
    + 'skew_high\t#6dded3\nskew_low\t#ad72e3\npairwise_match_min\t#FFE7E7\npairwise_match_max\t#FF7272\n';
  const savedDefaultColors = structuredClone(source);
  savedDefaultColors.resources['resource-0001'] = source.resources['colors-default-colors-file'];
  savedDefaultColors.resources['colors-default-colors-file'] = tsvResource('colors-default-colors-file',
    'colors-default-colors-file-default-colors.tsv', generatedDefaultColors);
  savedDefaultColors.webFiles.bindings = { schema: 2, linearSeqs, d_color: binding('resource-0001', 'colors-default-colors.tsv') };
  const sourceColors = source.renderRequest.diagramOptions.colors;
  const importedRules = `${tableText(source, sourceColors.colorTable || sourceColors.colorTableFile)}CDS\tgene\tneoU\t#000000\tRemoved after import\n`;
  const savedRules = structuredClone(source);
  savedRules.resources['resource-0002'] = tsvResource('colors-color-table-file', 'resource-0002-specific-colors.tsv', importedRules);
  savedRules.webFiles.bindings = { schema: 2, linearSeqs, t_color: binding('resource-0002', 'specific-colors.tsv') };
  for (const [name, saved] of Object.entries({ savedDefaultColors, savedRules })) {
    const { session, equivalence } = await prepareGallerySessionForPublication(saved);
    assert.deepEqual(equivalence.differences, [], name);
    const { colors } = session.renderRequest.diagramOptions;
    for (const [field, alias] of [['defaultColorsFile', 'defaultColors'], ['colorTable', 'colorTableFile']]) {
      const committed = saved.renderRequest.diagramOptions.colors;
      assert.equal(tableText(session, colors[field] || colors[alias]), tableText(saved, committed[field] || committed[alias]), `${name} ${field}`);
    }
  }
}

// A Session the Web app saves binds its Circular input as File cards. A card's uid keys the request's records and
// names them in the slice's record rows (Record Order). The published Session holds no bindings, so its Load makes
// the cards from the request: they keep the uids the records are keyed by, and the published Session rebuilds its
// own request. Before, Load named the card `record-1` and the publication refused the Session (OV-3022).
{
  const saved = await loadSession('ecoli_o157_sakai_rings.gbdraw');
  const [card] = saved.webFiles.bindings?.circularSources ?? [{ uid: saved.renderRequest.records[0].recordKey }];
  assert.match(card.uid, /^circular-source-/);
  assert.deepEqual(saved.modes.circular.config.adv.multi_record_positions, [{ selector: card.uid, row: 1 }]);
  const { session } = await prepareGallerySessionForPublication(saved);
  assert.equal(Object.hasOwn(session.webFiles, 'bindings'), false);
  const { projection } = await validateGalleryPublicationReadiness(session);
  assert.deepEqual(projection.files.circularSources.map(({ uid }) => uid), [card.uid]);
}

// A Session the Web app saves keeps each File card's own state in
// `webFiles.bindings.linearSeqs`: a record's typed name apart from the name
// inferred from its file. Publication keeps those rows, so Load and Reset
// Settings read the names as after that Save (OV-1701); the imported table
// files stay unpublished (OV-367).
{
  const source = await loadSession('BGC0000708-BGC0000713.gbdraw');
  const saved = structuredClone(source);
  const cardRows = saved.renderRequest.records.map((record, index) => ({
    uid: record.recordKey,
    gb: { resourceId: record.source.resourceId, name: `BGC-${index + 1}.gbk`, type: 'text/plain', lastModified: 0 },
    gff: null, fasta: null, depth: null, losat_gencode: 1,
    // The first record's name is typed over a different inferred name; the others draw their inferred name.
    definition: index === 0 ? record.presentation.label : '', record_subtitle: record.presentation.subtitle || '',
    file_definition: '', file_subtitle: '',
    inferred_definition: index === 0 ? 'Inferred name' : record.presentation.label,
    region_record_id: '', region_start: null, region_end: null,
    region_reverse: Boolean(record.presentation.reverseComplement)
  }));
  saved.webFiles.bindings = { schema: 2, linearSeqs: cardRows,
    d_color: { resourceId: 'record-1-genbank', name: 'imported-colors.tsv', type: 'text/plain', lastModified: 0 } };
  const { session, equivalence } = await prepareGallerySessionForPublication(saved);
  assert.deepEqual(equivalence.differences, []);
  assert.deepEqual(session.webFiles.bindings, { schema: 2, linearSeqs: cardRows });
}

// A Linear card that names the record of a file holding only that record (the Record select's choice) publishes: the
// committed request, as the CLI writes it, has no selector for the record, and the rebuild reads each card's records,
// as Generate does, and writes the same (`cliRecordSelector`, one rule with Circular).
{
  const source = await loadSession('BGC0000708-BGC0000713.gbdraw');
  const saved = structuredClone(source);
  const recordIdOf = (resourceId) => {
    const bytes = Buffer.from(saved.resources[resourceId].data, 'base64');
    const text = (bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes).toString('utf8');
    return /^VERSION\s+(\S+)/m.exec(text)[1];
  };
  const cardRows = saved.renderRequest.records.map((record, index) => ({
    uid: record.recordKey,
    gb: { resourceId: record.source.resourceId, name: `BGC-${index + 1}.gbk`, type: 'text/plain', lastModified: 0 },
    gff: null, fasta: null, depth: null, losat_gencode: 1,
    definition: '', record_subtitle: record.presentation.subtitle || '', file_definition: '', file_subtitle: '',
    inferred_definition: record.presentation.label,
    region_record_id: recordIdOf(record.source.resourceId), region_start: null, region_end: null,
    region_reverse: Boolean(record.presentation.reverseComplement)
  }));
  assert.ok(saved.renderRequest.records.every((record) => record.selector === null));
  saved.webFiles.bindings = { schema: 2, linearSeqs: cardRows };
  const { session, equivalence } = await prepareGallerySessionForPublication(saved);
  assert.deepEqual(equivalence.differences, []);
  assert.deepEqual(session.renderRequest.records.map((record) => record.selector), cardRows.map(() => null));
}

// A Session the CLI writes publishes as the Gallery entry of its command: a
// Circular request of one record is the single drawing the command saves
// (OV-3297) with its one output (OV-267), the resolved
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
      const file = path.join(directory, `${id}.gbdraw`);
      execFileSync('python', ['-m', 'gbdraw.cli', mode, ...args, '-o', path.join(directory, id),
        '-f', mode === 'circular' ? 'interactive_svg' : 'svg',
        '--session_output', file], { env: { ...process.env, PYTHONPATH: process.cwd() }, stdio: 'pipe', timeout: 1_800_000 });
      const source = await loadSessionFile(file);
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
        assert.deepEqual([source.renderRequest.grouping, session.renderRequest.grouping], ['single', 'single'], id);
        assert.equal(session.renderRequest.output.interactiveMetadataPolicy,
          source.renderRequest.output.interactiveMetadataPolicy, id);
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

// A Session that draws a Depth track publishes with Depth shown: the draft split
// keeps Show Depth only when the drawing has a Depth file, and it counts the
// files the projection bound, as Load binds them. The command states no Depth
// style, and the styles it leaves unstated compare as absent (OV-2819).
{
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-gallery-cli-depth-'));
  try {
    const depth = path.join(directory, 'lambda.depth.tsv');
    await writeFile(depth, Array.from({ length: 48502 }, (_, index) => `NC_001416.1\t${index + 1}\t${10 + (index % 50)}\n`).join(''));
    const file = path.join(directory, 'depth.gbdraw');
    execFileSync('python', ['-m', 'gbdraw.cli', 'circular', '--gbk', 'tests/test_inputs/NC_001416.gb',
      '--depth_track', depth, '-o', path.join(directory, 'depth'), '-f', 'interactive_svg', '--session_output', file],
    { env: { ...process.env, PYTHONPATH: process.cwd() }, stdio: 'pipe', timeout: 1_800_000 });
    const source = await loadSessionFile(file);
    assert.equal(source.renderRequest.diagramOptions.depthTracks.length, 1);
    const { session, equivalence } = await prepareGallerySessionForPublication(source);
    assert.deepEqual(equivalence.differences, [], 'CLI depth');
    assert.equal(session.modes.circular.config.form.show_depth, true);
    // A Session with the Web draft (a published or Web-saved one) publishes again with Depth shown.
    const again = await prepareGallerySessionForPublication(session);
    assert.equal(again.session.modes.circular.config.form.show_depth, true);
    const readiness = await validateGalleryPublicationReadiness(again.session);
    assert.deepEqual(readiness.equivalence.differences, [], 'published depth');
    // A style the command states still compares: a rebuild that draws another
    // color is refused (OV-2819 drops only the fields the command left unstated).
    const stated = path.join(directory, 'stated.gbdraw');
    execFileSync('python', ['-m', 'gbdraw.cli', 'circular', '--gbk', 'tests/test_inputs/NC_001416.gb',
      '--depth_track', depth, '--depth_track_color', '#123456', '-o', path.join(directory, 'stated'), '-f', 'interactive_svg',
      '--session_output', stated], { env: { ...process.env, PYTHONPATH: process.cwd() }, stdio: 'pipe', timeout: 1_800_000 });
    const recolored = createGallerySessionPublication({
      assertRequestsEquivalent: assertCanonicalRenderRequestsEquivalent,
      /** @param {Parameters<typeof buildCanonicalRenderRequest>} args */
      buildRequest: (...args) => {
        const built = buildCanonicalRenderRequest(...args);
        built.renderRequest.diagramOptions.depthTracks[0].color = '#654321';
        return built;
      },
      buildRequestState: buildCanonicalRequestState,
      readRequest: readPublishedRenderRequest,
      projectRequest: projectCanonicalSessionRequest,
      resolveComparisonPlan: resolveLinearComparisonPlan,
      planLinearUnits,
      splitDraft
    });
    await assert.rejects(recolored.prepareGallerySessionForPublication(await loadSessionFile(stated)),
      /depthTracks\[0\]\.color/);
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
}

// A CLI Session's wrap length publishes in the Web's request and the drawing's
// form (OV-2455). The equivalence leaves a CLI request's configuration out
// (OV-398), so it could not see the leaf the rebuild dropped.
{
  const { session } = await prepareGallerySessionForPublication(
    await loadReleasedSession('cli-linear-wrapped.v50.gbdraw-session.json.gz'));
  assert.equal(session.renderRequest.diagramOptions.configOverrides['canvas.linear.wrap_length'], 5000);
  assert.equal(session.modes.linear.config.form.linear_wrap_length, 5000);
}

// A CLI command's annotation slot params are text ("0", "false") and the CLI
// leaves an Axis index it was not given unstated; the published request compares
// them as Load reads them (OV-1299). A real difference in either still fails.
{
  const directory = await mkdtemp(path.join(tmpdir(), 'gbdraw-gallery-cli-slots-'));
  try {
    const data = 'gbdraw/web/tutorial-data/tobacco-plastome-regions';
    const file = path.join(directory, 'slots.gbdraw');
    execFileSync('python', ['-m', 'gbdraw.cli', 'circular', '--gbk', `${data}/NC_001879.gbk`,
      '--annotation_table', `${data}/nicotiana-tabacum-regions.tsv`,
      '--circular_track_slot', 'features:features@side=overlay,lane_direction=split',
      '--circular_track_slot', 'names:annotations@set_id=plastome_regions,side=inside,r=0.75,w=12px,'
        + 'padding_px=0,show_labels=false,overflow=compress',
      '--circular_track_slot', 'regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,'
        + 'show_labels=true,padding_px=1,overflow=compress',
      '--circular_track_slot', 'gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08',
      '-o', path.join(directory, 'slots'), '-f', 'interactive_svg', '--session_output', file],
    { env: { ...process.env, PYTHONPATH: process.cwd() }, stdio: 'pipe', timeout: 1_800_000 });
    const source = await loadSessionFile(file);
    const { tracks } = source.renderRequest.diagramOptions;
    assert.deepEqual([tracks.circularTrackSlots[1].params.padding_px, tracks.circularTrackSlots[1].params.show_labels,
      tracks.circularTrackSlots[2].params.padding_px, tracks.circularTrackAxisIndex ?? null], ['0', 'false', '1', null]);
    const { session, equivalence } = await prepareGallerySessionForPublication(source);
    assert.deepEqual(equivalence.differences, [], 'CLI slot params');
    assert.deepEqual(session.renderRequest.diagramOptions.tracks.circularTrackSlots[1].params.padding_px, 0);
    assert.equal(session.renderRequest.diagramOptions.tracks.circularTrackAxisIndex, 0);
    // A rebuild that reads a param or the Axis index differently is a difference.
    for (const [change, at] of [
      [(rebuilt) => { rebuilt.circularTrackSlots[2].params.padding_px = 3; }, /circularTrackSlots\[2\]\.params\.padding_px/],
      [(rebuilt) => { rebuilt.circularTrackAxisIndex = 1; }, /circularTrackAxisIndex/]
    ]) {
      const tampered = createGallerySessionPublication({
        assertRequestsEquivalent: assertCanonicalRenderRequestsEquivalent,
        buildRequest: (input) => {
          const built = buildCanonicalRenderRequest(input);
          change(built.renderRequest.diagramOptions.tracks);
          return built;
        },
        buildRequestState: buildCanonicalRequestState,
        readRequest: readPublishedRenderRequest,
        projectRequest: projectCanonicalSessionRequest,
        resolveComparisonPlan: resolveLinearComparisonPlan,
        planLinearUnits,
        splitDraft
      });
      await assert.rejects(tampered.prepareGallerySessionForPublication(structuredClone(source)),
        (error) => at.test(String(error?.message)) && error.comparison?.equivalent === false);
    }
  } finally {
    await rm(directory, { recursive: true, force: true });
  }
}
