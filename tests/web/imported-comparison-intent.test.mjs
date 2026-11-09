import assert from 'node:assert/strict';

import {
  IMPORTED_COMPARISON_ACTIONS,
  IMPORTED_COMPARISON_DISPOSITIONS,
  classifyImportedComparisonIntent,
  createImportedComparisonIntentState,
  importedComparisonExecution,
  inheritSavedComparison,
  savedComparisonSource,
  resolveImportedComparisonAction,
  restoreImportedComparisonIntent,
  serializeImportedComparisonResolution
} from '../../gbdraw/web/js/services/imported-comparison-intent.js';
import { normalizeCaughtError, normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';

const records = [
  {
    recordKey: 'record-a',
    presentation: { gridRow: 1, gridColumn: null }
  },
  {
    recordKey: 'record-b',
    presentation: { gridRow: 2, gridColumn: null }
  }
];
const resources = {
  comparison: {
    kind: 'nucleotide-blast',
    encoding: 'base64',
    data: 'Ymxhc3Q='
  },
  protein: {
    kind: 'canonical-tsv',
    encoding: 'base64',
    data: 'cHJvdGVpbg=='
  }
};
const request = (comparisons, overrides = {}) => ({
  schema: 6,
  mode: 'linear',
  records,
  diagramOptions: {
    bitscore: 50,
    evalue: 1e-2,
    identity: 0,
    alignmentLength: 0
  },
  comparisons,
  ...overrides
});
const pipelineSettings = {
  collinearityParams: {
    kind: 'lossless',
    parameters: {
      minAnchors: 1,
      maxUnitGap: 0,
      maxDiagonalDrift: 0,
      maxConflicts: 1,
      mergeOrientation: 'either'
    }
  },
  collinearityUnitMode: 'auto',
  collinearityAnchorMode: 'rbh',
  collinearitySearchScope: 'adjacent',
  collinearityColorMode: 'orientation',
  losatpBin: 'losat',
  ncbiBlastpBin: null,
  losatpThreads: null,
  proteinBlastpMaxHits: 5,
  proteinBlastpCandidateLimit: null,
  orthogroupMembershipMode: 'anchor_core_v1',
  orthogroupMemberMaxHits: 5,
  collinearMaxParalogLinksPerOrthogroup: 2,
  alignOrthogroupFeature: null
};
const currentPipelineSettings = { ...pipelineSettings };
delete currentPipelineSettings.alignOrthogroupFeature;

const editable = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'nucleotideBlast',
    resourceId: 'comparison',
    queryRecordIndex: 0,
    subjectRecordIndex: 1
  }]),
  resources
});
assert.equal(editable.disposition, IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE);
assert.equal(editable.hasCommittedComparison, true);

const editablePipeline = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'generatedProteinComparison',
    mode: 'pairwise',
    pairs: [{ queryRecordIndex: 0, subjectRecordIndex: 1 }],
    settings: pipelineSettings
  }]),
  resources
});
assert.equal(
  editablePipeline.disposition,
  IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
);
const currentEditablePipeline = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'generatedProteinComparison',
    mode: 'pairwise',
    pairs: [{ queryRecordIndex: 0, subjectRecordIndex: 1 }],
    settings: currentPipelineSettings
  }], { schema: 8 }),
  resources
});
assert.equal(
  currentEditablePipeline.disposition,
  IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE
);
const currentPipelineWithLegacyAlignment = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'generatedProteinComparison',
    mode: 'pairwise',
    pairs: [{ queryRecordIndex: 0, subjectRecordIndex: 1 }],
    settings: pipelineSettings
  }], { schema: 8 }),
  resources
});
assert.equal(
  currentPipelineWithLegacyAlignment.disposition,
  IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED
);

const preserved = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'precomputedProteinComparison',
    resourceId: 'protein',
    encoding: 'canonicalTsv',
    queryRecordIndex: 0,
    subjectRecordIndex: 1
  }]),
  resources
});
assert.equal(
  preserved.disposition,
  IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY
);
assert.match(preserved.message, /cannot edit/i);

const decision = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'nucleotideBlast',
    resourceId: 'missing',
    queryRecordIndex: 0,
    subjectRecordIndex: 1
  }]),
  resources
});
assert.equal(
  decision.disposition,
  IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED
);
assert.match(decision.message, /missing/i);
const incompletePipeline = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'generatedProteinComparison',
    mode: 'pairwise',
    pairs: [{ queryRecordIndex: 0, subjectRecordIndex: 1 }],
    settings: {}
  }]),
  resources
});
assert.equal(
  incompletePipeline.disposition,
  IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED
);
const nonExecutableEndpoints = classifyImportedComparisonIntent({
  renderRequest: request([{
    kind: 'nucleotideBlast',
    resourceId: 'comparison',
    queryRecordIndex: 0,
    subjectRecordIndex: 2
  }], {
    records: [
      ...records,
      { recordKey: 'record-c', presentation: { gridRow: 3, gridColumn: null } }
    ]
  }),
  resources
});
assert.equal(
  nonExecutableEndpoints.disposition,
  IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED
);
assert.match(nonExecutableEndpoints.message, /not on adjacent rows/i);

const inheritedState = createImportedComparisonIntentState();
restoreImportedComparisonIntent(inheritedState, preserved, {
  action: IMPORTED_COMPARISON_ACTIONS.INHERIT
});
assert.equal(inheritedState.action, IMPORTED_COMPARISON_ACTIONS.INHERIT);
assert.deepEqual(serializeImportedComparisonResolution(inheritedState), {
  action: IMPORTED_COMPARISON_ACTIONS.INHERIT
});
assert.deepEqual(importedComparisonExecution({
  intent: inheritedState,
  draftResolution: { valid: true, hasComparisonIntent: false }
}), { ok: true, mode: 'inherit' });

const decisionState = createImportedComparisonIntentState();
restoreImportedComparisonIntent(decisionState, decision);
// A refused Generate (empty draft), as the app reports it in Linear mode.
const generateFailure = (intent) => {
  const execution = importedComparisonExecution({
    intent,
    draftResolution: { valid: true, hasComparisonIntent: false }
  });
  assert.equal(execution.ok, false);
  return normalizeUserFacingError(execution.message, {
    operation: 'generate', stage: 'request-validation', code: execution.fallbackCode
  });
};
// A missing resource keeps its own diagnostic.
assert.equal(generateFailure(decisionState)?.code, 'COMPARISON_INPUT');
// OV-303: any other reason names the pending choice, not the UNKNOWN fallback.
for (const intent of [preserved, nonExecutableEndpoints]) {
  const error = generateFailure(intent);
  assert.deepEqual([error?.code, error?.actions], ['COMPARISON_CHOICE_REQUIRED', ['edit-comparison']]);
  assert.equal(error?.summary,
    'Choose how to handle the saved comparison in the Comparison panel, then Generate again.');
}
assert.equal(resolveImportedComparisonAction({
  intent: decisionState,
  action: IMPORTED_COMPARISON_ACTIONS.INHERIT,
  draftResolution: { valid: true, hasComparisonIntent: true }
}).ok, false);
assert.equal(resolveImportedComparisonAction({
  intent: decisionState,
  action: IMPORTED_COMPARISON_ACTIONS.REPLACE,
  draftResolution: { valid: false, hasComparisonIntent: true, error: 'Incomplete draft.' }
}).ok, false);
assert.deepEqual(resolveImportedComparisonAction({
  intent: decisionState,
  action: IMPORTED_COMPARISON_ACTIONS.REPLACE,
  draftResolution: { valid: true, hasComparisonIntent: true }
}), { ok: true, action: IMPORTED_COMPARISON_ACTIONS.REPLACE });
assert.deepEqual(importedComparisonExecution({
  intent: decisionState,
  draftResolution: { valid: true, hasComparisonIntent: true }
}), { ok: true, mode: 'draft' });
// OV-301: after Replace, an invalid draft runs and fails there with its own
// diagnostic; an empty draft names the two ways forward.
assert.deepEqual(importedComparisonExecution({
  intent: decisionState,
  draftResolution: { valid: false, hasComparisonIntent: true, error: 'Incomplete draft.' }
}), { ok: true, mode: 'draft' });
const emptyReplacement = generateFailure(decisionState);
assert.deepEqual([emptyReplacement?.code, emptyReplacement?.actions],
  ['COMPARISON_REPLACEMENT_EMPTY', ['edit-comparison']]);
assert.equal(emptyReplacement?.summary,
  'Set up a comparison in the Comparison panel, or choose Clear comparison, then Generate again.');
assert.equal(importedComparisonExecution({ intent: decisionState, draftResolution: null }).ok, false);

restoreImportedComparisonIntent(decisionState, decision);
assert.deepEqual(resolveImportedComparisonAction({
  intent: decisionState,
  action: IMPORTED_COMPARISON_ACTIONS.CLEAR,
  draftResolution: { valid: false, hasComparisonIntent: false }
}), { ok: true, action: IMPORTED_COMPARISON_ACTIONS.CLEAR });
assert.deepEqual(importedComparisonExecution({
  intent: decisionState,
  draftResolution: { valid: false, hasComparisonIntent: false }
}), { ok: true, mode: 'clear' });

const committed = {
  renderRequest: request([{
    kind: 'precomputedProteinComparison',
    resourceId: 'protein',
    encoding: 'canonicalTsv',
    queryRecordIndex: 0,
    subjectRecordIndex: 1
  }], { schema: 9 }),
  resources
};
const candidate = {
  renderRequest: request([], {
    records: [records[1], records[0]],
    diagramOptions: {
      bitscore: 99,
      evalue: 1e-20,
      identity: 90,
      alignmentLength: 200,
      pairwiseMatchStyle: 'curve'
    }
  }),
  resources: {}
};
inheritSavedComparison({ candidate, source: savedComparisonSource(committed) });
assert.deepEqual(candidate.renderRequest.comparisons, [{
  kind: 'precomputedProteinComparison',
  resourceId: 'protein',
  encoding: 'canonicalTsv',
  queryRecordIndex: 1,
  subjectRecordIndex: 0
}]);
assert.equal(candidate.resources.protein, resources.protein);
assert.deepEqual(
  {
    bitscore: candidate.renderRequest.diagramOptions.bitscore,
    evalue: candidate.renderRequest.diagramOptions.evalue,
    identity: candidate.renderRequest.diagramOptions.identity,
    alignmentLength: candidate.renderRequest.diagramOptions.alignmentLength
  },
  { bitscore: 50, evalue: 1e-2, identity: 0, alignmentLength: 0 }
);
assert.equal(candidate.renderRequest.diagramOptions.pairwiseMatchStyle, 'curve');

const circular = classifyImportedComparisonIntent({
  renderRequest: {
    schema: 6,
    mode: 'circular',
    records: [records[0]],
    comparisons: [],
    diagramOptions: {
      conservationBlastFiles: [{ resourceId: 'comparison' }]
    }
  },
  resources
});
assert.equal(circular.disposition, IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE);

// OV-796 (PR-D): a schema-1 request is classified as its reader reads it,
// through the lift, like the schema-2 re-save of the same Gallery entry.
{
  const { readFileSync } = await import('node:fs');
  const { gunzipSync } = await import('node:zlib');
  const session = JSON.parse(gunzipSync(readFileSync(new URL(
    '../fixtures/sessions/lambda_basic_linear.v31-schema1.gbdraw-session.json.gz', import.meta.url
  ))));
  assert.equal(session.renderRequest.schema, 1);
  assert.equal(classifyImportedComparisonIntent({
    renderRequest: session.renderRequest, resources: session.resources
  }).disposition, IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE);
}

// OV-2995: a request's comparison indexes count drawn units, one per window (D-01), as Python reads them; a record
// drawn in two windows draws two units. The classifier reads them through `unitRecordsOf`.
{
  const window = (start, end) => ({ start, end, reverseComplement: false, column: null });
  const windowedRecords = [
    { recordKey: 'record-a', presentation: { gridRow: 1, gridColumn: null }, windows: [window(1, 15000), window(25000, 40000)] },
    { recordKey: 'record-b', presentation: { gridRow: 2, gridColumn: null }, windows: [] }
  ];
  const windowedRequest = (comparisons) => request(comparisons, { schema: 10, records: windowedRecords });
  const unitPipeline = (pairs) => windowedRequest([{
    kind: 'generatedProteinComparison', mode: 'pairwise', settings: currentPipelineSettings,
    pairs: pairs.map(([queryRecordIndex, subjectRecordIndex]) => ({ queryRecordIndex, subjectRecordIndex }))
  }]);
  const classify = (renderRequest) => classifyImportedComparisonIntent({ renderRequest, resources });
  // Window 1 and window 2 of record-a against record-b: unit 2 is record-b, not a third record. The controls hold
  // record pairs, so the comparison runs as saved through Inherit.
  const intent = classify(unitPipeline([[0, 2], [1, 2]]));
  assert.equal(intent.disposition, IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY, JSON.stringify(intent));
  assert.equal(intent.hasCommittedComparison, true);
  // The preserved comparison holds what Inherit copies (OV-3001): Generate reads it from the drawing, not from a Result.
  assert.deepEqual(intent.source, savedComparisonSource({ renderRequest: unitPipeline([[0, 2], [1, 2]]), resources }));
  assert.deepEqual(intent.source.records.map(({ recordKey }) => recordKey), ['record-a', 'record-b']);
  assert.equal(classify(unitPipeline([[0, -1]])).source, null, 'a comparison that is not preserved holds no source');
  // An index that is no unit index, or a unit against itself, identifies no two records.
  assert.equal(classify(unitPipeline([[0, -1]])).disposition, IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED);
  assert.equal(classify(unitPipeline([[2, 2]])).disposition, IMPORTED_COMPARISON_DISPOSITIONS.DECISION_REQUIRED);
  // The CLI's LOSAT comparison of the same drawing: one protein table per unit pair (gbdraw linear --region twice).
  const cliRequest = windowedRequest([
    { kind: 'precomputedProteinComparison', resourceId: 'protein', encoding: 'canonicalTsv', queryRecordIndex: 0, subjectRecordIndex: 2 },
    { kind: 'precomputedProteinComparison', resourceId: 'protein', encoding: 'canonicalTsv', queryRecordIndex: 1, subjectRecordIndex: 2 },
    { kind: 'generatedProteinComparison', mode: 'none', settings: currentPipelineSettings, pairs: [] }
  ]);
  assert.equal(classify(cliRequest).disposition, IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY);
  // Inherit carries unit indexes into a request that draws the same units (the same records and windows), as saved.
  const committedWindowed = classify(cliRequest).source;
  const sameUnits = { renderRequest: request([], { schema: 10, records: structuredClone(windowedRecords) }), resources: {} };
  inheritSavedComparison({ candidate: sameUnits, source: committedWindowed });
  assert.deepEqual(sameUnits.renderRequest.comparisons.map((comparison) => [comparison.queryRecordIndex, comparison.subjectRecordIndex]),
    [[0, 2], [1, 2], [undefined, undefined]]);
  // A draft whose windows changed draws other units: the saved indexes name none of them.
  const otherUnits = { renderRequest: request([], { schema: 10, records: [{ ...windowedRecords[0], windows: [window(1, 15000)] }, windowedRecords[1]] }), resources: {} };
  // Generate shows the choice in the Comparison panel and why, not the UNKNOWN fallback.
  /** @param {() => unknown} inherit */
  const shownFailure = (inherit) => {
    try {
      inherit();
    } catch (error) {
      return normalizeCaughtError(error, { operation: 'generate', stage: 'request-validation' });
    }
    assert.fail('Inherit did not refuse');
  };
  const changed = shownFailure(() => inheritSavedComparison({ candidate: otherUnits, source: committedWindowed }));
  assert.equal(changed.code, 'COMPARISON_CHOICE_REQUIRED');
  assert.match(changed.summary, /The drawn records or windows changed since the comparison was saved\. Choose Replace or Clear\./);
  assert.deepEqual(changed.actions, ['edit-comparison']);
  // A drawing that holds no saved comparison gets the same choice.
  const missing = shownFailure(() => inheritSavedComparison({ candidate: structuredClone(sameUnits), source: null }));
  assert.equal(missing.code, 'COMPARISON_CHOICE_REQUIRED');
  assert.match(missing.summary, /This drawing holds no saved comparison to inherit\./);
}

console.log('imported comparison intent tests passed');
