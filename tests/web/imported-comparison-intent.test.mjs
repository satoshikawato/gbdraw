import assert from 'node:assert/strict';

import {
  IMPORTED_COMPARISON_ACTIONS,
  IMPORTED_COMPARISON_DISPOSITIONS,
  classifyImportedComparisonIntent,
  createImportedComparisonIntentState,
  importedComparisonExecution,
  inheritCommittedComparisonIntent,
  resolveImportedComparisonAction,
  restoreImportedComparisonIntent,
  serializeImportedComparisonResolution
} from '../../gbdraw/web/js/services/imported-comparison-intent.js';

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
assert.match(nonExecutableEndpoints.message, /not executable/i);

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
assert.equal(importedComparisonExecution({
  intent: decisionState,
  draftResolution: { valid: true, hasComparisonIntent: false }
}).ok, false);
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
  }]),
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
inheritCommittedComparisonIntent({ candidate, committed });
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

console.log('imported comparison intent tests passed');
