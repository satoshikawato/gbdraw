import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-linear-comparisons-'));
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'linear-comparisons.js'),
  join(tempRoot, 'linear-comparisons.js')
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const {
  LINEAR_COMPARISON_MODES,
  adjacentRowPairs,
  buildLinearComparisonTimeline,
  buildPairwiseLosatJobSpecs,
  createDefaultLinearComparisonPlan,
  createLinearComparisonEdge,
  linearComparisonEdgeKey,
  materializeResolvedEdgesAsSelectedPlan,
  normalizeLinearComparisonPlan,
  plainTextLinearRecordLabel,
  reconcileLinearComparisonPlan,
  resolveLinearComparisonPlan,
  validateLinearComparisonEdges
} = await import(pathToFileURL(join(tempRoot, 'linear-comparisons.js')));

assert.equal(plainTextLinearRecordLabel('Alpha <em>beta</em>'), 'Alpha beta');
assert.equal(plainTextLinearRecordLabel('Alpha <scr<script>ipt>'), 'Alpha');
assert.equal(plainTextLinearRecordLabel('Alpha <script'), 'Alpha script');
assert.equal(plainTextLinearRecordLabel('<b></b>'), 'Record');

const file = (name) => ({ name });
const sequences = [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }, { uid: 'd' }];
const twoRows = [
  { uid: 'a', row: 1 }, { uid: 'b', row: 1 },
  { uid: 'c', row: 2 }, { uid: 'd', row: 2 }
];
const sparseRows = [
  { uid: 'a', row: 2 }, { uid: 'b', row: 7 },
  { uid: 'c', row: 20 }, { uid: 'd', row: 20 }
];

assert.deepEqual(createDefaultLinearComparisonPlan(), {
  mode: 'none',
  defaultSource: 'losat',
  edges: []
});
assert.deepEqual(adjacentRowPairs(sequences, twoRows), [['a', 'c'], ['b', 'd']]);
assert.deepEqual(adjacentRowPairs(sequences, twoRows, true), [
  ['a', 'c'], ['a', 'd'], ['b', 'c'], ['b', 'd']
]);
assert.deepEqual(adjacentRowPairs(sequences, sparseRows), [['a', 'b'], ['b', 'c']]);
assert.deepEqual(adjacentRowPairs(sequences, sparseRows, true), [
  ['a', 'b'], ['b', 'c'], ['b', 'd']
]);

const retained = file('retained.tsv');
const normalized = normalizeLinearComparisonPlan({
  mode: 'unexpected',
  defaultSource: 'unexpected',
  edges: [{
    id: 'draft',
    queryUid: ' a ',
    subjectUid: ' b ',
    file: retained,
    losatFilename: 'retained.losat.tsv'
  }]
});
assert.equal(normalized.mode, 'adjacent');
assert.equal(normalized.defaultSource, 'losat');
assert.equal(
  normalizeLinearComparisonPlan({ defaultSource: 'upload', edges: [] }).mode,
  'adjacent',
  'a missing persisted mode keeps the legacy normalizer fallback'
);
assert.deepEqual(normalized.edges[0], {
  id: 'draft',
  queryUid: 'a',
  subjectUid: 'b',
  included: false,
  fileActive: false,
  losatFilenameActive: false,
  source: 'upload',
  file: retained,
  losatFilename: 'retained.losat.tsv'
});

const noneResolution = resolveLinearComparisonPlan({
  plan: {
    mode: 'none',
    defaultSource: 'upload',
    edges: [{
      ...normalized.edges[0],
      included: true
    }]
  },
  sequences: sequences.slice(0, 2)
});
assert.deepEqual(noneResolution.edges, []);
assert.equal(noneResolution.valid, true);
assert.equal(noneResolution.hasComparisonIntent, false);
assert.equal(noneResolution.hasLosatIntent, false);
assert.equal(noneResolution.hasUploadIntent, false);

assert.deepEqual(resolveLinearComparisonPlan({
  plan: createDefaultLinearComparisonPlan(),
  sequences: sequences.slice(0, 1)
}).edges, []);

const sparsePlan = {
  mode: 'adjacent',
  defaultSource: 'upload',
  edges: [
    createLinearComparisonEdge({
      id: 'a-b', queryUid: 'a', subjectUid: 'b', source: 'upload',
      file: file('a-b.tsv'), fileActive: true
    }),
    createLinearComparisonEdge({
      id: 'b-c', queryUid: 'b', subjectUid: 'c', source: 'upload',
      file: file('b-c-retained.tsv'), fileActive: false
    })
  ]
};
const sparseResolution = resolveLinearComparisonPlan({
  plan: sparsePlan,
  sequences: sequences.slice(0, 3)
});
assert.deepEqual(sparseResolution.edges.map((edge) => ({
  edgeKey: edge.edgeKey,
  ordinal: edge.ordinal,
  file: edge.file.name
})), [{ edgeKey: 'a->b', ordinal: 0, file: 'a-b.tsv' }]);
assert.equal(sparseResolution.valid, true);
assert.equal(sparseResolution.hasUploadIntent, true);

const independentActivationPlan = {
  mode: 'adjacent',
  defaultSource: 'losat',
  edges: [{
    id: 'activation', queryUid: 'a', subjectUid: 'b', included: false,
    source: 'upload', file: retained, fileActive: false,
    losatFilename: 'custom.tsv', losatFilenameActive: true
  }]
};
const activeName = resolveLinearComparisonPlan({
  plan: independentActivationPlan,
  sequences: sequences.slice(0, 2)
});
assert.equal(activeName.edges.length, 1);
assert.equal(activeName.edges[0].source, 'losat');
assert.equal(activeName.edges[0].fileActive, false);
assert.equal(activeName.edges[0].losatFilenameActive, true);
independentActivationPlan.defaultSource = 'upload';
assert.equal(resolveLinearComparisonPlan({
  plan: independentActivationPlan,
  sequences: sequences.slice(0, 2)
}).edges.length, 0);
independentActivationPlan.edges[0].fileActive = true;
const activeFile = resolveLinearComparisonPlan({
  plan: independentActivationPlan,
  sequences: sequences.slice(0, 2)
});
assert.equal(activeFile.edges.length, 1);
assert.equal(activeFile.edges[0].file, retained);
assert.equal(activeFile.edges[0].losatFilenameActive, true);

const selectedPlan = {
  mode: 'selected',
  defaultSource: 'losat',
  edges: [
    createLinearComparisonEdge({ id: 'selected', queryUid: 'a', subjectUid: 'b', source: 'losat' }),
    createLinearComparisonEdge({
      id: 'omitted', queryUid: 'b', subjectUid: 'c', source: 'upload',
      included: false, file: file('omitted.tsv'), fileActive: true
    })
  ]
};
const selectedResolution = resolveLinearComparisonPlan({
  plan: selectedPlan,
  sequences: sequences.slice(0, 3),
  losatProgram: 'blastn'
});
assert.deepEqual(selectedResolution.edges.map((edge) => edge.id), ['selected']);
assert.equal(selectedResolution.hasLosatIntent, true);
assert.equal(selectedResolution.hasUploadIntent, false);

const materialized = materializeResolvedEdgesAsSelectedPlan(sparsePlan, sparseResolution);
assert.equal(materialized.mode, LINEAR_COMPARISON_MODES.SELECTED);
assert.equal(materialized.edges.find((edge) => edge.id === 'a-b').included, true);
assert.equal(materialized.edges.find((edge) => edge.id === 'b-c').included, false);
assert.equal(materialized.edges.find((edge) => edge.id === 'b-c').file.name, 'b-c-retained.tsv');

const adjacentDuplicatePlan = {
  mode: 'adjacent',
  defaultSource: 'losat',
  edges: [
    createLinearComparisonEdge({
      id: 'adjacent-owner', queryUid: 'a', subjectUid: 'b', source: 'upload', included: false
    }),
    createLinearComparisonEdge({
      id: 'adjacent-duplicate', queryUid: 'a', subjectUid: 'b', source: 'losat', included: true
    })
  ]
};
const adjacentDuplicateResolution = resolveLinearComparisonPlan({
  plan: adjacentDuplicatePlan,
  sequences: sequences.slice(0, 2)
});
const materializedDuplicatePlan = materializeResolvedEdgesAsSelectedPlan(
  adjacentDuplicatePlan,
  adjacentDuplicateResolution
);
assert.deepEqual(materializedDuplicatePlan.edges.map(({ id, included }) => ({ id, included })), [
  { id: 'adjacent-owner', included: true },
  { id: 'adjacent-duplicate', included: false }
]);

const reordered = resolveLinearComparisonPlan({
  plan: selectedPlan,
  sequences: [{ uid: 'b' }, { uid: 'a' }, { uid: 'c' }],
  losatProgram: 'tblastx'
});
assert.equal(reordered.edges[0].edgeKey, linearComparisonEdgeKey('a', 'b'));
assert.equal(reordered.edges[0].queryIndex, 1);
assert.equal(reordered.edges[0].subjectIndex, 0);

const repairableFile = file('repairable.tsv');
const repairablePlan = {
  mode: 'selected',
  defaultSource: 'upload',
  edges: [createLinearComparisonEdge({
    id: 'repairable', queryUid: 'a', subjectUid: 'b', source: 'upload',
    file: repairableFile, fileActive: true
  })]
};
const missingUidPlan = reconcileLinearComparisonPlan(repairablePlan, [{ uid: 'b' }, { uid: 'c' }]);
assert.equal(missingUidPlan.edges.length, 1);
assert.equal(missingUidPlan.edges[0].id, 'repairable');
assert.equal(missingUidPlan.edges[0].file, repairableFile);
const missingUidResolution = resolveLinearComparisonPlan({
  plan: missingUidPlan,
  sequences: [{ uid: 'b' }, { uid: 'c' }]
});
assert.equal(missingUidResolution.valid, false);
assert(missingUidResolution.errors.some((issue) => issue.code === 'missing-uid'));
assert.throws(
  () => buildPairwiseLosatJobSpecs({ resolution: missingUidResolution }),
  /no longer available/
);
missingUidPlan.edges[0].queryUid = 'c';
const repairedPlan = reconcileLinearComparisonPlan(missingUidPlan, [{ uid: 'b' }, { uid: 'c' }]);
const repairedResolution = resolveLinearComparisonPlan({
  plan: repairedPlan,
  sequences: [{ uid: 'b' }, { uid: 'c' }]
});
assert.equal(repairedResolution.valid, true);
assert.equal(repairedResolution.edges[0].id, 'repairable');
assert.equal(repairedResolution.edges[0].edgeKey, 'c->b');
assert.equal(repairedResolution.edges[0].file, repairableFile);

const issueCodes = (edges, seqs = sequences, layout = []) => (
  validateLinearComparisonEdges({ edges, sequences: seqs, layout }).map((issue) => issue.code)
);
assert(issueCodes([
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'b', source: 'losat' }),
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'b', source: 'losat' })
]).includes('duplicate'));
assert(issueCodes([
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'a', source: 'losat' })
]).includes('self'));
assert(issueCodes([
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'missing', source: 'losat' })
]).includes('missing-uid'));
assert(issueCodes([
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'b', source: 'losat' })
], sequences, twoRows).includes('same-row'));
assert(!issueCodes([
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'b', source: 'losat' })
], sequences, sparseRows).includes('non-adjacent'));
assert(issueCodes([
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'c', source: 'losat' })
], sequences, sparseRows).includes('non-adjacent'));
assert(issueCodes([
  createLinearComparisonEdge({ queryUid: 'a', subjectUid: 'c', source: 'losat' })
], sequences, [
  { uid: 'a', row: 1 }, { uid: 'b', row: 2 },
  { uid: 'c', row: 3 }, { uid: 'd', row: 4 }
]).includes('non-adjacent'));
assert(issueCodes([
  createLinearComparisonEdge({
    queryUid: 'a', subjectUid: 'b', source: 'upload',
    file: retained, fileActive: false
  })
]).includes('missing-upload'));

const invalidSelectedUpload = resolveLinearComparisonPlan({
  plan: {
    mode: 'selected', defaultSource: 'upload',
    edges: [createLinearComparisonEdge({
      queryUid: 'a', subjectUid: 'b', source: 'upload', file: retained, fileActive: false
    })]
  },
  sequences: sequences.slice(0, 2)
});
assert.equal(invalidSelectedUpload.valid, false);
assert(invalidSelectedUpload.errors.some((issue) => issue.code === 'missing-upload'));

const selectedNonPairwise = resolveLinearComparisonPlan({
  plan: selectedPlan,
  sequences: sequences.slice(0, 3),
  losatProgram: 'blastp',
  blastpMode: 'orthogroup'
});
assert.equal(selectedNonPairwise.valid, false);
assert(selectedNonPairwise.errors.some((issue) => issue.code === 'selected-losat-requires-pairwise'));
const selectedUploadOnly = resolveLinearComparisonPlan({
  plan: {
    mode: 'selected', defaultSource: 'upload',
    edges: [createLinearComparisonEdge({
      queryUid: 'a', subjectUid: 'b', source: 'upload',
      file: retained, fileActive: true
    })]
  },
  sequences: sequences.slice(0, 2),
  losatProgram: 'blastp',
  blastpMode: 'collinear'
});
assert.equal(selectedUploadOnly.valid, true);

for (const program of ['blastn', 'tblastx', 'blastp']) {
  const resolution = resolveLinearComparisonPlan({
    plan: selectedPlan,
    sequences: sequences.slice(0, 3),
    losatProgram: program,
    blastpMode: 'pairwise'
  });
  assert.deepEqual(buildPairwiseLosatJobSpecs({
    resolution,
    program,
    blastpMode: 'pairwise'
  }), [{
    edgeKey: 'a->b',
    ordinal: 0,
    queryUid: 'a',
    subjectUid: 'b',
    queryIndex: 0,
    subjectIndex: 1,
    program
  }]);
}

const timelinePairs = (timeline) => timeline.rows.flatMap((row) => row.boundaryAfter?.pairs || []);
const adjacentLosatPlan = { mode: 'adjacent', defaultSource: 'losat', edges: [] };
const defaultTwoTimeline = buildLinearComparisonTimeline({
  sequences: sequences.slice(0, 2),
  plan: adjacentLosatPlan,
  resolution: resolveLinearComparisonPlan({
    plan: adjacentLosatPlan,
    sequences: sequences.slice(0, 2)
  })
});
assert.deepEqual(defaultTwoTimeline.rows.map((row) => ({
  row: row.row,
  records: row.records.map((record) => record.uid),
  boundary: row.boundaryAfter?.pairs.map((pair) => pair.edgeKey) || null
})), [
  { row: 1, records: ['a'], boundary: ['a->b'] },
  { row: 2, records: ['b'], boundary: null }
]);

const adjacentUploadPlan = { mode: 'adjacent', defaultSource: 'upload', edges: [] };
const adjacentUploadTimeline = buildLinearComparisonTimeline({
  sequences: sequences.slice(0, 2),
  plan: adjacentUploadPlan,
  resolution: resolveLinearComparisonPlan({
    plan: adjacentUploadPlan,
    sequences: sequences.slice(0, 2)
  })
});
assert.equal(timelinePairs(adjacentUploadTimeline)[0].source, 'upload');
assert.equal(timelinePairs(adjacentUploadTimeline)[0].active, false);

const fiveSequences = [...sequences, { uid: 'e' }];
const defaultFivePlan = adjacentLosatPlan;
const defaultFiveTimeline = buildLinearComparisonTimeline({
  sequences: fiveSequences,
  plan: defaultFivePlan,
  resolution: resolveLinearComparisonPlan({ plan: defaultFivePlan, sequences: fiveSequences })
});
assert.deepEqual(
  defaultFiveTimeline.rows.flatMap((row) => row.boundaryAfter?.pairs.map((pair) => pair.edgeKey) || []),
  ['a->b', 'b->c', 'c->d', 'd->e']
);

const zippedTimeline = buildLinearComparisonTimeline({
  sequences,
  layout: twoRows,
  plan: adjacentLosatPlan,
  resolution: resolveLinearComparisonPlan({
    plan: adjacentLosatPlan,
    sequences,
    layout: twoRows
  })
});
assert.equal(zippedTimeline.rows.length, 2);
assert.deepEqual(zippedTimeline.rows[0].records.map((record) => record.uid), ['a', 'b']);
assert.deepEqual(zippedTimeline.rows[0].boundaryAfter.pairs.map((pair) => pair.edgeKey), ['a->c', 'b->d']);

const sparseTimeline = buildLinearComparisonTimeline({
  sequences,
  layout: sparseRows,
  plan: adjacentLosatPlan,
  resolution: resolveLinearComparisonPlan({
    plan: adjacentLosatPlan,
    sequences,
    layout: sparseRows
  })
});
assert.deepEqual(sparseTimeline.rows.map((row) => ({
  row: row.row,
  lowerRow: row.boundaryAfter?.lowerRow || null,
  pairs: row.boundaryAfter?.pairs.map((pair) => pair.edgeKey) || []
})), [
  { row: 2, lowerRow: 7, pairs: ['a->b'] },
  { row: 7, lowerRow: 20, pairs: ['b->c'] },
  { row: 20, lowerRow: null, pairs: [] }
]);

const crossProductPlan = {
  mode: 'selected',
  defaultSource: 'losat',
  edges: [
    createLinearComparisonEdge({ id: 'a-c', queryUid: 'a', subjectUid: 'c', source: 'losat' }),
    createLinearComparisonEdge({ id: 'a-d', queryUid: 'a', subjectUid: 'd', source: 'losat' }),
    createLinearComparisonEdge({ id: 'b-c', queryUid: 'b', subjectUid: 'c', source: 'losat' }),
    createLinearComparisonEdge({ id: 'b-d', queryUid: 'b', subjectUid: 'd', source: 'losat' })
  ]
};
const crossProductTimeline = buildLinearComparisonTimeline({
  sequences,
  layout: twoRows,
  plan: crossProductPlan,
  resolution: resolveLinearComparisonPlan({ plan: crossProductPlan, sequences, layout: twoRows })
});
assert.deepEqual(
  crossProductTimeline.rows[0].boundaryAfter.pairs.map((pair) => pair.edgeKey),
  ['a->c', 'b->d', 'a->d', 'b->c']
);
assert.equal(new Set(timelinePairs(crossProductTimeline).map((pair) => pair.edgeKey)).size, 4);

const reversePlan = {
  mode: 'selected',
  defaultSource: 'losat',
  edges: [createLinearComparisonEdge({
    id: 'reverse', queryUid: 'd', subjectUid: 'b', source: 'losat'
  })]
};
const reverseTimeline = buildLinearComparisonTimeline({
  sequences,
  layout: twoRows,
  plan: reversePlan,
  resolution: resolveLinearComparisonPlan({ plan: reversePlan, sequences, layout: twoRows })
});
const reversePair = timelinePairs(reverseTimeline).find((pair) => pair.edgeKey === 'd->b');
assert.equal(reversePair.queryUid, 'd');
assert.equal(reversePair.subjectUid, 'b');
assert.equal(reversePair.queryIndex, 3);
assert.equal(reversePair.subjectIndex, 1);

const attachedFile = file('attached.tsv');
const attachedPlan = {
  mode: 'selected',
  defaultSource: 'upload',
  edges: [createLinearComparisonEdge({
    id: 'attached', queryUid: 'a', subjectUid: 'b', source: 'upload',
    file: attachedFile, fileActive: true,
    losatFilename: 'attached.raw.tsv', losatFilenameActive: true
  })]
};
const reorderedSequences = [{ uid: 'b' }, { uid: 'a' }];
const reorderedResolution = resolveLinearComparisonPlan({
  plan: attachedPlan,
  sequences: reorderedSequences
});
const reorderedTimeline = buildLinearComparisonTimeline({
  sequences: reorderedSequences,
  plan: attachedPlan,
  resolution: reorderedResolution
});
const attachedPair = timelinePairs(reorderedTimeline).find((pair) => pair.edgeKey === 'a->b');
assert.equal(attachedPair.queryIndex, 1);
assert.equal(attachedPair.subjectIndex, 0);
assert.equal(attachedPair.draft.id, 'attached');
assert.equal(attachedPair.draft.file, attachedFile);
assert.equal(attachedPair.draft.losatFilename, 'attached.raw.tsv');

const dormantPlan = {
  mode: 'none',
  defaultSource: 'losat',
  edges: [createLinearComparisonEdge({
    id: 'dormant', queryUid: 'a', subjectUid: 'b', source: 'losat',
    included: false,
    file: retained,
    fileActive: false,
    losatFilename: 'retained.raw.tsv',
    losatFilenameActive: false
  })]
};
const dormantTimeline = buildLinearComparisonTimeline({
  sequences: sequences.slice(0, 2),
  plan: dormantPlan,
  resolution: resolveLinearComparisonPlan({ plan: dormantPlan, sequences: sequences.slice(0, 2) })
});
const dormantPair = timelinePairs(dormantTimeline)[0];
assert.equal(dormantPair.source, 'none');
assert.equal(dormantPair.active, false);
assert.equal(dormantPair.draft.losatFilename, 'retained.raw.tsv');
assert.equal(dormantPair.draft.losatFilenameActive, false);
assert.equal(dormantPair.draft.file, retained);
assert.equal(dormantPair.draft.fileActive, false);

const duplicateOwnershipPlan = {
  mode: 'selected',
  defaultSource: 'losat',
  edges: [
    createLinearComparisonEdge({
      id: 'primary-inactive', queryUid: 'a', subjectUid: 'b', source: 'upload', included: false
    }),
    createLinearComparisonEdge({
      id: 'duplicate-active', queryUid: 'a', subjectUid: 'b', source: 'losat', included: true
    })
  ]
};
const duplicateOwnershipTimeline = buildLinearComparisonTimeline({
  sequences: sequences.slice(0, 2),
  plan: duplicateOwnershipPlan,
  resolution: resolveLinearComparisonPlan({
    plan: duplicateOwnershipPlan,
    sequences: sequences.slice(0, 2)
  })
});
const primaryPair = timelinePairs(duplicateOwnershipTimeline)[0];
assert.equal(primaryPair.edgeId, 'duplicate-active');
assert.equal(primaryPair.source, 'losat');
assert.equal(primaryPair.active, true);
assert.equal(primaryPair.draft.id, primaryPair.edgeId);
assert.equal(primaryPair.resolved.id, primaryPair.edgeId);
assert.equal(primaryPair.source, primaryPair.draft.source);
assert.equal(primaryPair.source, primaryPair.resolved.source);
assert.deepEqual(
  duplicateOwnershipTimeline.unplacedDrafts.map(({ draft }) => draft.id),
  ['primary-inactive']
);

const invalidLayout = [
  { uid: 'a', row: 1 }, { uid: 'b', row: 1 },
  { uid: 'c', row: 2 }, { uid: 'd', row: 3 }
];
const invalidDraftPlan = {
  mode: 'selected',
  defaultSource: 'losat',
  edges: [
    createLinearComparisonEdge({ id: 'placed', queryUid: 'a', subjectUid: 'c', source: 'losat' }),
    createLinearComparisonEdge({ id: 'duplicate', queryUid: 'a', subjectUid: 'c', source: 'losat' }),
    createLinearComparisonEdge({ id: 'missing', queryUid: 'a', subjectUid: 'missing', source: 'losat' }),
    createLinearComparisonEdge({ id: 'same-row', queryUid: 'a', subjectUid: 'b', source: 'losat' }),
    createLinearComparisonEdge({ id: 'non-adjacent', queryUid: 'a', subjectUid: 'd', source: 'losat' })
  ]
};
const invalidDraftTimeline = buildLinearComparisonTimeline({
  sequences,
  layout: invalidLayout,
  plan: invalidDraftPlan,
  resolution: resolveLinearComparisonPlan({
    plan: invalidDraftPlan,
    sequences,
    layout: invalidLayout
  })
});
assert.deepEqual(
  invalidDraftTimeline.unplacedDrafts.map(({ draft }) => draft.id),
  ['duplicate', 'missing', 'same-row', 'non-adjacent']
);
assert.deepEqual(
  invalidDraftTimeline.unplacedDrafts.map(({ issues }) => issues[0].code),
  ['duplicate', 'missing-uid', 'same-row', 'non-adjacent']
);
assert.equal(timelinePairs(invalidDraftTimeline).filter((pair) => pair.edgeKey === 'a->c').length, 1);

const oneRecordNoneTimeline = buildLinearComparisonTimeline({
  sequences: sequences.slice(0, 1),
  plan: { mode: 'none', defaultSource: 'losat', edges: [] },
  resolution: resolveLinearComparisonPlan({
    plan: { mode: 'none', defaultSource: 'losat', edges: [] },
    sequences: sequences.slice(0, 1)
  })
});
assert.equal(oneRecordNoneTimeline.rows.length, 1);
assert.equal(oneRecordNoneTimeline.rows[0].boundaryAfter, null);
assert.deepEqual(timelinePairs(oneRecordNoneTimeline), []);

console.log('linear-comparisons tests passed');
