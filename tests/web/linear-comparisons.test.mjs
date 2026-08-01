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
  buildPairwiseLosatJobSpecs,
  createDefaultLinearComparisonPlan,
  createLinearComparisonEdge,
  linearComparisonEdgeKey,
  materializeResolvedEdgesAsSelectedPlan,
  normalizeLinearComparisonPlan,
  reconcileLinearComparisonPlan,
  resolveLinearComparisonPlan,
  validateLinearComparisonEdges
} = await import(pathToFileURL(join(tempRoot, 'linear-comparisons.js')));

const file = (name) => ({ name });
const sequences = [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }, { uid: 'd' }];
const twoRows = [
  { uid: 'a', row: 1 }, { uid: 'b', row: 1 },
  { uid: 'c', row: 2 }, { uid: 'd', row: 2 }
];

assert.deepEqual(createDefaultLinearComparisonPlan(), {
  mode: 'adjacent',
  defaultSource: 'losat',
  edges: []
});
assert.deepEqual(adjacentRowPairs(sequences, twoRows), [['a', 'c'], ['b', 'd']]);
assert.deepEqual(adjacentRowPairs(sequences, twoRows, true), [
  ['a', 'c'], ['a', 'd'], ['b', 'c'], ['b', 'd']
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

const reordered = resolveLinearComparisonPlan({
  plan: selectedPlan,
  sequences: [{ uid: 'b' }, { uid: 'a' }, { uid: 'c' }],
  losatProgram: 'tblastx'
});
assert.equal(reordered.edges[0].edgeKey, linearComparisonEdgeKey('a', 'b'));
assert.equal(reordered.edges[0].queryIndex, 1);
assert.equal(reordered.edges[0].subjectIndex, 0);
assert.equal(reconcileLinearComparisonPlan(selectedPlan, [{ uid: 'b' }, { uid: 'c' }]).edges.length, 1);
assert.equal(reconcileLinearComparisonPlan(selectedPlan, [{ uid: 'c' }]).edges.length, 0);

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

console.log('linear-comparisons tests passed');
