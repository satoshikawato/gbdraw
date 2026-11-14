// Web Load keeps the comparison kind a CLI Session stored.
// OV-3573: a request of schema 5 to 9 stored a CLI comparisons-table BLAST table as a protein table
// (`precomputedProteinComparison` of the BLAST columns alone); Load reads it as the BLAST table it is.
// OV-3566: a LOSATP Pairwise search that compared no pair (every record in one row) is stored as the Pairwise
// pipeline without pairs; the controls hold it, and a redraw of it searches nothing.
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { COMPARISON_COLUMNS } from '../../gbdraw/web/js/contracts/schema.generated.js';
import {
  readStoredBlastTables, requestRunsProteinSearch
} from '../../gbdraw/web/js/services/session-request.js';
import {
  IMPORTED_COMPARISON_DISPOSITIONS, classifyImportedComparisonIntent
} from '../../gbdraw/web/js/services/imported-comparison-intent.js';

const base64 = (text) => Buffer.from(text, 'utf8').toString('base64');
const text = (data) => Buffer.from(data, 'base64').toString('utf8');
const BLAST_ROWS = 'rec0\trec1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n';
const table = (header, rows = BLAST_ROWS) => ({
  kind: 'canonical-tsv', name: 'comparison-explicit-1.tsv', type: 'text/tab-separated-values',
  size: 0, lastModified: 0, encoding: 'base64', data: base64(`${header}\n${rows}`)
});
const comparison = (resourceId, queryRecordIndex = 0, subjectRecordIndex = 1) => ({
  kind: 'precomputedProteinComparison', resourceId, encoding: 'canonicalTsv', queryRecordIndex, subjectRecordIndex
});

test('a schema 5-9 protein table of the BLAST columns alone loads as a BLAST table; any other stays as stored', () => {
  const marker = { kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: {} };
  const resources = {
    blast: table(COMPARISON_COLUMNS.join('\t')),
    crlf: table(`${COMPARISON_COLUMNS.join('\t')}\r`),
    protein: table([...COMPARISON_COLUMNS, 'query_protein_id', 'subject_protein_id'].join('\t'))
  };
  const stored = { schema: 9, mode: 'linear', comparisons: [comparison('blast'), comparison('protein', 1, 2), marker] };
  const crlf = { schema: 5, mode: 'linear', comparisons: [comparison('crlf')] };
  const current = { schema: 10, mode: 'linear', comparisons: [comparison('blast')] };
  const read = readStoredBlastTables([stored, crlf, current], resources);
  assert.deepEqual(read.requests[0].comparisons, [
    { kind: 'nucleotideBlast', resourceId: 'blast', queryRecordIndex: 0, subjectRecordIndex: 1 },
    stored.comparisons[1],
    marker
  ]);
  assert.equal(read.requests[0].comparisons[1], stored.comparisons[1]);
  assert.equal(read.requests[1].comparisons[0].kind, 'nucleotideBlast');
  // The request of the current schema is never rewritten.
  assert.equal(read.requests[2], current);
  for (const resourceId of ['blast', 'crlf']) {
    const resource = read.resources[resourceId];
    assert.equal(resource.kind, 'nucleotide-blast');
    assert.equal(text(resource.data), BLAST_ROWS);
    assert.equal(resource.size, Buffer.byteLength(BLAST_ROWS));
    assert.equal(resource.name, 'comparison-explicit-1.tsv');
  }
  assert.equal(read.resources.protein, resources.protein);
  // The file's own values are not edited.
  assert.equal(resources.blast.kind, 'canonical-tsv');
  assert.equal(stored.comparisons[0].kind, 'precomputedProteinComparison');

  const unchanged = readStoredBlastTables([current, { schema: 9, mode: 'linear', comparisons: [comparison('protein')] }],
    resources);
  assert.equal(unchanged.resources, resources);
  assert.equal(unchanged.requests[1].comparisons[0].kind, 'precomputedProteinComparison');
});

const pipelineSettings = {
  collinearityParams: { kind: 'lossless', parameters: {
    minAnchors: 1, maxUnitGap: 0, maxDiagonalDrift: 0, maxConflicts: 1, mergeOrientation: 'strand'
  } },
  collinearityUnitMode: 'auto', collinearityAnchorMode: 'rbh', collinearitySearchScope: 'adjacent',
  collinearityColorMode: 'orientation', losatpBin: 'losat', ncbiBlastpBin: null, losatpThreads: null,
  proteinBlastpMaxHits: 5, proteinBlastpCandidateLimit: null, orthogroupMembershipMode: 'anchor_core_v1',
  orthogroupMemberMaxHits: null, collinearInferOrthogroups: false, collinearMaxParalogLinksPerOrthogroup: 2
};
const pairwiseRequest = ({ pairs = [], positions = null, gridRows = [null, null] } = {}) => ({
  schema: 10, mode: 'linear',
  records: gridRows.map((gridRow, index) => ({ recordKey: `record-${index + 1}`, presentation: { gridRow } })),
  layout: { recordGapPx: 24, multiRecordPositions: positions },
  comparisons: [{ kind: 'generatedProteinComparison', mode: 'pairwise', pairs, settings: pipelineSettings }]
});

test('a Pairwise pipeline without pairs is the controls\' LOSATP Pairwise; over one row it searches nothing', () => {
  const oneRow = pairwiseRequest({ positions: ['#1@1', '#2@1'] });
  assert.equal(classifyImportedComparisonIntent({ renderRequest: oneRow, resources: {} }).disposition,
    IMPORTED_COMPARISON_DISPOSITIONS.EDITABLE);
  assert.equal(requestRunsProteinSearch(oneRow), false);
  // A record drawn in several windows: the search's pairs join drawn units, which only Inherit runs (OV-3003).
  const windowed = pairwiseRequest();
  windowed.records[0].windows = [{ start: 1, end: 100 }, { start: 200, end: 300 }];
  assert.equal(classifyImportedComparisonIntent({ renderRequest: windowed, resources: {} }).disposition,
    IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY);
  assert.equal(requestRunsProteinSearch(pairwiseRequest({ gridRows: [1, 1] })), false);
  // Two rows, the CLI's default of one row per record, or a listed pair: a search Generate runs.
  assert.equal(requestRunsProteinSearch(pairwiseRequest({ positions: ['#1@1', '#2@2'] })), true);
  assert.equal(requestRunsProteinSearch(pairwiseRequest({ gridRows: [1, 2] })), true);
  assert.equal(requestRunsProteinSearch(pairwiseRequest()), true);
  assert.equal(requestRunsProteinSearch(pairwiseRequest({
    positions: ['#1@1', '#2@1'], pairs: [{ queryRecordIndex: 0, subjectRecordIndex: 1 }]
  })), true);
});
