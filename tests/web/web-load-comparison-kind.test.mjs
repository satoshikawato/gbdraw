// Web Load keeps the comparison kind a CLI Session stored.
// OV-3573: a request of schema 2 or 5 to 9 stored a CLI comparisons-table BLAST table as a protein table
// (`precomputedProteinComparison` of the BLAST columns alone); Load reads it as the BLAST table it is.
// OV-3566: a LOSATP Pairwise search that compared no pair (every record in one row) is stored as the Pairwise
// pipeline without pairs; the controls hold it, and a redraw of it searches nothing.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';

import {
  COMPARISON_COLUMNS, LEGACY_GIVEN_COMPARISON_RESOURCE_PREFIX
} from '../../gbdraw/web/js/contracts/schema.generated.js';
import {
  readStoredBlastTables, requestRunsProteinSearch, withStoredTablesWithoutProteinMatches
} from '../../gbdraw/web/js/services/session-request.js';
import { pythonLegacySessionUpgrader } from './helpers/python-evaluator.mjs';
import {
  IMPORTED_COMPARISON_DISPOSITIONS, classifyImportedComparisonIntent, storedProteinResultMode
} from '../../gbdraw/web/js/services/imported-comparison-intent.js';

const base64 = (text) => Buffer.from(text, 'utf8').toString('base64');
const text = (data) => Buffer.from(data, 'base64').toString('utf8');
const BLAST_ROWS = 'rec0\trec1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n';
const table = (header, rows = BLAST_ROWS) => ({
  kind: 'canonical-tsv', name: 'comparison-explicit-1.tsv', type: 'text/tab-separated-values',
  size: 0, lastModified: 0, encoding: 'base64', data: base64(`${header}\n${rows}`)
});

const given = (number) => `${LEGACY_GIVEN_COMPARISON_RESOURCE_PREFIX}${number}`;

test('a schema 5-9 protein table of the BLAST columns alone under a given table\'s ID loads as a BLAST table', () => {
  assert.equal(given(1), 'comparison-explicit-1', 'the ID 0.14.x CLI and Python renders gave a table they were given');
  const marker = { kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: {} };
  const blastHeader = COMPARISON_COLUMNS.join('\t');
  const resources = {
    blast: table(blastHeader),
    crlf: table(`${blastHeader}\r`),
    protein: table([...COMPARISON_COLUMNS, 'query_protein_id', 'subject_protein_id'].join('\t')),
    // The IDs the Web and the CLI gave LOSATP tables: never read by the columns alone.
    webLosatp: table(blastHeader),
    cliLosatp: table(blastHeader),
    resolved: table(blastHeader),
    notNumbered: table(blastHeader)
  };
  const ids = {
    blast: given(1), crlf: given(2), protein: given(3), webLosatp: 'comparison-canonical-protein-1',
    cliLosatp: 'comparison-protein-1', resolved: 'comparison-resolved-1', notNumbered: given('1-copy')
  };
  Object.assign(resources, Object.fromEntries(Object.entries(ids).map(([key, id]) => [id, resources[key]])));
  for (const key of Object.keys(ids)) delete resources[key];
  const comparison = (key, queryRecordIndex = 0, subjectRecordIndex = 1) => ({
    kind: 'precomputedProteinComparison', resourceId: ids[key], encoding: 'canonicalTsv', queryRecordIndex,
    subjectRecordIndex
  });
  const stored = { schema: 9, mode: 'linear', comparisons: [comparison('blast'), comparison('protein', 1, 2), marker] };
  const crlf = { schema: 5, mode: 'linear', comparisons: [comparison('crlf')] };
  const current = { schema: 10, mode: 'linear', comparisons: [comparison('blast')] };
  const read = readStoredBlastTables([stored, crlf, current], resources);
  assert.deepEqual(read.requests[0].comparisons, [
    { kind: 'nucleotideBlast', resourceId: ids.blast, queryRecordIndex: 0, subjectRecordIndex: 1 },
    stored.comparisons[1],
    marker
  ]);
  assert.equal(read.requests[0].comparisons[1], stored.comparisons[1]);
  assert.equal(read.requests[1].comparisons[0].kind, 'nucleotideBlast');
  // The request of the current schema is never rewritten.
  assert.equal(read.requests[2], current);
  for (const key of ['blast', 'crlf']) {
    const resource = read.resources[ids[key]];
    assert.equal(resource.kind, 'nucleotide-blast');
    assert.equal(text(resource.data), BLAST_ROWS);
    assert.equal(resource.size, Buffer.byteLength(BLAST_ROWS));
    assert.equal(resource.name, 'comparison-explicit-1.tsv');
  }
  for (const key of ['protein', 'webLosatp', 'cliLosatp', 'resolved', 'notNumbered']) {
    assert.equal(read.resources[ids[key]], resources[ids[key]], key);
  }
  const losatp = readStoredBlastTables([{ schema: 9, mode: 'linear', comparisons: [
    comparison('webLosatp'), comparison('cliLosatp', 1, 2), comparison('resolved', 2, 3),
    comparison('notNumbered', 3, 4)
  ] }], resources);
  assert.equal(losatp.resources, resources);
  assert.ok(losatp.requests[0].comparisons.every((item) => item.kind === 'precomputedProteinComparison'));
  // The file's own values are not edited.
  assert.equal(resources[ids.blast].kind, 'canonical-tsv');
  assert.equal(stored.comparisons[0].kind, 'precomputedProteinComparison');

  // Schema 2 (main from 2026-07-18) wrote the given table under the same ID, as Python's decoder reads it
  // (`LEGACY_GIVEN_COMPARISON_SCHEMAS`); schema 1 had no record pairs.
  const schema2 = readStoredBlastTables([{ schema: 2, mode: 'linear', comparisons: [comparison('blast')] }], resources);
  assert.equal(schema2.requests[0].comparisons[0].kind, 'nucleotideBlast');
  assert.equal(schema2.resources[ids.blast].kind, 'nucleotide-blast');
  // A table bound to features has no BLAST form: the pre-pass leaves the file as stored, and Load's request reader
  // reads it with `proteinMatches: false`.
  const boundResources = { [given(1)]: table([...COMPARISON_COLUMNS, 'query_feature_index', 'subject_feature_index'].join('\t'),
    BLAST_ROWS.replace('\n', '\t1\t2\n')) };
  const boundRequest = { schema: 9, mode: 'linear', comparisons: [comparison('blast')] };
  const bound = readStoredBlastTables([boundRequest], boundResources);
  assert.equal(bound.resources, boundResources);
  assert.equal(bound.requests[0], boundRequest);
  assert.deepEqual(withStoredTablesWithoutProteinMatches(boundRequest, boundResources).comparisons,
    [{ ...comparison('blast'), proteinMatches: false }]);

  const protein = { schema: 9, mode: 'linear', comparisons: [comparison('protein')] };
  const unchanged = readStoredBlastTables([current, protein], resources);
  assert.equal(unchanged.resources, resources);
  assert.equal(unchanged.requests[1].comparisons[0].kind, 'precomputedProteinComparison');
});

// The vectors Python's decoder reads too (tests/test_cli_session_render.py): one rule on both surfaces.
const GIVEN_TABLE_VECTORS = JSON.parse(readFileSync(new URL('../fixtures/given-comparison-table-vectors.json', import.meta.url),
  'utf8')).cases;
for (const vector of GIVEN_TABLE_VECTORS) {
  test(`a stored table is read as Python reads it: ${vector.id}`, () => {
    const stored = {
      kind: 'precomputedProteinComparison', resourceId: vector.resourceId, encoding: 'canonicalTsv',
      queryRecordIndex: 0, subjectRecordIndex: 1
    };
    const tsv = [vector.columns, ...vector.rows].map((cells) => `${cells.join('\t')}\n`).join('');
    const resources = { [vector.resourceId]: { ...table('', ''), data: base64(tsv) } };
    const request = { schema: vector.schema, mode: 'linear', comparisons: [stored] };
    const read = readStoredBlastTables([request], resources);
    // The request Load's reader classifies and promotes.
    const readRequest = withStoredTablesWithoutProteinMatches(read.requests[0], read.resources);
    const [comparison] = readRequest.comparisons;
    if (vector.read === 'protein') {
      assert.equal(readRequest, request);
      assert.equal(read.resources, resources);
    } else if (vector.read === 'table') {
      // The file keeps the request as stored: its schema has no `proteinMatches`.
      assert.equal(read.requests[0], request);
      assert.deepEqual(comparison, { ...stored, proteinMatches: false });
      assert.equal(read.resources, resources);
    } else {
      assert.equal(vector.read, 'blast');
      assert.deepEqual(comparison, { kind: 'nucleotideBlast', resourceId: vector.resourceId, queryRecordIndex: 0,
        subjectRecordIndex: 1 });
      assert.equal(read.resources[vector.resourceId].kind, 'nucleotide-blast');
      assert.equal(text(read.resources[vector.resourceId].data), vector.blastRows);
    }
  });
}

// OV-3574: a table that holds no protein matches is no stored LOSATP result: beside the CLI's no-search pipeline
// the comparison is preserved and named by its tables, not opened as a LOSATP Pairwise comparison.
test('a table without protein matches is no stored protein result', () => {
  const marker = { kind: 'generatedProteinComparison', mode: 'none', pairs: [], settings: pipelineSettings };
  const bound = {
    kind: 'precomputedProteinComparison', resourceId: given(1), encoding: 'canonicalTsv', queryRecordIndex: 0,
    subjectRecordIndex: 1, proteinMatches: false
  };
  const { proteinMatches: _matches, ...protein } = bound;
  assert.equal(storedProteinResultMode([bound, marker]), '');
  assert.equal(storedProteinResultMode([protein, marker]), 'pairwise');
  const renderRequest = { ...pairwiseRequest({ positions: ['#1@1', '#2@2'] }), comparisons: [bound, marker] };
  const resources = { [given(1)]: table([...COMPARISON_COLUMNS, 'query_feature_index'].join('\t'), '') };
  const intent = classifyImportedComparisonIntent({ renderRequest, resources });
  assert.equal(intent.disposition, IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY);
  assert.equal(intent.message,
    'The saved comparison tables can run, but the current controls cannot represent them losslessly.');
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
  // A record entry of every record of its file (a Python API `RecordCardinality.ALL`): its rows are Python's, and
  // the search runs as on dev.
  const everyRecord = pairwiseRequest({ gridRows: [null] });
  everyRecord.records[0].cardinality = 'all';
  assert.equal(requestRunsProteinSearch(everyRecord), true);
  assert.equal(requestRunsProteinSearch(pairwiseRequest({
    positions: ['#1@1', '#2@1'], pairs: [{ queryRecordIndex: 0, subjectRecordIndex: 1 }]
  })), true);
});

// A Web LOSATP Session that main's Gallery shipped (Session 40, request schema 5; provenance in its
// .provenance.json): its protein tables start with the BLAST columns too. Load keeps them protein tables.
test('a released Web LOSATP Session keeps its protein tables and its LOSATP pipeline', () => {
  const fixture = new URL('../fixtures/performance/hepatoplasmataceae_collinear.2afba711.gbdraw-session.json.gz',
    import.meta.url);
  const session = JSON.parse(gunzipSync(readFileSync(fixture)).toString('utf8'));
  const stored = session.renderRequest;
  assert.equal(stored.schema, 5);
  const tables = stored.comparisons.filter((item) => item.kind === 'precomputedProteinComparison');
  assert.equal(tables.length, 4);
  for (const { resourceId } of tables) {
    assert.ok(text(session.resources[resourceId].data).startsWith(`${COMPARISON_COLUMNS.join('\t')}\t`));
  }
  // The pre-pass hands Load the very request and resources the file holds, so Load reads them as before.
  const read = readStoredBlastTables([stored], session.resources);
  assert.equal(read.requests[0], stored);
  assert.equal(read.resources, session.resources);
});

// OV-3589: a Session whose comparisons are BLAST tables (a CLI -b or comparisons-table render, which stores the
// protein pipeline as a no-search marker) is preserved, and its message names the tables, not a protein comparison.
test('a preserved BLAST table Session is named by its tables, a preserved protein comparison as one', () => {
  const fixture = new URL('../fixtures/sessions/BGC0000708-BGC0000709-comparisons-table.v50-cli0145.gbdraw-session.json.gz',
    import.meta.url);
  const session = JSON.parse(gunzipSync(readFileSync(fixture)).toString('utf8'));
  const read = readStoredBlastTables([session.renderRequest], session.resources);
  const tables = classifyImportedComparisonIntent({ renderRequest: read.requests[0], resources: read.resources });
  assert.equal(tables.disposition, IMPORTED_COMPARISON_DISPOSITIONS.PRESERVED_READ_ONLY);
  assert.equal(tables.message,
    'The saved comparison tables can run, but the current controls cannot represent them losslessly.');
  // A protein pipeline the controls cannot represent (its collinearity parameters are not lossless) keeps its message.
  const protein = classifyImportedComparisonIntent({ renderRequest: {
    ...pairwiseRequest({ positions: ['#1@1', '#2@2'] }),
    comparisons: [{ ...pairwiseRequest().comparisons[0], settings: {
      ...pipelineSettings,
      collinearityParams: { kind: 'lossless', parameters: { minAnchors: 1, maxUnitGap: 0, maxDiagonalDrift: 0, maxConflicts: 1 } }
    } }]
  }, resources: {} });
  assert.equal(protein.message, 'The saved protein comparison is executable, but current controls cannot represent it losslessly.');
});

// Python converts a Session 27-44 from the file Load read (`upgradeLegacySession`), whose requests keep their schema:
// the pre-pass writes no field that schema lacks. Session 40 (request schema 5, a Python LOSATP table under a given
// table's ID, its protein column removed: a table bound to features without protein matches).
test('a Session 27-44 with a table bound to features still converts after the pre-pass', async () => {
  const fixture = new URL('../fixtures/performance/hepatoplasmataceae_collinear.2afba711.gbdraw-session.json.gz',
    import.meta.url);
  const session = JSON.parse(gunzipSync(readFileSync(fixture)).toString('utf8'));
  const stored = session.renderRequest.comparisons.find((item) => item.kind === 'precomputedProteinComparison');
  const resource = session.resources[stored.resourceId];
  delete session.resources[stored.resourceId];
  const lines = text(resource.data).split('\n');
  const proteinColumn = lines[0].split('\t').indexOf('query_protein_id');
  const data = lines.map((line) => (line ? line.split('\t').filter((_cell, index) => index !== proteinColumn).join('\t') : line))
    .join('\n');
  const { checksum: _checksum, ...kept } = resource;
  session.resources[given(1)] = { ...kept, data: base64(data), size: Buffer.byteLength(data) };
  stored.resourceId = given(1);
  const read = readStoredBlastTables([session.renderRequest], session.resources);
  assert.equal(read.requests[0], session.renderRequest);
  assert.equal(withStoredTablesWithoutProteinMatches(session.renderRequest, session.resources).comparisons
    .find((item) => item.resourceId === given(1)).proteinMatches, false);
  const document = { ...session, renderRequest: read.requests[0], resources: read.resources };
  const { result } = await pythonLegacySessionUpgrader()('upgradeLegacySession', { document, storedRequestSchema: 5 });
  assert.equal(result.error, undefined);
});
