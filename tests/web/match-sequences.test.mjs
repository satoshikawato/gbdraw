import assert from 'node:assert/strict';
import test from 'node:test';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-match-sequences-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'feature-utils.js'),
  await readFile('gbdraw/web/js/app/feature-utils.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'feature-sequence-fasta.js'),
  await readFile('gbdraw/web/js/app/feature-sequence-fasta.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'match-sequences.js'),
  await readFile('gbdraw/web/js/app/match-sequences.js', 'utf8'),
  'utf8'
);

const {
  buildMatchSequenceBundle,
  createSequenceSourceRegistry,
  extractMatchedSpan,
  reverseComplementNucleotide
} = await import(pathToFileURL(join(tempDir, 'match-sequences.js')));

test('extracts 1-based inclusive forward and reverse IUPAC spans', () => {
  assert.deepEqual(extractMatchedSpan('AACCGGTT', 2, 5), {
    valid: true,
    start: 2,
    end: 5,
    orientation: '+',
    sequence: 'ACCG',
    sequenceLength: 4
  });
  assert.equal(extractMatchedSpan('AACCGGTT', 5, 2).sequence, 'CGGT');
  assert.equal(reverseComplementNucleotide('ARYKMBVDHN'), 'NDHBVKMRYT');
});

test('rejects invalid and out-of-bounds coordinates', () => {
  assert.equal(extractMatchedSpan('ACGT', 0, 2).valid, false);
  assert.equal(extractMatchedSpan('ACGT', 1.5, 2).valid, false);
  assert.match(extractMatchedSpan('ACGT', 1, 5).reason, /sequence length/);
});

test('registry refuses ambiguous record IDs and resolves stable keys first', () => {
  const registry = createSequenceSourceRegistry([
    { key: 'linear:record:0', recordId: 'dup', sequence: 'AAAA', origin: 'linear-record', recordIndex: 0 },
    { key: 'linear:record:1', recordId: 'dup', sequence: 'CCCC', origin: 'linear-record', recordIndex: 1 }
  ]);
  assert.equal(registry.resolve('linear:record:1', 'dup').source.sequence, 'CCCC');
  assert.match(registry.resolve('', 'dup', { origin: 'linear-record' }).reason, /ambiguous/);
});

test('builds deterministic single and combined FASTA in query-subject order', () => {
  const registry = createSequenceSourceRegistry([
    { key: 'linear:record:0', recordId: 'query', sequence: 'AACCGG', origin: 'linear-record', recordIndex: 0 },
    { key: 'linear:record:1', recordId: 'subject', sequence: 'TTGGCC', origin: 'linear-record', recordIndex: 1 }
  ]);
  const bundle = buildMatchSequenceBundle([
    { role: 'query', sourceKey: 'linear:record:0', recordId: 'query', start: 2, end: 4, displayRole: 'Query' },
    { role: 'subject', sourceKey: 'linear:record:1', recordId: 'subject', start: 5, end: 3, displayRole: 'Subject' }
  ], {
    matchId: 'pairwise_comparison1_match3',
    resolveSequenceSource: registry.resolve
  });
  assert.equal(bundle.entries[0].fasta, '>pairwise_comparison1_match3_query|record=query|coords=2..4|strand=+\nACC\n');
  assert.equal(bundle.entries[1].fasta, '>pairwise_comparison1_match3_subject|record=subject|coords=5..3|strand=-\nGCC\n');
  assert.equal(bundle.combinedFasta, `${bundle.entries[0].fasta}${bundle.entries[1].fasta}`);
  assert.equal(bundle.combinedFilename, 'pairwise_comparison1_match3_both.fna');
});
