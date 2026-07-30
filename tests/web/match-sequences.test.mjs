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
  join(tempDir, 'color-utils.js'),
  await readFile('gbdraw/web/js/app/color-utils.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'conservation-series.js'),
  await readFile('gbdraw/web/js/app/conservation-series.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'match-sequences.js'),
  await readFile('gbdraw/web/js/app/match-sequences.js', 'utf8'),
  'utf8'
);

const {
  buildRestoredMatchSequenceSources,
  buildMatchSequenceBundle,
  createSequenceSourceRegistry,
  extractMatchedSpan,
  reverseComplementNucleotide
} = await import(pathToFileURL(join(tempDir, 'match-sequences.js')));

const textFile = (value) => ({ text: async () => value });
const namedTextFile = (name, value) => ({ name, text: async () => value });

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

test('restores transformed linear record sources for pairwise session popups', async () => {
  const sources = await buildRestoredMatchSequenceSources({
    mode: 'linear',
    lInputType: 'gb',
    linearSeqs: [
      {
        gb: textFile(`LOCUS       query\nACCESSION   QUERY\nVERSION     QUERY.1\nORIGIN\n        1 aaccggtt\n//\n`),
        region_record_id: 'QUERY.1',
        region_start: 2,
        region_end: 5,
        region_reverse: true
      },
      {
        gb: textFile(`LOCUS       subject\nACCESSION   SUBJECT\nVERSION     SUBJECT.1\nORIGIN\n        1 ttggccaa\n//\n`)
      }
    ]
  });

  assert.deepEqual(
    sources.map(({ key, recordId, sequence, origin, recordIndex }) => ({
      key, recordId, sequence, origin, recordIndex
    })),
    [
      {
        key: 'linear:record:0',
        recordId: 'QUERY.1',
        sequence: 'CGGT',
        origin: 'linear-record',
        recordIndex: 0
      },
      {
        key: 'linear:record:1',
        recordId: 'SUBJECT.1',
        sequence: 'TTGGCCAA',
        origin: 'linear-record',
        recordIndex: 1
      }
    ]
  );
});

test('restores circular reference and display-ordered conservation sources', async () => {
  const sources = await buildRestoredMatchSequenceSources({
    mode: 'circular',
    cInputType: 'gb',
    files: {
      c_gb: textFile(`LOCUS       ref\nACCESSION   REF\nVERSION     REF.1\nORIGIN\n        1 aaccggtt\n//\n`),
      c_conservation_fastas: [
        textFile('>comparison_a\nAAAACCCC\n'),
        textFile('>comparison_b\nGGGGTTTT\n')
      ]
    },
    circularConservation: {
      source: 'losat',
      series: [{ sourceIndex: 1 }, { sourceIndex: 0 }]
    }
  });

  assert.deepEqual(
    sources.map(({ key, recordId, sequence, origin, recordIndex, sourceIndex }) => ({
      key, recordId, sequence, origin, recordIndex, sourceIndex
    })),
    [
      {
        key: 'circular:record:0',
        recordId: 'REF.1',
        sequence: 'AACCGGTT',
        origin: 'circular-reference',
        recordIndex: 0,
        sourceIndex: undefined
      },
      {
        key: 'homology:comparison:0:comparison_b',
        recordId: 'comparison_b',
        sequence: 'GGGGTTTT',
        origin: 'homology-comparison',
        recordIndex: undefined,
        sourceIndex: 0
      },
      {
        key: 'homology:comparison:1:comparison_a',
        recordId: 'comparison_a',
        sequence: 'AAAACCCC',
        origin: 'homology-comparison',
        recordIndex: undefined,
        sourceIndex: 1
      }
    ]
  );
});

test('keeps missing LOSAT comparison FASTA placeholders aligned with displayed series', async () => {
  const sources = await buildRestoredMatchSequenceSources({
    mode: 'circular',
    cInputType: 'gb',
    files: {
      c_gb: textFile(`LOCUS       ref\nACCESSION   REF\nVERSION     REF.1\nORIGIN\n        1 aaccggtt\n//\n`),
      c_conservation_fastas: [
        null,
        textFile('>comparison_b\nGGGGTTTT\n')
      ]
    },
    circularConservation: {
      source: 'losat',
      series: [{ sourceIndex: 0 }, { sourceIndex: 1 }]
    }
  });

  assert.deepEqual(
    sources
      .filter((source) => source.origin === 'homology-comparison')
      .map(({ recordId, sequence, sourceIndex }) => ({ recordId, sequence, sourceIndex })),
    [
      { recordId: 'comparison_b', sequence: 'GGGGTTTT', sourceIndex: 1 }
    ]
  );
});

test('uses BLAST identity rather than stale source indexes for upload companions', async () => {
  const sources = await buildRestoredMatchSequenceSources({
    mode: 'circular',
    cInputType: 'gb',
    files: {
      c_gb: textFile(`LOCUS       ref\nACCESSION   REF\nVERSION     REF.1\nORIGIN\n        1 aaccggtt\n//\n`),
      c_conservation_blasts: [
        namedTextFile('comparison-a.tsv', ''),
        namedTextFile('comparison-b.tsv', '')
      ],
      c_conservation_sequence_sources: [
        namedTextFile('comparison-a.fna', '>comparison_a\nAAAACCCC\n'),
        namedTextFile('comparison-b.fna', '>comparison_b\nGGGGTTTT\n')
      ]
    },
    circularConservation: {
      source: 'upload',
      series: [
        { fileName: 'comparison-b.tsv', sourceIndex: 0 },
        { fileName: 'comparison-a.tsv', sourceIndex: 1 }
      ]
    }
  });

  assert.deepEqual(
    sources
      .filter((source) => source.origin === 'homology-comparison')
      .map(({ recordId, sequence, sourceIndex }) => ({ recordId, sequence, sourceIndex })),
    [
      { recordId: 'comparison_b', sequence: 'GGGGTTTT', sourceIndex: 0 },
      { recordId: 'comparison_a', sequence: 'AAAACCCC', sourceIndex: 1 }
    ]
  );
});
