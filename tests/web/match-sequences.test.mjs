import assert from 'node:assert/strict';
import test from 'node:test';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-match-sequences-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'));
await mkdir(join(tempDir, 'services'));
await writeFile(
  join(tempDir, 'app', 'feature-utils.js'),
  await readFile('gbdraw/web/js/app/feature-utils.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'feature-sequence-fasta.js'),
  await readFile('gbdraw/web/js/app/feature-sequence-fasta.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'color-utils.js'),
  await readFile('gbdraw/web/js/app/color-utils.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'conservation-series.js'),
  await readFile('gbdraw/web/js/app/conservation-series.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'match-sequences.js'),
  await readFile('gbdraw/web/js/app/match-sequences.js', 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'file-content-cache.js'),
  'export const readFileText = (file) => file.text();\n',
  'utf8'
);

const {
  analyzeCatalogSequenceSourceCoverage,
  buildRestoredMatchSequenceSources,
  buildMatchSequenceBundle,
  createSequenceSourceRegistry,
  extractMatchedSpan,
  resolveCircularComparisonSequenceAvailability,
  reverseComplementNucleotide
} = await import(pathToFileURL(join(tempDir, 'app', 'match-sequences.js')));

const textFile = (value) => ({ text: async () => value });
const namedTextFile = (name, value) => ({ name, text: async () => value });

const analyzeCircularHomologyCoverage = ({
  comparisonSource,
  comparisonSourceAvailability
} = {}) => analyzeCatalogSequenceSourceCoverage({
  mode: 'circular',
  renderRequest: {
    records: [{ recordKey: 'record-reference' }]
  },
  catalogFeatureState: {
    items: [{
      recordKeys: ['record-reference'],
      sequenceSources: [{
        key: 'circular:record:0',
        recordId: 'reference',
        aliases: ['reference'],
        sequence: 'AACCGGTT',
        origin: 'circular-reference',
        recordIndex: 0
      }, ...(comparisonSource ? [comparisonSource] : [])],
      biologicalFeatures: [],
      comparisonMatches: [{
        id: 'homology-a',
        match_kind: 'homology',
        query_record_id: 'comparison',
        subject_record_id: 'reference',
        source_index: '0',
        reference_side: 'subject'
      }]
    }]
  },
  comparisonSourceAvailability
});

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
  assert.equal(
    registry.resolve('', 'wrong-record', { origin: 'linear-record', recordIndex: 0 }).source,
    null
  );
  assert.match(
    registry.resolve('', 'wrong-record', { origin: 'linear-record', recordIndex: 0 }).reason,
    /No sequence source matched/
  );
});

test('registry validates every supplied source identity field and duplicate key', () => {
  const registry = createSequenceSourceRegistry([
    {
      key: 'linear:record:0',
      recordId: 'record-a',
      aliases: ['record-a.1'],
      sequence: 'AAAA',
      origin: 'linear-record',
      recordIndex: 0,
      sourceIndex: 2
    },
    {
      key: 'linear:record:0',
      recordId: 'record-a',
      aliases: ['record-a.1'],
      sequence: 'AAAA',
      origin: 'linear-record',
      recordIndex: 0,
      sourceIndex: 2
    },
    {
      key: 'invalid:record',
      recordId: 'invalid',
      sequence: 'CCCC',
      recordIndex: -1
    },
    {
      key: 'duplicate:key',
      recordId: 'first',
      sequence: 'GGGG'
    },
    {
      key: 'duplicate:key',
      recordId: 'second',
      sequence: 'TTTT'
    }
  ]);

  assert.equal(
    registry.values().filter((source) => source.key === 'linear:record:0').length,
    1
  );
  assert.equal(registry.resolve('invalid:record', 'invalid').source, null);
  assert.match(registry.resolve('duplicate:key', '').reason, /ambiguous/);

  for (const context of [
    { recordIndex: 'bogus' },
    { recordIndex: -1 },
    { recordIndex: Number.MAX_SAFE_INTEGER + 1 },
    { sourceIndex: '0.5' },
    { sourceIndex: -1 },
    { sourceIndex: Number.MAX_SAFE_INTEGER + 1 }
  ]) {
    assert.equal(
      registry.resolve('linear:record:0', 'record-a', context).source,
      null
    );
  }

  for (const [recordId, context] of [
    ['wrong-record', {}],
    ['record-a', { origin: 'homology-comparison' }],
    ['record-a', { recordIndex: 1 }],
    ['record-a', { sourceIndex: 3 }]
  ]) {
    assert.equal(
      registry.resolve('linear:record:0', recordId, context).source,
      null
    );
  }

  assert.equal(
    registry.resolve('linear:record:0', 'record-a.1', {
      origin: 'linear-record',
      recordIndex: '0',
      sourceIndex: 2
    }).source?.sequence,
    'AAAA'
  );
});

test('coverage reports a match endpoint whose displayed record source is missing', () => {
  const source = {
    key: 'linear:record:0',
    recordId: 'record-a',
    aliases: ['record-a'],
    sequence: 'AACCGGTT',
    origin: 'linear-record',
    recordIndex: 0
  };
  const coverage = analyzeCatalogSequenceSourceCoverage({
    mode: 'linear',
    renderRequest: {
      records: [{ recordKey: 'record-a' }, { recordKey: 'record-b' }]
    },
    catalogFeatureState: {
      items: [{
        recordKeys: ['record-a', 'record-b'],
        sequenceSources: [source],
        biologicalFeatures: [{
          recordKey: 'record-a',
          biologicalFeatureId: 'feature-a',
          sequenceSourceIndex: 0
        }],
        comparisonMatches: [{
          id: 'match-a-b',
          match_kind: 'pairwise',
          query_record_id: 'record-a',
          query_record_index: 0,
          subject_record_id: 'record-b',
          subject_record_index: 1
        }]
      }]
    }
  });

  assert.equal(coverage.complete, false);
  assert.deepEqual(coverage.consumerCounts, {
    biologicalFeatures: 1,
    matchEndpoints: 2
  });
  assert.equal(coverage.resolvedConsumers.length, 1);
  assert.equal(coverage.missingConsumers.length, 1);
  assert.equal(
    coverage.missingConsumers[0].expectedSource.key,
    'linear:record:1'
  );
  assert.deepEqual(
    coverage.missingConsumers[0].exampleConsumers,
    [{ kind: 'match-endpoint', id: 'match-a-b/subject' }]
  );
  assert.match(coverage.missingConsumers[0].reasons[0], /unavailable/);
});

test('coverage reports a recoverable homology comparison source missing from the catalog', () => {
  const comparisonSourceAvailability = resolveCircularComparisonSequenceAvailability({
    files: {
      c_conservation_blasts: [namedTextFile('comparison.tsv', '')],
      c_conservation_sequence_sources: [namedTextFile('comparison.fna', '>comparison\nAAAA\n')]
    },
    circularConservation: {
      source: 'upload',
      series: [{ sourceIndex: 0 }]
    }
  });
  const coverage = analyzeCircularHomologyCoverage({ comparisonSourceAvailability });

  assert.equal(coverage.complete, false);
  assert.deepEqual(coverage.consumerCounts, {
    biologicalFeatures: 0,
    matchEndpoints: 2
  });
  assert.equal(coverage.resolvedConsumers.length, 1);
  assert.equal(coverage.missingConsumers.length, 1);
  assert.equal(
    coverage.missingConsumers[0].expectedSource.origin,
    'homology-comparison'
  );
  assert.equal(coverage.missingConsumers[0].expectedSource.sourceIndex, 0);
  assert.deepEqual(
    coverage.missingConsumers[0].exampleConsumers,
    [{ kind: 'match-endpoint', id: 'homology-a/query' }]
  );
  assert.deepEqual(coverage.optionalUnavailableConsumers, []);
});

test('coverage reports a never-supplied homology comparison FASTA as optionally unavailable', () => {
  const comparisonSourceAvailability = resolveCircularComparisonSequenceAvailability({
    files: {
      c_conservation_blasts: [namedTextFile('comparison.tsv', '')],
      c_conservation_sequence_sources: []
    },
    circularConservation: {
      source: 'upload',
      series: [{ sourceIndex: 0 }]
    }
  });
  const coverage = analyzeCircularHomologyCoverage({ comparisonSourceAvailability });

  assert.equal(coverage.complete, true);
  assert.deepEqual(coverage.missingConsumers, []);
  assert.equal(coverage.optionalUnavailableConsumers.length, 1);
  assert.equal(
    coverage.optionalUnavailableConsumers[0].expectedSource.sourceIndex,
    0
  );
  assert.deepEqual(
    coverage.optionalUnavailableConsumers[0].exampleConsumers,
    [{ kind: 'match-endpoint', id: 'homology-a/query' }]
  );
});

test('coverage resolves a catalog comparison source without consulting optional availability', () => {
  const coverage = analyzeCircularHomologyCoverage({
    comparisonSource: {
      key: 'homology:comparison:0:comparison',
      recordId: 'comparison',
      aliases: ['comparison'],
      sequence: 'AAAACCCC',
      origin: 'homology-comparison',
      sourceIndex: 0
    }
  });

  assert.equal(coverage.complete, true);
  assert.deepEqual(coverage.missingConsumers, []);
  assert.deepEqual(coverage.optionalUnavailableConsumers, []);
  assert.equal(coverage.resolvedConsumers.length, 2);
});

test('coverage fails closed when comparison source availability is unknown', () => {
  const coverage = analyzeCircularHomologyCoverage({
    comparisonSourceAvailability: []
  });

  assert.equal(coverage.complete, false);
  assert.equal(coverage.missingConsumers.length, 1);
  assert.deepEqual(coverage.optionalUnavailableConsumers, []);
  assert.match(coverage.missingConsumers[0].reasons[0], /availability is missing/);
});

test('comparison availability preserves nulls and duplicate-name source ordering', () => {
  const firstBlast = {
    name: 'duplicate.tsv',
    size: 12,
    lastModified: 34,
    text: async () => ''
  };
  const secondBlast = { ...firstBlast };
  const firstFasta = namedTextFile('first.fna', '>first\nAAAA\n');
  const availability = resolveCircularComparisonSequenceAvailability({
    files: {
      c_conservation_blasts: [firstBlast, secondBlast],
      c_conservation_sequence_sources: [firstFasta, null]
    },
    circularConservation: {
      source: 'upload',
      series: [
        { sourceKey: 'duplicate.tsv|12|34|1', sourceIndex: 1 },
        { sourceKey: 'duplicate.tsv|12|34|0', sourceIndex: 0 }
      ]
    }
  });

  assert.deepEqual(availability, [
    { sourceIndex: 0, file: null, availability: 'never-supplied' },
    { sourceIndex: 1, file: firstFasta, availability: 'recoverable' }
  ]);
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
