import assert from 'node:assert/strict';
import test from 'node:test';
import {
  analyzeCatalogSequenceSourceCoverage,
  buildRestoredMatchSequenceSources,
  buildMatchSequenceBundle,
  createSequenceSourceRegistry,
  extractMatchedSpan,
  resolveCircularComparisonSequenceAvailability,
  reverseComplementNucleotide
} from '../../gbdraw/web/js/services/match-sequences.js';
import { readRecordSourceSpan, recordSourceInterval } from '../../gbdraw/web/js/services/record-source-coordinates.js';

// The restore reads files through the shared browser sequence reader (`readParsedSequenceRecords`), which reads a
// File's bytes once per File object.
const textFile = (value) => new File([value], 'input.txt');
const namedTextFile = (name, value) => new File([value], name);

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

// OV-1097 (Owner decision 2026-10-18): the restored record is the region
// Generate draws, so an end past the record end is read to the record end, and
// a region that starts after the record ends has no sequence.
test('a restored region past the record end holds the record up to its end (OV-1097)', async () => {
  const gb = () => textFile(`LOCUS       query\nACCESSION   QUERY\nVERSION     QUERY.1\nORIGIN\n        1 aaccggtt\n//\n`);
  const sources = await buildRestoredMatchSequenceSources({
    mode: 'linear',
    lInputType: 'gb',
    linearSeqs: [
      { gb: gb(), region_record_id: 'QUERY.1', region_start: 5, region_end: 20 },
      { gb: gb(), region_record_id: 'QUERY.1', region_start: 9, region_end: 20 }
    ]
  });
  assert.deepEqual(sources.map(({ key, sequence }) => ({ key, sequence })), [
    { key: 'linear:record:0', sequence: 'GGTT' }
  ]);
});

test('a match on a restored cropped reverse record reports source coordinates (CO-10)', async () => {
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
  // The registry holds the displayed record: source 2..5 reverse-complemented.
  assert.deepEqual(
    sources.map(({ key, recordId, sequence, origin, recordIndex }) => ({
      key, recordId, sequence, origin, recordIndex
    })),
    [
      { key: 'linear:record:0', recordId: 'QUERY.1', sequence: 'CGGT', origin: 'linear-record', recordIndex: 0 },
      { key: 'linear:record:1', recordId: 'SUBJECT.1', sequence: 'TTGGCCAA', origin: 'linear-record', recordIndex: 1 }
    ]
  );

  // The renderer marks the cropped reverse record with its input-file span.
  const group = {
    getAttribute: (name) => ({
      'data-gbdraw-record-index': '0',
      'data-gbdraw-record-source-start': '2',
      'data-gbdraw-record-source-end': '5',
      'data-gbdraw-record-source-step': '-1'
    })[name] ?? null
  };
  const svg = {
    querySelectorAll: (selector) => (
      selector === '[data-gbdraw-record-index][data-gbdraw-record-source-step]' ? [group] : []
    )
  };
  const match = { ownerSVGElement: svg };
  const sourceSpan = readRecordSourceSpan(match, 0);
  assert.deepEqual(sourceSpan, { start: 2, end: 5, step: -1 });
  assert.equal(readRecordSourceSpan(match, 1), null);
  // Displayed bases 1..2 are source 5..4; the table (search frame) gives 4..3.
  assert.deepEqual(recordSourceInterval(sourceSpan, 1, 2), { start: 5, end: 4, table: { start: 4, end: 3 } });
  assert.deepEqual(recordSourceInterval({ start: 1, end: 8, step: -1 }, 1, 2), { start: 8, end: 7, table: null });

  const registry = createSequenceSourceRegistry(sources);
  const bundle = buildMatchSequenceBundle([
    { role: 'query', sourceKey: 'linear:record:0', recordId: 'QUERY.1', start: 1, end: 2, sourceSpan },
    { role: 'subject', sourceKey: 'linear:record:1', recordId: 'SUBJECT.1', start: 3, end: 4, sourceSpan: null }
  ], { matchId: 'pairwise_comparison1_match1', resolveSequenceSource: registry.resolve });
  // The sequence is unchanged; only the header coordinates are source based.
  assert.equal(bundle.entries[0].fasta, '>pairwise_comparison1_match1_query|record=QUERY.1|coords=5..4|strand=-\nCG\n');
  assert.equal(bundle.entries[0].filename, 'pairwise_comparison1_match1_query_QUERY.1_5-4.fna');
  assert.equal(bundle.entries[1].fasta, '>pairwise_comparison1_match1_subject|record=SUBJECT.1|coords=3..4|strand=+\nGG\n');
});

test('record source spans are read from one scan of the Result (OV-931)', () => {
  let scans = 0;
  const groupFor = (index, step) => ({
    getAttribute: (name) => ({
      'data-gbdraw-record-index': String(index),
      'data-gbdraw-record-source-start': '2',
      'data-gbdraw-record-source-end': '5',
      'data-gbdraw-record-source-step': String(step)
    })[name] ?? null
  });
  const groups = [groupFor(0, -1), groupFor(1, 1)];
  const svg = {
    querySelector: (selector) => {
      scans += 1;
      const index = /data-gbdraw-record-index="(\d+)"/.exec(selector)?.[1];
      return groups.find((group) => group.getAttribute('data-gbdraw-record-index') === index) || null;
    },
    querySelectorAll: (selector) => {
      scans += 1;
      return selector === '[data-gbdraw-record-index][data-gbdraw-record-source-step]' ? groups : [];
    }
  };
  const match = { ownerSVGElement: svg };
  for (let round = 0; round < 50; round += 1) {
    assert.deepEqual(readRecordSourceSpan(match, 0), { start: 2, end: 5, step: -1 });
    assert.deepEqual(readRecordSourceSpan(match, 1), { start: 2, end: 5, step: 1 });
    assert.equal(readRecordSourceSpan(match, 2), null);
  }
  assert.equal(scans, 1);
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

// OV-1680: a restored Linear card selects its record among every record the Python reader reads. A GenBank
// record without a sequence (a CON record) keeps its place, so `#2` names the record the drawing drew.
test('a restored record selector counts a record without a sequence (OV-1680)', async () => {
  const record = (name, division, body) => [
    `LOCUS       ${name}                       8 bp    DNA     linear   ${division} 01-JAN-1980`,
    'DEFINITION  x.', `ACCESSION   ${name}`, `VERSION     ${name}`, 'FEATURES             Location/Qualifiers',
    body, '//'
  ].join('\n');
  const gb = textFile([
    record('BBB', 'CON', 'CONTIG      join(AAA:1..8)'),
    record('DUP', 'UNK', 'ORIGIN\n        1 aaaaaaaa'),
    record('DUP', 'UNK', 'ORIGIN\n        1 cccccccc')
  ].join('\n') + '\n');
  const sources = await buildRestoredMatchSequenceSources({
    mode: 'linear', lInputType: 'gb', linearSeqs: [{ gb, region_record_id: '#2' }, { gb, region_record_id: '#1' }]
  });
  // Python: `#2` is the first DUP; `#1` is BBB, which has no sequence to register.
  assert.deepEqual(sources.map(({ key, recordId, sequence }) => ({ key, recordId, sequence })), [
    { key: 'linear:record:0', recordId: 'DUP', sequence: 'AAAAAAAA' }
  ]);
});

// OV-1681: a restored GFF3 card's `#N` counts the records the GFF3 names, as the drawing's loader does.
test('a restored GFF3 record selector counts the records the GFF3 names (OV-1681)', async () => {
  const gff = namedTextFile('named.gff',
    '##gff-version 3\nG1\t.\tCDS\t1\t9\t.\t+\t0\tID=c1\nG2\t.\tCDS\t1\t9\t.\t+\t0\tID=c2\n');
  const fasta = namedTextFile('extra.fasta', '>X\nCCCCCCCCC\n>G1\nATGAAATAA\n>G2\nATGGGGTAA\n');
  const sources = await buildRestoredMatchSequenceSources({
    mode: 'linear', lInputType: 'gff', linearSeqs: [{ gff, fasta, region_record_id: '#2' }]
  });
  assert.deepEqual(sources.map(({ recordId, sequence }) => ({ recordId, sequence })), [
    { recordId: 'G2', sequence: 'ATGGGGTAA' }
  ]);
});

// OV-1681: the automatic record of a restored GFF3 card is the first record the GFF3 names, and a card without
// both files keeps its selector, so the other cards keep their sources.
test('a restored GFF3 automatic record is the first record the GFF3 names (OV-1681)', async () => {
  const gff = namedTextFile('g1-only.gff', '##gff-version 3\nG1\t.\tCDS\t1\t9\t.\t+\t0\tID=c1\n');
  const fasta = namedTextFile('x-first.fasta', '>X\nCCCCCCCCC\n>G1\nATGAAATAA\n');
  const automatic = await buildRestoredMatchSequenceSources({
    mode: 'linear', lInputType: 'gff', linearSeqs: [{ gff, fasta, region_record_id: '' }]
  });
  assert.deepEqual(automatic.map(({ recordId, sequence }) => ({ recordId, sequence })), [
    { recordId: 'G1', sequence: 'ATGAAATAA' }
  ]);
  const missingGff = await buildRestoredMatchSequenceSources({
    mode: 'linear', lInputType: 'gff',
    linearSeqs: [{ gff, fasta, region_record_id: 'G1' }, { gff: null, fasta, region_record_id: '#1' }]
  });
  assert.deepEqual(missingGff.map(({ key, recordId }) => ({ key, recordId })), [
    { key: 'linear:record:0', recordId: 'G1' },
    { key: 'linear:record:1', recordId: 'X' }
  ]);
});
