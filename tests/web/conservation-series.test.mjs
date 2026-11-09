// D12: a Circular similarity ring row added from a GenBank or DDBJ file takes
// the label the file names itself (first record's DEFINITION, then organism),
// as read by the Python ring reader; FASTA rows keep the file-name default and a
// typed label always wins.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

const {
  applyComparisonSequenceRecordLabel,
  defaultConservationSeriesLabel,
  normalizeConservationSeriesColor,
  reconcileConservationSeries,
  withoutConservationSeries
} = await import('../../gbdraw/web/js/services/conservation-series.js');
const { ringGenomeFileLists } = await import('../../gbdraw/web/js/services/drawing/inputs.js');

const file = (name, size = 10) => ({ name, size, lastModified: 0 });

const rowsFor = (sourceFiles) => reconcileConservationSeries({
  sourceFiles,
  previousSeries: [],
  legacyLabels: []
});

test('a GenBank or DDBJ row without a typed label takes the record label', () => {
  const genbank = file('comparison-c.gbk');
  const ddbj = file('comparison-d.ddbj', 20);
  const sourceFiles = [file('comparison-b.fasta'), genbank, ddbj];
  const series = rowsFor(sourceFiles);
  assert.deepEqual(series.map(({ label }) => label), ['comparison-b', 'comparison-c', 'comparison-d']);

  assert.equal(applyComparisonSequenceRecordLabel({
    series, sourceFiles, file: genbank, recordLabel: 'synthetic comparison c'
  }), true);
  assert.equal(applyComparisonSequenceRecordLabel({
    series, sourceFiles, file: ddbj, recordLabel: 'Synthetic organism d'
  }), true);
  assert.deepEqual(series.map(({ label }) => label), [
    'comparison-b', 'synthetic comparison c', 'Synthetic organism d'
  ]);
});

test('a FASTA row or a file that names no label keeps the file-name default', () => {
  const fasta = file('comparison-b.fasta');
  const bare = file('bare.gbk');
  const sourceFiles = [fasta, bare];
  const series = rowsFor(sourceFiles);
  for (const [target, recordLabel] of [[fasta, null], [bare, ''], [bare, '   '], [bare, undefined]]) {
    assert.equal(applyComparisonSequenceRecordLabel({
      series, sourceFiles, file: target, recordLabel
    }), false);
  }
  assert.deepEqual(series.map(({ label }) => label), ['comparison-b', 'bare']);
});

test('the file-name default follows the shared vectors of the Python ring reader', async () => {
  // D-03: tests/test_comparison_sequences.py runs the same cases.
  const { cases } = JSON.parse(await readFile(
    new URL('../fixtures/comparison_ring_default_label_cases.json', import.meta.url), 'utf8'
  ));
  for (const { fileName, expected } of cases) {
    assert.equal(defaultConservationSeriesLabel(file(fileName), 0), expected, fileName);
  }
});

test('a typed label wins over the record label', () => {
  const genbank = file('comparison-c.gbk');
  const sourceFiles = [genbank];
  const series = rowsFor(sourceFiles);
  series[0].label = 'My ring';
  assert.equal(applyComparisonSequenceRecordLabel({
    series, sourceFiles, file: genbank, recordLabel: 'synthetic comparison c'
  }), false);
  assert.equal(series[0].label, 'My ring');
});

test('a row whose file was removed or replaced is not relabeled', () => {
  const genbank = file('comparison-c.gbk');
  const replacement = file('comparison-c.gbk');
  const sourceFiles = [replacement];
  const series = rowsFor(sourceFiles);
  assert.equal(applyComparisonSequenceRecordLabel({
    series, sourceFiles, file: genbank, recordLabel: 'synthetic comparison c'
  }), false);
  assert.equal(applyComparisonSequenceRecordLabel({
    series: [], sourceFiles: [], file: genbank, recordLabel: 'synthetic comparison c'
  }), false);
  assert.equal(series[0].label, 'comparison-c');
});

// OV-893: the one removal of a series (Remove, and Load for a series whose
// table it could not read) takes the place from every list; the series after
// it keep their own files and take the place before.
test('withoutConservationSeries leaves out one series with its files; later series keep theirs', () => {
  const series = ['b', 'c', 'd'].map((label, sourceIndex) => ({ label, sourceIndex, color: `#${sourceIndex}` }));
  const kept = withoutConservationSeries({
    blasts: ['b.tsv', null, 'd.tsv'], fastas: ['b.fa', 'c.fa', 'd.fa'], sequenceSources: ['b.fa', 'c.fa', null], series
  }, 1);
  assert.deepEqual(kept, {
    blasts: ['b.tsv', 'd.tsv'], fastas: ['b.fa', 'd.fa'], sequenceSources: ['b.fa', null],
    series: [{ label: 'b', sourceIndex: 0, color: '#0' }, { label: 'd', sourceIndex: 1, color: '#2' }],
    labels: 'b,d', enabled: true
  });
  assert.equal(withoutConservationSeries({ blasts: ['b.tsv'], fastas: [], sequenceSources: [], series: [series[0]] }, 0).enabled, false,
    'with no series left the comparison is off');
  assert.equal(series[2].sourceIndex, 2, 'the stored entries are not changed');
});

// A row without a colour takes its colour by position, repeating after ten. The
// list stays as it was when it was derived by rounding (OV-2393 moved the colour
// conversion to truncation); each colour is one of the CLI's `ten` scheme.
test('a row without a colour takes the colour of its position', () => {
  const colors = ['#81a5cc', '#e5a972', '#7cb075', '#d38282', '#8fc4c0',
    '#ead487', '#be92b2', '#ffafb7', '#ae8e7c', '#c6bebb'];
  assert.deepEqual(colors.map((_, index) => normalizeConservationSeriesColor(null, index)), colors);
  assert.equal(normalizeConservationSeriesColor('', 10), colors[0]);
  assert.equal(normalizeConservationSeriesColor('#ABC', 3), '#aabbcc');
});

// PR-B: a ring's place holds its comparison genome's files. A place that is one
// file or null (a Session 50 or older) reads as a genome of that one file, or of
// none; a ring of one file keeps the key it had before a ring took several.
test('a ring place is its genome\'s files; the one-file form of an older Session reads as one', () => {
  const [a, b, c] = [file('a.fasta'), file('b.fasta', 20), file('c.fasta', 30)];
  assert.deepEqual(ringGenomeFileLists([a, null, [b, null, c]]), [[a], [], [b, c]]);
  assert.deepEqual(ringGenomeFileLists(undefined), []);
  const [oneFlat] = rowsFor([a]);
  const [oneList] = rowsFor([[a]]);
  assert.equal(oneList.sourceKey, oneFlat.sourceKey, 'a one-file ring keeps its key');
  const series = rowsFor([[b, c], [], [a]]);
  assert.deepEqual(series.map(({ fileName, sourceIndex }) => [fileName, sourceIndex]), [['b.fasta', 0], ['a.fasta', 1]],
    'an empty place is no ring; a ring is named by its first file');
  assert.notEqual(series[0].sourceKey, rowsFor([[c, b]])[0].sourceKey, 'the order of a genome\'s files is part of its key');
});

// PR-B review 1, items 1 and 5: Save, the span export catalog and the ring rows
// read each ring's genome files by the identity Generate matches (source key),
// and a place a Session left empty on Load (a missing file) keeps the later
// rings at theirs instead of lending its ring the next genome's files.
test('a ring keeps its own genome files when an earlier ring\'s files are missing', async () => {
  const { ringGenomeFilesInRingOrder, ringGenomePlaces } = await import('../../gbdraw/web/js/services/conservation-series.js');
  const [a, b, plasmid] = [file('a.fasta', 8), file('b.fasta', 14), file('b-plasmid.fasta', 4)];
  const [entryA, entryB] = rowsFor([[a], [b, plasmid]]);
  const drawing = (fastas) => ({
    files: { c_conservation_blasts: [], c_conservation_blasts_source: null, c_conservation_fastas: fastas, c_conservation_sequence_sources: [] },
    circularConservation: { source: 'losat', series: [entryA, entryB] }
  });
  assert.deepEqual(ringGenomeFilesInRingOrder(drawing([[a], [b, plasmid]])), [[a], [b, plasmid]]);
  assert.deepEqual(ringGenomeFilesInRingOrder(drawing([[], [b, plasmid]])), [[], [b, plasmid]],
    'the ring whose file is missing lists none; the next ring keeps its files');
  assert.deepEqual(ringGenomePlaces(drawing([[], [b, plasmid]])).places, [null, 1]);
  // Series moved: the rows follow their keys, not their places.
  const moved = { ...drawing([[a], [b, plasmid]]), circularConservation: { source: 'losat', series: [entryB, entryA] } };
  assert.deepEqual(ringGenomeFilesInRingOrder(moved), [[b, plasmid], [a]]);
  assert.deepEqual(ringGenomePlaces(moved).places, [1, 0]);
});
