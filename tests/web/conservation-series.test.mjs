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
  reconcileConservationSeries
} = await import('../../gbdraw/web/js/services/conservation-series.js');

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
