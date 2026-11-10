// D12: a Circular similarity ring row added from a GenBank or DDBJ file takes
// the label the file names itself (first record's DEFINITION, then organism),
// as read by the Python ring reader; FASTA rows keep the file-name default and a
// typed label always wins.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

const {
  applyComparisonSequenceRecordLabel,
  conservationSourceDescriptors,
  defaultConservationSeriesLabel,
  effectiveRingColor,
  moveConservationSeriesEntry,
  orderedConservationSources,
  reassignDrawnRingColors,
  reassignRingColorsOnReorder,
  reconcileConservationSeries,
  ringDrawnPositions,
  ringKeysOfFiles,
  ringsOfDrawing,
  ringSwatchColor,
  withoutConservationSeries
} = await import('../../gbdraw/web/js/services/conservation-series.js');
const { ringGenomeFileLists } = await import('../../gbdraw/web/js/services/drawing/inputs.js');

const file = (name, size = 10) => ({ name, size, lastModified: 0 });

const orderedConservationSourcesOf = (sourceFiles, series) => orderedConservationSources(sourceFiles, { series })
  .map((source) => [source.label, source.files.reduce((total, member) => total + member.size, 0)]);
const rowsFor = (sourceFiles) => reconcileConservationSeries({
  sourceFiles,
  previousSeries: [],
  legacyLabels: []
}).series;

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
    enabled: true
  });
  assert.equal(withoutConservationSeries({ blasts: ['b.tsv'], fastas: [], sequenceSources: [], series: [series[0]] }, 0).enabled, false,
    'with no series left the comparison is off');
  assert.equal(series[2].sourceIndex, 2, 'the stored entries are not changed');
});

// PD-OI-110 (Ring colors): a ring without a colour of its own takes its
// scheme's colour at its scheme index, repeating; Same color leaves it to the
// configured colour (null, sent as ''); a hand-set colour wins on every scheme.
test('a ring takes its own colour on every scheme, else the scheme colour at its index', () => {
  const ten = ['#81a5cc', '#e5a972', '#8fc4c0', '#d38282', '#c6bebb', '#7cb075', '#be92b2', '#ead487', '#ae8e7c', '#ffafb7'];
  assert.deepEqual(ten.map((_, schemeIndex) => effectiveRingColor({ color: null, schemeIndex }, 'ten')), ten,
    'the ten colours in the receipt order (no Web order)');
  assert.equal(effectiveRingColor({ color: null, schemeIndex: 10 }, 'ten'), ten[0], 'ring 11 repeats ring 1');
  assert.equal(effectiveRingColor({ color: null, schemeIndex: 4 }, 'per_five'), '#8b9cc1');
  assert.equal(effectiveRingColor({ color: null, schemeIndex: 5 }, 'per_five'), '#c88890');
  assert.equal(effectiveRingColor({ color: null, schemeIndex: 30 }, 'per_five'), '#8b9cc1', 'ring 31 repeats ring 1');
  assert.equal(effectiveRingColor({ color: null, schemeIndex: 3 }, 'same'), null, 'Same color invents no colour');
  assert.equal(effectiveRingColor({ color: '', schemeIndex: 3 }, undefined), null, 'unset reads as Same color');
  for (const scheme of ['same', 'ten', 'per_five']) {
    assert.equal(effectiveRingColor({ color: '#ABC', schemeIndex: 0 }, scheme), '#aabbcc', `${scheme}: the own colour wins`);
  }
  const packaged = { unmanagedConfigOverrides: {} };
  assert.equal(ringSwatchColor(null, packaged), '#8b9cc1', 'a Same color ring shows the configured colour');
  assert.equal(ringSwatchColor('', packaged), '#8b9cc1');
  assert.equal(ringSwatchColor('#123456', packaged), '#123456');
  const configured = { unmanagedConfigOverrides: { 'objects.conservation.max_color': '#AB0' } };
  assert.equal(ringSwatchColor(null, configured), '#aabb00', 'the drawing\'s own max colour, as Python draws it');
  assert.equal(ringSwatchColor('#123456', configured), '#123456');
});

test('a new ring has no colour of its own and takes the first scheme index no ring uses (D-05)', () => {
  const [a, b, c, d] = ['a.fa', 'b.fa', 'c.fa', 'd.fa'].map((name) => file(name));
  const three = rowsFor([a, b, c]);
  assert.deepEqual(three.map(({ color, schemeIndex }) => [color, schemeIndex]), [[null, 0], [null, 1], [null, 2]]);
  const withoutB = reconcileConservationSeries({ sourceFiles: [a, c], previousSeries: three, legacyLabels: [] }).series;
  assert.deepEqual(withoutB.map(({ schemeIndex }) => schemeIndex), [0, 2], 'a removal keeps the others\' colours');
  const added = reconcileConservationSeries({ sourceFiles: [a, c, d], previousSeries: withoutB, legacyLabels: [] }).series;
  assert.deepEqual(added.map(({ schemeIndex }) => schemeIndex), [0, 2, 1], 'the added ring takes the freed index');
  assert.deepEqual(orderedConservationSources([a, c, d], { series: added, colorScheme: 'ten' }).map(({ color }) => color),
    ['#81a5cc', '#8fc4c0', '#e5a972']);
  assert.deepEqual(orderedConservationSources([a, c, d], { series: added }).map(({ color }) => color), ['', '', ''],
    'Same color sends no colour, as the CLI reads it');
});

test('a hand-set colour survives a scheme switch and a reorder; re-assignment changes scheme colours only', () => {
  const [a, b, c] = ['a.fa', 'b.fa', 'c.fa'].map((name) => file(name));
  const series = rowsFor([a, b, c]);
  series[1].color = '#123456';
  const colors = (scheme) => orderedConservationSources([a, b, c], { series, colorScheme: scheme }).map(({ color }) => color);
  assert.deepEqual(colors('ten'), ['#81a5cc', '#123456', '#8fc4c0']);
  assert.deepEqual(colors('same'), ['', '#123456', '']);
  assert.ok(moveConservationSeriesEntry(series, 2, -1) && moveConservationSeriesEntry(series, 1, -1), 'c moves to the outside');
  assert.deepEqual(series.map(({ fileName }) => fileName), ['c.fa', 'a.fa', 'b.fa']);
  assert.deepEqual(orderedConservationSources([a, b, c], { series, colorScheme: 'ten' }).map(({ color }) => color),
    ['#8fc4c0', '#81a5cc', '#123456'], 'with re-assignment off a colour stays with its genome');
  const drawing = { circularConservation: { series, colorScheme: 'ten' }, adv: { circular_track_slots_enabled: false } };
  assert.equal(reassignDrawnRingColors(drawing), true);
  assert.deepEqual(series.map(({ color, schemeIndex }) => [color, schemeIndex]), [[null, 0], [null, 1], ['#123456', 1]],
    'the button gives each scheme ring its position and keeps the hand-set colour');
  assert.equal(reassignDrawnRingColors(drawing), false, 'a second press changes nothing');
});

test('re-assignment on reorder runs only with the checkbox on, and not under Same color (D-21)', () => {
  const ring = (schemeIndex) => ({ color: null, schemeIndex });
  const drawing = (colorScheme, reassignOnReorder) => ({
    circularConservation: { series: [ring(1), ring(0)], colorScheme, reassignOnReorder },
    adv: { circular_track_slots_enabled: false }
  });
  const indexes = (shown) => shown.circularConservation.series.map(({ schemeIndex }) => schemeIndex);
  for (const [scheme, on, changed] of [['ten', true, true], ['per_five', true, true], ['ten', false, false],
    ['same', true, false], [undefined, true, false]]) {
    const shown = drawing(scheme, on);
    assert.equal(reassignRingColorsOnReorder(shown), changed, `${scheme}, checkbox ${on}`);
    assert.deepEqual(indexes(shown), changed ? [0, 1] : [1, 0]);
  }
  const same = drawing('same', false);
  assert.equal(reassignDrawnRingColors(same), false, 'the button changes nothing under Same color');
  assert.deepEqual(indexes(same), [1, 0]);
});

test('the drawn order is the stack order of the enabled ring slots, as Python counts it', () => {
  const slot = (trackIndex, enabled = true) => ({ renderer: 'sequence_conservation', enabled, params: { track_index: String(trackIndex) } });
  assert.deepEqual(ringDrawnPositions(3, null), [0, 1, 2]);
  assert.deepEqual(ringDrawnPositions(3, [{ renderer: 'features' }, slot(3), slot(1), slot(2)]), [1, 2, 0]);
  assert.deepEqual(ringDrawnPositions(3, [slot(2), slot(1, false)]), [1, 0, 2], 'a ring no slot draws follows the drawn ones');
});

// PD-OI-110 (D-38): Load pairs a draft entry with a ring of the stored request by the key the ring's file has among
// the request's ring files: the key the draft's entry for that file holds.
test('the ring keys of a request are its ring files\' source keys, none for a ring without a file', () => {
  const file = (name, size = 4) => ({ name, size, lastModified: 7 });
  const ringFiles = [file('a.fna'), null, file('a.gbk'), file('a.fna')];
  assert.deepEqual(ringKeysOfFiles(ringFiles), ['a.fna|4|7|0', null, 'a.gbk|4|7|0', 'a.fna|4|7|1']);
  assert.deepEqual(ringKeysOfFiles(ringFiles).filter(Boolean),
    conservationSourceDescriptors(ringFiles).map(({ sourceKey }) => sourceKey), 'the keys the draft\'s entries take');
});

// OV-3236: each ring takes the file its entry's `sourceKey` names, also among files of one name (by file name, two
// different files of one name took the files in file order whatever the rings' order).
test('the ring files in series order follow each entry\'s source key, also for files of one name', () => {
  const x = { name: 'genomic.fna', size: 8, lastModified: 1 };
  const y = { name: 'genomic.fna', size: 14, lastModified: 2 };
  const series = reconcileConservationSeries({ sourceFiles: [x, y], previousSeries: [], legacyLabels: [] }).series;
  assert.ok(moveConservationSeriesEntry(series, 1, -1));
  assert.deepEqual(orderedConservationSources([x, y], { series }).map(({ file }) => file), [y, x]);
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
  assert.deepEqual(series.map(({ fileName, sourceIndex }) => [fileName, sourceIndex]), [['b.fasta', 0], ['a.fasta', 2]],
    'an empty place is no ring; a ring is named by its first file and keeps its place');
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

// PR-B review 2, HIGH 1 (b): the one matcher works in passes over every entry (source keys first, then file names
// among the sources no entry claimed, then the stored place), so an entry whose files are missing does not take a
// later ring's genome whose first file has the same name.
test('a ring whose files are missing does not take a later ring of the same file name', async () => {
  const { orderedConservationSources, ringGenomeFilesInRingOrder, ringGenomePlaces } = await import('../../gbdraw/web/js/services/conservation-series.js');
  const [a, b] = [file('sequence.fasta', 8), file('sequence.fasta', 16)];
  const [entryA, entryB] = rowsFor([[a], [b]]).map((entry, index) => ({ ...entry, label: `Genome ${'AB'[index]}` }));
  const drawing = {
    files: { c_conservation_blasts: [], c_conservation_blasts_source: null, c_conservation_fastas: [[], [b]], c_conservation_sequence_sources: [] },
    circularConservation: { source: 'losat', series: [entryA, entryB] }
  };
  assert.deepEqual(ringGenomePlaces(drawing).places, [null, 1]);
  assert.deepEqual(ringGenomeFilesInRingOrder(drawing), [[], [b]], 'Save keeps the missing ring\'s place empty');
  assert.deepEqual(orderedConservationSources(drawing.files.c_conservation_fastas, drawing.circularConservation)
    .map((source) => [source.label, source.ownColor, source.schemeIndex, source.placeIndex, source.files]),
  [['Genome B', entryB.color, entryB.schemeIndex, 1, [b]]],
  'Generate draws Genome B on its own genome');
  assert.deepEqual(reconcileConservationSeries({ sourceFiles: [[b]], previousSeries: [entryA, entryB], legacyLabels: [] }).series
    .map((entry) => [entry.label, entry.color, entry.schemeIndex]), [['Genome B', entryB.color, entryB.schemeIndex]],
  'the reconcile keeps Genome B, not A on B\'s file');
});

// The series reconcile follows a change of the files: an entry whose files went leaves with them, so a re-picked set
// of uploaded tables starts with rows of their own, not the old rows' labels by place (as on dev). Generate, Save and
// the span export also read an entry by the place a Session or a LOSAT-cache replay stored (the last pass).
test('a re-picked set of tables does not keep the old rows\' labels by place', () => {
  const first = rowsFor([file('a.tsv', 5), file('b.tsv', 6)]);
  first[0].label = 'Alpha';
  const next = reconcileConservationSeries({ sourceFiles: [file('c.tsv', 7), file('d.tsv', 8)], previousSeries: first, legacyLabels: [] }).series;
  assert.deepEqual(next.map((entry) => [entry.fileName, entry.label]), [['c.tsv', 'c'], ['d.tsv', 'd']]);
});

// #1777 review 3, finding 1: the name pass claims a source only when its name is unambiguous among the unclaimed
// sources and the entries left; two rings whose tables (a LOSAT-cache replay: `<first FASTA stem>.circular_conservation.
// <suffix>.tsv`) or first files share a name are told apart by their stored places, which a move keeps.
test('two rings of one table name keep their own tables and genomes after a move', async () => {
  const { orderedConservationSources, ringGenomeFilesInRingOrder, ringGenomePlaces, moveConservationSeriesEntry } =
    await import('../../gbdraw/web/js/services/conservation-series.js');
  const shared = 'sequence.circular_conservation.losatn.tsv';
  const tables = [file(shared, 60), file(shared, 63), file('d.circular_conservation.losatn.tsv', 64)];
  const genomes = [[file('b.fasta', 1)], [file('c.fasta', 2)], [file('d.fasta', 3)]];
  // The entries name the tables and keep the keys of the genomes they searched.
  const series = rowsFor(genomes).map((entry, index) => ({ ...entry, fileName: tables[index].name, label: `Genome ${'BCD'[index]}` }));
  const drawing = {
    files: { c_conservation_blasts: tables, c_conservation_blasts_source: 'losat-cache', c_conservation_fastas: genomes, c_conservation_sequence_sources: [] },
    circularConservation: { source: 'losat', series }
  };
  moveConservationSeriesEntry(series, 0, 1);
  assert.deepEqual(orderedConservationSources(tables, drawing.circularConservation).map((source) => [source.label, source.file.size]),
    [['Genome C', 63], ['Genome B', 60], ['Genome D', 64]], 'Generate draws each ring on its own table');
  assert.deepEqual(ringGenomePlaces(drawing).places, [1, 0, 2]);
  assert.deepEqual(ringGenomeFilesInRingOrder(drawing), [genomes[1], genomes[0], genomes[2]], 'Save writes each ring\'s own genome');
});

test('two LOSAT rings whose first files share a name keep their own genomes when their keys no longer match', () => {
  const genomes = [[file('sequence.fasta', 8)], [file('sequence.fasta', 16)], [file('sequence.fasta', 24)]];
  // Keys of other dates (a Session whose bindings lacked them): only the names and the stored places are left.
  const series = rowsFor(genomes).map((entry, index) => ({ ...entry, sourceKey: `stale-${index}`, label: `Genome ${'ABC'[index]}` }));
  const moved = [series[1], series[0], series[2]];
  assert.deepEqual(orderedConservationSourcesOf(genomes, moved), [['Genome B', 16], ['Genome A', 8], ['Genome C', 24]]);
  assert.deepEqual(reconcileConservationSeries({ sourceFiles: genomes, previousSeries: moved, legacyLabels: [] }).series
    .map((entry) => [entry.label, entry.color, entry.sourceIndex]), moved.map((entry) => [entry.label, entry.color, entry.sourceIndex]),
  'the reconcile keeps each row on the place of its own name');
});

// #1777 review 4: a source key ends with the count of earlier sources of the same files, which shifts when one of them
// goes. Rings whose files share name, size and date (three NCBI `sequence.fasta` downloads of one size) are therefore
// told apart by their stored places, as rings of one name are: the key pass claims a key only when one source and one
// entry carry its files, so Remove series and a Load that lost one ring's file leave each other ring its own genome.
test('rings whose files share name, size and date keep their own genomes when one of them goes', () => {
  const genomes = [[file('sequence.fasta', 8)], [file('sequence.fasta', 8)], [file('sequence.fasta', 8)]];
  const series = rowsFor(genomes).map((entry, index) => ({ ...entry, label: `Genome ${'ABC'[index]}` }));
  assert.deepEqual(series.map((entry) => [entry.sourceKey, entry.sourceIndex]),
    [['sequence.fasta|8|0|0', 0], ['sequence.fasta|8|0|1', 1], ['sequence.fasta|8|0|2', 2]]);
  /** Which genome each drawn ring takes. */
  const drawn = (places, entries) => orderedConservationSources(places, { series: entries })
    .map((source) => [source.label, genomes.findIndex(([member]) => member === source.file)]);
  // Remove series of Genome A: the places after it move up one, and so do their entries' stored places.
  const removedFirst = [series[1], series[2]].map((entry) => ({ ...entry, sourceIndex: entry.sourceIndex - 1 }));
  assert.deepEqual(drawn(genomes.slice(1), removedFirst), [['Genome B', 1], ['Genome C', 2]]);
  // A Load that lost Genome A's file keeps its place empty: B and C keep theirs, and the reconcile keeps their places.
  assert.deepEqual(drawn([[], genomes[1], genomes[2]], series), [['Genome B', 1], ['Genome C', 2]]);
  assert.deepEqual(reconcileConservationSeries({ sourceFiles: [[], genomes[1], genomes[2]], previousSeries: series, legacyLabels: [] }).series
    .map((entry) => [entry.label, entry.sourceIndex]), [['Genome B', 1], ['Genome C', 2]]);
  // One source left of two entries of the same files: the entry whose file went takes none.
  assert.deepEqual(drawn([[], genomes[1]], series.slice(0, 2)), [['Genome B', 1]]);
});

// #1777 review 5: a Generate tells whether the draft's rings are still those its searches read by the rings' files, not
// by their keys alone: after the first of three rings of `sequence.fasta` files of one size and date goes and a new
// one comes, the two left and the new one carry the keys the three had.
test('the rings of a drawing differ when rings of files of one name, size and date change, though their keys do not', () => {
  const genome = (fileId) => ({ fileId, sourceId: fileId, name: 'sequence.fasta', type: 'text/plain', size: 8, lastModified: 0 });
  const drawingOf = (fileIds) => {
    const lists = fileIds.map((fileId) => [genome(fileId)]);
    const series = reconcileConservationSeries({ sourceFiles: lists, previousSeries: [], legacyLabels: [] }).series;
    return {
      files: { c_conservation_blasts: [], c_conservation_blasts_source: null, c_conservation_fastas: lists, c_conservation_sequence_sources: [] },
      circularConservation: { source: 'losat', series }
    };
  };
  const searched = drawingOf(['a', 'b', 'c']);
  const changed = drawingOf(['b', 'c', 'd']);
  assert.deepEqual(changed.circularConservation.series.map((entry) => entry.sourceKey), searched.circularConservation.series.map((entry) => entry.sourceKey));
  assert.notEqual(ringsOfDrawing(changed), ringsOfDrawing(searched));
  assert.equal(ringsOfDrawing(drawingOf(['a', 'b', 'c'])), ringsOfDrawing(searched));
});
