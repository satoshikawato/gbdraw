// Sequence search (PD-OI-117): a nucleotide query is searched with blastn and a protein query with tblastn against
// the same records. Python reads each query's program and reach (gbdraw/web_support/sequence_search.py); the panel
// plans one LOSAT job per program, and every LOSAT execution mode runs tblastn through the one program list of
// services/losat-runtime.js. The last case runs the vendored serial LOSAT, so it fails when the build lacks tblastn.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { test } from 'node:test';

const search = await import('../../gbdraw/web/js/services/sequence-search.js');
const runtime = await import('../../gbdraw/web/js/services/losat-runtime.js');
const generated = await import('../../gbdraw/web/js/utils/sequence-search-programs.generated.js');

const nucleotide = { name: 'n1', length: 40, program: 'blastn', reach: 40, fasta: '>n1\nACGT\n' };
const protein = { name: 'p1', length: 30, program: 'tblastn', reach: 90, fasta: '>p1\nMKV\n' };
const second = { name: 'n2', length: 20, program: 'blastn', reach: 20, fasta: '>n2\nTTTT\n' };

test('a search plans one LOSAT job per program, blastn with the panel task and tblastn with its default', () => {
  const plan = search.searchLosatPlan({ queries: [nucleotide, protein, second], subjectFasta: '>s\nACGT\n', task: 'megablast' });
  assert.deepEqual(plan.jobs, [
    { pairIndex: 0, cacheKey: 'sequence-search-blastn', program: 'blastn', querySequenceKey: 'search-query-blastn', subjectSequenceKey: 'search-subject', extraArgs: ['--task', 'megablast'] },
    { pairIndex: 1, cacheKey: 'sequence-search-tblastn', program: 'tblastn', querySequenceKey: 'search-query-tblastn', subjectSequenceKey: 'search-subject', extraArgs: [] }
  ]);
  assert.deepEqual([...plan.sequences], [
    ['search-subject', '>s\nACGT\n'],
    ['search-query-blastn', '>n1\nACGT\n>n2\nTTTT\n'],
    ['search-query-tblastn', '>p1\nMKV\n']
  ]);
  // A protein query pasted first changes nothing: blastn runs first, so its rows lead the raw table and its name.
  const pastedFirst = search.searchLosatPlan({ queries: [protein, nucleotide], subjectFasta: '>s\nA\n', task: 'blastn' });
  assert.deepEqual(pastedFirst.jobs.map(({ program }) => program), ['blastn', 'tblastn']);
  assert.deepEqual(search.searchPrograms([protein, second, nucleotide]), ['blastn', 'tblastn']);
  // Protein queries only: one tblastn job.
  assert.deepEqual(search.searchLosatPlan({ queries: [protein], subjectFasta: '>s\nA\n', task: 'blastn' }).jobs.map(({ program }) => program), ['tblastn']);
});

test('the origin window reaches as far as the query that covers the most bases: three per protein residue', () => {
  assert.deepEqual(search.searchSpans([{ recordKey: 'a' }, { recordKey: 'b' }], [nucleotide, protein]), [
    { recordKey: 'a', originWindow: 90 }, { recordKey: 'b', originWindow: 90 }
  ]);
});

test('the raw rows of the jobs are one table in job order', () => {
  assert.equal(search.searchRawText([{ text: 'a\t1\n' }, { text: '' }, { text: 'b\t2' }]), 'a\t1\nb\t2\n');
});

test('the programs and their order come from Python through the generated module', async () => {
  assert.deepEqual([...generated.SEARCH_PROGRAMS], ['blastn', 'tblastn']);
  const source = await readFile(new URL('../../gbdraw/web/js/services/sequence-search.js', import.meta.url), 'utf8');
  assert.doesNotMatch(source, /=== 'tblastn'/, 'the reply check keeps no program list of its own');
});

test('a query read without its program or reach is refused', () => {
  assert.deepEqual(search.readSearchQueriesReply({ queries: [nucleotide, protein] }), [nucleotide, protein]);
  for (const query of [{ ...protein, program: 'blastp' }, { ...protein, reach: undefined }, { name: 'q', length: 4, fasta: '>q\nACGT\n' }]) {
    assert.throws(() => search.readSearchQueriesReply({ queries: [query] }), (error) => error?.code === 'RESULT_INVALID');
  }
});

test('the panel names the programs of its queries, and the task applies to blastn only', async () => {
  globalThis.window ??= {
    Vue: {
      ref: (value) => ({ value }), reactive: (value) => value,
      computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
    }
  };
  const { createSequenceSearch } = await import('../../gbdraw/web/js/app/sequence-search.js');
  const sent = [];
  let reply = [nucleotide];
  const panel = createSequenceSearch({
    state: { activeDrawing: () => ({ id: 'linear', mode: 'linear', annotationSets: [] }) },
    readAxisRecords: () => [{ recordKey: 'r' }],
    readRecordSpans: async (spans) => {
      sent.push(spans);
      return { fasta: '>r:1-100\nACGT\n', spans: [] };
    },
    committedSets: () => [],
    readRecordBindings: () => ({ bindingKeyOf: () => 'r', selectorOf: () => null, selectorOfBinding: () => null }),
    runLosatJobs: async (plan) => {
      sent.push(plan.jobs.map(({ program }) => program));
      return { status: 'canceled' };
    },
    runHelper: async () => ({ result: { queries: reply } }),
    helperOperations: { readQueries: 'readSearchQueries', buildHitSet: 'searchHitAnnotationSet' },
    addSearchSet: () => {},
    shownRoot: () => null,
    closeRegionPopup: () => {},
    revealSection: async () => {}
  });
  assert.deepEqual([panel.programs.value, panel.taskApplies.value], [['blastn'], true]);
  reply = [protein];
  await panel.readQueries('>p1\nMKV\n');
  assert.deepEqual([panel.programs.value, panel.taskApplies.value, panel.programText.value], [['tblastn'], false, 'tblastn']);
  assert.equal(panel.queryLengthText(protein), '30 aa');
  reply = [protein, nucleotide];
  await panel.readQueries('>p1\nMKV\n>n1\nACGT\n');
  assert.deepEqual([panel.programs.value, panel.taskApplies.value, panel.programText.value], [['blastn', 'tblastn'], true, 'blastn, task megablast · tblastn']);
  await panel.search();
  assert.deepEqual(sent, [[{ recordKey: 'r', originWindow: 90 }], ['blastn', 'tblastn']]);
});

test('every LOSAT runner reads the one program list, which has tblastn', async () => {
  assert.deepEqual([...runtime.LOSAT_PROGRAMS], ['blastn', 'tblastx', 'blastp', 'tblastn']);
  assert.throws(() => runtime.requireLosatProgram('blastx'), /LOSAT program must be one of blastn, tblastx, blastp, tblastn\./);
  for (const path of ['services/losat.js', 'workers/losat-worker.js']) {
    const source = await readFile(new URL(`../../gbdraw/web/js/${path}`, import.meta.url), 'utf8');
    assert.match(source, /requireLosatProgram\(program\)/, path);
    assert.doesNotMatch(source, /new Set\(\['blastn'/, `${path} keeps no program list of its own`);
  }
});

test('the vendored serial LOSAT runs tblastn: a protein query hits its gene on the minus strand in genome coordinates', async () => {
  const web = new URL('../../gbdraw/web/', import.meta.url);
  const shim = await import(new URL('vendor/browser_wasi_shim/dist/index.js', web).href);
  const wasmModule = await WebAssembly.compile(await readFile(new URL('wasm/losat/losat.wasm', web)));
  // A 60-residue protein and its gene, reverse complemented, at bases 301-483 of a 900-base subject.
  const peptide = 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL';
  const CODONS = { M: 'ATG', S: 'TCT', K: 'AAA', G: 'GGT', E: 'GAA', L: 'CTG', F: 'TTC', T: 'ACC', V: 'GTT', P: 'CCG', I: 'ATC', D: 'GAT', N: 'AAC', H: 'CAT', A: 'GCT', Y: 'TAC', C: 'TGC', W: 'TGG' };
  const gene = [...peptide].map((residue) => CODONS[residue]).join('') + 'TAA';
  const complement = { A: 'T', C: 'G', G: 'C', T: 'A' };
  const reverse = [...gene].reverse().map((base) => complement[base]).join('');
  const filler = (length, seed) => Array.from({ length }, (_, index) => 'ACGT'[(index * seed + (index >> 3)) % 4]).join('');
  const subject = filler(300, 7) + reverse + filler(900 - 300 - reverse.length, 11);
  const text = await runtime.runLosatPairDirect({
    program: 'tblastn', queryFasta: `>prot\n${peptide}\n`, subjectFasta: `>rec:1-900\n${subject}\n`, wasiShim: shim, wasmModule
  });
  const [row] = text.trim().split('\n').map((line) => line.split('\t'));
  assert.deepEqual(row.slice(0, 2), ['prot', 'rec:1-900']);
  // qstart qend sstart send: the whole protein, read on the minus strand (sstart > send) over the gene's codons.
  assert.deepEqual(row.slice(6, 10).map(Number), [1, 60, 300 + reverse.length, 304]);
});

test('the hit table writes "aa" after the length of a protein query\'s hit; a nucleotide hit keeps its number', async () => {
  globalThis.window ??= {
    Vue: {
      ref: (value) => ({ value }), reactive: (value) => value,
      computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}
    }
  };
  const { createSearchHitTable } = await import('../../gbdraw/web/js/app/sequence-search/hit-table.js');
  const hit = (query, alignmentLength) => ({ annotationId: query, query, alignmentLength });
  const panel = {
    chip: '', selectedHit: '',
    outcome: {
      drawingId: 'linear', setId: null, rowsLeftOut: 0, recordIds: {}, rawText: '',
      hits: [hit('n1', 1234), hit('p1', 237)],
      queries: [{ ...nucleotide, color: '#000' }, { ...protein, color: '#111' }]
    }
  };
  const table = createSearchHitTable(panel, { activeDrawing: () => ({ id: 'linear', annotationSets: [] }), shownRoot: () => null });
  assert.deepEqual(table.hitRows.value.map((row) => table.hitLengthText(row)), ['1234', '237 aa']);
});
