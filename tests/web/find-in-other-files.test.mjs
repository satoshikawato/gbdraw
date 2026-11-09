// Find in other files (PD-OI-123): the dialog's owner over a fake Worker that answers as Python's steps answer
// (`gbdraw/web_support/locate.py`; each reply is checked against the generated description, so a reply here of another
// shape fails). The search runs under the LOSAT overlay by File; Cancel, Back and closing add nothing; Create adds the
// derived drawing through the drawing list as one step, after the source drawing.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createFindInOtherFiles } = await import('../../gbdraw/web/js/app/find-in-other-files.js');
const { referenceWindowOf, locateLosatJobs } = await import('../../gbdraw/web/js/services/find-in-other-files.js');
const { proposalView } = await import('../../gbdraw/web/js/app/find-in-other-files/proposal.js');

const record = (recordKey, recordId, length, circular = false) => ({ recordKey, recordId, length, circular });
const FILES = [
  { fileKey: 'ref', records: [record('ref', 'REF', 12000)] },
  { fileKey: 'a', records: [record('a', 'A', 9000)] },
  { fileKey: 'b', records: [record('b:1', 'B1', 5500), record('b:2', 'B2', 8500), record(null, 'B3', 700)] }
];
const sent = {
  envelope: {
    renderRequest: { mode: 'linear', records: [
      { recordKey: 'ref', source: { kind: 'genbank', resourceId: 'r0' } },
      { recordKey: 'a', source: { kind: 'genbank', resourceId: 'r1' } },
      { recordKey: 'b', source: { kind: 'gffFasta', gffResourceId: 'r2', fastaResourceId: 'r3' } }
    ] },
    resources: {}
  },
  filesByResource: new Map([['r0', { name: 'ref.gb' }], ['r1', { name: 'a.gb' }], ['r2', { name: 'b.gff' }]])
};

const cds = (index, start, end) => ({
  protein_id: `p${index}`, locus_tag: `g${index}`, gene: null, product: null, start, end, strand: '+', length_aa: 100, length_nt: 303, source_protein_id: null
});
const gene = (index) => ({ protein_id: `p${index}`, order: index, window_index: 0, start: 100 * index, end: 100 * index + 300, strand: '+', query_coverage: 1, bits: 200, status: 'complete' });
const planned = (recordKey, recordId, column, start, end, { selected = true, reason = null, reverse = false, genes = [1, 2], ends = [false, false] } = {}) => ({
  record_key: recordKey, record_id: recordId, record_index: 1, record_length: 9000, record_circular: false,
  window: { start, end, reverse_complement: reverse, column }, selected, reason, cluster: 0,
  touches_start: ends[0], touches_end: ends[1], across_origin: false, gap_to_next_same_record: null, evidence_bits: 400,
  genes: genes.map(gene), ref_span: [genes[0], genes.at(-1)]
});
const coverage = (complete, partial = [], absent = []) => ({ complete, partial, absent, reference_not_found: [] });
const tally = (read, kept) => ({ rows_read: read, rows_kept: kept, thresholds: [] });
const filePlan = (fileKey, windows, candidates, cover, warnings = [], counts = tally(10, 10)) => ({
  file_key: fileKey, windows, candidates, coverage: cover, tally: counts, warnings
});
const options = {
  program: 'tlosatn', thresholds: { evalue: 0.01, bitscore: 50, identity: 0, alignment_length: 0 }, best_ratio: 0.8,
  merge_distance: 10000, flank: 1000, min_genes: 2, min_share: 0.05, jump_min: 2, jump_share: 0.05, max_windows: 8
};
/** The search's plan: A has one window, B two and one not selected, the failed File none. */
const searchPlan = () => ({
  reference: {
    file_key: 'ref', record_key: 'ref', record_id: 'REF', record_length: 12000, circular: false,
    windows: [{ start: 1501, end: 8100, reverse_complement: false, column: 1 }],
    genes: [1, 2, 3].map((index) => ({ cds: cds(index, 100 * index, 100 * index + 300), window_index: 0, order: index })),
    cut_genes: [cds(7, 8000, 8400)]
  },
  options, runtime: null, warnings: [],
  files: [
    filePlan('a', [planned('a', 'A', 1, 1234, 7654, { reverse: true, genes: [1, 3], ends: [true, false] })], [], coverage(['p1', 'p3'], [], ['p2'])),
    filePlan('b', [planned('b:1', 'B1', 1, 1, 2000, { genes: [1] }), planned('b:2', 'B2', 2, 3000, 8000, { genes: [2, 3] })],
      [planned('b:2', 'B2', 3, 100, 900, { selected: false, reason: 'COVERED_BY_OTHER_WINDOW', genes: [3] })], coverage(['p1', 'p2', 'p3']))
  ]
});
const NAMES = { p1: 'g1', p2: 'g2', p3: 'g3', p7: 'g7' };
const drawingReply = (plan) => ({
  status: 'derived', plan, pairs: [],
  drawing: { status: 'derived', mode: 'linear', renderRequest: { mode: 'linear' }, resourceSources: {}, returnedToAuto: [],
    notCarried: [{ kind: 'LINEAR_COMPARISONS', count: 1, items: [] }] }
});
/** Python's `choose_locate_windows`, for the cases here: the kept windows by column, with the edited span. */
const chosenPlan = (plan, choices) => ({
  ...plan,
  files: plan.files.map((file) => {
    const all = [...file.windows, ...file.candidates];
    const kept = choices.filter((choice) => choice.fileKey === file.file_key).map((choice) => {
      const item = all.find((candidate) => candidate.window.column === choice.column);
      return { ...item, selected: true, reason: null,
        window: { start: choice.start, end: choice.end, reverse_complement: choice.reverseComplement, column: choice.column } };
    });
    return { ...file, windows: kept, candidates: all.filter((item) => !kept.some((k) => k.window.column === item.window.column)) };
  })
});

const harness = ({ losat = 'ok', files = FILES } = {}) => {
  const calls = { steps: [], losat: [], created: [] };
  const owner = createFindInOtherFiles({
    shownDrawingId: () => 'linear',
    projectShownRequest: async () => sent,
    runLocateStep: async ({ locate }) => {
      calls.steps.push(locate);
      const reply = {
        files: () => ({ files }),
        reference: () => (locate.reference.windows.some((window) => window.end > 12000)
          ? { status: 'refused', reason: 'WINDOW_PAST_END', start: locate.reference.windows[0].start, end: locate.reference.windows[0].end }
          : { status: 'ready', fileKey: 'ref', recordId: 'REF', genes: 3, cutGenes: ['g7'] }),
        search: () => ({ program: 'tlosatn', search: 'tblastn', cacheArgs: ['-evalue', '10'], jobs: locate.search.targets.map((fileKey) => ({
          fileKey, queryFasta: '>q\nMK\n', subjectFasta: `>${fileKey}\nACGT\n`, rawKey: `raw-${fileKey}` })) }),
        plan: () => ({ plan: searchPlan(), pairs: [], geneNames: NAMES }),
        drawing: () => {
          const edited = locate.drawing.choices.find((choice) => choice.start > choice.end);
          return edited ? { status: 'refused', reason: 'WINDOW_ACROSS_LINEAR_ORIGIN', start: edited.start, end: edited.end, fileKey: edited.fileKey, column: edited.column }
            : drawingReply(chosenPlan(locate.drawing.plan, locate.drawing.choices));
        }
      }[locate.step]();
      return { result: JSON.parse(JSON.stringify(reply)) };
    },
    runLosatJobs: async (plan) => {
      calls.losat.push(plan);
      if (losat === 'canceled') return { status: 'canceled' };
      return { status: 'ok', results: plan.jobs.map((job) => ({ ...job, text: 'rows' })) };
    },
    shownDrawing: () => ({ adv: { min_bitscore: '75', evalue: '', identity: null, alignment_length: 30 },
      linearSeqs: files.map(() => ({ gb: { name: 'x' }, gff: null, fasta: null })) }),
    newDrawingFromSearch: async (sourceId, derivation, kind) => {
      calls.created.push({ sourceId, derivation, kind });
      return { status: 'ok', drawingId: 'linear-2' };
    }
  });
  return { owner, calls };
};

test('the dialog opens on the drawing\'s Files with the picked span as the reference and counts its genes', async () => {
  const { owner, calls } = harness();
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  const { dialog } = owner;
  assert.equal(dialog.show, true);
  assert.deepEqual(dialog.rows.map(({ recordKey, start, end }) => [recordKey, start, end]), [['ref', '1501', '8100']]);
  assert.equal(dialog.genes, 3);
  // Search in: every other File, ticked at first, with its record count and length; the reference File is not offered.
  assert.deepEqual(owner.searchFiles().map((file) => [file.name, file.size, dialog.targets[file.fileKey]]), [
    ['a.gb', '1 record · 9,000 bp', true], ['b.gff', '3 records · 14,700 bp', true]
  ]);
  // Search settings start at Pairwise Comparisons' (a blank takes Linear's default).
  assert.deepEqual({ ...dialog.settings }, { program: 'tlosatn', bitscore: '75', evalue: '1e-2', identity: '0', alignmentLength: '30' });
  assert.equal(owner.searchReason(), '');
  assert.deepEqual(calls.steps.map((step) => step.step), ['files', 'reference']);
});

test('the Drawing menu opens one empty row, and Search says why it is disabled', async () => {
  const { owner } = harness();
  await owner.open();
  assert.deepEqual(owner.dialog.rows.map(({ recordKey, start, end }) => [recordKey, start, end]), [['', '', '']]);
  assert.equal(owner.searchReason(), 'Choose the record of each reference window.');
  const [row] = owner.dialog.rows;
  await owner.editRow(row.id, { recordKey: 'ref', start: '1,501', end: '99000' });
  assert.equal(owner.searchReason(), 'The window 1,501–99,000 reaches past the end of the record.');
  await owner.editRow(row.id, { end: '8100' });
  await owner.addRow();
  assert.equal(owner.dialog.rows[1].recordKey, 'ref');
  assert.equal(owner.searchReason(), 'Enter the Start and End of each reference window.');
  await owner.removeRow(owner.dialog.rows[1].id);
  for (const file of owner.searchFiles()) owner.setTarget(file.fileKey, false);
  assert.equal(owner.searchReason(), 'Tick at least one file to search in.');

  const single = harness({ files: FILES.slice(0, 1) });
  assert.equal(single.owner.entryReason(), 'Add another file to search in');
  await single.owner.open();
  assert.equal(single.owner.searchReason(), 'Add another file to search in');
});

test('a reference row reads as the record card reads its Region fields', () => {
  assert.deepEqual(referenceWindowOf({ start: '8,100', end: '1501', acrossOrigin: false }), { window: { start: 1501, end: 8100, reverseComplement: true } });
  assert.deepEqual(referenceWindowOf({ start: '8100', end: '1501', acrossOrigin: true }), { window: { start: 8100, end: 1501, reverseComplement: false } });
  assert.deepEqual(referenceWindowOf({ start: '100', end: '1501', acrossOrigin: true }), { reason: 'ORIGIN_NEEDS_START_AFTER_END' });
  assert.deepEqual(referenceWindowOf({ start: '1.5', end: '9', acrossOrigin: false }), { reason: 'NOT_WHOLE' });
});

test('a cancelled search adds nothing and keeps the dialog on the search', async () => {
  const { owner, calls } = harness({ losat: 'canceled' });
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  await owner.search();
  assert.equal(owner.dialog.screen, 'search');
  assert.equal(owner.dialog.status, 'ready');
  assert.equal(calls.losat.length, 1);
  assert.equal(calls.losat[0].byFile, true);
  assert.deepEqual(calls.steps.map((step) => step.step), ['files', 'reference', 'search']);
  assert.deepEqual(calls.created, []);
});

test('the search runs one job per File and proposes the windows Python planned, worded from its values', async () => {
  const { owner, calls } = harness();
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  owner.setSetting({ evalue: '0.001' });
  await owner.search();
  const [run] = calls.losat;
  assert.deepEqual(run.jobs.map((job) => [job.fileKey, job.program, job.cacheKey, job.pairIndex]), [
    ['a', 'tblastn', 'raw-a', 0], ['b', 'tblastn', 'raw-b', 1]
  ]);
  const planStep = calls.steps.find((step) => step.step === 'plan');
  assert.deepEqual(planStep.plan.options, { evalue: 0.001, bitscore: 75, identity: 0, alignmentLength: 30 });
  assert.deepEqual(planStep.plan.outcomes, [{ fileKey: 'a', text: 'rows', error: null }, { fileKey: 'b', text: 'rows', error: null }]);

  const { dialog } = owner;
  assert.equal(dialog.screen, 'proposal');
  assert.equal(dialog.reference, 'Reference: REF in ref.gb, 1 window, 3 genes.');
  assert.equal(dialog.cutGenes, 'Cut by a window edge and not searched: g7.');
  const [a, b] = dialog.proposal;
  assert.deepEqual([a.name, a.summary, b.name, b.summary], ['a.gb', '1 window · 2 of 3 genes', 'b.gff', '2 windows · 3 of 3 genes']);
  assert.deepEqual(a.windows.map(({ number, recordId, span, genes, ends }) => [number, recordId, span, genes, ends]), [
    ['①', 'A', '1,234–7,654 (−)', '2 · g1 to g3', ['contig start']]
  ]);
  assert.deepEqual([a.found, a.absent, a.partial], ['2 found in full, 0 in part, 1 not at all.', 'g2', '']);
  assert.deepEqual(b.windows.map((window) => window.number), ['①', '②']);
  assert.deepEqual(b.notSelected.map(({ ticked, number, reason }) => [ticked, number, reason]), [
    [false, '', 'Not selected: adds no gene that the selected windows do not cover.']
  ]);
  assert.equal(dialog.count, '3 windows in 2 files');
  assert.deepEqual(dialog.notCarried, ['Linear comparisons']);

  // Ticking the window that was not selected adds it in the order of the reference; unticking all disables Create.
  await owner.toggleWindow(b.notSelected[0].key, true);
  assert.deepEqual(owner.dialog.proposal[1].notSelected.map((window) => [window.ticked, window.number]), [[true, '③']]);
  assert.equal(owner.dialog.count, '4 windows in 2 files');
  for (const file of owner.dialog.proposal) {
    for (const window of [...file.windows, ...file.notSelected]) await owner.toggleWindow(window.key, false);
  }
  assert.equal(owner.dialog.createReason, 'Tick at least one window to create the drawing.');
  assert.equal(await owner.confirm(), null);
  assert.deepEqual(calls.created, []);
});

test('an edit changes only its window, and one Python refuses is undone with the reason', async () => {
  const { owner, calls } = harness();
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  await owner.search();
  const key = owner.dialog.proposal[1].windows[0].key;
  owner.startEdit(key);
  owner.editWindowField({ start: '200', reverseComplement: true });
  await owner.applyEdit();
  assert.equal(owner.dialog.edit, null);
  assert.equal(owner.dialog.proposal[1].windows[0].span, '200–2,000 (−)');
  assert.equal(owner.dialog.proposal[1].windows[1].span, '3,000–8,000');

  owner.startEdit(key);
  owner.editWindowField({ start: '3000', end: '100' });
  await owner.applyEdit();
  assert.equal(owner.dialog.edit.message, 'The window 3,000–100 crosses the origin, which a linear record does not have.');
  assert.equal(owner.dialog.proposal[1].windows[0].span, '200–2,000 (−)');
  const last = calls.steps.filter((step) => step.step === 'drawing').at(-1);
  assert.equal(last.drawing.choices.find((choice) => choice.fileKey === 'b' && choice.column === 1).start, 3000);
});

test('Back adds nothing; Create adds the derived drawing after the source drawing and closes', async () => {
  const { owner, calls } = harness();
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  await owner.search();
  owner.back();
  assert.equal(owner.dialog.screen, 'search');
  assert.deepEqual(owner.dialog.proposal, []);
  assert.deepEqual(calls.created, []);

  await owner.search();
  const drawingStep = calls.steps.filter((step) => step.step === 'drawing').at(-1);
  assert.deepEqual(drawingStep.drawing.labels, { ref: 'ref.gb', a: 'a.gb', b: 'b.gff' });
  const outcome = await owner.confirm();
  assert.deepEqual(outcome, { status: 'ok', drawingId: 'linear-2' });
  assert.equal(calls.created.length, 1);
  const [{ sourceId, derivation, kind }] = calls.created;
  assert.deepEqual([sourceId, kind, derivation.sent], ['linear', 'search', sent]);
  assert.equal(derivation.derived.mode, 'linear');
  assert.equal(owner.dialog.show, false);
});

test('a reply of another shape is refused, and a File whose search failed says why', async () => {
  const { owner } = harness({ files: [{ fileKey: 'ref' }] });
  await owner.open();
  assert.equal(owner.dialog.status, 'error');

  // A File whose search failed is listed with the reason and no windows; a File without a window says how many
  // alignments the thresholds removed.
  const plan = { ...searchPlan(), files: [
    filePlan('a', [], [], coverage([], [], ['p1', 'p2', 'p3']), [{ code: 'FILE_SEARCH_FAILED', args: { error: 'LOSAT exited with status 1' } }], tally(0, 0)),
    filePlan('b', [], [], coverage([], [], ['p1', 'p2', 'p3']), [{ code: 'FILE_NO_WINDOWS', args: { rows_read: 12, rows_kept: 0 } }], tally(12, 0))
  ] };
  const view = proposalView({ plan, pairs: [], geneNames: NAMES }, plan, [], { a: 'a.gb', b: 'b.gff' });
  assert.deepEqual(view.map((file) => [file.name, file.empty, file.windows.length]), [
    ['a.gb', 'The search of this file failed: LOSAT exited with status 1', 0],
    ['b.gff', 'No window found. The thresholds removed 12 alignments.', 0]
  ]);
  assert.deepEqual(locateLosatJobs({ program: 'tlosatn', search: 'tblastn', cacheArgs: [], jobs: [] }), { jobs: [], sequences: new Map() });
});
