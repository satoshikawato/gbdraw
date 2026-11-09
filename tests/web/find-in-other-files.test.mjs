// Find in other files (PD-OI-123): the dialog's owner over a fake Worker that answers as Python's steps answer
// (`gbdraw/web_support/locate.py`; each reply is checked against the generated description, so a reply here of another
// shape fails). The search runs under the full-window overlay from its first Python step to its last, LOSAT by File;
// Cancel, Back and closing add nothing; Create adds the derived drawing through the drawing list as one step, after
// the source drawing.
import assert from 'node:assert/strict';
import test from 'node:test';

const { createFindInOtherFiles } = await import('../../gbdraw/web/js/app/find-in-other-files.js');
const { referenceWindowOf, locateLosatJobs } = await import('../../gbdraw/web/js/services/find-in-other-files.js');
const { proposalView, windowDrafts } = await import('../../gbdraw/web/js/app/find-in-other-files/proposal.js');

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
/**
 * Python's `locate_file_counts` (gbdraw/comparisons/locate.py), step by step: the drawn row places a File's kept
 * windows, in plan order, at columns 1..k (`_drawn_row`), the drawing numbers a window by that column, and the reply
 * gives each number by the window's column in the plan (its identity there, windows not kept counted), in that order;
 * then the genes found and the rows the thresholds removed.
 */
const countsOf = (plan) => plan.files.map((file) => ({
  fileKey: file.file_key,
  windowNumbers: file.windows.map((item, index) => ({ column: item.window.column, number: index + 1 }))
    .sort((a, b) => a.column - b.column),
  genesFound: file.coverage.complete.length + file.coverage.partial.length,
  genesTotal: plan.reference.genes.length,
  alignmentsRemoved: file.tally.rows_read - file.tally.rows_kept
}));
const drawingReply = (plan, counts = countsOf(plan)) => ({
  status: 'derived', plan, counts, pairs: [],
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

/** A promise and its resolve, for a step a test answers later. */
const deferred = () => {
  let resolve;
  const promise = new Promise((done) => { resolve = done; });
  return { promise, resolve };
};

/**
 * The owner over fakes. The overlay fake is `runInLosatOverlay`'s contract: it is up while the task runs, `cancel()`
 * aborts it, `until` drops a reply that comes after Cancel, and a run it refuses (`busy`) never starts the task.
 * `held` names a Python step whose reply waits for `calls.release()`; `counts` replaces Python's counts.
 */
const harness = ({ losat = 'ok', files = FILES, busy = '', held = '', counts = null, drawingAdv = null } = {}) => {
  const calls = { steps: [], losat: [], created: [], status: [], overlay: false, cancel: () => {}, release: () => {} };
  const owner = createFindInOtherFiles({
    shownDrawingId: () => 'linear',
    projectShownRequest: async () => sent,
    runLocateStep: async ({ locate }) => {
      calls.steps.push({ ...locate, overlay: calls.overlay, dialogStatus: owner.dialog.status, dialogShown: owner.dialog.show });
      if (locate.step === held) {
        const wait = deferred();
        calls.release = wait.resolve;
        await wait.promise;
      }
      const reply = {
        files: () => ({ files }),
        reference: () => (locate.reference.windows.some((window) => window.end > 12000)
          ? { status: 'refused', reason: 'WINDOW_PAST_END', start: locate.reference.windows[0].start, end: locate.reference.windows[0].end }
          : { status: 'ready', fileKey: 'ref', recordId: 'REF', genes: 3 }),
        search: () => ({ program: 'tlosatn', search: 'tblastn', cacheArgs: ['-evalue', '10'], jobs: locate.search.targets.map((fileKey) => ({
          fileKey, queryFasta: '>q\nMK\n', subjectFasta: `>${fileKey}\nACGT\n`, rawKey: `raw-${fileKey}` })) }),
        plan: () => ({ plan: searchPlan(), pairs: [], geneNames: NAMES }),
        drawing: () => {
          const edited = locate.drawing.choices.find((choice) => choice.start > choice.end);
          return edited ? { status: 'refused', reason: 'WINDOW_ACROSS_LINEAR_ORIGIN', start: edited.start, end: edited.end, fileKey: edited.fileKey, column: edited.column }
            : drawingReply(...((chosen) => [chosen, counts ?? countsOf(chosen)])(chosenPlan(locate.drawing.plan, locate.drawing.choices)));
        }
      }[locate.step]();
      return { result: JSON.parse(JSON.stringify(reply)) };
    },
    runInOverlay: async (task) => {
      if (busy) return { status: 'busy', reason: busy };
      const controller = new AbortController();
      const canceled = new Promise((_resolve, reject) => controller.signal.addEventListener('abort', () => reject(controller.signal.reason)));
      canceled.catch(() => {});
      calls.cancel = () => controller.abort(new Error('canceled'));
      calls.overlay = true;
      try {
        return { status: 'ok', value: await task({
          signal: controller.signal,
          status: (text) => calls.status.push(text),
          until: (step) => Promise.race([step, canceled]),
          losat: async (plan) => {
            calls.losat.push(plan);
            if (losat === 'canceled') calls.cancel();
            await Promise.resolve();
            if (controller.signal.aborted) throw controller.signal.reason;
            return plan.jobs.map((job) => ({ ...job, text: 'rows' }));
          }
        }) };
      } catch (error) {
        if (controller.signal.aborted) return { status: 'canceled' };
        throw error;
      } finally {
        calls.overlay = false;
      }
    },
    shownDrawing: () => ({ adv: drawingAdv ?? { min_bitscore: '75', evalue: '', identity: null, alignment_length: 30 },
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

test('a span picked across the origin of a circular record opens as a window across the origin', async () => {
  // The Region popup gives the span in source coordinates on the source strand, so Start > End is a span across the
  // origin (a window of a record drawn in several windows included); the row reads it so, not as reverse-complemented.
  const { owner, calls } = harness({ files: [{ fileKey: 'ref', records: [record('ref', 'REF', 12000, true)] }, ...FILES.slice(1)] });
  await owner.open({ recordKey: 'ref', start: 11001, end: 500 });
  assert.deepEqual(owner.dialog.rows.map(({ recordKey, start, end, acrossOrigin }) => [recordKey, start, end, acrossOrigin]),
    [['ref', '11001', '500', true]]);
  assert.deepEqual(calls.steps.at(-1).reference.windows, [{ start: 11001, end: 500, reverseComplement: false }]);
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
  const view = proposalView({ plan, pairs: [], geneNames: NAMES }, { plan, counts: countsOf(plan) }, [], { a: 'a.gb', b: 'b.gff' });
  assert.deepEqual(view.map((file) => [file.name, file.empty, file.windows.length]), [
    ['a.gb', 'The search of this file failed: LOSAT exited with status 1', 0],
    ['b.gff', 'No window found. The thresholds removed 12 alignments.', 0]
  ]);
  assert.deepEqual(locateLosatJobs({ program: 'tlosatn', search: 'tblastn', cacheArgs: [], jobs: [] }), { jobs: [], sequences: new Map() });
});

test('the overlay is up from the first Python step of a search to its last, and the dialog steps aside only then', async () => {
  const { owner, calls } = harness();
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  assert.equal(calls.overlay, false);
  await owner.search();
  // Python reads the Files and plans under the overlay, with the dialog hidden ('searching') and still open.
  const during = calls.steps.filter((step) => step.step === 'search' || step.step === 'plan');
  assert.deepEqual(during.map(({ step, overlay, dialogStatus, dialogShown }) => [step, overlay, dialogStatus, dialogShown]), [
    ['search', true, 'searching', true], ['plan', true, 'searching', true]
  ]);
  // The overlay says what runs at each step; LOSAT's "File k of N" is the by-File runner's (services/losat.js).
  assert.deepEqual(calls.status, ['Reading the files of the drawing…', 'Finding the windows in each file…']);
  // The proposal is worked out after the overlay is down, in the dialog.
  const drawing = calls.steps.find((step) => step.step === 'drawing');
  assert.deepEqual([drawing.overlay, owner.dialog.screen], [false, 'proposal']);
});

test('Cancel during a Python step of the search discards its reply and keeps the dialog on the search', async () => {
  for (const held of ['search', 'plan']) {
    const { owner, calls } = harness({ held });
    await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
    const searching = owner.search();
    await new Promise((resolve) => setTimeout(resolve, 0));
    assert.equal(calls.overlay, true, held);
    calls.cancel();
    await searching;
    assert.equal(calls.overlay, false, held);
    assert.deepEqual([owner.dialog.screen, owner.dialog.status, owner.dialog.message], ['search', 'ready', ''], held);
    // The step's reply comes after Cancel and changes nothing; the search adds nothing.
    calls.release();
    await new Promise((resolve) => setTimeout(resolve, 0));
    assert.deepEqual([owner.dialog.screen, owner.dialog.status], ['search', 'ready'], held);
    assert.deepEqual(calls.steps.map((step) => step.step), held === 'search' ? ['files', 'reference', 'search'] : ['files', 'reference', 'search', 'plan'], held);
    assert.equal(calls.losat.length, held === 'search' ? 0 : 1, held);
    assert.deepEqual(calls.created, [], held);
  }
});

test('a search while a Generate or a Session operation runs does not start, and says why', async () => {
  const { owner, calls } = harness({ busy: 'Generating diagram. Retry after generation finishes.' });
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  await owner.search();
  assert.deepEqual([owner.dialog.screen, owner.dialog.status, owner.dialog.message],
    ['search', 'ready', 'Generating diagram. Retry after generation finishes.']);
  assert.deepEqual(calls.steps.map((step) => step.step), ['files', 'reference']);
});

test('the window numbers, the genes found and the alignments removed are the ones Python gives', async () => {
  // Python numbers the kept windows as the drawing numbers them; the dialog shows its numbers, not their positions.
  const counts = [
    { fileKey: 'b', windowNumbers: [{ column: 1, number: 2 }, { column: 2, number: 1 }], genesFound: 2, genesTotal: 5, alignmentsRemoved: 0 },
    { fileKey: 'a', windowNumbers: [{ column: 1, number: 1 }], genesFound: 1, genesTotal: 5, alignmentsRemoved: 7 }
  ];
  const { owner } = harness({ counts });
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  await owner.search();
  const [a, b] = owner.dialog.proposal;
  assert.deepEqual(b.windows.map((window) => [window.recordId, window.number]), [['B1', '②'], ['B2', '①']]);
  // Each File is found by its key (the reply lists them in another order here).
  assert.deepEqual([a.summary, b.summary], ['1 window · 1 of 5 genes', '2 windows · 2 of 5 genes']);

  const plan = { ...searchPlan(), files: [filePlan('a', [], [], coverage([], [], ['p1', 'p2', 'p3']), [], tally(12, 0))] };
  const view = proposalView({ plan, pairs: [], geneNames: NAMES },
    { plan, counts: [{ fileKey: 'a', windowNumbers: [], genesFound: 0, genesTotal: 3, alignmentsRemoved: 9 }] }, [], { a: 'a.gb' });
  assert.equal(view[0].empty, 'No window found. The thresholds removed 9 alignments.');

  // A kept window without its number in the reply is not the reply of this plan: refused, never shown unnumbered.
  const kept = searchPlan();
  const missing = countsOf(kept).map((item) => (item.fileKey === 'b' ? { ...item, windowNumbers: item.windowNumbers.slice(1) } : item));
  assert.throws(() => proposalView({ plan: kept, pairs: [], geneNames: NAMES }, { plan: kept, counts: missing }, windowDrafts(kept), {}),
    (error) => error.code === 'RESULT_INVALID' && error.operation === 'findInOtherFiles');
});

test('each search setting starts from its own Pairwise Comparisons value; an invalid one takes its default alone', async () => {
  const { owner } = harness({ drawingAdv: { min_bitscore: '100', evalue: 'abc', identity: '', alignment_length: '-5' } });
  await owner.open();
  assert.deepEqual({ ...owner.dialog.settings }, { program: 'tlosatn', bitscore: '100', evalue: '1e-2', identity: '0', alignmentLength: '0' });
});

test('a tick or an edit while Python answers the last one is not taken', async () => {
  const { owner, calls } = harness();
  await owner.open({ recordKey: 'ref', start: 1501, end: 8100 });
  await owner.search();
  const before = calls.steps.filter((step) => step.step === 'drawing').length;
  const [first, second] = owner.dialog.proposal[1].windows;
  const pending = owner.toggleWindow(first.key, false);
  assert.equal(owner.dialog.status, 'deriving');
  assert.equal(await owner.toggleWindow(second.key, false), '');
  await pending;
  assert.equal(calls.steps.filter((step) => step.step === 'drawing').length, before + 1);
  assert.deepEqual(owner.dialog.proposal[1].windows.map((window) => window.ticked), [false, true]);
  // The window after an unticked one moves up in the drawing: B2 is drawn first, numbered ①.
  assert.deepEqual(owner.dialog.proposal[1].windows.map((window) => window.number), ['', '①']);
});

test('typing in Start or End counts the genes once typing pauses', async () => {
  const { owner, calls } = harness();
  await owner.open();
  const [row] = owner.dialog.rows;
  await owner.editRow(row.id, { recordKey: 'ref' });
  const counted = calls.steps.filter((step) => step.step === 'reference').length;
  const typed = ['1', '15', '150', '1501'].map((start) => owner.editRow(row.id, { start }));
  assert.equal(owner.searchReason(), 'Counting the genes of the windows…');
  owner.editRow(row.id, { end: '8100' });
  await Promise.all(typed);
  assert.equal(calls.steps.filter((step) => step.step === 'reference').length, counted + 1);
  assert.deepEqual(calls.steps.at(-1).reference.windows, [{ start: 1501, end: 8100, reverseComplement: false }]);
  assert.equal(owner.dialog.genes, 3);
  // Closing drops a count that still waits.
  owner.editRow(row.id, { start: '1502' });
  owner.close();
  await new Promise((resolve) => setTimeout(resolve, 300));
  assert.equal(calls.steps.filter((step) => step.step === 'reference').length, counted + 1);
});
