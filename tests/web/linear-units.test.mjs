// The drawn units of a Linear Generate's comparison searches (OV-2738, D-01, D-37). Python plans the units and the
// unit pairs (tests/test_web_linear_units.py checks its plans); the Web asks for the plan with the records its request
// writes, then searches the units and writes their indexes. Each case here sends the Web's payload to the real Python
// helper and runs the Web over Python's reply.
import assert from 'node:assert/strict';
import test from 'node:test';

import { createLinearSeq } from '../../gbdraw/web/js/services/drawing/inputs.js';
import { buildLosatJobSpecs, linearComparisonEdgeKey } from '../../gbdraw/web/js/services/linear-comparisons.js';
import { losatSearchSequences } from '../../gbdraw/web/js/services/linear-sources.js';
import {
  linearDrawingUnitPlaces, linearWindowsPayload, normalizeWebGridColumnOrdering, requireLinearComparisonPlanSnapshot
} from '../../gbdraw/web/js/services/session-request.js';
import {
  linearUnitPlanPayload,
  planLinearUnitSearch,
  readLinearUnitPlan
} from '../../gbdraw/web/js/services/linear-units.js';
import {
  inheritSavedComparison,
  inheritUnitPlans,
  inheritedPipelineSearch,
  savedComparisonSource
} from '../../gbdraw/web/js/services/imported-comparison-intent.js';
import { pythonLinearUnitPlan } from './helpers/linear-unit-plan.mjs';

const pythonPlan = pythonLinearUnitPlan;
/** Python's `planLinearUnits`, as the diagram Worker runs it. */
const planLinearUnits = async (payload) => ({ result: pythonPlan(payload) });
/**
 * Generate's Inherit after its searches: Python plans the saved and the drawn units, and the comparison is carried by
 * them (`inheritUnitPlans`).
 * @param {Parameters<typeof inheritSavedComparison>[0]} input
 */
const inheritPlanned = async (input) => inheritSavedComparison({
  ...input, units: await inheritUnitPlans({ candidate: input.candidate, source: input.source, runHelper: planLinearUnits })
});
/**
 * The places the request writes for the drawn cards (`linearDrawingUnitPlaces`, the request's own owner), with the
 * Record Layout's rows when `rows` is given.
 * @param {{ uid: string, row: number }[] | null} [rows]
 */
const placesOf = (rows = null) => (cards) => linearDrawingUnitPlaces(
  { linearRecordLayoutEnabled: { value: rows !== null }, linearRecordRows: rows ?? [], linearSeqs: cards },
  () => ({ length: null, id: '', count: 1 }))(cards);

const window = (start, end, reverseComplement = false, column = null) => ({ start, end, reverseComplement, column });
const cases = [
  {
    id: 'two windows against a whole record',
    cards: [{ uid: 'a', windows: [window(1, 150), window(151, 300, true)] }, { uid: 'b' }],
    pairs: [['a', 'b', 'losat']],
    unitPairs: [[0, 2], [1, 2]]
  },
  {
    id: 'windows on both sides and a cropped record',
    cards: [
      { uid: 'a', windows: [window(1, 150, false, 1), window(151, 300, false, 2)] },
      { uid: 'b', windows: [window(10, 90), window(200, 280)] },
      { uid: 'c', region_start: 5, region_end: 250, region_reverse: true }
    ],
    pairs: [['a', 'b', 'losat'], ['b', 'c', 'losat']],
    unitPairs: [[0, 2], [0, 3], [1, 2], [1, 3], [2, 4], [3, 4]]
  },
  {
    id: 'a window across the origin',
    cards: [{ uid: 'a', windows: [window(16000, 700, true)] }, { uid: 'b' }],
    pairs: [['a', 'b', 'losat']],
    unitPairs: [[0, 1]]
  },
  {
    id: 'no windows',
    cards: [{ uid: 'a' }, { uid: 'b' }, { uid: 'c' }],
    pairs: [['a', 'b', 'losat'], ['c', 'b', 'upload']],
    unitPairs: [[0, 1], [2, 1]]
  }
].map((item) => ({
  ...item,
  cards: item.cards.map((card) => createLinearSeq(card)),
  pairs: item.pairs.map(([queryRecordKey, subjectRecordKey, source]) => ({ queryRecordKey, subjectRecordKey, source }))
}));

/** The comparison plan over the case's cards (one edge per record pair). */
const resolutionOf = (testCase, cards) => {
  const indexOf = new Map(cards.map((card, index) => [card.uid, index]));
  const edges = testCase.pairs.map(({ queryRecordKey, subjectRecordKey, source }, ordinal) => Object.freeze({
    id: `edge-${ordinal}`, queryUid: queryRecordKey, subjectUid: subjectRecordKey,
    queryIndex: indexOf.get(queryRecordKey), subjectIndex: indexOf.get(subjectRecordKey),
    edgeKey: linearComparisonEdgeKey(queryRecordKey, subjectRecordKey), ordinal, source, included: true,
    file: null, fileActive: false, losatFilename: '', losatFilenameActive: false
  }));
  return Object.freeze({
    mode: 'selected', defaultSource: 'losat', edges: Object.freeze(edges), errors: Object.freeze([]), error: '',
    valid: true, hasComparisonIntent: edges.length > 0, hasLosatIntent: edges.some((edge) => edge.source === 'losat'),
    hasUploadIntent: edges.some((edge) => edge.source === 'upload')
  });
};

for (const testCase of cases) {
  test(`the Web asks for and searches Python's unit plan: ${testCase.id}`, async () => {
    const { cards } = testCase;
    const resolution = resolutionOf(testCase, cards);
    const payload = linearUnitPlanPayload(cards, resolution, placesOf());
    assert.deepEqual(payload.pairs, testCase.pairs);
    const plan = pythonPlan(payload);
    const search = await planLinearUnitSearch(cards, resolution, placesOf(), async (asked) => {
      assert.deepEqual(asked, payload);
      return { result: plan };
    });
    const pairs = plan.pairs.map(({ queryUnit, subjectUnit }) => [queryUnit, subjectUnit]);
    assert.deepEqual(pairs, testCase.unitPairs);
    // Every request index is a unit index of Python's plan.
    assert.deepEqual(search.resolution.edges.map((edge) => [edge.queryIndex, edge.subjectIndex]), pairs);
    assert.deepEqual(search.resolution.edges.map((edge) => edge.ordinal), pairs.map((_, index) => index));
    assert.equal(new Set(search.resolution.edges.map((edge) => edge.edgeKey)).size, pairs.length);
    // Every unit pair of a LOSAT record pair is searched, in the plan's order.
    const losatPairs = plan.pairs
      .filter(({ pair }) => testCase.pairs[pair].source === 'losat')
      .map(({ queryUnit, subjectUnit }) => [queryUnit, subjectUnit]);
    const specs = buildLosatJobSpecs({
      resolution: search.resolution, recordCount: search.sequences.length,
      recordUids: search.sequences.map((seq) => seq.uid), program: 'blastn', collinearInferOrthogroups: false
    });
    assert.deepEqual(specs.map((spec) => [spec.queryIndex, spec.subjectIndex]), losatPairs);
    // Each unit searches the window of Python's plan, under its unit key when it is one of several windows; a
    // window across the origin is searched too, as Python crops it (OV-2739).
    plan.units.forEach((unit, index) => {
      const card = cards.find((item) => item.uid === unit.recordKey);
      assert.equal(search.cards[index], card);
      assert.equal(search.recordKeys[index], unit.recordKey);
      assert.equal(search.windowOf[index], unit.windowOf);
      assert.equal(search.windows[index], unit.window);
      assert.equal(search.gridPairs, plan.gridPairs);
      if (unit.windowOf === null) assert.equal(search.sequences[index], card);
      else assert.deepEqual(search.sequences[index], { ...card, uid: unit.unitKey });
    });
    const across = plan.units.filter(({ window: drawn }) => drawn && drawn.start > drawn.end);
    assert.deepEqual(across.map(({ window: drawn }) => drawn.searchRegion),
      testCase.id === 'a window across the origin' ? ['16000-700:origin'] : []);
    // The OFF cards follow the units, and no drawn card is searched twice.
    const off = createLinearSeq({ uid: 'off' });
    const searched = losatSearchSequences(search.sequences, [...cards, off], cards);
    assert.deepEqual(searched, { sequences: [...search.sequences, off], drawnCount: search.sequences.length });
  });
}

test('a reply that does not plan the records and pairs asked is refused', () => {
  const [testCase] = cases;
  const asked = linearUnitPlanPayload(testCase.cards, resolutionOf(testCase, testCase.cards), placesOf());
  const plan = pythonPlan(asked);
  const refused = (result) => assert.throws(() => readLinearUnitPlan(result, asked),
    (error) => error.code === 'RESULT_INVALID' && error.operation === 'planLinearUnits');
  refused(null);
  refused({ units: plan.units.slice(0, 2), pairs: plan.pairs });
  refused({ ...plan, pairs: [{ pair: 0, queryUnit: 2, subjectUnit: 0 }] });
  refused({ ...plan, pairs: [{ pair: 1, queryUnit: 0, subjectUnit: 2 }] });
  refused({ units: plan.units, pairs: plan.pairs });
  refused({ ...plan, gridPairs: [{ queryUnit: 2, subjectUnit: 0 }] });
  refused({ ...plan, gridPairs: [{ queryUnit: 0, subjectUnit: plan.units.length }] });
  assert.equal(readLinearUnitPlan(plan, asked), plan);
});

// OV-3248: which unit pairs a drawing displays is Python's (`grid_display_pairs`), whatever wrote its rows. The Web
// sends the rows its request writes (the Record Layout's, without columns) and reads the reply; it derives no pair.
test('the Web reads the unit pairs a grid displays from the unit plan of Python', async () => {
  const cards = ['a', 'b', 'c'].map((uid) => createLinearSeq({ uid }));
  const pairs = { pairs: [['a', 'c'], ['b', 'c']].map(([queryRecordKey, subjectRecordKey]) => (
    { queryRecordKey, subjectRecordKey, source: 'losat' })) };
  const resolution = resolutionOf(pairs, cards);
  const rows = [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }, { uid: 'c', row: 2 }];
  const searchOf = async (places) => planLinearUnitSearch(cards, resolution, places,
    async (asked) => ({ result: pythonPlan(asked) }));
  assert.deepEqual(linearUnitPlanPayload(cards, resolution, placesOf(rows)).records.map(({ gridRow }) => gridRow), [1, 1, 2]);
  assert.deepEqual((await searchOf(placesOf(rows))).gridPairs,
    [{ queryUnit: 0, subjectUnit: 2 }, { queryUnit: 1, subjectUnit: 2 }]);
  // Without the Record Layout every record takes its own row, in card order: no grid.
  assert.equal((await searchOf(placesOf())).gridPairs, null);
});

// OV-3003: Inherit of a saved LOSATP pipeline that names unit pairs (Find in other files keeps only the window pairs
// that share a reference gene) runs exactly those pairs, with the planner and runner the draft path uses; the draft
// plan of the same record pair compares every window pair.
test('Inherit plans the LOSAT jobs of exactly the saved unit pairs of a pipeline', async () => {
  const cards = [
    createLinearSeq({ uid: 'a', windows: [window(1, 150), window(151, 300)] }),
    createLinearSeq({ uid: 'b', windows: [window(10, 90), window(200, 280)] })
  ];
  const windowsOf = (card) => linearWindowsPayload(card, { arranged: false, ordinal: 1, record: { length: null, id: '' } });
  const settings = { proteinBlastpMaxHits: 3, proteinBlastpCandidateLimit: null };
  const saved = (pairs, records = cards, mode = 'pairwise') => savedComparisonSource({
    renderRequest: {
      schema: 10, mode: 'linear', diagramOptions: { bitscore: 40 },
      records: records.map((card) => ({ recordKey: card.uid, windows: windowsOf(card) })),
      comparisons: [{
        kind: 'generatedProteinComparison', mode, settings,
        pairs: pairs.map(([queryRecordIndex, subjectRecordIndex]) => ({ queryRecordIndex, subjectRecordIndex }))
      }]
    },
    resources: {}
  });
  /** @type {unknown[]} */
  const asked = [];
  const search = (source, layout = []) => inheritedPipelineSearch({
    source, cards, placesOf: placesOf(),
    runHelper: async (payload) => {
      asked.push(payload);
      return { result: pythonPlan(payload) };
    },
    requirePlan: requireLinearComparisonPlanSnapshot,
    plan: { sequences: cards, recordsOff: [], layout }
  });
  // Units 0 and 1 are a's windows, 2 and 3 b's: the pipeline keeps window 1 against window 1, window 2 against 2.
  const inherited = await search(saved([[0, 2], [1, 3]]));
  assert.ok(inherited);
  assert.deepEqual(inherited.pipeline.settings, settings);
  // Python plans the units of the saved records and of the drawn cards, which Inherit carries the saved indexes by,
  // then those of the searches; the Web derives no unit itself.
  const unitPlan = { records: cards.map((card) => ({ recordKey: card.uid, windows: windowsOf(card), gridRow: null })), pairs: [] };
  assert.deepEqual(asked, [unitPlan, unitPlan, unitPlan]);
  const { resolution, sequences } = inherited.unitSearch;
  assert.equal(resolution.valid, true, resolution.error);
  assert.equal(resolution.hasLosatIntent, true);
  assert.deepEqual(resolution.edges.map((edge) => [edge.queryIndex, edge.subjectIndex]), [[0, 2], [1, 3]]);
  const specs = buildLosatJobSpecs({
    resolution, recordCount: sequences.length, recordUids: sequences.map((seq) => seq.uid),
    program: 'blastp', blastpMode: 'pairwise', collinearInferOrthogroups: false
  });
  assert.deepEqual(specs.map((spec) => [spec.queryIndex, spec.subjectIndex]), [[0, 2], [1, 3]]);
  assert.deepEqual(specs.map((spec) => spec.edgeKey), resolution.edges.map((edge) => edge.edgeKey));
  // Generate writes what it searched, the tables and the marker, in place of the saved pipeline.
  const searched = [{ kind: 'precomputedProteinComparison', resourceId: 'searched', queryRecordIndex: 0, subjectRecordIndex: 2 }];
  const candidate = {
    renderRequest: {
      records: cards.map((card) => ({ recordKey: card.uid, windows: windowsOf(card) })), diagramOptions: {},
      comparisons: structuredClone(searched)
    },
    resources: { searched: { kind: 'canonical-tsv' } }
  };
  await inheritPlanned({ candidate, source: saved([[0, 2], [1, 3]]), searched: true });
  assert.deepEqual(candidate.renderRequest.comparisons, searched);
  assert.deepEqual(candidate.renderRequest.diagramOptions, { bitscore: 40 });
  // A pipeline that names no pairs (a CLI Session before 0.15.0) searches what Python's Pairwise searches without
  // pairs: every unit pair of the records of adjacent rows, as Python plans them for the adjacent record pair.
  asked.length = 0;
  const adjacent = await search(saved([]));
  assert.ok(adjacent);
  assert.deepEqual(asked.map((payload) => payload.pairs), [[], [], [{ queryRecordKey: 'a', subjectRecordKey: 'b', source: 'losat' }]]);
  const adjacentSpecs = buildLosatJobSpecs({
    resolution: adjacent.unitSearch.resolution, recordCount: adjacent.unitSearch.sequences.length,
    recordUids: adjacent.unitSearch.sequences.map((seq) => seq.uid),
    program: 'blastp', blastpMode: 'pairwise', collinearInferOrthogroups: false
  });
  assert.deepEqual(adjacentSpecs.map((spec) => [spec.queryIndex, spec.subjectIndex]), [[0, 2], [0, 3], [1, 2], [1, 3]]);
  // A saved comparison without a Pairwise pipeline is copied as saved; changed windows ask for the choice before a
  // search.
  assert.equal(await search(saved([], cards, 'none')), null);
  // Saved pairs whose records now share one row (moved in Record Layout) stop as the draft's would, with the pair
  // named, before any LOSAT job: Generate builds no job from a plan it refuses.
  await assert.rejects(search(saved([[0, 2], [1, 3]]), [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }]),
    (error) => error.code === 'COMPARISON_INPUT' && error.context?.reason === 'PAIR_TOPOLOGY'
      && error.context?.queryRecordIndex === 0 && error.context?.subjectRecordIndex === 1);
  asked.length = 0;
  await assert.rejects(search(saved([[0, 2]], [cards[0], createLinearSeq({ uid: 'b', windows: [window(10, 90)] })])),
    (error) => error.code === 'COMPARISON_CHOICE_REQUIRED' && error.context?.reason === 'INHERIT_UNITS_CHANGED');
  assert.deepEqual(asked.map((payload) => payload.pairs), [[], []], 'the two unit plans, and no search');
});

// OV-3365: Load orders the cards of a request Python wrote (`normalizeWebGridColumnOrdering`), and that order can differ
// from the order the request stores its records in, which the saved unit indexes count. Inherit carries each saved
// unit to the drawn unit of the same key and window in Python's plans, so the first Generate searches exactly the saved
// unit pairs. The cases: the drawing Find in other files writes for the Vibrio plasmids (a row of contig 80 in two
// windows at columns 1 and 4, contig 79 at column 2 and contig 86 at column 3), which Load orders 79, 80, 86; a row of
// A in two windows and B at column 2, which Load orders B, A; and A drawn twice in the same span (D-18, "+ Add window"
// on an uncropped card), whose order Load keeps.
test('Inherit carries each saved unit pair to the same windows when Load orders the cards differently', async () => {
  const settings = { proteinBlastpMaxHits: 3, proteinBlastpCandidateLimit: null };
  const cases = [
    {
      records: [
        { recordKey: 'reference', presentation: { gridRow: 1, gridColumn: 1 }, windows: [] },
        { recordKey: 'plasmid', presentation: { gridRow: 2, gridColumn: 1 }, windows: [window(32762, 32761, true)] },
        { recordKey: 'contig-80', presentation: { gridRow: 3, gridColumn: null },
          windows: [window(1, 13298, true, 1), window(13299, 30248, true, 4)] },
        { recordKey: 'contig-79', presentation: { gridRow: 3, gridColumn: 2 }, windows: [] },
        { recordKey: 'contig-86', presentation: { gridRow: 3, gridColumn: 3 }, windows: [window(1, 27527)] }
      ],
      loaded: ['reference', 'plasmid', 'contig-79', 'contig-80', 'contig-86'],
      // Units in stored order: reference 0, plasmid 1, contig 80's windows 2 and 3, contig 79 4, contig 86 5.
      pairs: [[0, 1], [1, 2], [1, 3], [1, 4], [1, 5]]
    },
    {
      records: [
        { recordKey: 'R', presentation: { gridRow: 1, gridColumn: 1 }, windows: [] },
        { recordKey: 'C', presentation: { gridRow: 3, gridColumn: 1 }, windows: [] },
        { recordKey: 'A', presentation: { gridRow: 2, gridColumn: null }, windows: [window(1, 100, false, 1), window(201, 300, false, 3)] },
        { recordKey: 'B', presentation: { gridRow: 2, gridColumn: 2 }, windows: [] }
      ],
      loaded: ['R', 'B', 'A', 'C'],
      // Units in stored order: R 0, C 1, A's windows 2 and 3, B 4.
      pairs: [[0, 2], [0, 3], [0, 4]]
    },
    {
      records: [
        { recordKey: 'A', presentation: { gridRow: 1, gridColumn: null }, windows: [window(1, 100, false, 1), window(1, 100, false, 2)] },
        { recordKey: 'B', presentation: { gridRow: 2, gridColumn: 1 }, windows: [] }
      ],
      loaded: ['A', 'B'],
      // Units in stored order: A's two windows of one span 0 and 1, B 2.
      pairs: [[0, 2], [1, 2]]
    }
  ];
  for (const { records, loaded, pairs } of cases) {
    const byKey = new Map(records.map((record) => [record.recordKey, record]));
    const cardOf = (record) => createLinearSeq({ uid: record.recordKey, windows: record.windows.map((entry) => ({ ...entry })) });
    const cards = normalizeWebGridColumnOrdering(records).records.map((record) => cardOf(byKey.get(record.recordKey)));
    assert.deepEqual(cards.map((card) => card.uid), loaded, 'the order Load gives the cards');
    const rows = records.map((record) => ({ uid: record.recordKey, row: record.presentation.gridRow }));
    // The identity of each unit, in a list of records' unit order: Python's plan of those records names it.
    const unitsOf = (keys) => pythonPlan({ records: keys.map((key) => ({ recordKey: key, windows: byKey.get(key).windows,
      gridRow: null })), pairs: [] }).units.map((unit) => unit.unitKey);
    const storedUnits = unitsOf(records.map((record) => record.recordKey));
    const cardUnits = unitsOf(cards.map((card) => card.uid));
    const source = savedComparisonSource({
      renderRequest: {
        schema: 10, mode: 'linear', diagramOptions: {},
        records: records.map(({ recordKey, windows }) => ({ recordKey, windows })),
        comparisons: [{ kind: 'generatedProteinComparison', mode: 'pairwise', settings,
          pairs: pairs.map(([queryRecordIndex, subjectRecordIndex]) => ({ queryRecordIndex, subjectRecordIndex })) }]
      },
      resources: {}
    });
    const inherited = await inheritedPipelineSearch({
      source, cards, placesOf: placesOf(rows), runHelper: planLinearUnits,
      requirePlan: requireLinearComparisonPlanSnapshot, plan: { sequences: cards, recordsOff: [], layout: rows }
    });
    assert.ok(inherited);
    const searched = inherited.unitSearch.resolution.edges.map((edge) => [cardUnits[edge.queryIndex], cardUnits[edge.subjectIndex]]);
    assert.deepEqual(searched.map((pair) => pair.join(' ')).sort(),
      pairs.map(([query, subject]) => `${storedUnits[query]} ${storedUnits[subject]}`).sort());
    // Generate then writes the searched tables in place of the pipeline: Inherit accepts the cards' units.
    const candidate = { renderRequest: { records: cards.map((card) => ({ recordKey: card.uid,
      windows: linearWindowsPayload(card, { arranged: true, ordinal: 1, record: { length: null, id: '' } }) })),
      diagramOptions: {}, comparisons: [] }, resources: {} };
    await inheritPlanned({ candidate, source, searched: true });
    // A saved comparison copied as saved names the same units by their place in this request.
    const copied = { renderRequest: { ...structuredClone(candidate.renderRequest), comparisons: [] }, resources: {} };
    const tables = savedComparisonSource({ renderRequest: { schema: 10, mode: 'linear', records: source.records, diagramOptions: {},
      comparisons: pairs.map(([queryRecordIndex, subjectRecordIndex]) => ({ kind: 'precomputedProteinComparison',
        resourceId: 'protein', encoding: 'canonicalTsv', queryRecordIndex, subjectRecordIndex })) },
    resources: { protein: { kind: 'canonical-tsv' } } });
    await inheritPlanned({ candidate: copied, source: tables });
    assert.deepEqual(copied.renderRequest.comparisons.map((table) => `${cardUnits[table.queryRecordIndex]} ${cardUnits[table.subjectRecordIndex]}`),
      pairs.map(([query, subject]) => `${storedUnits[query]} ${storedUnits[subject]}`));
  }
});

// A saved window that is no request window (a hand-edited or damaged Session) is not sent to Python: Inherit asks for
// the choice in the Comparison panel, as for any comparison it cannot carry.
test('Inherit asks for the choice when a saved window is malformed', async () => {
  const records = [
    { recordKey: 'A', windows: [window(1, 100, false, 1), { start: 'x', end: 200, reverseComplement: false, column: 2 }] },
    { recordKey: 'B', windows: [] }
  ];
  const source = savedComparisonSource({
    renderRequest: { schema: 10, mode: 'linear', diagramOptions: {}, records,
      comparisons: [{ kind: 'precomputedProteinComparison', resourceId: 'protein', encoding: 'canonicalTsv',
        queryRecordIndex: 0, subjectRecordIndex: 2 }] },
    resources: { protein: { kind: 'canonical-tsv' } }
  });
  const candidate = { renderRequest: { records: [{ recordKey: 'A', windows: [window(1, 100, false, 1), window(101, 200, false, 2)] },
    { recordKey: 'B', windows: [] }], diagramOptions: {}, comparisons: [] }, resources: {} };
  await assert.rejects(inheritPlanned({ candidate, source }),
    (error) => error.code === 'COMPARISON_CHOICE_REQUIRED' && error.context?.reason === 'INHERIT_UNITS_CHANGED');
});
