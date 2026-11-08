// The drawn units of a Linear Generate's comparison searches (OV-2738, D-01, D-37). Python plans the units and the
// unit pairs (tests/test_web_linear_units.py checks its plans); the Web asks for the plan with the records its request
// writes, then searches the units and writes their indexes. Each case here sends the Web's payload to the real Python
// helper and runs the Web over Python's reply.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import test from 'node:test';

import { createLinearSeq } from '../../gbdraw/web/js/services/drawing/inputs.js';
import { buildLosatJobSpecs, linearComparisonEdgeKey } from '../../gbdraw/web/js/services/linear-comparisons.js';
import { losatSearchSequences } from '../../gbdraw/web/js/services/linear-sources.js';
import { linearWindowsPayload } from '../../gbdraw/web/js/services/session-request.js';
import {
  linearUnitPlanPayload,
  planLinearUnitSearch,
  readLinearUnitPlan
} from '../../gbdraw/web/js/services/linear-units.js';

/** Python's reply to `planLinearUnits` for a payload. */
const pythonPlan = (payload) => JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json, sys
from gbdraw.web_support.linear_units import plan_linear_units_json
payload = json.load(sys.stdin)
print(plan_linear_units_json(json.dumps(payload["records"]), json.dumps(payload["pairs"])))
`], { input: JSON.stringify(payload), encoding: 'utf8' }));

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
    const windowsOf = (card) => linearWindowsPayload(card, { arranged: true });
    const payload = linearUnitPlanPayload(cards, resolution, windowsOf);
    assert.deepEqual(payload.pairs, testCase.pairs);
    const plan = pythonPlan(payload);
    const search = await planLinearUnitSearch(cards, resolution, windowsOf, async (asked) => {
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
  const asked = linearUnitPlanPayload(testCase.cards, resolutionOf(testCase, testCase.cards),
    (card) => linearWindowsPayload(card, { arranged: true }));
  const plan = pythonPlan(asked);
  const refused = (result) => assert.throws(() => readLinearUnitPlan(result, asked),
    (error) => error.code === 'RESULT_INVALID' && error.operation === 'planLinearUnits');
  refused(null);
  refused({ units: plan.units.slice(0, 2), pairs: plan.pairs });
  refused({ ...plan, pairs: [{ pair: 0, queryUnit: 2, subjectUnit: 0 }] });
  refused({ ...plan, pairs: [{ pair: 1, queryUnit: 0, subjectUnit: 2 }] });
  assert.equal(readLinearUnitPlan(plan, asked), plan);
});
