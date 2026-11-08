// @ts-check
// The drawn units of a Linear Generate's comparison searches, planned by Python (`planLinearUnits`,
// gbdraw/web_support/linear_units.py; D-01, D-37). A record drawn in several windows draws one unit per window, and a
// comparison of two records searches every window of one against every window of the other. The page sends the
// records its request will draw and the record pairs of its comparison plan, reads Python's plan, and runs its
// searches over the units: every index a search writes into the request is a unit index, as Python reads it. The Web
// derives no unit and no pair itself.
/** @import { LinearSeq } from './drawing/inputs.js' */
/** @import { LinearComparisonResolution, LinearComparisonResolvedEdge } from './linear-comparisons.js' */
/** @import { CanonicalRecordWindow } from '../contracts/request.generated.js' */
import { linearComparisonEdgeKey } from './linear-comparisons.js';
import { linearSeqRecordKey } from './record-draw-selection.js';
import { linearWindowsPayload } from './session-request.js';
import { diagnosticError } from '../utils/error-normalization.js';

const FAILURE = /** @type {const} */ ({ stage: 'helper', operation: 'planLinearUnits' });

/**
 * A record the request will draw (`LinearUnitRecordInput`): its key and its windows, as the request writes them.
 * @typedef {{ recordKey: string, windows: CanonicalRecordWindow[] }} LinearUnitRecordInput
 */
/**
 * A record pair of the comparison plan (`LinearUnitPairInput`), by record keys.
 * @typedef {{ queryRecordKey: string, subjectRecordKey: string, source: 'losat' | 'upload' }} LinearUnitPairInput
 */
/** The payload of `planLinearUnits`. @typedef {{ records: LinearUnitRecordInput[], pairs: LinearUnitPairInput[] }} LinearUnitPlanPayload */
/**
 * The span of its record a unit draws (`LinearUnitWindow`), in source coordinates, 1-based and inclusive (a start
 * after the end runs across the origin of a circular record), and `searchRegion`, the region string of the span a
 * comparison search reads (the crop before its reverse complement), which only the Worker's extraction reads.
 * @typedef {{ start: number, end: number, reverseComplement: boolean, searchRegion: string }} LinearUnitWindow
 */
/**
 * One drawn unit (`LinearDrawnUnit`): its key, its record's key, the record key when it is one of several windows of
 * it (else null), and its window (null for a record drawn whole).
 * @typedef {{ unitKey: string, recordKey: string, windowOf: string | null, window: LinearUnitWindow | null }} LinearDrawnUnit
 */
/** One unit pair of record pair `pair` (`LinearUnitPair`). @typedef {{ pair: number, queryUnit: number, subjectUnit: number }} LinearUnitPair */
/** The `planLinearUnits` reply (`LinearUnitPlan`). @typedef {{ units: LinearDrawnUnit[], pairs: LinearUnitPair[] }} LinearUnitPlan */

/**
 * The searches of one Generate over the drawn units: each unit's card, under the unit's key when it is one of several
 * windows (`sequences`), the card itself (`cards`), the unit's record key, `windowOf` and window (the span its search
 * reads, null for a record drawn whole), and the comparison plan whose edges join units (`resolution`).
 * @typedef {object} LinearUnitSearch
 * @property {LinearSeq[]} sequences
 * @property {LinearSeq[]} cards
 * @property {string[]} recordKeys
 * @property {(string | null)[]} windowOf
 * @property {(LinearUnitWindow | null)[]} windows
 * @property {LinearComparisonResolution} resolution
 */

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/** @param {unknown} value @returns {value is number} */
const isIndex = (value) => Number.isSafeInteger(value) && Number(value) >= 0;

/** @param {unknown} value @returns {value is LinearUnitWindow | null} */
const isWindow = (value) => value === null || (isRecord(value) && isIndex(value.start) && isIndex(value.end)
  && typeof value.reverseComplement === 'boolean' && typeof value.searchRegion === 'string' && value.searchRegion !== '');

/** @param {unknown} value @returns {value is LinearDrawnUnit} */
const isUnit = (value) => isRecord(value) && typeof value.unitKey === 'string' && typeof value.recordKey === 'string'
  && (value.windowOf === null || typeof value.windowOf === 'string') && isWindow(value.window);

/**
 * The helper's result, checked against `LinearUnitPlan`: every unit names a record asked, every record asked has a
 * unit, and every pair joins units of a pair asked.
 * @param {unknown} result
 * @param {LinearUnitPlanPayload} asked
 * @returns {LinearUnitPlan}
 */
export const readLinearUnitPlan = (result, asked) => {
  const keys = new Set(asked.records.map(({ recordKey }) => recordKey));
  const units = isRecord(result) && Array.isArray(result.units) ? result.units : null;
  const pairs = isRecord(result) && Array.isArray(result.pairs) ? result.pairs : null;
  const valid = units !== null && pairs !== null && units.every(isUnit)
    && units.every((unit) => keys.has(unit.recordKey))
    && new Set(units.map((unit) => unit.recordKey)).size === keys.size
    && pairs.every((pair) => isRecord(pair) && isIndex(pair.pair) && pair.pair < asked.pairs.length
      && isIndex(pair.queryUnit) && pair.queryUnit < units.length
      && isIndex(pair.subjectUnit) && pair.subjectUnit < units.length
      && units[pair.queryUnit].recordKey === asked.pairs[pair.pair].queryRecordKey
      && units[pair.subjectUnit].recordKey === asked.pairs[pair.pair].subjectRecordKey);
  if (!valid) throw diagnosticError('RESULT_INVALID', {}, FAILURE);
  return /** @type {LinearUnitPlan} */ (/** @type {unknown} */ (result));
};

/**
 * The card of each unit, found by its record key; a unit that is one of several windows of its record searches
 * under its unit key. The span a unit's search reads is its window (`LinearUnitSearch.windows`), not the card's crop.
 * @param {readonly LinearSeq[]} cards
 * @param {LinearUnitPlan} plan
 */
const unitSequences = (cards, plan) => {
  const cardByKey = new Map(cards.map((card, index) => [linearSeqRecordKey(card, index), card]));
  return plan.units.map((unit) => {
    const card = /** @type {LinearSeq} */ (cardByKey.get(unit.recordKey));
    return { card, sequence: unit.windowOf === null ? card : { ...card, uid: unit.unitKey } };
  });
};

/**
 * The comparison plan over units: each edge of `resolution` becomes the edges of its unit pairs, numbered in order;
 * an edge of a record drawn whole keeps its key, so a plan without windows is `resolution`'s.
 * @param {LinearComparisonResolution} resolution
 * @param {LinearUnitPlan} plan
 * @returns {LinearComparisonResolution}
 */
const unitResolution = (resolution, plan) => {
  const edges = plan.pairs.map(({ pair, queryUnit, subjectUnit }, ordinal) => {
    const edge = resolution.edges[pair];
    const queryUid = plan.units[queryUnit].unitKey;
    const subjectUid = plan.units[subjectUnit].unitKey;
    /** @type {LinearComparisonResolvedEdge} */
    const unitEdge = {
      ...edge, queryUid, subjectUid, queryIndex: queryUnit, subjectIndex: subjectUnit,
      edgeKey: linearComparisonEdgeKey(queryUid, subjectUid), ordinal
    };
    return Object.freeze(unitEdge);
  });
  return Object.freeze({ ...resolution, edges: Object.freeze(edges) });
};

/**
 * The searches of one Generate over the units of `plan` (`LinearUnitSearch`).
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {LinearComparisonResolution} resolution The comparison plan over the cards.
 * @param {LinearUnitPlan} plan
 * @returns {LinearUnitSearch}
 */
export const linearUnitSearch = (cards, resolution, plan) => {
  const units = unitSequences(cards, plan);
  return {
    sequences: units.map(({ sequence }) => sequence),
    cards: units.map(({ card }) => card),
    recordKeys: plan.units.map(({ recordKey }) => recordKey),
    windowOf: plan.units.map(({ windowOf }) => windowOf),
    windows: plan.units.map(({ window }) => window),
    resolution: unitResolution(resolution, plan)
  };
};

/**
 * The unit searches of a Linear request's comparison plan: Python's plan of the drawn units when `resolution`
 * compares records, else null (the request writes `resolution` itself). Generate and the draft request (Save without
 * a Result, New drawing from a region) read the plan here, so both write the same unit indexes (OV-2738).
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {LinearComparisonResolution} resolution
 * @param {boolean} arranged Whether the Record Layout arranges the rows (`linearWindowsPayload`).
 * @param {(payload: LinearUnitPlanPayload) => Promise<{ result?: unknown }>} runHelper Runs the Worker helper.
 * @returns {Promise<LinearUnitSearch | null>}
 */
export const planRequestUnitSearch = async (cards, resolution, arranged, runHelper) => {
  if (resolution.hasComparisonIntent !== true) return null;
  const payload = requestUnitPlanPayload(cards, resolution, arranged);
  const reply = await runHelper(payload);
  return linearUnitSearch(cards, resolution, readLinearUnitPlan(reply?.result, payload));
};

/**
 * What `planLinearUnits` is asked: the drawn cards' records, keyed and windowed as the request writes them, and the
 * record pairs of the plan's edges.
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {LinearComparisonResolution} resolution
 * @param {boolean} arranged Whether the Record Layout arranges the rows (`linearWindowsPayload`).
 * @returns {LinearUnitPlanPayload}
 */
export const requestUnitPlanPayload = (cards, resolution, arranged) => ({
  records: cards.map((card, index) => ({
    recordKey: linearSeqRecordKey(card, index), windows: linearWindowsPayload(card, { arranged })
  })),
  pairs: resolution.edges.map((edge) => ({ queryRecordKey: edge.queryUid, subjectRecordKey: edge.subjectUid, source: edge.source }))
});

