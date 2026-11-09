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
import { diagnosticError } from '../utils/error-normalization.js';

const FAILURE = /** @type {const} */ ({ stage: 'helper', operation: 'planLinearUnits' });

/**
 * A record the request will draw (`LinearUnitRecordInput`): its key and its windows, as the request writes them.
 * @typedef {{ recordKey: string, windows: CanonicalRecordWindow[] }} LinearUnitRecordInput
 */
/**
 * The windows the request writes for a card (`linearWindowsPayload` of `services/session-request.js`, which the
 * caller passes with its layout, so this module reads the request writer's windows without importing it).
 * @typedef {(card: LinearSeq) => CanonicalRecordWindow[]} LinearUnitWindowsOf
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
 * Whether drawn cards or request records draw several units: only one whose window list holds several windows does
 * (D-01). Without one, Python's unit plan is the records themselves and the plan's own edges, so a request needs no
 * plan to write them, and a unit index is a record index.
 * @param {readonly { windows: readonly unknown[] }[]} records Drawn cards, or request records.
 */
export const drawsSeveralUnits = (records) => records.some((record) => record.windows.length > 1);

/**
 * The record each drawn unit draws, in unit order: the record's key (`recordKeys`), and that key again when the unit
 * is one of several windows of the record, else null (`windowOf`). Python's unit plan gives it (`LinearUnitSearch`),
 * and so does a catalog item (`catalogUnitRecords`).
 * @typedef {Pick<LinearUnitSearch, 'recordKeys' | 'windowOf'>} LinearUnitRecords
 */
/**
 * The record each drawn unit of a catalog item draws (`LinearUnitRecords`): a unit `windowOf` names is one of
 * several windows of the record it maps to; any other unit is its own record (Python's `catalog_record_units`,
 * gbdraw/render/interactive_svg.py). A `windowOf` that is no object names none; admission validates it
 * (feature-catalog.js).
 * @param {readonly string[]} recordKeys The item's `recordKeys`: its units, in drawing order.
 * @param {unknown} windowOf The item's `windowOf`.
 * @returns {LinearUnitRecords}
 */
export const catalogUnitRecords = (recordKeys, windowOf) => {
  const windows = isRecord(windowOf) ? windowOf : {};
  const recordOf = recordKeys.map((unitKey) => (Object.hasOwn(windows, unitKey) ? String(windows[unitKey]) : null));
  return { recordKeys: recordKeys.map((unitKey, unit) => recordOf[unit] ?? unitKey), windowOf: recordOf };
};

/**
 * Where a drawn unit lies (`LinearUnitPlace`): the key of the record it draws (its card's), and, when it is one of
 * several windows of that record, the window's number among them (from 1), else null.
 * @typedef {{ recordKey: string, window: number | null }} LinearUnitPlace
 */
/**
 * The place of each unit a drawn Result's catalog names, by unit key (`LinearUnitPlace`): its record
 * (`catalogUnitRecords`) and its window number, counted in drawing order, which is the order of the record's window
 * list (Python's `record_drawn_windows`). Empty without a catalog item. Every item names the same units
 * (`catalogDrawnUnits`), so the first is read.
 * @param {{ readonly items: readonly { recordKeys: readonly string[], windowOf?: unknown }[] } | null} catalog
 * @returns {Map<string, LinearUnitPlace>}
 */
export const catalogUnitPlaces = (catalog) => {
  const [item] = catalog === null ? [] : catalog.items;
  if (!item) return new Map();
  const units = catalogUnitRecords(item.recordKeys, item.windowOf);
  /** @type {Map<string, number>} */
  const counted = new Map();
  return new Map(item.recordKeys.map((unitKey, unit) => {
    const of = units.windowOf[unit];
    const window = of === null ? null : (counted.get(of) ?? 0) + 1;
    if (of !== null && window !== null) counted.set(of, window);
    return [unitKey, { recordKey: units.recordKeys[unit], window }];
  }));
};

/**
 * Whether a LOSAT cache-info row is a window pair of the drawn Result (OV-2742, OV-2918): one of its two units is a
 * window of a record drawn in several (`places`, the Result's units), and each unit is such a window or a record drawn
 * whole, of a record `drawn` names. Its indexes count that Result's units, so a reorder or a Load keeps it as saved.
 * @param {{ queryUid?: unknown, subjectUid?: unknown }} entry
 * @param {ReadonlyMap<string, LinearUnitPlace>} places
 * @param {(recordKey: string) => boolean} drawn
 */
export const isDrawnWindowPairRow = (entry, places, drawn) => {
  /** @param {unknown} uid */
  const drawnWindow = (uid) => {
    const place = places.get(String(uid ?? ''));
    return place !== undefined && place.window !== null && drawn(place.recordKey);
  };
  /** @param {unknown} uid */
  const drawnUnit = (uid) => drawnWindow(uid) || drawn(String(uid ?? ''));
  return (drawnWindow(entry.queryUid) || drawnWindow(entry.subjectUid))
    && drawnUnit(entry.queryUid) && drawnUnit(entry.subjectUid);
};

/**
 * One window pair of a card pair (`LinearWindowPair`): the edge key of its LOSAT cache-info row (the pair of drawn
 * units Generate searched), the row's file name, and each side's window number (null for a card drawn whole).
 * @typedef {{ edgeKey: string, filename: string, queryWindow: number | null, subjectWindow: number | null }} LinearWindowPair
 */
/**
 * The window pairs of each card pair, by the card pair's edge key (its two records' keys): the cache-info rows whose
 * two units `places` knows and at least one of which is one of several windows, in row order. A row whose two units
 * are cards drawn whole is the card pair's own row, which the card reads by its edge key.
 * @param {readonly { edgeKey?: string, queryUid?: string, subjectUid?: string, filename: string }[]} entries
 * @param {ReadonlyMap<string, LinearUnitPlace>} places
 * @returns {Map<string, LinearWindowPair[]>}
 */
export const linearWindowPairsByCardEdge = (entries, places) => {
  /** @type {Map<string, LinearWindowPair[]>} */
  const byCard = new Map();
  entries.forEach((entry) => {
    const query = places.get(String(entry.queryUid ?? ''));
    const subject = places.get(String(entry.subjectUid ?? ''));
    if (!entry.edgeKey || !query || !subject || (query.window === null && subject.window === null)) return;
    const cardEdge = linearComparisonEdgeKey(query.recordKey, subject.recordKey);
    byCard.set(cardEdge, [...(byCard.get(cardEdge) ?? []), {
      edgeKey: entry.edgeKey, filename: entry.filename, queryWindow: query.window, subjectWindow: subject.window
    }]);
  });
  return byCard;
};

/**
 * The request record a drawn unit draws: its index among the request's records, and whether the unit is one of
 * several windows of it.
 * @typedef {{ recordIndex: number, oneOfSeveral: boolean }} LinearUnitRecord
 */

/**
 * The request record each drawn unit draws, found by its record key (null for a unit whose record `records` does not
 * hold). A request's comparison indexes and a catalog's sources and match endpoints count drawn units; the CLI's `#n`
 * and `records[]` count records: this is the one reader from the first to the second. Without `units`, the units are
 * the records when none is drawn in several windows (`drawsSeveralUnits`), and unknown (null) otherwise.
 * @param {readonly { recordKey: string, windows: readonly unknown[] }[]} records The request's records.
 * @param {LinearUnitRecords | null} units The request's drawn units, as Python planned or drew them.
 * @returns {(LinearUnitRecord | null)[] | null}
 */
export const unitRecordsOf = (records, units) => {
  if (!units) return drawsSeveralUnits(records) ? null : records.map((_, recordIndex) => ({ recordIndex, oneOfSeveral: false }));
  const indexByKey = new Map(records.map(({ recordKey }, index) => [recordKey, index]));
  return units.recordKeys.map((recordKey, unit) => {
    const recordIndex = indexByKey.get(recordKey);
    return recordIndex === undefined ? null : { recordIndex, oneOfSeveral: units.windowOf[unit] !== null };
  });
};

/**
 * The drawn cards' records, keyed and windowed as the request writes them: what Python plans the units of.
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {LinearUnitWindowsOf} windowsOf The windows the request writes for a card.
 * @returns {LinearUnitRecordInput[]}
 */
export const unitRecordInputs = (cards, windowsOf) => cards.map((card, index) => (
  { recordKey: linearSeqRecordKey(card, index), windows: windowsOf(card) }));

/**
 * What `planLinearUnits` is asked: the drawn cards' records (`unitRecordInputs`) and the record pairs of the plan's
 * edges.
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {LinearComparisonResolution} resolution
 * @param {LinearUnitWindowsOf} windowsOf The windows the request writes for a card.
 * @returns {LinearUnitPlanPayload}
 */
export const linearUnitPlanPayload = (cards, resolution, windowsOf) => ({
  records: unitRecordInputs(cards, windowsOf),
  pairs: resolution.edges.map((edge) => ({ queryRecordKey: edge.queryUid, subjectRecordKey: edge.subjectUid, source: edge.source }))
});

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
 * Plans the drawn units of a Linear request's comparison plan with Python: Generate, the draft request (Save without
 * a Result, New drawing from a region) and the LOSAT job estimate read it here, so all write or count unit indexes.
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {LinearComparisonResolution} resolution
 * @param {LinearUnitWindowsOf} windowsOf
 * @param {(payload: LinearUnitPlanPayload) => Promise<{ result?: unknown }>} runHelper Runs the Worker helper.
 * @returns {Promise<LinearUnitSearch>}
 */
export const planLinearUnitSearch = async (cards, resolution, windowsOf, runHelper) => {
  const payload = linearUnitPlanPayload(cards, resolution, windowsOf);
  const reply = await runHelper(payload);
  return linearUnitSearch(cards, resolution, readLinearUnitPlan(reply?.result, payload));
};

/**
 * The searches of saved unit pairs (a saved LOSATP pipeline that Inherit runs, OV-3003): Python plans the units of the
 * drawn cards, asked for no record pair, and the plan joins exactly `pairs`, each under the record pair its two units
 * draw. `resolve` gives the comparison plan of those record pairs, one edge per pair in their order.
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {readonly Omit<LinearUnitPair, 'pair'>[]} pairs Indexes of the units Python plans for `cards`.
 * @param {LinearUnitWindowsOf} windowsOf
 * @param {(payload: LinearUnitPlanPayload) => Promise<{ result?: unknown }>} runHelper Runs the Worker helper.
 * @param {(pairs: LinearUnitPairInput[]) => LinearComparisonResolution} resolve
 * @returns {Promise<LinearUnitSearch>}
 */
export const savedPairsUnitSearch = async (cards, pairs, windowsOf, runHelper, resolve) => {
  /** @type {LinearUnitPlanPayload} */
  const payload = { records: unitRecordInputs(cards, windowsOf), pairs: [] };
  const { units } = readLinearUnitPlan((await runHelper(payload))?.result, payload);
  if (pairs.some(({ queryUnit, subjectUnit }) => !units[queryUnit] || !units[subjectUnit])) {
    throw diagnosticError('RESULT_INVALID', {}, FAILURE);
  }
  /** @type {LinearUnitPairInput[]} */
  const recordPairs = [];
  /** @type {Map<string, number>} */
  const pairOf = new Map();
  const unitPairs = pairs.map(({ queryUnit, subjectUnit }) => {
    const [queryRecordKey, subjectRecordKey] = [units[queryUnit].recordKey, units[subjectUnit].recordKey];
    const key = linearComparisonEdgeKey(queryRecordKey, subjectRecordKey);
    if (!pairOf.has(key)) pairOf.set(key, recordPairs.push({ queryRecordKey, subjectRecordKey, source: 'losat' }) - 1);
    return { pair: /** @type {number} */ (pairOf.get(key)), queryUnit, subjectUnit };
  });
  return linearUnitSearch(cards, resolve(recordPairs), { units, pairs: unitPairs });
};

/**
 * The units a Generate searches: Python's plan when the drawing has a comparison; without one nothing is searched
 * and the cards are the units (no window, none a window of another).
 * @param {readonly LinearSeq[]} cards The drawn cards, in request order.
 * @param {LinearComparisonResolution} resolution
 * @param {LinearUnitWindowsOf} windowsOf
 * @param {(payload: LinearUnitPlanPayload) => Promise<{ result?: unknown }>} runHelper Runs the Worker helper.
 * @returns {Promise<LinearUnitSearch>}
 */
export const generateUnitSearch = async (cards, resolution, windowsOf, runHelper) => (resolution.hasComparisonIntent === true
  ? planLinearUnitSearch(cards, resolution, windowsOf, runHelper)
  : {
      sequences: [...cards], cards: [...cards], recordKeys: cards.map((card, index) => linearSeqRecordKey(card, index)),
      windowOf: cards.map(() => null), windows: cards.map(() => null), resolution
    });
