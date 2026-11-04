// @ts-check
import { diagnosticError } from '../utils/error-normalization.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ProjectedDraftRequest, RegionDerivation } from '../services/config.js' */
/** @import { DerivedRegionDrawing, RefusedRegionDrawing, RegionDrawingHelperPayload, RegionDrawingRegion, RegionDrawingReply } from '../services/diagram-worker-protocol.js' */
/** @import { DrawingListOutcome } from './drawing-list.js' */

/**
 * @typedef {object} RegionDrawingPorts
 * @property {() => DrawingId} shownDrawingId
 * @property {() => Promise<ProjectedDraftRequest>} projectShownRequest The shown drawing's draft request.
 * @property {(payload: RegionDrawingHelperPayload) => Promise<{ result?: unknown }>} deriveRegion
 *   Runs the Worker helper `deriveRegionDrawing`; the reply is read here.
 * @property {(sourceId: DrawingId, derivation: RegionDerivation) => Promise<DrawingListOutcome>} newDrawingFromRegion
 */

/**
 * What the creation screen shows: the refusal, or the derived drawing (what
 * returned to Auto, what is not carried) with the Create that adds it.
 * @typedef {RefusedRegionDrawing | (DerivedRegionDrawing & { create: () => Promise<DrawingListOutcome> })} RegionDrawingOffer
 */

const REPLY_INVALID = () => diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: 'deriveRegionDrawing' });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/**
 * The Worker's reply, checked against `RegionDrawingReply`.
 * @param {unknown} result
 * @returns {RegionDrawingReply}
 */
const readRegionDrawingReply = (result) => {
  if (!isRecord(result)) throw REPLY_INVALID();
  if (result.status === 'refused') {
    if (result.reason !== 'ORIGIN_CROSSING' || typeof result.recordKey !== 'string'
      || !Number.isSafeInteger(result.start) || !Number.isSafeInteger(result.end)) throw REPLY_INVALID();
    return /** @type {RefusedRegionDrawing} */ (result);
  }
  const { resourceSources } = result;
  if (result.status !== 'derived' || (result.mode !== 'circular' && result.mode !== 'linear')
    || !isRecord(result.renderRequest) || result.renderRequest.mode !== result.mode
    || !isRecord(resourceSources) || !Object.values(resourceSources).every((id) => typeof id === 'string')
    || !Array.isArray(result.returnedToAuto) || !Array.isArray(result.notCarried)) throw REPLY_INVALID();
  return /** @type {DerivedRegionDrawing} */ (/** @type {unknown} */ (result));
};

/**
 * New drawing from a region (PD-OI-090 〔作成〕(4)): Python derives the
 * drawing from the shown drawing's draft request; Create adds it with the
 * Web files of the request that was sent.
 * @param {RegionDrawingPorts} ports
 */
export const createRegionDrawing = ({ shownDrawingId, projectShownRequest, deriveRegion, newDrawingFromRegion }) => {
  /**
   * @param {readonly RegionDrawingRegion[]} regions Source coordinates, 1-based and inclusive.
   * @param {boolean} adaptSizes Whether hand-set sizes whose Auto differs return to Auto.
   * @returns {Promise<RegionDrawingOffer>}
   */
  const deriveFromShown = async (regions, adaptSizes) => {
    const sourceId = shownDrawingId();
    const sent = await projectShownRequest();
    const { renderRequest, resources } = sent.envelope;
    const reply = readRegionDrawingReply((await deriveRegion({
      canonicalRequest: renderRequest,
      resources: /** @type {RegionDrawingHelperPayload['resources']} */ (resources),
      regions,
      adaptSizes
    })).result);
    if (reply.status === 'refused') return reply;
    return { ...reply, create: () => newDrawingFromRegion(sourceId, { derived: reply, sent }) };
  };
  return { deriveFromShown };
};
