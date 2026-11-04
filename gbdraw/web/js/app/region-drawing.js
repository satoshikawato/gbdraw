// @ts-check
import { diagnosticError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { msg, t, tm } from '../utils/i18n.js';
import { reactive } from '../utils/reactivity.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ProjectedDraftRequest, RegionDerivation } from '../services/config.js' */
/** @import { DerivedRegionDrawing, RefusedRegionDrawing, RegionNotCarriedItem, RegionDrawingHelperPayload, RegionDrawingRegion, RegionDrawingReply } from '../services/diagram-worker-protocol.js' */
/** @import { DrawingListOutcome } from './drawing-list.js' */

/**
 * @typedef {object} RegionDrawingPorts
 * @property {() => DrawingId} shownDrawingId
 * @property {() => Promise<ProjectedDraftRequest>} projectShownRequest The shown drawing's draft request.
 * @property {(payload: RegionDrawingHelperPayload) => Promise<{ result?: unknown }>} deriveRegion
 *   Runs the Worker helper `deriveRegionDrawing`; the reply is read here.
 * @property {(sourceId: DrawingId, derivation: RegionDerivation) => Promise<DrawingListOutcome>} newDrawingFromRegion
 * @property {() => readonly RegionSourceRecord[]} shownRecords The shown drawing's records, in drawing order.
 * @property {(setting: string) => string} settingLabel A setting's label in the settings catalog.
 */

/**
 * A record of the source drawing as the dialog offers it. `features` are the
 * shown Result's admitted features of the record, or `null` without a Result.
 * @typedef {object} RegionSourceRecord
 * @property {string} recordKey The record key of the source drawing's request (the identity).
 * @property {string} caption Display only.
 * @property {number} length Source length in bases.
 * @property {readonly { featureId: string, label: string }[] | null} features
 */

/**
 * One record row of the dialog: display values and the user's input only.
 * @typedef {object} RegionDialogRow
 * @property {string} recordKey
 * @property {string} caption
 * @property {number} length
 * @property {boolean} checked
 * @property {'whole' | 'coordinates' | 'feature'} kind
 * @property {string} start
 * @property {string} end
 * @property {string} featureId
 * @property {readonly { featureId: string, label: string }[] | null} features
 */

/**
 * The dialog's view. `status`: 'idle' (no region yet), 'invalid' (an input
 * is not a region), 'deriving', 'derived', 'refused', or 'error'.
 * @typedef {object} RegionDialogView
 * @property {boolean} show
 * @property {RegionDialogRow[]} rows
 * @property {boolean} keepValues "Keep my values": derive with `adaptSizes: false`.
 * @property {'idle' | 'invalid' | 'deriving' | 'derived' | 'refused' | 'error'} status
 * @property {string} message The reason shown for 'idle', 'invalid', 'refused' and 'error'.
 * @property {{ setting: string, label: string, explicit: string, auto: string }[]} returnedToAuto
 * @property {string[]} notCarried
 */

/** @type {Readonly<Record<RegionNotCarriedItem['kind'], string>>} One English line per kind (wording owner). */
const NOT_CARRIED_TEXT = Object.freeze({
  OTHER_MODE_SETTINGS: msg('Settings of the other mode'),
  FEATURE_PLACEMENTS: msg('Feature placements'),
  FEATURE_PLACEMENT_TABLE: msg('The feature placement table'),
  SIMILARITY_ALIGNMENT: msg('The Similarity alignment'),
  CONSERVATION_RINGS: msg('Conservation rings'),
  CONSERVATION_RING_TABLES: msg('The conservation ring tables'),
  LINEAR_COMPARISONS: msg('Linear comparisons'),
  COMPARISON_TABLES: msg('The uploaded comparison tables'),
  LOSAT_NO_PAIR: msg('The LOSAT comparison (no pair of drawn records)'),
  ANNOTATION_TABLE: msg('The annotation table'),
  ANNOTATION_CROSSES_ORIGIN: msg('Annotations that cross the origin'),
  ANNOTATION_DRAWN_COORDINATES: msg('Annotations in drawn coordinates'),
  ANNOTATIONS_OUTSIDE: msg('Annotations outside the regions')
});

/** @param {unknown} value */
const shownValue = (value) => (value === null || value === undefined ? t('Auto') : String(value));

/** @param {string} text @returns {number | null} */
const wholeNumber = (text) => (/^\s*[1-9]\d*\s*$/.test(text) ? Number(text) : null);

/**
 * The regions of the checked rows, or the reason they are not regions yet.
 * @param {readonly RegionDialogRow[]} rows
 * @returns {{ regions: RegionDrawingRegion[] } | { reason: string, status: 'idle' | 'invalid' }}
 */
const regionsOf = (rows) => {
  const checked = rows.filter((row) => row.checked);
  if (checked.length === 0) return { status: 'idle', reason: t('Choose at least one record.') };
  /** @type {RegionDrawingRegion[]} */
  const regions = [];
  for (const row of checked) {
    if (row.kind === 'whole') regions.push({ recordKey: row.recordKey, start: 1, end: row.length });
    else if (row.kind === 'feature') {
      if (!row.featureId) return { status: 'invalid', reason: t('Choose a feature of {name}.', { name: row.caption }) };
      regions.push({ recordKey: row.recordKey, featureId: row.featureId });
    } else {
      const start = wholeNumber(row.start);
      const end = wholeNumber(row.end);
      if (start === null || end === null) {
        return { status: 'invalid', reason: t('Enter a start and an end of {name} as whole numbers.', { name: row.caption }) };
      }
      regions.push({ recordKey: row.recordKey, start, end });
    }
  }
  return { regions };
};

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
export const createRegionDrawing = ({
  shownDrawingId, projectShownRequest, deriveRegion, newDrawingFromRegion, shownRecords, settingLabel
}) => {
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

  // The dialog (the one owner of its open, close and reset).
  /** @type {RegionDialogView} */
  const dialog = reactive({
    show: false, rows: [], keepValues: false, status: 'idle', message: '', returnedToAuto: [], notCarried: []
  });
  /** @type {(() => Promise<DrawingListOutcome>) | null} */
  let create = null;
  let request = 0;

  /** @param {RegionDialogView['status']} status @param {string} message */
  const showStatus = (status, message) => {
    create = null;
    Object.assign(dialog, { status, message, returnedToAuto: [], notCarried: [] });
  };

  /** Derives again from the inputs; a reply that a newer request superseded is dropped. */
  const refresh = async () => {
    const mine = ++request;
    const choice = regionsOf(dialog.rows);
    if (!('regions' in choice)) return showStatus(choice.status, choice.reason);
    showStatus('deriving', '');
    try {
      const offer = await deriveFromShown(choice.regions, !dialog.keepValues);
      if (mine !== request) return undefined;
      if (offer.status === 'refused') {
        const refusal = normalizeUserFacingError(diagnosticError('REGION_INVALID', { field: 'region', reason: offer.reason }));
        return showStatus('refused', refusal ? errorText(refusal.summaryParts) : '');
      }
      Object.assign(dialog, {
        status: 'derived',
        message: '',
        returnedToAuto: offer.returnedToAuto.map((item) => ({
          setting: item.setting, label: settingLabel(item.setting), explicit: shownValue(item.explicit), auto: shownValue(item.auto)
        })),
        notCarried: offer.notCarried.map((item) => tm(NOT_CARRIED_TEXT[item.kind]))
      });
      create = offer.create;
    } catch (error) {
      if (mine !== request) return undefined;
      const normalized = normalizeUserFacingError(error);
      showStatus('error', normalized ? errorText(normalized.summaryParts) : String(error));
    }
    return undefined;
  };

  const close = () => {
    request += 1;
    create = null;
    Object.assign(dialog, {
      show: false, rows: [], keepValues: false, status: 'idle', message: '', returnedToAuto: [], notCarried: []
    });
  };

  const open = () => {
    close();
    dialog.rows = shownRecords().map((record) => ({
      recordKey: record.recordKey, caption: record.caption, length: record.length, checked: false, kind: 'whole',
      start: '1', end: String(record.length), featureId: '', features: record.features
    }));
    dialog.show = true;
    showStatus('idle', t('Choose at least one record.'));
  };

  /** Create: the derived drawing as one project step; the dialog closes. */
  const confirm = async () => {
    if (dialog.status !== 'derived' || !create) return null;
    const outcome = await create();
    close();
    return outcome;
  };

  return { deriveFromShown, dialog, open, refresh, confirm, close };
};
