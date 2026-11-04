// @ts-check
import { diagnosticError, normalizeUserFacingError, reasonText } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { msg, t, tm } from '../utils/i18n.js';
import { reactive } from '../utils/reactivity.js';
import { drawnRecordIdentities } from '../services/record-draw-selection.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ProjectedDraftRequest, RegionDerivation } from '../services/config.js' */
/** @import { DerivedRegionDrawing, RefusedRegionDrawing, RegionNotCarriedItem, RegionDrawingHelperPayload, RegionDrawingRegion, RegionDrawingReply } from '../services/diagram-worker-protocol.js' */
/** @import { DrawingListOutcome } from './drawing-list.js' */
/** @import { IdentityRequestRecord, RecordIdentitySources } from '../services/record-draw-selection.js' */

/**
 * @typedef {object} RegionDrawingPorts
 * @property {() => DrawingId} shownDrawingId
 * @property {() => Promise<ProjectedDraftRequest>} projectShownRequest The shown drawing's draft request.
 * @property {(payload: RegionDrawingHelperPayload) => Promise<{ result?: unknown }>} deriveRegion
 *   Runs the Worker helper `deriveRegionDrawing`; the reply is read here.
 * @property {(sourceId: DrawingId, derivation: RegionDerivation) => Promise<DrawingListOutcome>} newDrawingFromRegion
 * @property {() => Promise<readonly RegionSourceRecord[]>} shownRecords The shown drawing's records, in drawing
 *   order (`regionSourceRecords`).
 */

/**
 * A record of the source drawing as the dialog offers it. `features` are the
 * shown Result's admitted features of the record, or `null` without a Result.
 * @typedef {object} RegionSourceRecord
 * @property {string} recordKey The record key of the source drawing's request (the identity).
 * @property {string} caption Display only.
 * @property {number} length Source length in bases.
 * @property {readonly RegionSourceFeature[] | null} features
 */

/** @typedef {{ featureId: string, label: string }} RegionSourceFeature A feature by its biological feature ID. */

/**
 * The records the dialog offers: the drawn records of the shown drawing's
 * request, keyed as that request keys them, from the one owner of a request
 * record's source entry (`drawnRecordIdentities`), with the discovered length.
 * A record of unknown length is not offered (its whole span is not known).
 * @param {readonly IdentityRequestRecord[]} requestRecords The projected request's records.
 * @param {RecordIdentitySources} sources The discovered record lists of the drawing.
 * @param {(recordKey: string) => readonly RegionSourceFeature[] | null} featuresOf A record's features in the
 *   shown Result, or null without a Result that drew it.
 * @returns {RegionSourceRecord[]}
 */
export const regionSourceRecords = (requestRecords, sources, featuresOf) => drawnRecordIdentities(requestRecords, sources)
  .flatMap(({ requestKey, entry }) => {
    const length = Number(entry.recordLength);
    if (!Number.isSafeInteger(length) || length <= 0) return [];
    const caption = String(entry.recordId ?? entry.record_id ?? requestKey);
    return [{ recordKey: requestKey, caption, length, features: featuresOf(requestKey) }];
  });

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
 * @property {readonly RegionSourceFeature[] | null} features
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

/**
 * The English name of each setting that can return to Auto, by setting ID
 * (`autoSizes.adaptableSettings` of `mode-profiles.generated.js`, from
 * `ADAPTABLE_SETTINGS` in `gbdraw/auto_sizes.py`): the text of the settings
 * pane's control for it, or a name of its own where no control sets it. A
 * setting that each mode names differently has one text per mode.
 * @type {Readonly<Record<string, string | Readonly<Record<DerivedRegionDrawing['mode'], string>>>>}
 */
const ADAPTABLE_SETTING_TEXT = Object.freeze({
  'objects.features.block_stroke_width': msg('Block Stroke Width'),
  'objects.features.line_stroke_width': msg('Line Stroke Width'),
  'objects.legends.color_rect_size': msg('Legend swatch size'),
  'objects.legends.font_size': msg('Legend font size'),
  'labels.font_size': msg('Label Font Size'),
  'labels.font_size.linear': msg('Label Font Size'),
  'objects.axis.circular.stroke_width': msg('Axis Stroke Width'),
  'objects.axis.linear.stroke_width': msg('Axis Stroke Width'),
  'canvas.linear.arrow_length_parameter': msg('Linear arrow length'),
  'canvas.linear.default_cds_height': msg('Feature Height'),
  'objects.definition.linear.font_size': msg('Default font size'),
  'objects.scale.font_size': msg('Scale Font Size'),
  'objects.scale.ruler_label_font_size': msg('Ruler Label Font Size'),
  window: msg('Window'),
  step: msg('Step'),
  depth_window: msg('Depth Window'),
  depth_step: msg('Depth Step'),
  'circularSlot.features.width': msg('Feature Track Width'),
  'circularSlot.gc_content.width': msg('GC Content Width'),
  'circularSlot.gc_skew.width': msg('GC Skew Width'),
  'circularSlot.depth.width': msg('Depth Width'),
  'objects.scale.interval': Object.freeze({
    circular: msg('Circular scale interval'), linear: msg('Linear ruler interval')
  }),
  'canvas.linear.default_gc_height': msg('GC Height'),
  'canvas.linear.depth_height': msg('Default Depth Height'),
  'objects.gc_content.tick_font_size': msg('Tick Font Size'),
  'objects.depth.tick_font_size': msg('Depth tick font size')
});

/**
 * A setting's name in the UI language, as a drawing of `mode` names it; an
 * unknown setting shows its ID.
 * @param {DerivedRegionDrawing['mode']} mode
 * @param {string} setting A setting ID of `autoSizes.adaptableSettings`.
 * @returns {string}
 */
export const adaptableSettingLabel = (mode, setting) => {
  const text = Object.hasOwn(ADAPTABLE_SETTING_TEXT, setting) ? ADAPTABLE_SETTING_TEXT[setting] : null;
  if (text === null) return setting;
  return tm(typeof text === 'string' ? text : text[mode]);
};

/**
 * One English line per kind (wording owner); a kind worded per mode is keyed,
 * as `ADAPTABLE_SETTING_TEXT` is, by the new drawing's mode. The other mode's
 * settings are those of the mode the new drawing is not.
 * @type {Readonly<Record<RegionNotCarriedItem['kind'], string | Readonly<Record<DerivedRegionDrawing['mode'], string>>>>}
 */
const NOT_CARRIED_TEXT = Object.freeze({
  OTHER_MODE_SETTINGS: Object.freeze({
    linear: msg('Settings that apply only to Circular drawings'), circular: msg('Settings that apply only to Linear drawings')
  }),
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

/**
 * A not-carried line in the UI language, as a drawing of `mode` words it.
 * @param {DerivedRegionDrawing['mode']} mode
 * @param {RegionNotCarriedItem['kind']} kind
 * @returns {string}
 */
const notCarriedLine = (mode, kind) => {
  const text = NOT_CARRIED_TEXT[kind];
  return tm(typeof text === 'string' ? text : text[mode]);
};

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
    if ((result.reason !== 'ORIGIN_CROSSING' && result.reason !== 'START_AFTER_END') || typeof result.recordKey !== 'string'
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
  shownDrawingId, projectShownRequest, deriveRegion, newDrawingFromRegion, shownRecords
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
      // The reason alone, from the wording owner of diagnostic reasons: the dialog's rows are the region.
      if (offer.status === 'refused') return showStatus('refused', tm(/** @type {string} */ (reasonText(offer.reason))));
      Object.assign(dialog, {
        status: 'derived',
        message: '',
        returnedToAuto: offer.returnedToAuto.map((item) => ({
          setting: item.setting, label: adaptableSettingLabel(offer.mode, item.setting), explicit: shownValue(item.explicit), auto: shownValue(item.auto)
        })),
        notCarried: offer.notCarried.map((item) => notCarriedLine(offer.mode, item.kind))
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

  /** Opens on the shown drawing's records; a close before they are read keeps it closed. */
  const open = async () => {
    close();
    const mine = request;
    dialog.show = true;
    showStatus('deriving', '');
    let records;
    try {
      records = await shownRecords();
    } catch (error) {
      if (mine !== request) return;
      const normalized = normalizeUserFacingError(error);
      showStatus('error', normalized ? errorText(normalized.summaryParts) : String(error));
      return;
    }
    if (mine !== request) return;
    dialog.rows = records.map((record) => ({
      recordKey: record.recordKey, caption: record.caption, length: record.length, checked: false, kind: 'whole',
      start: '1', end: String(record.length), featureId: '', features: record.features
    }));
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
