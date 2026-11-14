// @ts-check
import { diagnosticError, normalizeUserFacingError, reasonText } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { msg, numberText, t, tm } from '../utils/i18n.js';
import { reactive } from '../utils/reactivity.js';
import { drawnRecordIdentities } from '../services/record-draw-selection.js';
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';
import { copyRefusalLine } from './region-drawing/copy-text.js';
import { createLastTarget, pasteLinesOf, pastePayloadOf, pasteTargetsOf } from './region-drawing/paste.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ProjectedDraftRequest, RegionDerivation } from '../services/config.js' */
/** @import { DerivedRegionDrawing, NewDrawingCopy, NotCarriedItem, PastedRecords, RefusedCopy } from '../contracts/locate.generated.js' */
/** @import { CopyTarget, PasteTargetPorts } from './region-drawing/paste.js' */
/** @import { PageExistingDrawingCopy, RefusedRegionDrawing, RegionDrawingHelperPayload, RegionDrawingRegion, RegionDrawingReply } from '../services/diagram-worker-protocol.js' */
/** @import { DrawingListOutcome } from './drawing-list.js' */
/** @import { IdentityRequestRecord, RecordIdentitySources } from '../services/record-draw-selection.js' */

/**
 * @typedef {object} RegionDrawingPorts
 * @property {() => DrawingId} shownDrawingId
 * @property {() => Promise<ProjectedDraftRequest>} projectShownRequest The shown drawing's draft request.
 * @property {(payload: RegionDrawingHelperPayload) => Promise<{ result?: unknown }>} deriveRegion
 *   Runs the Worker helper `deriveRegionDrawing`; the reply is read here.
 * @property {(sourceId: DrawingId, derivation: RegionDerivation) => Promise<DrawingListOutcome>} newDrawingFromRegion
 * @property {() => Promise<RecordIdentitySources | null>} readRecordSources The shown drawing's discovered record
 *   lists, read first when a Load left them deferred; null without the record discovery owner.
 * @property {(recordKey: string) => readonly RegionSourceFeature[] | null} resultFeatures The features the shown
 *   Result drew under a record key (the committed request's), or null when it drew none under that key.
 * @property {(drawn: DrawnRecordChoice) => string | null} committedRecordKeyOf The committed request's record key of a
 *   drawn record (its definition group's index and ID), found as the record popup finds it; null when none.
 * @property {PasteTargetPorts} pasteTargets The existing drawings a copy can paste into.
 */

/** @typedef {{ recordIndex: number, recordId: string }} DrawnRecordChoice A drawn record (the record popup's identity). */

/**
 * A record of the source drawing as the dialog offers it. `features` are the
 * shown Result's admitted features of the record, or `null` without a Result.
 * @typedef {object} RegionSourceRecord
 * @property {string} recordKey The record key of the source drawing's request (the identity).
 * @property {readonly string[]} aliases The other keys the record has (`drawnRecordIdentities`): a committed request
 *   may key it so.
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
 * @param {(recordKey: string) => readonly RegionSourceFeature[] | null} featuresOf The features the shown Result
 *   drew under a record key, or null when it drew no record under that key; asked for the record's request key and
 *   then its aliases (`drawnRecordIdentities`).
 * @returns {RegionSourceRecord[]}
 */
export const regionSourceRecords = (requestRecords, sources, featuresOf) => drawnRecordIdentities(requestRecords, sources)
  .flatMap(({ requestKey, aliases, entry }) => {
    const length = Number(entry.recordLength);
    if (!Number.isSafeInteger(length) || length <= 0) return [];
    const caption = String(entry.recordId ?? entry.record_id ?? requestKey);
    // The Result keys the record as the committed request did: its draft key, or a key it has under the other
    // canvas (its aliases).
    const features = [...new Set([requestKey, ...aliases])].map(featuresOf).find((found) => found !== null) ?? null;
    return [{ recordKey: requestKey, aliases, caption, length, features }];
  });

/**
 * One record row of the dialog: display values and the user's input only.
 * @typedef {object} RegionDialogRow
 * @property {string} recordKey
 * @property {readonly string[]} aliases
 * @property {string} caption
 * @property {number} length
 * @property {boolean} checked
 * @property {'whole' | 'coordinates' | 'feature'} kind
 * @property {string} start
 * @property {string} end
 * @property {string} featureId
 * @property {readonly RegionSourceFeature[] | null} features
 * @property {string} featureReason Why A feature cannot be chosen, or '' when it can.
 */

/**
 * The user's input to one row; the other row fields are display values.
 * @typedef {Partial<Pick<RegionDialogRow, 'checked' | 'kind' | 'start' | 'end' | 'featureId'>>} RegionRowEdit
 */

/**
 * The dialog's view. `status`: 'idle' (no region yet), 'invalid' (an input
 * is not a region), 'deriving', 'derived', 'creating' (Create runs; the derived
 * view stays shown), 'refused', or 'error'.
 * @typedef {object} RegionDialogView
 * @property {boolean} show
 * @property {RegionDialogRow[]} rows
 * @property {boolean} keepValues "Keep my values": derive with `adaptSizes: false`.
 * @property {NewDrawingCopy['mode']} newMode The mode the user chose for the new drawing; null takes Python's default
 *   (the source mode for whole records only, else Linear: PD-OI-090 (4)).
 * @property {NewDrawingCopy['settings']} settings "Keep this drawing's settings and edits" or "Start from the mode defaults".
 * @property {DerivedRegionDrawing['mode'] | null} mode The mode of the derived drawing, shown in the mode choice.
 * @property {'idle' | 'invalid' | 'deriving' | 'derived' | 'creating' | 'refused' | 'error'} status
 * @property {string} message The reason shown for 'idle', 'invalid', 'refused' and 'error'.
 * @property {{ setting: string, label: string, explicit: string, auto: string }[]} returnedToAuto
 * @property {string[]} notCarried
 * @property {CopyTarget} target A new drawing, or the existing drawing the copy pastes into.
 * @property {{ id: DrawingId, name: string }[]} targets The existing drawings offered: every drawing but the shown one.
 * @property {PageExistingDrawingCopy['edits']} edits "Bring the records' edits" or "Records only" (an existing drawing).
 * @property {PageExistingDrawingCopy['duplicates']} duplicates "Skip it" or "Add it again" for what the drawing has already.
 * @property {string[]} alreadyPresent The lines of what the existing drawing has already.
 * @property {string[]} leftOut The lines of what the existing drawing leaves out, with the reason.
 */

/**
 * The English name of each setting that can return to Auto, by setting ID
 * (`autoSizes.adaptableSettings` of `mode-profiles.generated.js`, from
 * `ADAPTABLE_SETTINGS` in `gbdraw/auto_sizes.py`): the accessible-name key of
 * the settings pane's control that writes it, or a name of its own where no
 * control sets it or where the control's text needs its section to be read
 * (Window, Tick Font Size). A setting that each mode names differently has one
 * text per mode. A twin of the pane labels: region-drawing-dialog.test.mjs
 * checks each against the control that writes the setting.
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
  'objects.definition.linear.font_size': msg('Definition font size'),
  'objects.scale.font_size': msg('Linear scale font size'),
  'objects.scale.ruler_label_font_size': msg('Linear ruler label font size'),
  window: msg('GC content/skew window'),
  step: msg('GC content/skew step'),
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
  'objects.gc_content.tick_font_size': msg('GC content tick font size'),
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
 * One English line per kind of `regionDrawing.notCarriedKinds` (from Python's
 * `NotCarriedKind`; wording owner); a kind worded per mode is keyed,
 * as `ADAPTABLE_SETTING_TEXT` is, by the new drawing's mode. The other mode's
 * settings are those of the mode the new drawing is not.
 * @type {Readonly<Record<NotCarriedItem['kind'], string | Readonly<Record<DerivedRegionDrawing['mode'], string>>>>}
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
 * @param {NotCarriedItem['kind']} kind
 * @returns {string}
 */
export const notCarriedLine = (mode, kind) => {
  const text = NOT_CARRIED_TEXT[kind];
  return tm(typeof text === 'string' ? text : text[mode]);
};

/** A returned size's value: Auto, or a number grouped as the row's length is. @param {unknown} value */
const shownValue = (value) => {
  if (value === null || value === undefined) return t('Auto');
  return typeof value === 'number' ? numberText(value) : String(value);
};

/**
 * What a derived drawing's dialog lists: the sizes that return to Auto and the lines of what is not carried.
 * Copy to another drawing and Find in other files word them alike.
 * @param {Pick<DerivedRegionDrawing, 'mode' | 'returnedToAuto' | 'notCarried'>} derived
 * @returns {Pick<RegionDialogView, 'returnedToAuto' | 'notCarried'>}
 */
export const derivedDrawingNotes = ({ mode, returnedToAuto, notCarried }) => ({
  returnedToAuto: returnedToAuto.map((item) => ({
    setting: item.setting, label: adaptableSettingLabel(mode, item.setting), explicit: shownValue(item.explicit), auto: shownValue(item.auto)
  })),
  notCarried: notCarried.map((item) => notCarriedLine(mode, item.kind))
});

/** The lines of a derived drawing or a paste, none shown. */
const NO_LINES = Object.freeze({ returnedToAuto: [], notCarried: [], alreadyPresent: [], leftOut: [] });
/** The dialog as it opens and after Close. @returns {Omit<RegionDialogView, 'show' | 'rows'>} */
const startingChoices = () => ({
  keepValues: false, newMode: null, settings: 'keep', mode: null, status: 'idle', message: '', ...NO_LINES,
  target: 'new', targets: [], edits: 'edits', duplicates: 'skip'
});

/**
 * The user's choices of the dialog that are not rows: Keep my values, the target, the new drawing's mode and settings,
 * and an existing drawing's edits and duplicates.
 * @typedef {Partial<Pick<RegionDialogView, 'keepValues' | 'newMode' | 'settings' | 'target' | 'edits' | 'duplicates'>>} RegionChoiceEdit
 */

/** A copy as the dialog's choices make it: a new drawing, or a paste into the drawing of `drawingId`. */
/** @typedef {NewDrawingCopy | ({ into: 'existing', drawingId: DrawingId } & Pick<PageExistingDrawingCopy, 'edits' | 'duplicates'>)} CopyChoice */

/** The copy the dialog's choices make. @param {RegionDialogView} dialog @returns {CopyChoice} */
const copyOf = ({ target, settings, newMode, edits, duplicates }) => (target === 'new'
  ? { into: 'new', settings, mode: newMode }
  : { into: 'existing', drawingId: target, edits, duplicates });

/** The reason a refusal shows, from the wording owner of its vocabulary. @param {RefusedRegionDrawing | RefusedCopy} offer */
const refusalText = (offer) => (offer.reason === 'START_AFTER_END'
  ? tm(/** @type {string} */ (reasonText(offer.reason)))
  : copyRefusalLine(offer.reason));

// A feature is chosen from the features the shown Result drew of the record.
/** @param {readonly RegionSourceFeature[] | null} features */
const featureReasonOf = (features) => {
  if (features === null) return t('Generate the drawing to choose a feature.');
  return features.length === 0 ? t('This record has no features to choose.') : '';
};

/**
 * The dialog's rows of `records`; a record that `chosen` names by its key or an alias takes the chosen edit.
 * @param {readonly RegionSourceRecord[]} records @param {{ keys: readonly string[], edit: RegionRowEdit }} chosen
 * @returns {RegionDialogRow[]}
 */
const dialogRowsOf = (records, chosen) => records.map((record) => ({
  recordKey: record.recordKey, aliases: record.aliases, caption: record.caption, length: record.length, checked: false,
  kind: /** @type {RegionDialogRow['kind']} */ ('whole'), start: '1', end: String(record.length), featureId: '',
  features: record.features, featureReason: featureReasonOf(record.features),
  ...([record.recordKey, ...record.aliases].some((key) => chosen.keys.includes(key)) ? chosen.edit : {})
}));

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
 * @typedef {RefusedRegionDrawing | RefusedCopy | ((DerivedRegionDrawing | PastedRecords) & { create: () => Promise<unknown> })} RegionDrawingOffer
 */

const REPLY_INVALID = () => diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: 'deriveRegionDrawing' });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/** The reply's vocabularies, from Python's `StartAfterEndReason`, `CopyRefusalReason` and `NotCarriedKind`. */
const {
  refusalReasons: REFUSAL_REASONS, copyRefusalReasons: COPY_REFUSAL_REASONS, notCarriedKinds: NOT_CARRIED_KINDS,
  pasteNotCarriedKinds: PASTE_NOT_CARRIED_KINDS
} = MODE_PROFILE_DATA.regionDrawing;

/** @param {unknown} value @returns {value is Record<string, string>} */
const isIdTable = (value) => isRecord(value) && Object.values(value).every((id) => typeof id === 'string');

/** A paste's reply, checked against `PastedRecords`. @param {Record<string, unknown>} result @returns {PastedRecords} */
const readPastedRecords = (result) => {
  const lists = [result.pastedRecordKeys, result.alreadyPresent, result.leftOut, result.notCarried, result.addedResources];
  if ((result.mode !== 'circular' && result.mode !== 'linear') || !isRecord(result.renderRequest)
    || !isIdTable(result.resourceSources) || !isIdTable(result.targetResourceSources) || !lists.every(Array.isArray)
    || !/** @type {unknown[]} */ (result.leftOut).every((item) => isRecord(item) && COPY_REFUSAL_REASONS.includes(String(item.reason)))
    || !/** @type {unknown[]} */ (result.notCarried).every((item) => isRecord(item) && PASTE_NOT_CARRIED_KINDS.includes(String(item.kind)))) {
    throw REPLY_INVALID();
  }
  return /** @type {PastedRecords} */ (/** @type {unknown} */ (result));
};

/**
 * The Worker's reply, checked against `RegionDrawingReply` or a refused copy (`RefusedCopy`).
 * @param {unknown} result
 * @returns {RegionDrawingReply | RefusedCopy | PastedRecords}
 */
const readRegionDrawingReply = (result) => {
  if (!isRecord(result)) throw REPLY_INVALID();
  if (result.status === 'pasted') return readPastedRecords(result);
  if (result.status === 'refused' && Object.keys(result).length === 2) {
    if (typeof result.reason !== 'string' || !COPY_REFUSAL_REASONS.includes(result.reason)) throw REPLY_INVALID();
    return /** @type {RefusedCopy} */ (result);
  }
  if (result.status === 'refused') {
    if (typeof result.reason !== 'string' || !REFUSAL_REASONS.includes(result.reason) || typeof result.recordKey !== 'string'
      || !Number.isSafeInteger(result.start) || !Number.isSafeInteger(result.end)) throw REPLY_INVALID();
    return /** @type {RefusedRegionDrawing} */ (result);
  }
  const { resourceSources } = result;
  if (result.status !== 'derived' || (result.mode !== 'circular' && result.mode !== 'linear')
    || !isRecord(result.renderRequest) || result.renderRequest.mode !== result.mode
    || !isRecord(resourceSources) || !Object.values(resourceSources).every((id) => typeof id === 'string')
    || !Array.isArray(result.returnedToAuto) || !Array.isArray(result.notCarried)
    || !result.notCarried.every((item) => isRecord(item) && typeof item.kind === 'string' && NOT_CARRIED_KINDS.includes(item.kind))) {
    throw REPLY_INVALID();
  }
  return /** @type {DerivedRegionDrawing} */ (/** @type {unknown} */ (result));
};

/**
 * The offer of a copy of `regions` of the shown drawing: Python derives the new drawing, or pastes into an existing
 * drawing, whose draft request travels with the source's (`pastePayloadOf`).
 * @param {RegionDrawingPorts} ports
 * @param {readonly RegionDrawingRegion[]} regions Source coordinates, 1-based and inclusive.
 * @param {boolean} adaptSizes Whether hand-set sizes whose Auto differs return to Auto.
 * @param {CopyChoice} choice Where the regions go and with which settings or edits.
 * @returns {Promise<RegionDrawingOffer>}
 */
const deriveCopy = async (ports, regions, adaptSizes, choice) => {
  const { shownDrawingId, projectShownRequest, deriveRegion, newDrawingFromRegion, pasteTargets } = ports;
  const sourceId = shownDrawingId();
  const sent = await projectShownRequest();
  const resources = /** @type {RegionDrawingHelperPayload['resources']} */ (sent.envelope.resources);
  const into = choice.into === 'existing' ? { choice, targetSent: await pasteTargets.projectDrawing(choice.drawingId) } : null;
  const payload = into
    ? pastePayloadOf(resources, into.targetSent, { edits: into.choice.edits, duplicates: into.choice.duplicates })
    : { resources, copy: /** @type {NewDrawingCopy} */ (choice) };
  const reply = readRegionDrawingReply((await deriveRegion({ canonicalRequest: sent.envelope.renderRequest, regions, adaptSizes, ...payload })).result);
  if (reply.status === 'refused') return reply;
  if (reply.status !== 'pasted') return { ...reply, create: () => newDrawingFromRegion(sourceId, { derived: reply, sent }) };
  if (!into) throw REPLY_INVALID();
  return { ...reply, create: () => pasteTargets.pasteInto(into.choice.drawingId, { pasted: reply, sent, targetSent: into.targetSent }) };
};

/**
 * Copy to another drawing (PD-OI-090 〔作成〕(4)): Python derives the
 * drawing from the shown drawing's draft request; Copy adds it with the
 * Web files of the request that was sent, or pastes it into an existing drawing.
 * @param {RegionDrawingPorts} ports
 */
export const createRegionDrawing = (ports) => {
  const { shownDrawingId, projectShownRequest, readRecordSources, resultFeatures, committedRecordKeyOf, pasteTargets } = ports;
  /** @param {readonly RegionDrawingRegion[]} regions @param {boolean} adaptSizes @param {CopyChoice} choice */
  const deriveFromShown = (regions, adaptSizes, choice) => deriveCopy(ports, regions, adaptSizes, choice);
  const lastTarget = createLastTarget(pasteTargets.drawingOf);

  // The dialog (the one owner of its open, close, reset and input edits; the template assigns none of its state).
  /** @type {RegionDialogView} */
  const dialog = reactive({ show: false, rows: [], ...startingChoices() });
  /** @type {(() => Promise<unknown>) | null} */
  let create = null;
  let request = 0;

  /** @param {RegionDialogView['status']} status @param {string} message */
  const showStatus = (status, message) => {
    create = null;
    Object.assign(dialog, { status, message, ...NO_LINES });
  };

  /** Derives again from the inputs; a reply that a newer request superseded is dropped. None while Create runs. */
  const refresh = async () => {
    if (dialog.status === 'creating') return undefined;
    const mine = ++request;
    const choice = regionsOf(dialog.rows);
    if (!('regions' in choice)) return showStatus(choice.status, choice.reason);
    showStatus('deriving', '');
    try {
      const offer = await deriveFromShown(choice.regions, !dialog.keepValues, copyOf(dialog));
      if (mine !== request) return undefined;
      // The reason alone, from the wording owner of its vocabulary: the dialog's rows are the region.
      if (offer.status === 'refused') return showStatus('refused', refusalText(offer));
      Object.assign(dialog, { status: 'derived', message: '', mode: offer.mode, ...NO_LINES, ...(offer.status === 'pasted'
        ? pasteLinesOf(offer, captionOf, dialog.targets.find(({ id }) => id === dialog.target)?.name ?? '')
        : derivedDrawingNotes(offer)) });
      create = offer.create;
    } catch (error) {
      if (mine !== request) return undefined;
      const normalized = normalizeUserFacingError(error);
      showStatus('error', normalized ? errorText(normalized.summaryParts) : String(error));
    }
    return undefined;
  };

  /**
   * A row edit: writes the user's input into the row found by its record key and derives again; an unknown key
   * changes nothing.
   * @param {string} recordKey
   * @param {RegionRowEdit} patch
   */
  const editRow = (recordKey, patch) => {
    const row = dialog.rows.find((candidate) => candidate.recordKey === recordKey);
    if (!row) return Promise.resolve(undefined);
    Object.assign(row, patch);
    return refresh();
  };

  /** A source record's caption, by its record key; the key when no row names it. @param {string} recordKey */
  const captionOf = (recordKey) => dialog.rows.find((row) => row.recordKey === recordKey)?.caption ?? recordKey;

  /** A choice that is not a row (Keep my values, the target, the mode, the settings): derives again. @param {RegionChoiceEdit} edit */
  const setChoice = (edit) => {
    Object.assign(dialog, edit);
    return refresh();
  };

  const close = () => {
    request += 1;
    create = null;
    Object.assign(dialog, { show: false, rows: [], ...startingChoices() });
  };

  /**
   * Opens on the shown drawing's records, those `chosen` names (by record key) with its edit, and derives when one is
   * chosen; a close before they are read keeps it closed.
   * @param {{ keys: readonly string[], edit: RegionRowEdit }} [chosen]
   */
  const open = async (chosen = { keys: [], edit: {} }) => {
    close();
    const mine = request;
    dialog.show = true;
    showStatus('deriving', '');
    let records;
    try {
      // The read writes the record list that the projection reads, so it finishes first (as for Generate).
      const sources = await readRecordSources();
      const sent = sources ? await projectShownRequest() : null;
      records = sources && sent ? regionSourceRecords(sent.envelope.renderRequest.records, sources, resultFeatures) : [];
    } catch (error) {
      if (mine !== request) return;
      const normalized = normalizeUserFacingError(error);
      showStatus('error', normalized ? errorText(normalized.summaryParts) : String(error));
      return;
    }
    if (mine !== request) return undefined;
    dialog.rows = dialogRowsOf(records, chosen);
    dialog.targets = pasteTargetsOf(pasteTargets.tabs(), shownDrawingId());
    dialog.target = lastTarget.read(dialog.targets);
    return refresh();
  };

  /**
   * Opens the dialog on a region chosen on the figure (Select region): its record checked with its coordinates, so
   * Copy makes it as the dialog makes any region. A record the dialog does not offer is not chosen.
   * @param {{ recordKey: string, start: number, end: number }} region Source coordinates, 1-based and inclusive.
   */
  const openWithRegion = ({ recordKey, start, end }) => open({
    keys: [recordKey], edit: { checked: true, kind: 'coordinates', start: String(start), end: String(end) }
  });

  /**
   * Opens the dialog with drawn records chosen whole (the record popup, the records of a Shift+drag), each found by
   * its committed record key, never by its caption. @param {readonly DrawnRecordChoice[]} drawn
   */
  const openWithRecords = (drawn) => open({
    keys: drawn.map(committedRecordKeyOf).filter(/** @returns {key is string} */ (key) => key !== null), edit: { checked: true, kind: 'whole' }
  });

  /**
   * Create: the derived drawing as one project step; the dialog closes. While it runs the dialog is 'creating',
   * so a second Create adds nothing; a failure shows its reason and leaves the dialog open.
   */
  const confirm = async () => {
    if (dialog.status !== 'derived' || !create) return null;
    const run = create;
    const { target } = dialog;
    create = null;
    dialog.status = 'creating';
    try {
      const outcome = await run();
      lastTarget.remember(target);
      close();
      return outcome;
    } catch (error) {
      const normalized = normalizeUserFacingError(error);
      showStatus('error', normalized ? errorText(normalized.summaryParts) : String(error));
      return null;
    }
  };

  return { deriveFromShown, dialog, open, openWithRegion, openWithRecords, refresh, editRow, setChoice, confirm, close };
};
