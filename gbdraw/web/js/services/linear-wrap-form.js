// @ts-check
// The Linear wrap length of a drawing's form (`form.linear_wrap_length`): what
// the form holds and how a stored value reads into it. The request leaf is
// `canvas.linear.wrap_length` (services/session-request.js compiles it);
// Python decides the row length and whether a record is cut
// (gbdraw/layout/linear_wrap.py), and parses a length typed with a unit.
/** @import { InputSummary, RecordSummary } from './input-summary.js' */
/** @import { CanonicalRenderRequest } from '../contracts/request.generated.js' */

/** The request leaf of the wrap length (a config override path). */
export const LINEAR_WRAP_LENGTH_PATH = 'canvas.linear.wrap_length';

/**
 * `null`: off; `'auto'`: Python's automatic row length; a number: a row length in bp, a positive safe integer.
 * @typedef {null | 'auto' | number} LinearWrapFormValue
 */
/**
 * A stored value read into the form: `lost` when neither the form nor Python can read it, and then the value is off.
 * @typedef {{ value: LinearWrapFormValue, lost: boolean }} LinearWrapReading
 */
/**
 * Python's readings of the stored wrap lengths the form cannot hold as stored (a length typed with a unit such as
 * `"5kb"`, `"Auto"`, `" 20000 "`), by `wrapReadingKey` of the stored value. Load asks Python
 * (`linear_wrap_length_from_config`, through the Worker's configuration check) before its readers run; a value
 * Python refuses has no reading.
 * @typedef {ReadonlyMap<string, LinearWrapFormValue>} StoredWrapReadings
 */
/** @type {StoredWrapReadings} */
export const NO_WRAP_READINGS = new Map();
/** @param {unknown} raw */
export const wrapReadingKey = (raw) => JSON.stringify(raw);

/**
 * Reads a stored wrap length (a Session's form value, a request or config leaf, a preserved override) as the form
 * holds it. The form's own values and Python's canonical ones read as they are: absent, `null` and `"none"` are
 * off, `"auto"` is Auto, a positive safe integer is a length in bp. Any other value is Python's to read: its
 * reading in `readings`, else lost.
 * @param {unknown} raw
 * @param {StoredWrapReadings} [readings]
 * @returns {LinearWrapReading}
 */
export const wrapFormValue = (raw, readings = NO_WRAP_READINGS) => {
  if (raw === undefined || raw === null || raw === 'none') return { value: null, lost: false };
  if (raw === 'auto') return { value: 'auto', lost: false };
  if (typeof raw === 'number' && Number.isSafeInteger(raw) && raw > 0) return { value: raw, lost: false };
  const reading = readings.get(wrapReadingKey(raw));
  return reading === undefined ? { value: null, lost: true } : { value: reading, lost: false };
};

/**
 * The wrap length a committed Linear request stores, as stored: its override leaf, or else its full config's (a CLI
 * request); null when it stores none. Python's validator reads it where the stored text decides (the export check).
 * @param {CanonicalRenderRequest} request
 * @returns {unknown}
 */
export const requestLinearWrapLeaf = (request) => {
  const { configOverrides, config } = /** @type {{ configOverrides?: Record<string, unknown> | null,
    config?: { canvas?: { linear?: { wrap_length?: unknown } } } | null }} */ (request.diagramOptions);
  const leaf = configOverrides && Object.hasOwn(configOverrides, LINEAR_WRAP_LENGTH_PATH)
    ? configOverrides[LINEAR_WRAP_LENGTH_PATH]
    : config?.canvas?.linear?.wrap_length;
  return leaf ?? null;
};

/**
 * The Wrap Rows value a committed Linear request asked for, as the form holds it (`requestLinearWrapLeaf`). A value
 * the form cannot hold reads as null.
 * @param {CanonicalRenderRequest} request
 * @returns {LinearWrapFormValue}
 */
export const requestLinearWrapLength = (request) => wrapFormValue(requestLinearWrapLeaf(request)).value;

// Turning Wrap Rows on switches the settings rows cannot draw with, and the one rows decide (receipt PD-OI-113:
// each row draws its ruler on its axis), in the same edit: Features on axis to Above, Bar (Simple) to Ruler
// (Ticks), Normalize Record Lengths off, Ruler on Axis on. It does so also when no record is cut, so the figure is
// the form's visible values. Turning it off writes only the length; Undo restores the others.
/**
 * @typedef {{ linear_wrap_length: LinearWrapFormValue | number, linear_track_layout?: 'above' | 'below',
 *   scale_style?: 'ruler', normalize_length?: false, linear_ruler_on_axis?: true }} LinearWrapFormPatch
 *   A length that is no positive safe integer (a typed `1.5`) is sent as typed, and Python refuses it.
 */
/** The form fields a Wrap Rows edit writes; the layout transition saves and restores them. */
export const LINEAR_WRAP_FORM_FIELDS = Object.freeze(
  /** @type {const} */ (['linear_wrap_length', 'linear_track_layout', 'scale_style', 'normalize_length', 'linear_ruler_on_axis'])
);

/**
 * The fields one Wrap Rows edit writes.
 * @param {{ linear_track_layout: unknown }} form The drawing's form before the edit.
 * @param {LinearWrapFormValue | number} value
 * @returns {LinearWrapFormPatch}
 */
export const linearWrapFormPatch = (form, value) => (value === null
  ? { linear_wrap_length: null }
  : {
      linear_wrap_length: value,
      linear_track_layout: form.linear_track_layout === 'below' ? 'below' : 'above',
      scale_style: 'ruler',
      normalize_length: false,
      linear_ruler_on_axis: true
    });

/**
 * What the shown Linear drawing holds that rows cannot draw with, each read from the owner its request reads.
 * @typedef {object} LinearWrapBlockerInputs
 * @property {boolean} comparisons BLAST, comparison pairs, LOSAT, Similarity groups or a Similarity alignment.
 * @property {boolean} sharedRecordRow Record Layout puts several records in one row.
 * @property {boolean} featureAlignment Records where Align at feature put them (its `applied`), also when it moved none.
 * @property {boolean} movedRecords A drawn record moved sideways (`movesARecordSideways`), whatever moved it.
 * @property {boolean} featureTrackOnAxis Custom Track Slots put the feature track on the axis (side overlay).
 */
/**
 * Python's refusal reason for each (gbdraw/layout/linear_wrap.py), so the field's hint and the Generate refusal
 * show one text (`reasonText`).
 * @typedef {'WRAP' | 'WRAP_GRID_ROW' | 'WRAP_RECORD_POSITION' | 'WRAP_TRACK_LAYOUT'} LinearWrapBlocker
 */
/** @type {ReadonlyArray<[keyof LinearWrapBlockerInputs, LinearWrapBlocker]>} */
// Align at feature and a moved record share one reason and one hint (D-65).
const BLOCKERS = Object.freeze([
  ['comparisons', 'WRAP'], ['sharedRecordRow', 'WRAP_GRID_ROW'], ['featureAlignment', 'WRAP_RECORD_POSITION'],
  ['movedRecords', 'WRAP_RECORD_POSITION'], ['featureTrackOnAxis', 'WRAP_TRACK_LAYOUT']
]);
/**
 * What keeps Wrap Rows from being turned on or its length changed. Turning it off is never blocked, so a Session
 * that holds both can follow the refusal's "turn Wrap Rows off". Whether a record is cut is Python's.
 * @param {LinearWrapBlockerInputs} inputs
 * @returns {LinearWrapBlocker[]}
 */
export const linearWrapBlockers = (inputs) => [...new Set(BLOCKERS.filter(([input]) => inputs[input]).map(([, blocker]) => blocker))];

// The length's unit is view state (D-27): the form and the Session hold bp, and the Web never sends unit text,
// so these factors are the field's own, not a copy of Python's length parser (D-57). Each is 10 to its exponent.
const UNIT_EXPONENTS = Object.freeze({ bp: 0, kbp: 3, Mbp: 6 });
export const LINEAR_WRAP_UNITS = Object.freeze({ bp: 1, kbp: 1_000, Mbp: 1_000_000 });
/** @typedef {keyof typeof LINEAR_WRAP_UNITS} LinearWrapUnit */
/** @type {readonly LinearWrapUnit[]} */
const LARGEST_FIRST = Object.freeze(['Mbp', 'kbp', 'bp']);

/** @param {unknown} unit @returns {unit is LinearWrapUnit} */
export const isLinearWrapUnit = (unit) => typeof unit === 'string' && Object.hasOwn(LINEAR_WRAP_UNITS, unit);

/**
 * The unit a length shows in when it comes from outside the field (Load, Undo, Reset): the largest one whose
 * amount is a whole number; Auto and Off keep `current`.
 * @param {unknown} value
 * @param {LinearWrapUnit} current
 * @returns {LinearWrapUnit}
 */
export const linearWrapUnitOf = (value, current) => (typeof value === 'number' && Number.isSafeInteger(value) && value > 0
  ? /** @type {LinearWrapUnit} */ (LARGEST_FIRST.find((unit) => value % LINEAR_WRAP_UNITS[unit] === 0))
  : current);

/**
 * The amount the length field shows: empty for Auto (D-28) and Off.
 * @param {unknown} value
 * @param {LinearWrapUnit} unit
 */
export const linearWrapAmountText = (value, unit) => (typeof value === 'number' ? String(value / LINEAR_WRAP_UNITS[unit]) : '');

/**
 * The value a typed amount writes: empty is Auto (D-28); otherwise amount x unit in bp. A decimal amount is
 * shifted by the unit's exponent digit by digit, so it writes a whole number exactly when it names one (`1.1` kbp
 * is 1100 bp, `1.0000000001` kbp is not whole). The number field also gives a sign or an exponent (`+5`, `1e3`):
 * the browser reads it as a number, which is multiplied by the unit (`1e3` kbp writes 1,000,000 bp, accepted when
 * whole, although the CLI refuses `--wrap_length 1e3kb`). A bp value that is not whole is sent as it is, and Python's
 * validator refuses it, as the CLI refuses `--wrap_length 1.0000000001kb` (D-51, D-57).
 * @param {string} text
 * @param {LinearWrapUnit} unit
 * @returns {'auto' | number}
 */
export const linearWrapInputValue = (text, unit) => {
  const trimmed = String(text).trim();
  const amount = Number(trimmed);
  if (!trimmed || !Number.isFinite(amount)) return 'auto';
  const decimal = /^(\d*)(?:\.(\d*))?$/.exec(trimmed);
  if (!decimal) return amount * LINEAR_WRAP_UNITS[unit];
  const [, whole, fraction = ''] = decimal;
  const shift = UNIT_EXPONENTS[unit];
  const digits = `${whole}${fraction.padEnd(shift, '0')}`;
  // The digits past the unit's exponent are the bp's fraction.
  return Number(`${digits.slice(0, whole.length + shift) || '0'}.${digits.slice(whole.length + shift) || '0'}`);
};

/**
 * What the shown Result says about the drawing's Wrap Rows value. `off`: wrap is off. `rows`: the Result cut
 * records with this value over these records (its input summary's `linearWrap`). `one-row`: the Result asked for
 * this value over these records and cut nothing, so the axis is the longest record. `unknown`: no Result drew
 * this draft.
 * @typedef {{ kind: 'off' } | { kind: 'unknown' } | { kind: 'one-row' }
 *   | { kind: 'rows', rowLength: number, rowCount: number, drawnAxisLength: number }} ShownLinearWrap
 */
/**
 * A record the next Linear Generate draws, as its card asks for it: the bases of its source record and the region of
 * it the card states (1-based, both ends included, the smaller first), null without one. Python's input summary
 * names each drawn record by the same two facts (`RecordSummary`), so the Web never predicts a cropped length.
 * @typedef {RecordSummary} LinearDraftRecord
 */
/**
 * @typedef {object} ShownLinearWrapInputs
 * @property {unknown} value The drawing's `form.linear_wrap_length`.
 * @property {LinearWrapFormValue | undefined} requested The shown Result's committed request's wrap length (read
 *   with `wrapFormValue`); undefined when no Result of this drawing is shown.
 * @property {Pick<InputSummary, 'records' | 'linearWrap'> | undefined} summary The shown Result's input
 *   summary; undefined when the render sent none for it (a loaded Result).
 * @property {readonly LinearDraftRecord[]} records The records the next Generate draws, [] while unknown.
 */
/** @param {LinearDraftRecord['region']} a @param {LinearDraftRecord['region']} b */
const sameRegion = (a, b) => (a === null || b === null ? a === b : a.start === b.start && a.end === b.end);
/**
 * Reads Python's summary only; the Auto rule, the cut test and the drawn axis stay Python's.
 * @param {ShownLinearWrapInputs} inputs
 * @returns {ShownLinearWrap}
 */
export const shownLinearWrap = ({ value, requested, summary, records }) => {
  if (value === null) return { kind: 'off' };
  if (requested === undefined || requested !== value || summary === undefined) return { kind: 'unknown' };
  // Rows and one row alike describe the next Generate only over the records the Result drew: the same source
  // lengths, each cropped to the same region (review 3 F2, D-53).
  const sameRecords = summary.records.length === records.length
    && summary.records.every(({ sourceLength, region }, index) => sourceLength === records[index].sourceLength
      && sameRegion(region, records[index].region));
  const wrap = summary.linearWrap;
  if (!sameRecords || (wrap !== null && wrap.requested !== value)) return { kind: 'unknown' };
  return wrap === null
    ? { kind: 'one-row' }
    : { kind: 'rows', rowLength: wrap.rowLength, rowCount: wrap.rowCount, drawnAxisLength: wrap.drawnAxisLength };
};

/**
 * The exported Result's wrap as the PNG/PDF size check reads it (gbdraw/web_support/export_limit.py `_WrapPayload`):
 * the committed request's wrap leaf as stored, and each record the Result drew, by its source length and region.
 * @typedef {{ requested: unknown, records: LinearDraftRecord[] }} ExportLinearWrap
 */
/**
 * What the size check of a PNG or PDF export reads about the exported Result's wrap (D-52, D-55, D-70): the
 * committed request's wrap leaf as stored, and the records the Result drew, from its input summary, or, for a Result
 * without one (a loaded Session, a Gallery Session before Generate), the records the drawing's discovery reads,
 * awaited as Generate awaits it. Python reads the leaf with its config validator (`"5kb"` and `"Auto"` included),
 * crops each record with the drawing's rule, and decides whether a record is cut and what would fit.
 * @param {{ leaf: unknown, summary: Pick<InputSummary, 'records'> | undefined,
 *   discoverRecords: () => Promise<readonly LinearDraftRecord[]> }} inputs `leaf` as `requestLinearWrapLeaf` reads it,
 *   undefined when no Linear Result of the drawing is shown; `discoverRecords` the records the drawing reads once its
 *   discovery has run, [] while unknown.
 * @returns {Promise<ExportLinearWrap | null>} null with wrap off or the records unknown.
 */
export const exportLinearWrap = async ({ leaf, summary, discoverRecords }) => {
  if (leaf === undefined || leaf === null) return null;
  const records = summary === undefined ? await discoverRecords() : summary.records;
  return records.length === 0 ? null
    : { requested: leaf, records: records.map(({ sourceLength, region }) => ({ sourceLength, region })) };
};
