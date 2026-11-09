// @ts-check
// One circle: a Circular drawing whose records are drawn one after another on one circle
// (`layout.recordLayout` `one_circle`, `gbdraw circular --one_circle`). This module owns the Circular record
// layout of a draft (separate figures, the Multi-Record Canvas grid, or one circle: `circularRecordLayoutOf`),
// the choice between them (`chooseRecordLayout`), and the draft's one-circle settings
// (`adv.circular_record_layout`, `one_circle_gap_degrees`, `one_circle_min_length`, `one_circle_record_band`,
// `one_circle_left`): their fresh values and limits, which are Python's (`CircularMultiRecordOptions`, through
// the generated request contract), what a request's layout stores for them, and what a stored layout reads
// back as. A drawing draws one circle when its Multi-Record Canvas switch is on (all records in one figure) and
// its record layout is `one_circle`.
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { ONE_CIRCLE_PROFILE, readGcContentMode } from '../mode-profiles.js';
import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { projectOptionalNumber } from '../utils/optional-positive-number.js';
import { isPlainObject } from '../utils/contract-shape.js';
/** @import { ObjectShape } from '../utils/contract-shape.js' */
/** @import { CanonicalCircularLayout } from '../contracts/request.generated.js' */
/** @import { GcContentMode } from '../mode-profiles.js' */

/** @typedef {'grid' | 'one_circle'} CircularRecordLayout */
/** @typedef {'separate' | 'grid' | 'one_circle'} CircularRecordLayoutChoice */
/** @typedef {Exclude<CircularRecordLayoutChoice, 'one_circle'>} LeftRecordLayout */
/**
 * What one circle left when it was chosen, which leaving it restores: the layout, and the two values one circle
 * draws as the command does (GC content mode, Hide GC Skew).
 * @typedef {object} OneCircleLeft
 * @property {LeftRecordLayout} layout
 * @property {GcContentMode} gc_content_mode
 * @property {boolean} suppress_skew
 */
/**
 * The one-circle settings of a Circular draft.
 * @typedef {object} OneCircleAdv
 * @property {CircularRecordLayout} circular_record_layout
 * @property {number} one_circle_gap_degrees Degrees between two records.
 * @property {number} one_circle_min_length The shortest record drawn, in bp; shorter ones are counted, not drawn.
 * @property {boolean} one_circle_record_band Whether a band outside the tracks names each record.
 * @property {OneCircleLeft | null} one_circle_left What one circle left; null unless one circle is drawn.
 */
/**
 * The fields of a Circular draft the record layout reads.
 * @typedef {{ form: { multi_record_canvas: boolean }, adv: Pick<OneCircleAdv, 'circular_record_layout'> }} CircularLayoutDraft
 */
/**
 * The fields of a Circular draft the layout choice writes.
 * @typedef {{
 *   form: { multi_record_canvas: boolean, suppress_skew: boolean },
 *   adv: Pick<OneCircleAdv, 'circular_record_layout' | 'one_circle_left'> & { gc_content_mode: string }
 * }} CircularLayoutChoiceDraft
 */
/**
 * The one-circle fields of a request layout, all present.
 * @typedef {Required<Pick<CanonicalCircularLayout, 'recordLayout' | 'oneCircleGapDegrees' | 'oneCircleMinLength'
 *   | 'oneCircleRecordBand'>>} OneCircleLayout
 */
/** @typedef {{ kind: string, values?: readonly string[], minimum?: number, maximum?: number }} LayoutFieldShape */

const LAYOUT_SHAPE = /** @type {ObjectShape} */ (CANONICAL_REQUEST_DESCRIPTION.types.CanonicalCircularLayout);
const LAYOUT_DEFAULTS = /** @type {OneCircleLayout} */ (/** @type {unknown} */ (LAYOUT_SHAPE.defaults));
const LAYOUT_FIELDS = /** @type {Readonly<Record<keyof OneCircleLayout, LayoutFieldShape>>} */ (
  /** @type {unknown} */ (LAYOUT_SHAPE.fields));

/** The record layouts a Circular draft may hold: Python's (`CircularMultiRecordOptions.record_layout`). */
export const CIRCULAR_RECORD_LAYOUTS = /** @type {readonly CircularRecordLayout[]} */ (LAYOUT_FIELDS.recordLayout.values);

/** The ranges of the two numbers, as Python checks them: the inputs' limits and the Session readers'. */
export const ONE_CIRCLE_LIMITS = Object.freeze({
  gapDegrees: Object.freeze({
    minimum: /** @type {number} */ (LAYOUT_FIELDS.oneCircleGapDegrees.minimum),
    maximum: /** @type {number} */ (LAYOUT_FIELDS.oneCircleGapDegrees.maximum)
  }),
  minLength: Object.freeze({ minimum: /** @type {number} */ (LAYOUT_FIELDS.oneCircleMinLength.minimum) })
});

/** @param {number} value */
export const oneCircleGapInRange = (value) => (
  value >= ONE_CIRCLE_LIMITS.gapDegrees.minimum && value <= ONE_CIRCLE_LIMITS.gapDegrees.maximum
);
/** @param {number} value */
export const oneCircleMinLengthInRange = (value) => Number.isInteger(value) && value >= ONE_CIRCLE_LIMITS.minLength.minimum;

/** A fresh draft's one-circle settings. @type {Readonly<OneCircleAdv>} */
export const ONE_CIRCLE_ADV_DEFAULTS = Object.freeze({
  circular_record_layout: LAYOUT_DEFAULTS.recordLayout,
  one_circle_gap_degrees: LAYOUT_DEFAULTS.oneCircleGapDegrees,
  one_circle_min_length: LAYOUT_DEFAULTS.oneCircleMinLength,
  one_circle_record_band: LAYOUT_DEFAULTS.oneCircleRecordBand,
  one_circle_left: null
});

/**
 * The record layout a Circular draft draws: separate figures (Multi-Record Canvas off), the Multi-Record Canvas
 * grid, or one circle. Every reader that names the layout reads it here.
 * @param {CircularLayoutDraft} drawing
 * @returns {CircularRecordLayoutChoice}
 */
export const circularRecordLayoutOf = (drawing) => {
  if (!drawing.form.multi_record_canvas) return 'separate';
  return drawing.adv.circular_record_layout === 'one_circle' ? 'one_circle' : 'grid';
};

/**
 * Whether a Circular draft draws its records on one circle.
 * @param {CircularLayoutDraft} drawing
 */
export const drawsOneCircle = (drawing) => circularRecordLayoutOf(drawing) === 'one_circle';

/** The layout a fresh drawing draws, which a one-circle Session without a left layout returns to. */
const FRESH_LAYOUT = /** @type {LeftRecordLayout} */ (WEB_UX_PROFILE.circular.gridByDefault ? 'grid' : 'separate');

/**
 * Chooses a Circular record layout (PD-OI-019) and returns the Hide GC Skew value the caller applies through the
 * track editor, or null when it stays. Choosing one circle keeps what it leaves (`one_circle_left`) and draws GC
 * content and GC skew as the command does (`ONE_CIRCLE_PROFILE`, Python's `ONE_CIRCLE_CONFIG_OVERRIDES`); a
 * value the user then changes stays theirs while one circle is drawn. Leaving one circle (`'left'`: back to the
 * layout it left, or another layout) puts back the layout's GC content mode and Hide GC Skew. Any other choice
 * writes only the layout.
 * @param {CircularLayoutChoiceDraft} drawing
 * @param {CircularRecordLayoutChoice | 'left'} choice
 * @returns {boolean | null}
 */
export const chooseRecordLayout = (drawing, choice) => {
  const leaving = circularRecordLayoutOf(drawing);
  const left = drawing.adv.one_circle_left;
  /** @param {CircularRecordLayoutChoice} layout */
  const writeLayout = (layout) => {
    drawing.form.multi_record_canvas = layout !== 'separate';
    drawing.adv.circular_record_layout = layout === 'one_circle' ? 'one_circle' : 'grid';
  };
  if (choice === 'one_circle') {
    if (leaving === 'one_circle') return null;
    drawing.adv.one_circle_left = {
      layout: leaving, gc_content_mode: readGcContentMode(drawing.adv.gc_content_mode).value,
      suppress_skew: Boolean(drawing.form.suppress_skew)
    };
    writeLayout('one_circle');
    drawing.adv.gc_content_mode = ONE_CIRCLE_PROFILE.gcContentMode;
    return !ONE_CIRCLE_PROFILE.showSkew;
  }
  if (leaving !== 'one_circle') {
    if (choice !== 'left') writeLayout(choice);
    return null;
  }
  drawing.adv.one_circle_left = null;
  writeLayout(choice === 'left' ? (left?.layout ?? FRESH_LAYOUT) : choice);
  if (!left) return null;
  drawing.adv.gc_content_mode = left.gc_content_mode;
  return left.suppress_skew;
};

/**
 * What a one-circle request's Load reads one circle as having left: the fresh layout and, for the GC content
 * mode and Hide GC Skew, what the request draws, unless it is the one-circle profile's, which reads as unset
 * (Python's rule: an unset value is the profile's), so leaving one circle draws the fresh value.
 * @param {{ gc_content_mode: string, suppress_skew: boolean }} drawn What the request draws.
 * @param {{ gc_content_mode: string, suppress_skew: boolean }} fresh A fresh drawing's values.
 * @returns {OneCircleLeft}
 */
export const oneCircleLeftOfRequest = (drawn, fresh) => ({
  layout: FRESH_LAYOUT,
  gc_content_mode: readGcContentMode(drawn.gc_content_mode === ONE_CIRCLE_PROFILE.gcContentMode ? fresh.gc_content_mode
    : drawn.gc_content_mode).value,
  suppress_skew: drawn.suppress_skew === !ONE_CIRCLE_PROFILE.showSkew ? fresh.suppress_skew : drawn.suppress_skew
});

/**
 * A stored `one_circle_left`: the record it names, else null (a lost value when something was stored). Its GC
 * content mode is read as the Advanced setting reads one (`readGcContentMode`), so leaving one circle never
 * writes back a mode the setting would not hold.
 * @param {unknown} value
 * @returns {{ value: OneCircleLeft | null, lost: boolean }}
 */
export const readOneCircleLeft = (value) => {
  if (value === null || value === undefined) return { value: null, lost: false };
  const record = isPlainObject(value) ? /** @type {Record<string, unknown>} */ (value) : null;
  const gcContentMode = readGcContentMode(record?.gc_content_mode);
  if (!record || !['separate', 'grid'].includes(String(record.layout)) || typeof record.gc_content_mode !== 'string'
    || gcContentMode.lost || typeof record.suppress_skew !== 'boolean') return { value: null, lost: true };
  return {
    value: {
      layout: /** @type {LeftRecordLayout} */ (record.layout),
      gc_content_mode: gcContentMode.value,
      suppress_skew: record.suppress_skew
    },
    lost: false
  };
};

/**
 * The one-circle fields a draft's request layout stores: none for the grid, as Python writes none. A blank
 * number is the default; any other number is sent literally (the Python layer checks its range).
 * @param {OneCircleAdv} adv
 * @returns {OneCircleLayout | {}}
 */
export const oneCircleRequestLayout = (adv) => (adv.circular_record_layout !== 'one_circle' ? {} : {
  recordLayout: 'one_circle',
  oneCircleGapDegrees: projectOptionalNumber(adv.one_circle_gap_degrees, { field: 'one_circle_gap_degrees' })
    ?? LAYOUT_DEFAULTS.oneCircleGapDegrees,
  oneCircleMinLength: projectOptionalNumber(adv.one_circle_min_length, { field: 'one_circle_min_length' })
    ?? LAYOUT_DEFAULTS.oneCircleMinLength,
  oneCircleRecordBand: adv.one_circle_record_band !== false
});

/**
 * A stored Circular layout's one-circle settings as the draft holds them; an absent field is the default.
 * @param {Partial<CanonicalCircularLayout> | null} layout
 * @param {OneCircleLeft | null} left What one circle left (`oneCircleLeftOfRequest`), for a one-circle layout.
 * @returns {OneCircleAdv}
 */
export const oneCircleAdvOfLayout = (layout, left) => {
  const recordLayout = layout?.recordLayout ?? LAYOUT_DEFAULTS.recordLayout;
  return {
    circular_record_layout: recordLayout,
    one_circle_gap_degrees: layout?.oneCircleGapDegrees ?? LAYOUT_DEFAULTS.oneCircleGapDegrees,
    one_circle_min_length: layout?.oneCircleMinLength ?? LAYOUT_DEFAULTS.oneCircleMinLength,
    one_circle_record_band: layout?.oneCircleRecordBand ?? LAYOUT_DEFAULTS.oneCircleRecordBand,
    one_circle_left: recordLayout === 'one_circle' ? left : null
  };
};
