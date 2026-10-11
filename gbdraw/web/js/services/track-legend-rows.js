// @ts-check
import { uniqueLegendKey } from './specific-color-rules.js';

/** @import { PythonLegendKey } from './legend-svg.js' */

// OV-484: the Legend rows Python draws for the content and skew tracks of a
// request, each with the track slot that draws it and the color that fills
// it. Python allocates the keys (gbdraw/legend/table.py `prepare_legend_table`
// without explicit slots; `_sync_legend_table_for_circular_slots` and
// `_sync_legend_table_for_linear_slots` in gbdraw/diagrams/*/assemble.py with
// them); tests/fixtures/legend_track_slot_rows.json holds Python's answers.
// The Legend panel, the live compile and the palette's track painter read a
// track row's color here, so a row is never colored from the shape of its
// caption.

/**
 * What the track rows read of a track slot of the canonical request
 * (`buildCircularTrackSlotPayload`, `buildLinearTrackSlotPayload`).
 * @typedef {object} TrackRowSlot
 * @property {string} id
 * @property {string} renderer
 * @property {boolean} [enabled] Absent: enabled. The request writer keeps enabled slots only.
 * @property {{ nt?: string, legend_label?: string, positive_color?: string, negative_color?: string }} [params]
 *   Absent values: the request's dinucleotide, Python's label, the palette's skew colors.
 */

/**
 * What the track rows read of the request that drew a Result
 * (`CanonicalRenderRequest`, session-request.js).
 * @typedef {object} TrackRowRequest
 * @property {string} mode
 * @property {{
 *   dinucleotide?: string,
 *   tracks?: { circularTrackSlots?: readonly TrackRowSlot[] | null, linearTrackSlots?: readonly TrackRowSlot[] | null }
 * }} diagramOptions `dinucleotide` absent: GC. `tracks` or its mode's slots absent or null: the
 *   request has no explicit slots, and Python draws the default tracks.
 */

/** @typedef {'gc_content' | 'skew_high' | 'skew_low'} TrackPaletteKey The palette key of a track's fill. */

/**
 * A Legend row Python draws for a track.
 * @typedef {object} TrackLegendRow
 * @property {PythonLegendKey} key
 * @property {string | null} slotId The track slot that draws the row; null: a default track (no explicit slot draws it).
 * @property {TrackPaletteKey} paletteKey The palette key that fills the row unless the slot has its own color.
 * @property {string | null} ownColor The slot's own color (`positive_color`, `negative_color`), which wins
 *   over the palette; null: the palette's.
 */

/** @type {Readonly<Record<string, 'content' | 'skew'>>} */
const TRACK_KIND = Object.freeze({ dinucleotide_content: 'content', dinucleotide_skew: 'skew' });

/** @param {string | undefined} value @returns {string | null} */
const textOrNull = (value) => (value === undefined ? '' : String(value)).trim() || null;

// Python draws one row under a track's label when its two fills are equal
// (the content track always; a skew track of one color), else "<label> (+)"
// and "<label> (-)". The fills are those of the last Generate, which a live
// palette edit may make equal or distinct, so every form is listed: the rows
// of a Result are found by their keys.
/** @type {ReadonlyArray<readonly [suffix: string, side: 'positive' | 'negative']>} */
const ROW_FORMS = Object.freeze([['', 'positive'], [' (+)', 'positive'], [' (-)', 'negative']]);

/**
 * The explicit track slots of a request's mode that Python draws; null when it has none.
 * @param {TrackRowRequest} request
 * @returns {readonly TrackRowSlot[] | null}
 */
const drawnTrackSlots = ({ mode, diagramOptions }) => {
  const slots = mode === 'linear' ? diagramOptions.tracks?.linearTrackSlots : diagramOptions.tracks?.circularTrackSlots;
  return slots ? slots.filter((slot) => slot.enabled !== false) : null;
};

/**
 * The Legend rows of the content and skew tracks of the request that drew a
 * Result, by Python's key; none without a request.
 * @param {TrackRowRequest | null} request
 * @returns {Map<PythonLegendKey, TrackLegendRow>}
 */
export const trackLegendRows = (request) => {
  /** @type {Map<PythonLegendKey, TrackLegendRow>} */
  const rows = new Map();
  if (!request) return rows;
  const defaultNt = (textOrNull(request.diagramOptions.dinucleotide) || 'GC').toUpperCase();
  /** @type {Set<string>} */
  const reserved = new Set();
  /** @param {'content' | 'skew'} kind @param {TrackRowSlot | null} slot */
  const add = (kind, slot) => {
    const label = textOrNull(slot?.params?.legend_label) || `${(textOrNull(slot?.params?.nt) || defaultNt).toUpperCase()} ${kind}`;
    ROW_FORMS.forEach(([suffix, side]) => {
      const key = /** @type {PythonLegendKey} */ (uniqueLegendKey(reserved, `${label}${suffix}`));
      reserved.add(key);
      rows.set(key, kind === 'content'
        ? { key, slotId: slot?.id ?? null, paletteKey: 'gc_content', ownColor: null }
        : {
            key, slotId: slot?.id ?? null, paletteKey: side === 'positive' ? 'skew_high' : 'skew_low',
            ownColor: textOrNull(side === 'positive' ? slot?.params?.positive_color : slot?.params?.negative_color)
          });
    });
  };
  const slots = drawnTrackSlots(request);
  // Linear slots replace only the skew rows: the content row stays the default one.
  if (!slots || request.mode === 'linear') add('content', null);
  if (!slots) add('skew', null);
  (slots || []).forEach((slot) => {
    const kind = TRACK_KIND[slot.renderer];
    if (kind && !(kind === 'content' && request.mode === 'linear')) add(kind, slot);
  });
  return rows;
};

/**
 * The fill of a track's Legend row and of the part of the track it names: the
 * slot's own color, else the palette's ('' when the palette has none).
 * @param {Pick<TrackLegendRow, 'ownColor' | 'paletteKey'>} row
 * @param {Readonly<Record<string, string>>} palette
 * @returns {string}
 */
export const trackRowColor = (row, palette) => row.ownColor || palette[row.paletteKey] || '';

/**
 * The fill of one side of a skew track of a Result, by the slot that draws it
 * (`data-gbdraw-slot-id`; null: the default skew track).
 * @param {Iterable<TrackLegendRow>} rows
 * @param {string | null} slotId
 * @param {'skew_high' | 'skew_low'} paletteKey
 * @param {Readonly<Record<string, string>>} palette
 * @returns {string}
 */
export const skewTrackColor = (rows, slotId, paletteKey, palette) => {
  for (const row of rows) {
    if (row.slotId === slotId && row.paletteKey === paletteKey) return trackRowColor(row, palette);
  }
  return trackRowColor({ ownColor: null, paletteKey }, palette);
};
