// @ts-check
import { normalizeCircularTrackSlot, normalizeNt, parseCircularTrackSlotSpec } from './circular-track-slot-model.js';
import { parseLinearTrackSlotSpecs } from './linear-track-slot-model.js';
import { normalizeOptionalText } from './track-slot-display.js';

/** @import { PythonLegendKey } from './legend-svg.js' */

// OV-484: the content and skew tracks a request draws, and the Legend rows
// Python drew for them in a Result, each with the track slot that draws it and
// the color that fills it. The Legend panel, the live compile and the feature
// popup read a track row's color here (`trackLegendRows`), and the palette's
// track painter a skew track's (`skewTrackColor`), so neither is read from the
// shape of a row's caption. A content track has no color of its own: every
// one takes the palette's gc_content, which the painter reads directly.
//
// The slots are read through the slot model's parser and normalizer
// (`parseCircularTrackSlotSpec`, `normalizeCircularTrackSlot`,
// `parseLinearTrackSlotSpecs`), so a spec string of a Python API request and a
// slot object give the same track.
// Python allocates the row keys (gbdraw/legend/table.py `prepare_legend_table`
// without explicit slots; `_sync_legend_table_for_circular_slots` and
// `_sync_legend_table_for_linear_slots` in gbdraw/diagrams/*/assemble.py with
// them); tests/fixtures/legend_track_slot_rows.json holds Python's answers.

/**
 * What this module reads of a track slot the slot parsers return
 * (`parseCircularTrackSlotSpec`, `normalizeCircularTrackSlot`, `parseLinearTrackSlotSpecs`).
 * @typedef {object} TrackRowSlot
 * @property {string} id
 * @property {string} renderer
 * @property {boolean} enabled
 * @property {{ nt?: string, legend_label?: unknown, label?: unknown, positive_color?: string, negative_color?: string }} params
 *   `nt` absent: the request's dinucleotide. `legend_label`, else `label` where the slot has no
 *   `legend_label` key: the rows' label (blank or null: Python's). The colors are resolved, the
 *   `high_color` and `low_color` aliases included; absent: the palette's.
 */

/**
 * A content or skew track a request draws.
 * @typedef {object} DrawnTrack
 * @property {string | null} slotId The track slot that draws it; null: a default track, which no explicit slot
 *   draws (Python tags its group with a preset slot id of its own).
 * @property {'content' | 'skew'} kind
 * @property {string} label The label of its Legend rows.
 * @property {string | null} positiveColor A skew slot's own `positive_color`; null: the palette's.
 * @property {string | null} negativeColor A skew slot's own `negative_color`; null: the palette's.
 */

/**
 * The track settings a Result was drawn with, where no request states them: the Web shows a
 * Session 27 to 30 (gbdraw 0.13.0 and earlier) from its saved settings and its saved preview.
 * @typedef {object} LoadedTrackSettings
 * @property {'circular' | 'linear'} mode
 * @property {string} dinucleotide
 * @property {readonly TrackRowSlot[] | null} slots The mode's Custom Track Slots as loaded; null (the slots
 *   are off): the default tracks.
 */

/** @typedef {'gc_content' | 'skew_high' | 'skew_low'} TrackPaletteKey The palette key of a track's fill. */

/**
 * A Legend row Python drew for a track.
 * @typedef {object} TrackLegendRow
 * @property {PythonLegendKey} key
 * @property {string | null} slotId The track slot that draws the row; null: a default track.
 * @property {TrackPaletteKey} paletteKey The palette key that fills the row unless the slot has its own color.
 * @property {string | null} ownColor The slot's own color, which wins over the palette; null: the palette's.
 */

/** @typedef {ReadonlyMap<PythonLegendKey, TrackLegendRow>} TrackLegendRows The track rows of a Result, by Python's key. */

/** @type {Readonly<Record<string, DrawnTrack['kind']>>} */
const TRACK_KIND = Object.freeze({ dinucleotide_content: 'content', dinucleotide_skew: 'skew' });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => typeof value === 'object' && value !== null && !Array.isArray(value);

/**
 * What the tracks read of the request that drew a Result, checked: its mode,
 * its dinucleotide and the slot entries of its mode (null: no explicit slots,
 * so Python draws the default tracks). Null for a value that is no such request.
 * @param {unknown} request
 * @returns {{ mode: 'circular' | 'linear', defaultNt: string, entries: ReadonlyArray<string | Record<string, unknown>> | null } | null}
 */
const trackRequestView = (request) => {
  if (!isRecord(request)) return null;
  const { mode, diagramOptions } = request;
  if ((mode !== 'circular' && mode !== 'linear') || !isRecord(diagramOptions)) return null;
  const dinucleotide = diagramOptions.dinucleotide ?? 'GC';
  const tracks = diagramOptions.tracks ?? null;
  if (typeof dinucleotide !== 'string' || (tracks !== null && !isRecord(tracks))) return null;
  const entries = tracks?.[mode === 'linear' ? 'linearTrackSlots' : 'circularTrackSlots'] ?? null;
  if (entries !== null && !(Array.isArray(entries) && entries.every((entry) => typeof entry === 'string' || isRecord(entry)))) return null;
  return { mode, defaultNt: normalizeNt(dinucleotide), entries };
};

/**
 * The content and skew tracks a mode draws from its track slots, in Python's order.
 * @param {'circular' | 'linear'} mode
 * @param {string} defaultNt The dinucleotide of a slot that names none.
 * @param {readonly TrackRowSlot[] | null} slots null: the default tracks.
 * @returns {readonly DrawnTrack[]}
 */
const slotTracks = (mode, defaultNt, slots) => {
  /** @param {DrawnTrack['kind']} kind @param {TrackRowSlot | null} slot @returns {DrawnTrack} */
  const track = (kind, slot) => {
    const params = slot?.params ?? {};
    return {
      slotId: slot?.id ?? null,
      kind,
      label: normalizeOptionalText(Object.hasOwn(params, 'legend_label') ? params.legend_label : params.label)
        || `${normalizeNt(params.nt, defaultNt)} ${kind}`,
      positiveColor: kind === 'skew' ? normalizeOptionalText(params.positive_color) : null,
      negativeColor: kind === 'skew' ? normalizeOptionalText(params.negative_color) : null
    };
  };
  if (!slots) return [track('content', null), track('skew', null)];
  // Linear slots replace only the skew rows: the content row stays the default one.
  const drawn = slots.filter((slot) => slot.enabled && TRACK_KIND[slot.renderer] && !(mode === 'linear' && TRACK_KIND[slot.renderer] === 'content'));
  return [...(mode === 'linear' ? [track('content', null)] : []), ...drawn.map((slot) => track(TRACK_KIND[slot.renderer], slot))];
};
/**
 * The content and skew tracks the request that drew a Result draws. None for a
 * value that is no request (a Result loaded without one is read with
 * `loadedTracks`) and for slot entries the slot parsers reject (as Session
 * load would): no track is named, so no row and no skew track takes a color here.
 * @param {unknown} request
 * @returns {readonly DrawnTrack[]}
 */
export const drawnTracks = (request) => {
  const view = trackRequestView(request);
  if (!view) return [];
  const { mode, defaultNt, entries } = view;
  try {
    return slotTracks(mode, defaultNt, entries && (mode === 'linear'
      ? parseLinearTrackSlotSpecs(entries, defaultNt)
      : entries.map((entry, index) => (typeof entry === 'string'
          ? parseCircularTrackSlotSpec(entry, index, defaultNt) : normalizeCircularTrackSlot(entry, index, defaultNt)))));
  } catch {
    return [];
  }
};

/**
 * The content and skew tracks of a Result that was loaded without a request, from the track
 * settings Load recorded for it. No settings name no track.
 * @param {LoadedTrackSettings | null} settings
 * @returns {readonly DrawnTrack[]}
 */
export const loadedTracks = (settings) => (
  settings ? slotTracks(settings.mode, normalizeNt(settings.dinucleotide), settings.slots) : []
);

/**
 * The palette key and the own color of one side of a track.
 * @param {DrawnTrack} track
 * @param {'positive' | 'negative'} side
 * @returns {Pick<TrackLegendRow, 'paletteKey' | 'ownColor'>}
 */
const trackSideFill = (track, side) => {
  if (track.kind === 'content') return { paletteKey: 'gc_content', ownColor: null };
  return side === 'positive'
    ? { paletteKey: 'skew_high', ownColor: track.positiveColor }
    : { paletteKey: 'skew_low', ownColor: track.negativeColor };
};

/**
 * The Legend rows Python drew for the tracks in a Result, by Python's key.
 *
 * Python draws one row under a track's label when its two fills were equal at
 * the Generate, else "<label> (+)" and "<label> (-)", and numbers a key the
 * table already holds ("<key> (2)"). A skew slot with both colors of its own
 * says which; any other track takes the form whose rows the Result has left
 * over. A key is given to a track only where the Result has exactly as many
 * rows of that key as tracks that draw one, in their order: a row the tracks
 * cannot be told apart for (a feature row of the same caption, tracks whose
 * forms the Result leaves open) is no track's row here and keeps the color
 * Python drew.
 * @param {readonly DrawnTrack[]} tracks `drawnTracks` of the request that drew the Result.
 * @param {Iterable<string>} pythonKeys The keys of the Legend rows Python drew in the Result.
 * @returns {TrackLegendRows}
 */
export const trackLegendRows = (tracks, pythonKeys) => {
  const drawnKeys = new Set(pythonKeys);
  /** @param {string} preferred @returns {PythonLegendKey[]} The drawn keys Python numbered from one preferred key. */
  const numberedKeys = (preferred) => {
    const keys = [];
    for (let key = preferred, number = 2; drawnKeys.has(key); key = `${preferred} (${number++})`) keys.push(/** @type {PythonLegendKey} */ (key));
    return keys;
  };
  /** @type {Map<string, DrawnTrack[]>} */
  const tracksByLabel = new Map();
  tracks.forEach((track) => tracksByLabel.set(track.label, [...(tracksByLabel.get(track.label) || []), track]));
  /** @type {Map<PythonLegendKey, TrackLegendRow>} */
  const rows = new Map();
  /** @param {PythonLegendKey[]} keys @param {DrawnTrack[]} drawers @param {'positive' | 'negative'} side */
  const give = (keys, drawers, side) => {
    if (keys.length !== drawers.length) return;
    keys.forEach((key, index) => rows.set(key, { key, slotId: drawers[index].slotId, ...trackSideFill(drawers[index], side) }));
  };
  tracksByLabel.forEach((group, label) => {
    const oneRowKeys = numberedKeys(label);
    const positiveKeys = numberedKeys(`${label} (+)`);
    const negativeKeys = numberedKeys(`${label} (-)`);
    const stated = group.map((track) => (track.positiveColor && track.negativeColor ? track.positiveColor === track.negativeColor : null));
    const leftOver = (/** @type {number} */ keyCount, /** @type {boolean} */ oneRow) => keyCount - stated.filter((form) => form === oneRow).length;
    const open = leftOver(oneRowKeys.length, true) <= 0 ? false
      : leftOver(Math.max(positiveKeys.length, negativeKeys.length), false) <= 0 ? true : null;
    const oneRow = stated.map((form) => form ?? open);
    if (oneRow.includes(null)) return;
    const twoRowTracks = group.filter((_, index) => !oneRow[index]);
    give(oneRowKeys, group.filter((_, index) => oneRow[index]), 'positive');
    give(positiveKeys, twoRowTracks, 'positive');
    give(negativeKeys, twoRowTracks, 'negative');
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
 * The fill of one side of a skew track of a Result, by the slot id of its
 * group (`data-gbdraw-slot-id`). The default skew track of a request without
 * explicit slots is the only skew track of its Result, so it answers for every
 * skew group, whatever preset id Python tagged it with. '' for a group the
 * request names no slot for, which keeps its drawn color.
 * @param {readonly DrawnTrack[]} tracks `drawnTracks` of the request that drew the Result.
 * @param {string | null} slotId
 * @param {'positive' | 'negative'} side
 * @param {Readonly<Record<string, string>>} palette
 * @returns {string}
 */
export const skewTrackColor = (tracks, slotId, side, palette) => {
  const track = tracks.find((candidate) => candidate.kind === 'skew' && (candidate.slotId === null || candidate.slotId === slotId));
  return track ? trackRowColor(trackSideFill(track, side), palette) : '';
};
