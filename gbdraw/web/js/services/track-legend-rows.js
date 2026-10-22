// @ts-check
import { normalizeCircularTrackSlot, normalizeNt, parseCircularTrackSlotSpec } from './circular-track-slot-model.js';
import { parseLinearTrackSlotSpecs } from './linear-track-slot-model.js';
import { defaultDepthSeriesLabel } from './depth-track-state.js';
import { normalizeOptionalText } from './track-slot-display.js';
import { samePaint } from '../utils/color-utils.js';

/** @import { PythonLegendKey, PythonLegendRow } from './legend-svg.js' */

// OV-484: the content, skew and Depth tracks a request draws, and the Legend
// rows Python drew for them in a Result, each with the track slot that draws it
// and the color that fills it. The live compile and the Legend list also read
// here which rows leave with a track the draft hides (`hiddenTrackLegendKeys`,
// OV-1198, OV-1215), by the predicate the track hide reads (`drawnTrackHidden`). The Legend panel, the live compile and the feature
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
// them; the Depth rows: gbdraw/analysis/depth_tracks.py
// `sync_depth_track_legend_entries`, and the Depth slots of the Circular sync);
// tests/fixtures/legend_track_slot_rows.json holds Python's answers.

/**
 * What this module reads of a track slot the slot parsers return
 * (`parseCircularTrackSlotSpec`, `normalizeCircularTrackSlot`, `parseLinearTrackSlotSpecs`).
 * @typedef {object} TrackRowSlot
 * @property {string} id
 * @property {string} renderer
 * @property {boolean} enabled
 * @property {{ nt?: string, legend_label?: unknown, label?: unknown, positive_color?: string, negative_color?: string,
 *   track_index?: unknown }} params
 *   `nt` absent: the request's dinucleotide. `legend_label`, else `label` where the slot has no
 *   `legend_label` key: the rows' label (blank or null: Python's). The colors are resolved, the
 *   `high_color` and `low_color` aliases included; absent: the palette's. `track_index`: the Depth
 *   series a Depth slot draws (absent: the first).
 */

/**
 * A content, skew or Depth track a request draws.
 * @typedef {object} DrawnTrack
 * @property {string | null} slotId The track slot that draws it; null: a default track, which no explicit slot
 *   draws (Python tags its group with a preset slot id of its own), or a Depth series without Depth slots.
 * @property {'content' | 'skew' | 'depth'} kind A Depth track has one row, in its series' own fill.
 * @property {string} label The label of its Legend rows.
 * @property {string | null} positiveColor A skew slot's own `positive_color`; null: the palette's.
 * @property {string | null} negativeColor A skew slot's own `negative_color`; null: the palette's.
 * @property {true} [saved] The track is read from a Session's saved settings (`loadedTracks`), which the
 *   Result need not have been drawn with: it has a row in its own color only where Python drew the row in it.
 */

/**
 * The track settings a Session 27 to 30 (gbdraw 0.13.0 and earlier) was saved with. Such a Session has
 * no request, and the Web shows its saved preview, so they stand in for the request that drew the
 * Result. They are the saved draft: a track setting changed after the last Generate did not draw the
 * Result (a slot's own color is checked against the drawn row; where the saved dinucleotide is not the
 * drawn one, the rows have no track until Generate).
 * @typedef {object} LoadedTrackSettings
 * @property {'circular' | 'linear'} mode
 * @property {string} dinucleotide
 * @property {readonly TrackRowSlot[] | null} slots The mode's Custom Track Slots as loaded; null (the slots
 *   are off): the default tracks.
 */

/** @typedef {'gc_content' | 'skew_high' | 'skew_low'} TrackPaletteKey The palette key of a track's fill. */

/**
 * The Depth series of a request: each series' label (Python's
 * `_default_depth_track_labels` for one without), by its index.
 * @typedef {readonly string[]} DepthSeriesLabels
 */

/**
 * A Legend row Python drew for a track.
 * @typedef {object} TrackLegendRow
 * @property {PythonLegendKey} key
 * @property {string | null} slotId The track slot that draws the row; null: a default track.
 * @property {DrawnTrack['kind']} kind
 * @property {TrackPaletteKey | null} paletteKey The palette key that fills the row unless the slot has its own
 *   color; null: a Depth row, filled by its series.
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
 * @returns {{ mode: 'circular' | 'linear', defaultNt: string, entries: ReadonlyArray<string | Record<string, unknown>> | null,
 *   depth: DepthSeriesLabels } | null}
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
  const depthTracks = diagramOptions.depthTracks ?? [];
  if (!Array.isArray(depthTracks) || !depthTracks.every(isRecord)) return null;
  const depth = depthTracks.map(({ label }, index) => (typeof label === 'string' && label
    ? label : defaultDepthSeriesLabel(index, depthTracks.length)));
  return { mode, defaultNt: normalizeNt(dinucleotide), entries, depth };
};

/** @type {Readonly<Record<DrawnTrack['kind'], number>>} The order of a Linear request's rows by kind. */
const LINEAR_KIND_ORDER = Object.freeze({ content: 0, depth: 1, skew: 2 });

/**
 * The content, skew and Depth tracks drawn from track slots, in Python's
 * order: Circular in slot order, Linear content rows, then Depth rows, then
 * skew rows. A Depth slot draws the series of its `track_index`, if the
 * request has it. Without slots, each mode draws each Depth series under its
 * label (Linear through `default_linear_track_slots`).
 * @param {'circular' | 'linear'} mode
 * @param {string} defaultNt The dinucleotide of a slot that names none.
 * @param {readonly TrackRowSlot[] | null} slots null: the default tracks.
 * @param {DepthSeriesLabels} [depth] The request's Depth series.
 * @returns {readonly DrawnTrack[]}
 */
const slotTracks = (mode, defaultNt, slots, depth = []) => {
  /** @param {DrawnTrack['kind']} kind @param {TrackRowSlot | null} slot @param {string} [fallback] @returns {DrawnTrack} */
  const track = (kind, slot, fallback = `${normalizeNt(slot?.params?.nt, defaultNt)} ${kind}`) => {
    const params = slot?.params ?? {};
    return {
      slotId: slot?.id ?? null,
      kind,
      label: normalizeOptionalText(Object.hasOwn(params, 'legend_label') ? params.legend_label : params.label) || fallback,
      positiveColor: kind === 'skew' ? normalizeOptionalText(params.positive_color) : null,
      negativeColor: kind === 'skew' ? normalizeOptionalText(params.negative_color) : null
    };
  };
  if (!slots) {
    return [track('content', null), track('skew', null), ...depth.map((label) => track('depth', null, label))];
  }
  /** @param {TrackRowSlot} slot @returns {DrawnTrack | null} */
  const slotTrack = (slot) => {
    if (slot.renderer !== 'depth') return TRACK_KIND[slot.renderer] ? track(TRACK_KIND[slot.renderer], slot) : null;
    const index = Number(slot.params.track_index ?? 0);
    return Number.isInteger(index) && index >= 0 && index < depth.length ? track('depth', slot, depth[index]) : null;
  };
  // Each enabled content, skew and Depth slot draws its own rows, in both modes (OV-846).
  const tracks = slots.filter((slot) => slot.enabled).map(slotTrack).filter((each) => each !== null);
  return mode === 'linear' ? tracks.sort((left, right) => LINEAR_KIND_ORDER[left.kind] - LINEAR_KIND_ORDER[right.kind]) : tracks;
};
/**
 * The content, skew and Depth tracks the request that drew a Result draws. None for a
 * value that is no request (a Result loaded without one is read with
 * `loadedTracks`) and for slot entries the slot parsers reject (as Session
 * load would): no track is named, so no row and no skew track takes a color here.
 * @param {unknown} request
 * @returns {readonly DrawnTrack[]}
 */
export const drawnTracks = (request) => {
  const view = trackRequestView(request);
  if (!view) return [];
  const { mode, defaultNt, entries, depth } = view;
  try {
    return slotTracks(mode, defaultNt, entries && (mode === 'linear'
      ? parseLinearTrackSlotSpecs(entries, defaultNt)
      : entries.map((entry, index) => (typeof entry === 'string'
          ? parseCircularTrackSlotSpec(entry, index, defaultNt) : normalizeCircularTrackSlot(entry, index, defaultNt)))), depth);
  } catch {
    return [];
  }
};

/**
 * The content and skew tracks of a Result that was loaded without a request, from the track
 * settings its Session was saved with (`DrawnTrack.saved`). No settings name no track.
 * @param {LoadedTrackSettings | null} settings
 * @returns {readonly DrawnTrack[]}
 */
export const loadedTracks = (settings) => (
  settings ? slotTracks(settings.mode, normalizeNt(settings.dinucleotide), settings.slots).map((track) => ({ ...track, saved: true })) : []
);

/**
 * The palette key and the own color of one side of a track.
 * @param {DrawnTrack} track
 * @param {'positive' | 'negative'} side
 * @returns {Pick<TrackLegendRow, 'paletteKey' | 'ownColor'>}
 */
const trackSideFill = (track, side) => {
  if (track.kind === 'depth') return { paletteKey: null, ownColor: null };
  if (track.kind === 'content') return { paletteKey: 'gc_content', ownColor: null };
  return side === 'positive'
    ? { paletteKey: 'skew_high', ownColor: track.positiveColor }
    : { paletteKey: 'skew_low', ownColor: track.negativeColor };
};

/**
 * The Legend rows Python drew for the tracks in a Result, by Python's key.
 *
 * Python gives a track one key under its label when its two fills were equal
 * at the Generate, else "<label> (+)" and "<label> (-)", and numbers a key
 * taken before it ("<key> (2)"), a feature row's included. A track whose row
 * has the name and the swatch of an earlier row draws no row and still takes
 * its key (Python's shared keys), so the tracks of a key take the numbered keys
 * after the one a row of another source holds, whether or not they are drawn
 * (tests/fixtures/legend_track_slot_rows.json, `shared`). A track claims only
 * a key Python drew (D-15-65): `CDS (2)` of a track labelled like the feature
 * row `CDS`, never `CDS` itself. A skew slot with both colors of its own says
 * which form it has; any other track takes the form whose keys the Result has
 * left over. A key is given to its track only where the Result has at most
 * one key more of that label than tracks that take one: a row the tracks
 * cannot be told apart for (tracks whose forms the Result leaves open) is no
 * track's row here and keeps the color Python drew.
 * @param {readonly DrawnTrack[]} tracks `drawnTracks` of the request that drew the Result.
 * @param {ReadonlyMap<PythonLegendKey, PythonLegendRow>} pythonRows The Legend rows Python drew in the Result.
 * @param {ReadonlySet<PythonLegendKey>} [sharedLegendKeys] Python's keys of the Result that no row draws, which
 *   hold the shared keys (its Legend row facts' `suppressed`); none: a Result without facts.
 * @returns {TrackLegendRows}
 */
export const trackLegendRows = (tracks, pythonRows, sharedLegendKeys = new Set()) => {
  /** @type {Set<string>} */
  const drawnKeys = new Set(pythonRows.keys());
  /** @param {string} key */
  const known = (key) => drawnKeys.has(key) || sharedLegendKeys.has(/** @type {PythonLegendKey} */ (key));
  /** @param {string} preferred @param {number} index @returns {PythonLegendKey} The key of the n-th source of one preferred key. */
  const numberedKey = (preferred, index) => /** @type {PythonLegendKey} */ (index ? `${preferred} (${index + 1})` : preferred);
  /** @param {string} preferred @param {number} count @returns {number} How many of the first `count` keys of one preferred key Python took. */
  const knownAmong = (preferred, count) => Array.from({ length: count }, (_, index) => numberedKey(preferred, index)).filter(known).length;
  /** @type {Map<string, DrawnTrack[]>} */
  const tracksByLabel = new Map();
  tracks.forEach((track) => tracksByLabel.set(track.label, [...(tracksByLabel.get(track.label) || []), track]));
  /** @type {Map<PythonLegendKey, TrackLegendRow>} */
  const rows = new Map();
  /** @param {string} preferred @param {DrawnTrack[]} drawers @param {'positive' | 'negative'} side */
  const give = (preferred, drawers, side) => {
    // 1: a row of another source (a feature row) holds the unnumbered key.
    const offset = knownAmong(preferred, drawers.length + 1) - drawers.length;
    if (offset !== 0 && offset !== 1) return;
    drawers.forEach((drawer, index) => {
      const key = numberedKey(preferred, offset + index);
      if (!drawnKeys.has(key)) return;
      const row = { key, slotId: drawer.slotId, kind: drawer.kind, ...trackSideFill(drawer, side) };
      if (drawer.saved && row.ownColor && !samePaint(row.ownColor, pythonRows.get(key)?.color)) return;
      rows.set(key, row);
    });
  };
  tracksByLabel.forEach((group, label) => {
    // One key more than the tracks could take tells a row of another source.
    const oneRowKeys = knownAmong(label, group.length + 1);
    const twoRowKeys = Math.max(knownAmong(`${label} (+)`, group.length + 1), knownAmong(`${label} (-)`, group.length + 1));
    // A Depth track and a skew slot with both colors of its own say their form.
    const stated = group.map((track) => (track.kind === 'depth' ? true
      : track.positiveColor && track.negativeColor ? track.positiveColor === track.negativeColor : null));
    const leftOver = (/** @type {number} */ keyCount, /** @type {boolean} */ oneRow) => keyCount - stated.filter((form) => form === oneRow).length;
    const open = leftOver(oneRowKeys, true) <= 0 ? false
      : leftOver(twoRowKeys, false) <= 0 ? true : null;
    const oneRow = stated.map((form) => form ?? open);
    if (oneRow.includes(null)) return;
    const twoRowTracks = group.filter((_, index) => !oneRow[index]);
    give(label, group.filter((_, index) => oneRow[index]), 'positive');
    give(`${label} (+)`, twoRowTracks, 'positive');
    give(`${label} (-)`, twoRowTracks, 'negative');
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
export const trackRowColor = (row, palette) => row.ownColor || (row.paletteKey && palette[row.paletteKey]) || '';

/**
 * The draft toggles of the drawn tracks: Hide GC Content and Hide GC Skew in Circular, Show GC Content and
 * Show GC Skew in Linear, and Show Depth.
 * @typedef {{ suppress_gc?: unknown, suppress_skew?: unknown, show_gc?: unknown, show_skew?: unknown, show_depth?: unknown }} TrackVisibilityForm
 */

/**
 * The draft's Custom Track Slots of each mode. While a mode's stack is on, its enabled slots name the tracks
 * the request draws instead of the toggles that the stack disables.
 * @typedef {{ circular_track_slots_enabled?: unknown, linear_track_slots_enabled?: unknown,
 *   circular_track_slots?: unknown, linear_track_slots?: unknown }} TrackStackDraft
 */

/** @typedef {{ form: TrackVisibilityForm, adv: TrackStackDraft }} TrackVisibilityDrawing */

/**
 * Whether the draft's request draws Depth: an enabled Depth slot while the mode's custom stack is on, else
 * Show Depth. The request's serializer and the Generate input check read it too.
 * @param {'circular' | 'linear'} mode
 * @param {TrackVisibilityDrawing} drawing
 * @returns {boolean}
 */
export const draftDrawsDepth = (mode, { form, adv }) => {
  const linear = mode === 'linear';
  if (!(linear ? adv.linear_track_slots_enabled : adv.circular_track_slots_enabled)) return Boolean(form.show_depth);
  const slots = linear ? adv.linear_track_slots : adv.circular_track_slots;
  return Array.isArray(slots) && slots.some((slot) => isRecord(slot) && slot.enabled !== false && slot.renderer === 'depth');
};

/**
 * Whether the draft hides the drawn tracks of a kind (`app/svg-styles.js` `applyTrackVisibility`): a track
 * the draft's request does not draw. Depth follows `draftDrawsDepth`. A Linear custom stack draws its
 * content and skew slots whatever Show GC Content and Show GC Skew say (Python's
 * `default_linear_track_slots` is only the fallback); Circular reads Hide GC Content and Hide GC Skew,
 * which its slot editor keeps in step with the slots.
 * @param {DrawnTrack['kind']} kind
 * @param {'circular' | 'linear'} mode
 * @param {TrackVisibilityDrawing} drawing
 * @returns {boolean}
 */
export const drawnTrackHidden = (kind, mode, drawing) => {
  if (kind === 'depth') return !draftDrawsDepth(mode, drawing);
  const { form, adv } = drawing;
  if (mode === 'circular') return Boolean(kind === 'content' ? form.suppress_gc : form.suppress_skew);
  return !adv.linear_track_slots_enabled && !(kind === 'content' ? form.show_gc : form.show_skew);
};

/**
 * The keys of a Result's Legend rows whose drawn track the draft hides: each
 * leaves the live Result and the Legend list with its track and returns with
 * it (OV-1198, OV-1215, PD-OI-104). Generate draws none of them.
 * @param {TrackLegendRows} rows `trackLegendRows` of the displayed Result.
 * @param {'circular' | 'linear'} mode
 * @param {TrackVisibilityDrawing} drawing
 * @returns {PythonLegendKey[]}
 */
export const hiddenTrackLegendKeys = (rows, mode, drawing) => (
  [...rows.values()].filter((row) => drawnTrackHidden(row.kind, mode, drawing)).map((row) => row.key)
);

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
