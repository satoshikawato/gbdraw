// @ts-check
// The feature density tracks of the shown drawing (`adv.feature_density_tracks`) as its Feature density section edits
// them: Add, Duplicate, Move, Remove, and the fields of a track. A button or a field the user changes is one History
// step (app/history-inputs.js); a Remove that asks runs its chosen removal as its own step. A track is found by its
// id, never by its place or caption. A track's rows in the drawing's stack follow it: Add and Duplicate put one row
// for the new track into a stack that draws (Use custom stack on) or that a Session or the user saved, and leave a
// stack no Session saved, which Use custom stack builds again with every track (as PD-OI-058 does for Depth rows);
// Remove takes the track's rows out, asking first when other rows start at or extend to them (PD-OI-104). Python
// owns the track's values (R7): the fields write the request's shape, and Generate refuses a wrong value in Python's
// words.
/** @import { CanonicalFeatureDensityTrack } from '../contracts/request.generated.js' */
/** @import { DrawingState } from '../state.js' */
/** @import { AskingStackPorts, SwitchedTrackRow } from './track-slot-edits.js' */
/** @import { FeatureDensityTrackField } from '../services/feature-density-tracks.js' */
import { circularFeatureDensitySlot, circularStacksNoSessionSaved } from '../services/circular-track-slot-model.js';
import { linearFeatureDensitySlot, linearStacksNoSessionSaved } from '../services/linear-track-slot-model.js';
import {
  FEATURE_DENSITY_TRACK_CHOICES, copiedFeatureDensityTrack, featureDensitySlotId, insertFeatureDensityRow, newFeatureDensityTrack, removeFeatureDensityRows,
  writeFeatureDensityTrackField
} from '../services/feature-density-tracks.js';
import { isUnsavedDefaultStack, nextNumberedId } from '../services/track-slot-validation.js';
import { msg } from '../utils/i18n.js';
import { askTrackSourceRemoval } from './track-slot-edits.js';

/**
 * The words the editor shows for Python's metric and style choices (English, marked; shown with `tm`). A test checks
 * that they name exactly those choices.
 * @type {Readonly<Record<string, string>>}
 */
export const FEATURE_DENSITY_CHOICE_WORDS = Object.freeze({
  count: msg('Count'), coverage: msg('Coverage (%)'), heatmap: msg('Heatmap'), histogram: msg('Histogram')
});

/**
 * A stack row as this owner reads and writes it.
 * @typedef {{ id: string, renderer: string, enabled?: boolean, side?: string | null, params?: Record<string, unknown> }} DensityStackRow
 */

/**
 * @typedef {object} FeatureDensityTrackEditorPorts
 * @property {() => DrawingState} shownDrawing
 * @property {() => unknown} busy The session's busy answer; a value stops an edit.
 * @property {(slotsKey: 'circular_track_slots' | 'linear_track_slots') => AskingStackPorts<SwitchedTrackRow>} askingStackPorts
 * @property {(drawing: DrawingState) => void} normalizeStack The mode's stack normalizer, after a row change.
 */

/**
 * The drawing's tracks, as the request writes them (`adv` is untyped).
 * @param {DrawingState} drawing @returns {CanonicalFeatureDensityTrack[]}
 */
const tracksOf = (drawing) => drawing.adv.feature_density_tracks;

/** @param {DrawingState} drawing */
const stackOf = (drawing) => (drawing.mode === 'circular'
  ? { slotsKey: /** @type {const} */ ('circular_track_slots'), axisKey: /** @type {const} */ ('circular_track_slots_axis_index'),
      enabled: drawing.adv.circular_track_slots_enabled, makeRow: circularFeatureDensitySlot,
      unsaved: () => circularStacksNoSessionSaved(drawing.adv.nt, drawing.form.track_type) }
  : { slotsKey: /** @type {const} */ ('linear_track_slots'), axisKey: /** @type {const} */ ('linear_track_slots_axis_index'),
      enabled: drawing.adv.linear_track_slots_enabled, makeRow: linearFeatureDensitySlot,
      unsaved: () => linearStacksNoSessionSaved(drawing.adv.nt, drawing.form.linear_track_layout) });

/**
 * @param {FeatureDensityTrackEditorPorts} ports
 */
export const createFeatureDensityTrackEditor = ({ shownDrawing, busy, askingStackPorts, normalizeStack }) => {
  /** @param {DrawingState} drawing @param {string} trackId */
  const trackOf = (drawing, trackId) => tracksOf(drawing).find(({ id }) => id === trackId) ?? null;

  /**
   * Puts `track` into the drawing at `at` and its row into the stack that should draw it.
   * @param {DrawingState} drawing @param {CanonicalFeatureDensityTrack} track @param {number} at
   */
  const putTrack = (drawing, track, at) => {
    tracksOf(drawing).splice(at, 0, track);
    const stack = stackOf(drawing);
    const rows = /** @type {DensityStackRow[]} */ (drawing.adv[stack.slotsKey]);
    if (!stack.enabled && isUnsavedDefaultStack(rows, stack.unsaved())) return;
    const row = /** @type {DensityStackRow} */ (stack.makeRow(track.id));
    row.id = nextNumberedId(featureDensitySlotId(track.id), rows);
    drawing.adv[stack.axisKey] = insertFeatureDensityRow(rows, drawing.adv[stack.axisKey], row);
    normalizeStack(drawing);
  };

  /** @param {DrawingState} drawing @param {string} trackId */
  const takeTrack = (drawing, trackId) => {
    const tracks = tracksOf(drawing);
    const at = tracks.findIndex(({ id }) => id === trackId);
    if (at < 0) return;
    tracks.splice(at, 1);
    const stack = stackOf(drawing);
    drawing.adv[stack.axisKey] = removeFeatureDensityRows(
      /** @type {DensityStackRow[]} */ (drawing.adv[stack.slotsKey]), drawing.adv[stack.axisKey], trackId
    );
    normalizeStack(drawing);
  };

  /** Adds a track that counts the first feature type the drawing draws (its Feature Types list; every type when empty). */
  const addFeatureDensityTrack = () => {
    const stop = busy();
    if (stop) return stop;
    const drawing = shownDrawing();
    const tracks = tracksOf(drawing);
    const featureType = String(drawing.adv.features[0] ?? '') || FEATURE_DENSITY_TRACK_CHOICES.everyFeatureType;
    return putTrack(drawing, newFeatureDensityTrack(tracks, featureType), tracks.length);
  };

  /** @param {string} trackId */
  const duplicateFeatureDensityTrack = (trackId) => {
    const stop = busy();
    if (stop) return stop;
    const drawing = shownDrawing();
    const track = trackOf(drawing, trackId);
    if (!track) return undefined;
    const tracks = tracksOf(drawing);
    return putTrack(drawing, copiedFeatureDensityTrack(track, tracks), tracks.indexOf(track) + 1);
  };

  /**
   * Moves a track one place up (-1) or down (1): the simple stack draws the tracks in this order.
   * @param {string} trackId @param {-1 | 1} offset
   */
  const moveFeatureDensityTrack = (trackId, offset) => {
    const stop = busy();
    if (stop) return stop;
    const tracks = tracksOf(shownDrawing());
    const from = tracks.findIndex(({ id }) => id === trackId);
    const to = from + offset;
    if (from < 0 || to < 0 || to >= tracks.length) return undefined;
    tracks.splice(to, 0, ...tracks.splice(from, 1));
    return undefined;
  };

  /**
   * Removes a track and its rows; asks first when other rows start at or extend to its rows (PD-OI-104).
   * @param {string} trackId
   */
  const removeFeatureDensityTrack = (trackId) => {
    const stop = busy();
    if (stop) return stop;
    const drawing = shownDrawing();
    if (!trackOf(drawing, trackId)) return undefined;
    const stack = stackOf(drawing);
    const rows = /** @type {DensityStackRow[]} */ (drawing.adv[stack.slotsKey]);
    const leaving = new Set(rows
      .filter((row) => row.renderer === 'feature_density' && String(row.params?.track ?? '') === trackId)
      .map(({ id }) => id));
    return askTrackSourceRemoval(askingStackPorts(stack.slotsKey), {
      kind: 'feature-density-track',
      track: trackId,
      loss: { leaving, stayingOn: new Set(), emptySeries: () => [] },
      sources: () => tracksOf(drawing),
      label: msg('Remove feature density track'),
      apply: (turnOffReaching) => {
        turnOffReaching();
        takeTrack(drawing, trackId);
      }
    }) ?? takeTrack(drawing, trackId);
  };

  /**
   * Writes one field of a track as the request holds it (`writeFeatureDensityTrackField`).
   * @param {string} trackId @param {FeatureDensityTrackField} field @param {unknown} value
   */
  const setFeatureDensityTrackField = (trackId, field, value) => {
    const stop = busy();
    if (stop) return stop;
    const track = trackOf(shownDrawing(), trackId);
    return track ? writeFeatureDensityTrackField(track, field, value) : false;
  };

  return {
    // The template's choices: Python's, in its order, with their words.
    featureDensityTrackChoices: FEATURE_DENSITY_TRACK_CHOICES,
    featureDensityChoiceWords: FEATURE_DENSITY_CHOICE_WORDS,
    addFeatureDensityTrack,
    duplicateFeatureDensityTrack,
    moveFeatureDensityTrack,
    removeFeatureDensityTrack,
    setFeatureDensityTrackField
  };
};
