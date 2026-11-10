// @ts-check
// The feature density tracks of a drawing's draft (`adv.feature_density_tracks`): the tracks a Session's request
// gave the drawing (`diagramOptions.featureDensityTracks`). Load reads them into the draft
// (`projectCanonicalSessionRequest`), and Generate and Save write them back into the request
// (`projectCanonicalRenderInput`), so a Session the CLI or the Python API wrote keeps its tracks. Python owns a
// track (gbdraw/analysis/feature_density.py) and its request shape (`CanonicalFeatureDensityTrack`, carried by
// tools/generate_web_contracts.py); the Web checks a stored track against that shape and never derives one. The
// copies are JSON copies: a drawing's draft is a Vue proxy, which `structuredClone` refuses.
/** @import { CanonicalFeatureDensityTrack, CanonicalFeatureQuery } from '../contracts/request.generated.js' */
/** @import { FieldReading } from '../utils/field-reading.js' */
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { shapeViolations } from '../utils/contract-shape.js';
import { isUnsetValue } from '../utils/field-reading.js';
import { cloneJsonData } from './json-clone.js';
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';
import { classifyOptionalNumber } from '../utils/optional-positive-number.js';
import { nextCopyId, nextNumberedId } from './track-slot-validation.js';

const TRACK_SHAPE = Object.freeze({ kind: /** @type {const} */ ('ref'), name: 'CanonicalFeatureDensityTrack' });

/**
 * Whether a value is one track as the request writes it (every field, no other key).
 * @param {unknown} value
 * @returns {value is CanonicalFeatureDensityTrack}
 */
const isFeatureDensityTrack = (value) => shapeViolations(value, TRACK_SHAPE, CANONICAL_REQUEST_DESCRIPTION.types, '').length === 0;

/**
 * A stored `adv.feature_density_tracks` read into the draft: the tracks of the request's shape. Anything but a list
 * reads as none and is lost; a list entry of another shape is dropped and counted.
 * @param {unknown} value
 * @returns {FieldReading<CanonicalFeatureDensityTrack[]>}
 */
export const readFeatureDensityTracks = (value) => {
  if (!Array.isArray(value)) return { value: [], lost: !isUnsetValue(value) };
  const tracks = value.filter(isFeatureDensityTrack);
  const dropped = value.length - tracks.length;
  return { value: cloneJsonData(tracks), lost: dropped > 0, ...(dropped > 0 ? { dropped } : {}) };
};

/**
 * The draft's tracks of a loaded request's `diagramOptions.featureDensityTracks`, which Load has checked against the
 * request shape; a request without tracks (absent or null) gives none.
 * @param {readonly CanonicalFeatureDensityTrack[] | null | undefined} requested
 * @returns {CanonicalFeatureDensityTrack[]}
 */
export const projectFeatureDensityTracks = (requested) => cloneJsonData([...(requested ?? [])]);

/**
 * The request's `featureDensityTracks` of a draft: the key only when the drawing has tracks, so the request of a
 * drawing without them is the request it was before the tracks existed. A drawing's draft always holds the list
 * (`createDefaultAdv`).
 * @param {{ feature_density_tracks: readonly CanonicalFeatureDensityTrack[] }} adv
 * @returns {{ featureDensityTracks?: CanonicalFeatureDensityTrack[] }}
 */
export const requestFeatureDensityTracks = ({ feature_density_tracks: tracks }) => (
  tracks.length > 0 ? { featureDensityTracks: cloneJsonData([...tracks]) } : {}
);

/**
 * The id of the slot that draws track `trackId` in a stack the Web builds from the simple settings, as Python names
 * it in a stack the reader did not write (gbdraw/tracks/parsing.py `feature_density_slot_id`).
 * @param {string} trackId
 * @returns {string}
 */
export const featureDensitySlotId = (trackId) => `${MODE_PROFILE_DATA.featureDensitySlotIdPrefix}${trackId}`;

/**
 * A `feature_density` row's params as both stacks read them, edited in place: the track id trimmed, as Python reads
 * it (gbdraw/tracks/parsing.py `feature_density_slot_track_id`).
 * @param {Record<string, unknown>} params
 */
export const normalizeFeatureDensitySlotParams = (params) => {
  if (Object.hasOwn(params, 'track')) params.track = String(params.track ?? '').trim();
};

/**
 * The ids of a draft's tracks, in order.
 * @param {{ feature_density_tracks: readonly CanonicalFeatureDensityTrack[] }} adv
 * @returns {string[]}
 */
export const featureDensityTrackIds = ({ feature_density_tracks: tracks }) => tracks.map(({ id }) => id);

/**
 * The track a new `feature_density` row draws: the first of the drawing's tracks no enabled row draws yet, else the
 * first track; '' when the drawing has none.
 * @param {readonly CanonicalFeatureDensityTrack[]} tracks
 * @param {readonly { renderer?: unknown, enabled?: unknown, params?: { track?: unknown } }[]} slots
 * @returns {string}
 */
export const nextFeatureDensitySlotTrack = (tracks, slots) => {
  const drawn = new Set(slots
    .filter((slot) => slot.enabled !== false && slot.renderer === 'feature_density')
    .map((slot) => String(slot.params?.track ?? '')));
  return (tracks.find(({ id }) => !drawn.has(id)) ?? tracks[0])?.id ?? '';
};

/**
 * Puts `row` after the last feature density row of a stack, else after its Features row, else at the top: where
 * Python puts a track's row in a stack the reader did not write (`insertFeatureDensitySlots`), so the rows stay in
 * track order. Edits `rows` in place and returns the Axis index that keeps every other row on its side of the Axis.
 * @template {{ renderer?: unknown }} T
 * @param {T[]} rows
 * @param {number | null} axisIndex
 * @param {T} row
 * @returns {number | null}
 */
export const insertFeatureDensityRow = (rows, axisIndex, row) => {
  const last = rows.findLastIndex(({ renderer }) => renderer === 'feature_density');
  const at = (last >= 0 ? last : rows.findIndex(({ renderer }) => renderer === 'features')) + 1;
  rows.splice(at, 0, row);
  return axisIndex !== null && at < axisIndex ? axisIndex + 1 : axisIndex;
};

/**
 * Puts one row per track (`makeSlot`) after the Features row of a stack that draws no feature density track (at the
 * top without one), as Python does to a stack the reader did not write (gbdraw/tracks/parsing.py
 * `with_feature_density_slots`): the stacks the Web builds from the simple settings, and the stack Use custom stack
 * turns on, draw the tracks the simple stack drew. Edits `rows` in place and returns the Axis index that keeps every
 * other row on its side of the Axis.
 * @template {{ renderer?: unknown }} T
 * @param {T[]} rows
 * @param {number | null} axisIndex
 * @param {readonly string[]} trackIds
 * @param {(trackId: string) => T} makeSlot
 * @returns {number | null}
 */
export const insertFeatureDensitySlots = (rows, axisIndex, trackIds, makeSlot) => (
  rows.some(({ renderer }) => renderer === 'feature_density') ? axisIndex
    : trackIds.reduce((axis, trackId) => insertFeatureDensityRow(rows, axis, makeSlot(trackId)), axisIndex)
);

/**
 * Takes the rows that draw track `trackId` (on or off) out of a stack, in place, and returns the Axis index that
 * keeps every other row on its side of the Axis.
 * @template {{ renderer?: unknown, params?: { track?: unknown } }} T
 * @param {T[]} rows
 * @param {number | null} axisIndex
 * @param {string} trackId
 * @returns {number | null}
 */
export const removeFeatureDensityRows = (rows, axisIndex, trackId) => {
  const drawsTrack = (/** @type {T} */ row) => row.renderer === 'feature_density' && String(row.params?.track ?? '') === trackId;
  const before = axisIndex === null ? 0 : rows.slice(0, axisIndex).filter(drawsTrack).length;
  rows.splice(0, rows.length, ...rows.filter((row) => !drawsTrack(row)));
  return axisIndex === null ? null : axisIndex - before;
};

/**
 * The choices of a track's metric and style, and what a track the Web adds starts with (Python's).
 * @typedef {object} FeatureDensityTrackChoices
 * @property {readonly CanonicalFeatureDensityTrack['metric'][]} metrics
 * @property {readonly CanonicalFeatureDensityTrack['style'][]} styles
 * @property {CanonicalFeatureDensityTrack['metric']} defaultMetric
 * @property {CanonicalFeatureDensityTrack['style']} defaultStyle
 * @property {string} everyFeatureType The feature type that counts every type (`*`).
 */
/** @type {FeatureDensityTrackChoices} */
export const FEATURE_DENSITY_TRACK_CHOICES = /** @type {FeatureDensityTrackChoices} */ (MODE_PROFILE_DATA.featureDensityTrack);

/**
 * A track the Web adds: the first unused id of `track`, `track_2`, ... (its stack row is `density_track`); a query of `featureType` alone; Python's
 * metric and style; Auto window and step (the records' GC content window); the track color (`color` null) and the
 * Legend name of its query (`caption` null).
 * @param {readonly CanonicalFeatureDensityTrack[]} tracks
 * @param {string} featureType
 * @returns {CanonicalFeatureDensityTrack}
 */
export const newFeatureDensityTrack = (tracks, featureType) => ({
  id: nextNumberedId('track', tracks),
  query: { featureType, qualifier: null, value: null, minLength: null, maxLength: null },
  metric: FEATURE_DENSITY_TRACK_CHOICES.defaultMetric,
  style: FEATURE_DENSITY_TRACK_CHOICES.defaultStyle,
  window: null,
  step: null,
  color: null,
  caption: null
});

/**
 * A copy of `track` under the next copy id (`nextCopyId`: `genes` -> `genes_2`).
 * @param {CanonicalFeatureDensityTrack} track
 * @param {readonly CanonicalFeatureDensityTrack[]} tracks
 * @returns {CanonicalFeatureDensityTrack}
 */
export const copiedFeatureDensityTrack = (track, tracks) => ({ ...cloneJsonData(track), id: nextCopyId(track.id, tracks) });

/**
 * A field of a track the editor writes: one of the query's, or the track's own.
 * @typedef {keyof CanonicalFeatureQuery | 'metric' | 'style' | 'window' | 'step' | 'color' | 'caption'} FeatureDensityTrackField
 */
/** @type {ReadonlySet<FeatureDensityTrackField>} */
const QUERY_FIELDS = new Set(['featureType', 'qualifier', 'value', 'minLength', 'maxLength']);
/** @type {ReadonlySet<FeatureDensityTrackField>} */
const NUMBER_FIELDS = new Set(['minLength', 'maxLength', 'window', 'step']);

/**
 * Writes an edited field into `track` as the request holds it: a number field's blank as null and its number as that
 * number (a text that is no number changes nothing; Python checks the range, R7); a text field trimmed, its blank as
 * null; the feature type trimmed, its blank as every type; the metric and style only as one of Python's choices.
 * Returns whether the track changed.
 * @param {CanonicalFeatureDensityTrack} track
 * @param {FeatureDensityTrackField} field
 * @param {unknown} raw
 * @returns {boolean}
 */
export const writeFeatureDensityTrackField = (track, field, raw) => {
  /** @type {unknown} */
  let value;
  if (NUMBER_FIELDS.has(field)) {
    const number = classifyOptionalNumber(raw);
    if (number.status === 'invalid') return false;
    value = number.value;
  } else if (field === 'metric' || field === 'style') {
    const choice = String(raw);
    /** @type {readonly string[]} */
    const choices = field === 'metric' ? FEATURE_DENSITY_TRACK_CHOICES.metrics : FEATURE_DENSITY_TRACK_CHOICES.styles;
    if (!choices.includes(choice)) return false;
    value = choice;
  } else {
    const text = String(raw ?? '').trim();
    value = field === 'featureType' ? text || FEATURE_DENSITY_TRACK_CHOICES.everyFeatureType : text || null;
  }
  /** @type {Record<string, unknown>} */
  const owner = QUERY_FIELDS.has(field) ? track.query : track;
  if (owner[field] === value) return false;
  owner[field] = value;
  return true;
};
