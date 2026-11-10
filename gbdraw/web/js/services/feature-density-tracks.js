// @ts-check
// The feature density tracks of a drawing's draft (`adv.feature_density_tracks`): the tracks a Session's request
// gave the drawing (`diagramOptions.featureDensityTracks`). Load reads them into the draft
// (`projectCanonicalSessionRequest`), and Generate and Save write them back into the request
// (`projectCanonicalRenderInput`), so a Session the CLI or the Python API wrote keeps its tracks. Python owns a
// track (gbdraw/analysis/feature_density.py) and its request shape (`CanonicalFeatureDensityTrack`, carried by
// tools/generate_web_contracts.py); the Web checks a stored track against that shape and never derives one. The
// copies are JSON copies: a drawing's draft is a Vue proxy, which `structuredClone` refuses.
/** @import { CanonicalFeatureDensityTrack } from '../contracts/request.generated.js' */
/** @import { FieldReading } from '../utils/field-reading.js' */
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { shapeViolations } from '../utils/contract-shape.js';
import { isUnsetValue } from '../utils/field-reading.js';
import { cloneJsonData } from './json-clone.js';
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';

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
export const insertFeatureDensitySlots = (rows, axisIndex, trackIds, makeSlot) => {
  if (trackIds.length === 0 || rows.some(({ renderer }) => renderer === 'feature_density')) return axisIndex;
  const at = rows.findIndex(({ renderer }) => renderer === 'features') + 1;
  rows.splice(at, 0, ...trackIds.map(makeSlot));
  return axisIndex !== null && at < axisIndex ? axisIndex + trackIds.length : axisIndex;
};
