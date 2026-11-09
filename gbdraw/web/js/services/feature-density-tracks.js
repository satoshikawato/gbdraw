// @ts-check
// The feature density tracks of a drawing's draft (`adv.feature_density_tracks`): the tracks a Session's request
// gave the drawing (`diagramOptions.featureDensityTracks`). Load reads them into the draft
// (`projectCanonicalSessionRequest`), and Generate and Save write them back into the request
// (`projectCanonicalRenderInput`), so a Session the CLI or the Python API wrote keeps its tracks. Python owns a
// track (gbdraw/analysis/feature_density.py) and its request shape (`CanonicalFeatureDensityTrack`, carried by
// tools/generate_web_contracts.py); the Web checks a stored track against that shape and never derives one.
/** @import { CanonicalFeatureDensityTrack } from '../contracts/request.generated.js' */
/** @import { FieldReading } from '../utils/field-reading.js' */
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { shapeViolations } from '../utils/contract-shape.js';
import { isUnsetValue } from '../utils/field-reading.js';

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
  return { value: structuredClone(tracks), lost: dropped > 0, ...(dropped > 0 ? { dropped } : {}) };
};

/**
 * The draft's tracks of a loaded request's `diagramOptions.featureDensityTracks`, which Load has checked against the
 * request shape; a request without tracks (absent or null) gives none.
 * @param {readonly CanonicalFeatureDensityTrack[] | null | undefined} requested
 * @returns {CanonicalFeatureDensityTrack[]}
 */
export const projectFeatureDensityTracks = (requested) => structuredClone([...(requested ?? [])]);

/**
 * The request's `featureDensityTracks` of a draft: the key only when the drawing has tracks, so the request of a
 * drawing without them is the request it was before the tracks existed. A drawing's draft always holds the list
 * (`createDefaultAdv`); a draft object built without it, as the request tests build theirs, has no tracks.
 * @param {{ feature_density_tracks?: readonly CanonicalFeatureDensityTrack[] }} adv
 * @returns {{ featureDensityTracks?: CanonicalFeatureDensityTrack[] }}
 */
export const requestFeatureDensityTracks = ({ feature_density_tracks: tracks = [] }) => (
  tracks.length > 0 ? { featureDensityTracks: structuredClone([...tracks]) } : {}
);
