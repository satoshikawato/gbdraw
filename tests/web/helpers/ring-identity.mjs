// The identity invariant of a Circular drawing's comparison rings (#1777 round 8), as the ring tests and the ring
// identity walk check it after every ring transition:
// (a) each series entry's `sourceKey`, `fileName` and `sourceIndex` name one and the same source of the drawing's
//     current ring sources (`circularRingSources`), and no two entries name one source;
// (b) each managed stack row's `series_key` is the key the one matcher (`conservationSourcesOfSeries`) gives its
//     entry over those sources: with the stack on, the managed rows are one per drawn ring, in no other key; with the
//     stack off, a row that still carries a key carries a drawn ring's key (the rows follow the rings again when the
//     stack turns on).
import {
  circularRingSources, conservationSourceDescriptors, isManagedConservationSlot, orderedConservationSources
} from '../../../gbdraw/web/js/services/conservation-series.js';

/**
 * Every way the drawing's rings break the invariant, as readable lines; none when it holds.
 * @param {import('../../../gbdraw/web/js/state.js').DrawingState} drawing
 * @returns {string[]}
 */
export const ringIdentityDefects = (drawing) => {
  /** @type {string[]} */
  const defects = [];
  const sources = circularRingSources(drawing);
  const descriptors = conservationSourceDescriptors(sources);
  /** @type {Set<string>} */
  const named = new Set();
  drawing.circularConservation.series.forEach((entry, index) => {
    const source = descriptors.find((descriptor) => descriptor.sourceKey === entry.sourceKey);
    const shown = JSON.stringify([entry.sourceKey, entry.fileName, entry.sourceIndex]);
    if (!source) {
      defects.push(`(a) entry ${index} ${shown}: its key names no source`);
    } else if (source.fileName !== entry.fileName || source.placeIndex !== entry.sourceIndex) {
      defects.push(`(a) entry ${index} ${shown}: its key names ${JSON.stringify([source.sourceKey, source.fileName, source.placeIndex])}`);
    } else if (named.has(source.sourceKey)) {
      defects.push(`(a) entry ${index} ${shown}: another entry names that source`);
    }
    named.add(entry.sourceKey);
  });
  if (drawing.circularConservation.enabled !== true) return defects;
  const drawn = orderedConservationSources(sources, drawing.circularConservation).map((source) => source.sourceKey);
  const rowKeys = (Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [])
    .filter((slot) => isManagedConservationSlot(slot))
    .map((slot) => String(slot.params?.series_key || '').trim());
  const sorted = (keys) => JSON.stringify([...keys].sort());
  if (drawing.adv.circular_track_slots_enabled === true) {
    if (sorted(rowKeys) !== sorted(drawn)) defects.push(`(b) row keys ${sorted(rowKeys)} are not the drawn rings' ${sorted(drawn)}`);
  } else {
    const carried = rowKeys.filter(Boolean);
    if (carried.some((key) => !drawn.includes(key)) || new Set(carried).size !== carried.length) {
      defects.push(`(b) with the stack off, row keys ${sorted(carried)} are not distinct keys of the drawn rings ${sorted(drawn)}`);
    }
  }
  return defects;
};
