// The identity invariant of a Circular drawing's comparison rings (#1777 round 8, strengthened in review 7), as the
// ring tests and the ring identity walk check it after every ring transition:
// (a) each series entry's `sourceKey`, `fileName` and `sourceIndex` name one and the same source of the drawing's
//     current ring sources (`circularRingSources`), and no two entries name one source; in a LOSAT-cache replay, whose
//     rings are its tables, each table sits beside its ring's genome list (one genome list per table, as the next
//     Generate, a ring row and `keepRingsAtPlaces` read them), and the genome the next Generate searches for an entry
//     is the one at that entry's place;
// (b) each drawn ring has exactly one managed stack row, the one carrying the key the one matcher
//     (`conservationSourcesOfSeries`) gives its entry, and with the stack on no managed row carries another key; with
//     the stack off, a row that still carries a key carries a drawn ring's key (the rows follow the rings again when the
//     stack turns on). Given `rows` (each drawn ring's row id, in drawn order, from the caller's own model), the row of
//     each drawn ring is that row: a row that moved to another ring is a defect even when the keys stay one per ring.
import {
  circularRingSources, conservationSourceDescriptors, conservationSourcesOfSeries, isManagedConservationSlot, orderedConservationSources
} from '../../../gbdraw/web/js/services/conservation-series.js';
import { ringGenomeFileLists } from '../../../gbdraw/web/js/services/drawing/inputs.js';

/**
 * Every way the drawing's rings break the invariant, as readable lines; none when it holds.
 * @param {import('../../../gbdraw/web/js/state.js').DrawingState} drawing
 * @param {{ rows?: readonly (string | null)[] }} [expected] Each drawn ring's own row id, in drawn order (null: no
 *   row expected yet, any one row passes).
 * @returns {string[]}
 */
export const ringIdentityDefects = (drawing, { rows } = {}) => {
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
  const replay = String(drawing.circularConservation.source || '').trim().toLowerCase() !== 'upload'
    && drawing.files.c_conservation_blasts_source === 'losat-cache' && sources.length > 0;
  if (replay) {
    const genomes = ringGenomeFileLists(drawing.files.c_conservation_fastas);
    if (genomes.length !== sources.length) {
      defects.push(`(a) ${sources.length} tables beside ${genomes.length} genome lists`);
    }
    const searched = conservationSourcesOfSeries(genomes, drawing.circularConservation).matched;
    drawing.circularConservation.series.forEach((entry, index) => {
      const genome = searched[index];
      if (genome && genome.placeIndex !== entry.sourceIndex) {
        defects.push(`(a) entry ${index}: the next Generate searches the genome at place ${genome.placeIndex}, its table is at ${entry.sourceIndex}`);
      }
    });
  }
  if (drawing.circularConservation.enabled !== true) return defects;
  const drawn = orderedConservationSources(sources, drawing.circularConservation);
  const managed = (Array.isArray(drawing.adv.circular_track_slots) ? drawing.adv.circular_track_slots : [])
    .filter((slot) => isManagedConservationSlot(slot));
  /** @param {{ params?: { series_key?: unknown } }} slot */
  const keyOf = (slot) => String(slot.params?.series_key || '').trim();
  const drawnKeys = new Set(drawn.map((ring) => ring.sourceKey));
  if (drawing.adv.circular_track_slots_enabled === true) {
    drawn.forEach((ring, index) => {
      const own = managed.filter((slot) => keyOf(slot) === ring.sourceKey);
      if (own.length !== 1) {
        defects.push(`(b) drawn ring ${index} (${ring.label}) has ${own.length} rows`);
      } else if (rows && rows[index] !== null && rows[index] !== undefined && String(own[0].id) !== rows[index]) {
        defects.push(`(b) drawn ring ${index} (${ring.label}) has row ${own[0].id}, its own is ${rows[index]}`);
      }
    });
    managed.filter((slot) => !drawnKeys.has(keyOf(slot)))
      .forEach((slot) => defects.push(`(b) row ${slot.id} carries ${JSON.stringify(keyOf(slot))}, no drawn ring's key`));
  } else {
    const carried = managed.map(keyOf).filter(Boolean);
    if (carried.some((key) => !drawnKeys.has(key)) || new Set(carried).size !== carried.length) {
      defects.push(`(b) with the stack off, row keys ${JSON.stringify(carried)} are not distinct keys of the drawn rings`);
    }
  }
  return defects;
};
