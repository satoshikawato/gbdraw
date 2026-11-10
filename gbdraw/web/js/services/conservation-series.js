// @ts-check
import {
  CONSERVATION_COLOR_SCHEMES, CONSERVATION_DEFAULT_COLOR_SCHEME, CONSERVATION_MAX_COLOR_PATH, CONSERVATION_SAME_COLOR
} from '../utils/conservation-colors.generated.js';
import { normalizeOptionalHexColor } from '../utils/color-utils.js';
import { ringGenomeFileLists } from './drawing/inputs.js';
/** @import { InputFileRef } from './history-files.js' */
/** @import { ConservationColorScheme } from '../utils/conservation-colors.generated.js' */

/**
 * What a comparison source slot holds, as these helpers read it (its name,
 * size and date): an `InputFileRef` in the project's bindings, a File or
 * Session resource view while a Session is decoded.
 * @typedef {{ readonly name?: string, readonly size?: number, readonly lastModified?: number }} ConservationFile
 */

/**
 * One comparison ring: the source file it reads and how it is drawn. `color` is the ring's own colour (hand-set,
 * lowercase #rrggbb), or null when it takes the colour of the drawing's Ring colors scheme at `schemeIndex`
 * (`effectiveRingColor`, PD-OI-110).
 * @typedef {{
 *   sourceKey: string, fileName: string, sourceIndex: number, label: string, color: string | null,
 *   schemeIndex: number, losat_gencode: number
 * }} ConservationSeriesEntry
 */
/**
 * A series entry as a Session or an older draft holds it: every field is unvalidated.
 * @typedef {{
 *   sourceKey?: unknown, fileName?: unknown, sourceIndex?: unknown, label?: unknown, name?: unknown,
 *   color?: unknown, schemeIndex?: unknown, losat_gencode?: unknown
 * }} RawConservationSeriesEntry
 */
/**
 * The comparison rings of a drawing as these helpers read them: the series, the legacy names, and the Ring colors
 * scheme (`colorScheme`; unset or unknown reads as `same`).
 * @typedef {{ series?: RawConservationSeriesEntry[] | null, labels?: unknown, colorScheme?: unknown }} CircularConservationLike
 */
/**
 * A ring as the request, the managed stack row, and the panel read it: `color` is the colour it draws
 * (`effectiveRingColor`), or '' when it draws the configured colour (`same`), as Python reads ''.
 * @typedef {Omit<ConservationSeriesEntry, 'sourceKey' | 'fileName' | 'sourceIndex' | 'color' | 'schemeIndex'>
 *   & { color: string }} DrawnRingFields
 */
/**
 * One comparison source: `files` are what it reads, `file` the first of them, which names it (its file name, its
 * default label). A BLAST table is a source of one file; a ring's comparison genome is a source of its files, in the
 * listed order (a chromosome and its plasmids).
 * @template {ConservationFile} [F=ConservationFile]
 * `sourceIndex` counts the sources; `placeIndex` is the source's place in the given list, empty places counted.
 * @typedef {{ file: F, files: readonly F[], sourceKey: string, sourceIndex: number, placeIndex: number,
 *   fileName: string, defaultLabel: string }} ConservationSourceDescriptor
 */
/**
 * The places of a source list as the helpers read them: a file (a BLAST table), a ring's genome files, or nothing.
 * @template {ConservationFile} F
 * @typedef {F | readonly (F | null)[] | null | undefined} ConservationSourcePlace
 */
/**
 * A ring in drawn order: its source, the fields it draws, its place, and the colour fields of the series entry it
 * reads (`ownColor`, the entry's own colour or null, and `schemeIndex`), which a Generate's adoption gives back to the
 * ring it searched (the request holds only the drawn colour).
 * @template {ConservationFile} [F=ConservationFile]
 * @typedef {ConservationSourceDescriptor<F> & DrawnRingFields
 *   & { orderIndex: number, ownColor: string | null, schemeIndex: number }} OrderedConservationSource
 */

export const CONSERVATION_SLOT_MANAGER = 'circular_conservation';

/**
 * The Ring colors scheme a stored value names, else the default (`CONSERVATION_DEFAULT_COLOR_SCHEME`, Python's).
 * @param {unknown} value
 * @returns {ConservationColorScheme}
 */
export const readRingColorScheme = (value) => (
  typeof value === 'string' && Object.hasOwn(CONSERVATION_COLOR_SCHEMES, value)
    ? /** @type {ConservationColorScheme} */ (value) : CONSERVATION_DEFAULT_COLOR_SCHEME
);

/**
 * The one owner of the colour a comparison ring draws (PD-OI-110): its own colour (hand-set) on every scheme, else
 * the scheme's colour at its `schemeIndex`, repeating, else null under `same` (the configured colour). The panel,
 * the managed stack row, and the request read it through `orderedConservationSources`.
 * @param {RawConservationSeriesEntry | null | undefined} entry
 * @param {unknown} scheme
 * @returns {string | null}
 */
export const effectiveRingColor = (entry, scheme) => {
  const own = normalizeOptionalHexColor(entry?.color);
  if (own) return own;
  const cycle = CONSERVATION_COLOR_SCHEMES[readRingColorScheme(scheme)];
  if (cycle.length === 0) return null;
  return cycle[readSchemeIndex(entry?.schemeIndex, 0) % cycle.length];
};

/**
 * What the Same color colour reads of a drawing: the configuration values its Session or request stored that the
 * Web does not manage, by dotted config path (the drawing's preserved settings).
 * @typedef {{ readonly unmanagedConfigOverrides: Readonly<Record<string, unknown>> }} SameRingColorDrawing
 */
/**
 * The colour a swatch shows for a ring colour: for a ring that draws the configured colour ('' or null), the
 * drawing's `objects.conservation.max_color` (a preserved setting), else the packaged one, as Python draws it.
 * @param {string | null | undefined} color
 * @param {SameRingColorDrawing | null} drawing
 */
export const ringSwatchColor = (color, drawing) => color
  || normalizeOptionalHexColor(drawing?.unmanagedConfigOverrides[CONSERVATION_MAX_COLOR_PATH])
  || CONSERVATION_SAME_COLOR;

/**
 * A stored scheme index: an integer from 0, else `fallback` (the ring's place).
 * @param {unknown} value
 * @param {number} fallback
 */
export const readSchemeIndex = (value, fallback) => (
  typeof value === 'number' && Number.isInteger(value) && value >= 0 ? value : fallback
);

/**
 * The scheme index a new ring takes (D-05): the first that no ring of the series uses, so a ring added after a
 * removal takes the freed colour.
 * @param {readonly { schemeIndex?: unknown }[]} series
 */
const firstUnusedSchemeIndex = (series) => {
  const used = new Set(series.map((entry) => entry.schemeIndex));
  let index = 0;
  while (used.has(index)) index += 1;
  return index;
};

/**
 * A track slot as the drawn order reads it: a Web stack row or a request's slot.
 * @typedef {{ renderer?: unknown, enabled?: unknown, params?: { track_index?: unknown } | null }} RingSlot
 */
/**
 * Each series place's drawn position (outermost first), as Python counts it for the scheme colours
 * (`_drawn_ring_colors`): the enabled `sequence_conservation` slots in stack order, each naming its ring by
 * `params.track_index` (the 1-based series place); without a stack (`slots` null), the series order. A ring no slot
 * draws follows the drawn ones, in series order.
 * @param {number} count The number of series entries.
 * @param {readonly RingSlot[] | null} slots
 * @returns {number[]}
 */
export const ringDrawnPositions = (count, slots) => {
  /** @type {number[]} */
  const drawn = [];
  (slots || []).forEach((slot) => {
    if (slot?.enabled === false || slot?.renderer !== 'sequence_conservation') return;
    const place = Number(slot.params?.track_index ?? 0) - 1;
    if (Number.isInteger(place) && place >= 0 && place < count && !drawn.includes(place)) drawn.push(place);
  });
  if (slots === null) drawn.splice(0, drawn.length, ...Array.from({ length: count }, (_, place) => place));
  for (let place = 0; place < count; place += 1) if (!drawn.includes(place)) drawn.push(place);
  /** @type {number[]} */
  const positions = new Array(count).fill(0);
  drawn.forEach((place, position) => { positions[place] = position; });
  return positions;
};

/**
 * A ring's own colour (hand-set, lowercase #rrggbb), or null when it takes the scheme's colour.
 * @param {{ color?: unknown } | null | undefined} entry
 */
export const ringOwnColor = (entry) => normalizeOptionalHexColor(entry?.color);

/**
 * Re-assign the scheme colours to the current order (PD-OI-110): each ring without its own colour takes the scheme
 * index of its drawn position (`positions`, by series place); a hand-set colour stays.
 * @param {{ color?: unknown, schemeIndex?: unknown }[]} series
 * @param {readonly number[]} positions
 * @returns {boolean} Whether an index changed.
 */
export const reassignRingSchemeIndexes = (series, positions) => {
  let changed = false;
  series.forEach((entry, place) => {
    if (ringOwnColor(entry) || entry.schemeIndex === positions[place]) return;
    entry.schemeIndex = positions[place];
    changed = true;
  });
  return changed;
};

/**
 * What `reassignDrawnRingColors` reads of a drawing: its rings, their scheme, and its slot stack.
 * @typedef {{
 *   circularConservation: { series: { color?: unknown, schemeIndex?: unknown }[], colorScheme: unknown },
 *   adv: { circular_track_slots_enabled?: unknown, circular_track_slots?: readonly RingSlot[] }
 * }} RingColorDrawing
 */
/**
 * Re-assign a drawing's scheme colours to the order its rings are drawn in (`ringDrawnPositions` over its stack
 * while the custom stack is in use, else the series order): the button 「今の並びで色を振り直す」, and a reorder
 * while 「トラックを並び替えたら色を振り直す」 is on. Under Same color both are disabled, and disabled has no
 * effect (PD-OI-110 「同じ色のときはどちらも無効」): the indexes stay as they were.
 * @param {RingColorDrawing} drawing
 * @returns {boolean} Whether an index changed.
 */
export const reassignDrawnRingColors = (drawing) => {
  const { series, colorScheme } = drawing.circularConservation;
  if (readRingColorScheme(colorScheme) === 'same') return false;
  const slots = drawing.adv.circular_track_slots_enabled === true ? drawing.adv.circular_track_slots ?? [] : null;
  return reassignRingSchemeIndexes(series, ringDrawnPositions(series.length, slots));
};

/**
 * After a reorder (the series arrows or the custom stack), re-assign the scheme colours while
 * 「トラックを並び替えたら色を振り直す」 is on, in the move's own step (PD-OI-110).
 * @param {RingColorDrawing & { circularConservation: { reassignOnReorder?: unknown } }} drawing
 * @returns {boolean} Whether an index changed.
 */
export const reassignRingColorsOnReorder = (drawing) => (
  drawing.circularConservation.reassignOnReorder === true && reassignDrawnRingColors(drawing)
);

/**
 * @template T
 * @param {T | (T | null)[] | null | undefined} files a null place is dropped
 * @returns {T[]}
 */
export const normalizeFileList = (files) => (Array.isArray(files) ? files.filter(/** @returns {file is T} */ (file) => Boolean(file)) : (files ? [files] : []));

/**
 * What `circularRingSources` reads of a Circular drawing.
 * @template {ConservationFile} [F=InputFileRef]
 * @typedef {{
 *   files: {
 *     readonly c_conservation_blasts: readonly F[], readonly c_conservation_blasts_source: string | null,
 *     readonly c_conservation_fastas: readonly (readonly F[])[]
 *   },
 *   circularConservation: { readonly source?: unknown }
 * }} RingSourceDrawing
 */
/** @param {RingSourceDrawing<ConservationFile>} drawing */
const isUploadedRingDrawing = (drawing) => String(drawing.circularConservation.source || '').trim().toLowerCase() === 'upload';
/**
 * Whether a drawing's rings follow BLAST tables: uploaded ones, or those of a LOSAT-cache replay.
 * @param {RingSourceDrawing<ConservationFile>} drawing
 */
const ringsFollowTables = (drawing) => isUploadedRingDrawing(drawing) || (
  drawing.files.c_conservation_blasts_source === 'losat-cache' && normalizeFileList([...drawing.files.c_conservation_blasts]).length > 0
);
/**
 * The sources of a Circular drawing's rings, the one reading of which files
 * its ring series follow: the BLAST tables of uploaded rings and of a
 * LOSAT-cache replay, else each LOSAT ring's comparison genome files at its
 * place, an empty place included (a Session loaded without a ring's files
 * keeps its place, as Generate and Save read it), so a series entry's stored
 * place names the same ring for every reader.
 * @param {RingSourceDrawing} drawing
 * @returns {readonly (InputFileRef | readonly InputFileRef[])[]}
 */
export const circularRingSources = (drawing) => {
  if (ringsFollowTables(drawing)) return normalizeFileList([...drawing.files.c_conservation_blasts]);
  return ringGenomeFileLists(drawing.files.c_conservation_fastas);
};

/**
 * The files of a source place, in order: one for a file, a ring genome's files for a list.
 * @template {ConservationFile} F
 * @param {ConservationSourcePlace<F>} place
 * @returns {F[]}
 */
const filesOfPlace = (place) => {
  if (Array.isArray(place)) {
    return /** @type {readonly (F | null)[]} */ (place).filter(/** @returns {file is F} */ (file) => Boolean(file));
  }
  return place ? [/** @type {F} */ (place)] : [];
};

/**
 * The one removal of a comparison series with its files: Remove in the
 * Circular comparison panel, and Load for a series whose BLAST table it could
 * not read (OV-893, E-OD30). Each list holds one place per series, in series
 * order; the place `index` goes from every list, and the series after it keep
 * their own files and take the place before (`sourceIndex`). The comparison
 * stays on while a series keeps its BLAST table (`enabled`).
 * @template B The place of a BLAST table.
 * @template G The place of a ring genome's files.
 * @template {object} S
 * @param {{ blasts: B[], fastas: G[], sequenceSources: G[], series: S[] }} lists
 * @param {number} index
 */
export const withoutConservationSeries = ({ blasts, fastas, sequenceSources, series }, index) => {
  /** @template T @param {T[]} list */
  const without = (list) => list.filter((_, place) => place !== index);
  const kept = without(series).map((entry, sourceIndex) => ({ ...entry, sourceIndex }));
  const keptBlasts = without(blasts);
  return {
    blasts: keptBlasts,
    fastas: without(fastas),
    sequenceSources: without(sequenceSources),
    series: kept,
    enabled: keptBlasts.length > 0
  };
};

/**
 * @param {{ name?: unknown } | null | undefined} file A File, or a series entry's file name; only `name` is read.
 * @param {unknown} index
 */
export const defaultConservationSeriesLabel = (file, index) => {
  const raw = String(file?.name || `Series ${Number(index) + 1}`).trim();
  const withoutExtension = raw.replace(/\.[^.]+$/, '').trim();
  return withoutExtension || `Series ${Number(index) + 1}`;
};

/**
 * @param {unknown} value
 * @returns {string[]}
 */
export const parseConservationLabelText = (value) =>
  String(value || '')
    .split(/[\n,]+/)
    .map((label) => label.trim())
    .filter(Boolean);

/**
 * @param {unknown} value
 * @param {number} [fallback]
 */
const normalizePositiveInteger = (value, fallback = 1) => {
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : fallback;
};

/**
 * @param {ConservationFile} file `normalizeFileList` drops falsy entries
 * @param {unknown} index
 */
const fileSignature = (file, index) => {
  const fileName = String(file?.name || `source_${Number(index) + 1}`).trim();
  const size = Number.isFinite(Number(file?.size)) ? Number(file.size) : 0;
  const lastModified = Number.isFinite(Number(file?.lastModified)) ? Number(file.lastModified) : 0;
  return `${fileName}|${size}|${lastModified}`;
};

/**
 * The sources of a source list, in order; a place without a file is left out,
 * so `sourceIndex` counts the sources and `placeIndex` names the place. A source's key is its files'
 * signatures (one file's is that file's, as before a ring took several), with
 * the count of earlier sources of the same files, which shifts when one of them
 * goes (`filesOfSourceKey`).
 * @template {ConservationFile} F
 * @param {ConservationSourcePlace<F> | readonly ConservationSourcePlace<F>[]} sourceFiles
 * @returns {ConservationSourceDescriptor<F>[]}
 */
export const conservationSourceDescriptors = (sourceFiles) => {
  /** @type {Map<string, number>} */
  const seen = new Map();
  /** @type {readonly ConservationSourcePlace<F>[]} */
  const places = Array.isArray(sourceFiles) ? sourceFiles : [sourceFiles];
  return places.flatMap((place, placeIndex) => {
    const files = filesOfPlace(place);
    return files.length > 0 ? [{ files, placeIndex }] : [];
  }).map(({ files, placeIndex }, index) => {
    const [file] = files;
    const signature = files.map((member) => fileSignature(member, index)).join(' + ');
    const duplicateIndex = seen.get(signature) || 0;
    seen.set(signature, duplicateIndex + 1);
    const fileName = String(file?.name || `source_${Number(index) + 1}`).trim();
    const sourceKey = `${signature}|${duplicateIndex}`;
    return {
      file,
      files,
      sourceKey,
      sourceIndex: index,
      placeIndex,
      fileName,
      defaultLabel: defaultConservationSeriesLabel(file, index)
    };
  });
};

/**
 * A source key's files: the key less its count of earlier sources of the same files.
 * @param {string} sourceKey
 */
const filesOfSourceKey = (sourceKey) => sourceKey.replace(/\|\d+$/, '');

/**
 * The source key of each ring, in ring order: the `sourceKey` `conservationSourceDescriptors` gives the ring's file
 * among the ring files (a BLAST table of a stored request; the file's name, size and date, copies of one file by their
 * order). A ring whose file is not held has no key: no entry is guessed for it.
 * @param {readonly (ConservationFile | null | undefined)[]} ringFiles The ring files, one place per ring, in ring
 *   order.
 * @returns {(string | null)[]}
 */
export const ringKeysOfFiles = (ringFiles) => {
  const keys = conservationSourceDescriptors(ringFiles.filter(/** @returns {file is ConservationFile} */ (file) => Boolean(file)))
    .map(({ sourceKey }) => sourceKey);
  return ringFiles.map((file) => (file ? /** @type {string} */ (keys.shift()) : null));
};

/**
 * A series entry of a source: a kept entry keeps its colour and scheme index (a missing index reads as its place),
 * a new one (`previousEntry` null) has no colour of its own and the first scheme index the series leaves (D-05).
 * @param {ConservationSourceDescriptor} descriptor
 * @param {RawConservationSeriesEntry | null} previousEntry
 * @param {readonly { schemeIndex?: unknown }[]} series The entries before it.
 * @param {string | null} [legacyLabel]
 * @returns {ConservationSeriesEntry}
 */
const mergeSeriesEntry = (descriptor, previousEntry, series, legacyLabel = null) => {
  const rawLabel = previousEntry?.label ?? previousEntry?.name ?? legacyLabel ?? descriptor.defaultLabel;
  return {
    sourceKey: descriptor.sourceKey,
    fileName: descriptor.fileName,
    // The entry's stored place (`conservationSourcesOfSeries` reads it), which an empty place before it keeps.
    sourceIndex: descriptor.placeIndex,
    label: String(rawLabel || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
    color: normalizeOptionalHexColor(previousEntry?.color),
    schemeIndex: previousEntry ? readSchemeIndex(previousEntry.schemeIndex, series.length) : firstUnusedSchemeIndex(series),
    losat_gencode: normalizePositiveInteger(previousEntry?.losat_gencode, 1)
  };
};

/**
 * How a change of the rings moves their source keys, by which the track stack's managed ring rows name their rings
 * (`series_key`): the key a ring's row carries before the change to the key of that ring's source after it, `''` for a
 * ring that goes. A source key names files, and its count of earlier sources of the same files shifts when one of them
 * goes, so a row follows the matcher's result (`conservationSourcesOfSeries`), never the key it carried; a key not in
 * the map stays.
 * @typedef {ReadonlyMap<string, string>} RingKeyChanges
 */
/** @type {RingKeyChanges} */
export const NO_RING_KEY_CHANGES = new Map();
/**
 * The key a ring row carries after a change (`RingKeyChanges`).
 * @param {RingKeyChanges} keyChanges @param {string} key
 */
export const ringKeyAfter = (keyChanges, key) => (keyChanges.has(key) ? String(keyChanges.get(key)) : key);
/**
 * The key changes of two changes made in turn: a row's key goes through `first`, then through `second`.
 * @param {RingKeyChanges} first @param {RingKeyChanges} second
 * @returns {RingKeyChanges}
 */
export const thenRingKeyChanges = (first, second) => new Map([
  ...second,
  ...[...first].map(([before, between]) => /** @type {[string, string]} */ ([before, ringKeyAfter(second, between)]))
]);
/**
 * The key changes of a change that keeps the series entries in their order, less the one at `removed`: the key of each
 * entry's source before to the key of the source of the entry it becomes, `''` for the removed entry and for one left
 * without a source.
 * @param {readonly ({ sourceKey?: unknown } | null | undefined)[]} before each entry's source before (or the entry, whose
 *   key names it)
 * @param {readonly ({ sourceKey?: unknown } | null | undefined)[]} after each entry's source after, the removed one left out
 * @param {number} [removed] the entry the change removes; none by default
 * @returns {RingKeyChanges}
 */
export const ringKeyChangesOf = (before, after, removed = -1) => new Map(before.flatMap((source, index) => {
  const key = String(source?.sourceKey || '').trim();
  if (!key) return [];
  const next = index === removed ? null : after[removed >= 0 && index > removed ? index - 1 : index];
  return [/** @type {[string, string]} */ ([key, String(next?.sourceKey || '').trim()])];
}));

/**
 * The series with each entry naming the source at its index of `sources`, all three fields of its identity from that
 * one source: its key, its first file's name and its place (`sourceIndex`); an entry given no source stays as it is.
 *
 * The rings' identity invariant, which every ring transition keeps and `tests/web/ring-identity-walk.test.mjs`
 * checks after each step: (a) each series entry's `sourceKey`, `fileName` and `sourceIndex` name one and the same
 * source of the drawing's current ring sources (`circularRingSources`), and no two entries name one source; (b) each
 * managed stack row's `series_key` is the key the one matcher (`conservationSourcesOfSeries`) gives its entry over
 * those sources. With (a) the matcher gives each entry the source it names, so an entry's stored key is the key its row
 * carries, which the reconcile's key changes read (`reconcileConservationSeries`).
 * @param {readonly RawConservationSeriesEntry[]} series
 * @param {readonly (ConservationSourceDescriptor | null | undefined)[]} sources
 * @returns {RawConservationSeriesEntry[]}
 */
export const seriesNamingSources = (series, sources) => series.map((entry, index) => {
  const source = sources[index];
  return source ? { ...entry, sourceKey: source.sourceKey, fileName: source.fileName, sourceIndex: source.placeIndex } : entry;
});

/**
 * The label the stored label list (`circularConservation.labels`) gives a source that no series entry claims. The list
 * names the rings by their series index, as the writers write it from the series (Generate's adoption, the reconcile)
 * and as an older Session kept it for places past its series (OV-898). Its label at an index the series covers is
 * that entry's ring, never a new source, so such a source takes its file's name (OV-3166).
 * @param {readonly string[]} labels @param {number} sourceIndex @param {number} seriesLength
 * @returns {string | null}
 */
const listLabelOfUnclaimedSource = (labels, sourceIndex, seriesLength) => (
  sourceIndex >= seriesLength ? labels[sourceIndex] ?? null : null
);
/**
 * The entries the series gains for the sources no entry of it claims, after its `seriesLength` entries, each labelled
 * by the stored label list where an older Session named it there (`listLabelOfUnclaimedSource`), else by its file. Both
 * reconciles, of a change of the files and of a LOSAT-cache replay, fold the list into the series this way, so the list
 * names no ring after them.
 * @param {readonly ConservationSourceDescriptor[]} unclaimed
 * @param {readonly string[]} labels The stored label list, parsed (`parseConservationLabelText`).
 * @param {number} seriesLength The entries of the series the list was stored with.
 * @param {readonly { schemeIndex?: unknown }[]} kept The entries kept before them, whose scheme indexes a gained entry
 *   does not take (D-05).
 * @returns {ConservationSeriesEntry[]}
 */
export const seriesEntriesOfUnclaimedSources = (unclaimed, labels, seriesLength, kept) => unclaimed.reduce(
  (/** @type {ConservationSeriesEntry[]} */ gained, descriptor) => [...gained, mergeSeriesEntry(descriptor, null, [...kept, ...gained],
    listLabelOfUnclaimedSource(labels, descriptor.sourceIndex, seriesLength))],
  []
);

/**
 * The series a change of the ring files leaves, and how the change moves the rings' keys: each previous entry's key
 * (the key its row carries: the rings' identity invariant, `seriesNamingSources`, held before the change, and a
 * transition that re-keyed the rows first wrote its entries the same keys) to the key of the source the matcher gives
 * it, `''` when it has none.
 * @param {{
 *   sourceFiles: ConservationSourcePlace<ConservationFile> | readonly ConservationSourcePlace<ConservationFile>[],
 *   previousSeries: RawConservationSeriesEntry[] | null | undefined,
 *   legacyLabels: string[] | null | undefined
 * }} options
 * @returns {{ series: ConservationSeriesEntry[], keyChanges: RingKeyChanges }}
 */
export const reconcileConservationSeries = ({
  sourceFiles,
  previousSeries,
  legacyLabels
}) => {
  const previous = Array.isArray(previousSeries) ? previousSeries : [];
  // A change of the files drops an entry whose files went: a stored place lends it no file of another name.
  const { matched, unclaimed } = conservationSourcesOfSeries(sourceFiles, { series: previous }, { storedPlace: 'same name' });
  const fallbackLabels = Array.isArray(legacyLabels) ? legacyLabels : [];
  /** @type {ConservationSeriesEntry[]} */
  const nextSeries = [];
  matched.forEach((descriptor, seriesIndex) => {
    if (!descriptor) return;
    nextSeries.push(mergeSeriesEntry(descriptor, previous[seriesIndex], nextSeries, fallbackLabels[descriptor.sourceIndex]));
  });
  nextSeries.push(...seriesEntriesOfUnclaimedSources(unclaimed, fallbackLabels, previous.length, nextSeries));
  return { series: nextSeries, keyChanges: ringKeyChangesOf(previous, matched) };
};

// D12: a row added from a GenBank or DDBJ file takes the label the file names
// itself (first record's DEFINITION, then organism), read by the Python ring
// reader. Only the row's file-name default is replaced, so a typed label wins;
// a FASTA file names no label and keeps the file-name default.
/**
 * @param {{
 *   series: ConservationSeriesEntry[] | null | undefined,
 *   sourceFiles: ConservationSourcePlace<ConservationFile> | readonly ConservationSourcePlace<ConservationFile>[],
 *   file: ConservationFile | null | undefined,
 *   recordLabel: unknown
 * }} options
 */
export const applyComparisonSequenceRecordLabel = ({ series, sourceFiles, file, recordLabel }) => {
  const label = String(recordLabel ?? '').trim();
  if (!label || !file || !Array.isArray(series)) return false;
  // The row of the file's source as Generate matches it (`conservationSourcesOfSeries`), not by an equal key.
  const { matched } = conservationSourcesOfSeries(sourceFiles, { series });
  const seriesIndex = matched.findIndex((candidate) => candidate?.file === file);
  const entry = seriesIndex >= 0 ? series[seriesIndex] : null;
  const descriptor = seriesIndex >= 0 ? matched[seriesIndex] : null;
  if (!entry || !descriptor || String(entry.label ?? '').trim() !== descriptor.defaultLabel) return false;
  entry.label = label;
  return true;
};

/**
 * How the matcher reads a stored place (`conservationSourcesOfSeries`). Generate, the rows, Save and the span export
 * take the source at an entry's stored place (`storedPlace: 'any'`): a Session writes one place per ring, and a
 * LOSAT-cache replay's tables sit at their rings' places under names of other files. The series reconcile, which
 * follows a change of the files, takes it only when its first file still has the entry's file name
 * (`'same name'`), so an entry whose files went leaves with them. Load's repair of the tables it could not read lets
 * an entry left without a source take the empty place it names (`emptyPlaces`), so that entry is the one reported.
 * @typedef {{ readonly storedPlace?: 'any' | 'same name', readonly emptyPlaces?: boolean }} ConservationMatchPasses
 */
/**
 * @template {ConservationFile} F
 * @typedef {{ descriptors: ConservationSourceDescriptor<F>[], matched: (ConservationSourceDescriptor<F> | null)[],
 *   unclaimed: ConservationSourceDescriptor<F>[], places: (number | null)[], unplaced: number[] }} ConservationSeriesMatch
 */

/**
 * The source of each series entry, in series order: the one matching of a ring
 * or a BLAST row to its files, which Generate, the track slots, the ring rows,
 * Save, the span export catalog, the series reconcile and Load's repair of
 * unread tables read. It works in passes over every entry, so that no entry
 * takes a source a later entry names more closely: each entry first takes the
 * source of its source key when that key is unambiguous, its files carried by
 * one source and one entry (three rings of NCBI `sequence.fasta` files of one
 * size and date are not: their keys differ only in a count of earlier ones,
 * which Remove series or a file Load could not read shifts); then an entry
 * left takes the source of its file name when that name is unambiguous, carried by one source no entry claimed
 * and by one entry left (two rings whose first files are both NCBI's
 * `sequence.fasta`, or whose replayed tables are named after them, are not);
 * then the source at the place its source index names, which a move keeps
 * (`passes`). `matched[i]` is entry i's source, `null` when none is left, and
 * `places[i]` its place; `unclaimed` are the sources no entry takes and
 * `unplaced` the places no entry takes, in order.
 * @template {ConservationFile} F
 * @param {ConservationSourcePlace<F> | readonly ConservationSourcePlace<F>[]} sourceFiles
 * @param {CircularConservationLike | null | undefined} circularConservation
 * @param {ConservationMatchPasses} [passes]
 * @returns {ConservationSeriesMatch<F>}
 */
export const conservationSourcesOfSeries = (sourceFiles, circularConservation, { storedPlace = 'any', emptyPlaces = false } = {}) => {
  // The places' union hides F from the call's inference; the descriptors hold the given files.
  const descriptors = /** @type {ConservationSourceDescriptor<F>[]} */ (conservationSourceDescriptors(sourceFiles));
  /** @type {readonly ConservationSourcePlace<F>[]} */
  const allPlaces = Array.isArray(sourceFiles) ? sourceFiles : [sourceFiles];
  const series = Array.isArray(circularConservation?.series) ? circularConservation.series : [];
  /** @type {Set<string>} */
  const claimed = new Set();
  /** @type {(ConservationSourceDescriptor<F> | null)[]} */
  const matched = series.map(() => null);
  /** @param {unknown} entry @returns {entry is RawConservationSeriesEntry} */
  const isEntry = (entry) => Boolean(entry) && typeof entry === 'object';
  /** @param {RawConservationSeriesEntry} entry */
  const nameOf = (entry) => String(entry.fileName || '').trim();
  /** @param {RawConservationSeriesEntry} entry @param {number} seriesIndex */
  const storedPlaceOf = (entry, seriesIndex) => {
    const sourceIndex = Number(entry.sourceIndex);
    return Number.isInteger(sourceIndex) ? sourceIndex : seriesIndex;
  };
  /**
   * One pass: each entry without a source yet takes the unclaimed source `fits` names.
   * @param {(entry: RawConservationSeriesEntry, seriesIndex: number) => (source: ConservationSourceDescriptor<F>) => boolean} fits
   */
  const claim = (fits) => series.forEach((entry, seriesIndex) => {
    if (matched[seriesIndex] || !isEntry(entry)) return;
    const fitsEntry = fits(entry, seriesIndex);
    const source = descriptors.find((candidate) => !claimed.has(candidate.sourceKey) && fitsEntry(candidate)) || null;
    if (!source) return;
    claimed.add(source.sourceKey);
    matched[seriesIndex] = source;
  });
  /** @param {string[]} names @returns {Map<string, number>} */
  const tally = (names) => names.reduce((counts, name) => counts.set(name, (counts.get(name) || 0) + 1), new Map());
  /** @param {RawConservationSeriesEntry} entry */
  const keyOf = (entry) => String(entry.sourceKey || '').trim();
  const sourcesOfFiles = tally(descriptors.map((source) => filesOfSourceKey(source.sourceKey)));
  const entriesOfFiles = tally(series.filter(isEntry).map(keyOf).filter(Boolean).map(filesOfSourceKey));
  claim((entry) => {
    const sourceKey = keyOf(entry);
    const files = filesOfSourceKey(sourceKey);
    const unambiguous = sourceKey !== '' && sourcesOfFiles.get(files) === 1 && entriesOfFiles.get(files) === 1;
    return (source) => unambiguous && source.sourceKey === sourceKey;
  });
  const sourcesOfName = tally(descriptors.filter((source) => !claimed.has(source.sourceKey)).map((source) => source.fileName));
  const entriesOfName = tally(series.filter((entry, seriesIndex) => !matched[seriesIndex] && isEntry(entry)).map(nameOf));
  claim((entry) => {
    const fileName = nameOf(entry);
    const unambiguous = fileName !== '' && sourcesOfName.get(fileName) === 1 && entriesOfName.get(fileName) === 1;
    return (source) => unambiguous && source.fileName === fileName;
  });
  claim((entry, seriesIndex) => {
    const place = storedPlaceOf(entry, seriesIndex);
    return (source) => source.placeIndex === place && (storedPlace === 'any' || source.fileName === nameOf(entry));
  });
  const places = matched.map((source) => (source ? source.placeIndex : null));
  const placeHoldsNothing = (/** @type {number} */ place) => place >= 0 && place < allPlaces.length && filesOfPlace(allPlaces[place]).length === 0;
  if (emptyPlaces) series.forEach((entry, seriesIndex) => {
    if (matched[seriesIndex] || !isEntry(entry)) return;
    const place = storedPlaceOf(entry, seriesIndex);
    if (placeHoldsNothing(place) && !places.includes(place)) places[seriesIndex] = place;
  });
  const unclaimed = descriptors.filter((source) => !claimed.has(source.sourceKey));
  const unplaced = emptyPlaces
    ? allPlaces.map((_, place) => place).filter((place) => !places.includes(place))
    : unclaimed.map((source) => source.placeIndex);
  return { descriptors, matched, unclaimed, places, unplaced };
};

/**
 * @template {ConservationFile} F
 * @param {ConservationSourcePlace<F> | readonly ConservationSourcePlace<F>[]} sourceFiles
 * @param {CircularConservationLike | null | undefined} circularConservation
 * @returns {OrderedConservationSource<F>[]}
 */
export const orderedConservationSources = (sourceFiles, circularConservation) => {
  // The places' union hides F from the call's inference, as in `conservationSourcesOfSeries`.
  const { matched, unclaimed } = /** @type {ConservationSeriesMatch<F>} */ (conservationSourcesOfSeries(sourceFiles, circularConservation));
  const legacyLabels = parseConservationLabelText(circularConservation?.labels);
  const series = Array.isArray(circularConservation?.series) ? circularConservation.series : [];
  const scheme = circularConservation?.colorScheme;
  // Each entry's scheme index as `reconcileConservationSeries` reads it: a stored one, else its place; a source no
  // entry claims takes the first one left (D-05).
  /** @type {OrderedConservationSource<F>[]} */
  const entries = [];
  matched.forEach((descriptor, seriesIndex) => {
    if (!descriptor) return;
    const entry = series[seriesIndex];
    const schemeIndex = readSchemeIndex(entry.schemeIndex, entries.length);
    const ownColor = normalizeOptionalHexColor(entry.color);
    entries.push({
      ...descriptor,
      label: String(entry.label || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
      color: effectiveRingColor({ color: ownColor, schemeIndex }, scheme) ?? '',
      ownColor,
      schemeIndex,
      losat_gencode: normalizePositiveInteger(entry.losat_gencode, 1),
      orderIndex: entries.length
    });
  });
  unclaimed.forEach((descriptor) => {
    const fallbackLabel = listLabelOfUnclaimedSource(legacyLabels, descriptor.sourceIndex, series.length) || descriptor.defaultLabel;
    const schemeIndex = firstUnusedSchemeIndex(entries);
    entries.push({
      ...descriptor,
      label: String(fallbackLabel || descriptor.defaultLabel).trim() || descriptor.defaultLabel,
      color: effectiveRingColor({ schemeIndex }, scheme) ?? '',
      ownColor: null,
      schemeIndex,
      losat_gencode: 1,
      orderIndex: entries.length
    });
  });
  return entries;
};

/**
 * What `ringGenomePlaces` reads of a Circular drawing.
 * @template {ConservationFile} [F=InputFileRef]
 * @typedef {RingSourceDrawing<F> & {
 *   files: RingSourceDrawing<F>['files'] & { readonly c_conservation_sequence_sources: readonly (readonly F[])[] },
 *   circularConservation: CircularConservationLike & { readonly source?: unknown }
 * }} RingGenomeDrawing
 */
/** @typedef {'c_conservation_fastas' | 'c_conservation_sequence_sources'} RingGenomeKey */
/**
 * Where each ring's comparison genome files are, in series order: the list
 * that holds them (a LOSAT ring's own files, or the comparison sequence files
 * beside an uploaded or replayed BLAST table) and each series entry's place in
 * it, or `null` for an entry without a source. Entries are matched to sources
 * as Generate matches them (`conservationSourcesOfSeries`, in passes over every
 * entry: source keys, then unambiguous file names, then the place a Session wrote), so a
 * place left empty on Load keeps the later rings at theirs. `unplaced` are the places of LOSAT rings
 * that no entry names, which Generate draws after the others.
 * @template {ConservationFile} F
 * @param {RingGenomeDrawing<F>} drawing
 * @returns {{ key: RingGenomeKey, lists: F[][], places: (number | null)[], unplaced: number[] }}
 */
export const ringGenomePlaces = (drawing) => {
  const tables = ringsFollowTables(drawing);
  /** @type {RingGenomeKey} */
  const key = isUploadedRingDrawing(drawing) ? 'c_conservation_sequence_sources' : 'c_conservation_fastas';
  /** @type {F[][]} */
  const lists = ringGenomeFileLists(drawing.files[key]);
  // A table's place is its genome files' place, beside it.
  const { places, unplaced } = conservationSourcesOfSeries(
    tables ? normalizeFileList([...drawing.files.c_conservation_blasts]) : lists, drawing.circularConservation
  );
  return { key, lists, places, unplaced: tables ? [] : unplaced };
};

/**
 * Each ring's comparison genome files in the order Generate draws the rings,
 * which Save and the span export catalog read: a LOSAT ring's own files, one
 * list per series entry (empty for an entry whose files are missing, as a
 * Session loaded without them shows its stored rings), then the rings no entry
 * names; or the comparison sequence files at each uploaded or replayed table's
 * place, in the tables' ring order (`orderedConservationSources`).
 * @template {ConservationFile} F
 * @param {RingGenomeDrawing<F>} drawing
 * @returns {F[][]}
 */
export const ringGenomeFilesInRingOrder = (drawing) => {
  const { lists, places, unplaced } = ringGenomePlaces(drawing);
  if (ringsFollowTables(drawing)) {
    const tables = orderedConservationSources(normalizeFileList([...drawing.files.c_conservation_blasts]), drawing.circularConservation);
    return tables.map((table) => [...(lists[table.placeIndex] || [])]);
  }
  return [...places, ...unplaced].map((place) => (place === null ? [] : [...(lists[place] || [])]));
};

/**
 * The rings of a drawing as one string, equal only for the same rings in the same order: each by the key its series
 * entry holds and by the ids of its genome files at the place Generate matches it to (`ringGenomeFilesInRingOrder`).
 * A Generate reads it to tell whether the draft's rings are still those its searches read. The keys alone do not tell
 * rings apart: a key names files by name, size and date, so after one of three rings of such files goes, the two left
 * and a new one carry the keys the three had.
 * @param {RingGenomeDrawing} drawing
 */
export const ringsOfDrawing = (drawing) => JSON.stringify([
  (Array.isArray(drawing.circularConservation.series) ? drawing.circularConservation.series : []).map((entry) => String(entry.sourceKey)),
  ringGenomeFilesInRingOrder(drawing).map((files) => files.map((file) => file.fileId))
]);

/**
 * @template T
 * @param {T[]} series
 * @param {unknown} index
 * @param {unknown} direction
 */
export const moveConservationSeriesEntry = (series, index, direction) => {
  const idx = Number(index);
  const step = Math.sign(Number(direction));
  if (!Array.isArray(series) || !Number.isInteger(idx) || step === 0) return false;
  const target = idx + step;
  if (idx < 0 || idx >= series.length || target < 0 || target >= series.length) return false;
  const [entry] = series.splice(idx, 1);
  series.splice(target, 0, entry);
  return true;
};

/**
 * @param {{ label?: unknown, fileName?: unknown } | null | undefined} entry
 * @param {unknown} index
 * @param {Set<string>} [existingIds]
 */
export const safeConservationSlotId = (entry, index, existingIds = new Set()) => {
  const source = String(entry?.label || entry?.fileName || `comparison_${Number(index) + 1}`).trim();
  const cleaned = source
    .replace(/\.[^.]+$/, '')
    .replace(/[^A-Za-z0-9_]+/g, '_')
    .replace(/^_+|_+$/g, '')
    .slice(0, 40) || `comparison_${Number(index) + 1}`;
  let id = `conservation_${cleaned}`;
  let suffix = 2;
  while (existingIds.has(id)) {
    id = `conservation_${cleaned}_${suffix}`;
    suffix += 1;
  }
  existingIds.add(id);
  return id;
};

/** @param {{ renderer?: unknown, params?: { managed?: unknown } | null } | null | undefined} slot */
export const isManagedConservationSlot = (slot) => (
  slot?.renderer === 'sequence_conservation' &&
  String(slot?.params?.managed || '').trim() === CONSERVATION_SLOT_MANAGER
);
