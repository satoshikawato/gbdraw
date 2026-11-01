// @ts-check
// E1: one drawing's generated artifact, which the composition root keeps, by
// drawing ID, while another drawing is shown. The History snapshot service
// captures and installs it (`captureArtifactSlot`, `installArtifactSlot`);
// Session Save writes it and Session Load builds it (`services/config.js`).
/** @import { DrawingId, DrawingIdentity } from './drawing-view.js' */

/**
 * A drawing's feature-editor status, a runtime field of its artifact slot (never saved in a Session).
 * `state.setFeatureEditorStatus` writes it; the Editor header reads it (`featureEditorStatusText`).
 * @typedef {object} FeatureEditorStatus
 * @property {'idle' | 'summary-ready' | 'failed'} status `failed`: Load could not recover the feature metadata.
 * @property {string | null} error The failure, or the warning of a recovered Load (the Editor header's title).
 */
/** @type {Readonly<FeatureEditorStatus>} */
export const FEATURE_EDITOR_STATUS_IDLE = Object.freeze({ status: 'idle', error: null });

/**
 * E1: the generated artifact of one drawing, which each drawing keeps while
 * another is shown. The state ref of every key and the value an empty slot
 * installs, given the mode of the slot's drawing and the current value: the displayed-Result
 * part of the owner set, the selected Result, the feature-editor status, the generated Legend colors and
 * stroke defaults, and the alignment Reset receipt (it binds the
 * committed Linear Session). Draft and project state stays shared: the specific
 * color rules, the file Legend captions, the Linear record orientations, the
 * pending palette, the LOSAT caches and manifests, the alignment plan and
 * record translations, and every editor intent (R2).
 * @satisfies {Readonly<Record<string, (mode: 'circular' | 'linear', current: any) => any>>}
 */
export const ARTIFACT_SLOT_EMPTY = Object.freeze({
  results: () => [],
  selectedResultIndex: () => 0,
  featureCatalog: () => null,
  extractedFeatures: () => [],
  featureEditorStatus: () => FEATURE_EDITOR_STATUS_IDLE,
  biologicalFeatures: () => [],
  featureRecordIds: () => [],
  orthogroups: () => [],
  featureOrthogroupIndex: () => new Map(),
  collinearGroups: () => [],
  trackSlotResolvedGeometry: () => null,
  annotationWarnings: () => [],
  featureIdentityNotices: () => [],
  comparisonWarnings: () => [],
  labelFitReports: () => [],
  inputSummaries: () => [],
  lastRunInfo: () => null,
  pairwiseMatchFactors: () => ({}),
  editableLabels: () => [],
  generatedLegendPosition: () => 'left',
  generatedMultiRecordCanvas: () => false,
  generatedCircularPlotTitlePosition: () => 'none',
  // An empty slot has no Result to project a palette onto: it keeps the applied one,
  // as a slot whose applied palette is null does (`installArtifactSlot`).
  appliedPaletteName: (_mode, current) => current,
  appliedPaletteColors: (_mode, current) => current,
  similarityAlignmentResetReceipt: () => null,
  originalLegendColors: () => ({}),
  originalSvgStroke: () => ({ color: null, width: null })
});
// The cast: `ARTIFACT_SLOT_EMPTY` is the frozen literal above, so its own keys are exactly the keys of its type.
export const ARTIFACT_SLOT_KEYS = Object.freeze(/** @type {(keyof typeof ARTIFACT_SLOT_EMPTY)[]} */ (Object.keys(ARTIFACT_SLOT_EMPTY)));

/**
 * One drawing's generated artifact, held by reference (no clone).
 * @typedef {object} ArtifactSlot
 * @property {DrawingId} drawingId The drawing whose artifact it is (its mode is the drawing's).
 * @property {Readonly<Record<string, any>>} values The value of every `ARTIFACT_SLOT_KEYS` key.
 * @property {readonly string[]} legendInventory The selected Result's generated Legend inventory
 *   (`originalLegendOrder`). The Legend owner keeps each Result's inventory (OV-47), so an install
 *   does not write it: the mode transition hands it to the Legend owner, which adopts it when the
 *   Result is displayed.
 * @property {any} matchSequenceOwner The match-sequence registry's trusted owner, or null.
 * @property {Record<string, any> | null} runtimeState The runtime owner's state (committed Session, CLI helper files).
 * @property {Readonly<Record<string, any>> | null} transportIdentity The History transport identity of the artifact.
 * @property {number} retainedBytes
 */

/**
 * The artifact slots a Session Load installs before its preview mounts:
 * `opening`, the artifact the Session opens on when it is not the committed
 * set's, and `stashed`, the ones kept for the drawings not shown.
 * @typedef {{ opening: Readonly<ArtifactSlot> | null, stashed: readonly Readonly<ArtifactSlot>[] }} LoadedArtifactSlots
 */

/**
 * A slot of `drawing`'s artifact; a value it is not given is the empty
 * slot's for the drawing's mode.
 * @param {{ drawing: DrawingIdentity, values: ArtifactSlot['values'], legendInventory?: ArtifactSlot['legendInventory'],
 *   matchSequenceOwner?: ArtifactSlot['matchSequenceOwner'], runtimeState?: ArtifactSlot['runtimeState'],
 *   transportIdentity?: ArtifactSlot['transportIdentity'], retainedBytes?: number }} slot
 * @returns {Readonly<ArtifactSlot>}
 */
export const createArtifactSlot = ({
  drawing, values, legendInventory = [], matchSequenceOwner = null, runtimeState = null,
  transportIdentity = null, retainedBytes = 0
}) => {
  // The cast widens the key list, so that a field of any name can be looked up in it.
  const unknown = Object.keys(values).filter((key) => !(/** @type {readonly string[]} */ (ARTIFACT_SLOT_KEYS)).includes(key));
  if (unknown.length > 0) throw new Error(`Unknown artifact slot field(s): ${unknown.join(', ')}.`);
  return Object.freeze({
    drawingId: drawing.id,
    values: Object.freeze(Object.fromEntries(ARTIFACT_SLOT_KEYS.map((key) => [
      key, Object.hasOwn(values, key) ? values[key] : ARTIFACT_SLOT_EMPTY[key](drawing.mode, null)
    ]))),
    legendInventory: Object.freeze([...legendInventory]),
    matchSequenceOwner,
    runtimeState,
    transportIdentity,
    retainedBytes: Math.max(0, Number(retainedBytes) || 0)
  });
};

/**
 * `slot` with `values` in place of its own, or, with no slot, `drawing`'s
 * empty slot with `values`.
 * @param {Readonly<ArtifactSlot> | null} slot
 * @param {DrawingIdentity} drawing
 * @param {Readonly<Record<string, unknown>>} values
 * @returns {Readonly<ArtifactSlot>}
 */
const withArtifactSlotValues = (slot, drawing, values) => createArtifactSlot({
  drawing,
  values: { ...slot?.values, ...values },
  legendInventory: slot?.legendInventory,
  matchSequenceOwner: slot?.matchSequenceOwner ?? null,
  runtimeState: slot?.runtimeState ?? null,
  transportIdentity: slot?.transportIdentity ?? null,
  retainedBytes: slot?.retainedBytes ?? 0
});

/**
 * `slot`'s artifact as `drawing`'s (Duplicate: the copy shows its source's
 * Result without a Generate). The artifact's values are shared, as a slot
 * keeps them by reference; only the drawing it belongs to differs.
 * @param {Readonly<ArtifactSlot>} slot
 * @param {DrawingIdentity} drawing
 * @returns {Readonly<ArtifactSlot>}
 */
export const artifactSlotAs = (slot, drawing) => withArtifactSlotValues(slot, drawing, {});

/**
 * What `ArtifactSlotStore.snapshot` returns for `restore`.
 * @typedef {readonly Readonly<ArtifactSlot>[]} ArtifactSlotSnapshot
 */

/**
 * The kept artifacts (E1): the artifact of each drawing that is not shown,
 * by its drawing ID and by reference, until the mode transition installs it.
 * A drawing has no entry while it is shown or before it has an artifact, and
 * an entry's key is its slot's `drawingId`. Only these operations change it.
 * @typedef {object} ArtifactSlotStore
 * @property {() => readonly Readonly<ArtifactSlot>[]} kept Every kept slot.
 * @property {<T>(drawingId: DrawingId, shownDrawingId: DrawingId, readShown: () => T,
 *   readKept: (slot: Readonly<ArtifactSlot>) => T) => T | undefined} readOfDrawing
 *   A value of drawing `drawingId`'s artifact: read from the displayed artifact while that drawing
 *   is shown, otherwise from the slot that keeps it (undefined for an empty slot). Every read of
 *   "this drawing's artifact, shown or kept" goes through here. A value kept in a slot is the same
 *   object when the drawing is shown again.
 * @property {(departing: Readonly<ArtifactSlot>, arrivingId: DrawingId) => Readonly<ArtifactSlot> | null} swap
 *   The mode transition: keeps the departing drawing's artifact and takes out the arriving one's
 *   (null: it has none yet).
 * @property {(stashed: readonly Readonly<ArtifactSlot>[]) => void} replaceLoaded
 *   Session Load: the loaded Session's kept artifacts replace every slot.
 * @property {() => ArtifactSlotSnapshot} snapshot
 * @property {(snapshot: ArtifactSlotSnapshot) => void} restore The Load rollback.
 * @property {(drawing: DrawingIdentity, values: Readonly<Record<string, unknown>>) => void} writeValues
 *   Writes slot values of a drawing that is not shown (History and the Instant Preview toggle;
 *   OV-1168, OV-1195); with no slot, into the drawing's empty slot.
 * @property {(drawingId: DrawingId) => Readonly<ArtifactSlot> | null} takeSlot
 *   A removed drawing's slot goes with it (`removeDrawing` in state.js); the
 *   step that removed the drawing holds it (null: it has none) for its Undo.
 * @property {(slot: Readonly<ArtifactSlot>) => void} putSlot
 *   A held slot back, with its drawing (the Undo of a Delete, the Redo of a New).
 */

/** @returns {ArtifactSlotStore} */
export const createArtifactSlotStore = () => {
  /** @type {Map<DrawingId, Readonly<ArtifactSlot>>} */
  const slots = new Map();
  /** @param {Readonly<ArtifactSlot>} slot */
  const keep = (slot) => { slots.set(slot.drawingId, slot); };
  const kept = () => Object.freeze([...slots.values()]);
  return {
    kept,
    readOfDrawing: (drawingId, shownDrawingId, readShown, readKept) => {
      if (shownDrawingId === drawingId) return readShown();
      const kept = slots.get(drawingId);
      return kept ? readKept(kept) : undefined;
    },
    swap: (departing, arrivingId) => {
      keep(departing);
      const arriving = slots.get(arrivingId) ?? null;
      slots.delete(arrivingId);
      return arriving;
    },
    replaceLoaded: (stashed) => {
      slots.clear();
      stashed.forEach(keep);
    },
    snapshot: kept,
    restore: (snapshot) => {
      slots.clear();
      snapshot.forEach(keep);
    },
    takeSlot: (drawingId) => {
      const taken = slots.get(drawingId) ?? null;
      slots.delete(drawingId);
      return taken;
    },
    putSlot: keep,
    writeValues: (drawing, values) => keep(withArtifactSlotValues(slots.get(drawing.id) ?? null, drawing, values))
  };
};
