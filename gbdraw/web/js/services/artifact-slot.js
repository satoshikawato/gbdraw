// @ts-check
// E1: one diagram mode's generated artifact, which the composition root keeps
// while the other mode is shown. The History snapshot service captures and
// installs it (`captureArtifactSlot`, `installArtifactSlot`); Session Save
// writes it and Session Load builds it (`services/config.js`).

/**
 * E1: the generated artifact of one diagram mode, which each mode keeps while
 * the other is shown. The state ref of every key and the value an empty slot
 * installs, given the slot's mode and the current value: the displayed-Result
 * part of the owner set, the selected Result, the generated Legend colors and
 * stroke defaults, and the alignment Reset receipt (it binds the
 * committed Linear Session). Draft and project state stays shared: the specific
 * color rules, the file Legend captions, the Linear record orientations, the
 * pending palette, the LOSAT caches and manifests, the alignment plan and
 * record translations, and every editor intent (R2).
 * @type {Readonly<Record<string, (mode: 'circular' | 'linear', current: any) => any>>}
 */
export const ARTIFACT_SLOT_EMPTY = Object.freeze({
  results: () => [],
  selectedResultIndex: () => 0,
  featureCatalog: () => null,
  extractedFeatures: () => [],
  biologicalFeatures: () => [],
  featureRecordIds: () => [],
  orthogroups: () => [],
  featureOrthogroupIndex: () => new Map(),
  collinearGroups: () => [],
  trackSlotResolvedGeometry: () => null,
  annotationWarnings: () => [],
  featureIdentityNotices: () => [],
  comparisonWarnings: () => [],
  lastRunInfo: () => null,
  pairwiseMatchFactors: () => ({}),
  editableLabels: () => [],
  generatedLegendPosition: () => 'left',
  generatedMode: (mode) => mode,
  generatedMultiRecordCanvas: () => false,
  generatedCircularPlotTitlePosition: () => 'none',
  // An empty slot has no Result to project a palette onto: it keeps the applied one.
  appliedPaletteName: (_mode, current) => current,
  appliedPaletteColors: (_mode, current) => current,
  similarityAlignmentResetReceipt: () => null,
  originalLegendColors: () => ({}),
  originalSvgStroke: () => ({ color: null, width: null })
});
export const ARTIFACT_SLOT_KEYS = Object.freeze(Object.keys(ARTIFACT_SLOT_EMPTY));

/**
 * One mode's generated artifact, held by reference (no clone).
 * @typedef {object} ArtifactSlot
 * @property {'circular' | 'linear'} mode
 * @property {Readonly<Record<string, any>>} values The value of every `ARTIFACT_SLOT_KEYS` key.
 * @property {readonly string[]} legendInventory The selected Result's generated Legend inventory
 *   (`originalLegendOrder`). The Legend owner keeps each Result's inventory (OV-47), so an install
 *   does not write it: the mode transition hands it to the Legend owner, which adopts it when the
 *   Result is displayed.
 * @property {readonly Record<string, any>[] | null} legendRows The Legend rows shown with the selected
 *   Result (`legendEntries`), or null when they were never shown (a Result loaded from a Session's
 *   `otherModeResult`). The Legend edits stay shared between the modes until per-drawing Legend
 *   edits: a switch rebuilds the arriving Result's rows from these and the shared edits made since.
 * @property {any} matchSequenceOwner The match-sequence registry's trusted owner, or null.
 * @property {Record<string, any> | null} runtimeState The runtime owner's state (committed Session, CLI helper files).
 * @property {Readonly<Record<string, any>> | null} transportIdentity The History transport identity of the artifact.
 * @property {number} retainedBytes
 */

/**
 * @param {{ mode: 'circular' | 'linear', values: Record<string, any>, legendInventory?: readonly string[],
 *   legendRows?: readonly Record<string, any>[] | null, matchSequenceOwner?: any,
 *   runtimeState?: Record<string, any> | null, transportIdentity?: Readonly<Record<string, any>> | null,
 *   retainedBytes?: number }} slot
 * @returns {Readonly<ArtifactSlot>}
 */
export const createArtifactSlot = ({
  mode, values, legendInventory = [], legendRows = null, matchSequenceOwner = null, runtimeState = null,
  transportIdentity = null, retainedBytes = 0
}) => {
  const unknown = Object.keys(values).filter((key) => !ARTIFACT_SLOT_KEYS.includes(key));
  if (unknown.length > 0) throw new Error(`Unknown artifact slot field(s): ${unknown.join(', ')}.`);
  const slotMode = mode === 'linear' ? 'linear' : 'circular';
  return Object.freeze({
    mode: slotMode,
    values: Object.freeze(Object.fromEntries(ARTIFACT_SLOT_KEYS.map((key) => [
      key, Object.hasOwn(values, key) ? values[key] : ARTIFACT_SLOT_EMPTY[key](slotMode, null)
    ]))),
    legendInventory: Object.freeze([...legendInventory]),
    legendRows: legendRows ? Object.freeze([...legendRows]) : null,
    matchSequenceOwner,
    runtimeState,
    transportIdentity,
    retainedBytes: Math.max(0, Number(retainedBytes) || 0)
  });
};
