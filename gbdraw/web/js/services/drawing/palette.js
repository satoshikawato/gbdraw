// @ts-check
// The palette slice of a drawing: its palette (`selectedPalette`, Session 46
// `modes.<m>.config.palette`), its default colors (`currentColors`,
// `config.colors`), and the palette queued for the next Generate while Instant
// Preview is off (`pendingPaletteName`, `pendingPaletteColors`, `ui`). The
// drawing shows them read-only to the type checker; only this module writes
// them, through `drawing.slices.palette`: the palette controls' writers below
// (app/results.js, app/palettes.js, the color table and preset edits), and
// Load, a History restore, and Reset (`install`, `reset`). The palette the
// shown Result was drawn with (`appliedPaletteName`) is app state, not here.
import { ref } from '../../utils/reactivity.js';
import { createSliceHandle, sliceMembers } from './slice-handle.js';
/** @import { ReadonlyRef, Ref } from '../../utils/reactivity.js' */
/** @import { SliceHandle } from './slice-handle.js' */

/** @typedef {Record<string, string>} PaletteColors A color by feature or track key. */
/**
 * @typedef {{
 *   selectedPalette: Ref<string>, currentColors: Ref<PaletteColors>,
 *   pendingPaletteName: Ref<string>, pendingPaletteColors: Ref<PaletteColors>
 * }} PaletteMembers
 */
/**
 * What `capture` returns and `install` writes.
 * @typedef {{ selectedPalette: string, currentColors: PaletteColors, pendingPaletteName: string, pendingPaletteColors: PaletteColors }} PaletteData
 */
/**
 * @typedef {{
 *   readonly selectedPalette: ReadonlyRef<string>, readonly currentColors: ReadonlyRef<Readonly<PaletteColors>>,
 *   readonly pendingPaletteName: ReadonlyRef<string>, readonly pendingPaletteColors: ReadonlyRef<Readonly<PaletteColors>>
 * }} PaletteView
 */
/** @typedef {SliceHandle<'palette'>} PaletteSlice The handle a writer takes. */
/** @typedef {PaletteView & { readonly slices: { readonly palette: PaletteSlice } }} PaletteDrawing A drawing's palette, and the handle its writers take. */

/** @returns {PaletteData} */
export const createDefaultPaletteData = () => ({ selectedPalette: 'default', currentColors: {}, pendingPaletteName: '', pendingPaletteColors: {} });

// A palette handle exists only through `create`, which gives it these members.
/** @param {PaletteSlice} slice */
const membersOf = (slice) => /** @type {PaletteMembers} */ (sliceMembers(slice, 'palette'));

/**
 * A new slice, or one over a drawing's existing members (a test fixture).
 * @param {Partial<PaletteMembers>} [members]
 * @returns {{ handle: PaletteSlice, view: PaletteView }}
 */
const create = (members = {}) => {
  const defaults = createDefaultPaletteData();
  /** @type {PaletteMembers} */
  const data = {
    selectedPalette: members.selectedPalette ?? ref(defaults.selectedPalette),
    currentColors: members.currentColors ?? ref(defaults.currentColors),
    pendingPaletteName: members.pendingPaletteName ?? ref(defaults.pendingPaletteName),
    pendingPaletteColors: members.pendingPaletteColors ?? ref(defaults.pendingPaletteColors)
  };
  return { handle: createSliceHandle('palette', data), view: Object.freeze({ ...data }) };
};

/**
 * @param {PaletteSlice} slice
 * @returns {PaletteData}
 */
const capture = (slice) => {
  const { selectedPalette, currentColors, pendingPaletteName, pendingPaletteColors } = membersOf(slice);
  return {
    selectedPalette: selectedPalette.value,
    currentColors: { ...currentColors.value },
    pendingPaletteName: pendingPaletteName.value,
    pendingPaletteColors: { ...pendingPaletteColors.value }
  };
};

// Each member `data` holds becomes a copy of it; the others are kept.
/**
 * @param {PaletteSlice} slice
 * @param {Partial<Readonly<PaletteData>>} data
 */
const install = (slice, { selectedPalette, currentColors, pendingPaletteName, pendingPaletteColors }) => {
  const members = membersOf(slice);
  if (selectedPalette !== undefined) members.selectedPalette.value = String(selectedPalette);
  if (currentColors) members.currentColors.value = { ...currentColors };
  if (pendingPaletteName !== undefined) members.pendingPaletteName.value = String(pendingPaletteName);
  if (pendingPaletteColors) members.pendingPaletteColors.value = { ...pendingPaletteColors };
};

/**
 * The queued palette as Save and a History step hold it (`ui.pendingPaletteName`,
 * `ui.pendingPaletteColors`).
 * @param {PaletteSlice} slice
 */
export const queuedPaletteData = (slice) => {
  const { pendingPaletteName, pendingPaletteColors } = capture(slice);
  return { pendingPaletteName, pendingPaletteColors };
};

/** @param {PaletteSlice} slice */
const reset = (slice) => install(slice, createDefaultPaletteData());

/**
 * The palette, and the default colors when given (a palette switch, or the
 * first palette of a drawing that has no colors yet).
 * @param {PaletteSlice} slice
 * @param {string} name
 * @param {Readonly<PaletteColors>} [colors]
 */
export const selectPalette = (slice, name, colors) => install(slice, { selectedPalette: name, currentColors: colors });

/**
 * The default colors become `colors` (Default colors Reset).
 * @param {PaletteSlice} slice
 * @param {Readonly<PaletteColors>} colors
 */
export const replaceDefaultColors = (slice, colors) => install(slice, { currentColors: colors });

/**
 * Sets the default color of each key `colors` holds (an edit, a color table,
 * a preset); the other keys keep theirs.
 * @param {PaletteSlice} slice
 * @param {Readonly<PaletteColors>} colors
 */
export const setDefaultColors = (slice, colors) => {
  const { currentColors } = membersOf(slice);
  currentColors.value = { ...currentColors.value, ...colors };
};

/**
 * Queues `name` with `colors` for the next Generate.
 * @param {PaletteSlice} slice
 * @param {string} name
 * @param {Readonly<PaletteColors>} colors
 */
export const queuePalette = (slice, name, colors) => install(slice, { pendingPaletteName: name, pendingPaletteColors: colors });

/**
 * Restores the queued palette of a generated artifact's owner set. The slice
 * keeps the colors object: the owner set holds it by reference, as it holds
 * every artifact value (History counts an artifact once per identity).
 * @param {PaletteSlice} slice
 * @param {string} name
 * @param {PaletteColors} colors
 */
export const restoreArtifactQueuedPalette = (slice, name, colors) => {
  const members = membersOf(slice);
  members.pendingPaletteName.value = name;
  members.pendingPaletteColors.value = colors;
};

/** @param {PaletteSlice} slice */
export const clearQueuedPalette = (slice) => install(slice, { pendingPaletteName: '', pendingPaletteColors: {} });

export const PALETTE_SLICE = Object.freeze({
  members: Object.freeze(['selectedPalette', 'currentColors', 'pendingPaletteName', 'pendingPaletteColors']),
  create, capture, install, reset
});
