// @ts-check
// The drawing a state-free service takes (R13: such a service never imports
// state.js). state.js's `DrawingState` must satisfy it: the History service
// is generic over the drawing its ports take, so the composition root's
// ports check that. The members a slice owns are that slice's read-only
// view, so the type checker rejects a write that does not go through the
// slice's writers; `slices` holds the handles the writers take. A setting no
// slice owns yet has the type its readers rely on, or `unknown` (state.js
// builds the settings with Vue's untyped `reactive`).
/** @import { Ref } from '../utils/reactivity.js' */
/** @import { ModeProfileName } from '../mode-profiles.js' */
/** @import { LayoutPreferences } from './layout-preferences.js' */
/** @import { LinearComparisonPlan } from './linear-comparisons.js' */
/** @import { AnnotationSetsSlice, AnnotationSetsView } from './drawing/annotation-sets.js' */
/** @import { CanvasSlice, CanvasView } from './drawing/canvas.js' */
/** @import { ColorRulesSlice, ColorRulesView } from './drawing/color-rules.js' */
/** @import { FeatureEditsSlice, FeatureEditsView } from './drawing/feature-edits.js' */
/** @import { GroupNamesSlice, GroupNamesView } from './drawing/group-names.js' */
/** @import { InputType, InputsSlice, InputsView } from './drawing/inputs.js' */
/** @import { LabelTablesSlice, LabelTablesView } from './drawing/label-tables.js' */
/** @import { PaintSlice, PaintView } from './drawing/paint.js' */
/** @import { PaletteSlice, PaletteView } from './drawing/palette.js' */
/** @import { SettingsSlice } from './drawing/settings.js' */
/** @import { VisibilityRulesSlice, VisibilityRulesView } from './drawing/visibility-rules.js' */

/**
 * The identity of a drawing, its key in `state.drawings`. A drawing's mode is
 * its `mode` member: code never reads a mode from an ID. Until a project holds
 * several drawings of one mode, each mode has one drawing, whose ID is the
 * mode's name, so a writer that knows only a mode (a mode button, a Session's
 * `ui.mode`, a Result's request) names the drawing of that mode through
 * `drawingIdOfMode`, the one constructor of an ID. The brand keeps a mode from
 * standing in for an ID unnoticed.
 * @typedef {string & { readonly __drawingId: true }} DrawingId
 */
/**
 * A drawing's identity members, constant: its ID and its diagram mode.
 * @typedef {{ readonly id: DrawingId, readonly mode: ModeProfileName }} DrawingIdentity
 */

/**
 * The handles of a drawing's slices.
 * @typedef {{
 *   readonly visibilityRules: VisibilityRulesSlice,
 *   readonly colorRules: ColorRulesSlice,
 *   readonly palette: PaletteSlice,
 *   readonly labelTables: LabelTablesSlice,
 *   readonly annotationSets: AnnotationSetsSlice,
 *   readonly groupNames: GroupNamesSlice,
 *   readonly canvas: CanvasSlice,
 *   readonly featureEdits: FeatureEditsSlice,
 *   readonly paint: PaintSlice,
 *   readonly settings: SettingsSlice,
 *   readonly inputs: InputsSlice
 * }} DrawingSlices
 */

/**
 * The settings members of a drawing that typed modules read; the others are
 * not typed here yet. The settings slice (drawing/settings.js) installs,
 * captures, and resets all of its members; their controls write them
 * (plan-r2 §13 item 3).
 * @typedef {object} DrawingSettings
 * @property {unknown} form
 * @property {unknown} adv
 * @property {unknown} losat
 * @property {Ref<string>} losatProgram
 * @property {LinearComparisonPlan} linearComparisonPlan
 * @property {LayoutPreferences} layoutPreferences
 * @property {Ref<boolean>} linearTypographyLinked
 */

/**
 * @typedef {DrawingIdentity & VisibilityRulesView & ColorRulesView & PaletteView & LabelTablesView & AnnotationSetsView
 *   & GroupNamesView & CanvasView & FeatureEditsView & PaintView & DrawingSettings & InputsView
 *   & { readonly slices: DrawingSlices }} DrawingView
 */

/**
 * The ID of the drawing of `mode`: each mode has one drawing, named by the
 * mode. The one cast to a `DrawingId`.
 * @param {ModeProfileName} mode
 * @returns {DrawingId}
 */
export const drawingIdOfMode = (mode) => /** @type {DrawingId} */ (mode);

/**
 * The ID a new drawing of `mode` gets beside `takenIds`: the mode when free,
 * else `<mode>-<n>` with the smallest free n >= 2 (B2-D4). Python's
 * `gbdraw.session_drawings.allocate_drawing_id` takes the same rule; both read
 * the `ids` vectors of `tests/fixtures/sessions/drawings-vectors.json`.
 * @param {ModeProfileName} mode
 * @param {Iterable<string>} takenIds
 * @returns {DrawingId}
 */
export const allocateDrawingId = (mode, takenIds) => {
  const taken = new Set(takenIds);
  let id = String(mode);
  for (let number = 2; taken.has(id); number += 1) id = `${mode}-${number}`;
  return /** @type {DrawingId} */ (id);
};

// The name of each mode's first drawing (Python's `DRAWING_NAMES`).
/** @type {Readonly<Record<ModeProfileName, string>>} */
const MODE_DRAWING_NAMES = Object.freeze({ circular: 'Circular', linear: 'Linear' });

/**
 * The name a new drawing gets: its mode's name, followed by its ID's number
 * for a second drawing of the mode (`linear-2` is "Linear 2"). Python's
 * `gbdraw.session_drawings.allocate_session_drawing` takes the same rule; both
 * read the `names.defaults` vectors of
 * `tests/fixtures/sessions/drawings-vectors.json`.
 * @param {ModeProfileName} mode
 * @param {DrawingId} drawingId
 * @returns {string}
 */
export const defaultDrawingName = (mode, drawingId) => {
  const number = String(drawingId).slice(String(mode).length + 1);
  return MODE_DRAWING_NAMES[mode] + (number ? ` ${number}` : '');
};

/**
 * The first drawing of a mode, its ID and tab name: a Session 27-44's
 * drawings until the Worker names them (B2d-4), and the unused drawing Load
 * keeps of a mode the Session has none of (D-15-93).
 * @param {ModeProfileName} mode
 * @returns {Readonly<DrawingIdentity & { name: string }>}
 */
export const firstDrawingOfMode = (mode) => {
  const id = drawingIdOfMode(mode);
  return Object.freeze({ id, mode, name: defaultDrawingName(mode, id) });
};

/**
 * Why `names` cannot name the drawings of a project, if they cannot: a blank
 * name, or two names that differ only in case. Python's
 * `gbdraw.session_drawings._identity_problem` takes the same rule (case
 * folding as `str.casefold`); both read the `names.identities` vectors of
 * `tests/fixtures/sessions/drawings-vectors.json`.
 * @param {readonly string[]} names
 * @returns {'BLANK_NAME' | 'REPEATED_NAME' | null}
 */
export const drawingNamesProblem = (names) => {
  if (names.some((name) => name.trim() === '')) return 'BLANK_NAME';
  const folded = names.map((name) => name.toUpperCase().toLowerCase());
  return new Set(folded).size === folded.length ? null : 'REPEATED_NAME';
};

/**
 * The Session 48 and History `ui` key of a drawing's input type: the field of
 * its mode. A History restore of one drawing writes only its own (OIC-030).
 * @param {DrawingIdentity} drawing
 * @returns {'cInputType' | 'lInputType'}
 */
export const inputTypeKeyOf = (drawing) => (drawing.mode === 'linear' ? 'lInputType' : 'cInputType');

/**
 * The Session 48 input type fields (`ui.cInputType`, `ui.lInputType`), each
 * the input type of the drawing of its mode: Save, the Load rollback, and
 * History write them from the drawings until the flip saves each drawing's own.
 * @param {readonly (DrawingIdentity & { readonly inputType: Ref<InputType> })[]} drawings
 * @returns {{ cInputType: InputType, lInputType: InputType }}
 */
export const sessionInputTypes = (drawings) => {
  /** @type {{ cInputType: InputType, lInputType: InputType }} */
  const types = { cInputType: 'gb', lInputType: 'gb' };
  drawings.forEach((drawing) => { types[inputTypeKeyOf(drawing)] = drawing.inputType.value; });
  return types;
};
