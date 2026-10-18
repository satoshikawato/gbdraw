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
 *   readonly settings: SettingsSlice
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
 *   & GroupNamesView & CanvasView & FeatureEditsView & PaintView & DrawingSettings
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
 * The input type a drawing's slots read: the app-level ref of its mode. A
 * History restore of one drawing writes only this one (OIC-030).
 * @param {DrawingIdentity} drawing
 * @returns {'cInputType' | 'lInputType'}
 */
export const inputTypeKeyOf = (drawing) => (drawing.mode === 'linear' ? 'lInputType' : 'cInputType');
