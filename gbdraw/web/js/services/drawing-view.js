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
 * @typedef {VisibilityRulesView & ColorRulesView & PaletteView & LabelTablesView & AnnotationSetsView
 *   & GroupNamesView & CanvasView & FeatureEditsView & PaintView & DrawingSettings
 *   & { readonly slices: DrawingSlices }} DrawingView
 */

export {};
