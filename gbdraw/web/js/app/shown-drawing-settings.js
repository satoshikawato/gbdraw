// @ts-check
// What the settings and the feature popup offer for the shown drawing:
// whether the scale styling and the ruler on the axis apply, the clicked
// feature's color-rule Legend groups and the feature types a visibility rule
// suggests, whether the Label settings show, the template's names for the
// members of the shown drawing, and the default color grid's groups. The
// composition root passes the state and the owners' readers as ports.
import { linearAxisRulerApplies } from '../mode-profiles.js';
import { SETTINGS_SLICE } from '../services/drawing/settings.js';
import { labelSettingsVisible } from '../services/feature-placement.js';
import { displayedLegendRowContext } from '../services/feature-visibility.js';
import { colorRuleLegendGroups, groupLegendRowKeys } from '../services/specific-color-rules.js';
import { COMPARISON_COLOR_KEYS } from '../utils/color-utils.js';
import { toRaw } from '../utils/vue-raw.js';
/** @import { DrawingState, WebState } from '../state.js' */
/** @import { CircularTrackSlot } from '../services/circular-track-slot-model.js' */

const { computed, isRef } = window.Vue;

/**
 * @typedef {object} ScaleStyleAvailabilityPorts
 * @property {WebState} state
 * @property {(slot: CircularTrackSlot | null | undefined) => boolean} circularTrackSlotEffectiveEnabled
 *   Whether a Custom Track Slot draws.
 */

/**
 * Registers the `canUseLinearRulerOnAxis` and `canUseCircularScaleStyling`
 * computeds; the root calls it where they stood.
 * @param {ScaleStyleAvailabilityPorts} ports
 */
export const createScaleStyleAvailability = ({ state, circularTrackSlotEffectiveEnabled }) => {
  const canUseLinearRulerOnAxis = computed(
    () => linearAxisRulerApplies(state.activeDrawing().form)
  );

  const canUseCircularScaleStyling = computed(
    () => {
      const drawing = state.activeDrawing();
      return (
        drawing.adv.circular_track_slots_enabled
          ? drawing.adv.circular_track_slots.some((/** @type {CircularTrackSlot} */ slot) => (
              slot?.renderer === 'ticks' &&
              circularTrackSlotEffectiveEnabled(slot)
            ))
          : drawing.form.show_scale !== false
      );
    }
  );
  return { canUseLinearRulerOnAxis, canUseCircularScaleStyling };
};

/**
 * @typedef {object} FeatureChoiceReaderPorts
 * @property {WebState} state
 * @property {(caption: unknown) => string} shownLegendName The name the Legend shows for a row.
 */

/**
 * Registers the `popupLegendGroups` and `featureVisibilityFeatureSuggestions`
 * computeds; the root calls it where they stood.
 * @param {FeatureChoiceReaderPorts} ports
 */
export const createFeatureChoiceReaders = ({ state, shownLegendName }) => {
  const { clickedFeature, featureKeys, extractedFeatures } = state;
  // The list under the popup's Legend name: the color-rule Legend groups of
  // the clicked feature's type, by name (PD-OI-061 revision 4). Each is shown
  // as the Legend shows its row: by the row's drawn name, and a group that
  // draws no row by its name as written, never by a key.
  const popupLegendGroups = computed(() => {
    const type = clickedFeature.value?.feat.type;
    if (type === undefined) return [];
    const drawing = state.activeDrawing();
    const rowKeyOf = groupLegendRowKeys(displayedLegendRowContext(state, drawing));
    return colorRuleLegendGroups(drawing.manualSpecificRules, String(type))
      .map((group) => ({ group, shownName: shownLegendName(rowKeyOf(group)) }))
      .sort((left, right) => left.group.name.localeCompare(right.group.name) || left.shownName.localeCompare(right.shownName));
  });

  const featureVisibilityFeatureSuggestions = computed(() => {
    const values = new Set(['*', ...featureKeys]);
    for (const feat of extractedFeatures.value) {
      const type = String(feat?.type || '').trim();
      if (type) values.add(type);
    }
    return Array.from(values);
  });
  return { popupLegendGroups, featureVisibilityFeatureSuggestions };
};

/**
 * @typedef {object} ShownDrawingMemberPorts
 * @property {WebState} state
 * @property {{ readonly value: string }} linearComparisonGlobalAction The comparison every Linear pair runs.
 * @property {(drawing: DrawingState, colors?: Record<string, string | null>) => [string, string][]} userDefaultColorRows
 *   The user's own default colors of the drawing, as key and color.
 */

/**
 * Registers the `labelSettingsShown` and `defaultColorGroups` computeds; the
 * root calls it where they stood. Each call of `drawingMember` creates one
 * computed, so the root creates this before its return calls it.
 * @param {ShownDrawingMemberPorts} ports
 */
export const createShownDrawingMembers = ({ state, linearComparisonGlobalAction, userDefaultColorRows }) => {
  const { mode } = state;
  const labelSettingsShown = computed(() => {
    const drawing = state.activeDrawing();
    const scope = mode.value === 'linear' ? drawing.form.show_labels_linear : drawing.form.labels_mode;
    return labelSettingsVisible(scope, drawing.featureOverrides);
  });

  // The template's names for the members of the shown mode's drawing: a ref
  // member is its value, any other member the object itself. Only a setting's
  // ref takes a write, its control's (plan-r2 §13 item 3); an editor member's
  // writers are ports, and a derived value is read-only.
  /** @param {keyof DrawingState} key */
  const drawingMember = (key) => {
    if (!isRef(state.activeDrawing()[key])) return computed(() => state.activeDrawing()[key]);
    if (!(/** @type {readonly string[]} */ (SETTINGS_SLICE.members)).includes(key)) return computed(() => state.activeDrawing()[key].value);
    return computed({
      get: () => state.activeDrawing()[key].value,
      /** @param {unknown} value */
      set: (value) => SETTINGS_SLICE.install(state.activeDrawing().slices.settings, { [key]: value })
    });
  };

  // The default colors of the feature types the shown figure draws,
  // and the user's own default colors, come first; the others wait under More.
  // The values are read untracked, so editing a color never moves it to the
  // other grid; the grouping follows the drawn types, the palette, the keys,
  // and a replaced color table (Reset, palette change, Session load, an
  // imported Default colors file).
  const defaultColorGroups = computed(() => {
    const drawing = state.activeDrawing();
    /** @type {Record<string, string>} */
    const colors = drawing.currentColors.value;
    const keys = Object.keys(colors);
    const drawn = new Set(state.extractedFeatures.value
      .map((/** @type {{ type?: string }} */ feature) => String(feature?.type || '')));
    if (!drawn.size) return { primary: keys, more: /** @type {string[]} */ ([]) };
    const circular = mode.value === 'circular';
    if (circular ? !drawing.form.suppress_gc : drawing.form.show_gc) drawn.add('gc_content');
    if (circular ? !drawing.form.suppress_skew : drawing.form.show_skew) ['skew_high', 'skew_low'].forEach((key) => drawn.add(key));
    if (!circular && drawing.linearSeqs.length > 1 && linearComparisonGlobalAction.value !== 'none') {
      COMPARISON_COLOR_KEYS.forEach((key) => drawn.add(key));
    }
    userDefaultColorRows(drawing, toRaw(colors)).forEach(([key]) => drawn.add(key));
    const primary = keys.filter((key) => key === 'default' || drawn.has(key));
    return { primary, more: keys.filter((key) => !drawn.has(key) && key !== 'default') };
  });
  return { labelSettingsShown, drawingMember, defaultColorGroups };
};
