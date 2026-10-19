// @ts-check
/** @import { DiagramMode, DrawingState } from '../../state.js' */
/** @import { FeaturePlacementDraftRow, FeaturePlacementTarget } from '../../services/feature-placement.js' */
/** @import { ChangeTrackLayout, LayoutControl } from '../track-slot-edits.js' */
/** @import { RecordKeyedFeature } from '../../services/feature-identity.js' */
import {
  canonicalFeaturePlacements,
  featureIdentityKeyOf,
  parseFeatureIdentityKey
} from '../../services/feature-placement.js';
import { resolveCircularTrackFeaturePlacement } from '../../services/circular-track-slot-model.js';
import { createDefaultLinearTrackSlots, effectiveLinearSlotPlacement } from '../../services/linear-track-slot-model.js';
import { validateCustomTrackPlan } from '../../services/track-slot-validation.js';
import { msg, t } from '../../utils/i18n.js';

/**
 * A feature the placement control lists: a row of the Features list or a clicked feature.
 * @typedef {Record<string, unknown>} PlacedFeature
 */

/**
 * The inputs of a feature slot, as the draft of one drawing holds them.
 * @typedef {{ mode: DiagramMode, form: DrawingState['form'], adv: DrawingState['adv'] }} PlacementDraft
 */

/**
 * What `saveTrackLayout` keeps: the layout form fields and, per mode, the stack and a copy of each row.
 * @typedef {object} SavedTrackLayout
 * @property {Record<string, unknown>} form
 * @property {{ mode: string, slots: unknown, enabled: unknown, axis: unknown, rows: [Record<string, unknown>, unknown][] }[]} stacks
 */

const SIDES = { circular: ['outward', 'inward'], linear: ['above', 'below'] };
/** The shown name of a side's first lane, one whole key each. @type {Record<string, () => string>} */
const LANE_ONE_LABEL = {
  outward: () => t('Outward lane 1'),
  inward: () => t('Inward lane 1'),
  above: () => t('Above lane 1'),
  below: () => t('Below lane 1')
};

// The placement targets of a draft's feature slot: the one availability
// predicate of the popup choices and of a layout change (R3, R10). It resolves
// the draft through the same slot owner used by request projection; Result
// geometry continues to describe the last successful Generate.
/**
 * @param {PlacementDraft} draft
 * @returns {{ kind: string, side?: string, level?: number }[]} The main track and, when the slot splits, one lane per side.
 */
export const draftPlacementTargets = ({ mode, form, adv }) => {
  const trackType = mode === 'circular' ? form.track_type : form.linear_track_layout;
  const slot = adv[`${mode}_track_slots_enabled`]
    ? validateCustomTrackPlan({
      mode, slots: adv[`${mode}_track_slots`],
      axisIndex: adv[`${mode}_track_slots_axis_index`], trackType
    }).enabledSlots.find((entry) => entry.renderer === 'features')
    : (mode === 'circular' ? {} : createDefaultLinearTrackSlots({ trackLayout: trackType })[0]);
  if (!slot) return [];
  const bidirectional = mode === 'circular'
    ? resolveCircularTrackFeaturePlacement(slot, trackType).laneDirection === 'split'
    : effectiveLinearSlotPlacement(slot) === 'overlay' && !form.separate_strands;
  return [{ kind: 'main' }, ...(bidirectional
    ? SIDES[mode].map((side) => ({ kind: 'lane', side, level: 1 })) : [])];
};

// The draft inputs of draftPlacementTargets, saved and restored in place so
// the stack rows keep their identity (R10).
const LAYOUT_FORM_FIELDS = ['track_type', 'linear_track_layout', 'separate_strands'];
/**
 * @param {unknown} value
 * @returns {unknown}
 */
const copy = (value) => (Array.isArray(value) ? value.map(copy) : value && typeof value === 'object'
  ? Object.fromEntries(Object.entries(value).map(([key, entry]) => [key, copy(entry)])) : value);
/**
 * @param {Record<string, any>} draft The draft's `form` and `adv`.
 * @returns {SavedTrackLayout}
 */
export const saveTrackLayout = ({ form, adv }) => ({
  form: Object.fromEntries(LAYOUT_FORM_FIELDS.map((field) => [field, form[field]])),
  stacks: Object.keys(SIDES).map((mode) => {
    const slots = adv[`${mode}_track_slots`];
    return { mode, slots, enabled: adv[`${mode}_track_slots_enabled`], axis: adv[`${mode}_track_slots_axis_index`],
      rows: /** @type {[Record<string, unknown>, unknown][]} */ ((Array.isArray(slots) ? slots : []).map((slot) => [slot, copy(slot)])) };
  })
});
/**
 * @param {Record<string, any>} draft The draft's `form` and `adv`.
 * @param {SavedTrackLayout} saved The result of `saveTrackLayout`.
 */
export const restoreTrackLayout = ({ form, adv }, saved) => {
  Object.assign(form, saved.form);
  for (const { mode, slots, enabled, axis, rows } of saved.stacks) {
    Object.assign(adv, { [`${mode}_track_slots_enabled`]: enabled, [`${mode}_track_slots_axis_index`]: axis,
      [`${mode}_track_slots`]: slots });
    for (const [slot, content] of rows) {
      Object.keys(slot).forEach((key) => delete slot[key]);
      Object.assign(slot, content);
    }
    if (Array.isArray(slots)) slots.splice(0, slots.length, ...rows.map(([slot]) => slot));
  }
};

/**
 * @typedef {object} FeaturePlacementActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoable History's undoable step.
 * @property {() => ({ mode?: string } | null)} getCommittedRequest The committed canonical request.
 * @property {(feature: RecordKeyedFeature) => boolean} isCurrentFeature Whether the feature belongs to the displayed Result.
 * @property {<T extends object>(value: T) => T} [reactive] Vue `reactive`
 * @property {() => Promise<any>} [nextTick] Vue `nextTick`
 */

/** @param {FeaturePlacementActionsOptions} options */
export const createFeaturePlacementActions = ({
  state, runUndoable, getCommittedRequest, isCurrentFeature,
  reactive = (value) => value, nextTick = () => Promise.resolve()
}) => {
  /** @param {DrawingState} drawing */
  const draft = (drawing) => ({ mode: state.mode.value, form: drawing.form, adv: drawing.adv });
  /**
   * @param {DrawingState} drawing
   * @param {PlacedFeature[]} features
   */
  const targetsFor = (drawing, features) => {
    const mode = state.mode.value;
    if (!features.length || getCommittedRequest()?.mode !== mode) return [];
    /** @type {{ recordKeys: unknown[] }[]} */
    const items = state.featureCatalog.value?.items || [];
    if (!features.every((feature) => items.some((item) => item.recordKeys.includes(feature?.record_key)))) return [];
    return draftPlacementTargets(draft(drawing));
  };
  /** @param {PlacedFeature[]} features */
  const choices = (features) => {
    const drawing = state.activeDrawing();
    const targets = targetsFor(drawing, features);
    return ['auto', 'main', ...SIDES[/** @type {DiagramMode} */ (state.mode.value)]].map((value) => {
      const enabled = features.length > 0 && features.every((feature) => {
        if (!featureIdentityKeyOf(feature) || !isCurrentFeature(feature)) return false;
        if (value === 'auto') return true;
        return targets.some((target) => value === 'main' ? target.kind === 'main' : target.side === value);
      });
      return { value, enabled, label: value === 'auto' ? t('Auto') : value === 'main' ? t('Main') : LANE_ONE_LABEL[value](),
      reason: enabled ? '' : t('Unavailable in the current draft feature slot.') };
    });
  };
  /**
   * @param {PlacedFeature[]} features
   * @param {string} value
   */
  const setPlacement = (features, value) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    const choice = choices(features).find((entry) => entry.value === value);
    if (!choice?.enabled) throw new Error(choice?.reason || 'Unknown feature placement.');
    return runUndoable(features.length === 1 ? msg('Change feature placement') : msg('Change selected feature placements'), () => {
      for (const feature of features) {
        const key = featureIdentityKeyOf(feature);
        if (value === 'auto') delete drawing.featurePlacementOverrides[key];
        else {
          // `key` comes from featureIdentityKeyOf, so parseFeatureIdentityKey reads it back; canonicalFeaturePlacements validates the row.
          const row = /** @type {FeaturePlacementDraftRow} */ ({ ...parseFeatureIdentityKey(key), placement: /** @type {FeaturePlacementTarget} */ (value === 'main'
            ? { kind: 'main' } : { kind: 'lane', side: value, level: 1 }) });
          canonicalFeaturePlacements({ [key]: row }, state.mode.value);
          drawing.featurePlacementOverrides[key] = row;
        }
      }
    });
  };

  // Q3 (Owner, 2026-10-04) and R10: the one transition for an edit of a draft
  // feature-slot input. It applies the edit and, when a lane placement of the
  // drawing loses its lane, restores the inputs and asks. A drawing's
  // placements are of its own mode (PD-OI-086). Reset applies the
  // edit and removes those rows as one History step; Cancel keeps the control's
  // value and records none. The control's own History adapter records an edit
  // that loses nothing (R11). Restores (Undo/Redo, Session load, Import, Reset
  // Settings) install state as is; Generate names what they leave (R6).
  // The composition root passes `changeTrackLayout` to the track stack editors
  // as their one port (R13); each routes its feature-slot edits through it.
  const layoutChange = reactive({ open: false, count: 0, setting: '', value: '' });
  /** @type {{ control: LayoutControl | null, apply: () => any } | null} */
  let pendingLayout = null;
  /** @param {DrawingState} drawing */
  // `filter(Boolean)` drops the undefined side of the main target.
  const laneSides = (drawing) => /** @type {string[]} */ (draftPlacementTargets({ mode: state.mode.value, form: drawing.form, adv: drawing.adv })
    .map((target) => target.side).filter(Boolean));
  /**
   * @param {DrawingState} drawing
   * @param {string[]} before
   */
  const lostLaneKeys = (drawing, before) => {
    const after = laneSides(drawing);
    return Object.keys(drawing.featurePlacementOverrides).filter((key) => {
      const { placement } = drawing.featurePlacementOverrides[key] || {};
      return before.includes(placement?.side) && !after.includes(placement?.side);
    });
  };
  /** @param {() => HTMLElement | null | undefined} find */
  const focusAfterRender = (find) => nextTick().then(() => find()?.focus?.());
  /**
   * The `ChangeTrackLayout` port of the track stack editors (app/track-slot-edits.js).
   * @type {ChangeTrackLayout}
   */
  const changeTrackLayout = (apply, control = null) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!Object.values(drawing.featurePlacementOverrides).some((row) => row?.placement?.kind === 'lane')) return apply();
    const before = laneSides(drawing);
    const saved = saveTrackLayout(drawing);
    const result = apply();
    const lost = lostLaneKeys(drawing, before);
    if (!lost.length) return result;
    restoreTrackLayout(drawing, saved);
    const checkbox = control?.type === 'checkbox';
    const value = checkbox ? (control.checked ? 'On' : 'Off') : control?.selectedOptions?.[0]?.text || '';
    // The checkbox shows the kept value again; a bound select re-renders to it.
    if (checkbox) control.checked = !control.checked;
    pendingLayout = { control, apply };
    Object.assign(layoutChange, { open: true, count: lost.length, value,
      setting: control?.getAttribute?.('aria-label') || control?.labels?.[0]?.textContent.trim()
        || (control?.tagName === 'BUTTON' && control.textContent.trim()) || control?.title || t('This change') });
    // The primary choice of the layout dialog is a button.
    void focusAfterRender(() => /** @type {HTMLElement | null | undefined} */ (globalThis.document?.querySelector('[data-placement-layout-primary]')));
    return false;
  };
  /** @param {'reset' | 'cancel'} choice */
  const resolveLayoutChange = async (choice) => {
    const drawing = state.activeDrawing();
    const change = pendingLayout;
    pendingLayout = null;
    layoutChange.open = false;
    if (!change) return false;
    const result = choice === 'reset' && await runUndoable(msg('Change setting and reset Feature placements'), () => {
      const before = laneSides(drawing);
      change.apply();
      for (const key of lostLaneKeys(drawing, before)) delete drawing.featurePlacementOverrides[key];
    });
    // A Row actions item closed its menu when it ran: focus returns to the menu's toggle.
    void focusAfterRender(() => /** @type {HTMLElement | null | undefined} */ (
      change.control?.closest?.('details:not([open])')?.querySelector(':scope > summary')
    ) ?? change.control);
    return result;
  };
  /**
   * @param {Event} event A change of the select or checkbox that holds the setting.
   * @param {string} field
   */
  const changeLayoutSetting = (event, field) => {
    const drawing = state.activeDrawing();
    // The listeners sit on a select or a checkbox.
    const control = /** @type {LayoutControl & { value?: string }} */ (event.target);
    const value = control.type === 'checkbox' ? control.checked : control.value;
    return changeTrackLayout(() => { drawing.form[field] = value; }, control);
  };

  return { choices, setPlacement, layoutChange, resolveLayoutChange, changeLayoutSetting, changeTrackLayout,
    valueFor: (/** @type {PlacedFeature} */ feature) => {
      const drawing = state.activeDrawing();
      // The control lists the shown mode's drawing's placements (R2).
      const target = drawing.featurePlacementOverrides[featureIdentityKeyOf(feature)]?.placement;
      return target?.kind === 'main' ? 'main' : target?.side || 'auto';
    } };
};
