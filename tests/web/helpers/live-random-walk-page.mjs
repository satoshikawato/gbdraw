// The identity by which a random-walk step names a feature
// (tests/web/helpers/live-random-walk.mjs). The walk records it and a replay
// finds the feature by it after the fixture has booted again, so it must not
// change between boots: a feature's `id` holds its record key, which for a
// Linear record is a UUID made at each boot (OV-2236). The record's position
// in the displayed Result (`record_idx`) and the biological feature ID (the
// key Python gives the feature in its record) do not change. The page imports
// this module (`/tests/web/helpers/live-random-walk-page.mjs`), so every
// reader of the identity uses this one derivation. It also finds the controls
// a `setting` step sets (`walkSettingControls`) and the value it gives one.
import { ROWS_DRIVEN_BY_KIND, isSettableControl, liveRequestRows, settingControlAction } from './live-setting-rows.mjs';

/**
 * A feature row of the displayed Result as the identity reads it.
 * @typedef {{ record_idx: unknown, biological_feature_id: unknown }} WalkFeatureRow
 */

/** @param {WalkFeatureRow} feature @returns {string} */
export const walkFeatureKey = (feature) => {
  const recordIndex = feature.record_idx;
  const featureId = String(feature.biological_feature_id ?? '');
  if (!Number.isInteger(recordIndex) || !featureId) {
    throw new Error(`a feature row without a record position or biological feature ID: ${JSON.stringify({ recordIndex, featureId })}`);
  }
  return `record ${recordIndex}, feature ${featureId}`;
};

/**
 * The feature of `features` that `key` names, or undefined.
 * @template {WalkFeatureRow} T
 * @param {readonly T[]} features
 * @param {string} key
 * @returns {T | undefined}
 */
export const findWalkFeature = (features, key) => features.find((feature) => walkFeatureKey(feature) === key);

/** @import { SettingControlAction, SettingRowLike } from './live-setting-rows.mjs' */
/**
 * A control a `setting` step of the walk can set now: the live row it belongs
 * to (its nearest `data-setting` scope), its name as the History list shows it
 * (`describeControl`; a second control of one name in a row is `<name> #2`),
 * and how the step sets it.
 * @typedef {{ row: string, control: string, action: SettingControlAction, element: HTMLElement }} WalkSettingControl
 */
/**
 * Whether the page shows `element`: no ancestor is hidden. A closed `<details>`
 * hides nothing here: the step opens it, as a user does.
 * @param {Element} element
 */
const displayed = (element) => {
  for (let node = /** @type {Element | null} */ (element); node; node = node.parentElement) {
    const style = getComputedStyle(node);
    if (style.display === 'none' || style.visibility === 'hidden') return false;
  }
  return true;
};
/**
 * The controls a `setting` step can set now, in document order: of the live
 * rows the automatic rerender applies (`liveRequestRows`) that a control
 * drives, each enabled, shown control of a type the step sets that records a
 * History step and has a name. A select needs another option to pick.
 * The rows are `SETTING_APPLICATION`, read from its module now; whether a
 * control records a History step is `isIgnoredTarget` of app/history-inputs.js,
 * and its name the one `describeControl` gives.
 * @returns {Promise<WalkSettingControl[]>}
 */
export const walkSettingControls = async () => {
  const [{ SETTING_APPLICATION }, { describeControl }, { isIgnoredTarget }] = await Promise.all([
    import('/gbdraw/web/js/services/setting-application.js'),
    import('/gbdraw/web/js/utils/control-description.js'),
    import('/gbdraw/web/js/app/history-inputs.js')
  ]);
  /** @type {readonly SettingRowLike[]} */
  const rows = SETTING_APPLICATION;
  const isIgnored = (/** @type {Element} */ element) => isIgnoredTarget(element);
  const nameOf = (/** @type {Element} */ element) => describeControl(element)?.name ?? '';
  const ids = new Set(liveRequestRows(rows).map(({ id }) => id).filter((id) => !(id in ROWS_DRIVEN_BY_KIND)));
  /** @type {WalkSettingControl[]} */
  const controls = [];
  for (const element of /** @type {HTMLElement[]} */ ([...document.querySelectorAll('[data-setting] input, [data-setting] select, [data-setting] textarea')])) {
    const row = element.closest('[data-setting]')?.getAttribute('data-setting') ?? '';
    const action = settingControlAction(element.tagName, element.getAttribute('type'));
    if (!ids.has(row) || element.matches(':disabled') || isIgnored(element) || !displayed(element)) continue;
    const name = nameOf(element);
    const otherOption = element.tagName === 'SELECT'
      && [.../** @type {HTMLSelectElement} */ (element).options].some((option) => !option.disabled && !option.selected);
    if (!isSettableControl({ action, name, otherOption }) || !action) continue;
    const same = controls.filter((control) => control.row === row && control.control.replace(/ #\d+$/, '') === name).length;
    controls.push({ row, control: same ? `${name} #${same + 1}` : name, action, element });
  }
  return controls;
};

const WALK_TEXTS = Object.freeze(['Walk setting text', 'Second walk text']);
const WALK_COLORS = Object.freeze(['#c83366', '#2a9d8f', '#264653']);
/**
 * The value a `setting` step gives `element`, drawn by `pick` from values the
 * control does not hold: a number scaled from the one it holds (10 for an empty
 * Auto field) within its `min`, `max` and `step`; a walk text; a walk color;
 * another enabled option of a select. Null for a checkbox, which the step clicks.
 * @param {HTMLElement} element
 * @param {SettingControlAction} action
 * @param {number} pick
 * @returns {string | null}
 */
export const walkSettingValue = (element, action, pick) => {
  const field = /** @type {HTMLInputElement} */ (element);
  if (action === 'toggle') return null;
  const other = (/** @type {readonly string[]} */ values) => {
    const left = values.filter((value) => value !== field.value);
    return left[pick % left.length];
  };
  if (action === 'choose') return other([.../** @type {HTMLSelectElement} */ (element).options].filter((option) => !option.disabled).map((option) => option.value));
  if (field.type === 'color') return other(WALK_COLORS);
  if (field.type !== 'number') return other(WALK_TEXTS);
  const step = Number(field.step) > 0 ? Number(field.step) : 1;
  const min = field.min === '' ? -Infinity : Number(field.min);
  const max = field.max === '' ? Infinity : Number(field.max);
  const held = field.value === '' ? Number.NaN : Number(field.value);
  const base = held > 0 ? held : 10;
  const decimals = (String(step).split('.')[1] || '').length;
  const fit = (/** @type {number} */ value) => Number((Math.min(max, Math.max(min, Math.round(value / step) * step))).toFixed(decimals));
  const candidates = [1.5, 0.75, 2].map((factor) => fit(base * factor)).filter((value) => value !== held);
  return String(candidates.length ? candidates[pick % candidates.length] : fit(base + step));
};
