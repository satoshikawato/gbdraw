// The live setting rows the random walk's `setting` step drives
// (tests/web/helpers/live-random-walk.mjs): the rows of `SETTING_APPLICATION`
// (gbdraw/web/js/services/setting-application.js) that the automatic
// rerender applies from the draft (R1 (c)): `application` `live` with
// `requestPaths`. The walk reads them from the module at run time, so a new
// live row joins the walk without a list here. A pure module: the page imports
// it (`/tests/web/helpers/live-setting-rows.mjs`) to find the controls, and the
// node test tests/web/live-random-walk-settings.test.mjs checks in index.html
// that every row has a control the step can drive.
//
// A control belongs to the row of its nearest `data-setting` scope (the rule of
// tests/web/setting-application.test.mjs). The step drives a control of a type
// it can set (below) that records a History step: not one inside
// `data-history-ignore`, `data-history-managed`, or `data-history-scope="transient"`
// (app/history-inputs.js skips those).

/**
 * A row of `SETTING_APPLICATION` as this module reads it.
 * @typedef {{ readonly id: string, readonly application: string, readonly requestPaths?: Readonly<Record<string, string | readonly string[]>> }} SettingRowLike
 */

/**
 * The rows the automatic rerender applies from the draft.
 * @template {SettingRowLike} T
 * @param {readonly T[]} rows `SETTING_APPLICATION`
 * @returns {T[]}
 */
export const liveRequestRows = (rows) => rows.filter((row) => row.application === 'live' && Boolean(row.requestPaths));

// The live rows whose draft only an action writes, not a control: the walk
// kind that runs that action. A Sequence search's set (PD-OI-117) is the
// `sequenceSearch` kind; the search panel's own fields are transient.
export const ROWS_DRIVEN_BY_KIND = Object.freeze({ 'sequence-search': 'sequenceSearch' });

// The controls the step sets, by tag and (for an input) type, and how: `fill`
// types a value and leaves the field, `toggle` clicks it, `choose` picks
// another option.
/** @typedef {'fill' | 'toggle' | 'choose'} SettingControlAction */
/** @type {Readonly<Record<string, SettingControlAction>>} */
const INPUT_ACTIONS = Object.freeze({ number: 'fill', text: 'fill', color: 'fill', checkbox: 'toggle' });

/**
 * How the step sets a control, or null for one it does not set.
 * @param {string} tag the element's tag name, any case
 * @param {string | null | undefined} type an input's `type` (default `text`)
 * @returns {SettingControlAction | null}
 */
export const settingControlAction = (tag, type) => {
  const name = tag.toLowerCase();
  if (name === 'select') return 'choose';
  if (name === 'textarea') return 'fill';
  if (name === 'input') return INPUT_ACTIONS[String(type || 'text').toLowerCase()] ?? null;
  return null;
};

/**
 * Whether a `setting` step can set a control: a type it sets (`action`), a
 * name, and for a select another option to pick. The one requirement that
 * the page reader (`walkSettingControls`, the name `describeControl` gives)
 * and the index.html scan (tests/web/live-random-walk-settings.test.mjs, the
 * name the markup declares) both apply.
 * @param {{ action: SettingControlAction | null, name: string, otherOption: boolean }} control
 */
export const isSettableControl = ({ action, name, otherOption }) => action !== null && name !== '' && (action !== 'choose' || otherOption);
