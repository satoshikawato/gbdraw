// @ts-check
// R10, Q3, R13: a track stack editor routes each of its edits of a feature-slot
// input through `changeTrackLayout`, the feature placement owner's transition
// that the composition root injects as a port. A DOM event the template passes
// last names the control for Cancel and focus (its listener's element, not an
// icon inside a button); the edit itself runs without the event.
/** @import { AskChoice } from './app-messages.js' */
import { overlayRowsReaching } from '../services/track-slot-validation.js';
import { joinSentences, listText, msg, t, tn } from '../utils/i18n.js';

/**
 * The edited checkbox, select, or button, with the members the placement owner reads.
 * @typedef {HTMLElement & {
 *   type?: string,
 *   checked?: boolean,
 *   selectedOptions?: ArrayLike<{ text: string }>,
 *   labels?: ArrayLike<{ textContent: string }> | null
 * }} LayoutControl
 */

/**
 * R10: the feature placement owner's transition for an edit of a draft
 * feature-slot input. It runs `apply` (and may ask before keeping it) and
 * returns what `apply` returns, a busy outcome, or `false` while the user is asked.
 * @typedef {(apply: () => any, control?: LayoutControl | null) => any} ChangeTrackLayout
 */

/**
 * The Legend port a stack edit runs in: a drawn track the edit hides or shows
 * again takes its Legend rows with it (R2, OV-1198).
 * @typedef {<T>(change: () => T) => T} ChangeTrackVisibility
 */

/**
 * The port a track stack editor passes to `featureSlotEdits`: each stack edit
 * (add, renderer change, enable, disable, remove, stack switch) runs in the
 * Legend port once. The port wraps the whole layout transition, so an edit
 * whose lane check restores the inputs and asks retires nothing; the edit that
 * the Reset choice of that dialog applies later runs in the port on its own,
 * inside the Reset's History step.
 * @param {ChangeTrackLayout} changeTrackLayout
 * @param {ChangeTrackVisibility} changeTrackVisibility
 * @returns {ChangeTrackLayout}
 */
export const trackLayoutWithVisibility = (changeTrackLayout, changeTrackVisibility) => (apply, control) => {
  let inTransition = true;
  try {
    return changeTrackVisibility(() => changeTrackLayout(
      () => (inTransition ? apply() : changeTrackVisibility(apply)), control
    ));
  } finally {
    inTransition = false;
  }
};

/**
 * @template {Record<string, (...args: any[]) => any>} Actions
 * @param {ChangeTrackLayout} changeTrackLayout
 * @param {Actions} actions
 * @returns {Actions} The same names; each takes the optional trailing DOM event as well.
 */
export const featureSlotEdits = (changeTrackLayout, actions) => /** @type {Actions} */ (Object.fromEntries(
  Object.entries(actions).map(([name, action]) => [name, (...args) => {
    const event = typeof Event === 'function' && args.at(-1) instanceof Event ? args.pop() : null;
    return changeTrackLayout(() => action(...args), event?.currentTarget);
  }])
));

/**
 * A stack row as a hide reads and turns it off: its slot id and its switch.
 * @typedef {{ id?: unknown, enabled?: unknown }} SwitchedTrackRow
 */
/**
 * The question a hide asks on its own (Hide GC content?): its texts and its
 * hide and keep buttons.
 * @typedef {{ title: string, description: string, hide: string, keep: string }} TrackHideQuestion
 */
/**
 * @typedef {object} TrackHide
 * @property {() => readonly SwitchedTrackRow[]} rows The draft stack's rows, read at the question and at the answer.
 * @property {ReadonlySet<string>} hiding The slot ids of the drawn rows the hide turns off.
 * @property {TrackHideQuestion | null} question The hide's own question, or null for a row's own switch.
 * @property {string} label The History step of the chosen hide.
 * @property {() => boolean} unchanged Whether the draft is still as asked when the answer comes.
 * @property {(turnOffReaching: () => void) => unknown} apply The hide, with the turning off of the rows that
 *   reach it, which it runs inside its own transition.
 */
/**
 * PD-OI-104 scenario revision 6: a hide of drawn stack rows that other rows
 * start at or extend to (`overlayRowsReaching`) does not hide them alone (a
 * stack Generate refuses). It asks first: the question names those rows as a
 * refused Generate names a row ("Track row N", with the row's slot id) and
 * offers to hide and turn them off as one History step; Cancel changes
 * nothing. Rows are found by slot id at the answer.
 * @param {{ askChoice: AskChoice, runUndoable: (label: string, fn: () => unknown) => unknown }} ports
 * @param {TrackHide} hide
 * @returns {Promise<unknown> | null} null when nothing asks: the caller applies the hide itself.
 */
export const askTrackHide = ({ askChoice, runUndoable }, { rows, hiding, question, label, unchanged, apply }) => {
  const stack = rows();
  const reaching = overlayRowsReaching(stack, hiding);
  if (!question && reaching.length === 0) return null;
  /** @param {number} index */
  const rowName = (index) => t('Track row {row}', { row: index + 1 });
  const ids = new Set(reaching.map((index) => String(stack[index].id)));
  const switched = question ? -1 : stack.findIndex((row) => hiding.has(String(row.id)));
  const own = question ?? {
    title: t('Turn off {row}?', { row: rowName(switched) }),
    description: '',
    hide: t('Turn off {rows}', { rows: listText([switched, ...reaching].map(rowName)) }),
    keep: t('Cancel')
  };
  const reach = reaching.length === 0 ? '' : tn(reaching.length,
    'This row starts at or extends to a track being turned off. It is turned off too, as Generate cannot draw it without that track.',
    'These rows start at or extend to a track being turned off. They are turned off too, as Generate cannot draw them without that track.');
  return askChoice({
    title: own.title,
    description: joinSentences(...[own.description, reach].filter(Boolean)),
    items: reaching.map((index) => t('{row} ({id})', { row: rowName(index), id: String(stack[index].id) })),
    actions: [{
      id: 'hide',
      label: question && reaching.length > 0 ? t('{action} and turn off {rows}', { action: own.hide, rows: listText(reaching.map(rowName)) }) : own.hide,
      tone: 'primary'
    }],
    cancelLabel: own.keep
  }).then((choice) => (choice === 'hide' && unchanged()
    ? runUndoable(label, () => apply(() => rows().forEach((row) => { if (ids.has(String(row.id))) row.enabled = false; })))
    : false));
};

/**
 * A stack row's own switch that asks before it turns off a drawn row other
 * rows start at or extend to, and turns them off with it (`askTrackHide`).
 * The composition root wraps each stack editor's switch with it.
 * @template {SwitchedTrackRow} Row
 * @param {{ askChoice: AskChoice, runUndoable: (label: string, fn: () => unknown) => unknown,
 *   busy: () => unknown, rows: () => readonly Row[] }} ports `rows`: the shown drawing's stack.
 * @param {(slot: Row | null | undefined, enabled: unknown, turnOffReaching?: (() => void) | Event) => unknown} setEnabled
 *   The editor's switch, run in its stack edit transition.
 * @returns {(slot: Row | null | undefined, enabled: unknown, event?: Event) => unknown}
 */
export const askingTrackRowSwitch = ({ askChoice, runUndoable, busy, rows }, setEnabled) => (slot, enabled, event) => {
  const drawn = () => Boolean(slot) && slot?.enabled !== false && !busy() && rows().includes(/** @type {Row} */ (slot));
  const asked = !enabled && drawn() ? askTrackHide({ askChoice, runUndoable }, {
    rows,
    hiding: new Set([String(slot?.id)]),
    question: null,
    label: msg('Turn off track rows'),
    unchanged: drawn,
    apply: (turnOffReaching) => setEnabled(slot, false, turnOffReaching)
  }) : null;
  if (!asked) return setEnabled(slot, enabled, event);
  // The checkbox shows the row on until the choice.
  const box = /** @type {{ checked?: boolean } | null} */ (event?.target ?? null);
  if (box) box.checked = true;
  return asked;
};
