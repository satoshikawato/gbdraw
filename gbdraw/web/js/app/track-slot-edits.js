// @ts-check
// R10, Q3, R13: a track stack editor routes each of its edits of a feature-slot
// input through `changeTrackLayout`, the feature placement owner's transition
// that the composition root injects as a port. A DOM event the template passes
// last names the control for Cancel and focus (its listener's element, not an
// icon inside a button); the edit itself runs without the event.
/** @import { AskChoice } from './app-messages.js' */
/** @import { DepthSourceLoss } from '../services/depth-track-state.js' */
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
 * The texts of a question: its title and action, one whole sentence each, say
 * what happens to the track and to the rows that reach it.
 * @typedef {{ title: string, description?: string, action: string, keep: string }} TrackHideQuestionText
 */
/**
 * @typedef {object} TrackHide
 * @property {() => readonly SwitchedTrackRow[]} rows The draft stack's rows, read at the question and at the answer.
 * @property {'circular' | 'linear'} mode The stack's mode.
 * @property {ReadonlySet<string>} hiding The slot ids of the rows the change removes or turns off.
 * @property {ReadonlySet<string>} [alsoOff] The slot ids of rows that turn off with the rows that reach them.
 * @property {(rows: readonly string[], turnedOff: ReadonlySet<string>) => TrackHideQuestionText} question The texts,
 *   from the names and the slot ids of the rows the choice turns off.
 * @property {boolean} [asksAlone] Whether it asks also when no row reaches them (default: it does not).
 * @property {string} label The History step of the chosen change.
 * @property {() => boolean} unchanged Whether the draft is still as asked when the answer comes, beyond the stack's rows.
 * @property {(turnOffReaching: () => void) => unknown} apply The change, with the turning off of the rows that
 *   reach it, which it runs inside its own transition.
 */
/** @param {number} index A row's index in its stack. */
const trackRowName = (index) => t('Track row {row}', { row: index + 1 });
/**
 * PD-OI-104 scenario revisions 6 and 7: a change that removes or turns off
 * stack rows that other rows start at or extend to (`overlayRowsReaching`)
 * does not go ahead alone (a stack Generate refuses). It asks first: the
 * question names those rows as a refused Generate names a row ("Track row N",
 * with the row's slot id) and offers to make the change and turn them off as
 * one History step; Cancel changes nothing. Rows are found by slot id at the
 * answer, and the answer does nothing when the stack's rows, or the rows that
 * reach, are no longer those it asked about.
 * @param {{ askChoice: AskChoice, runUndoable: (label: string, fn: () => unknown) => unknown }} ports
 * @param {TrackHide} hide
 * @returns {Promise<unknown> | null} null when nothing asks: the caller makes the change itself.
 */
export const askTrackHide = ({ askChoice, runUndoable }, { rows, mode, hiding, alsoOff = new Set(), question, asksAlone = false, label, unchanged, apply }) => {
  const stack = [...rows()];
  const reaching = overlayRowsReaching(stack, mode, hiding);
  if (!asksAlone && reaching.length === 0) return null;
  const reachingIds = new Set(reaching.map((index) => String(stack[index].id)));
  const ids = new Set([...reachingIds, ...alsoOff]);
  const asked = () => {
    const now = rows();
    const reachingNow = overlayRowsReaching(now, mode, hiding);
    return now.length === stack.length && now.every((row, index) => row === stack[index])
      && reachingNow.length === reachingIds.size && reachingNow.every((index) => reachingIds.has(String(now[index].id)));
  };
  const own = question(stack.flatMap((row, index) => (ids.has(String(row.id)) ? [trackRowName(index)] : [])), ids);
  const reach = reaching.length === 0 ? '' : tn(reaching.length,
    'This row starts at or extends to a track that this change removes or turns off. It is turned off too, as it cannot be drawn without that track.',
    'These rows start at or extend to a track that this change removes or turns off. They are turned off too, as they cannot be drawn without that track.');
  return askChoice({
    title: own.title,
    description: joinSentences(...[own.description ?? '', reach].filter(Boolean)),
    items: reaching.map((index) => trackRowItem(index, stack[index].id)),
    actions: [{ id: 'hide', label: own.action, tone: 'primary' }],
    cancelLabel: own.keep
  }).then((choice) => (choice === 'hide' && asked() && unchanged()
    ? runUndoable(label, () => apply(() => rows().forEach((row) => { if (ids.has(String(row.id))) row.enabled = false; })))
    : false));
};

/**
 * What the composition root gives a stack edit that asks: `rows` is the shown drawing's stack.
 * @template {SwitchedTrackRow} Row
 * @typedef {{ askChoice: AskChoice, runUndoable: (label: string, fn: () => unknown) => unknown,
 *   busy: () => unknown, rows: () => readonly Row[], mode: 'circular' | 'linear' }} AskingStackPorts
 */

/**
 * A stack row's own switch that asks before it turns off a drawn row other
 * rows start at or extend to, and turns them off with it (`askTrackHide`).
 * The composition root wraps each stack editor's switch with it.
 * @template {SwitchedTrackRow} Row
 * @param {AskingStackPorts<Row>} ports
 * @param {(slot: Row | null | undefined, enabled: unknown, turnOffReaching?: (() => void) | Event) => unknown} setEnabled
 *   The editor's switch, run in its stack edit transition.
 * @returns {(slot: Row | null | undefined, enabled: unknown, event?: Event) => unknown}
 */
export const askingTrackRowSwitch = ({ askChoice, runUndoable, busy, rows, mode }, setEnabled) => (slot, enabled, event) => {
  const drawn = () => Boolean(slot) && slot?.enabled !== false && !busy() && rows().includes(/** @type {Row} */ (slot));
  const asked = !enabled && drawn() ? askTrackHide({ askChoice, runUndoable }, {
    rows,
    mode,
    hiding: new Set([String(slot?.id)]),
    question: (names) => {
      const row = trackRowName(rows().indexOf(/** @type {Row} */ (slot)));
      return { title: t('Turn off {row}?', { row }), action: t('Turn off {rows}', { rows: listText([row, ...names]) }), keep: t('Cancel') };
    },
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

/**
 * A stack row's Remove that asks before it removes a row other rows start at
 * or extend to, and turns them off with the removal (`askTrackHide`). The
 * composition root wraps each stack editor's Remove with it.
 * @template {SwitchedTrackRow} Row
 * @param {AskingStackPorts<Row>} ports
 * @param {(index: unknown, turnOffReaching?: (() => void) | Event) => unknown} removeRow
 *   The editor's Remove, run in its stack edit transition.
 * @returns {(index: unknown, event?: Event) => unknown}
 */
export const askingTrackRowRemoval = ({ askChoice, runUndoable, busy, rows, mode }, removeRow) => (index, event) => {
  const row = rows()[Number(index)];
  const asked = row && !busy() ? askTrackHide({ askChoice, runUndoable }, {
    rows,
    mode,
    hiding: new Set([String(row.id)]),
    question: (names) => {
      const texts = { row: trackRowName(Number(index)), rows: listText(names) };
      return { title: t('Remove {row} and turn off {rows}?', texts), action: t('Remove {row} and turn off {rows}', texts), keep: t('Cancel') };
    },
    label: msg('Remove track row'),
    unchanged: () => !busy(),
    apply: (turnOffReaching) => removeRow(rows().indexOf(row), turnOffReaching)
  }) : null;
  return asked ?? removeRow(index, event);
};

/** @param {number} index A row's index in its stack. @param {unknown} id Its slot id. */
const trackRowItem = (index, id) => t('{row} ({id})', { row: trackRowName(index), id: String(id) });

/**
 * @typedef {object} TrackSourceRemoval
 * @property {'depth-file' | 'depth-series' | 'linear-file' | 'comparison'} kind What goes: a Depth file, a whole
 *   Depth series, a Linear File with the Depth files it holds, or a comparison source.
 * @property {number} [series] The Depth series' index, of a `depth-series`.
 * @property {DepthSourceLoss} loss What it does to the stack (depth-track-state.js).
 * @property {() => readonly unknown[]} sources What the removal reads, compared by identity at the answer.
 * @property {string} label The History step of the chosen removal.
 * @property {(turnOffReaching: () => void) => unknown} apply The removal, with the turning off of the rows that reach what leaves.
 */
/**
 * The title and action of a removal's question.
 * @param {TrackSourceRemoval['kind']} kind @param {{ rows: string, series: number }} texts
 */
const trackSourceRemovalTexts = (kind, texts) => {
  if (kind === 'depth-file') return { title: t('Remove the Depth file and turn off {rows}?', texts), action: t('Remove the Depth file and turn off {rows}', texts) };
  if (kind === 'depth-series') return { title: t('Remove Depth series #{series} and turn off {rows}?', texts), action: t('Remove Depth series #{series} and turn off {rows}', texts) };
  if (kind === 'linear-file') return { title: t('Remove the File and turn off {rows}?', texts), action: t('Remove the File and turn off {rows}', texts) };
  return { title: t('Remove the comparison source and turn off {rows}?', texts), action: t('Remove the comparison source and turn off {rows}', texts) };
};

/**
 * The removal of a source that asks when a row other rows start at or extend
 * to goes with it (`askTrackHide`). A Depth row that stays in the stack with no
 * file and that another row reaches turns off with those rows, so Generate
 * draws; one that no row reaches stays on, with its row issue. It writes
 * nothing before the answer, and the answer does nothing when the sources
 * changed meanwhile.
 * @template {SwitchedTrackRow} Row
 * @param {AskingStackPorts<Row>} ports
 * @param {TrackSourceRemoval} removal
 * @returns {Promise<unknown> | null} null when nothing asks: the caller removes it itself.
 */
export const askTrackSourceRemoval = ({ askChoice, runUndoable, busy, rows, mode }, { kind, series = 0, loss, sources, label, apply }) => {
  const asked = [...sources()];
  const stack = rows();
  const alsoOff = new Set([...loss.stayingOn].filter((id) => overlayRowsReaching(stack, mode, new Set([id])).length > 0));
  return askTrackHide({ askChoice, runUndoable }, {
    rows,
    mode,
    hiding: new Set([...loss.leaving, ...alsoOff]),
    alsoOff,
    question: (names, turnedOff) => {
      const staying = stack.flatMap((row, index) => (alsoOff.has(String(row.id)) ? [trackRowItem(index, row.id)] : []));
      const emptySeries = loss.emptySeries(turnedOff);
      return {
        ...trackSourceRemovalTexts(kind, { rows: listText(names), series: series + 1 }),
        description: joinSentences(...[
          emptySeries.length === 0 ? '' : tn(emptySeries.length,
            'Depth series {series} stays empty, and Generate will ask for a file or its removal.',
            'Depth series {series} stay empty, and Generate will ask for a file or their removal.',
            { series: listText(emptySeries.map((index) => `#${index + 1}`)) }),
          staying.length === 0 ? '' : tn(staying.length,
            '{rows} has no Depth file left. It stays in the stack, turned off.',
            '{rows} have no Depth file left. They stay in the stack, turned off.',
            { rows: listText(staying) })
        ].filter(Boolean)),
        keep: t('Cancel')
      };
    },
    label,
    unchanged: () => {
      const now = sources();
      return !busy() && now.length === asked.length && now.every((source, index) => source === asked[index]);
    },
    apply
  });
};
