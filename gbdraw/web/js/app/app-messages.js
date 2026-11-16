// @ts-check
/** @import { Wording } from '../utils/i18n.js' */
import { t, tm, tx } from '../utils/i18n.js';

// The texts of a question are shown as given: the caller words them with `t`
// or `tn` where it asks (utils/i18n.js); the choice dialog is modal, so the UI
// language cannot change while it shows. A notice stays shown while the UI
// language may change, so its caller gives how to word its text, its lines
// and its link (`Wording`), and the banner words them where it shows them: a
// language change rewords a shown notice. The status region speaks a notice
// once, in the language it was shown in. The notices of Load Session and Save
// Session, and the warnings kept as `msg` text, are worded here.

/**
 * The style of a choice button: one primary choice, secondary choices, or a
 * choice that removes work (danger).
 * @typedef {'primary' | 'secondary' | 'danger'} ChoiceTone
 *
 * @typedef {object} ChoiceAction
 * @property {string} id What the promise of `ask` resolves to.
 * @property {string} label The button text; it names the outcome.
 * @property {ChoiceTone} [tone]
 * @property {string} [detail] A line under the label that says what the choice does. When an action has one, the
 *   actions show as a list of options above Cancel.
 * @property {string} [icon] The Phosphor class of the icon that leads an option of that list, such as `ph-copy`.
 *
 * A short line of a choice dialog, led by an icon (a Phosphor class such as `ph-warning`).
 * @typedef {object} ChoiceLine
 * @property {string} icon
 * @property {string} text
 *
 * @typedef {object} ChoiceRequest
 * @property {string} title
 * @property {string} [description] One or two sentences: what the choice is about.
 * @property {ChoiceLine[]} [lines] Short icon-led lines in place of a description.
 * @property {string[]} [items] Short lines listed under the description.
 * @property {string} [note] A closing sentence under the list.
 * @property {ChoiceAction[]} actions
 * @property {string} [cancelLabel] The text of the button that changes nothing.
 *
 * Asks with the app's choice dialog; resolves to the chosen action's id, or
 * `cancel` for Cancel, Escape, a backdrop click, or a later question.
 * @typedef {(request: ChoiceRequest) => Promise<string>} AskChoice
 *
 * @typedef {object} TextField
 * @property {string} label The field's visible label, which is also its accessible name.
 * @property {string} value
 *
 * @typedef {object} TextRequest
 * @property {string} title
 * @property {string} [description]
 * @property {string} label The text field's label.
 * @property {string} [value] The text the field starts with, selected.
 * @property {string} actionLabel The button that keeps the text; Enter in the field chooses it.
 *
 * Asks for one line of text with the app's choice dialog; resolves to the
 * text as typed, or null for Cancel, Escape, a backdrop click, or a later question.
 * @typedef {(request: TextRequest) => Promise<string | null>} AskText
 *
 * @typedef {'info' | 'warning'} NoticeTone
 *
 * @typedef {object} NoticeRequest
 * @property {Wording} text
 * @property {Wording[]} [items] Warnings listed under the text.
 * @property {NoticeTone} [tone]
 *
 * A line under the notice's text with a link that opens in a new tab: `text`
 * says what the link offers, `label` names it.
 * @typedef {object} NoticeLink
 * @property {Wording} text
 * @property {Wording} label
 * @property {string} href
 *
 * The wordings of the shown notice.
 * @typedef {object} NoticeWordings
 * @property {Wording} text
 * @property {Wording[]} items
 * @property {NoticeLink | null} link
 *
 * The notice one `ShowNotice` call showed: `dismiss` ends it while it still
 * shows, and leaves a later notice.
 * @typedef {{ dismiss: () => void }} NoticeHandle
 *
 * Shows the outcome of an operation in the notice banner, replacing the
 * notice it shows. A caller without a banner returns nothing.
 * @typedef {(request: NoticeRequest) => NoticeHandle | void} ShowNotice
 */

/**
 * A confirmation asks with the app's choice dialog, never a native
 * `confirm()`, and a question for text asks with the same dialog and a text
 * field, never a native `prompt()`. One question shows at a time: a new
 * question cancels the open one.
 * @param {{ reactive: <T extends object>(value: T) => T }} options Vue's `reactive`.
 */
export const createChoicePrompt = ({ reactive }) => {
  const prompt = reactive({
    show: false,
    title: '',
    description: '',
    lines: /** @type {ChoiceLine[]} */ ([]),
    items: /** @type {string[]} */ ([]),
    note: '',
    actions: /** @type {ChoiceAction[]} */ ([]),
    cancelLabel: '',
    field: /** @type {TextField | null} */ (null)
  });
  /** @type {((choice: string) => void) | null} */
  let answer = null;
  /** @param {string} choice */
  const resolve = (choice) => {
    const settle = answer;
    answer = null;
    prompt.show = false;
    settle?.(choice);
  };
  /** @type {AskChoice} */
  const ask = (request) => {
    if (answer) resolve('cancel');
    Object.assign(prompt, {
      title: request.title,
      description: request.description || '',
      lines: [...(request.lines || [])],
      items: [...(request.items || [])],
      note: request.note || '',
      actions: request.actions.map((action) => ({ tone: 'secondary', ...action })),
      cancelLabel: request.cancelLabel || t('Cancel'),
      field: null,
      show: true
    });
    return new Promise((settle) => { answer = settle; });
  };
  /** @type {AskText} */
  const askText = async ({ title, description, label, value = '', actionLabel }) => {
    const choice = ask({ title, description, actions: [{ id: 'ok', label: actionLabel, tone: 'primary' }] });
    // The dialog's field writes the text through `prompt.field` into `field`.
    const field = { label, value };
    prompt.field = field;
    return (await choice) === 'ok' ? field.value : null;
  };
  // Enter in the text field chooses the action; an Enter that ends an IME composition does not.
  /** @param {KeyboardEvent} event */
  const submitField = (event) => {
    if (!prompt.field || event.isComposing || event.keyCode === 229) return;
    event.preventDefault();
    resolve('ok');
  };
  return { prompt, ask, askText, submitField, resolve };
};

/**
 * The outcome of an operation that has no place of its own
 * (Load Session, Save Session, a table import or export) shows in one notice
 * banner above the Result, never in a native `alert()`. Warnings show as a list.
 * @param {{ reactive: <T extends object>(value: T) => T }} options Vue's `reactive`.
 */
export const createNoticeBanner = ({ reactive }) => {
  // `text`, `items` and `link` word the shown notice each time they are read,
  // so the banner's render depends on the UI language.
  const notice = reactive({
    wordings: /** @type {NoticeWordings | null} */ (null),
    tone: /** @type {NoticeTone} */ ('info'),
    announcement: '',
    /** @returns {string} */
    get text() { return this.wordings ? this.wordings.text() : ''; },
    /** @returns {string[]} */
    get items() { return this.wordings ? this.wordings.items.map((item) => item()) : []; },
    /** @returns {{ text: string, label: string, href: string } | null} */
    get link() {
      const link = this.wordings ? this.wordings.link : null;
      return link ? { text: link.text(), label: link.label(), href: link.href } : null;
    }
  });
  let shownCount = 0;
  const dismiss = () => { Object.assign(notice, { wordings: null, tone: 'info', announcement: '' }); };
  /** @type {(request: NoticeRequest) => NoticeHandle} */
  const show = ({ text, items = [], tone }) => {
    shownCount += 1;
    const id = shownCount;
    const shown = items.filter((item) => item() !== '');
    const spoken = [text(), ...shown.map((item) => item())].join(' ');
    Object.assign(notice, {
      wordings: { text, items: shown, link: null },
      tone: tone || (shown.length ? 'warning' : 'info'),
      // The status region speaks a change of its text only: a notice that
      // repeats the shown one differs from it by a trailing no-break space.
      announcement: notice.announcement === spoken ? `${spoken}\u00A0` : spoken
    });
    return { dismiss: () => { if (shownCount === id) dismiss(); } };
  };
  // #1383 (PD-OI-108): the notice of a loaded Gallery example offers its
  // tutorial. The link joins the notice that shows and goes with it.
  /** @param {NoticeLink} link */
  const addLink = (link) => {
    const wordings = notice.wordings;
    if (!wordings) return;
    Object.assign(notice, { wordings: { ...wordings, link }, announcement: `${notice.announcement} ${link.text()}` });
  };
  return { notice, show, dismiss, addLink };
};

/**
 * The warning notice of a text kept until it is shown: English marked with
 * `msg` (utils/i18n.js), worded with `tm` where the notice shows.
 * @param {string} text
 * @returns {NoticeRequest}
 */
export const keptTextWarning = (text) => ({ text: () => tm(text), tone: 'warning' });

/**
 * The notice of a loaded Session: what the load changed without a loss is
 * listed under it.
 * @param {string} file The loaded file's name; empty when the load had none.
 * @param {Wording[]} items
 * @returns {NoticeRequest}
 */
export const sessionLoadedNotice = (file, items) => ({
  text: () => (file ? t('Session loaded: {file}.', { file }) : t('Session loaded.')), items
});

/** The notice of a loaded configuration file of an older gbdraw. @type {NoticeRequest} */
export const LEGACY_CONFIGURATION_NOTICE = Object.freeze({
  text: () => t('Configuration loaded from an older gbdraw. Save Session to keep it in the current format.')
});

/**
 * The notice of a saved Session, with its warnings.
 * @param {string} file The downloaded file's name.
 * @param {Wording[]} items
 * @returns {NoticeRequest}
 */
export const sessionSavedNotice = (file, items) => ({ text: () => t('Session saved: {file}.', { file }), items });

/**
 * #1383 (PD-OI-108): the line of a loaded Gallery example's notice that
 * offers its tutorial.
 * @param {string} href
 * @returns {NoticeLink}
 */
export const galleryTutorialLink = (href) => ({
  text: () => t('This diagram has a step-by-step tutorial.'),
  label: () => tx('English-only Gallery tutorial', 'Open the tutorial'),
  href
});
