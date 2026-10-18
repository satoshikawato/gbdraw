// @ts-check

/**
 * The style of a choice button: one primary choice, secondary choices, or a
 * choice that removes work (danger).
 * @typedef {'primary' | 'secondary' | 'danger'} ChoiceTone
 *
 * @typedef {object} ChoiceAction
 * @property {string} id What the promise of `ask` resolves to.
 * @property {string} label The button text; it names the outcome.
 * @property {ChoiceTone} [tone]
 *
 * @typedef {object} ChoiceRequest
 * @property {string} title
 * @property {string} [description] One or two sentences: what the choice is about.
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
 * @property {string} text
 * @property {string[]} [items] Warnings listed under the text.
 * @property {NoticeTone} [tone]
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
    items: /** @type {string[]} */ ([]),
    note: '',
    actions: /** @type {ChoiceAction[]} */ ([]),
    cancelLabel: 'Cancel',
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
      items: [...(request.items || [])],
      note: request.note || '',
      actions: request.actions.map((action) => ({ tone: 'secondary', ...action })),
      cancelLabel: request.cancelLabel || 'Cancel',
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
  const notice = reactive({
    text: '', items: /** @type {string[]} */ ([]), tone: /** @type {NoticeTone} */ ('info'), announcement: ''
  });
  let shownCount = 0;
  const dismiss = () => { Object.assign(notice, { text: '', items: [], tone: 'info', announcement: '' }); };
  /** @type {(request: NoticeRequest) => NoticeHandle} */
  const show = ({ text, items = [], tone }) => {
    shownCount += 1;
    const id = shownCount;
    const shown = items.filter(Boolean);
    const spoken = [text, ...shown].join(' ');
    Object.assign(notice, {
      text,
      items: shown,
      tone: tone || (shown.length ? 'warning' : 'info'),
      // The status region speaks a change of its text only: a notice that
      // repeats the shown one differs from it by a trailing no-break space.
      announcement: notice.announcement === spoken ? `${spoken}\u00A0` : spoken
    });
    return { dismiss: () => { if (shownCount === id) dismiss(); } };
  };
  return { notice, show, dismiss };
};
