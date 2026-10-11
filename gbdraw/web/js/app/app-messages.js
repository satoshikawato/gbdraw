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
 * @typedef {'info' | 'warning'} NoticeTone
 *
 * @typedef {object} NoticeRequest
 * @property {string} text
 * @property {string[]} [items] Warnings listed under the text.
 * @property {NoticeTone} [tone]
 *
 * Shows the outcome of an operation in the notice banner, replacing the
 * notice it shows.
 * @typedef {(request: NoticeRequest) => void} ShowNotice
 */

/**
 * UX-A07: a confirmation asks with the app's choice dialog, never a native
 * `confirm()`. One question shows at a time: a new question cancels the open one.
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
    cancelLabel: 'Cancel'
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
      show: true
    });
    return new Promise((settle) => { answer = settle; });
  };
  return { prompt, ask, resolve };
};

/**
 * UX-A07, UX-C12: the outcome of an operation that has no place of its own
 * (Load Session, Save Session, a table import or export) shows in one notice
 * banner above the Result, never in a native `alert()`. Warnings show as a list.
 * @param {{ reactive: <T extends object>(value: T) => T }} options Vue's `reactive`.
 */
export const createNoticeBanner = ({ reactive }) => {
  const notice = reactive({ text: '', items: /** @type {string[]} */ ([]), tone: /** @type {NoticeTone} */ ('info') });
  /** @type {ShowNotice} */
  const show = ({ text, items = [], tone }) => {
    Object.assign(notice, {
      text,
      items: items.filter(Boolean),
      tone: tone || (items.some(Boolean) ? 'warning' : 'info')
    });
  };
  const dismiss = () => { Object.assign(notice, { text: '', items: [], tone: 'info' }); };
  return { notice, show, dismiss };
};
