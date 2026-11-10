// @ts-check
// `v-keeps-list-focus` on the element that holds a keyed list: a row the list
// moves keeps the focus it held (OV-3070). A keyed list moves a row's element
// when the row changes place (Vue takes it out and puts it back), and Chromium
// then blurs the focused control inside it and ends its edit: `change` (a
// typed but uncommitted name would be committed), `blur` and `focusout`. While
// the list is patched, this directive holds those events back on the element
// that had the focus, and after the patch it focuses that element again when
// the move left the focus nowhere (on the page body). The focus then stays on
// the control it was on, which a keyed list keeps for its row. The rows'
// directives run their `updated` before the list's, so a typed field asks
// `listHoldsFocus` whether it still counts as focused.
//
// A `change` held back is still owed: Chromium takes the text it carried as
// committed and fires no `change` for it again. The field's next commit end
// (leaving it, or an Enter no handler took) fires it then, unless the text
// changed again since, in which case the browser fires its own
// (`oweChange`).

/** @type {WeakMap<Element, { focused: HTMLElement, stop: () => void }>} */
const heldFocus = new WeakMap();
/** @type {WeakSet<Element>} The elements whose focus a list holds while it is patched. */
const held = new WeakSet();
const ENDS = /** @type {const} */ (['change', 'blur', 'focusout', 'focus', 'focusin']);

/**
 * The `change` held back from `field` while it carried `text`, fired at its
 * next commit end unless the browser fires one itself (the text changed).
 * @param {HTMLElement} field
 * @param {string} text
 */
const oweChange = (field, text) => {
  const settle = () => {
    field.removeEventListener('change', settle, true);
    field.removeEventListener('blur', end, true);
    field.removeEventListener('keypress', end);
  };
  /** @param {Event} event */
  const end = (event) => {
    if (event.type === 'keypress' && (/** @type {KeyboardEvent} */ (event).key !== 'Enter' || event.defaultPrevented)) return;
    settle();
    if (/** @type {HTMLInputElement} */ (field).value === text) field.dispatchEvent(new Event('change', { bubbles: true }));
  };
  field.addEventListener('change', settle, true);
  field.addEventListener('blur', end, true);
  field.addEventListener('keypress', end);
};

/** @param {Element} list */
const release = (list) => {
  const holding = heldFocus.get(list);
  if (!holding) return;
  holding.stop();
  held.delete(holding.focused);
  heldFocus.delete(list);
};

/**
 * Whether a list holds the focus of `element` while it is patched: the element
 * had the focus before the patch and gets it back after it. A field's own
 * directive counts it as focused then (app/committed-value.js), so the move
 * does not replace its typed text.
 * @param {Element} element
 */
export const listHoldsFocus = (element) => held.has(element);

export const keepsListFocus = {
  /** @param {Element} list */
  beforeUpdate(list) {
    release(list);
    const focused = list.ownerDocument.activeElement;
    if (!(focused instanceof HTMLElement) || !list.contains(focused)) return;
    /** @param {Event} event */
    const hold = (event) => {
      if (event.target !== focused) return;
      event.stopImmediatePropagation();
      if (event.type === 'change') oweChange(focused, /** @type {HTMLInputElement} */ (focused).value);
    };
    // Capture listeners on the list run before the focused element's own.
    ENDS.forEach((type) => list.addEventListener(type, hold, true));
    heldFocus.set(list, { focused, stop: () => ENDS.forEach((type) => list.removeEventListener(type, hold, true)) });
    held.add(focused);
  },
  /** @param {Element} list */
  updated(list) {
    const holding = heldFocus.get(list);
    if (!holding) return;
    const { focused } = holding;
    const active = list.ownerDocument.activeElement;
    if (focused.isConnected && list.contains(focused) && (!active || active === list.ownerDocument.body)) {
      focused.focus({ preventScroll: true });
    }
    release(list);
  },
  /** @param {Element} list */
  beforeUnmount: release
};
