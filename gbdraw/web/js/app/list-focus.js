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
//
// When the focus does not come back to the element after the patch, its edit
// ended during the patch: another element took the focus (a dialog the same
// render mounts), or its row left the list. The events held back then end it
// after the patch (`endHeld`).

/**
 * A `change` held back from a field: `pay` fires it now unless the text changed since, `drop` forgets it.
 * @typedef {{ pay: () => void, drop: () => void }} OwedChange
 */
/**
 * The focus a list holds while it is patched.
 * @typedef {object} HeldFocus
 * @property {HTMLElement} focused The element that had the focus when the patch began.
 * @property {Set<string>} ended The types of the events held back from it in this patch.
 * @property {() => void} stop Removes the hold's listeners.
 */

/** @type {WeakMap<Element, HeldFocus>} */
const heldFocus = new WeakMap();
/** @type {WeakSet<Element>} The elements whose focus a list holds while it is patched. */
const held = new WeakSet();
const ENDS = /** @type {const} */ (['change', 'blur', 'focusout', 'focus', 'focusin']);
/** @type {WeakMap<HTMLElement, OwedChange>} The `change` a field is owed: the last one held back from it. */
const owedChanges = new WeakMap();

/**
 * The `change` held back from `field` while it carried `text`, fired at its
 * next commit end unless the browser fires one itself (the text changed).
 * @param {HTMLElement} field
 * @param {string} text
 */
const oweChange = (field, text) => {
  owedChanges.get(field)?.drop();
  const drop = () => {
    if (owedChanges.get(field) === owed) owedChanges.delete(field);
    field.removeEventListener('change', drop, true);
    field.removeEventListener('blur', end, true);
    field.removeEventListener('keypress', end);
  };
  const pay = () => {
    drop();
    if (/** @type {HTMLInputElement} */ (field).value === text) field.dispatchEvent(new Event('change', { bubbles: true }));
  };
  /** @param {Event} event */
  const end = (event) => {
    if (event.type === 'keypress' && (/** @type {KeyboardEvent} */ (event).key !== 'Enter' || event.defaultPrevented)) return;
    pay();
  };
  /** @type {OwedChange} */
  const owed = { pay, drop };
  owedChanges.set(field, owed);
  field.addEventListener('change', drop, true);
  field.addEventListener('blur', end, true);
  field.addEventListener('keypress', end);
};

/**
 * Ends the edit of an element whose focus did not come back after the patch:
 * the held `change`, then the held `blur` and `focusout`, which History's
 * input adapter ends the element's step on. An element whose row left the
 * list commits nothing: its handler would name the row that took its place.
 * @param {HeldFocus} holding
 * @param {Element} list
 */
const endHeld = ({ focused, ended }, list) => {
  const relatedTarget = list.ownerDocument.activeElement;
  const owed = owedChanges.get(focused);
  if (owed && list.contains(focused)) owed.pay();
  else if (owed) owed.drop();
  if (ended.has('blur')) focused.dispatchEvent(new FocusEvent('blur', { relatedTarget }));
  if (ended.has('focusout')) focused.dispatchEvent(new FocusEvent('focusout', { bubbles: true, relatedTarget }));
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
    /** @type {HeldFocus} */
    const holding = { focused, ended: new Set(), stop: () => ENDS.forEach((type) => list.removeEventListener(type, hold, true)) };
    /** @param {Event} event */
    const hold = (event) => {
      if (event.target !== focused) return;
      event.stopImmediatePropagation();
      holding.ended.add(event.type);
      if (event.type === 'change') oweChange(focused, /** @type {HTMLInputElement} */ (focused).value);
    };
    // Capture listeners on the list run before the focused element's own.
    ENDS.forEach((type) => list.addEventListener(type, hold, true));
    heldFocus.set(list, holding);
    held.add(focused);
  },
  /** @param {Element} list */
  updated(list) {
    const holding = heldFocus.get(list);
    if (!holding) return;
    const { focused } = holding;
    const active = list.ownerDocument.activeElement;
    if (list.contains(focused) && (!active || active === list.ownerDocument.body)) {
      focused.focus({ preventScroll: true });
    }
    release(list);
    if (list.ownerDocument.activeElement !== focused) endHeld(holding, list);
  },
  /** @param {Element} list */
  beforeUnmount: release
};
