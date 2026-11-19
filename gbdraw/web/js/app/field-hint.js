// @ts-check
// A field's hint (its placeholder): the one owner of the room it takes and of
// showing it whole.
// - `v-hint-width` on the input of a `.hint-width-field` writes the input's
//   hint, and its kind (number or text), into the hidden copy beside it
//   (`[data-hint-width]`), whose width the field takes (css/tailwind.input.css).
//   The template writes the hint once, in `:placeholder`.
// - A number field's spin buttons take room from its text. Their width is
//   measured once and published as `--number-spin-width` on the root element,
//   which the hidden copy of a number field's hint reads; the cut check below
//   uses the same value.
// - A hint wider than its field ends with an ellipsis (`input:placeholder-shown`
//   in css/tailwind.input.css). While such a hint is shown and the field is
//   under a hovering pointer or holds the focus, one house tooltip
//   (`.ui-tooltip`, placed by utils/tooltip-placement.js) shows it whole.
//   Escape hides it until the field is left. Screen readers read the hint from
//   the placeholder, so the tooltip is `aria-hidden`.
// app.js registers the directive and installs `FIELD_HINT_LISTENERS` on the window.

import { placeTooltip } from '../utils/tooltip-placement.js';

/** The custom property that holds the spin buttons' width. */
export const NUMBER_SPIN_WIDTH_PROPERTY = '--number-spin-width';

/** @type {number | null} */
let spinWidth = null;
/** @type {CanvasRenderingContext2D | null} */
let context = null;

/**
 * The width a number field's spin buttons take from its text: the scroll
 * width a long value gains in a number field over a text field of the same
 * width (0 where the browser reserves none). Measured once, with the classes
 * of the first field that needs it, and published on the root element.
 * @param {HTMLInputElement} field
 * @returns {number}
 */
const numberSpinWidth = (field) => {
  if (spinWidth !== null) return spinWidth;
  const widths = ['number', 'text'].map((type) => {
    const probe = document.createElement('input');
    probe.type = type;
    probe.className = field.className;
    Object.assign(probe.style, { position: 'absolute', visibility: 'hidden', width: '64px' });
    probe.value = '1'.repeat(40);
    document.body.append(probe);
    const width = probe.scrollWidth;
    probe.remove();
    return width;
  });
  spinWidth = Math.max(0, widths[0] - widths[1]);
  document.documentElement.style.setProperty(NUMBER_SPIN_WIDTH_PROPERTY, `${spinWidth}px`);
  return spinWidth;
};

/**
 * Write the field's hint and kind into its hidden copy, where they differ.
 * @param {HTMLInputElement} field
 */
const writeHintCopy = (field) => {
  const cell = field.parentElement;
  const copy = cell && cell.classList.contains('hint-width-field') ? cell.querySelector(':scope > [data-hint-width]') : null;
  if (!(copy instanceof HTMLElement)) throw new Error('v-hint-width needs its field in a .hint-width-field beside a [data-hint-width] copy.');
  const kind = field.type === 'number' ? 'number' : 'text';
  if (kind === 'number') numberSpinWidth(field);
  if (copy.dataset.hintWidth !== kind) copy.dataset.hintWidth = kind;
  if (copy.textContent !== field.placeholder) copy.textContent = field.placeholder;
};

// `v-hint-width` on the input of a `.hint-width-field`.
export const hintWidth = {
  /** @param {HTMLInputElement} field */
  mounted: writeHintCopy,
  /** @param {HTMLInputElement} field */
  updated: writeHintCopy
};

/**
 * Whether the field shows its hint and the hint is wider than the room for
 * its text.
 * @param {HTMLInputElement} field
 * @returns {boolean}
 */
const hintCut = (field) => {
  if (!field.placeholder || !field.isConnected || field.getClientRects().length === 0 || !field.matches(':placeholder-shown')) return false;
  context ??= document.createElement('canvas').getContext('2d');
  if (!context) return false;
  const style = getComputedStyle(field);
  context.font = `${style.fontStyle} ${style.fontWeight} ${style.fontSize} ${style.fontFamily}`;
  const room = field.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight)
    - (field.type === 'number' ? numberSpinWidth(field) : 0);
  return context.measureText(field.placeholder).width > room + 0.5;
};

// The field under a hovering (not touch) pointer, the field holding the focus,
// and the field whose tooltip Escape hid.
/** @type {HTMLInputElement | null} */
let pointed = null;
/** @type {HTMLInputElement | null} */
let focused = null;
/** @type {HTMLInputElement | null} */
let dismissed = null;
// The tooltip and the field whose hint it shows, while it is shown.
/** @type {{ element: HTMLDivElement, field: HTMLInputElement } | null} */
let shown = null;

const hide = () => {
  shown?.element.remove();
  shown = null;
};

/**
 * Show the whole hint of the pointed or focused field while its hint is cut,
 * placed at the field; hide it otherwise.
 */
const showCutHint = () => {
  if (dismissed && dismissed !== pointed && dismissed !== focused) dismissed = null;
  const field = [pointed, focused].find((candidate) => candidate && candidate !== dismissed && hintCut(candidate)) ?? null;
  if (!field) {
    hide();
    return;
  }
  if (!shown) {
    const element = document.createElement('div');
    element.className = 'ui-tooltip';
    element.setAttribute('role', 'tooltip');
    element.setAttribute('aria-hidden', 'true');
    document.body.append(element);
    shown = { element, field };
  }
  shown.field = field;
  if (shown.element.textContent !== field.placeholder) shown.element.textContent = field.placeholder;
  Object.assign(shown.element.style, placeTooltip(field.getBoundingClientRect(), shown.element.offsetHeight, window.innerWidth));
};

/**
 * The input an event reached, or null.
 * @param {Event} event
 * @returns {HTMLInputElement | null}
 */
const inputOf = (event) => (event.target instanceof HTMLInputElement ? event.target : null);

/** @param {Event} event */
const onPointerOver = (event) => {
  pointed = /** @type {PointerEvent} */ (event).pointerType === 'touch' ? null : inputOf(event);
  showCutHint();
};
/** @param {Event} event */
const onPointerOut = (event) => {
  if (inputOf(event) === pointed) pointed = null;
  showCutHint();
};
/** @param {Event} event */
const onFocusIn = (event) => {
  focused = inputOf(event);
  showCutHint();
};
/** @param {Event} event */
const onFocusOut = (event) => {
  if (inputOf(event) === focused) focused = null;
  showCutHint();
};
// After a key or a press, the page may have removed the field (a dialog that
// closed) without a pointerout or focusout: check again once it has rendered.
const recheckLater = () => { requestAnimationFrame(showCutHint); };
/** @param {Event} event */
const onKeyDown = (event) => {
  const key = /** @type {KeyboardEvent} */ (event);
  if (key.key === 'Escape' && !key.isComposing && shown) {
    // Escape closes the shown tooltip and goes no further, as a help tip's does.
    event.stopPropagation();
    dismissed = shown.field;
    hide();
    return;
  }
  recheckLater();
};
const onPlaceChange = () => { if (shown) showCutHint(); };

/**
 * The window listeners of the cut hint's tooltip: [type, listener, capture].
 * Escape is taken in the capture phase, before a dialog around the field.
 * Scrolls of any element are seen in the capture phase.
 * @type {ReadonlyArray<[string, (event: Event) => void, boolean]>}
 */
export const FIELD_HINT_LISTENERS = Object.freeze([
  ['pointerover', onPointerOver, false],
  ['pointerout', onPointerOut, false],
  ['focusin', onFocusIn, false],
  ['focusout', onFocusOut, false],
  ['input', showCutHint, false],
  ['keydown', onKeyDown, true],
  ['pointerdown', recheckLater, false],
  ['scroll', onPlaceChange, true],
  ['resize', onPlaceChange, false]
]);
