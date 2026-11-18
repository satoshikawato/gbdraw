// @ts-check
// The whole hint of a field that cannot show it. A field's hint (its
// placeholder) that is wider than the field ends with an ellipsis
// (css/tailwind.input.css, `input:placeholder-shown`); this module is the one
// owner that keeps such a hint readable. When the pointer or the focus reaches
// an empty field whose hint does not fit, the field's title is the hint; when
// the hint fits again, or the field has a value, that title goes. A title the
// template gives the field is never replaced. app.js installs the listener on
// the window for `pointerover` and `focusin`.

/** The attribute that marks a title this module set. */
const HINT_TITLE = 'data-hint-title';

/** @type {CanvasRenderingContext2D | null} */
let context = null;
/** @type {number | null} */
let spinWidth = null;

/**
 * The width a number field's spin buttons take from its text: the scroll
 * width a long value gains in a number field over a text field of the same
 * width (0 where the browser reserves none). Measured once.
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
  return spinWidth;
};

/**
 * Whether an empty field's hint is wider than the room for its text.
 * @param {HTMLInputElement} field
 * @returns {boolean}
 */
const hintOverflows = (field) => {
  context ??= document.createElement('canvas').getContext('2d');
  if (!context) return false;
  const style = getComputedStyle(field);
  context.font = `${style.fontStyle} ${style.fontWeight} ${style.fontSize} ${style.fontFamily}`;
  const room = field.clientWidth - parseFloat(style.paddingLeft) - parseFloat(style.paddingRight)
    - (field.type === 'number' ? numberSpinWidth(field) : 0);
  return context.measureText(field.placeholder).width > room + 0.5;
};

/**
 * Give the field under the pointer or the focus its hint as its title while
 * the hint is cut, and take that title away when it is not.
 * @param {Event} event
 */
export const titleCutHint = (event) => {
  const field = event.target;
  if (!(field instanceof HTMLInputElement) || !field.placeholder) return;
  const ownTitle = field.hasAttribute('title') && !field.hasAttribute(HINT_TITLE);
  if (ownTitle) return;
  if (field.value === '' && hintOverflows(field)) {
    field.setAttribute('title', field.placeholder);
    field.setAttribute(HINT_TITLE, '');
  } else if (field.hasAttribute(HINT_TITLE)) {
    field.removeAttribute('title');
    field.removeAttribute(HINT_TITLE);
  }
};
