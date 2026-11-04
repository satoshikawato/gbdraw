// @ts-check
// A control's name and the value it holds now, as the page shows them: the
// one owner that the History list (a step's detail, PD-OI-112) and the
// Feature placement layout dialog (Q3) read.
import { t } from './i18n.js';

/** The input types whose value is text the user typed. */
export const TEXT_INPUT_TYPES = new Set([
  '',
  'date',
  'datetime-local',
  'email',
  'month',
  'number',
  'password',
  'search',
  'tel',
  'text',
  'time',
  'url',
  'week'
]);

/**
 * A control's name, as shown when it was read, and the value it held: a
 * text, a checkbox's state, or a color. A description without `name` names
 * no control.
 * @typedef {{ name?: string, value?: string, checked?: boolean, color?: string }} ControlDescription
 */
/**
 * An input, textarea, select, button, or contenteditable element; `type` exists on some of them only.
 * @typedef {Element & { type?: string }} DescribedControl
 */

/** @param {Element | null | undefined} element */
const shownText = (element) => String(/** @type {HTMLElement | null | undefined} */ (element)?.innerText ?? element?.textContent ?? '')
  .replace(/\s+/g, ' ').trim();

/**
 * The control's accessible name as shown now: its `aria-label`, the elements
 * its `aria-labelledby` names, its labels, or a button's text.
 * @param {DescribedControl} element
 */
const accessibleName = (element) => {
  const aria = element.getAttribute?.('aria-label')?.trim();
  if (aria) return aria;
  const labelledBy = String(element.getAttribute?.('aria-labelledby') || '').split(/\s+/).filter(Boolean)
    .map((id) => shownText(element.ownerDocument?.getElementById(id))).join(' ').trim();
  if (labelledBy) return labelledBy;
  const labels = [.../** @type {HTMLInputElement} */ (element).labels ?? []].map(shownText).join(' ').trim();
  if (labels) return labels;
  return String(element.tagName || '').toLowerCase() === 'button' ? shownText(element) : '';
};

/**
 * A radio's group name: its fieldset's legend, or its radio group's name.
 * @param {DescribedControl} element
 */
const radioGroupName = (element) => {
  const group = element.closest?.('fieldset, [role="radiogroup"]');
  if (!group) return '';
  const legend = String(group.tagName || '').toLowerCase() === 'fieldset' ? group.querySelector?.('legend') : null;
  return legend ? shownText(legend) : accessibleName(group);
};

/**
 * The control's name and the value it holds now: a select's option text, a
 * checkbox's state, a color, a text field's value, a radio's label under its
 * group's name (empty when the group has none); a button, a file input, or an
 * editable element gives its name alone. A control without a name has none.
 * @param {DescribedControl} element
 * @returns {ControlDescription | null}
 */
export const describeControl = (element) => {
  const name = accessibleName(element);
  if (!name) return null;
  const tag = String(element.tagName || '').toLowerCase();
  const type = String(element.type || '').toLowerCase();
  if (type === 'radio') return { name: radioGroupName(element), value: name };
  if (tag === 'select') {
    const option = /** @type {HTMLSelectElement} */ (element).selectedOptions?.[0];
    return { name, value: option ? shownText(option) : '' };
  }
  const value = String(/** @type {HTMLInputElement} */ (element).value ?? '');
  if (type === 'checkbox') return { name, checked: Boolean(/** @type {HTMLInputElement} */ (element).checked) };
  if (type === 'color') return { name, value, color: value };
  if (tag === 'textarea' || (tag === 'input' && type !== 'password' && TEXT_INPUT_TYPES.has(type))) return { name, value };
  return { name };
};

/**
 * The value a description holds, as shown: a checkbox's state as On or Off
 * in the shown language, else its text.
 * @param {ControlDescription} description
 */
export const controlValueText = ({ value, checked }) => (checked === undefined ? value ?? '' : checked ? t('On') : t('Off'));
