// @ts-check

// A number field reads Arabic-Indic, Persian, Devanagari and Bengali digits as 0-9 and the
// Arabic decimal separator as '.', in every UI language (OV-2780). Chromium drops these
// characters from a number field and Firefox rejects the separator, so the text a key or a
// paste is about to insert goes in mapped instead. Text that an input method composes cannot
// be held back this way (its beforeinput cannot be cancelled) and stays the browser's.

/** The zero of each digit set read as 0-9: Arabic-Indic, Persian, Devanagari, Bengali. */
const OTHER_DIGIT_ZEROS = Object.freeze([0x0660, 0x06F0, 0x0966, 0x09E6]);
const ARABIC_DECIMAL_SEPARATOR = '٫';

/**
 * `text` with each Arabic-Indic, Persian, Devanagari and Bengali digit as 0-9 and the Arabic
 * decimal separator as '.'; every other character unchanged.
 * @param {string} text
 * @returns {string}
 */
export const asciiNumberText = (text) => Array.from(text, (character) => {
  if (character === ARABIC_DECIMAL_SEPARATOR) return '.';
  const code = character.codePointAt(0) ?? 0;
  const zero = OTHER_DIGIT_ZEROS.find((digitZero) => code >= digitZero && code <= digitZero + 9);
  return zero === undefined ? character : String(code - zero);
}).join('');

/**
 * A field whose text is a number: `type="number"`, or a text field whose `inputmode` is
 * numeric or decimal.
 * @param {EventTarget | null} target
 * @returns {target is HTMLInputElement}
 */
const isNumberField = (target) => target instanceof HTMLInputElement
  && (target.type === 'number' || target.inputMode === 'numeric' || target.inputMode === 'decimal');

/**
 * The page's `beforeinput` listener (capture phase): text typed or pasted into a number field
 * goes in with its digits as 0-9 and its decimal separator as '.', at the caret or over the
 * selection, through the browser's own insert so the field's input event and undo follow.
 * @param {InputEvent} event
 */
export const readOtherDigitsInNumberFields = (event) => {
  if (event.inputType !== 'insertText' && event.inputType !== 'insertFromPaste') return;
  if (!isNumberField(event.target)) return;
  const text = event.data ?? event.dataTransfer?.getData('text/plain') ?? '';
  const ascii = asciiNumberText(text);
  if (ascii === text) return;
  event.preventDefault();
  document.execCommand('insertText', false, ascii);
};
