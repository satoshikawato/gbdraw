// @ts-check
// The Result executor records Python's value of each paint attribute the
// first time it changes it on an element, in `data-gbdraw-base-<attribute>`
// (empty: Python drew none), so a later reconcile can return the element to
// it. The records travel with the Result through History and Sessions;
// exports strip them.
const RESULT_BASE_ATTRIBUTES = Object.freeze(['fill', 'stroke', 'stroke-width', 'display']);
/** @param {string} attribute */
export const resultBaseAttribute = (attribute) => `data-gbdraw-base-${attribute}`;
export const RESULT_BASE_SELECTOR = RESULT_BASE_ATTRIBUTES
  .map((attribute) => `[${resultBaseAttribute(attribute)}]`).join(', ');

// The value Python drew for `attribute` on an element: its record, else the
// attribute as it is (null: none).
/**
 * @param {Element | null | undefined} element
 * @param {string} attribute
 * @returns {string | null}
 */
export const pythonDrawnAttribute = (element, attribute) => {
  const base = element?.getAttribute?.(resultBaseAttribute(attribute));
  if (base !== null && base !== undefined) return base === '' ? null : base;
  return element?.getAttribute?.(attribute) ?? null;
};

/** @param {Element | null | undefined} root */
export const stripResultBaseAttributes = (root) => {
  if (!root) return;
  [root, ...Array.from(root.querySelectorAll(RESULT_BASE_SELECTOR))].forEach((element) => {
    RESULT_BASE_ATTRIBUTES.forEach((attribute) => element.removeAttribute(resultBaseAttribute(attribute)));
  });
};
