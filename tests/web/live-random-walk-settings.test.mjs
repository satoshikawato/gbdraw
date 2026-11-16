// The live setting rows the random walk's `setting` step drives
// (tests/web/helpers/live-setting-rows.mjs): every row of `SETTING_APPLICATION`
// that the automatic rerender applies from the draft has, in index.html, a
// control the step can set (a control of a type it sets, whose nearest
// `data-setting` scope is the row, outside the subtrees that record no History
// step), or is named with the walk kind that writes its draft. A new live row
// without either fails here with its id instead of never being drawn.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { SETTING_APPLICATION } from '../../gbdraw/web/js/services/setting-application.js';
import { isIgnoredTarget } from '../../gbdraw/web/js/app/history-inputs.js';
import { ROWS_DRIVEN_BY_KIND, liveRequestRows, settingControlAction } from './helpers/live-setting-rows.mjs';
import { WALK_KINDS } from './helpers/live-random-walk-plan.mjs';
import { descendants, parseElementTree } from './helpers/html-element-tree.mjs';

const INDEX_HTML = new URL('../../gbdraw/web/index.html', import.meta.url);

/** @typedef {import('./helpers/html-element-tree.mjs').HtmlElement} ScannedElement */

// An attribute selector list such as `[a], [b="c"]`, the form isIgnoredTarget passes to `closest`.
/** @param {string} selector @returns {(element: ScannedElement) => boolean} */
const attributeSelector = (selector) => {
  const parts = selector.split(',').map((part) => {
    const match = /^\[([\w-]+)(?:="([^"]*)")?\]$/.exec(part.trim());
    if (!match) throw new Error(`the scanned elements match attribute selectors only; got ${selector}`);
    return /** @param {ScannedElement} element */ (element) => element.attributes.has(match[1])
      && (match[2] === undefined || element.attributes.get(match[1]) === match[2]);
  });
  return (element) => parts.some((matches) => matches(element));
};
/** The element as `closest` sees it. @param {ScannedElement} element */
const asClosestTarget = (element) => ({
  closest: (/** @type {string} */ selector) => {
    const matches = attributeSelector(selector);
    for (let node = /** @type {ScannedElement | null} */ (element); node; node = node.parent) if (matches(node)) return node;
    return null;
  }
});

// The elements of an HTML text, by the scanner tests/web/setting-application.test.mjs reads index.html with.
/** @param {string} html @returns {ScannedElement[]} */
const scanElements = (html) => [...descendants(parseElementTree(html))];

/** The static `data-setting` of the element's nearest scope, or ''. @param {ScannedElement} element */
const nearestSetting = (element) => {
  for (let node = /** @type {ScannedElement | null} */ (element); node; node = node.parent) {
    if (node.attributes.has('data-setting')) return node.attributes.get('data-setting') ?? '';
    if (node.attributes.has(':data-setting')) return '';
  }
  return '';
};

/** The controls of each row that a `setting` step can set, by row id. */
const settableControls = () => {
  /** @type {Map<string, ScannedElement[]>} */
  const byRow = new Map();
  for (const element of scanElements(readFileSync(INDEX_HTML, 'utf8'))) {
    if (!settingControlAction(element.tag, element.attributes.get('type'))) continue;
    if (isIgnoredTarget(/** @type {Element} */ (/** @type {unknown} */ (asClosestTarget(element))))) continue;
    const row = nearestSetting(element);
    if (row) byRow.set(row, [...(byRow.get(row) ?? []), element]);
  }
  return byRow;
};

test('every live setting row has a control the setting step sets, or the walk kind that writes its draft', () => {
  const controls = settableControls();
  const rows = liveRequestRows(SETTING_APPLICATION);
  assert.ok(rows.length > 0, 'SETTING_APPLICATION lists live rows with request paths');
  const undriven = rows.filter(({ id }) => !(id in ROWS_DRIVEN_BY_KIND) && !controls.get(id)?.length).map(({ id }) => id);
  assert.deepEqual(undriven, [], 'live setting rows no setting step can drive: give each a control in its data-setting scope '
    + 'that records a History step, or name the walk kind that writes its draft in ROWS_DRIVEN_BY_KIND');
  for (const [id, kind] of Object.entries(ROWS_DRIVEN_BY_KIND)) {
    assert.ok(rows.some((row) => row.id === id), `${id} is a live setting row with request paths`);
    assert.ok(WALK_KINDS.includes(kind), `${id}: ${kind} is a walk kind`);
  }
});

test('the scan finds the controls by their nearest scope and leaves out those that record no History step', () => {
  const [plain, nested, ignored, typed] = [
    '<div data-setting="legend"><input type="number" aria-label="a"></div>',
    '<div data-setting="legend"><div data-setting="layout"><select aria-label="b"><option value="x"></option></select></div></div>',
    '<div data-setting="legend"><fieldset data-history-ignore><input type="text" :disabled="x > 1" aria-label="c"></fieldset></div>',
    '<div data-setting="legend"><input type="file"><input type="range"><textarea></textarea></div>'
  ].map((html) => scanElements(html).filter((element) => settingControlAction(element.tag, element.attributes.get('type'))));
  assert.deepEqual(plain.map(nearestSetting), ['legend']);
  assert.deepEqual(nested.map(nearestSetting), ['layout']);
  assert.deepEqual(ignored.map((element) => isIgnoredTarget(/** @type {Element} */ (/** @type {unknown} */ (asClosestTarget(element))))), [true]);
  assert.deepEqual(typed.map(({ tag }) => tag), ['textarea']);
  assert.equal(settingControlAction('INPUT', null), 'fill');
  assert.equal(settingControlAction('input', 'checkbox'), 'toggle');
  assert.equal(settingControlAction('select', undefined), 'choose');
});
