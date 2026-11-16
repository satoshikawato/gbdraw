// A template expression of gbdraw/web/index.html, read from the page and evaluated with the names it reads, so a node
// test words a stored value the way the shown page does (`tm`, `errorText` of a value the state keeps as data).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';

const page = readFileSync(new URL('../../../gbdraw/web/index.html', import.meta.url), 'utf8');

/**
 * The expression that `pattern`'s first group captures in index.html, as a function of the names it reads.
 * @param {RegExp} pattern Matches the binding once; its first group is the expression.
 * @returns {(scope: Record<string, unknown>) => unknown}
 */
export const templateBinding = (pattern) => {
  const matches = [...page.matchAll(new RegExp(pattern.source, `${pattern.flags.replace('g', '')}g`))];
  assert.equal(matches.length, 1, `index.html binds ${pattern} once`);
  const expression = matches[0][1].replaceAll('&quot;', '"').replaceAll('&amp;', '&');
  return (scope) => new Function(...Object.keys(scope), `return (${expression});`)(...Object.values(scope));
};
