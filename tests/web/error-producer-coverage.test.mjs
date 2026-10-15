// R6 / G-G(2): Web JavaScript raises user-correctable failures through
// diagnosticError(code, context), and the English message classifier does not
// grow. Both baselines are registered design-rule allowlists
// (tools/web-design-rule-guards.json, R6): an observed count may not exceed
// its baseline, lowering a baseline is optional, and raising one is an
// authority-only change (the Python half is
// tests/test_web_error_producer_coverage.py).
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';
import { sep } from 'node:path';
import test from 'node:test';
import { normalizeUserFacingError } from '../../gbdraw/web/js/utils/error-normalization.js';
import { compareCountBaseline, reportLowerable } from './helpers/count-baseline.mjs';

const WEB_ROOT = new URL('../../gbdraw/web/js/', import.meta.url);
const UNCLASSIFIED = new Set(['UNKNOWN', 'VALIDATION_UNCLASSIFIED']);
const FILLS = ['', '0', '1', 'PRIVATE', '#1'];

// Every `.js` file under gbdraw/web/js/, workers included. The throw sites are
// counted in all of them, as the Python half counts the whole package: a
// hand-kept list of owners would free its sites for the other modules whenever
// a file left it.
const WEB_MODULES = readdirSync(WEB_ROOT, { recursive: true })
  .map((file) => file.split(sep).join('/'))
  .filter((file) => file.endsWith('.js'))
  .sort();

// The throw sites of all modules whose message normalizes to UNKNOWN or
// VALIDATION_UNCLASSIFIED (mostly internal invariants today). One total (kind
// `count`), so a site that moves to another module with its code does not
// change it.
const UNCLASSIFIED_THROW_BASELINE = 584;

// Literal parts of the first argument of `throw new Error(`; null is a hole.
const literalAt = (source, index) => {
  const quote = source[index];
  if (!['\'', '"', '`'].includes(quote)) return null;
  const parts = [];
  let current = '';
  for (let cursor = index + 1; cursor < source.length;) {
    const char = source[cursor];
    if (char === '\\') { current += source[cursor + 1]; cursor += 2; continue; }
    if (char === quote) { parts.push(current); return parts; }
    if (quote === '`' && char === '$' && source[cursor + 1] === '{') {
      parts.push(current, null);
      current = '';
      let depth = 1;
      cursor += 2;
      while (cursor < source.length && depth) {
        if (source[cursor] === '{') depth += 1;
        else if (source[cursor] === '}') depth -= 1;
        cursor += 1;
      }
      continue;
    }
    current += char;
    cursor += 1;
  }
  return null;
};
const renderings = (parts) => {
  const holes = parts.filter((part) => part === null).length;
  const choices = holes <= 2
    ? (function* combos(count) {
        if (!count) { yield []; return; }
        for (const fill of FILLS) for (const rest of combos(count - 1)) yield [fill, ...rest];
      })(holes)
    : FILLS.map((fill) => Array(holes).fill(fill));
  return [...choices].map((choice) => {
    let next = 0;
    return parts.map((part) => (part === null ? choice[next++] : part)).join('');
  });
};
const unclassifiedThrowSites = (file) => {
  const source = readFileSync(new URL(file, WEB_ROOT), 'utf8');
  const sites = [];
  for (const match of source.matchAll(/throw new (?:Error|TypeError|RangeError)\(\s*/g)) {
    const parts = literalAt(source, match.index + match[0].length);
    const classified = parts && renderings(parts).some((message) => (
      !UNCLASSIFIED.has(normalizeUserFacingError(new Error(message)).code)
    ));
    if (!classified) sites.push(source.slice(match.index, match.index + 100).split('\n')[0]);
  }
  return sites;
};

test('JS throw sites normalize to a recognized diagnostic or do not grow in number (G-G(2))', (t) => {
  const sitesByModule = WEB_MODULES.map((file) => [file, unclassifiedThrowSites(file).length]);
  const total = sitesByModule.reduce((sum, [, sites]) => sum + sites, 0);
  if (total < UNCLASSIFIED_THROW_BASELINE) {
    t.diagnostic(`UNCLASSIFIED_THROW_BASELINE may be lowered: ${UNCLASSIFIED_THROW_BASELINE} -> ${total}`);
  }
  assert.ok(total <= UNCLASSIFIED_THROW_BASELINE,
    `${total} unclassified throw sites under gbdraw/web/js/, above UNCLASSIFIED_THROW_BASELINE ${UNCLASSIFIED_THROW_BASELINE}; `
      + 'raise the new ones with diagnosticError(code, context). Sites per module (compare with dev):\n'
      + sitesByModule.filter(([, sites]) => sites).map(([file, sites]) => `${file} ${sites}`).join('\n'));
});

// The message-classification tables in the JS wording owner (R6 ratchet).
const NATIVE_VALIDATION_BASELINE = { exactMessages: 92, patterns: 19 };

test('the JS message-classification tables only shrink (R6 ratchet)', async (t) => {
  const source = readFileSync(new URL('utils/error-normalization.js', WEB_ROOT), 'utf8');
  // The module has no imports, so a copy can expose its private table size.
  const probe = await import(`data:text/javascript,${encodeURIComponent(
    `${source}\nexport const nativeValidationSize = NATIVE_VALIDATIONS.size;`
  )}`);
  const body = source.slice(source.indexOf('const nativeValidation = '), source.indexOf('// The native track validator'));
  const actual = {
    exactMessages: probe.nativeValidationSize,
    patterns: (body.match(/\/\^/g) || []).length
  };
  const { above, lowerable } = compareCountBaseline(actual, NATIVE_VALIDATION_BASELINE);
  reportLowerable(t, 'NATIVE_VALIDATION_BASELINE', lowerable);
  assert.deepEqual(above.map(({ key, observed, recorded }) => `${key}: ${observed}, above ${recorded}`), [],
    'Raise new failures with diagnosticError; raising NATIVE_VALIDATION_BASELINE is an authority-only change.');
});
