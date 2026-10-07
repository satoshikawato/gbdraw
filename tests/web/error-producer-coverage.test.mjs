// R6 / G-G(2): JS validation owners raise user-correctable failures through
// diagnosticError(code, context); the English message classifier may only
// shrink. Both counts are shrink-only baselines: a file may never gain an
// unclassified throw site, and a baseline is lowered as soon as a site is
// migrated (the Python half is tests/test_web_error_producer_coverage.py).
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import test from 'node:test';
import { normalizeUserFacingError } from '../../gbdraw/web/js/services/error-normalization.js';

const WEB_ROOT = new URL('../../gbdraw/web/js/', import.meta.url);
const UNCLASSIFIED = new Set(['UNKNOWN', 'VALIDATION_UNCLASSIFIED']);
const FILLS = ['', '0', '1', 'PRIVATE', '#1'];

// Shrink-only: validation owner -> throw sites whose message normalizes to
// UNKNOWN or VALIDATION_UNCLASSIFIED (mostly internal invariants today). The
// Result admission and preview binding owners are counted too, because a
// failure there ends Generate (OV-13).
const UNCLASSIFIED_THROW_BASELINE = {
  'app/annotations/table-codec.js': 13,
  'app/annotations/validation.js': 0,
  'app/annotations/record-catalog.js': 0,
  'app/candidate-render.js': 3,
  'app/circular-track-slots.js': 2,
  'app/current-option-values.js': 7,
  'app/feature-editor/label-actions.js': 0,
  'app/feature-metadata-extraction.js': 2,
  'app/file-imports.js': 0,
  'app/legend-layout/decoration-continuity.js': 0,
  'app/linear-track-slots.js': 34,
  'app/preview-runtime.js': 7,
  'app/record-display-options.js': 23,
  'app/run-analysis.js': 18,
  'app/track-slot-validation.js': 31,
  'mode-profiles.js': 9,
  'services/config.js': 15,
  'services/current-option-values.js': 7,
  'services/feature-metadata-extraction.js': 2,
  'services/file-imports.js': 0,
  'services/session-file.js': 6,
  'services/session-import-client.js': 0,
  'services/session-request.js': 94,
  'services/svg-result-ingestion.js': 22,
  'services/track-slot-validation.js': 31,
  'utils/feature-rendering.js': 3,
  'utils/optional-positive-number.js': 0
};

// Web layering slices S3, S4 and S6 move four of these modules from app/ to
// services/. A move PR is a runtime change and may not rename a key of this
// registered baseline, so each new key sits beside the old one (same count when
// written) and names the file it replaces. Exactly one of the pair exists; the
// move PR deletes the old key, and the cleanup after the slices deletes this
// table.
const PENDING_MOVES = {
  'services/current-option-values.js': 'app/current-option-values.js',
  'services/feature-metadata-extraction.js': 'app/feature-metadata-extraction.js',
  'services/file-imports.js': 'app/file-imports.js',
  'services/track-slot-validation.js': 'app/track-slot-validation.js'
};
const fileExists = (file) => existsSync(new URL(file, WEB_ROOT));

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

test('JS validation throw sites normalize to a recognized diagnostic or shrink (G-G(2))', () => {
  for (const [file, baseline] of Object.entries(UNCLASSIFIED_THROW_BASELINE)) {
    if (!fileExists(file)) {
      const old = PENDING_MOVES[file];
      assert.ok(old && fileExists(old) && old in UNCLASSIFIED_THROW_BASELINE,
        `${file}: no such file; delete the key (a move PR deletes the old key, not this one)`);
      continue;
    }
    const sites = unclassifiedThrowSites(file);
    assert.ok(sites.length <= baseline,
      `${file}: new unclassified throw sites; raise them with diagnosticError(code, context):\n${sites.join('\n')}`);
    assert.equal(sites.length, baseline, `${file}: lower UNCLASSIFIED_THROW_BASELINE to ${sites.length}`);
  }
});

// The message-classification tables in the JS wording owner (R6 ratchet).
const NATIVE_VALIDATION_BASELINE = { exactMessages: 92, patterns: 20 };

test('the JS message-classification tables only shrink (R6 ratchet)', async () => {
  const source = readFileSync(new URL('services/error-normalization.js', WEB_ROOT), 'utf8');
  // The module has no imports, so a copy can expose its private table size.
  const probe = await import(`data:text/javascript,${encodeURIComponent(
    `${source}\nexport const nativeValidationSize = NATIVE_VALIDATIONS.size;`
  )}`);
  const body = source.slice(source.indexOf('const nativeValidation = '), source.indexOf('// The native track validator'));
  const actual = {
    exactMessages: probe.nativeValidationSize,
    patterns: (body.match(/\/\^/g) || []).length
  };
  assert.deepEqual(actual, NATIVE_VALIDATION_BASELINE,
    'Raise new failures with diagnosticError and lower NATIVE_VALIDATION_BASELINE when a row or pattern is removed.');
});
