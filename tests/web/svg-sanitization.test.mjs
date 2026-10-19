import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import {
  SVG_SANITIZE_OPTIONS,
  sanitizeSvgContent
} from '../../gbdraw/web/js/services/svg-sanitization.js';

const trustedRendererSemanticAttributes = [
  'data-edge-kind',
  'data-gbdraw-annotation-clipped',
  'data-gbdraw-definition-part',
  'data-gbdraw-composition',
  'data-gbdraw-composition-role',
  'data-gbdraw-composition-schema',
  'data-gbdraw-feature-offset-y',
  'data-gbdraw-label-binding-schema',
  'data-gbdraw-label-visibility-preview',
  'data-gbdraw-slot-id',
  'data-gbdraw-slot-renderer',
  'data-ortholog-path-id',
  'data-query-orthogroup-assignment-reason',
  'data-query-orthogroup-confidence',
  'data-query-orthogroup-member-count',
  'data-query-orthogroup-representative',
  'data-query-orthogroup-role',
  'data-query-row',
  'data-query-stable-feature-svg-id',
  'data-rbh-orthogroup-id',
  'data-record-column',
  'data-record-index',
  'data-record-key',
  'data-record-row',
  'data-render-role',
  'data-subject-orthogroup-assignment-reason',
  'data-subject-orthogroup-confidence',
  'data-subject-orthogroup-member-count',
  'data-subject-orthogroup-representative',
  'data-subject-orthogroup-role',
  'data-subject-row',
  'data-subject-stable-feature-svg-id'
];

test('trusted renderer semantic hooks survive the SVG sanitizer contract', () => {
  trustedRendererSemanticAttributes.forEach((attribute) => {
    assert.ok(
      SVG_SANITIZE_OPTIONS.ADD_ATTR.includes(attribute),
      `${attribute} is missing from the sanitizer allowlist`
    );
  });
  const rawSvg = [
    '<svg xmlns="http://www.w3.org/2000/svg">',
    '<g data-query-row="0" data-subject-row="1" ',
    'data-gbdraw-slot-id="comparison-1" data-record-key="record-a"></g>',
    '</svg>'
  ].join('');
  let receivedOptions = null;
  const sanitized = sanitizeSvgContent(rawSvg, {
    sanitize(source, options) {
      receivedOptions = options;
      return source;
    }
  });

  assert.equal(sanitized, rawSvg);
  assert.equal(receivedOptions, SVG_SANITIZE_OPTIONS);
});

test('standalone-only metadata is not admitted to plain preview sanitization', () => {
  [
    'data-gbdraw-interactive-feature',
    'data-gbdraw-interactive-match',
    'data-gbdraw-interactive-svg',
    'data-gbdraw-original-height',
    'data-gbdraw-original-viewbox',
    'data-gbdraw-original-width',
    'data-popup-mode',
    'data-result-index',
    'data-result-name',
    'data-schema'
  ].forEach((attribute) => {
    assert.equal(SVG_SANITIZE_OPTIONS.ADD_ATTR.includes(attribute), false);
  });
});

test('only the three schema 1 composition attributes are allowlisted', () => {
  assert.deepEqual(
    SVG_SANITIZE_OPTIONS.ADD_ATTR
      .filter((attribute) => attribute.startsWith('data-gbdraw-composition'))
      .sort(),
    [
      'data-gbdraw-composition',
      'data-gbdraw-composition-role',
      'data-gbdraw-composition-schema'
    ]
  );
});

test('retired alternate viewBox metadata is not allowlisted', () => {
  [
    'data-horizontal-viewbox',
    'data-vertical-viewbox'
  ].forEach((attribute) => {
    assert.equal(SVG_SANITIZE_OPTIONS.ADD_ATTR.includes(attribute), false);
  });
});

// The vendored DOMPurify keeps an attribute only when its built-in SVG lists,
// its aria-/data- rules, or ADD_ATTR admit it. This stub applies that rule to
// the attributes of the elements gbdraw writes, so the test needs no DOM.
const VENDORED_DOMPURIFY = readFileSync(
  new URL('../../gbdraw/web/vendor/dompurify/purify.min.js', import.meta.url),
  'utf8'
);
const domPurifyAttributeFilter = {
  /**
   * @param {string} source
   * @param {Readonly<Record<string, unknown>>} options
   */
  sanitize(source, options) {
    const added = /** @type {readonly string[]} */ (options.ADD_ATTR);
    return source.replace(/\s([\w:-]+)="[^"]*"/g, (attribute, name) => (
      added.includes(name)
      || /^(?:aria|data)-/.test(name)
      || VENDORED_DOMPURIFY.includes(`"${name}"`)
        ? attribute
        : ''
    ));
  }
};

test('label halos, underline casings and boxes keep pointer-events none (LP-D105, LP-D106)', () => {
  assert.equal(VENDORED_DOMPURIFY.includes('"pointer-events"'), false);
  const rawSvg = [
    '<svg xmlns="http://www.w3.org/2000/svg">',
    '<use href="#label-1" fill="none" stroke="#ffffff" stroke-width="1.92" ',
    'stroke-linejoin="round" stroke-linecap="round" aria-hidden="true" pointer-events="none"/>',
    '<line x1="0" y1="10" x2="40" y2="10" stroke="#ffffff" stroke-width="5.52" ',
    'stroke-linecap="butt" aria-hidden="true" pointer-events="none" data-label-feature-id="f1"/>',
    '<rect x="-4" y="-6" width="48" height="16" rx="4.32" ry="4.32" fill="#ffffff" fill-opacity="0.8" stroke="none" ',
    'aria-hidden="true" pointer-events="none" data-label-feature-id="f1"/>',
    '<text id="label-1" x="0" y="8" font-size="24">gene</text>',
    '</svg>'
  ].join('');
  const sanitized = sanitizeSvgContent(rawSvg, domPurifyAttributeFilter);

  assert.match(sanitized, /<use [^>]*pointer-events="none"/);
  assert.match(sanitized, /<line [^>]*pointer-events="none"/);
  assert.match(sanitized, /<rect [^>]*pointer-events="none"/);
});
