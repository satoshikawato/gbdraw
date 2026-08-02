import assert from 'node:assert/strict';
import test from 'node:test';

import {
  SVG_SANITIZE_OPTIONS,
  sanitizeSvgContent
} from '../../gbdraw/web/js/services/svg-sanitization.js';

const trustedRendererSemanticAttributes = [
  'data-edge-kind',
  'data-gbdraw-annotation-clipped',
  'data-gbdraw-definition-part',
  'data-gbdraw-feature-offset-y',
  'data-gbdraw-slot-id',
  'data-gbdraw-slot-renderer',
  'data-horizontal-viewbox',
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
  'data-subject-stable-feature-svg-id',
  'data-vertical-viewbox'
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
