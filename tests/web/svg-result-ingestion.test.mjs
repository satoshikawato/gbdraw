import assert from 'node:assert/strict';
import test from 'node:test';

import {
  getCommittedSvgContent,
  getCommittedSvgResultMetadata,
  ingestCatalogBackedSvgResults,
  ingestSvgResult,
  ingestSvgResults,
  isCommittedSvgResult,
  markCommittedSvgResultMounted,
  markCommittedSvgResultUnmounted
} from '../../gbdraw/web/js/services/svg-result-ingestion.js';
import { normalizeLogicalResults } from '../../gbdraw/web/js/services/result-normalization.js';
class FakeDomParser {
  static calls = 0;

  parseFromString(content, mediaType) {
    FakeDomParser.calls += 1;
    assert.equal(mediaType, 'image/svg+xml');
    const feature = {
      getAttribute(name) {
        const attributes = {
          id: 'rendered-a',
          'data-gbdraw-feature-id': 'rendered-a',
          'data-gbdraw-stable-feature-id': 'stable-a',
          'data-gbdraw-record-index': '2',
          'data-gbdraw-record-id': 'record-a'
        };
        return attributes[name] ?? null;
      }
    };
    return {
      documentElement: {
        localName: content.includes('<svg') ? 'svg' : 'html',
        content,
        cloneNode() { return this; },
        getAttribute(name) {
          return name === 'xmlns' || name === 'xmlns:xlink' ? 'present' : null;
        },
        setAttribute() {},
        removeAttribute() {},
        querySelectorAll(selector) {
          return selector.includes('data-gbdraw-feature-id') ? [feature] : [];
        }
      },
      querySelector: () => null
    };
  }
}

globalThis.XMLSerializer = class {
  serializeToString(node) {
    return node.content;
  }
};

test('one ingestion transaction owns sanitization, parsing, and runtime trust', () => {
  let sanitizeCalls = 0;
  const sanitizer = {
    sanitize(content) {
      sanitizeCalls += 1;
      return content.replace(/<script>[\s\S]*?<\/script>/gi, '');
    }
  };
  const raw = {
    name: 'unsafe.svg',
    content: '<svg><script>unsafe()</script><path data-gbdraw-feature-id="rendered-a"/></svg>',
    trusted: true
  };

  const committed = ingestSvgResult(raw, {
    sanitizer,
    parser: FakeDomParser
  });

  assert.equal(sanitizeCalls, 1);
  assert.equal(FakeDomParser.calls, 1);
  assert.equal(isCommittedSvgResult(raw), false);
  assert.equal(isCommittedSvgResult(committed), true);
  assert.equal(committed.content, '<svg><path data-gbdraw-feature-id="rendered-a"/></svg>');
  assert.equal(getCommittedSvgContent(committed), committed.content);
  const identities = getCommittedSvgResultMetadata(committed).renderedFeatureIdentities;
  assert.deepEqual(identities.byRenderedId.get('rendered-a'), {
    renderedId: 'rendered-a',
    stableId: 'stable-a',
    recordIndex: 2,
    recordId: 'record-a',
    elementId: 'rendered-a'
  });

  assert.equal(markCommittedSvgResultMounted(committed), true);
  const persistedEdit = { ...committed, content: '<svg><path fill="#abcdef"/></svg>' };
  assert.equal(getCommittedSvgResultMetadata(persistedEdit), getCommittedSvgResultMetadata(committed));
  const normalizedCommitted = normalizeLogicalResults([committed])[0];
  assert.equal(getCommittedSvgResultMetadata(normalizedCommitted), getCommittedSvgResultMetadata(committed));
  assert.equal(
    getCommittedSvgContent(persistedEdit),
    '<svg><path data-gbdraw-feature-id="rendered-a"/></svg>',
    'persisting an active live-DOM edit must not replace the mounted root'
  );
  assert.equal(markCommittedSvgResultUnmounted(persistedEdit), true);
  assert.equal(getCommittedSvgContent(persistedEdit), persistedEdit.content);

  const persistedRoundTrip = JSON.parse(JSON.stringify(committed));
  assert.equal(isCommittedSvgResult(persistedRoundTrip), false);
  assert.equal(getCommittedSvgContent(persistedRoundTrip), null);
  assert.equal(getCommittedSvgResultMetadata(persistedRoundTrip), null);

  assert.equal(ingestSvgResults([committed], { sanitizer, parser: FakeDomParser })[0], committed);
  assert.equal(sanitizeCalls, 1);
  assert.equal(FakeDomParser.calls, 1);
});

test('a malformed sanitized document never becomes committed', () => {
  assert.throws(
    () => ingestSvgResult(
      { name: 'bad.svg', content: 'not svg' },
      {
        sanitizer: { sanitize: (content) => content },
        parser: FakeDomParser
      }
    ),
    /malformed SVG/
  );
});

test('a validated catalog commits sanitized saved SVG without a detached parse or scan', () => {
  const catalogItem = {
    recordKeys: ['record-a'],
    biologicalFeatures: [{
      recordKey: 'record-a',
      biologicalFeatureId: 'stable-a',
      record_id: 'public-record-a'
    }],
    features: [{
      recordKey: 'record-a',
      biologicalFeatureId: 'stable-a',
      svgId: 'rendered-a'
    }]
  };
  let sanitizeCalls = 0;
  const beforeParserCalls = FakeDomParser.calls;
  const committed = ingestCatalogBackedSvgResults([{
    name: 'saved.svg',
    content: '<?xml version="1.0"?><svg><script>bad()</script><path id="rendered-a"/></svg>'
  }], {
    catalogItems: [catalogItem],
    sanitizer: {
      sanitize(content) {
        sanitizeCalls += 1;
        return content.replace(/^<\?xml[^>]*>/, '').replace(/<script>[\s\S]*?<\/script>/, '');
      }
    }
  })[0];

  assert.equal(sanitizeCalls, 1);
  assert.equal(FakeDomParser.calls, beforeParserCalls);
  assert.equal(isCommittedSvgResult(committed), true);
  assert.equal(committed.content, '<svg><path id="rendered-a"/></svg>');
  assert.deepEqual(
    getCommittedSvgResultMetadata(committed).renderedFeatureIdentities.byRenderedId
      .get('rendered-a'),
    {
      renderedId: 'rendered-a',
      stableId: 'stable-a',
      recordIndex: 0,
      recordId: 'public-record-a',
      elementId: 'rendered-a'
    }
  );

  assert.throws(
    () => ingestCatalogBackedSvgResults(
      [{ name: 'bad.svg', content: 'not svg' }],
      {
        catalogItems: [catalogItem],
        sanitizer: { sanitize: (content) => content }
      }
    ),
    /malformed SVG/
  );
});
