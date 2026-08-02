import assert from 'node:assert/strict';
import {
  normalizeLogicalResults
} from '../../gbdraw/web/js/services/result-normalization.js';

assert.deepEqual(
  normalizeLogicalResults([
    { name: 'layout.svg', content: 'plain' },
    { name: 'layout.interactive.svg', content: 'interactive' }
  ]),
  [{ name: 'layout.svg', content: 'plain' }]
);

assert.deepEqual(
  normalizeLogicalResults([
    { name: 'layout.interactive.svg', content: 'interactive' },
    { name: 'layout.svg', content: 'plain' },
    { name: 'layout_2.svg', content: 'batch-2' }
  ]),
  [
    { name: 'layout.svg', content: 'plain' },
    { name: 'layout_2.svg', content: 'batch-2' }
  ]
);

assert.deepEqual(
  normalizeLogicalResults([
    { name: 'interactive-map.svg', content: 'distinct' },
    { name: 'interactive-map.interactive.svg', content: 'paired' },
    { name: 'only.interactive.svg', content: 'legacy-only' }
  ]),
  [
    { name: 'interactive-map.svg', content: 'distinct' },
    { name: 'only.interactive.svg', content: 'legacy-only' }
  ]
);

assert.deepEqual(normalizeLogicalResults(null), []);
