import assert from 'node:assert/strict';

import {
  EXCLUDED_GROUP_SELECTOR
} from '../../gbdraw/web/js/app/feature-editor/label-actions.js';

assert.match(
  EXCLUDED_GROUP_SELECTOR,
  /g\[data-gbdraw-slot-renderer="ticks"\]/,
  'custom Circular tick-slot TextPaths must not be exposed as editable feature labels'
);
assert.match(EXCLUDED_GROUP_SELECTOR, /g\[id="tick"\]/);
assert.match(EXCLUDED_GROUP_SELECTOR, /g\[id\^="tick_"\]/);
