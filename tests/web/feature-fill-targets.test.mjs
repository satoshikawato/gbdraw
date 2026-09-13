import assert from 'node:assert/strict';
import { test } from 'node:test';
import {
  FEATURE_PART_BLOCK, filterFeatureFillTargets, getFeaturePart
} from '../../gbdraw/web/js/app/feature-dom.js';

const element = (id, part, fill) => ({
  getAttribute: key => ({ id, 'data-gbdraw-feature-part': part, fill })[key] ?? null
});

test('fill edits exclude outline-only and connector paths while preserving transparent blocks', () => {
  const block = element('f1__part1', 'block', 'none');
  const outline = element('f1__part1__outline', 'block', 'none');
  const connector = element('f1__line1', 'connector', 'none');
  assert.deepEqual(filterFeatureFillTargets([block, outline, connector]), [block]);
  assert.equal(getFeaturePart(outline), FEATURE_PART_BLOCK, 'outlines remain targets for block stroke edits');
});
