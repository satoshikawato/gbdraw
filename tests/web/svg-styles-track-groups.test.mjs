import assert from 'node:assert/strict';

import { getGroupsByBaseIds } from '../../gbdraw/web/js/app/svg-styles.js';

const makeGroup = (id, renderer = null) => ({
  getAttribute(name) {
    if (name === 'id') return id;
    if (name === 'data-gbdraw-slot-renderer') return renderer;
    return null;
  }
});

const semanticSkew = makeGroup('track_slot_arbitrary_hash', 'dinucleotide_skew');
const semanticLegacySkew = makeGroup('gc_skew_semantic', 'dinucleotide_skew');
const legacySkew = makeGroup('gc_skew_2');
const semanticContent = makeGroup('track_slot_content_hash', 'dinucleotide_content');
const semanticDepth = makeGroup('track_slot_depth_hash', 'depth');
const semanticTicks = makeGroup('track_slot_ticks_hash', 'ticks');
const groups = [
  semanticSkew,
  semanticLegacySkew,
  legacySkew,
  semanticContent,
  semanticDepth,
  semanticTicks
];

const svg = {
  querySelectorAll(selector) {
    const rendererMatch = selector.match(
      /^g\[data-gbdraw-slot-renderer="([^"]+)"\]$/
    );
    if (rendererMatch) {
      return groups.filter(
        (group) => group.getAttribute('data-gbdraw-slot-renderer') === rendererMatch[1]
      );
    }

    const idMatchers = Array.from(selector.matchAll(/g\[id(\^)?="([^"]+)"\]/g));
    return groups.filter((group) => {
      const id = group.getAttribute('id');
      return idMatchers.some(([, prefix, value]) => (
        prefix ? id.startsWith(value) : id === value
      ));
    });
  }
};

assert.deepEqual(
  getGroupsByBaseIds(
    svg,
    ['skew', 'gc_skew'],
    ['dinucleotide_skew']
  ),
  [semanticSkew, semanticLegacySkew],
  'semantic slot groups should suppress persisted-ID fallback matches'
);

assert.deepEqual(
  getGroupsByBaseIds(svg, ['gc_content'], ['dinucleotide_content']),
  [semanticContent]
);
assert.deepEqual(
  getGroupsByBaseIds(svg, ['depth'], ['depth']),
  [semanticDepth]
);
assert.deepEqual(
  getGroupsByBaseIds(svg, ['tick'], ['ticks']),
  [semanticTicks]
);

const persistedSvg = {
  querySelectorAll(selector) {
    if (selector.startsWith('g[data-gbdraw-slot-renderer=')) return [];
    return selector.includes('gc_skew') ? [legacySkew] : [];
  }
};
assert.deepEqual(
  getGroupsByBaseIds(
    persistedSvg,
    ['skew', 'gc_skew'],
    ['dinucleotide_skew']
  ),
  [legacySkew],
  'legacy IDs should remain a fallback for persisted SVGs without semantic hooks'
);
