import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const sourceUrl = new URL('../../gbdraw/web/js/app/depth-track-state.js', import.meta.url);
const displaySourceUrl = new URL('../../gbdraw/web/js/app/track-slot-display.js', import.meta.url);
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-depth-track-state-'));
const tempModulePath = join(tempDir, 'depth-track-state.mjs');
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n');
await writeFile(tempModulePath, await readFile(sourceUrl, 'utf8'));
await writeFile(join(tempDir, 'track-slot-display.js'), await readFile(displaySourceUrl, 'utf8'));

const {
  depthFileSlotsFromValue,
  reindexDepthSlots,
  reconcileDepthTracksToFiles,
  removeDepthTrackAt,
  syncDepthSlotLabels
} = await import(pathToFileURL(tempModulePath));

const file = (name) => ({ name });
const labels = (tracks) => tracks.map((track) => track.label);
const depthIndexes = (slots) => slots
  .filter((slot) => slot.renderer === 'depth' && slot.enabled !== false)
  .map((slot) => slot.params.track_index);
const depthLabels = (slots) => slots
  .filter((slot) => slot.renderer === 'depth' && slot.enabled !== false)
  .map((slot) => slot.params.legend_label);
const depthSlots = (count) => Array.from({ length: count }, (_, index) => ({
  id: `depth_${index + 1}`,
  renderer: 'depth',
  params: {
    track_index: index,
    legend_label: ['24 hpi', '12 hpi', '6 hpi'][index]
  }
}));

{
  const removal = removeDepthTrackAt({
    files: [file('24 hpi.tsv'), file('12 hpi.tsv'), file('6 hpi.tsv')],
    depthTracks: [{ label: '24 hpi' }, { label: '12 hpi' }, { label: '6 hpi' }],
    index: 1
  });
  const slots = reindexDepthSlots({
    slots: depthSlots(3),
    removedIndex: 1,
    activeCount: depthFileSlotsFromValue(removal.files).length
  });
  syncDepthSlotLabels({ slots, depthTracks: removal.depthTracks, activeCount: 2 });

  assert.deepEqual(removal.files.map((entry) => entry.name), ['24 hpi.tsv', '6 hpi.tsv']);
  assert.deepEqual(labels(removal.depthTracks), ['24 hpi', '6 hpi']);
  assert.deepEqual(depthIndexes(slots), [0, 1]);
  assert.deepEqual(depthLabels(slots), ['24 hpi', '6 hpi']);
}

{
  const removal = removeDepthTrackAt({
    files: [file('24 hpi.tsv'), file('12 hpi.tsv'), file('6 hpi.tsv')],
    depthTracks: [{ label: '24 hpi' }, { label: '12 hpi' }, { label: '6 hpi' }],
    index: 0
  });
  const slots = reindexDepthSlots({
    slots: depthSlots(3),
    removedIndex: 0,
    activeCount: 2
  });
  syncDepthSlotLabels({ slots, depthTracks: removal.depthTracks, activeCount: 2 });

  assert.deepEqual(removal.files.map((entry) => entry.name), ['12 hpi.tsv', '6 hpi.tsv']);
  assert.deepEqual(labels(removal.depthTracks), ['12 hpi', '6 hpi']);
  assert.deepEqual(depthIndexes(slots), [0, 1]);
  assert.deepEqual(depthLabels(slots), ['12 hpi', '6 hpi']);
}

{
  const removal = removeDepthTrackAt({
    files: [file('24 hpi.tsv'), file('12 hpi.tsv'), file('6 hpi.tsv')],
    depthTracks: [{ label: '24 hpi' }, { label: '12 hpi' }, { label: '6 hpi' }],
    index: 2
  });
  const slots = reindexDepthSlots({
    slots: depthSlots(3),
    removedIndex: 2,
    activeCount: 2
  });
  syncDepthSlotLabels({ slots, depthTracks: removal.depthTracks, activeCount: 2 });

  assert.deepEqual(removal.files.map((entry) => entry.name), ['24 hpi.tsv', '12 hpi.tsv']);
  assert.deepEqual(labels(removal.depthTracks), ['24 hpi', '12 hpi']);
  assert.deepEqual(depthIndexes(slots), [0, 1]);
}

{
  const removal = removeDepthTrackAt({
    files: [file('24 hpi.tsv')],
    depthTracks: [{ label: '24 hpi' }],
    index: 0
  });
  const slots = reindexDepthSlots({
    slots: depthSlots(1),
    removedIndex: 0,
    activeCount: 0
  });

  assert.deepEqual(removal.files, []);
  assert.deepEqual(depthIndexes(slots), []);
  assert.equal(removal.depthTracks.length, 1);
}

{
  const slots = reindexDepthSlots({
    slots: [
      { id: 'depth_1', renderer: 'depth', params: { track_index: 0 } },
      { id: 'custom_removed', renderer: 'depth', width: '10px', params: { track_index: 1 } },
      { id: 'depth_3', renderer: 'depth', params: { track_index: 2 } }
    ],
    removedIndex: 1,
    activeCount: 2
  });

  const custom = slots.find((slot) => slot.id === 'custom_removed');
  assert.equal(custom.enabled, false);
  assert.equal(Object.hasOwn(custom.params, 'track_index'), false);
  assert.deepEqual(depthIndexes(slots), [0, 1]);
}

{
  const reconciled = reconcileDepthTracksToFiles({
    files: [
      file('UAZ00-173B.SRR14027705.minimap2.depth.tsv'),
      file('UAZ00-173B.SRR14027721.minimap2.depth.tsv')
    ],
    depthTracks: [
      { label: '24 hpi SRR14027705)' },
      { label: '6 hpi (SRR14027712)' },
      { label: '6 hpi (SRR14027721)' }
    ],
    targetCount: 2
  });

  assert.deepEqual(labels(reconciled), ['24 hpi SRR14027705)', '6 hpi (SRR14027721)']);
}

console.log(`depth-track-state tests passed (${fileURLToPath(sourceUrl)})`);
