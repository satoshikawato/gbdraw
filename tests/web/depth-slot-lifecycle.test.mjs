// Managed Depth row lifecycle (PD-OI-058, TR-02/TR-03, rule R10): one table
// runs through the Circular and the Linear editor. A Depth source change is the
// only transition that adds or removes managed rows; manual rows are never
// changed and report a row issue when their series has no source (PD-OI-083).
import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-depth-slot-lifecycle-'));
for (const directory of ['app', 'utils']) {
  await cp(join(repoRoot, 'gbdraw', 'web', 'js', directory), join(tempRoot, directory), { recursive: true });
}
// circular-track-slots.js reads the conservation-series row helpers.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'conservation-series.js'),
  join(tempRoot, 'services', 'conservation-series.js')
);
// linear-track-slots.js reads the current option values.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'current-option-values.js'),
  join(tempRoot, 'services', 'current-option-values.js')
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
// The track-slot leaves the slot editors import.
for (const leaf of ['depth-track-state.js', 'track-slot-display.js', 'track-slot-validation.js',
  'circular-track-slot-model.js', 'linear-track-slot-model.js']) {
  await cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'services', leaf),
    join(tempRoot, 'services', leaf)
  );
}

const { createCircularTrackSlotEditor } = await import(pathToFileURL(join(tempRoot, 'app', 'circular-track-slots.js')));
const { createLinearTrackSlotEditor } = await import(pathToFileURL(join(tempRoot, 'app', 'linear-track-slots.js')));

const file = (trackIndex) => ({ name: `series-${trackIndex}.depth.tsv`, size: 1, lastModified: trackIndex });
// One record; `sourced` lists the logical series indexes that have a file.
const sourceRow = (sourced, width = Math.max(0, ...sourced.map((index) => index + 1))) => (
  Array.from({ length: width }, (_, index) => (sourced.includes(index) ? file(index) : null))
);
// `manual` rows carry a user parameter, so they are not default managed rows.
const depthRow = (side, { id, track, enabled = true, manual = false }) => ({
  id,
  renderer: 'depth',
  enabled,
  side,
  z: 0,
  params: { track_index: track, ...(manual ? { custom: 'keep' } : {}) }
});

const MODES = {
  circular: {
    create: (rows, sourced, width) => {
      const state = {
        mode: { value: 'circular' },
        form: { track_type: 'tuckin', show_depth: true, suppress_gc: false, suppress_skew: true, show_scale: true, separate_strands: true },
        adv: {
          circular_track_slots_enabled: true,
          circular_track_slots_axis_index: 1,
          circular_track_slots: [
            { id: 'features', renderer: 'features', enabled: true, side: 'outside', z: 0, params: { lane_direction: 'outside' } },
            { id: 'ticks', renderer: 'ticks', enabled: true, side: 'inside', z: 0, params: { tick_label_layout: 'label_out_tick_in' } },
            ...rows.map((row) => depthRow('inside', row)),
            { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'inside', z: 0, params: { nt: 'GC' } }
          ],
          nt: 'GC',
          features: ['CDS'],
          feature_shapes: { CDS: 'arrow' },
          depth_tracks: [],
          feature_width_circular: null,
          depth_width_circular: null,
          gc_content_width_circular: null,
          gc_content_radius_circular: null,
          gc_skew_width_circular: null,
          gc_skew_radius_circular: null
        },
        files: { c_depth: sourceRow(sourced, width), c_conservation_blasts: [], c_conservation_fastas: [] },
        circularConservation: { enabled: false, source: 'upload', labels: '', series: [] },
        annotationSets: [],
        circularRecordList: { value: [] }
      };
      const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
      editor.normalizeCircularTrackSlots();
      return {
        state,
        slots: () => state.adv.circular_track_slots,
        change: (nextSourced, nextWidth) => editor.changeCircularDepthSources(() => {
          state.files.c_depth = sourceRow(nextSourced, nextWidth);
        }),
        issue: (index) => editor.circularTrackSlotIssue(state.adv.circular_track_slots[index], index)
      };
    }
  },
  linear: {
    create: (rows, sourced, width) => {
      const state = {
        form: { linear_track_layout: 'middle', show_depth: true, show_gc: true, show_skew: false },
        adv: {
          nt: 'GC',
          linear_track_slots_enabled: true,
          linear_track_slots_axis_index: 0,
          linear_track_slots: [
            { id: 'features', renderer: 'features', enabled: true, side: 'overlay', height: '', spacing: '', z: 0, params: {} },
            ...rows.map((row) => ({ ...depthRow('below', row), height: '', spacing: '' })),
            { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'below', height: '', spacing: '', z: 0, params: { nt: 'GC' } }
          ],
          depth_tracks: [],
          depth_height: null,
          gc_height: null
        },
        linearSeqs: [{ depth: sourceRow(sourced, width) }],
        annotationSets: []
      };
      const editor = createLinearTrackSlotEditor({ state: withDrawings(state) });
      editor.normalizeLinearTrackSlots();
      return {
        state,
        slots: () => state.adv.linear_track_slots,
        change: (nextSourced, nextWidth) => editor.changeLinearDepthSources(() => {
          // PD-OI-025: clearing a source keeps the logical series width.
          state.linearSeqs[0].depth = sourceRow(nextSourced, nextWidth);
        }),
        issue: (index) => editor.linearTrackSlotIssue(state.adv.linear_track_slots[index], index)
      };
    }
  }
};

const depthRows = (slots) => slots
  .filter((slot) => slot.renderer === 'depth')
  .map((slot) => `${slot.id}:${slot.params.track_index}:${slot.enabled === false ? 'off' : 'on'}${slot.params.custom ? ':manual' : ''}`);
const otherRows = (slots) => slots.filter((slot) => slot.renderer !== 'depth').map((slot) => JSON.stringify(slot));

const CASES = [
  {
    name: 'a series that gains its first source gets one managed row',
    rows: [], before: [], after: [0],
    expected: ['depth:0:on']
  },
  {
    name: 'a disabled row on the series blocks the managed row and stays disabled',
    rows: [{ id: 'depth', track: 0, enabled: false }], before: [], after: [0],
    expected: ['depth:0:off']
  },
  {
    name: 'a manual row on the series blocks the managed row',
    rows: [{ id: 'my_depth', track: 0, manual: true }], before: [], after: [0],
    expected: ['my_depth:0:on:manual']
  },
  {
    name: 'replacing a source does not restore a deleted row',
    rows: [], before: [0], after: [0],
    expected: []
  },
  {
    name: 'a series that loses its last source loses its managed rows only',
    rows: [{ id: 'depth', track: 0 }, { id: 'my_depth', track: 0, manual: true }], before: [0], after: [], width: 1,
    expected: ['my_depth:0:on:manual'],
    issueFor: 'my_depth'
  },
  {
    name: 'only a gained series without a referencing row gets a row',
    rows: [{ id: 'depth_1', track: 0 }, { id: 'depth_3', track: 2, enabled: false }], before: [0], after: [0, 1, 2],
    expected: ['depth_1:0:on', 'depth_3:2:off', 'depth_2:1:on']
  },
  {
    name: 'only the series that loses its last source loses its managed row',
    rows: [{ id: 'depth_1', track: 0 }, { id: 'depth_2', track: 1 }], before: [0, 1], after: [1], width: 2,
    expected: ['depth_2:1:on']
  }
];

for (const [mode, adapter] of Object.entries(MODES)) {
  for (const testCase of CASES) {
    test(`${mode}: ${testCase.name}`, () => {
      const width = testCase.width ?? null;
      const editor = adapter.create(testCase.rows, testCase.before, width ?? undefined);
      const before = otherRows(editor.slots());
      editor.change(testCase.after, width ?? undefined);
      assert.deepEqual(depthRows(editor.slots()), testCase.expected);
      assert.deepEqual(otherRows(editor.slots()), before, 'non-Depth rows keep their order and contents');
      if (testCase.issueFor) {
        const index = editor.slots().findIndex((slot) => slot.id === testCase.issueFor);
        assert.match(editor.issue(index), /Depth track 'my_depth' has no logical Depth source\./);
      }
    });
  }
}
