import assert from 'node:assert/strict';
import test from 'node:test';
import {
  buildLegendStyleRetirement,
  trackDataLegendCaptions
} from '../../gbdraw/web/js/app/legend/track-data-styles.js';

test('track data names the legend labels of annotation sets and the sourced Depth series (OV-65)', () => {
  const captions = trackDataLegendCaptions({
    annotationSets: [
      { legendLabel: ' Set ', annotations: [{ legendLabel: 'Region X' }, { legendLabel: null }] },
      { legendLabel: '', annotations: [] }
    ],
    depthTracks: [{ label: 'sample' }, { label: 'other' }],
    depthSlots: [
      { renderer: 'depth', params: { track_index: 0, legend_label: 'Coverage' } },
      { renderer: 'depth', params: { track_index: 1, legend_label: 'Unsourced' } },
      { renderer: 'features', params: { legend_label: 'Not Depth' } }
    ],
    sourcedDepthTrackIndexes: [0]
  });
  assert.deepEqual([...captions].sort(), ['Coverage', 'Region X', 'Set', 'sample']);
});

test('a Depth series without a label is named by the fallback label', () => {
  assert.deepEqual([...trackDataLegendCaptions({ sourcedDepthTrackIndexes: [0, 1] })], ['Depth', 'Depth 2']);
});

test('styles follow the caption: a data change retires only the captions the data stops naming', () => {
  const sets = [{ legendLabel: '', annotations: [{ legendLabel: 'Keep' }, { legendLabel: 'Drop' }] }];
  const legendColorOverrides = { Keep: '#111111', Drop: '#222222', Feature: '#333333' };
  const legendStrokeOverrides = { Drop: { strokeColor: '#000000' }, Feature: { strokeWidth: 2 } };
  const retire = buildLegendStyleRetirement({
    legendColorOverrides,
    legendStrokeOverrides,
    namedCaptions: () => trackDataLegendCaptions({ annotationSets: sets })
  });
  const result = retire(() => {
    sets[0].annotations.pop();
    return 'changed';
  });
  assert.equal(result, 'changed');
  assert.deepEqual(legendColorOverrides, { Keep: '#111111', Feature: '#333333' });
  assert.deepEqual(legendStrokeOverrides, { Feature: { strokeWidth: 2 } });
});

test('a change that keeps every caption retires nothing', () => {
  const legendColorOverrides = { Keep: '#111111' };
  const retire = buildLegendStyleRetirement({
    legendColorOverrides,
    legendStrokeOverrides: {},
    namedCaptions: () => new Set(['Keep'])
  });
  retire(() => undefined);
  assert.deepEqual(legendColorOverrides, { Keep: '#111111' });
});
