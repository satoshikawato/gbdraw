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
    legendEntries: { value: [] },
    projectLegendEntries: () => assert.fail('no rename to retire'),
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
  const legendEntries = { value: [{ caption: 'Kept name', originalCaption: 'Keep' }] };
  const retire = buildLegendStyleRetirement({
    legendColorOverrides,
    legendStrokeOverrides: {},
    legendEntries,
    projectLegendEntries: () => assert.fail('no rename to retire'),
    namedCaptions: () => new Set(['Keep'])
  });
  retire(() => undefined);
  assert.deepEqual(legendColorOverrides, { Keep: '#111111' });
  assert.deepEqual(legendEntries.value, [{ caption: 'Kept name', originalCaption: 'Keep' }]);
});

// OV-87: names follow the caption as styles do. The rename of a row whose caption
// the data stops naming is retired with the styles stored under the new name.
test('a data change retires the Legend rename of a caption it stops naming, and the styles of the new name', () => {
  const named = new Set(['depth', 'Region X']);
  const legendColorOverrides = { Coverage: '#7b2cbf', 'Region Y': '#111111', CDS: '#222222' };
  const legendStrokeOverrides = { Coverage: { strokeWidth: 2 } };
  const entries = [
    { caption: 'CDS', originalCaption: 'CDS', color: '#222222' },
    { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf' },
    { caption: 'Region Y', originalCaption: 'Region X', color: '#111111' }
  ];
  const legendEntries = { value: entries };
  let projections = 0;
  const retire = buildLegendStyleRetirement({
    legendColorOverrides,
    legendStrokeOverrides,
    legendEntries,
    projectLegendEntries: () => { projections += 1; },
    namedCaptions: () => new Set(named)
  });
  retire(() => named.delete('depth'));
  assert.deepEqual(legendEntries.value.map(({ originalCaption, caption }) => `${originalCaption}=>${caption}`),
    ['CDS=>CDS', 'depth=>depth', 'Region X=>Region Y']);
  assert.equal(legendEntries.value[1].color, '#7b2cbf', 'the entry keeps what the displayed row shows');
  assert.deepEqual(legendColorOverrides, { 'Region Y': '#111111', CDS: '#222222' });
  assert.deepEqual(legendStrokeOverrides, {});
  assert.equal(projections, 1, 'the displayed Result shows the retired rename');
  assert.equal(entries[1].caption, 'Coverage', 'the History step keeps the entries it captured');
});

test('a new name that the new data names keeps its styles when the rename is retired', () => {
  const named = new Set(['depth']);
  const legendColorOverrides = { Coverage: '#7b2cbf' };
  const legendEntries = { value: [{ caption: 'Coverage', originalCaption: 'depth' }] };
  const retire = buildLegendStyleRetirement({
    legendColorOverrides,
    legendStrokeOverrides: {},
    legendEntries,
    projectLegendEntries: () => {},
    namedCaptions: () => new Set(named)
  });
  retire(() => { named.delete('depth'); named.add('Coverage'); });
  assert.deepEqual(legendEntries.value, [{ caption: 'depth', originalCaption: 'depth' }]);
  assert.deepEqual(legendColorOverrides, { Coverage: '#7b2cbf' });
});

// OV-120: a renamed row a Generate hid (Show Depth off) waits in the drawing;
// removing its data retires the waiting rename and its styles like a shown one.
test('a data change retires a waiting Legend rename of a caption it stops naming (OV-120)', () => {
  const named = new Set(['depth', 'Region X']);
  const legendColorOverrides = { Coverage: '#7b2cbf', 'Region Y': '#111111' };
  const legendStrokeOverrides = { Coverage: { strokeWidth: 2 } };
  const dormantLegendEntries = {
    value: [
      { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf' },
      { caption: 'Region Y', originalCaption: 'Region X', color: '#111111' }
    ]
  };
  const retire = buildLegendStyleRetirement({
    legendColorOverrides,
    legendStrokeOverrides,
    legendEntries: { value: [] },
    dormantLegendEntries,
    projectLegendEntries: () => assert.fail('no shown rename to retire'),
    namedCaptions: () => new Set(named)
  });
  retire(() => named.delete('depth'));
  assert.deepEqual(dormantLegendEntries.value, [{ caption: 'Region Y', originalCaption: 'Region X', color: '#111111' }]);
  assert.deepEqual(legendColorOverrides, { 'Region Y': '#111111' });
  assert.deepEqual(legendStrokeOverrides, {});
});
