import assert from 'node:assert/strict';

const {
  createDefaultLayoutPreferences,
  migrateLegacyLayoutPreferences,
  normalizeLayoutPreferences,
  resolveActiveLayoutPreference,
  resolveCircularLayoutPreference,
  requestLayoutPreferences,
  updateActiveLayoutPreference
} = await import('../../gbdraw/web/js/services/layout-preferences.js');
const { ABSENT_CIRCULAR_OUTPUT, ABSENT_LINEAR_OUTPUT } = await import('../../gbdraw/web/js/contracts/schema.generated.js');

// A request without `diagramOptions.output` (or with a null title side) reads
// with the sides Python draws it with, generated from its output options: a
// Circular one on one record and on a multi-record canvas alike.
for (const [mode, grouping, slot, absent] of [
  ['circular', 'single', (layout) => layout.circular.single, ABSENT_CIRCULAR_OUTPUT],
  ['circular', 'grid', (layout) => layout.circular.multi, ABSENT_CIRCULAR_OUTPUT],
  ['linear', 'single', (layout) => layout.linear, ABSENT_LINEAR_OUTPUT]
]) {
  assert.ok(absent, `${mode} absent output is generated`);
  assert.deepEqual(slot(requestLayoutPreferences({ mode, grouping, diagramOptions: {} })), absent, `${mode} ${grouping}`);
  assert.deepEqual(
    slot(requestLayoutPreferences({ mode, grouping, diagramOptions: { output: { legend: 'left', plotTitlePosition: null } } })),
    { legend: 'left', plotTitlePosition: absent.plotTitlePosition },
    `${mode} ${grouping}, null title side`
  );
}

const defaults = createDefaultLayoutPreferences();
assert.deepEqual(resolveCircularLayoutPreference(defaults, false), {
  legend: 'left',
  plotTitlePosition: 'none'
});
assert.deepEqual(resolveCircularLayoutPreference(defaults, true), {
  legend: 'left',
  plotTitlePosition: 'none'
});
assert.deepEqual(resolveActiveLayoutPreference(defaults, 'linear'), {
  legend: 'bottom',
  plotTitlePosition: 'bottom'
});
updateActiveLayoutPreference(defaults, 'circular', false, {
  legend: 'right',
  plotTitlePosition: 'top'
});
updateActiveLayoutPreference(defaults, 'linear', false, {
  legend: 'top',
  plotTitlePosition: 'center'
});
assert.deepEqual(defaults.circular.single, {
  legend: 'right',
  plotTitlePosition: 'top'
});
assert.deepEqual(defaults.linear, {
  legend: 'top',
  plotTitlePosition: 'center'
});

const normalized = normalizeLayoutPreferences({
  circular: {
    single: { legend: 'right', plotTitlePosition: 'top' },
    multi: { legend: '', plotTitlePosition: 'bottom' }
  },
  linear: { legend: 'top', plotTitlePosition: 'center' }
});
assert.deepEqual(normalized, {
  circular: {
    single: { legend: 'right', plotTitlePosition: 'top' },
    multi: { legend: null, plotTitlePosition: 'bottom' }
  },
  linear: { legend: 'top', plotTitlePosition: 'center' }
});
assert.deepEqual(resolveCircularLayoutPreference(normalized, true), {
  legend: 'right',
  plotTitlePosition: 'bottom'
});

const migrated = migrateLegacyLayoutPreferences(
  {
    circularSingleRecordLegendPosition: 'right',
    circularSingleRecordPlotTitlePosition: 'top',
    circularMultiRecordLegendPosition: 'bottom',
    circularMultiRecordPlotTitlePosition: 'bottom',
    linearLegendPosition: 'left',
    linearPlotTitlePosition: 'center'
  }
);
assert.deepEqual(migrated, {
  circular: {
    single: { legend: 'right', plotTitlePosition: 'top' },
    multi: { legend: 'bottom', plotTitlePosition: 'bottom' }
  },
  linear: { legend: 'left', plotTitlePosition: 'center' }
});

// A drawing's settings slice installs normalized preferences into its own objects.
const { SETTINGS_SLICE } = await import('../../gbdraw/web/js/services/drawing/settings.js');
const target = createDefaultLayoutPreferences();
const { single } = target.circular;
const settings = SETTINGS_SLICE.create({ layoutPreferences: target });
SETTINGS_SLICE.install(settings.handle, { layoutPreferences: normalizeLayoutPreferences(migrated) });
assert.deepEqual(target, migrated);
assert.equal(target.circular.single, single, 'the nested preferences stay the same objects');

// The Multi-Record Canvas toggle's smart defaults (A2 W1): a left Legend and
// no title move below the grid, other sides are copied, and a stored
// multi-record side means the toggle stored them already.
const { circularMultiRecordDefaults } = await import('../../gbdraw/web/js/services/layout-preferences.js');
assert.deepEqual(circularMultiRecordDefaults(createDefaultLayoutPreferences()), { legend: 'bottom', plotTitlePosition: 'bottom' });
const sided = createDefaultLayoutPreferences();
Object.assign(sided.circular.single, { legend: 'right', plotTitlePosition: 'top' });
assert.deepEqual(circularMultiRecordDefaults(sided), { legend: 'right', plotTitlePosition: 'top' });
sided.circular.multi.plotTitlePosition = 'none';
assert.equal(circularMultiRecordDefaults(sided), null);

console.log('layout preference tests passed');
