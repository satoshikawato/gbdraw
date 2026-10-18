import assert from 'node:assert/strict';

const {
  createDefaultLayoutPreferences,
  migrateLegacyLayoutPreferences,
  normalizeLayoutPreferences,
  resolveActiveLayoutPreference,
  resolveCircularLayoutPreference,
  updateActiveLayoutPreference
} = await import('../../gbdraw/web/js/services/layout-preferences.js');

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

console.log('layout preference tests passed');
