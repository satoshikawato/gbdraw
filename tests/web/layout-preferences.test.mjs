import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-layout-preferences-'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempRoot, 'app'));
for (const filename of ['layout-preferences.js', 'plot-title-position.js']) {
  await writeFile(
    join(tempRoot, 'app', filename),
    await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'app', filename), 'utf8'),
    'utf8'
  );
}
await writeFile(
  join(tempRoot, 'web-ux-profile.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'web-ux-profile.js'), 'utf8'),
  'utf8'
);

const {
  createDefaultLayoutPreferences,
  migrateLegacyLayoutPreferences,
  normalizeLayoutPreferences,
  replaceLayoutPreferences,
  resolveActiveLayoutPreference,
  resolveCircularLayoutPreference,
  updateActiveLayoutPreference
} = await import(pathToFileURL(join(tempRoot, 'app', 'layout-preferences.js')));

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

const target = createDefaultLayoutPreferences();
replaceLayoutPreferences(target, migrated);
assert.deepEqual(target, migrated);

console.log('layout preference tests passed');
