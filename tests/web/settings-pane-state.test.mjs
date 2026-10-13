// The display of the settings pane (services/settings-pane-state.js): the
// stacked settings or the icon rail, the rail entry whose sections are docked
// beside it, and the section in the popup beside the Result. Every transition
// is here; the display is 'stacked' whenever nothing is docked by the rail, a
// docked entry is always one the shown mode has, and a section is never in
// the dock and in the popup.
import assert from 'node:assert/strict';
import test from 'node:test';

import {
  SETTINGS_PANE_DISPLAYS, closeDockedSettingsEntry, dockSettingsSection, floatSettingsSection, moveSettingsPopup,
  reconcileSettingsPaneState, returnSettingsSection, showSettingsRail, showSettingsStacked
} from '../../gbdraw/web/js/services/settings-pane-state.js';

const paneState = (display = 'stacked', docked = null) => ({
  settingsPaneDisplay: { value: display },
  settingsPaneDocked: { value: docked },
  settingsPopup: { value: null }
});
const read = (state) => [state.settingsPaneDisplay.value, state.settingsPaneDocked.value];

test('the displays are stacked, the default, and rail', () => {
  assert.deepEqual([...SETTINGS_PANE_DISPLAYS], ['stacked', 'rail']);
});

test('switching to the rail docks the given entry, or Input when the mode lacks it', () => {
  const state = paneState();
  assert.equal(showSettingsRail(state, 'linear', 'labels'), 'labels');
  assert.deepEqual(read(state), ['rail', 'labels']);
  assert.equal(showSettingsRail(state, 'circular', 'comparison'), 'input-genomes');
  assert.deepEqual(read(state), ['rail', 'input-genomes']);
  // A page that opens in the remembered rail names no entry.
  const opened = paneState();
  assert.equal(showSettingsRail(opened, 'circular'), 'input-genomes');
  assert.deepEqual(read(opened), ['rail', 'input-genomes']);
});

test('switching to the stacked display docks nothing', () => {
  const state = paneState('rail', 'labels');
  showSettingsStacked(state);
  assert.deepEqual(read(state), ['stacked', null]);
});

test('docking a section docks the rail entry that shows it and keeps a docked entry docked', () => {
  const state = paneState('rail', null);
  assert.equal(dockSettingsSection(state, 'linear', 'comparison'), 'comparison');
  assert.equal(dockSettingsSection(state, 'linear', 'comparison'), 'comparison');
  assert.deepEqual(read(state), ['rail', 'comparison']);
  assert.equal(dockSettingsSection(state, 'linear', 'basic'), 'input-genomes');
  assert.equal(dockSettingsSection(state, 'circular', 'preserved-settings'), 'input-genomes');
  assert.equal(dockSettingsSection(state, 'circular', 'about'), 'about');
  assert.deepEqual(read(state), ['rail', 'about']);
  // An entry of the other mode resolves to Input.
  assert.equal(dockSettingsSection(state, 'circular', 'advanced-comparison'), 'input-genomes');
  assert.equal(dockSettingsSection(state, 'linear', 'pairwise-comparisons'), 'input-genomes');
});

test('the stacked display docks nothing', () => {
  const state = paneState('stacked', null);
  assert.equal(dockSettingsSection(state, 'linear', 'labels'), null);
  assert.deepEqual(read(state), ['stacked', null]);
});

test('closing leaves the rail alone', () => {
  const state = paneState('rail', 'layout');
  closeDockedSettingsEntry(state);
  assert.deepEqual(read(state), ['rail', null]);
});

// Axis & Scale is one card per mode with one title: its entry follows the mode.
// Comparison and Pairwise Comparisons are different sections under one short name.
test('a mode change docks the entry of the same title, or Input, in place of an entry the mode lacks, and leaves the rail alone as it is', () => {
  for (const [mode, docked, expected] of [
    ['circular', 'comparison', 'input-genomes'],
    ['circular', 'advanced-comparison', 'input-genomes'],
    ['circular', 'linear-axis', 'circular-axis'],
    ['linear', 'pairwise-comparisons', 'input-genomes'],
    ['linear', 'circular-axis', 'linear-axis'],
    ['linear', 'comparison', 'comparison'],
    ['circular', 'labels', 'labels'],
    ['linear', 'about', 'about'],
    ['circular', null, null],
    ['linear', null, null]
  ]) {
    const state = paneState('rail', docked);
    assert.equal(reconcileSettingsPaneState(state, mode), expected, `${mode} ${docked}`);
    assert.deepEqual(read(state), ['rail', expected], `${mode} ${docked}`);
  }
  const stacked = paneState('stacked', null);
  reconcileSettingsPaneState(stacked, 'circular');
  assert.deepEqual(read(stacked), ['stacked', null]);
});

test('the popup holds one section until it returns, its entry is docked, or the mode changes', () => {
  const at = { x: 500, y: 200, maxHeight: 600 };
  const state = paneState('rail', 'layout');
  floatSettingsSection(state, 'dinucleotide', at);
  assert.deepEqual(state.settingsPopup.value, { section: 'dinucleotide', ...at });
  floatSettingsSection(state, 'titles', { ...at, x: 40 });
  assert.deepEqual(state.settingsPopup.value, { section: 'titles', x: 40, y: 200, maxHeight: 600 }, 'one popup: the section is swapped');
  // Docking another entry, closing the dock, and the stacked display leave it open.
  dockSettingsSection(state, 'circular', 'labels');
  closeDockedSettingsEntry(state);
  assert.equal(state.settingsPopup.value.section, 'titles');
  dockSettingsSection(state, 'circular', 'titles');
  assert.deepEqual([...read(state), state.settingsPopup.value], ['rail', 'titles', null], 'the dock shows the section');
  floatSettingsSection(state, 'circular-axis', at);
  showSettingsRail(state, 'circular', 'circular-axis');
  assert.equal(state.settingsPopup.value, null);
  floatSettingsSection(state, 'dinucleotide', at);
  showSettingsStacked(state);
  assert.equal(state.settingsPopup.value.section, 'dinucleotide');
  reconcileSettingsPaneState(state, 'linear');
  assert.equal(state.settingsPopup.value, null, 'a mode change closes the popup');
  floatSettingsSection(state, 'titles', at);
  returnSettingsSection(state);
  assert.equal(state.settingsPopup.value, null);
});

// A drag, the placement by the popup's rendered size, and a window that changed size move the popup.
test('moving the popup keeps its section, and moves nothing while no section is in it', () => {
  const state = paneState();
  moveSettingsPopup(state, { x: 10, y: 20, maxHeight: 0 });
  assert.equal(state.settingsPopup.value, null);
  floatSettingsSection(state, 'titles', { x: 500, y: 200, maxHeight: 600 });
  moveSettingsPopup(state, { x: 10, y: 20, maxHeight: 0 });
  assert.deepEqual(state.settingsPopup.value, { section: 'titles', x: 10, y: 20, maxHeight: 0 });
});
