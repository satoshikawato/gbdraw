// @ts-check
// When each setting reaches the Result (R1; PD-OI-037 rev 2 / OIC-024): the one
// owner of the timing notes. `components.js` (`SettingNote`) and the template
// helpers in `app.js` render these rows; a later live apply reads the same rows.
// A state-free leaf: the request owner (`services/session-request.js`) may
// import it, and the page can import it without importing the request owner.

/**
 * When a setting's change reaches the Result (R1): `live` through an editor
 * edit or an automatic rerender, `generate` on the next Generate, `apply` by its
 * own Apply, `export` only in exported files.
 * @typedef {'live' | 'generate' | 'apply' | 'export'} SettingApplicationClass
 */
/**
 * One setting, or one group of settings that share an application, keyed by a
 * stable id. index.html marks the row's scope with `data-setting="<id>"`; a
 * control belongs to the row of its nearest scope, and `<setting-note>` renders
 * the row's note. `note` replaces the "Applies on Generate" text where part of
 * the change already shows before Generate.
 * @typedef {{ readonly id: SettingApplicationId, readonly application: SettingApplicationClass, readonly note?: string }} SettingApplicationRow
 */
/** @typedef {string} SettingApplicationId A row key such as `legend`; never a caption. */
/**
 * @typedef {{ application: SettingApplicationClass, text: string, qualified: boolean }} SettingApplicationNote
 *   `qualified`: the text is the row's own `note`.
 */

/** @type {readonly Readonly<SettingApplicationRow>[]} */
export const SETTING_APPLICATION = Object.freeze(/** @type {SettingApplicationRow[]} */ ([
  { id: 'inputs', application: 'generate' },
  { id: 'comparison', application: 'generate' },
  { id: 'output-prefix', application: 'export' },
  { id: 'species-strain', application: 'generate' },
  { id: 'advanced-comparison', application: 'generate' },
  { id: 'region-annotations', application: 'generate' },
  { id: 'pairwise-comparisons', application: 'generate' },
  { id: 'preserved-settings', application: 'generate' },
  { id: 'layout', application: 'generate' },
  // OV-1036: `applyTrackVisibility` (app/svg-styles.js) hides or shows a drawn
  // group at once; a track the Result lacks and the layout wait for Generate.
  {
    id: 'track-visibility',
    application: 'generate',
    note: 'Hiding or showing a drawn track applies at once; a new track and the layout apply on Generate'
  },
  { id: 'track-slots', application: 'generate' },
  { id: 'font', application: 'generate' },
  { id: 'titles', application: 'generate' },
  { id: 'legend', application: 'generate' },
  { id: 'colors', application: 'live' },
  { id: 'features', application: 'generate' },
  { id: 'feature-visibility', application: 'live' },
  { id: 'stroke', application: 'generate' },
  { id: 'labels', application: 'generate' },
  { id: 'axis-scale', application: 'generate' },
  { id: 'dinucleotide', application: 'generate' },
  { id: 'feature-placement', application: 'generate' }
]).map((row) => Object.freeze(row)));

/** @type {Readonly<Record<SettingApplicationClass, string>>} */
const SETTING_APPLICATION_TEXT = Object.freeze({
  live: 'Live edit',
  generate: 'Applies on Generate',
  apply: '',
  export: ''
});

/** @param {SettingApplicationId} id */
const settingApplicationRow = (id) => {
  const row = SETTING_APPLICATION.find((candidate) => candidate.id === id);
  if (!row) throw new Error(`Unknown setting application row: ${id}`);
  return row;
};

/**
 * The note a `<setting-note>` shows for a row; empty text shows none.
 * @param {SettingApplicationId} id
 * @returns {SettingApplicationNote}
 */
export const settingApplicationNote = (id) => {
  const row = settingApplicationRow(id);
  return { application: row.application, text: row.note || SETTING_APPLICATION_TEXT[row.application], qualified: Boolean(row.note) };
};

/**
 * The note text of an application class, for a note whose class follows state
 * (Palette Instant Preview, a pending palette draft).
 * @param {SettingApplicationClass} application
 */
export const settingApplicationText = (application) => SETTING_APPLICATION_TEXT[application];

/**
 * The row's timing as a sentence inside a help text or a dialog description.
 * @param {SettingApplicationId} id
 */
export const settingApplicationSentence = (id) => (
  settingApplicationRow(id).application === 'generate' ? 'Changes apply on Generate.' : ''
);
