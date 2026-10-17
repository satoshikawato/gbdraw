// @ts-check
// When each setting reaches the Result (R1; PD-OI-037 rev 2 / OIC-024): the one
// owner of the timing notes. index.html marks a row's scope with
// `data-setting="<id>"`, a control belongs to the row of its nearest scope, and
// `<setting-note>` (components.js) renders the row's note. A state-free leaf, so
// the page reads it without importing the request owner.

/**
 * `live`: an editor edit or an automatic rerender; `generate`: the next
 * Generate; `apply`: its own Apply; `export`: exported files only.
 * @typedef {'live' | 'generate' | 'apply' | 'export'} SettingApplicationClass
 */
/** @typedef {string} SettingApplicationId A row key such as `legend`; never a caption. */
/**
 * One setting, or a group of settings with one application. `note` replaces
 * the class text where parts of the change apply at different times.
 * @typedef {{ readonly id: SettingApplicationId, readonly application: SettingApplicationClass, readonly note?: string }} SettingApplicationRow
 */
/** @typedef {{ application: SettingApplicationClass, text: string, qualified: boolean }} SettingApplicationNote */

/** @param {SettingApplicationClass} application @param {SettingApplicationId[]} ids @returns {SettingApplicationRow[]} */
const rows = (application, ids) => ids.map((id) => ({ id, application }));

/** @type {readonly Readonly<SettingApplicationRow>[]} */
export const SETTING_APPLICATION = Object.freeze([
  // OV-1042: `output-prefix` names the Results at Generate, and exports take
  // their file names from the Result (export.js `getDownloadName`).
  ...rows('generate', ['output-prefix', 'inputs', 'comparison', 'species-strain', 'advanced-comparison', 'region-annotations',
    'pairwise-comparisons', 'preserved-settings', 'layout', 'track-slots', 'font', 'titles', 'legend', 'features',
    'stroke', 'labels', 'axis-scale', 'dinucleotide', 'feature-placement']),
  // OV-1036: `applyTrackVisibility` (app/svg-styles.js) hides or shows a drawn
  // group at once; a track the Result lacks and the layout wait for Generate.
  { id: 'track-visibility', application: /** @type {const} */ ('generate'),
    note: 'Hiding or showing a drawn track applies at once; a new track and the layout apply on Generate' },
  // Show Depth follows the Depth files: adding one turns it on, removing the
  // last one turns it off (OV-82), and `applyTrackVisibility` shows or hides
  // the drawn Depth track at once.
  { id: 'depth-sources', application: /** @type {const} */ ('generate'),
    note: 'Adding a Depth file shows and removing the last one hides the drawn Depth track at once; the rest applies on Generate' },
  ...rows('live', ['colors']),
  // Show and Off redraw at once; matching reads the rules at Generate
  // (`should_include_feature_in_analysis`).
  { id: 'feature-visibility', application: /** @type {const} */ ('live'),
    note: 'Showing and hiding apply at once; leaving a feature out of protein search applies on Generate' },
  { id: 'font-embedding', application: /** @type {const} */ ('export'),
    note: 'Applies to SVG and Interactive SVG exports, not to Generate.' }
].map((row) => Object.freeze(row)));

/** @type {Readonly<Record<SettingApplicationClass, string>>} */
const TEXT = Object.freeze({ live: 'Live edit', generate: 'Applies on Generate', apply: '', export: '' });

/** @param {SettingApplicationId} id */
const rowOf = (id) => {
  const row = SETTING_APPLICATION.find((candidate) => candidate.id === id);
  if (!row) throw new Error(`Unknown setting application row: ${id}`);
  return row;
};

/** The note a `<setting-note>` shows; empty text shows none. @param {SettingApplicationId} id @returns {SettingApplicationNote} */
export const settingApplicationNote = (id) => {
  const row = rowOf(id);
  return { application: row.application, text: row.note || TEXT[row.application], qualified: Boolean(row.note) };
};

/** The text of a note whose class follows state (Palette Instant Preview). @param {SettingApplicationClass} application */
export const settingApplicationText = (application) => TEXT[application];

/** The row's timing inside a help text or a dialog description. @param {SettingApplicationId} id */
export const settingApplicationSentence = (id) => (rowOf(id).application === 'generate' ? 'Changes apply on Generate.' : '');
