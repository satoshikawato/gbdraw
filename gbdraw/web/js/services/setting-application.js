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
 * A live row's draft fields and the request path each one applies to: a
 * dotted path of the drawing's History intent (`config.adv.legend_font_size`;
 * `*` matches one key) -> a request path that `services/session-request.js`
 * projects. Only the rows with these paths are the live setting rows the
 * automatic rerender applies from the draft (R1 (c), PD-OI-102/103); the
 * other live rows are editor intent.
 * @typedef {Readonly<Record<string, string>>} SettingRequestPaths
 */
/**
 * One setting, or a group of settings with one application. `note` replaces
 * the class text where parts of the change apply at different times.
 * @typedef {{ readonly id: SettingApplicationId, readonly application: SettingApplicationClass, readonly note?: string, readonly requestPaths?: SettingRequestPaths }} SettingApplicationRow
 */
/** @typedef {{ application: SettingApplicationClass, text: string, qualified: boolean }} SettingApplicationNote */

/** @param {SettingApplicationClass} application @param {SettingApplicationId[]} ids @returns {SettingApplicationRow[]} */
const rows = (application, ids) => ids.map((id) => ({ id, application }));

/** @type {readonly Readonly<SettingApplicationRow>[]} */
export const SETTING_APPLICATION = Object.freeze([
  // OV-1042: `output-prefix` names the Results at Generate, and exports take
  // their file names from the Result (export.js `getDownloadName`).
  ...rows('generate', ['output-prefix', 'inputs', 'comparison', 'species-strain', 'advanced-comparison', 'region-annotations',
    'pairwise-comparisons', 'preserved-settings', 'layout', 'track-slots', 'font', 'titles', 'features',
    'stroke', 'labels', 'axis-scale', 'dinucleotide', 'feature-placement']),
  // G1 (PD-OI-103): a committed value redraws the shown Result through the
  // automatic rerender (PD-OI-102). Legend position and Plot Title position
  // are the drawing's layout preferences (`ui.layoutPreferences`, per mode and
  // for Circular per single or multi-record canvas).
  { id: 'legend', application: /** @type {const} */ ('live'), requestPaths: {
    'ui.layoutPreferences.linear.legend': 'diagramOptions.output.legend',
    'ui.layoutPreferences.circular.*.legend': 'diagramOptions.output.legend',
    'config.adv.legend_box_size': 'diagramOptions.configOverrides.objects.legends.color_rect_size',
    'config.adv.legend_font_size': 'diagramOptions.configOverrides.objects.legends.font_size'
  } },
  { id: 'plot-title', application: /** @type {const} */ ('live'), requestPaths: {
    'config.form.plot_title': 'diagramOptions.plotTitle',
    'ui.layoutPreferences.linear.plotTitlePosition': 'diagramOptions.output.plotTitlePosition',
    'ui.layoutPreferences.circular.*.plotTitlePosition': 'diagramOptions.output.plotTitlePosition',
    'config.adv.plot_title_font_size': 'diagramOptions.plotTitleFontSize'
  } },
  // Circular multi-record canvas layout; Multi-Record Canvas itself and
  // Reload Records change the record set and apply on Generate.
  { id: 'multi-record-layout', application: /** @type {const} */ ('live'), requestPaths: {
    'config.adv.multi_record_size_mode': 'layout.multiRecordSizeMode',
    'config.adv.multi_record_min_radius_ratio': 'layout.multiRecordMinRadiusRatio',
    'config.adv.multi_record_column_gap_ratio': 'layout.multiRecordColumnGapRatio',
    'config.adv.multi_record_row_gap_ratio': 'layout.multiRecordRowGapRatio',
    'config.adv.multi_record_positions': 'layout.multiRecordPositions'
  } },
  // OV-1036: `applyTrackVisibility` (app/svg-styles.js) hides or shows a drawn
  // group at once; a track the Result lacks and the layout wait for Generate.
  { id: 'track-visibility', application: /** @type {const} */ ('generate'),
    note: 'Hiding or showing a drawn track applies at once; a new track and the layout apply on Generate' },
  // Show Depth follows the Depth files: adding one, or a series while one is
  // loaded (`addCircularDepthTrack`, `addLinearDepthTrack`), turns it on;
  // removing the last one, also with its Linear File, turns it off (OV-82,
  // OV-1041); `applyTrackVisibility` shows or hides the drawn Depth track at
  // once. A Legend title change retires the old caption, which clears a Legend
  // tab rename of it at once (`buildLegendStyleRetirement`, OV-1043).
  { id: 'depth-sources', application: /** @type {const} */ ('generate'),
    note: 'Adding a Depth file, or a series while one is loaded, turns Show Depth on, and removing the last Depth file or the File that holds it turns it off; the drawn Depth track shows or hides at once. A changed Legend title clears a Legend tab rename of the old title at once. The rest applies on Generate' },
  ...rows('live', ['colors']),
  // Editor intent the Result shows at once (the Legend and Features tabs of
  // the editor, the feature popup's Appearance).
  { id: 'legend-editor', application: /** @type {const} */ ('live'),
    note: 'Live edit: legend text, order, stroke, and removal update the current Result.' },
  { id: 'features-editor', application: /** @type {const} */ ('live'),
    note: 'Live edit: feature color, label text, and visibility update the current Result; geometry may rerender automatically.' },
  { id: 'feature-appearance', application: /** @type {const} */ ('live'), note: 'Updates the current Result' },
  // Show and Off redraw at once; matching reads the rules at Generate
  // (`should_include_feature_in_analysis`).
  { id: 'feature-visibility', application: /** @type {const} */ ('live'),
    note: 'Showing and hiding apply at once; leaving a feature out of protein search applies on Generate' },
  { id: 'font-embedding', application: /** @type {const} */ ('export'),
    note: 'Applies to SVG and Interactive SVG exports, not to Generate.' }
].map((row) => Object.freeze(row.requestPaths ? { ...row, requestPaths: Object.freeze({ ...row.requestPaths }) } : row)));

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

/**
 * The request paths of the live setting rows that a History step changed:
 * its change paths, `modes.<mode>.`-prefixed as a step records them or
 * already relative to the drawing as a restore passes them. A change path and
 * a row's draft path match when one is a prefix of the other, so a step that
 * replaced a whole object names every live field under it; a path whose
 * projected value equals the committed one applies nothing.
 * @param {readonly { path?: readonly (string | number)[] }[] | null | undefined} changes
 * @param {'circular' | 'linear'} mode The drawing whose fields count.
 * @returns {string[]}
 */
export const liveSettingRequestPaths = (changes, mode) => {
  /** @type {Set<string>} */
  const paths = new Set();
  for (const change of Array.isArray(changes) ? changes : []) {
    const raw = Array.isArray(change?.path) ? change.path.map(String) : [];
    if (raw[0] === 'modes' && raw.length > 1 && raw[1] !== mode) continue;
    const path = raw[0] === 'modes' ? raw.slice(2) : raw;
    for (const row of SETTING_APPLICATION) {
      Object.entries(row.requestPaths || {}).forEach(([draftPath, requestPath]) => {
        const keys = draftPath.split('.');
        const overlap = Math.min(keys.length, path.length);
        if (keys.slice(0, overlap).every((key, index) => key === '*' || key === path[index])) paths.add(requestPath);
      });
    }
  }
  return [...paths];
};

/** Every request path a live setting row applies, as Reset Settings and an artifact restore request them. */
export const allLiveSettingRequestPaths = () => [...new Set(SETTING_APPLICATION.flatMap((row) => Object.values(row.requestPaths || {})))];
