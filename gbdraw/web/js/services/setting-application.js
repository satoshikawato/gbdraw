// @ts-check
// When each setting reaches the Result (R1; PD-OI-037 rev 2 / OIC-024): the one
// owner of the timing notes. index.html marks a row's scope with
// `data-setting="<id>"`, a control belongs to the row of its nearest scope, and
// `<setting-note>` (components.js) renders the row's note. A state-free leaf, so
// the page reads it without importing the request owner. The notes are English
// marked with `msg`; the functions that return a shown text translate them.
import { msg, t, tm } from '../utils/i18n.js';

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
 * projects. A field that decides whether a request path applies to the shown
 * Result (`liveSettingApplies`) names it too, so a change of that field
 * requests it again (OV-1078). A draft path may request several request paths
 * that one projection writes together (the GT set). Only the rows with these
 * paths are the live setting rows the automatic rerender applies from the
 * draft (R1 (c), PD-OI-102/103/104); the other live rows are editor intent.
 * @typedef {Readonly<Record<string, string | readonly string[]>>} SettingRequestPaths
 */
/**
 * One setting, or a group of settings with one application. `note` replaces
 * the class text where parts of the change apply at different times;
 * `partialNote` is a live row's note while the shown Result can draw only some
 * of its request paths.
 * @typedef {{ readonly id: SettingApplicationId, readonly application: SettingApplicationClass, readonly note?: string, readonly partialNote?: string, readonly requestPaths?: SettingRequestPaths }} SettingApplicationRow
 */
/** @typedef {{ application: SettingApplicationClass, text: string, qualified: boolean }} SettingApplicationNote */

/** @param {SettingApplicationClass} application @param {SettingApplicationId[]} ids @returns {SettingApplicationRow[]} */
const rows = (application, ids) => ids.map((id) => ({ id, application }));
const CANVAS_CHANGE_NOTE = msg('Live edit; the position applies on Generate after the Multi-Record Canvas change');
// The GT set (PD-OI-104 rev 5): one projection (`projectTrackPlan`,
// services/session-request.js) writes these request paths together, so every
// draft input of a drawn GC, skew or Depth track requests all of them.
const TRACK_VISIBILITY_REQUEST_PATHS = Object.freeze([
  'diagramOptions.configOverrides.canvas.show_gc',
  'diagramOptions.configOverrides.canvas.show_skew',
  'diagramOptions.depthTracks',
  'diagramOptions.tracks.circularTrackSlots',
  'diagramOptions.tracks.linearTrackSlots'
]);
/** @param {readonly string[]} draftPaths @returns {SettingRequestPaths} */
const trackVisibilityPaths = (draftPaths) => Object.fromEntries(draftPaths.map((path) => [path, TRACK_VISIBILITY_REQUEST_PATHS]));

/** @type {readonly Readonly<SettingApplicationRow>[]} */
export const SETTING_APPLICATION = Object.freeze([
  // OV-1042: `output-prefix` names the Results at Generate, and exports take
  // their file names from the Result (export.js `getDownloadName`).
  ...rows('generate', ['output-prefix', 'inputs', 'comparison', 'species-strain', 'advanced-comparison', 'region-annotations',
    'pairwise-comparisons', 'preserved-settings', 'layout', 'font', 'titles', 'features',
    'stroke', 'labels', 'axis-scale', 'dinucleotide', 'feature-placement']),
  // D-201, PD-OI-103 rev 2: the stack's structure edits apply on Generate; a
  // switch that hides or shows a GC, skew or Depth row is the track-visibility row's.
  { id: 'track-slots', application: /** @type {const} */ ('generate'),
    note: msg('Applies on Generate, except turning a GC content, GC skew or Depth row on or off, removing one, and Use custom stack, which apply at once; turning on a row added since the last Generate, a row that was off when the Session was loaded, or a Depth row after a Depth file change apply on Generate, and so does Use custom stack after a Session load (except turned off on a Linear Result drawn with it) or while the Region Annotation sets differ from those of the Result, with every track it would hide or show') },
  // G1 (PD-OI-103): a committed value redraws the shown Result through the
  // automatic rerender (PD-OI-102). Legend position and Plot Title position
  // are the drawing's layout preferences (`ui.layoutPreferences`, per mode and
  // for Circular per single or multi-record canvas): a pending Multi-Record
  // Canvas change draws the other canvas's position on Generate (OV-1077), and
  // the position applies when the draft draws the shown canvas again (OV-1078).
  { id: 'legend', application: /** @type {const} */ ('live'), partialNote: CANVAS_CHANGE_NOTE, requestPaths: {
    'ui.layoutPreferences.linear.legend': 'diagramOptions.output.legend',
    'ui.layoutPreferences.circular.*.legend': 'diagramOptions.output.legend',
    'config.form.multi_record_canvas': 'diagramOptions.output.legend',
    'config.adv.legend_box_size': 'diagramOptions.configOverrides.objects.legends.color_rect_size',
    'config.adv.legend_font_size': 'diagramOptions.configOverrides.objects.legends.font_size'
  } },
  { id: 'plot-title', application: /** @type {const} */ ('live'), partialNote: CANVAS_CHANGE_NOTE, requestPaths: {
    'config.form.plot_title': 'diagramOptions.plotTitle',
    'ui.layoutPreferences.linear.plotTitlePosition': 'diagramOptions.output.plotTitlePosition',
    'ui.layoutPreferences.circular.*.plotTitlePosition': 'diagramOptions.output.plotTitlePosition',
    'config.form.multi_record_canvas': 'diagramOptions.output.plotTitlePosition',
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
  // Linear Lock Definition Column (PD-OI-103 G1).
  { id: 'definition-column', application: /** @type {const} */ ('live'), requestPaths: {
    'config.form.keep_definition_left_aligned': 'diagramOptions.configOverrides.canvas.linear.keep_definition_left_aligned'
  } },
  // G2 (PD-OI-103 rev 2, PD-OI-059 rev 3, D-150, D-151): each row a control
  // group the evidence rule admitted; the other controls of their sections
  // stay in their `generate` rows. A row with a path of one mode shows its
  // note by the shown mode's paths (`liveSettingNoteApplies`). An annotation
  // item's colors and a Linear row's Height show a note per item, which
  // states whether that item's value is live (`liveSettingItemApplies`).
  { id: 'definition-text', application: /** @type {const} */ ('live'), requestPaths: {
    'config.adv.def_font_size': 'diagramOptions.configOverrides.objects.definition.font_size',
    'config.adv.linear_definition_line_styles': 'diagramOptions.configOverrides.objects.definition.linear.line_styles'
  } },
  { id: 'species', application: /** @type {const} */ ('live'), requestPaths: {
    'config.form.species': 'diagramOptions.species'
  } },
  { id: 'dinucleotide-window', application: /** @type {const} */ ('live'), requestPaths: {
    'config.adv.window_size': 'diagramOptions.window'
  } },
  { id: 'annotation-style', application: /** @type {const} */ ('live'), requestPaths: {
    'config.annotationSets.*.annotations.*.style': 'diagramOptions.annotations.sets'
  } },
  { id: 'comparison-height', application: /** @type {const} */ ('live'), requestPaths: {
    'config.adv.comparison_height': 'diagramOptions.configOverrides.canvas.linear.comparison_height'
  } },
  { id: 'track-geometry', application: /** @type {const} */ ('live'), requestPaths: {
    'config.adv.center_reserved_radius': 'diagramOptions.tracks.centerReservedRadius',
    'config.adv.linear_track_slots.*.height': 'diagramOptions.tracks.linearTrackSlots'
  } },
  // The label fields only; the Linear card order decides whether they apply
  // (`liveSettingApplies`), so a change of the cards requests them again.
  { id: 'record-names', application: /** @type {const} */ ('live'), requestPaths: {
    'files.linearSeqs': 'records',
    'files.linearSeqs.*.definition': 'records',
    'config.form.circular_record_label': 'records',
    'config.form.circular_record_subtitle': 'records'
  } },
  // PD-OI-104 rev 5, D-199..D-202, D-207: the GC content, GC skew and Depth
  // toggles, the slot switches that hide or show one of those tracks, and Use
  // custom stack redraw the Result through the automatic rerender. A show
  // returns the parts the last Generate built (the kept track parts); it waits
  // for Generate where they cannot give it (`partialNote`).
  { id: 'track-visibility', application: /** @type {const} */ ('live'),
    partialNote: msg('Live edit. Showing a track applies on Generate if the Result did not draw it when its Session was loaded; after a Depth file or series change, or a change of the records the Depth files map to; for a custom stack row added since the last Generate or a stack the last Generate could not build; and, with a Circular custom track stack or Circular track widths or radii, for Hide GC Content and Hide GC Skew turned off'),
    requestPaths: trackVisibilityPaths([
      'config.form.show_gc', 'config.form.show_skew', 'config.form.suppress_gc', 'config.form.suppress_skew',
      'config.form.show_depth', 'files.c_depth', 'files.linearSeqs.*.depth', 'config.adv.depth_tracks',
      'config.adv.circular_track_slots', 'config.adv.circular_track_slots_enabled',
      'config.adv.linear_track_slots', 'config.adv.linear_track_slots_enabled'
    ]) },
  // Show Depth follows the Depth files: adding one, or a series while one is
  // loaded (`addCircularDepthTrack`, `addLinearDepthTrack`), turns it on;
  // removing the last one, also with its Linear File, turns it off (OV-82,
  // OV-1041), and the track-visibility row redraws the Result. A Legend title
  // change retires the old caption, which clears a Legend tab rename of it at
  // once (`buildLegendStyleRetirement`, OV-1043).
  { id: 'depth-sources', application: /** @type {const} */ ('generate'),
    note: msg('Adding a Depth file, or a series while one is loaded, turns Show Depth on, and removing the last Depth file or the File that holds it turns it off; turning it off hides the drawn Depth track at once. A changed Legend title clears a Legend tab rename of the old title at once. The rest, and a Depth track from changed files, applies on Generate') },
  ...rows('live', ['colors']),
  // Editor intent the Result shows at once (the Legend and Features tabs of
  // the editor, the feature popup's Appearance).
  { id: 'legend-editor', application: /** @type {const} */ ('live'),
    note: msg('Live edit: legend text, order, stroke, and removal update the current Result.') },
  { id: 'features-editor', application: /** @type {const} */ ('live'),
    note: msg('Live edit: feature color, label text, and visibility update the current Result; geometry may rerender automatically.') },
  { id: 'feature-appearance', application: /** @type {const} */ ('live'), note: msg('Updates the current Result') },
  // Show and Off redraw at once; matching reads the rules at Generate
  // (`should_include_feature_in_analysis`).
  { id: 'feature-visibility', application: /** @type {const} */ ('live'),
    note: msg('Showing and hiding apply at once; leaving a feature out of protein search applies on Generate') },
  { id: 'font-embedding', application: /** @type {const} */ ('export'),
    note: msg('Applies to SVG and Interactive SVG exports, not to Generate.') }
].map((row) => Object.freeze(row.requestPaths ? { ...row, requestPaths: Object.freeze({ ...row.requestPaths }) } : row)));

/** @type {Readonly<Record<SettingApplicationClass, string>>} */
const TEXT = Object.freeze({ live: msg('Live edit'), generate: msg('Applies on Generate'), apply: '', export: '' });

/** @param {SettingApplicationId} id */
const rowOf = (id) => {
  const row = SETTING_APPLICATION.find((candidate) => candidate.id === id);
  if (!row) throw new Error(`Unknown setting application row: ${id}`);
  return row;
};

/**
 * The Vue injection key of the maker of the `applies` predicate `SettingNote`
 * (components.js) passes to `settingApplicationNote`, one per render: the app
 * setup provides `committedSettingNotePredicate` (services/config.js).
 */
export const SETTING_NOTE_APPLIES = Symbol('setting-note-applies');

/**
 * The note a `<setting-note>` shows; empty text shows none. A live setting
 * row applies on Generate while the shown Result can draw none of its request
 * paths (`applies`, `liveSettingApplies` of the committed request: a Circular
 * layout field without a committed grid), and shows its `partialNote` while
 * it can draw some (a Circular position during a Multi-Record Canvas change).
 * @param {SettingApplicationId} id
 * @param {(requestPath: string) => boolean} [applies]
 * @returns {SettingApplicationNote}
 */
export const settingApplicationNote = (id, applies = () => true) => {
  const row = rowOf(id);
  const paths = [...new Set(Object.values(row.requestPaths || {}).flat())];
  const drawn = paths.filter((path) => applies(path)).length;
  if (drawn === paths.length) {
    return { application: row.application, text: tm(row.note || TEXT[row.application]), qualified: Boolean(row.note) };
  }
  return drawn > 0 && row.partialNote
    ? { application: row.application, text: tm(row.partialNote), qualified: true }
    : { application: 'generate', text: tm(TEXT.generate), qualified: false };
};

/** The text of a note whose class follows state (Palette Instant Preview). @param {SettingApplicationClass} application */
export const settingApplicationText = (application) => tm(TEXT[application]);

/** The row's timing inside a help text or a dialog description. @param {SettingApplicationId} id */
export const settingApplicationSentence = (id) => (rowOf(id).application === 'generate' ? t('Changes apply on Generate.') : '');

/**
 * The request paths of the live setting rows that a History step of the shown
 * drawing changed: its change paths, relative to the drawing (each drawing has
 * its own History, PD-OI-090). A change path and a row's draft path match when
 * one is a prefix of the other, so a step that replaced a whole object names
 * every live field under it; a path whose projected value equals the committed
 * one applies nothing.
 * @param {readonly { path?: readonly (string | number)[] }[] | null | undefined} changes
 * @returns {string[]}
 */
export const liveSettingRequestPaths = (changes) => {
  /** @type {Set<string>} */
  const paths = new Set();
  for (const change of Array.isArray(changes) ? changes : []) {
    const path = Array.isArray(change?.path) ? change.path.map(String) : [];
    for (const row of SETTING_APPLICATION) {
      Object.entries(row.requestPaths || {}).forEach(([draftPath, requestPaths]) => {
        const keys = draftPath.split('.');
        const overlap = Math.min(keys.length, path.length);
        if (keys.slice(0, overlap).every((key, index) => key === '*' || key === path[index])) {
          [requestPaths].flat().forEach((requestPath) => paths.add(requestPath));
        }
      });
    }
  }
  return [...paths];
};

/**
 * The kind of track whose Circular hide control an Undo took back.
 * @typedef {'content' | 'skew'} UndoneTrackHide
 */
/** @type {readonly { readonly control: 'suppress_gc' | 'suppress_skew', readonly track: UndoneTrackHide }[]} */
const HIDE_CONTROLS = Object.freeze([{ control: 'suppress_gc', track: 'content' }, { control: 'suppress_skew', track: 'skew' }]);
/**
 * The hides a History step changed. Its Undo takes them back: their rows stop
 * waiting for Generate and return from the kept parts of the last Generate, as
 * the track projection in `services/session-request.js` gives them. A row whose
 * hide is still on stays hidden: the projection reads the hide first.
 * @param {readonly { readonly path: readonly unknown[] }[]} changes The step's changes.
 * @returns {Set<UndoneTrackHide>}
 */
export const changedTrackHides = (changes) => {
  const paths = new Set(changes.map(({ path }) => path.join('.')));
  return new Set(HIDE_CONTROLS.filter(({ control }) => paths.has(`config.form.${control}`)).map(({ track }) => track));
};
/**
 * The hides that are on: an Undo's name lives only while its hide is off.
 * @param {Readonly<Record<'suppress_gc' | 'suppress_skew', boolean>>} form
 * @returns {Set<UndoneTrackHide>}
 */
export const trackHidesOn = (form) => new Set(HIDE_CONTROLS.filter(({ control }) => form[control]).map(({ track }) => track));

/** Every request path a live setting row applies, as Reset Settings and an artifact restore request them. */
export const allLiveSettingRequestPaths = () => [...new Set(SETTING_APPLICATION.flatMap((row) => Object.values(row.requestPaths || {}).flat()))];
