// @ts-check
// When each setting reaches the Result (R1; PD-OI-037 rev 3 / OIC-024): the one
// owner of the timing notes and of how they show. index.html marks a row's scope with
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
 * the class text where parts of the change apply at different times (on a
 * live row it names the part that waits for Generate); `partialNote` is a live
 * row's note while the shown Result can draw only some of its request paths,
 * and `oneCirclePartialNote` replaces it while the change that waits enters or
 * leaves One circle (`canvasChangeIsOneCircle`, services/session-request.js).
 * `help` is the opening sentence of an editor tab's help, which says what the
 * tab's edits change at once (`settingApplicationHelp`). `refitPaths` name
 * the draft fields of an editor row that Python draws and sizes, by their
 * dotted path in the drawing's History intent (as in `SettingRequestPaths`):
 * on a figure drawn to a requested size, a History step that changes what
 * they draw asks the automatic rerender to fit the figure again (Figure size,
 * PD-OI-114 rev 2; `refitFieldsBeforeStep`), and the help says so in
 * `refitHelp`, which names the edits that do (review S4, D-61).
 * @typedef {{ readonly id: SettingApplicationId, readonly application: SettingApplicationClass, readonly note?: string, readonly partialNote?: string, readonly oneCirclePartialNote?: string, readonly help?: string, readonly requestPaths?: SettingRequestPaths, readonly refitPaths?: Readonly<Record<string, RefitField>>, readonly refitHelp?: string }} SettingApplicationRow
 */
/**
 * A drawing's draft field that Python draws and sizes: the Legend's listed and
 * deleted rows (the paint slice), the color rules, whose captions name the
 * rows Python draws from them (a Legend rename of a feature type row can write
 * them; the color rules slice), and the Canvas padding (the canvas slice).
 * @typedef {'legendEntries' | 'deletedLegendEntries' | 'rules' | 'canvasPadding'} RefitField
 */
/**
 * How a `<setting-note>` shows its note (PD-OI-037 rev 3): `clock`, a clock
 * icon named "Applies on Generate" whose tooltip and description are the text,
 * for a change that applies on Generate wholly or in part; `hint`, a line of
 * text (an export-only row); `none`, no mark: the change applies at once.
 * @typedef {'clock' | 'hint' | 'none'} SettingNoteDisplay
 */
/** @typedef {{ application: SettingApplicationClass, text: string, display: SettingNoteDisplay }} SettingApplicationNote */

/** @param {SettingApplicationClass} application @param {SettingApplicationId[]} ids @returns {SettingApplicationRow[]} */
const rows = (application, ids) => ids.map((id) => ({ id, application }));
const CANVAS_CHANGE_NOTE = msg('Live edit; the position applies on Generate after the Multi-Record Canvas change.');
const ONE_CIRCLE_CHANGE_NOTE = msg('Live edit; the position applies on Generate after the One circle change.');
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
    'feature-density', 'pairwise-comparisons', 'preserved-settings', 'layout', 'font', 'titles', 'features',
    'stroke', 'labels', 'axis-scale', 'dinucleotide', 'feature-placement']),
  // D-201, PD-OI-103 rev 2: the stack's structure edits apply on Generate; a
  // switch that hides or shows a GC, skew or Depth row is the track-visibility row's.
  { id: 'track-slots', application: /** @type {const} */ ('generate'),
    note: msg('Applies on Generate, except turning a GC content, GC skew or Depth row on or off, removing one, and Use custom stack, which apply at once. These still apply on Generate: turning on a row added since the last Generate, a row that was off when the Session was loaded, or a Depth row after a Depth file change; and Use custom stack after a Session load (except when turned off on a Linear Result drawn with it) or while the Region Annotation sets differ from those of the Result, with every track it would hide or show.') },
  // G1 (PD-OI-103): a committed value redraws the shown Result through the
  // automatic rerender (PD-OI-102). Legend position and Plot Title position
  // are the drawing's layout preferences (`ui.layoutPreferences`, per mode and
  // for Circular per single or multi-record canvas): a pending Multi-Record
  // Canvas change draws the other canvas's position on Generate (OV-1077), and
  // the position applies when the draft draws the shown canvas again (OV-1078).
  { id: 'legend', application: /** @type {const} */ ('live'), partialNote: CANVAS_CHANGE_NOTE, oneCirclePartialNote: ONE_CIRCLE_CHANGE_NOTE,
    requestPaths: {
    'ui.layoutPreferences.linear.legend': 'diagramOptions.output.legend',
    'ui.layoutPreferences.circular.*.legend': 'diagramOptions.output.legend',
    'config.form.multi_record_canvas': 'diagramOptions.output.legend',
    'config.adv.legend_box_size': 'diagramOptions.configOverrides.objects.legends.color_rect_size',
    'config.adv.legend_font_size': 'diagramOptions.configOverrides.objects.legends.font_size'
  } },
  { id: 'plot-title', application: /** @type {const} */ ('live'), partialNote: CANVAS_CHANGE_NOTE, oneCirclePartialNote: ONE_CIRCLE_CHANGE_NOTE,
    requestPaths: {
    'config.form.plot_title': 'diagramOptions.plotTitle',
    'ui.layoutPreferences.linear.plotTitlePosition': 'diagramOptions.output.plotTitlePosition',
    'ui.layoutPreferences.circular.*.plotTitlePosition': 'diagramOptions.output.plotTitlePosition',
    'config.form.multi_record_canvas': 'diagramOptions.output.plotTitlePosition',
    'config.adv.plot_title_font_size': 'diagramOptions.plotTitleFontSize'
  } },
  // Circular multi-record canvas layout; Multi-Record Canvas itself and
  // Reload Records change the record set and apply on Generate. The grid's
  // fields and One circle's are two rows: each applies live to a Result drawn
  // in its own layout (`liveSettingApplies`; PD-OI-103 G1, PD-OI-019 rev 2).
  { id: 'multi-record-layout', application: /** @type {const} */ ('live'), requestPaths: {
    'config.adv.multi_record_size_mode': 'layout.multiRecordSizeMode',
    'config.adv.multi_record_min_radius_ratio': 'layout.multiRecordMinRadiusRatio',
    'config.adv.multi_record_column_gap_ratio': 'layout.multiRecordColumnGapRatio',
    'config.adv.multi_record_row_gap_ratio': 'layout.multiRecordRowGapRatio',
    'config.adv.multi_record_positions': 'layout.multiRecordPositions'
  } },
  { id: 'one-circle-layout', application: /** @type {const} */ ('live'), requestPaths: {
    'config.adv.one_circle_gap_degrees': 'layout.oneCircleGapDegrees',
    'config.adv.one_circle_min_length': 'layout.oneCircleMinLength',
    'config.adv.one_circle_record_band': 'layout.oneCircleRecordBand'
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
  // PD-OI-117: the set a Sequence search adds is drawn at once, and its Undo
  // and Redo apply at once. The projection adds or removes only the sets a
  // search step added (`searchSetIdsOfStep`, by set id); the other changes of
  // the annotation sets keep the Region Annotations timing.
  { id: 'sequence-search', application: /** @type {const} */ ('live'), requestPaths: {
    'config.annotationSets': 'diagramOptions.annotations.sets'
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
    partialNote: msg('Live edit. Showing a track applies on Generate if the Result did not draw it when its Session was loaded; after a Depth file or series change, or a change of the records the Depth files map to; for a custom stack row added since the last Generate or a stack the last Generate could not build; and, with a Circular custom track stack or Circular track widths or radii, when Hide GC Content or Hide GC Skew is turned off.'),
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
    note: msg('Adding a Depth file, or a series while a Depth file is loaded, turns Show Depth on, and removing the last Depth file or its genome file turns Show Depth off; turning Show Depth off hides the drawn Depth track at once. A changed Legend title clears a Legend tab rename of the old title at once. The rest, and a Depth track from changed files, apply on Generate.') },
  ...rows('live', ['colors']),
  // Editor intent the Result shows at once (the Legend and Features tabs of
  // the editor, the feature popup's Appearance, Canvas padding): no mark. The
  // Legend names, deletions and order, and the Canvas padding, are drawn by
  // Python too (the request's config leaves, `legendEditConfigLeaves`): a
  // figure drawn to a requested size is fitted again after them (PD-OI-114
  // rev 2); the colors and strokes stay the Web's.
  /** @type {SettingApplicationRow} */ ({ id: 'legend-editor', application: 'live',
    help: msg('Live edit: legend text, order, stroke, and removal update the current Result.'),
    refitPaths: Object.freeze({
      'editorState.legend.entries': 'legendEntries', 'editorState.legend.deletedEntries': 'deletedLegendEntries', 'config.rules': 'rules'
    }),
    refitHelp: msg('With a Figure size set, a rename, removal, or reorder then redraws the figure to fit it.') }),
  /** @type {SettingApplicationRow} */ ({ id: 'canvas-padding', application: 'live',
    help: msg('Live edit: the padding updates the current Result.'),
    refitPaths: Object.freeze({ 'ui.canvasPadding': 'canvasPadding' }),
    refitHelp: msg('With a Figure size set, the figure is then redrawn to fit it.') }),
  { id: 'features-editor', application: /** @type {const} */ ('live'),
    help: msg('Live edit: feature color, label text, and visibility update the current Result; geometry may rerender automatically.') },
  { id: 'feature-appearance', application: /** @type {const} */ ('live') },
  // Show and Off redraw at once; matching reads the rules at Generate
  // (`should_include_feature_in_analysis`).
  { id: 'feature-visibility', application: /** @type {const} */ ('live'),
    note: msg('Showing and hiding apply at once; leaving a feature out of protein search applies on Generate.') },
  { id: 'font-embedding', application: /** @type {const} */ ('export'),
    note: msg('Applies to SVG and Interactive SVG exports, not to Generate.') }
].map((row) => Object.freeze(row.requestPaths ? { ...row, requestPaths: Object.freeze({ ...row.requestPaths }) } : row)));

/** @type {Readonly<Record<SettingApplicationClass, string>>} */
const TEXT = Object.freeze({ live: '', generate: msg('Applies on Generate'), apply: '', export: '' });
/** @type {Readonly<Record<SettingApplicationClass, SettingNoteDisplay>>} */
const DISPLAY = Object.freeze({ live: 'none', generate: 'clock', apply: 'none', export: 'hint' });
/**
 * @param {SettingApplicationClass} application
 * @param {string} text English, marked with `msg`.
 * @param {boolean} waits Whether a part of a live change waits for Generate.
 * @returns {SettingApplicationNote}
 */
const noteOf = (application, text, waits) => ({ application, text: tm(text), display: waits ? 'clock' : DISPLAY[application] });

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
 * The Vue injection key of whether the change that waits for Generate enters
 * or leaves One circle, which `SettingNote` passes to `settingApplicationNote`:
 * the app setup provides `committedOneCircleChange` (services/config.js).
 */
export const SETTING_NOTE_ONE_CIRCLE_CHANGE = Symbol('setting-note-one-circle-change');

/**
 * The note a `<setting-note>` shows; empty text shows none. A live setting
 * row applies on Generate while the shown Result can draw none of its request
 * paths (`applies`, `liveSettingApplies` of the committed request: a Circular
 * layout field without a committed grid), and shows its `partialNote` while
 * it can draw some (a Circular position during a Multi-Record Canvas change,
 * or its `oneCirclePartialNote` during a change that enters or leaves One
 * circle, `oneCircleChange`).
 * @param {SettingApplicationId} id
 * @param {(requestPath: string) => boolean} [applies]
 * @param {boolean} [oneCircleChange]
 * @returns {SettingApplicationNote}
 */
export const settingApplicationNote = (id, applies = () => true, oneCircleChange = false) => {
  const row = rowOf(id);
  const paths = [...new Set(Object.values(row.requestPaths || {}).flat())];
  const drawn = paths.filter((path) => applies(path)).length;
  if (drawn === paths.length) return noteOf(row.application, row.note || TEXT[row.application], row.application === 'live' && Boolean(row.note));
  const partialNote = (oneCircleChange && row.oneCirclePartialNote) || row.partialNote;
  return drawn > 0 && partialNote ? noteOf(row.application, partialNote, true) : noteOf('generate', TEXT.generate, false);
};

/**
 * The note of a class that follows state instead of a row's (Palette Instant
 * Preview, a queued palette).
 * @param {SettingApplicationClass} application
 * @returns {SettingApplicationNote}
 */
export const settingClassNote = (application) => noteOf(application, TEXT[application], false);

/** The accessible name of the clock: "Applies on Generate". */
export const settingClockLabel = () => tm(TEXT.generate);

/**
 * The opening sentences of an editor tab's help: what its edits change at
 * once, and, for a row with `refitPaths`, which of them fit a sized figure again (`refitHelp`).
 * @param {SettingApplicationId} id
 */
export const settingApplicationHelp = (id) => {
  const row = rowOf(id);
  return [row.help || '', row.refitHelp || ''].filter((text) => text !== '').map((text) => tm(text)).join(' ');
};

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

/**
 * One changed value of a History step (`HistoryIntentChange` in services/history.js).
 * @typedef {{ readonly path?: readonly (string | number)[], readonly before?: unknown, readonly after?: unknown,
 *   readonly beforeHas?: boolean, readonly afterHas?: boolean }} HistoryChangeLike
 */
/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/**
 * `root` with the value at `keys` replaced (removed when `has` is false); each object on the way is copied.
 * @param {unknown} root @param {readonly string[]} keys @param {boolean} has @param {unknown} value
 * @returns {unknown}
 */
const withValueAt = (root, keys, has, value) => {
  if (keys.length === 0) return has ? value : undefined;
  const copy = isRecord(root) ? { ...root } : {};
  const [key, ...rest] = keys;
  const next = withValueAt(copy[key], rest, has, value);
  if (next === undefined) delete copy[key];
  else copy[key] = next;
  return copy;
};
/**
 * The refit fields (`refitPaths`) on the other side of a History step of the
 * shown drawing: `current` holds them as the step left them (after its commit
 * or Redo; before it, for its Undo), and each one the step changed takes the
 * value it had on the other side. Null when the step changed none of them.
 * @param {readonly HistoryChangeLike[] | null | undefined} changes
 * @param {Readonly<Record<RefitField, unknown>>} current
 * @param {boolean} undone The step is being undone.
 * @returns {Record<RefitField, unknown> | null}
 */
export const refitFieldsBeforeStep = (changes, current, undone) => {
  const side = undone ? 'after' : 'before';
  const fields = SETTING_APPLICATION.flatMap((row) => Object.entries(row.refitPaths || {}))
    .map(([path, field]) => /** @type {const} */ ([path.split('.'), field]));
  /** @type {Record<RefitField, unknown>} */
  const other = { ...current };
  let changed = false;
  for (const change of Array.isArray(changes) ? changes : []) {
    const path = Array.isArray(change?.path) ? change.path.map(String) : [];
    const has = Boolean(change[`${side}Has`]);
    for (const [keys, field] of fields) {
      const shared = Math.min(keys.length, path.length);
      if (path.length === 0 || !keys.slice(0, shared).every((key, index) => key === path[index])) continue;
      changed = true;
      // A change at the field or above it holds the whole value; one below it, a part of the value.
      other[field] = path.length <= keys.length
        ? keys.slice(path.length).reduce((/** @type {unknown} */ node, key) => (isRecord(node) ? node[key] : undefined), has ? change[side] : undefined)
        : withValueAt(other[field], path.slice(keys.length), has, change[side]);
    }
  }
  return changed ? other : null;
};

/** Every request path a live setting row applies, as Reset Settings and an artifact restore request them. */
export const allLiveSettingRequestPaths = () => [...new Set(SETTING_APPLICATION.flatMap((row) => Object.values(row.requestPaths || {}).flat()))];
