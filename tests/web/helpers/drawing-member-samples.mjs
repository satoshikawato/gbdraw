// One non-default sample per member of a drawing (the members of
// `DRAWING_SLICES`, helpers/drawing-state.mjs), in the form the member's slice
// `capture` returns and `install` takes. A sample is a value the product
// itself writes and keeps: a fixed point of Session Load, so a test that saves
// a drawing holding it and loads the Session finds it unchanged. Where Load
// normalizes a value, the sample is already in Load's form, and its comment
// names the normalization.
//
// A sample is a function of (mode, drawing, context):
// - `mode` is the drawing's mode ('circular' or 'linear');
// - `drawing` is the drawing the sample goes into, read for the member's
//   current value (a sample of an object member changes every leaf of it,
//   `variedLeaves`) and for the records, rows, and ring tables it binds;
// - `context` holds `registerFile(file)`, which registers a File in the
//   project's file store (`state.inputFiles`) and returns its binding.
// A member whose sample needs more than a new drawing holds lists it in
// `needs`: 'records' (the drawing binds a GenBank file of two records, whose
// cards or rows the sample names), 'rings' (a Circular drawing binds a ring
// table, whose ring row the sample changes), 'rows' (a Linear drawing has two
// record rows, each with its own file, which the sample's comparison edge
// pairs), 'groups' (the drawing holds a Result with
// Similarity groups, whose ids the sample names), or 'comparison' (the
// drawing holds a Result with a comparison read from a CLI table). A context
// inside the drawing's own slices is `prepare(drawing)` (the queued palette a
// queued color belongs to). The test that installs the sample builds that
// context first, as the user would.
//
// A member that only an action writes (`byContext`: the records a file read
// discovers, the answer to a saved comparison) has no value of its own: the
// test takes the action named by `needs` as the edit, and its result is the
// sample. `value` gives the member as that action leaves it in `drawing`.
//
// `modes` lists the modes in which the member alone is a change Save writes
// (both by default): Save counts a setting of one mode, such as the Linear
// comparison plan, only in a drawing of that mode (the `modes` of its row in
// the mode-scoped registry, gbdraw/web_support/mode_scoped_settings.py), and
// the record rows and File cards only in the mode that draws them. In a
// drawing of the other mode, written for another change, the member is kept
// through every path, Save then Load included, unless `otherModeLoad` is
// 'default': that drawing's Session does not hold it, and its Load gives the
// mode's default. tests/web/drawing-member-paths.test.mjs runs both.
//
// A member that no Session holds, which Save writes into other members
// instead (the records a File card names of its file, `recordsNamed`, which
// Save places into the records OFF), declares `drawn(drawing)`: what the
// drawing draws of it. Save then Load gives the member a new drawing's value
// and keeps `drawn`; every other path keeps the member itself. A slice's own
// Save then Load (tests/web/drawing-slices.test.mjs, `sliceSample`) compares
// members, so there the member keeps the drawing's value.
//
// Import this module after the test has set up `window`.
import { DRAWING_SLICES } from './drawing-state.mjs';
import { normalizeAnnotationSets } from '../../../gbdraw/web/js/services/annotation-state.js';
import { normalizeFeatureVisibilityRule } from '../../../gbdraw/web/js/services/feature-visibility.js';
import { normalizeSpecificRule } from '../../../gbdraw/web/js/services/specific-color-rules.js';
import { normalizePaletteColors } from '../../../gbdraw/web/js/utils/color-utils.js';
import { createCircularSource, createLinearSeq } from '../../../gbdraw/web/js/services/drawing/inputs.js';
import { queuePalette } from '../../../gbdraw/web/js/services/drawing/palette.js';
import { circularRecordsOff, expandedRecordKey } from '../../../gbdraw/web/js/services/record-draw-selection.js';
import { circularRecordEntries } from '../../../gbdraw/web/js/services/record-options.js';
import { defaultSettingsContainers, installStoredSettings, readLinearComparisonPlan, resolvedTrackStacks } from '../../../gbdraw/web/js/services/config-field-readers.js';
import { createDefaultLayoutPreferences, normalizeLayoutPreferences } from '../../../gbdraw/web/js/services/layout-preferences.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../../../gbdraw/web/js/contracts/request-shape.generated.js';
import { LEGEND_POSITIONS, SCALAR_AXIS_LABELS } from '../../../gbdraw/web/js/contracts/schema.generated.js';
import { CIRCULAR_MULTI_RECORD_SIZE_MODES, LINEAR_LABEL_PLACEMENTS, LINEAR_TRACK_LAYOUTS, PROTEIN_BLASTP_MODES } from '../../../gbdraw/web/js/services/current-option-values.js';
import { FIGURE_LENGTH_UNITS } from '../../../gbdraw/web/js/services/figure-size-form.js';
import { FEATURE_RENDERING_VALUES } from '../../../gbdraw/web/js/utils/feature-rendering.js';
import { LINEAR_LABEL_VISIBILITY_MODES } from '../../../gbdraw/web/js/services/linear-label-visibility.js';
import { GC_CONTENT_MODES } from '../../../gbdraw/web/js/mode-profiles.js';
import { CIRCULAR_RECORD_LAYOUTS } from '../../../gbdraw/web/js/services/one-circle-layout.js';
import { CONSERVATION_COLOR_SCHEMES } from '../../../gbdraw/web/js/utils/conservation-colors.generated.js';
import { LINEAR_COMPARISON_MODES, LINEAR_COMPARISON_SOURCES, createLinearComparisonEdge } from '../../../gbdraw/web/js/services/linear-comparisons.js';

/** @typedef {'circular' | 'linear'} Mode */
/** @typedef {import('../../../gbdraw/web/js/state.js').DrawingState} DrawingState */
/** @typedef {{ registerFile: (file: File) => unknown }} SampleContext */
/** @typedef {'records' | 'rings' | 'rows' | 'groups' | 'comparison'} SampleNeed */
/**
 * @typedef {object} MemberSample
 * @property {(mode: Mode, drawing: DrawingState, context: SampleContext) => unknown} value
 * @property {readonly SampleNeed[]} [needs]
 * @property {(drawing: DrawingState) => void} [prepare]
 * @property {readonly Mode[]} [modes]
 * @property {'default'} [otherModeLoad]
 * @property {boolean} [byContext]
 * @property {(drawing: DrawingState) => DrawnRecords} [drawn]
 */
/**
 * What a Circular drawing draws of its records: the records it leaves OFF
 * (`circularRecordsOff`, over its read records) and its Record Order rows.
 * @typedef {{ recordsOff: string[], recordOrder: unknown }} DrawnRecords
 */

// The slice that holds `member`, and the member's current value in `drawing`.
const SLICE_OF = new Map(Object.entries(DRAWING_SLICES).flatMap(([name, slice]) => slice.members.map((member) => [member, name])));
export const sliceOfMember = (member) => SLICE_OF.get(member);
export const captureMember = (drawing, member) => {
  const name = SLICE_OF.get(member);
  return DRAWING_SLICES[name].capture(drawing.slices[name])[member];
};

// The records a 'records' sample names: the drawing's bound cards (Circular)
// or record rows (Linear), by uid.
const boundUids = (drawing, mode) => (mode === 'linear'
  ? captureMember(drawing, 'linearSeqs').filter((row) => row.gb).map((row) => row.uid)
  : captureMember(drawing, 'circularSources').filter((card) => card.gb).map((card) => card.uid));

// A GenBank file of `ids`, one 40 bp record each, with one CDS.
export const genbankFile = (name, ids) => new File([ids.map((id) => [
  `LOCUS       ${id} 40 bp DNA circular`,
  'FEATURES             Location/Qualifiers',
  '     CDS             1..30',
  `                     /locus_tag="${id}_1"`,
  `ORIGIN\n        1 ${'a'.repeat(40)}\n//\n`
].join('\n')).join('')], name, { lastModified: 1 });

// An object member's sample changes every leaf (a leaf is a value that is no
// non-empty plain object), so a path that drops one leaf on its own fails.
// The leaves are the member's own, read from the drawing; a leaf takes the
// first of its candidates that the member's Load reader, given that change
// alone, keeps as a change, in the form the reader gives it (Load resolves a
// track stack's sides, for one). The candidates: the values of the app's own
// lists (`APP_LEAF_CHOICES`), then the values `OTHER_LEAF_CHOICES` names, then
// a boolean's negation, an integer one more or one less, another number's
// half, 7 for null, and for a text the other mode's default or 'sample'. The
// reader must keep every value of an app list (`variedLeaves` reports each one
// it does not, `refused`), so a reader that starts to refuse one valid value
// fails even when another value still varies the leaf. `UNVARIED_LEAVES`
// lists, with the reason, each leaf that keeps the drawing's value;
// tests/web/drawing-member-paths.test.mjs checks that no other leaf does and
// that the reader keeps the whole sample.
const isPlainObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/** The paths of the leaves of `value`. @param {unknown} value @param {string[]} [path] @returns {string[][]} */
export const leafPaths = (value, path = []) => (isPlainObject(value) && Object.keys(value).length
  ? Object.entries(value).flatMap(([key, entry]) => leafPaths(entry, [...path, key])) : [path]);
export const leafAt = (value, path) => path.reduce((entry, key) => (isPlainObject(entry) ? entry[key] : undefined), value);
const withLeaf = (value, [key, ...rest], leaf) => (key === undefined ? leaf : { ...value, [key]: withLeaf(value[key], rest, leaf) });
const sameJson = (left, right) => JSON.stringify(left) === JSON.stringify(right);
const otherMode = (mode) => (mode === 'linear' ? 'circular' : 'linear');
// What Load makes of a stored member, null when it reports a loss: the
// settings merge and its rows (`installStoredSettings`), and for `adv` the
// track stacks as Load resolves them (`applyConfigData`); the layout
// preferences' and the comparison plan's readers.
const settingsRead = (container) => (mode, value) => {
  const targets = defaultSettingsContainers(mode);
  const { losses } = installStoredSettings(targets, { [container]: value }, mode);
  if (losses.length) return null;
  if (container !== 'adv') return targets[container];
  const { circularAxis, circularSlots, linearAxis, linearSlots } = resolvedTrackStacks(targets.adv, targets.form);
  return { ...targets.adv, circular_track_slots_axis_index: circularAxis, circular_track_slots: circularSlots,
    linear_track_slots_axis_index: linearAxis, linear_track_slots: linearSlots };
};
const LOAD_READERS = {
  form: settingsRead('form'), adv: settingsRead('adv'), losat: settingsRead('losat'),
  circularConservation: settingsRead('circularConservation'),
  layoutPreferences: (mode, value) => normalizeLayoutPreferences(value),
  linearComparisonPlan: (mode, value) => {
    const { value: plan, lost } = readLinearComparisonPlan(value);
    return lost ? null : plan;
  }
};
/** The object members whose sample changes every leaf. */
export const LEAF_VARIED_MEMBERS = Object.freeze(Object.keys(LOAD_READERS));
// The choices of the leaves whose reader takes a value only from a list, or
// whose type the default (null, a list) does not show; `*` matches any one key.
const COMPARISON_SETTINGS = CANONICAL_REQUEST_DESCRIPTION.types.CanonicalProteinComparisonSettings.fields;
// A request literal's values, but its alias spellings (a reader reads an
// alias as the value it spells: `either` is `strand`).
const literalValues = (field) => field.values.filter((value) => !Object.hasOwn(field.spellings?.aliases ?? {}, value));
// The values of a list the app exports (its option lists and generated contracts).
const APP_LEAF_CHOICES = {
  form: {
    linear_track_layout: LINEAR_TRACK_LAYOUTS,
    figure_width: FIGURE_LENGTH_UNITS.map((unit) => `180${unit}`), figure_height: FIGURE_LENGTH_UNITS.map((unit) => `120${unit}`)
  },
  adv: {
    'feature_shapes.*': FEATURE_RENDERING_VALUES, label_placement: LINEAR_LABEL_PLACEMENTS,
    linear_accession_visibility: LINEAR_LABEL_VISIBILITY_MODES, linear_length_visibility: LINEAR_LABEL_VISIBILITY_MODES,
    depth_axis_labels: SCALAR_AXIS_LABELS.values, gc_content_axis_labels: SCALAR_AXIS_LABELS.values, gc_content_mode: GC_CONTENT_MODES,
    multi_record_size_mode: CIRCULAR_MULTI_RECORD_SIZE_MODES, circular_record_layout: CIRCULAR_RECORD_LAYOUTS
  },
  losat: {
    'blastp.mode': PROTEIN_BLASTP_MODES,
    'blastp.orthogroupMembershipMode': literalValues(COMPARISON_SETTINGS.orthogroupMembershipMode),
    'blastp.collinearColorMode': literalValues(COMPARISON_SETTINGS.collinearityColorMode),
    'blastp.collinearUnitMode': literalValues(COMPARISON_SETTINGS.collinearityUnitMode),
    'blastp.collinearAnchorMode': literalValues(COMPARISON_SETTINGS.collinearityAnchorMode),
    'blastp.collinearMergeOrientation': literalValues(CANONICAL_REQUEST_DESCRIPTION.types.CanonicalLosslessCollinearityParameters.fields.mergeOrientation),
    'blastp.collinearSearchScope': literalValues(COMPARISON_SETTINGS.collinearitySearchScope)
  },
  layoutPreferences: { 'circular.*.legend': LEGEND_POSITIONS.circular, 'linear.legend': LEGEND_POSITIONS.linear },
  circularConservation: { colorScheme: Object.keys(CONSERVATION_COLOR_SCHEMES) },
  linearComparisonPlan: { mode: Object.values(LINEAR_COMPARISON_MODES), defaultSource: Object.values(LINEAR_COMPARISON_SOURCES) }
};
// Literals where the app exports no list, and values built from the leaf or
// the drawing (a function of `(leaf, drawing)`).
const SAMPLE_WINDOW = { start: 5, end: 25, reverseComplement: false, column: null, acrossOrigin: false };
const COLOR = ['#123456'];
const OTHER_LEAF_CHOICES = {
  form: {
    // literals: the app exports no list of these
    track_type: ['spreadout', 'middle'], scale_style: ['ruler'], labels_mode: ['out'], show_labels_linear: ['all'],
    circular_record_selector: ['#2'], circular_range_fields_record: ['REC_A'],
    circular_region_start: [10], circular_region_end: [30], linear_wrap_length: [5000],
    // a window as `readLinearSeqWindows` gives it
    circular_windows: [[SAMPLE_WINDOW]],
    circular_record_ranges: [[{ source: 'sample-card', record: 'REC_A', windows: [SAMPLE_WINDOW] }]]
  },
  adv: {
    features: (leaf) => [leaf.slice(0, -1)], arrow_shaft_width_ratio: [0.5],
    // literals: the app exports no list of these
    font_family: ['serif'], nt: ['AT'], evalue: ['1e-10'], label_rendering: ['embedded_only'], circular_label_placement: ['radial'],
    circular_grouping_intent: ['circular'],
    'linear_definition_line_styles.*.font_weight': ['bold'], 'linear_definition_line_styles.*.fill': COLOR,
    block_stroke_color: COLOR, line_stroke_color: COLOR, axis_stroke_color: COLOR, scale_stroke_color: COLOR,
    ruler_label_color: COLOR, depth_color: COLOR,
    depth_tracks: [[{ label: 'Sample depth', color: '#123456', height: 9, large_tick_interval: 9, small_tick_interval: 3, tick_font_size: 9 }]],
    feature_density_tracks: [[{ id: 'density-1', query: { featureType: 'CDS', qualifier: null, value: null, minLength: null, maxLength: null },
      metric: 'count', style: 'heatmap', window: 1000, step: 500, color: null, caption: 'Sample density' }]],
    identity_scale: [[40, 90]], multi_record_positions: [[{ selector: '#1', row: 2 }]],
    // What choosing one circle keeps of the layout it leaves (`chooseRecordLayout`).
    one_circle_left: () => [{ layout: 'separate', gc_content_mode: GC_CONTENT_MODES.at(-1), suppress_skew: true }],
    circular_track_slots_axis_index: [1], linear_track_slots_axis_index: [0],
    // A row of each stack changed: the Features row's size.
    linear_track_slots: (slots) => [slots.map((slot) => (slot.id === 'features' ? { ...slot, height: '12' } : slot))],
    circular_track_slots: (slots) => [slots.map((slot) => (slot.id === 'features' ? { ...slot, width: 12 } : slot))]
  },
  // literals: the app exports no list of these
  losat: { outfmt: ['7'], 'blastn.task': ['blastn'] },
  // literals: the app exports no list of the title positions
  layoutPreferences: { 'circular.*.plotTitlePosition': ['top', 'bottom', 'none'], 'linear.plotTitlePosition': ['center', 'top', 'bottom'] },
  circularConservation: {
    // literals: the app exports no list of these (the reader's `textChoice` lists, config-field-readers.js)
    source: ['upload'], losat_program: ['tblastx'], reference: ['query', 'subject'],
    // Each ring row of the drawing's ring tables (the 'rings' context) renamed and colored.
    series: (rows) => (rows.length ? [rows.map((row) => ({ ...row, label: 'Sample ring', color: '#123456' }))] : [])
  },
  linearComparisonPlan: {
    // An edge that pairs the drawing's first two record rows (the 'rows' context), as the pair card makes one.
    edges: (edges, drawing) => {
      const [query, subject] = boundUids(drawing, 'linear');
      return subject ? [[createLinearComparisonEdge({ id: 'sample-edge', queryUid: query, subjectUid: subject, included: true,
        source: LINEAR_COMPARISON_SOURCES.LOSAT, losatFilename: 'sample-pair.tsv', losatFilenameActive: true })]] : [];
    }
  }
};
/** The key of `table` whose choices a leaf of `member` at `path` takes. @param {Record<string, object>} table @param {string} member @param {string[]} path */
const choiceKeyAt = (table, member, path) => Object.keys(table[member] ?? {}).find((name) => {
  const parts = name.split('.');
  return parts.length === path.length && parts.every((part, index) => part === '*' || part === path[index]);
});
/**
 * Every key of the choice tables, as `<member> <key>`, and the keys a leaf of
 * `member` at `path` takes: a key no leaf takes in either mode names a leaf
 * the member no longer has (tests/web/drawing-member-paths.test.mjs).
 */
export const LEAF_CHOICE_KEYS = Object.freeze([APP_LEAF_CHOICES, OTHER_LEAF_CHOICES]
  .flatMap((table) => Object.entries(table).flatMap(([member, keys]) => Object.keys(keys).map((key) => `${member} ${key}`))));
/** @param {string} member @param {string[]} path @returns {string[]} */
export const leafChoiceKeysAt = (member, path) => [APP_LEAF_CHOICES, OTHER_LEAF_CHOICES]
  .map((table) => choiceKeyAt(table, member, path)).filter((key) => key !== undefined).map((key) => `${member} ${key}`);
const choicesAt = (table, member, path, leaf, drawing) => {
  const key = choiceKeyAt(table, member, path);
  const choices = key === undefined ? [] : table[member][key];
  return typeof choices === 'function' ? choices(leaf, drawing) : [...choices];
};
const leafCandidates = (member, mode, path, leaf, drawing) => {
  const named = [...choicesAt(APP_LEAF_CHOICES, member, path, leaf, drawing), ...choicesAt(OTHER_LEAF_CHOICES, member, path, leaf, drawing)];
  if (typeof leaf === 'boolean') return [...named, !leaf];
  if (typeof leaf === 'number') return [...named, ...(Number.isInteger(leaf) ? [leaf + 1, leaf - 1] : [leaf / 2])];
  if (leaf === null) return [...named, 7];
  if (typeof leaf === 'string') return [...named, leafAt(defaultSettingsContainers(otherMode(mode))[member], path), 'sample'];
  return named;
};

/**
 * The leaves of an object member that keep the drawing's value, by member, with
 * their reason; a function of the mode lists the leaves of that mode. A leaf is
 * here only when it cannot vary alone.
 * @type {Record<string, Record<string, string> | ((mode: Mode) => Record<string, string>)>}
 */
export const UNVARIED_LEAVES = {
  adv: {
    // Load writes the current version of each stack's rows (`applyConfigData`).
    linear_track_slots_schema_version: 'the current version of the stack rows',
    circular_track_slots_schema_version: 'the current version of the stack rows'
  },
  // The request schema's membership modes: while it has one, the leaf has no other value.
  losat: literalValues(COMPARISON_SETTINGS.orthogroupMembershipMode).length === 1
    ? { 'blastp.orthogroupMembershipMode': 'the request schema has one membership mode' } : {},
  circularConservation: (mode) => ({
    // The label list an older Session stored for rings past its series: Load
    // reads it into the series (`withoutUnreadConservationSeries`), every ring
    // reconcile clears it (OV-3166), and the settings reader keeps none.
    labels: 'an older Session\'s label list, which Load reads into the series',
    // A ring row names a table of the rings' source: a change of the source
    // alone derives the rows again (`changeCircularConservationSource`), so
    // the sample keeps the source of its ring table ('upload', no default).
    ...(mode === 'circular' ? { source: 'the source of the ring table the rows name' }
      : { series: 'a Linear drawing has no comparison rings' })
  }),
  linearComparisonPlan: (mode) => (mode === 'circular' ? { edges: 'a Circular drawing has no Linear record rows to pair' } : {}),
  // A drawing's mode never changes, and its Session holds the layout slot of
  // its own mode only.
  layoutPreferences: (mode) => Object.fromEntries(leafPaths(createDefaultLayoutPreferences())
    .filter(([slot]) => slot !== mode).map((path) => [path.join('.'), 'the other mode\'s slot']))
};
export const unvariedOf = (member, mode) => {
  const entries = UNVARIED_LEAVES[member] ?? {};
  return typeof entries === 'function' ? entries(mode) : entries;
};

/**
 * `value` (an object member of `drawing`, of `mode`) with every leaf changed as
 * above; the paths of the leaves no candidate changes; the values of an app
 * list the Load reader does not keep (`refused`); and what the reader makes of
 * the result.
 * @param {string} member
 * @param {Mode} mode
 * @param {unknown} value
 * @param {DrawingState} drawing
 */
export const variedLeaves = (member, mode, value, drawing) => {
  const read = LOAD_READERS[member];
  const kept = unvariedOf(member, mode);
  const readingAt = (path, candidate) => leafAt(read(mode, withLeaf(value, path, candidate)), path);
  let varied = value;
  const unchanged = [];
  const refused = [];
  for (const path of leafPaths(value)) {
    const leaf = leafAt(value, path);
    choicesAt(APP_LEAF_CHOICES, member, path, leaf, drawing).filter((candidate) => !sameJson(readingAt(path, candidate), candidate))
      .forEach((candidate) => refused.push(`${path.join('.')} = ${JSON.stringify(candidate)}`));
    const changed = Object.hasOwn(kept, path.join('.')) ? undefined : leafCandidates(member, mode, path, leaf, drawing)
      .map((candidate) => readingAt(path, candidate))
      .find((reading) => reading !== undefined && !sameJson(reading, leaf));
    if (changed === undefined) unchanged.push(path.join('.'));
    else varied = withLeaf(varied, path, changed);
  }
  // A leaf whose reading follows others (a stack row's side follows the
  // stack's axis index) takes its reading of the whole sample, while that is
  // still a change.
  for (let round = 0; round < 3; round += 1) {
    const whole = read(mode, varied);
    const follows = whole ? leafPaths(varied).filter((path) => !sameJson(leafAt(whole, path), leafAt(varied, path))
      && !sameJson(leafAt(whole, path), leafAt(value, path))) : [];
    if (!follows.length) break;
    follows.forEach((path) => { varied = withLeaf(varied, path, leafAt(whole, path)); });
  }
  return { varied, unchanged, refused, read: read(mode, varied) };
};

// Every file slot of a drawing (the slots its `files` capture holds) is
// bound: a slot that holds a file or nothing (null) a new file, a list one
// file, unless `FILE_SLOT_SAMPLES` gives the slot's form. A Circular-only
// slot keeps its value in a Linear drawing: Save writes a Linear drawing's
// record rows, not these slots (`serializeActiveRenderFiles`, config.js), and
// a Linear drawing has no control that binds them.
const FILE_SLOT_SAMPLES = {
  // A Depth slot holds one record's series; each ring place a list of genome files.
  c_depth: (file) => [file()],
  c_conservation_fastas: (file) => [[file()]],
  c_conservation_sequence_sources: (file) => [[file()]],
  // The flag of a LOSAT-cache replay: its BLAST rows are the cache's (whose
  // History intent holds them, B23).
  c_conservation_blasts_source: () => 'losat-cache',
  // What a Generate wrote (the Linear comparisons, which an artifact
  // checkpoint holds, not the History intent) keeps its value; so do the
  // cards the rows of a record-major Depth matrix stand for, which only a
  // matrix has. The cards an older Session's flat record positions count over
  // (OV-2867) are marked by the Load of a Session 27-50 only: a Save writes
  // none, so the sample holds none, as a Load of a Session 51 gives.
  linearCanonicalComparisons: null,
  c_depth_rows: null,
  c_legacy_record_cards: () => null
};
const CIRCULAR_FILE_SLOTS = Object.freeze(['c_depth', 'c_conservation_blasts', 'c_conservation_blasts_source', 'c_conservation_fastas',
  'c_conservation_sequence_sources']);
const filesSample = (mode, drawing, { registerFile }) => Object.fromEntries(Object.entries(captureMember(drawing, 'files')).map(([slot, bound]) => {
  const file = () => registerFile(new File(['CDS\t#123456\n'], `sample-${mode}-${slot}.tsv`, { lastModified: 1 }));
  if ((mode === 'linear' && CIRCULAR_FILE_SLOTS.includes(slot)) || FILE_SLOT_SAMPLES[slot] === null) return [slot, bound];
  if (FILE_SLOT_SAMPLES[slot]) return [slot, FILE_SLOT_SAMPLES[slot](file)];
  return [slot, Array.isArray(bound) ? [file()] : file()];
}));

/** @type {Record<string, MemberSample>} */
export const MEMBER_SAMPLES = {
  // visibility rules
  featureVisibilityManualRules: { value: () => [normalizeFeatureVisibilityRule({
    id: 'sample-rule', featureType: 'CDS', qualifier: 'locus_tag', value: '^A$', action: 'off'
  })] },
  // color rules: one rule a table file loaded and one the user added, which
  // holds no `fromFile` (`normalizeSpecificRule`'s shape, #984).
  manualSpecificRules: { value: () => [
    normalizeSpecificRule({ feat: 'CDS', qual: 'product', val: 'kinase', color: '#123456', cap: 'Kinase', fromFile: true }),
    normalizeSpecificRule({ feat: 'tRNA', qual: 'product', val: 'Leu', color: '#654321', cap: 'Leucine' })
  ] },
  // palette: Load gives a palette's colors the comparison colors (`normalizePaletteColors`).
  selectedPalette: { value: () => 'mint' },
  currentColors: { value: () => normalizePaletteColors({ CDS: '#123456', tRNA: '#654321' }) },
  pendingPaletteName: { value: () => 'queued' },
  // A palette queued while Instant Preview is off holds its name and the
  // colors it was queued with (`queuePalette`); the sample changes one.
  pendingPaletteColors: { prepare: (drawing) => queuePalette(drawing.slices.palette, 'queued', drawing.currentColors.value),
    value: (mode, drawing) => ({ ...captureMember(drawing, 'pendingPaletteColors'), CDS: '#abcdef' }) },
  // label tables
  labelTextBulkOverrides: { value: () => ({ 'hypothetical protein': 'HP' }) },
  canonicalLabelOverrideRows: { value: () => [
    { recordId: '*', featureType: 'CDS', qualifier: 'product', valueRegex: '^kinase$', labelText: 'Kinase' }
  ] },
  manualPriorityRules: { value: () => [{ feat: 'CDS', order: 'gene,product' }] },
  filterMode: { value: () => 'Whitelist' },
  manualBlacklist: { value: () => 'putative' },
  manualWhitelist: { value: () => [{ feat: 'CDS', qual: 'product', key: 'kinase' }] },
  // annotation sets, as `normalizeAnnotationSets` gives them
  annotationSets: { value: () => normalizeAnnotationSets([{ id: 'sample-set', legendLabel: 'Sample', annotations: [
    { id: 'region_1', target: { kind: 'coordinateSpan', start: 10, end: 20 }, label: 'Sample region' }
  ] }]) },
  // group names (Linear only): Load keeps a name or description only for a
  // group of the drawing's saved Result (`orthogroupEditsOf`, B47); a dormant
  // name, by member set, needs none.
  orthogroupNameOverrides: { modes: ['linear'], otherModeLoad: 'default', needs: ['groups'], value: () => ({ og_1: 'Sample group' }) },
  orthogroupDescriptionOverrides: { modes: ['linear'], otherModeLoad: 'default', needs: ['groups'], value: () => ({ og_1: 'Sample description' }) },
  orthogroupDormantOverrides: { modes: ['linear'], value: () => ({ '["sample-member"]': { name: 'Dormant name', description: 'Dormant description' } }) },
  // canvas
  canvasPadding: { value: () => ({ top: 12, right: 3, bottom: 7, left: 5 }) },
  // feature edits: Load keeps the rows `canonicalFeatureOverrides` and the
  // placement reader of the drawing's mode accept.
  featureOverrides: { value: () => ({ '["record-1","CDS_1"]': {
    recordKey: 'record-1', biologicalFeatureId: 'CDS_1', featureVisibility: 'off', labelVisibility: null, labelText: 'Sample', labelSourceText: null
  } }) },
  featurePlacementOverrides: { value: (mode) => ({ '["record-1","CDS_2"]': {
    recordKey: 'record-1', biologicalFeatureId: 'CDS_2', placement: { kind: 'lane', side: mode === 'linear' ? 'above' : 'outward', level: 1 }
  } }) },
  // paint: Load reads the Legend rows as `normalizeEditorStateData` gives them.
  featureColorOverrides: { value: () => ({ '["record-1","CDS_1"]': { color: '#123456', caption: 'Kinase' } }) },
  featureStrokeOverrides: { value: () => ({ 'record-1\u0000CDS_1': { strokeColor: '#654321', strokeWidth: 2 } }) },
  legendEntries: { value: () => [{ caption: 'Renamed CDS', originalCaption: 'CDS', color: '#abcdef', featureIds: [] }] },
  deletedLegendEntries: { value: () => [{ caption: 'tRNA', originalCaption: 'tRNA', color: '#fedcba', featureIds: [] }] },
  // A dormant row is Python's: it holds no color (V-OD22).
  dormantLegendEntries: { value: () => [{ caption: 'Coverage', originalCaption: 'depth', featureIds: [] }] },
  legendColorOverrides: { value: () => ({ CDS: '#abcdef' }) },
  legendStrokeOverrides: { value: () => ({ CDS: { strokeColor: '#222222', strokeWidth: 1.5 } }) },
  addedLegendCaptions: { value: () => ['Kinase'] },
  // settings: an object member changes every leaf (`variedLeaves`)
  form: { value: (mode, drawing) => variedLeaves('form', mode, captureMember(drawing, 'form'), drawing).varied },
  adv: { value: (mode, drawing) => variedLeaves('adv', mode, captureMember(drawing, 'adv'), drawing).varied },
  losat: { value: (mode, drawing) => variedLeaves('losat', mode, captureMember(drawing, 'losat'), drawing).varied },
  losatProgram: { modes: ['linear'], value: () => 'tblastx' },
  // The ring rows of the drawing's ring table change with the other leaves.
  circularConservation: { modes: ['circular'], needs: ['rings'],
    value: (mode, drawing) => variedLeaves('circularConservation', mode, captureMember(drawing, 'circularConservation'), drawing).varied },
  // An edge pairs the drawing's two record rows.
  linearComparisonPlan: { modes: ['linear'], needs: ['rows'],
    value: (mode, drawing) => variedLeaves('linearComparisonPlan', mode, captureMember(drawing, 'linearComparisonPlan'), drawing).varied },
  linearRecordLayoutEnabled: { modes: ['linear'], value: (mode, drawing) => !captureMember(drawing, 'linearRecordLayoutEnabled') },
  linearRecordGap: { modes: ['linear'], value: () => 40 },
  linearRecordRows: { modes: ['linear'], needs: ['records'], value: (mode, drawing) => boundUids(drawing, mode).slice(0, 1).map((uid) => ({ uid, row: 2 })) },
  recordDisplayDrafts: { needs: ['records'], value: (mode, drawing) => boundUids(drawing, mode).slice(0, 1).map((uid) => ({
    sourceUid: uid, selector: '#1', recordId: 'REC_A', topologyOverride: null, startCoordinate: 7,
    reverseComplementOverride: null, anchorIntent: null
  })) },
  // A Linear record row is turned off by its ID, a Circular record by its
  // card's record key (a flat `#N` is the form of a Session 46-50, OV-2869).
  recordsOff: { needs: ['records'], value: (mode, drawing) => boundUids(drawing, mode).slice(0, 1)
    .map((uid) => (mode === 'linear' ? uid : expandedRecordKey(uid, 2))) },
  // The first record of the first card's file, named as the request a Load
  // or a copy installs names it (its request key and selector, OV-3530): the
  // card draws that record alone. The card is read, as after an Undo across
  // the read that placed the named records (`circularRecordsOff`); Save
  // places them into the records OFF (`placeSavedNamedRecords`, config.js).
  recordsNamed: { modes: ['circular'], otherModeLoad: 'default', needs: ['records'],
    value: (mode, drawing) => boundUids(drawing, mode).slice(0, 1).map((uid) => ({
      uid, named: [{ recordKey: 'record-1', selector: { kind: 'recordIndex', index: 0 } }]
    })),
    drawn: (drawing) => ({
      recordsOff: circularRecordsOff(drawing, circularRecordEntries(drawing.circularRecordList.value)),
      recordOrder: captureMember(drawing, 'adv').multi_record_positions
    }) },
  // A configuration path the Python evaluator accepts.
  unmanagedConfigOverrides: { value: () => ({ 'objects.gc_content.percent_background_opacity': 0.42 }) },
  // The answer to a saved comparison, which only the drawing of a Linear
  // Result holds: Inherit saved comparison.
  importedComparisonIntent: { modes: ['linear'], needs: ['comparison'], byContext: true,
    value: (mode, drawing) => captureMember(drawing, 'importedComparisonIntent') },
  // A drawing keeps the layout preferences of its own mode.
  layoutPreferences: { value: (mode, drawing) => variedLeaves('layoutPreferences', mode, captureMember(drawing, 'layoutPreferences'), drawing).varied },
  linearTypographyLinked: { modes: ['linear'], value: (mode, drawing) => !captureMember(drawing, 'linearTypographyLinked') },
  fileImports: { value: () => ({ blacklist: { words: ['sample'], before: 'putative' } }) },
  // inputs: a file bound in the drawing's own input fields, which the store
  // registers. A Circular drawing binds its GenBank file in a File card, a
  // Linear drawing in its one record row.
  // A drawing whose only change is a card file is saved (OV-3701).
  files: { value: filesSample },
  // Save writes a drawing's record rows (Linear) or File cards (Circular),
  // not both (`serializeActiveRenderFiles` in services/config.js).
  linearSeqs: { modes: ['linear'], otherModeLoad: 'default', value: (mode, drawing, { registerFile }) => [createLinearSeq({
    uid: 'sample-row', gb: registerFile(genbankFile(`sample-${mode}.gb`, [`sample_${mode}`]))
  })] },
  // The card holds a file for each input type, so it is a complete input
  // whichever type the drawing reads (`inputType`'s sample is GFF3).
  circularSources: { modes: ['circular'], otherModeLoad: 'default', value: (mode, drawing, { registerFile }) => [createCircularSource({
    uid: 'sample-card', gb: registerFile(genbankFile(`sample-${mode}.gb`, [`sample_${mode}`])),
    gff: registerFile(new File([`##gff-version 3\nsample_${mode}\tsample\tCDS\t1\t30\t.\t+\t0\tID=cds1\n`], `sample-${mode}.gff`, { lastModified: 1 })),
    fasta: registerFile(new File([`>sample_${mode}\n${'a'.repeat(40)}\n`], `sample-${mode}.fasta`, { lastModified: 1 }))
  })] },
  // A drawing whose only change is its input type is saved (OV-3701).
  inputType: { value: () => 'gff' },
  // What the record read of the drawing's File cards finds.
  circularRecordList: { modes: ['circular'], needs: ['records'], byContext: true, value: (mode, drawing) => captureMember(drawing, 'circularRecordList') },
  circularRecordDiscovery: { modes: ['circular'], needs: ['records'], byContext: true,
    value: (mode, drawing) => captureMember(drawing, 'circularRecordDiscovery') }
};

// The contexts a test builds in the drawing itself; a sample of 'rings' or
// 'rows' changes the leaves the drawing's tables and rows allow.
export const DRAWING_NEEDS = Object.freeze(['records', 'rings', 'rows']);
// A slice's sample: the merge of its members' samples, after their contexts
// (`prepare`). A member keeps the drawing's value when its sample needs a
// fixture's Result, or bound records the drawing does not have, when only an
// action writes it, when a Session of `mode` does not hold it
// (`otherModeLoad`), or when Save writes it into other members (`drawn`).
export const sliceSample = (name, mode, drawing, context) => Object.fromEntries(DRAWING_SLICES[name].members.map((member) => {
  const sample = MEMBER_SAMPLES[member];
  const needs = sample.needs ?? [];
  const kept = sample.byContext || Boolean(sample.drawn) || needs.some((need) => !DRAWING_NEEDS.includes(need))
    || (needs.includes('records') && boundUids(drawing, mode).length === 0)
    || (!(sample.modes ?? ['circular', 'linear']).includes(mode) && sample.otherModeLoad === 'default');
  if (!kept) sample.prepare?.(drawing);
  return [member, kept ? captureMember(drawing, member) : sample.value(mode, drawing, context)];
}));

// A capture as JSON without the store's file IDs: a Load binds each input
// file again, under a new ref of the same file. A record row without a file
// is a placeholder (a Circular drawing saves no rows), new at each Load. An
// object's keys are sorted: an object member is a map read by key (a color by
// feature key, a record row's fields), and a Load or History restore may
// build it in another key order. No member shows its key order (what the user
// orders, such as the Legend rows or the track slots, is a list), so the sort
// hides no loss; a list keeps its order.
export const refFree = (value) => JSON.stringify(value, (key, entry) => {
  if (entry && typeof entry === 'object' && typeof entry.fileId === 'string' && typeof entry.sourceId === 'string') {
    return { name: entry.name, type: entry.type, size: entry.size, lastModified: entry.lastModified };
  }
  if (key === 'linearSeqs') return entry.map((row) => (row.gb || row.gff || row.fasta ? row : { ...row, uid: '' }));
  if (entry && typeof entry === 'object' && !Array.isArray(entry)) {
    return Object.fromEntries(Object.keys(entry).sort().map((name) => [name, entry[name]]));
  }
  return entry;
});
