// @ts-check
/** @import { FileImports } from './file-imports.js' */
/** @import { FeaturePlacementDraft } from './feature-placement.js' */
/** @import { RawSessionResourceDescriptor } from './session-resource-backing.js' */
/** @import { ProteinBlastpMode } from './current-option-values.js' */
/** @import { CollinearMergeOrientation } from './losat-normalization.js' */
/** @import { CanonicalProteinComparisonSettings } from '../contracts/request.generated.js' */
/** @import { CanonicalFeatureDensityTrack, CanonicalLinearLayout, CanonicalLinearRenderRequest, CanonicalSimilarityAlignmentPlan } from '../contracts/request.generated.js' */
/** @import { FieldLeaf, FieldReading } from '../utils/field-reading.js' */
/** @import { LinearWrapFormValue } from './linear-wrap-form.js' */
/** @import { FigureLengthFormValue } from './figure-size-form.js' */
import { adoptCurrentSessionResources, createSessionResourceFileView, readSessionResourceBytes } from './session-resource-backing.js';
import { sha256Hex, textToBytes } from './byte-utils.js';
import { canonicalRecordReverseComplement } from './record-display-model.js';
import { createDefaultLinearDefinitionLineStyles } from './definition-line-style-state.js'; import { CIRCULAR_TRACK_RENDERERS, createDefaultCircularTrackSlots } from './circular-track-slot-model.js';
import { LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_RENDERERS, LINEAR_TRACK_SLOT_SCHEMA_VERSION, createDefaultLinearTrackSlots } from './linear-track-slot-model.js'; import { trackSlotLeafFaults, validateTrackSlotBindingInvariants } from './track-slot-validation.js';
import { requireCurrentCollinearAnchorMode, requireCurrentCollinearColorMode, requireCurrentCollinearMaxConflicts, requireCurrentCollinearMaxDiagonalDrift, requireCurrentCollinearMaxParalogLinks, requireCurrentCollinearMaxUnitGap, requireCurrentCollinearMergeOrientation, requireCurrentCollinearMinAnchors, requireCurrentCollinearInferOrthogroups, requireCurrentCollinearSearchScope, requireCurrentCollinearUnitMode, requireCurrentOrthogroupMemberMaxHits, requireCurrentOrthogroupMembershipMode, requireCurrentProteinBlastpCandidateLimit, requireCurrentProteinBlastpMaxHits, requireCurrentProteinBlastpMode, requireCurrentConfigFields } from './current-option-values.js'; import { DEFAULT_ARROW_SHAFT_WIDTH_RATIO, createDefaultFeatureRenderings } from '../utils/feature-rendering.js';
import { MODE_DEFAULT_FEATURE_TYPES, comparisonStateForMode, labelStyleForMode, managedAdvStateForMode, trackDefaultsForMode } from '../mode-profiles.js'; import { WEB_UX_PROFILE } from '../web-ux-profile.js';
import { assertSafeObjectKeys } from './safe-object-keys.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { INPUT_TYPES } from './input-types.js';
import { CONSERVATION_DEFAULT_COLOR_SCHEME } from '../utils/conservation-colors.generated.js';
/** @import { InputType, LinearSeqWindow } from './drawing/inputs.js' */
/** @import { CircularRecordRange } from './record-options.js' */
import { isUnsetValue } from '../utils/field-reading.js';
import { depthTrackEditDefaults, ensureDepthTrackConfigCount } from './depth-track-state.js';
import { validateRecordDisplayDrafts } from './record-display-model.js';
import { isRecordDrawKey } from './record-draw-selection.js';
import { canonicalFeaturePlacements } from './feature-placement.js';
import { validateScopedFeaturePlacements } from './feature-edit-migration.js';
import { FRESH_PROTEIN_COMPARISON_DEFAULTS, ORTHOGROUP_INFERENCE_VERSION, SCALAR_AXIS_LABELS, SESSION_PNG_DPI_CHOICES } from '../contracts/schema.generated.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../contracts/request-shape.generated.js';
import { ONE_CIRCLE_ADV_DEFAULTS } from './one-circle-layout.js';
/** @import { ObjectShape } from '../utils/contract-shape.js' */
import { DEFAULT_COLLINEAR_ANCHOR_MODE, DEFAULT_MERGE_ORIENTATION } from './losat-normalization.js';
const circularTracks = trackDefaultsForMode('circular'), linearTracks = trackDefaultsForMode('linear');
// The Feature overlap tolerance of a new drawing, in bp.
export const DEFAULT_FEATURE_OVERLAP_TOLERANCE_BP = 0;
/**
 * The Feature overlap tolerance a draft draws and saves with: a blank (a
 * cleared field holds null) or an unset value is the default, as R7
 * reads a blank (OV-826, BA-OD112). Any other value is the draft's own, which
 * `validateCurrentWriterActiveConfig` and Python check.
 * @param {{ feature_overlap_tolerance_bp?: unknown }} adv
 * @returns {unknown}
 */
export const draftFeatureOverlapTolerance = ({ feature_overlap_tolerance_bp: value }) => (
  value === undefined || value === null || value === '' ? DEFAULT_FEATURE_OVERLAP_TOLERANCE_BP : value
);
export const CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 4, LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 3;
/**
 * @template T
 * @param {T} value
 * @returns {value is T & Record<string, unknown>}
 */
const isObject = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value), has = /** @type {(value: unknown, key: string) => boolean} */ ((value, key) => Object.prototype.hasOwnProperty.call(value, key));
/**
 * The Session's active Web configuration (`config` of the stored Session) in
 * its current writer format: the domains of `CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS`
 * plus the compatibility field `colorsAreOverrides`. `form` and `adv` are
 * required; the option fields inside them belong to Python (R7) and stay
 * `Record<string, any>`. A reader takes unvalidated data and returns this type.
 * @typedef {object} ActiveWebConfig
 * @property {Record<string, any>} form
 * @property {Record<string, any>} adv
 * @property {Record<string, any>} [losat]
 * @property {Record<string, any>} [cliOptions]
 * @property {Record<string, any>} [colors]
 * @property {string} [palette]
 * @property {boolean} [paletteInstantPreviewEnabled]
 * @property {Array<Record<string, any>>} [rules]
 * @property {Array<Record<string, any>>} [qualifierPriorityRules]
 * @property {string} [filterMode] `None`, `Whitelist`, or `Blacklist`.
 * @property {Array<Record<string, any>>} [whitelist]
 * @property {string} [blacklistText]
 * @property {FileImports} [fileImports] What the Default colors and Blacklist files set (OV-1280).
 * @property {string} [losatProgram]
 * @property {Record<string, any>} [circularConservation]
 * @property {Array<Record<string, any>>} [annotationSets]
 * @property {Array<Record<string, any>>} [recordDisplayDrafts]
 * @property {FeaturePlacementDraft} [featurePlacementOverrides]
 * @property {Record<string, any>} [modeProfiles]
 * @property {Record<string, any>} [unmanagedConfigOverrides]
 * @property {Record<string, any>} [linearRecordLayout]
 * @property {Record<string, any>} [linearComparisonPlan]
 * @property {Record<string, any>} [importedComparisonResolution]
 * @property {Record<string, any>} [webEdits]
 * @property {boolean} [colorsAreOverrides]
 */

export const createDefaultForm = () => ({
  prefix: '', species: '', strain: '', plot_title: '', track_type: 'tuckin', linear_track_layout: 'middle', show_scale: true, scale_style: 'bar',
  linear_ruler_on_axis: false, labels_mode: 'none', show_labels_linear: 'none', multi_record_canvas: WEB_UX_PROFILE.circular.gridByDefault,
  circular_record_selector: '', circular_region_start: null, circular_region_end: null, circular_reverse: false,
  // The window across the origin a CLI or Python Session gave the record, which the Start/End fields cannot hold
  // (`cardCropWindow`); shown read-only.
  circular_windows: /** @type {LinearSeqWindow[]} */ ([]),
  // The crop a CLI or Python Session gave each record of a Circular grid or batch, drawn read-only there and put
  // into the fields above when Record select chooses the record alone (PD-OI-044 rev 5).
  circular_record_ranges: /** @type {CircularRecordRange[]} */ ([]), circular_range_fields_record: '',
  circular_record_label: '', circular_record_subtitle: '',
  separate_strands: WEB_UX_PROFILE.separateStrands, suppress_gc: !circularTracks.gc, suppress_skew: !circularTracks.skew, align_center: false,
  keep_definition_left_aligned: true, show_gc: linearTracks.gc, show_skew: linearTracks.skew, show_depth: false, normalize_length: false,
  // Off. Wrap Rows writes it through the layout transition (services/linear-wrap-form.js `linearWrapFormPatch`).
  linear_wrap_length: /** @type {LinearWrapFormValue} */ (null),
  // Auto. The Figure size field of the Layout card writes each side (app/figure-size-field.js).
  figure_width: /** @type {FigureLengthFormValue} */ (null),
  figure_height: /** @type {FigureLengthFormValue} */ (null)
});
// A new drawing's Depth settings. It holds one Depth series, the default
// series of these settings, as 0.14.0 stored it in every mode (OV-1901,
// D-15-122); the drawing, Reset Settings and the generated mode defaults read
// it here.
const DEFAULT_DEPTH_SETTINGS = Object.freeze({ depth_height: null, depth_color: '#4A90E2' });
/**
 * @param {string} [mode] `circular` or `linear`.
 * @returns {Record<string, any>}
 */
export const createDefaultAdv = (mode = 'circular') => ({
  features: [...MODE_DEFAULT_FEATURE_TYPES], feature_shapes: createDefaultFeatureRenderings(), arrow_head_length_ratio: null,
  arrow_shaft_width_ratio: DEFAULT_ARROW_SHAFT_WIDTH_RATIO, font_family: null, window_size: null, step_size: null, nt: 'GC', def_font_size: null,
  circular_definition_interval: null, label_font_size: null, circular_label_spacing: null, linear_label_spacing: null, label_rendering: 'auto',
  ...labelStyleForMode(mode), circular_label_placement: 'horizontal', label_placement: 'auto', label_rotation: null, block_stroke_width: null, block_stroke_color: null,
  line_stroke_width: null, line_stroke_color: null, axis_stroke_width: null, axis_stroke_color: managedAdvStateForMode(mode).axis_stroke_color,
  legend_box_size: null, legend_font_size: null, resolve_overlaps: false, feature_overlap_tolerance_bp: DEFAULT_FEATURE_OVERLAP_TOLERANCE_BP, feature_height: null, track_axis_gap: null, linear_show_replicon: false,
  linear_accession_visibility: 'auto', linear_length_visibility: 'auto', linear_definition_line_styles: createDefaultLinearDefinitionLineStyles(), gc_height: null,
  ...DEFAULT_DEPTH_SETTINGS, depth_tracks: ensureDepthTrackConfigCount([], 1, depthTrackEditDefaults(DEFAULT_DEPTH_SETTINGS)), depth_window_size: null, depth_step_size: null, depth_share_axis: false,
  depth_min: null, depth_max: null, depth_normalize: false, depth_show_axis: true, depth_show_ticks: true, depth_axis_labels: SCALAR_AXIS_LABELS.default, depth_large_tick_interval: null,
  depth_small_tick_interval: null, depth_tick_font_size: null,
  // The feature density tracks a Session's request gave the drawing (services/feature-density-tracks.js).
  feature_density_tracks: /** @type {CanonicalFeatureDensityTrack[]} */ ([]), linear_track_slots_enabled: false, linear_track_slots_schema_version: LINEAR_TRACK_SLOT_SCHEMA_VERSION,
  linear_track_slots_axis_index: null, linear_track_slots: createDefaultLinearTrackSlots(), gc_content_mode: 'deviation', gc_content_min_percent: 0,
  gc_content_max_percent: 100, gc_content_show_axis: true, gc_content_show_ticks: true, gc_content_axis_labels: SCALAR_AXIS_LABELS.default, gc_content_tick_interval: 20, gc_content_small_tick_interval: null,
  gc_content_tick_font_size: null, comparison_height: null, pairwise_match_style: mode === 'linear' ? 'curve' : 'ribbon', ...comparisonStateForMode(mode), identity_scale: null, scale_interval: null,
  scale_font_size: null, ruler_label_font_size: null, scale_stroke_width: null, scale_stroke_color: null, ruler_label_color: null, circular_grouping_intent: 'auto',
  multi_record_size_mode: 'auto', multi_record_min_radius_ratio: 0.55, multi_record_column_gap_ratio: 0.10, multi_record_row_gap_ratio: 0.05,
  multi_record_positions: [], ...ONE_CIRCLE_ADV_DEFAULTS, tick_label_font_size: null, plot_title_font_size: null, keep_full_definition_with_plot_title: false, show_definition: true,
  center_reserved_radius: null, feature_width_circular: null, depth_width_circular: null, gc_content_width_circular: null, gc_content_radius_circular: null,
  gc_skew_width_circular: null, gc_skew_radius_circular: null, circular_track_slots_enabled: false,
  circular_track_slots_schema_version: CIRCULAR_TRACK_SLOT_SCHEMA_VERSION, circular_track_slots_axis_index: null, circular_track_slots: createDefaultCircularTrackSlots(),
  outer_label_x_offset: null, outer_label_y_offset: null, inner_label_x_offset: null, inner_label_y_offset: null
});
/** @typedef {{ candidateLimit: number | null, orthogroupMemberMaxHits: number | null }} LosatpHitLimits `null` is Unbounded */
/** @returns {Record<ProteinBlastpMode, LosatpHitLimits>} */
export const createDefaultLosatpHitLimits = () => ({
  pairwise: { ...FRESH_PROTEIN_COMPARISON_DEFAULTS.hitLimits },
  orthogroup: { ...FRESH_PROTEIN_COMPARISON_DEFAULTS.hitLimits },
  collinear: { ...FRESH_PROTEIN_COMPARISON_DEFAULTS.hitLimits }
});
/**
 * The values of the execution selects (option values are text; Parallel runs' Auto is `undefined`).
 * @typedef {{ executionMode: string, totalThreadBudget: string, threadsPerJob: string, parallelWorkers: string | undefined }} LosatExecution
 */
// How LOSAT runs: one app-level setting for both drawings (Session `ui.losatExecution`).
/** @returns {LosatExecution} */
export const createDefaultLosatExecution = () => ({
  executionMode: 'threaded', totalThreadBudget: 'safe', threadsPerJob: 'auto', parallelWorkers: undefined
});
/**
 * A drawing's LOSAT search settings. A count is `null` when its field is cleared (`v-blank-null`), which reads as its
 * default; a cleared limit is Unbounded.
 * `hitLimitsByMode` keeps each mode's two limits while another mode is shown. A choice holds a value Python reads.
 * @typedef {{ mode: ProteinBlastpMode, maxHits: number | null, candidateLimit: number | null,
 *   hitLimitsByMode: Record<ProteinBlastpMode, Pick<LosatBlastpSettings, keyof LosatpHitLimits>>,
 *   orthogroupMembershipMode: CanonicalProteinComparisonSettings['orthogroupMembershipMode'],
 *   orthogroupMemberMaxHits: number | null, collinearInferOrthogroups: boolean, collinearMinAnchors: number | null,
 *   collinearMaxUnitGap: number | null, collinearMaxDiagonalDrift: number | null, collinearMaxConflictsInMergeGap: number | null,
 *   collinearMaxParalogLinksPerOrthogroup: number | null,
 *   collinearColorMode: CanonicalProteinComparisonSettings['collinearityColorMode'],
 *   collinearUnitMode: CanonicalProteinComparisonSettings['collinearityUnitMode'],
 *   collinearAnchorMode: CanonicalProteinComparisonSettings['collinearityAnchorMode'],
 *   collinearMergeOrientation: CollinearMergeOrientation,
 *   collinearSearchScope: CanonicalProteinComparisonSettings['collinearitySearchScope'] }} LosatBlastpSettings
 * @typedef {{ outfmt: string, blastn: { task: string }, blastp: LosatBlastpSettings }} LosatSettings
 */
/** @returns {LosatSettings} */
export const createDefaultLosat = () => ({
  outfmt: '6',
  blastn: { task: 'megablast' }, blastp: { mode: 'orthogroup', hitLimitsByMode: createDefaultLosatpHitLimits(), maxHits: 5, candidateLimit: FRESH_PROTEIN_COMPARISON_DEFAULTS.hitLimits.candidateLimit, orthogroupMembershipMode: ORTHOGROUP_INFERENCE_VERSION,
    orthogroupMemberMaxHits: FRESH_PROTEIN_COMPARISON_DEFAULTS.hitLimits.orthogroupMemberMaxHits,
    collinearInferOrthogroups: /** @type {boolean} */ (FRESH_PROTEIN_COMPARISON_DEFAULTS.collinearInferOrthogroups), collinearMinAnchors: 1, collinearMaxUnitGap: 0, collinearMaxDiagonalDrift: 0,
    collinearMaxConflictsInMergeGap: 1, collinearMaxParalogLinksPerOrthogroup: 2, collinearColorMode: 'orientation',
    collinearUnitMode: 'auto', collinearAnchorMode: DEFAULT_COLLINEAR_ANCHOR_MODE, collinearMergeOrientation: DEFAULT_MERGE_ORIENTATION, collinearSearchScope: 'adjacent' }
});
// `colorScheme` and `reassignOnReorder` are the Ring colors setting (PD-OI-110): the scheme of the rings without a
// colour of their own, and whether a reorder re-assigns their scheme colours.
export const createDefaultCircularConservation = () => ({ enabled: false, source: 'losat', losat_program: 'blastn',
  subject_gencode: 1, reference: 'auto', labels: '', series: [], ring_width: null, ring_gap: null,
  colorScheme: CONSERVATION_DEFAULT_COLOR_SCHEME, reassignOnReorder: false });
export const CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS = Object.freeze([
  'form', 'adv', 'losat', 'cliOptions', 'colors', 'palette', 'paletteInstantPreviewEnabled', 'rules',
  'qualifierPriorityRules', 'filterMode', 'whitelist', 'blacklistText', 'fileImports', 'losatProgram', 'circularConservation',
  'annotationSets', 'recordDisplayDrafts', 'recordsOff', 'featurePlacementOverrides', 'modeProfiles', 'unmanagedConfigOverrides',
  'linearRecordLayout', 'linearComparisonPlan', 'importedComparisonResolution', 'webEdits'
]);
// The `config` keys the CLI writer of Sessions 40 and 41 derived from its
// options; a Web draft of those Sessions always has `annotationSets`,
// `linearComparisonPlan` and `webEdits` too (OV-269, `_CLI_WRITER_CONFIG_DOMAINS`
// in gbdraw/session_io.py).
const CLI_WRITER_CONFIG_DOMAINS = new Set(['form', 'adv', 'losat', 'cliOptions', 'colors', 'palette', 'rules',
  'qualifierPriorityRules', 'filterMode', 'whitelist', 'blacklistText', 'losatProgram', 'circularConservation']);
/**
 * A Session 40 or 41 the CLI wrote holds a `config` derived from its options,
 * no Web draft: Load reads its tables from the request.
 * @param {Record<string, any>} session
 * @returns {boolean}
 */
export const holdsCliWriterConfig = (session) => [40, 41].includes(session?.version)
  && session.cliInvocation?.generatedBy === 'gbdraw' && isObject(session.config)
  && Object.keys(session.config).every((key) => CLI_WRITER_CONFIG_DOMAINS.has(key));
/** The configuration domains a draft stores, with the compatibility field `colorsAreOverrides`. */
export const STORED_CONFIG_DOMAINS = Object.freeze([...CURRENT_WRITER_ACTIVE_CONFIG_DOMAINS, 'colorsAreOverrides']);
/** The input types of a Session's `ui.cInputType` and `ui.lInputType` (Save refuses another; Load reads it, OV-594). */
export const SESSION_INPUT_TYPES = Object.freeze(/** @type {InputType[]} */ (Object.keys(INPUT_TYPES)));
/**
 * The word after each PNG DPI in the DPI menu; the DPIs are Python's (`SESSION_PNG_DPI_CHOICES`), and so is
 * the one the menu starts at (`DEFAULT_PNG_DPI`, which the command line and Python draw a PNG at too).
 */
export const PNG_DPI_CHOICE_WORDS = Object.freeze(/** @type {Readonly<Record<number, string>>} */ ({
  72: 'Web', 96: 'Screen', 150: 'Draft', 300: 'Print', 600: 'High'
}));
/**
 * The PNG DPI choices the DPI menu lists; Load reads a stored DPI among them (OV-595), as the
 * command line's replay does (`read_session_png_dpi` in gbdraw/session.py).
 */
export const PNG_DPI_CHOICES = Object.freeze(SESSION_PNG_DPI_CHOICES.map((dpi) => ({
  value: String(dpi), label: `${dpi} (${PNG_DPI_CHOICE_WORDS[dpi]})`
})));

/**
 * Stored settings a gbdraw writer wrote and later retired or moved (E-OD03):
 * Load reads them where they moved, or drops them, with no dialog row (TB19,
 * TB32), and Save accepts the `form` and `adv` ones in a Session 40-44 draft.
 * `rich_feature_popup` is read from a Session 40-44 draft only: Session 46
 * keeps the app preference in `ui.richFeaturePopup`.
 */
export const RETIRED_STORED_FIELDS = Object.freeze([
  'form.legend', 'adv.plot_title_position', 'adv.losatProgram', 'adv.rich_feature_popup', 'blastSource', 'priorityRules'
]);
/** @param {string} domain */
const retiredFieldsOf = (domain) => RETIRED_STORED_FIELDS.filter((at) => at.startsWith(`${domain}.`)).map((at) => at.slice(domain.length + 1));
export const CURRENT_WRITER_FORM_FIELDS = Object.freeze([...Object.keys(createDefaultForm()), ...retiredFieldsOf('form')]);
export const CURRENT_WRITER_ADV_FIELDS = Object.freeze([...Object.keys(createDefaultAdv()), ...retiredFieldsOf('adv')]);
const DOMAIN_SHAPES = Object.freeze({ form: 'object', adv: 'object', losat: 'object', cliOptions: 'object', colors: 'object',
  circularConservation: 'object', modeProfiles: 'object', linearRecordLayout: 'object', linearComparisonPlan: 'object', webEdits: 'object',
  importedComparisonResolution: 'object', unmanagedConfigOverrides: 'object',
  recordDisplayDrafts: 'array', recordsOff: 'array', featurePlacementOverrides: 'object', rules: 'array', qualifierPriorityRules: 'array', whitelist: 'array', annotationSets: 'array', palette: 'string', filterMode: 'string', blacklistText: 'string', fileImports: 'object',
  losatProgram: 'string', paletteInstantPreviewEnabled: 'boolean' });
/** The fields of a stored rule row (Save refuses another; Load drops it with a row). */
// A stored annotation set and its rows are the request's (Python's `AnnotationSet` and `RegionAnnotation`, through
// the generated request shape), so a field Python adds is saved without a second list here (OV-2580).
/** @param {string} type @returns {string[]} */
const requestFields = (type) => Object.keys(/** @type {ObjectShape} */ (CANONICAL_REQUEST_DESCRIPTION.types[type]).fields);
export const STORED_ROW_FIELDS = Object.freeze({ rules: ['feat', 'qual', 'val', 'color', 'cap', 'fromFile'],
  qualifierPriorityRules: ['feat', 'order', 'importedFrom', 'handOrder'], whitelist: ['feat', 'qual', 'key', 'importedFrom'],
  annotationSets: requestFields('CanonicalAnnotationSet') });
const ANNOTATION_FIELDS = requestFields('CanonicalRegionAnnotation');
const OBSOLETE_SLOT_KEYS = new Set(['gapAfter', 'gap_after', 'innerRadius', 'inner_radius', 'outerRadius',
  'outer_radius', 'placement', 'spacing', 'strict', 'compress', 'reserve']);
const OBSOLETE_SLOT_PARAM_KEYS = new Set(['side', 'radius', 'width', 'spacing', 'inner_gap_px', 'outer_gap_px', 'strict', 'compress', 'reserve']);
/** @param {Record<string, unknown>} value @param {Set<string>} fields @param {string} path */
const assertFields = (value, fields, path) => {
  const unknown = Object.keys(value).filter((key) => !fields.has(key));
  if (unknown.length) throw new Error(`Current session active configuration contains unknown ${path} field(s): ${unknown.join(', ')}.`);
};
/** @param {Record<string, unknown>} config */
const validateDomainShapes = (config) => {
  for (const [domain, shape] of Object.entries(DOMAIN_SHAPES)) {
    if (!has(config, domain) || (config[domain] === undefined && domain === 'cliOptions')) continue;
    const value = config[domain], valid = shape === 'object' ? isObject(value) : shape === 'array' ? Array.isArray(value) : typeof value === shape;
    if (!valid) throw new Error(`Current session active configuration config.${domain} must be ${shape}.`);
  }
  if (has(config, 'colorsAreOverrides') && typeof config.colorsAreOverrides !== 'boolean')
    throw new Error('Current session compatibility field config.colorsAreOverrides must be a boolean.');
};
/** @param {unknown[]} rows @param {string[]} fields @param {string} path */
const assertRows = (rows, fields, path) => rows.forEach((row, index) => {
  if (!isObject(row)) throw new Error(`Current session active configuration ${path}[${index}] must be an object.`);
  assertFields(row, new Set(fields), `${path}[${index}]`);
});
/**
 * The first loop checks that every annotation set is an object, which the
 * second loop reads.
 * @param {{ annotationSets?: { annotations?: unknown }[], [domain: string]: unknown }} config
 */
const validateCollections = (config) => {
  for (const [domain, fields] of Object.entries(STORED_ROW_FIELDS))
    if (Array.isArray(config[domain])) assertRows(config[domain], fields, `config.${domain}`);
  (config.annotationSets || []).forEach((set, index) => {
    if (!Array.isArray(set.annotations)) throw new Error(`Current session active configuration config.annotationSets[${index}].annotations must be an array.`);
    assertRows(set.annotations, ANNOTATION_FIELDS, `config.annotationSets[${index}].annotations`);
  });
};
/** @param {Record<string, unknown>} slot The obsolete fields of a Circular slot row, at its top level and in its params. */
const obsoleteCircularSlotFields = (slot) => [...OBSOLETE_SLOT_KEYS].filter((key) => has(slot, key))
  .concat(isObject(slot.params) ? [...OBSOLETE_SLOT_PARAM_KEYS].filter((key) => has(slot.params, key)).map((key) => `params.${key}`) : []);
/** @param {unknown} slots */
const obsoleteCircularSlotField = (slots) => {
  if (!Array.isArray(slots)) return null;
  for (const [slotIndex, slot] of slots.entries()) {
    if (!isObject(slot)) continue;
    const [field] = obsoleteCircularSlotFields(slot);
    if (field) return { field: field.replace(/^params\./, ''), slotIndex };
  }
  return null;
};
/** What each mode's slot check reads (`validateTrackSlotBindingInvariants`). */
const TRACK_SLOT_CHECKS = Object.freeze({
  circular: { modeLabel: 'Circular', layoutKind: 'circular', supportedRenderers: CIRCULAR_TRACK_RENDERERS,
    supportedSides: ['inside', 'outside', 'overlay'], anchorlessRenderers: ['ticks', 'spacer'] },
  linear: { modeLabel: 'Linear', layoutKind: 'linear', supportedRenderers: LINEAR_TRACK_RENDERERS,
    supportedSides: ['above', 'below', 'overlay'], anchorlessRenderers: ['spacer'] }
});
/**
 * @param {Record<string, any>} [config] The active configuration.
 * @param {{ depthTrackCount?: number | null }} [options]
 * @returns {void}
 */
export const validateImportedCircularTrackSlots = (config = {}, { depthTrackCount = null } = {}) => {
  const adv = isObject(config) ? config.adv : null;
  if (!isObject(adv) || !has(adv, 'circular_track_slots')) return;
  if (adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION)
    throw new Error(`Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${CIRCULAR_TRACK_SLOT_SCHEMA_VERSION}.`);
  // Any value and either enablement: readers drop only the lossless legacy null.
  const obsolete = obsoleteCircularSlotField(adv.circular_track_slots);
  if (obsolete) throw diagnosticError('TRACK_INVALID', { ...obsolete, reason: 'OBSOLETE_TRACK_FIELD' });
  validateTrackSlotBindingInvariants(adv.circular_track_slots, { ...TRACK_SLOT_CHECKS.circular, depthTrackCount });
};
/**
 * @param {Record<string, any>} [config] The active configuration.
 * @param {{ depthTrackCount?: number | null }} [options]
 * @returns {void}
 */
export const validateImportedLinearTrackSlots = (config = {}, { depthTrackCount = null } = {}) => {
  const adv = isObject(config) ? config.adv : null;
  if (!isObject(adv) || !has(adv, 'linear_track_slots')) return;
  if (![LEGACY_LINEAR_TRACK_SLOT_SCHEMA_VERSION, LINEAR_TRACK_SLOT_SCHEMA_VERSION].includes(adv.linear_track_slots_schema_version))
    throw new Error(`Custom Track Slots use an obsolete schema. Recreate the slots with schema version ${LINEAR_TRACK_SLOT_SCHEMA_VERSION}.`);
  if (!Array.isArray(adv.linear_track_slots)) throw new Error('Custom Track Slots must be an array.');
  validateTrackSlotBindingInvariants(adv.linear_track_slots, { ...TRACK_SLOT_CHECKS.linear, depthTrackCount });
};
/**
 * Load's reading of a drawing's track stack (OV-758, A10), where Save checks
 * it (`validateImported*TrackSlots`, E-OD02). A slot value the check refuses
 * is dropped and its default used (`params` that is no object reads as {});
 * a stack the check still refuses then (a slot that is no object, a missing
 * or repeated id, an unknown renderer, slots that contradict each other)
 * reads as the drawing's default stack, which `reset` reports: the stack's
 * axis index, a place in it, is unset with it.
 * @param {'circular' | 'linear'} mode
 * @param {unknown} slots A list (the settings' merge keeps a stack of another kind out).
 * @returns {FieldReading<unknown> & { reset: boolean }}
 */
export const readImportedTrackSlots = (mode, slots) => {
  const check = TRACK_SLOT_CHECKS[mode];
  const defaultStack = () => (mode === 'circular' ? createDefaultCircularTrackSlots() : createDefaultLinearTrackSlots());
  if (!Array.isArray(slots)) return { value: defaultStack(), lost: false, reset: true };
  /** @type {FieldLeaf[]} */
  const leaves = [];
  const read = slots.map((slot, index) => {
    if (!isObject(slot) || typeof slot.id !== 'string') return slot;
    /** @type {Record<string, unknown>} */
    const row = { ...slot, ...(isObject(slot.params) ? { params: { ...slot.params } } : {}) };
    const faults = trackSlotLeafFaults(row, slot.id.trim(), { ...check, sideSet: new Set(check.supportedSides) }).map(({ field }) => field);
    new Set([...(mode === 'circular' ? obsoleteCircularSlotFields(row) : []), ...faults]).forEach((field) => {
      const [container, key] = field === 'params' ? [null, field] : field.startsWith('params.') ? [/** @type {Record<string, unknown>} */ (row.params), field.slice(7)] : [row, field];
      const stored = container ? container[key] : row.params;
      if (container) delete container[key];
      else row.params = {};
      if (!isUnsetValue(stored)) leaves.push({ at: `[${index}].${field}`, stored, used: container ? undefined : {} });
    });
    return row;
  });
  try {
    validateTrackSlotBindingInvariants(read, check);
  } catch {
    return { value: defaultStack(), lost: slots.length > 0, reset: true };
  }
  return { value: read, lost: leaves.length > 0, ...(leaves.length ? { leaves } : {}), reset: false };
};
/**
 * The containers Load cannot read a drawing's settings without: the mode, the
 * configuration object and its `form` and `adv`. Load reads the values inside
 * them, the track stacks included (`readImportedTrackSlots`), with the reader
 * table (`config-field-readers.js`), which lists what it cannot restore.
 * @param {{ mode: unknown, storedConfig: unknown }} input
 * @returns {Record<string, unknown> & { form: Record<string, unknown>, adv: Record<string, unknown> }} The configuration.
 */
export const requireActiveConfigContainers = ({ mode, storedConfig: config }) => {
  if (!(/** @type {readonly unknown[]} */ (['circular', 'linear'])).includes(mode)) throw new Error(`Current session active configuration has unsupported mode: ${String(mode)}.`);
  if (!isObject(config)) throw new Error('Current session is missing its active Web configuration.');
  assertSafeObjectKeys(config);
  if (!isObject(config.form) || !isObject(config.adv)) throw new Error('Current session is missing its active form or advanced settings.');
  // The cast: the two tests above.
  return /** @type {Record<string, unknown> & { form: Record<string, unknown>, adv: Record<string, unknown> }} */ (config);
};
/**
 * Throws when `storedConfig` is not a current-writer active configuration
 * (Save and Gallery publication): a Session 46 slice of `mode`, or
 * (`scopedDrafts`) the flat draft of a Session 40-44, whose record display and
 * placement rows name their mode.
 * @param {{ mode: unknown, storedConfig: unknown, scopedDrafts?: boolean }} input
 * @returns {void}
 */
export const validateCurrentWriterActiveConfig = ({ mode, storedConfig, scopedDrafts = false }) => {
  const config = requireActiveConfigContainers({ mode, storedConfig });
  validateImportedCircularTrackSlots(config); validateImportedLinearTrackSlots(config);
  const domains = new Set(STORED_CONFIG_DOMAINS);
  const unknownDomains = Object.keys(config).filter((domain) => !domains.has(domain));
  if (unknownDomains.length)
    throw new Error(`Current session active configuration contains unknown domain(s): ${unknownDomains.join(', ')}.`);
  // A list of Linear card uids or Circular source selectors `#N`, once each, as
  // Python reads it; that a Linear uid names a bound card is checked with the
  // bindings (session-authority.js).
  if (has(config, 'recordsOff') && (!Array.isArray(config.recordsOff)
    || new Set(config.recordsOff).size !== config.recordsOff.length
    || !config.recordsOff.every((/** @type {unknown} */ key) => isRecordDrawKey(/** @type {'circular' | 'linear'} */ (mode), key))))
    throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  validateDomainShapes(config); validateCollections(config); requireCurrentConfigFields(config);
  if (has(config, 'recordDisplayDrafts')) validateRecordDisplayDrafts(config.recordDisplayDrafts, { scoped: scopedDrafts });
  if (has(config, 'featurePlacementOverrides')) {
    if (scopedDrafts) validateScopedFeaturePlacements(config.featurePlacementOverrides);
    else canonicalFeaturePlacements(config.featurePlacementOverrides, /** @type {'circular' | 'linear'} */ (mode)); // the first check
  }
  if (has(config.adv, 'feature_overlap_tolerance_bp') && (!Number.isSafeInteger(config.adv.feature_overlap_tolerance_bp)
    || /** @type {number} */ (config.adv.feature_overlap_tolerance_bp) < 0)) throw new Error('Feature overlap tolerance must be a non-negative integer.'); // a safe integer here
  assertFields(config.form, new Set(CURRENT_WRITER_FORM_FIELDS), 'config.form'); assertFields(config.adv, new Set(CURRENT_WRITER_ADV_FIELDS), 'config.adv');
  if (has(config.form, 'keep_definition_left_aligned') && typeof config.form.keep_definition_left_aligned !== 'boolean')
    throw new Error('Current session active configuration config.form.keep_definition_left_aligned must be a boolean.');
  for (const [path, value] of [['config.adv.losatProgram', config.adv.losatProgram], ['config.losatProgram', config.losatProgram]]) {
    if (value !== undefined && !(/** @type {readonly unknown[]} */ (['blastn', 'tblastx', 'blastp'])).includes(value))
      throw new Error(`Current session active configuration ${path} is invalid.`);
  }
  // `losat` is an object when present (`validateDomainShapes`); the test proves its `blastp` an object.
  if (isObject(/** @type {{ blastp?: unknown } | undefined} */ (config.losat)?.blastp)) {
    const blastp = /** @type {{ blastp: Record<string, unknown> }} */ (config.losat).blastp;
    if (has(blastp, 'mode')) requireCurrentProteinBlastpMode(blastp.mode);
    if (has(blastp, 'hitLimitsByMode')) {
      if (!isObject(blastp.hitLimitsByMode)) throw new Error('LOSATP mode hit limits must be an object.');
      for (const [mode, limits] of Object.entries(blastp.hitLimitsByMode)) {
        requireCurrentProteinBlastpMode(mode);
        if (!isObject(limits)) throw new Error('LOSATP mode hit limits must be an object.');
        assertFields(limits, new Set(['candidateLimit', 'orthogroupMemberMaxHits']), `losat.blastp.hitLimitsByMode.${mode}`);
        requireCurrentProteinBlastpCandidateLimit(limits.candidateLimit);
        requireCurrentOrthogroupMemberMaxHits(limits.orthogroupMemberMaxHits);
      }
    }
    if (has(blastp, 'candidateLimit')) {
      requireCurrentProteinBlastpCandidateLimit(blastp.candidateLimit);
    }
    requireCurrentProteinBlastpMaxHits(blastp.maxHits);
    requireCurrentOrthogroupMembershipMode(blastp.orthogroupMembershipMode);
    requireCurrentOrthogroupMemberMaxHits(blastp.orthogroupMemberMaxHits);
    requireCurrentCollinearMinAnchors(blastp.collinearMinAnchors);
    requireCurrentCollinearMaxUnitGap(blastp.collinearMaxUnitGap);
    requireCurrentCollinearMaxDiagonalDrift(blastp.collinearMaxDiagonalDrift);
    requireCurrentCollinearMaxConflicts(blastp.collinearMaxConflictsInMergeGap);
    requireCurrentCollinearMaxParalogLinks(
      blastp.collinearMaxParalogLinksPerOrthogroup
    );
    requireCurrentCollinearUnitMode(blastp.collinearUnitMode);
    requireCurrentCollinearAnchorMode(blastp.collinearAnchorMode);
    requireCurrentCollinearMergeOrientation(blastp.collinearMergeOrientation);
    requireCurrentCollinearColorMode(blastp.collinearColorMode);
    requireCurrentCollinearSearchScope(blastp.collinearSearchScope);
    requireCurrentCollinearInferOrthogroups(blastp.collinearInferOrthogroups);
  }
  if (has(config, 'filterMode') && !(/** @type {readonly unknown[]} */ (['None', 'Whitelist', 'Blacklist'])).includes(config.filterMode)) throw new Error('Current session active configuration config.filterMode is invalid.');
  // `palette` is a string when present (`validateDomainShapes`).
  if (has(config, 'palette') && !/** @type {string} */ (config.palette).trim()) throw new Error('Current session active configuration config.palette cannot be empty.');
  if (isObject(config.importedComparisonResolution)) {
    assertFields(
      config.importedComparisonResolution,
      new Set(['action']),
      'config.importedComparisonResolution'
    );
    if (
      config.importedComparisonResolution.action !== null
      && !(/** @type {readonly unknown[]} */ (['INHERIT', 'REPLACE', 'CLEAR'])).includes(config.importedComparisonResolution.action)
    ) {
      throw new Error('Current session active configuration config.importedComparisonResolution.action is invalid.');
    }
  }
};

/** @param {{ scale_font_size?: unknown, ruler_label_font_size?: unknown }} [adv] */
const linearTypographyValuesMatch = (adv = {}) => (
  Object.is(adv.scale_font_size, adv.ruler_label_font_size)
);

/**
 * A stored Linear typography link: linked when omitted, as stored when a
 * boolean, and unlinked for anything else.
 * @param {unknown} value
 * @returns {{ value: boolean, lost: boolean }}
 */
export const readLinearTypographyLinked = (value) => (value === undefined || value === null ? { value: true, lost: false }
  : typeof value === 'boolean' ? { value, lost: false } : { value: false, lost: true });
/**
 * @param {{
 *   adv: Record<string, any>,
 *   linked: any,
 *   ui?: { readonly linearTypographyLinked?: unknown }
 * }} options
 */
export const reconcileImportedLinearTypographyLink = ({ adv, linked, ui = {} }) => {
  if (!linked || typeof linked !== 'object' || !('value' in linked)) return false;
  // Omission takes the fresh linked default; unequal values still open unlinked.
  linked.value = readLinearTypographyLinked(ui.linearTypographyLinked).value && linearTypographyValuesMatch(adv);
  return linked.value;
};

// Artifact metadata admission is shared by Session, History and Align/Reset.
// The binding deliberately excludes current orientation, translations, labels and
// row order: manual Reverse, style and stable reorder do not rewrite history.
/** @param {unknown} value @returns {unknown} */
const sortedValue = (value) => Array.isArray(value)
  ? value.map(sortedValue)
  : isObject(value) ? Object.fromEntries(Object.keys(value).sort()
    .map((key) => [key, sortedValue(value[key])])) : value;

// The parts of a canonical request that the Alignment reset receipt reads. The
// writer (`canonicalSimilarityAlignment` in session-request.js) owns the plan.
/**
 * @typedef {{
 *   records: { recordKey: string, status?: string }[],
 *   reference: { recordKey: string }
 * }} AlignmentPlanView
 */
/**
 * @typedef {{
 *   recordKey: string, source: Record<string, string>, selector?: unknown,
 *   windows: { start: number, end: number, reverseComplement: boolean }[],
 *   presentation: { reverseComplement: boolean }, display?: unknown
 * }} AlignmentResetRecord
 */
/**
 * @typedef {{
 *   mode: string, records: unknown[],
 *   layout: { similarityAlignment?: AlignmentPlanView | null, recordTranslations?: { recordKey: string, x?: number }[], [field: string]: unknown }
 * }} AlignmentResetRequest
 */
/**
 * @typedef {{
 *   renderRequest?: AlignmentResetRequest | null,
 *   resources?: Record<string, RawSessionResourceDescriptor> | null
 * }} AlignmentResetCanonical
 */
/**
 * The canonical session an Align started from (`before` of a reset), as `projectCommittedSimilarityAlignment`
 * returns it: the generated request types with a Linear layout. Its plan is null: Align materializes its
 * baseline without a plan. Not for a Session's request, which the `AlignmentResetCanonical` view above reads.
 * @typedef {{
 *   renderRequest: Omit<CanonicalLinearRenderRequest, 'layout'> & { layout: CanonicalLinearLayout },
 *   resources?: Record<string, RawSessionResourceDescriptor> | null
 * }} AlignmentResetBaselineCanonical
 */
/**
 * The canonical session an Align produced (`after` of a reset): the same, with the plan present.
 * @typedef {AlignmentResetBaselineCanonical & {
 *   renderRequest: { layout: { similarityAlignment: CanonicalSimilarityAlignmentPlan } }
 * }} AlignmentResetPlannedCanonical
 */
/**
 * What a reset writes beside the plan: the binding of the plan and the sources,
 * the direction changes, and the reference record's shift (or null).
 * @typedef {{
 *   binding: string,
 *   directions: { recordKey: string, before: boolean, after: boolean }[],
 *   referenceDeltaX: { recordKey: string, deltaX: number } | null
 * }} AlignmentResetReceipt
 */

/**
 * The one shape check of a stored alignment reset receipt (OV-115, OV-130):
 * whether `receipt` is absent or is the receipt of `request`'s plan.
 * @param {unknown} receipt The stored receipt, unvalidated.
 * @param {AlignmentResetRequest | null | undefined} request
 * @returns {boolean}
 */
export const alignmentResetReceiptShapeHolds = (receipt, request) => {
  if (receipt === null || receipt === undefined) return true;
  if (!isObject(receipt)
    || Object.keys(receipt).sort().join(',') !== 'binding,directions,referenceDeltaX'
    || typeof receipt.binding !== 'string' || !/^[0-9a-f]{64}$/.test(receipt.binding)
    || !Array.isArray(receipt.directions)
    || request?.mode !== 'linear' || !request.layout?.similarityAlignment) return false;
  const plan = request.layout.similarityAlignment;
  const eligible = new Set(plan.records.filter(({ status }) => status !== 'skipped')
    .map(({ recordKey }) => recordKey));
  /** @type {Set<unknown>} */
  const keys = new Set();
  const directionsHold = receipt.directions.every((delta) => {
    if (!isObject(delta) || Object.keys(delta).sort().join(',') !== 'after,before,recordKey'
      || !eligible.has(delta.recordKey) || keys.has(delta.recordKey)
      || typeof delta.before !== 'boolean' || typeof delta.after !== 'boolean'
      || delta.before === delta.after) return false;
    keys.add(delta.recordKey);
    return true;
  });
  const delta = receipt.referenceDeltaX;
  return directionsHold && (delta === null || (isObject(delta)
    && Object.keys(delta).sort().join(',') === 'deltaX,recordKey'
    && delta.recordKey === plan.reference.recordKey
    && typeof delta.deltaX === 'number' && Number.isFinite(delta.deltaX)
    && delta.deltaX !== 0));
};

/**
 * @param {unknown} receipt The stored receipt, unvalidated: this function checks it.
 * @param {AlignmentResetRequest | null | undefined} request
 * @returns {AlignmentResetReceipt | null}
 */
export const validateAlignmentResetReceiptShape = (receipt, request) => {
  // A malformed or stale receipt (OV-115, OV-130) is one user-facing failure.
  if (!alignmentResetReceiptShapeHolds(receipt, request)) throw diagnosticError('ALIGNMENT_RESET_EVIDENCE');
  // Absent, or every field checked above.
  return /** @type {AlignmentResetReceipt | null} */ (receipt ?? null);
};

/**
 * @param {AlignmentResetCanonical} canonical
 * @param {(record: AlignmentResetRecord) => unknown} [selectorOf] Each record's selector as the binding reads it.
 */
const resetBinding = async (canonical, selectorOf = (record) => record.selector) => {
  const request = canonical?.renderRequest;
  if (request?.mode !== 'linear' || !request.layout?.similarityAlignment) {
    throw new Error('Alignment reset receipt requires an active canonical plan.');
  }
  const table = adoptCurrentSessionResources(canonical.resources);
  /** @type {Map<string, Promise<string>>} */
  const fingerprints = new Map();
  /** @param {Record<string, string>} source */
  const sourceIdentity = async (source) => {
    /** @type {Record<string, string | undefined>} */
    const identity = {};
    for (const [key, value] of Object.entries(source)) {
      if (key === 'kind') identity[key] = value;
      else {
        if (!fingerprints.has(value)) fingerprints.set(value, sha256Hex(
          await readSessionResourceBytes(createSessionResourceFileView(table, value))
        ));
        identity[key] = await fingerprints.get(value);
      }
    }
    return identity;
  };
  // One window binds as the region it was before windows (request schema 10), so stored receipts still match
  // (Python's `region_binding` in `session_io.py`).
  /** @param {AlignmentResetRecord['windows']} windows */
  const regionBinding = (windows) => {
    const spans = windows.map(({ start, end }) => ({ start, end }));
    return spans.length === 1 ? spans[0] : spans.length ? spans : null;
  };
  // The mode check above admits a Linear request only; the request writer gives each of its records a `recordKey` and a `source`.
  const records = await Promise.all(/** @type {AlignmentResetRecord[]} */ (request.records).map(async (record) => ({
    recordKey: record.recordKey, source: await sourceIdentity(record.source),
    selector: selectorOf(record),
    region: regionBinding(record.windows),
    display: record.display
  })));
  const plan = { ...request.layout.similarityAlignment,
    records: [...request.layout.similarityAlignment.records].sort(
      (a, b) => a.recordKey < b.recordKey ? -1 : a.recordKey > b.recordKey ? 1 : 0) };
  records.sort((a, b) => a.recordKey < b.recordKey ? -1 : a.recordKey > b.recordKey ? 1 : 0);
  return sha256Hex(textToBytes(JSON.stringify(sortedValue({ plan, records }))));
};

/**
 * The one reader of a stored alignment reset receipt: the receipt when it is
 * the receipt of `canonical`'s plan and sources, else none, with `stale` set.
 * Load reads a stale receipt as none, with a row (A19); History, Align and
 * Save refuse it (`validateSimilarityAlignmentResetReceipt`).
 * @param {unknown} receipt The stored receipt, unvalidated.
 * @param {AlignmentResetCanonical | null | undefined} canonical
 * @returns {Promise<{ value: AlignmentResetReceipt | null, stale: boolean }>}
 */
export const readSimilarityAlignmentResetReceipt = async (receipt, canonical) => {
  if (!alignmentResetReceiptShapeHolds(receipt, canonical?.renderRequest)) return { value: null, stale: true };
  const value = /** @type {AlignmentResetReceipt | null} */ (receipt ?? null); // checked above
  // A receipt that passed the shape check needs `request.mode === 'linear'`, so `canonical` is set.
  if (value && value.binding !== await resetBinding(/** @type {AlignmentResetCanonical} */ (canonical))) return { value: null, stale: true };
  return { value, stale: false };
};

/**
 * The record selectors of a request the file stored before windows (its records hold a `region`), by record key;
 * null for any other request. Values are not validated.
 * @param {unknown} request The stored request, unvalidated.
 * @returns {Map<unknown, unknown> | null}
 */
const storedRegionRecordSelectors = (request) => {
  const records = isObject(request) && Array.isArray(request.records) ? request.records.filter(isObject) : [];
  if (!records.some((record) => Object.hasOwn(record, 'region'))) return null;
  return new Map(records.map((record) => [record.recordKey, record.selector ?? null]));
};

/**
 * Load: a stored receipt read against the request Load promoted (`readSimilarityAlignmentResetReceipt`). A request
 * the file stored before windows bound each record's selector where the file stored it: the v0.14.5 Web wrote a
 * cropped Linear card's selector in its region and bound the record's own, `null`, and promotion moves the selector
 * onto the record (`promoteRecordRegionsToWindows`, OV-2772). Such a receipt is read with the stored selectors, as
 * Python's verifier reads the stored request (`_validate_alignment_reset_receipt`), and is rebound to the promoted
 * request, which History, Align and Save read.
 * @param {unknown} receipt The stored receipt, unvalidated.
 * @param {AlignmentResetCanonical | null | undefined} canonical The promoted request and the resources.
 * @param {unknown} storedRequest The request as the file stored it, unvalidated.
 * @returns {Promise<{ value: AlignmentResetReceipt | null, stale: boolean }>}
 */
export const readStoredSimilarityAlignmentResetReceipt = async (receipt, canonical, storedRequest) => {
  const read = await readSimilarityAlignmentResetReceipt(receipt, canonical);
  const storedSelectors = storedRegionRecordSelectors(storedRequest);
  if (!read.stale || storedSelectors === null || !receipt
    || !alignmentResetReceiptShapeHolds(receipt, canonical?.renderRequest)) return read;
  // The shape check passed: a receipt of a Linear request with a plan, so `canonical` is set.
  const value = /** @type {AlignmentResetReceipt} */ (receipt);
  const current = /** @type {AlignmentResetCanonical} */ (canonical);
  const asStored = await resetBinding(current, (record) => (
    storedSelectors.has(record.recordKey) ? storedSelectors.get(record.recordKey) : record.selector));
  if (value.binding !== asStored) return read;
  return { value: { ...value, binding: await resetBinding(current) }, stale: false };
};

/**
 * @param {unknown} receipt The stored receipt, unvalidated.
 * @param {AlignmentResetCanonical | null | undefined} canonical
 */
export const validateSimilarityAlignmentResetReceipt = async (receipt, canonical) => {
  const read = await readSimilarityAlignmentResetReceipt(receipt, canonical);
  if (read.stale) throw diagnosticError('ALIGNMENT_RESET_EVIDENCE');
  return read.value;
};

/**
 * @param {{ before: AlignmentResetBaselineCanonical, after: AlignmentResetPlannedCanonical }} sessions
 */
export const buildSimilarityAlignmentResetReceipt = async ({ before, after }) => {
  const previous = new Map(before.renderRequest.records.map((record) => [record.recordKey, record]));
  const directions = after.renderRequest.records.flatMap((record) => {
    const old = previous.get(record.recordKey);
    if (!old) throw new Error('Alignment reset receipt record binding changed.');
    const beforeReverse = canonicalRecordReverseComplement(old);
    const afterReverse = canonicalRecordReverseComplement(record);
    return beforeReverse === afterReverse ? []
      : [{ recordKey: record.recordKey, before: beforeReverse, after: afterReverse }];
  });
  const referenceKey = after.renderRequest.layout.similarityAlignment.reference.recordKey;
  /** @param {AlignmentResetBaselineCanonical} canonical */
  const x = (canonical) => canonical.renderRequest.layout.recordTranslations
    .find(({ recordKey }) => recordKey === referenceKey)?.x ?? 0;
  const deltaX = x(after) - x(before);
  const receipt = { binding: await resetBinding(after), directions,
    referenceDeltaX: deltaX === 0 ? null : { recordKey: referenceKey, deltaX } };
  return validateSimilarityAlignmentResetReceipt(receipt, after);
};
