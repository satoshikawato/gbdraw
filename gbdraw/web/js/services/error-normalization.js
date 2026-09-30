const OPERATIONS = new Set(['unknown', 'generate', 'align', 'feature-extraction', 'export-svg', 'export-png', 'export-pdf', 'evaluateRules', 'readPdfFont', 'buildProteinLosatCacheKeys', 'convertLosatNucleotideToDisplayTsv',
  'convertLosatpPairsToGenomicPayload', 'extractCdsProteinFasta', 'extractFirstFasta', 'generateLegendEntrySvg',
  'hydrateProteinLosatTsv', 'listGffFastaRecords', 'listSequenceRecords', 'measureLegendText',
  'promoteLegacyLosatpCache', 'regenerateDefinitionSvgs', 'resolveLegacyProteinReferences',
  'resolveSimilarityAlignment', 'validateConfigOverrides']);
const STAGES = new Set(['unknown', 'initialization', 'resource-staging', 'request-validation',
  'helper', 'rule-validation', 'render', 'result-admission', 'cleanup', 'export-capture', 'export-conversion', 'font-validation', 'transport', 'read']);
const FIELDS = new Set(`legend title scale decorations pattern record_selector region start end sourceStart sourceEnd recordLength
recordIndex queryIndex subjectIndex depth min_depth max_depth window step tick font_size plot_title_font_size height large_tick_interval small_tick_interval tick_font_size
inner_gap_px outer_gap_px radius width spacing arrow_head_length_ratio arrow_shaft_width_ratio keep_definition_left_aligned color action feature_type qualifier value record_id label_text
config configOverrides records anchors schema recordKey groupId direction sourceStrand role blast files
input comparison protein_blastp_max_hits orthogroup_member_max_hits bitscore evalue identity
alignment_length collinear_min_anchors collinear_max_gene_gap collinear_block_merge_gap
collinear_singleton_merge_gap collinear_max_diagonal_drift collinear_gap_penalty
collinear_nearby_duplicate_window collinear_constant_anchor_score collinear_infer_orthogroups
collinear_min_score collinear_min_block_span record_gap_px record_axis_height depth_window depth_step
depth_min depth_max min_gc max_gc gc_tick_interval gc_axis_font_size depth_tick_interval
depth_axis_font_size protein_blastp_mode protein_blastp_candidate_limit collinear_search_scope collinear_unit_mode collinear_anchor_mode collinear_merge_orientation collinear_color_mode orthogroup_membership_mode collinear_max_unit_gap collinear_max_conflicts collinear_max_paralog_links_per_orthogroup circular_multi_record_size_mode linear_track_layout linear_label_placement set_id anchor_slot side renderer lane_gap_px padding_px cover_anchor overflow layer z axis match_height source fasta gff annotations featurePlacements output_prefix`.split(/\s+/));
const REASONS = Object.freeze({
  DECORATION_IDENTITY: 'Source, region, mode, grouping or record identity changed, is unknown, or is ambiguous.',
  DECORATION_TARGET: 'The decoration target is missing or ambiguous.',
  DECORATION_METADATA: 'Composition metadata is invalid or unavailable.',
  ANNOTATION_SET: 'Select an existing annotation set.', ANNOTATION_MARKS: 'Use supported annotation mark filters.',
  OVERFLOW: 'Choose error, compress, or clip.', LAYER: 'Choose foreground or underlay.',
  OVERLAY_ANCHOR: 'Select an eligible enabled anchor for an Overlay annotation track.',
  OVERLAY_SIDE: 'Use Overlay when an anchor is selected.', DRAW_ORDER: 'Place underlay annotations below the anchor and foreground annotations above it.',
  TRACK_RENDERER: 'Choose a renderer supported by the current diagram mode.', TRACK_SIDE: 'Choose a side supported by the current diagram mode and renderer.',
  TRACK_PARAMS: 'Use the parameters supported by the selected renderer.', TRACK_AXIS: 'Select one eligible enabled feature track as the axis.',
  TRACK_LANE: 'Use a supported feature lane and matching side.', CONSERVATION_SOURCE: 'Use the supported conservation workflow and select a comparison source.',
  TICK_PARAMETER: 'Use large_tick_interval instead of tick_interval.', FEATURES_COUNT: 'Use the supported number of feature tracks for this layout.',

  CIRCULAR_MULTI_RECORD_SIZE_MODE: 'Choose auto, linear, or equal.',
  LINEAR_TRACK_LAYOUT: 'Choose above, middle, or below.',
  LINEAR_LABEL_PLACEMENT: 'Choose auto or above_feature.',
  PROTEIN_BLASTP_MODE: 'Choose pairwise, orthogroup, or collinear.',
  COLLINEAR_SEARCH_SCOPE: 'Choose adjacent or all.',
  COLLINEAR_UNIT_MODE: 'Choose auto, cds, or locus.',
  COLLINEAR_ANCHOR_MODE: 'Choose all, one_to_one, or rbh.',
  COLLINEAR_MERGE_ORIENTATION: 'Choose strand, order, or either.',
  COLLINEAR_COLOR_MODE: 'Choose average_identity, orientation, or orientation_identity.',
  ORTHOGROUP_MEMBERSHIP_MODE: 'Choose anchor_core_v1.',

  UNKNOWN_CONFIG_PATH: 'Remove unknown configuration overrides or use a supported setting.',
  RESOURCE_SIZE: 'Session resource byte sizes must be non-negative safe integers.',
  TRACK_SCHEMA: 'Recreate Custom Track Slots with the current schema.',
  SESSION_FIELDS: 'Remove unsupported top-level Session fields.',
  RECORDS_REQUIRED: 'Session render requests require at least one record.',
  POSITIVE_UNIT_INTERVAL: 'Use a finite number greater than zero and at most one.',
  PIXEL_POSITIVE: 'Use a finite number of pixels greater than zero (px optional).',
  PIXEL_NONNEGATIVE: 'Use a finite number of pixels of zero or greater (px optional).',
  POSITIVE_SCALAR: 'Use a positive finite px or factor scalar.',
  CIRCULAR_GAPS: 'Use inner_gap_px and outer_gap_px for physical gaps.',
  SEPARATE_LINEAR_ROWS: 'Turn Normalize Record Lengths off or assign each record to a separate Linear row.',
  GFF_FASTA_MATCH: 'Ensure every GFF3 record has a matching FASTA entry.',
  UNTERMINATED_SET: 'Close the character set with ].',
  UNTERMINATED_GROUP: 'Close the group with ).', UNBALANCED_GROUP: 'Check the matching parentheses.',
  NOTHING_TO_REPEAT: 'Place a value before the repetition operator.', MULTIPLE_REPEAT: 'Remove the repeated repetition operator.',
  FLAGS_POSITION: 'Place global flags at the start of the pattern.', LOOKBEHIND_WIDTH: 'Use a fixed-width lookbehind.',
  INVALID_ESCAPE: 'Check the escape sequence.', UNKNOWN_EXTENSION: 'Check the Python group syntax.',
  CHARACTER_RANGE: 'Check the character range.', GROUP_REFERENCE: 'Check the group reference.',
  GROUP_NAME: 'Check the group name.', SYNTAX_ERROR: 'Correct the Python regular expression.',
  ADJACENT_ALL: 'Choose adjacent or all.', BLASTP_MODE: 'Choose pairwise, orthogroup, or collinear.',
  UNIQUE_IDS: 'Use distinct record or protein identifiers.', MATCH_IDS: 'Use matching protein FASTA and metadata identifiers.',
  JSON_FORMAT: 'Use valid JSON.', VISIBILITY_ACTION: 'Use show, off, or exclude_matching; accepted aliases are on, hide, false, and 0.',
  TARGET_RECORD: 'Choose an available target record.', TARGET_MODE: 'Clear the target record or enable Multi-record canvas.',
  DEPTH_SERIES: 'Select an existing Depth TSV or remove the slot.',
  BOTH_ENDPOINTS: 'Supply both region endpoints or leave both empty.', SPECIFIC_COLUMNS: 'Supply four or five tab-separated columns.',
  COLOR: 'Use none, a supported named color, or a hex color with 3 or 6 digits.',
  BOOLEAN: 'Use true or false.', INTEGER: 'Use an integer.', NONNEGATIVE: 'Use a finite value of zero or greater.',
  FINITE: 'Use a finite number.', POSITIVE_INTEGER: 'Use an integer greater than zero.',
  POSITIVE: 'Use a finite value greater than zero.', POSITIVE_OR_AUTO: 'Use Auto or a finite value greater than zero.',
  NONNEGATIVE_INTEGER: 'Use an integer of zero or greater.', PERCENT: 'Use a finite value between 0 and 100.',
  ARRAY: 'Use a list.', OBJECT: 'Use an object.', FIELDS: 'Check the required fields.', REQUIRED: 'Supply the required value.',
  STRICT_ORDER: 'The start must be less than the end.',
  ORDER: 'The start or minimum must not exceed the end or maximum.', RECORD_BOUNDS: 'Keep the region within the record length.',
  STRAND: 'Use -1, 1, or no strand.', CROP_START_CONFLICT: 'Choose a crop or an explicit display start.',
  REFERENCE_REQUIRED: 'Supply the depth reference column.', REFERENCE_MISMATCH: 'Match depth references to the selected record.',
  WORKSPACE: 'Retry the operation.', OUT_OF_RANGE: 'Choose a record within the loaded range.',
  NO_MATCH: 'Choose an available record.', AMBIGUOUS: 'Use #index to distinguish records with the same ID.',
  SELECTOR_FORMAT: 'Use #<number> or a record ID.', SELECT_ONE: 'Select exactly one record.',
  REGION_FORMAT: 'Use record_id:start-end[:rc] or #index:start-end[:rc].',
  CANNOT_FIT: 'Move the track, reduce widths, disable conflicting labels, or place it outside.',
  THREE_COLUMNS: 'Supply at least three tab-separated columns.', READ: 'Replace or reselect the input.',
  FORMAT: 'Use GenBank or the required GFF3 and FASTA inputs.', NO_PROTEINS: 'Choose input containing CDS proteins.',
  EMPTY_ENDPOINT: 'Check the comparison endpoints.', INDEX_ALIGNMENT: 'Check the comparison endpoints.',
  SOURCE_INDEX: 'Check the comparison endpoints.', SOURCE_VIEW_CONFLICT: 'Check the comparison inputs and display transforms.'
});
const DEFINITIONS = Object.freeze({
  SESSION_SIZE_LIMIT: ['The Session exceeds the browser size limit. Use a smaller Session file.', ['select-input', 'retry']],
  SESSION_BROWSER_UNSUPPORTED: ['This browser cannot complete this Session operation. Use a browser with the required gzip and Worker support.', ['retry']],
  SESSION_IMPORT_UNAVAILABLE: ['This browser cannot load Sessions because Session import Workers are unavailable. Use a browser with Worker support.', ['retry']],
  DECORATION_CONTINUITY: ['Could not preserve decoration placement. Reset the affected position or use Reset Layout on the previous Result, or restore matching settings, then Generate again.', ['edit-input', 'retry', 'save-session']],
  GENERATION_BUSY: ['A diagram generation request is already running. Wait for it to finish before retrying.', ['retry', 'save-session']],
  UNKNOWN: ['The operation failed without recognized diagnostic information. Retry; if it continues, save a Session for investigation.', ['retry', 'save-session']],
  VALIDATION_UNCLASSIFIED: ['Input validation failed. Review the inputs before retrying.', ['review-input', 'retry']],
  INPUT_INVALID: ['An input value is invalid.', ['edit-input', 'retry']],
  INPUT_REQUIRED: ['Supply GenBank input or matching GFF3 and FASTA inputs.', ['select-input', 'retry']],
  FASTA_REQUIRED: ['Supply a matching FASTA input for each GFF3 input.', ['select-input', 'retry']],
  INPUT_UNREADABLE: ['An input could not be read. Replace or reselect it.', ['select-input', 'retry']],
  NO_RECORDS: ['No records were found. Choose input containing records.', ['select-input', 'retry']],
  RECORD_SELECTION: ['The record selection is invalid.', ['select-record', 'retry']],
  REGION_INVALID: ['The region is invalid.', ['edit-region', 'retry']],
  DEPTH_INVALID: ['The depth input or settings are invalid.', ['edit-depth', 'disable-track', 'retry']],
  TABLE_INVALID: ['The table is invalid.', ['edit-table', 'retry']],
  COMPARISON_INPUT: ['The comparison input is invalid. Supply comparison FASTA or BLAST outfmt 6/7 as required.', ['edit-comparison', 'retry']],
  COMPARISON_IDENTITY: ['Comparison endpoints disagree with the displayed features. Review the comparison inputs and display transforms; save a Session if it continues.', ['edit-comparison', 'retry', 'save-session']],
  ANNOTATION_TARGET: ['The region annotation target is invalid.', ['edit-annotation', 'retry']],
  TRACK_INVALID: ['The track settings are invalid.', ['edit-track', 'retry']],
  TRACK_LAYOUT: ['A circular track does not fit inside.', ['edit-track', 'retry']],
  REGEX_SYNTAX: ['The Python regular expression is invalid.', ['edit-pattern', 'retry']],
  RESOURCE_INVALID: ['Input resource preparation failed. Reselect the input and retry.', ['select-input', 'retry']],
  HELPER_PROTOCOL: ['The helper request is invalid. Retry the operation.', ['retry', 'save-session']],
  RUNTIME_INCOMPATIBLE: ['The diagram engine is incompatible with this Web app. Reload the page and retry; contact the site administrator if it continues.', ['reload', 'retry']],
  WORKER_INIT: ['The diagram runtime could not start. Save a Session and reload before retrying.', ['save-session', 'reload', 'retry']],
  FEATURE_METADATA: ['The diagram engine returned incompatible feature metadata. Reload the page and Generate again.', ['reload', 'retry']],
  RESULT_INVALID: ['The generated result could not be accepted. Retry the operation.', ['retry', 'save-session']],
  EXPORT_INPUT: ['The current diagram could not be prepared for export. Generate a diagram and retry.', ['generate', 'retry']],
  EXPORT_DIMENSIONS: ['The diagram has no usable export dimensions. Generate again or use SVG.', ['generate', 'use-svg']],
  PNG_DPI: ['The PNG DPI is invalid. Use a positive finite DPI and retry.', ['edit-dpi', 'retry']],
  EXPORT_CONVERSION: ['The browser could not convert the diagram for export. Retry or use SVG.', ['retry', 'use-svg']],
  PDF_LIBRARY: ['The PDF libraries could not start. Reload and retry or use SVG.', ['reload', 'retry', 'use-svg']],
  PDF_GLYPH: ['PDF fonts cannot represent this text. Use SVG to retain this text.', ['use-svg']],
  CLEANUP_FAILED: ['Temporary resource cleanup failed. Save a Session and reload before retrying.', ['save-session', 'reload', 'retry']]
});

const FIELD_LABELS = Object.freeze({ protein_blastp_max_hits: 'Protein BLASTP Pairwise max hits',
  arrow_head_length_ratio: 'Arrow head length ratio', arrow_shaft_width_ratio: 'Arrow shaft width ratio',
  keep_definition_left_aligned: 'Lock Definition Column' });
const CODES = new Set(Object.keys(DEFINITIONS));
// Exact native JS validation messages with no document interpolation.
const NATIVE_VALIDATIONS = new Map([
  ['The diagram engine returned incompatible feature metadata. Reload the page and Generate again.', { code: 'FEATURE_METADATA', stage: 'result-admission' }],
  ['A File-like object with arrayBuffer() or text() is required.', { code: 'RESOURCE_INVALID', stage: 'resource-staging' }],
  ['A File-like object with arrayBuffer() is required.', { code: 'RESOURCE_INVALID', stage: 'resource-staging' }],
  ['Session file is too large.', { code: 'SESSION_SIZE_LIMIT', stage: 'read' }],
  ['Expanded session file is too large.', { code: 'SESSION_SIZE_LIMIT', stage: 'read' }],
  ['This browser does not support gzip session export.', { code: 'SESSION_BROWSER_UNSUPPORTED', stage: 'transport' }],
  ['This browser does not support gzip session import.', { code: 'SESSION_BROWSER_UNSUPPORTED', stage: 'read' }]
]);
for (const [message, field, reason] of [
  ['Depth minimum must be less than or equal to depth maximum.', 'min_depth', 'ORDER'],
  ['Depth window must be greater than 0.', 'window', 'POSITIVE'],
  ['Depth step must be greater than 0.', 'step', 'POSITIVE'],
  ['Depth large tick interval must be greater than 0.', 'large_tick_interval', 'POSITIVE'],
  ['Depth small tick interval must be greater than 0.', 'small_tick_interval', 'POSITIVE'],
  ['Depth tick font size must be greater than 0.', 'tick_font_size', 'POSITIVE'],
  ['GC content minimum percent must be a finite number.', 'min_gc', 'FINITE'],
  ['GC content maximum percent must be a finite number.', 'max_gc', 'FINITE'],
  ['GC content minimum percent must be less than or equal to maximum percent.', 'min_gc', 'ORDER'],
  ['GC content large tick interval must be greater than 0.', 'large_tick_interval', 'POSITIVE'],
  ['GC content small tick interval must be greater than 0.', 'small_tick_interval', 'POSITIVE'],
  ['GC content tick font size must be greater than 0.', 'tick_font_size', 'POSITIVE'],
  ['Pairwise Match Height must be Auto or a positive finite number.', 'match_height', 'POSITIVE_OR_AUTO']
]) NATIVE_VALIDATIONS.set(message, { code: 'INPUT_INVALID', stage: 'request-validation', context: { field, reason } });
for (const [message, code] of [
  ['Please upload a GenBank file.', 'INPUT_REQUIRED'], ['GFF3 and FASTA are required.', 'FASTA_REQUIRED'],
  ['Please upload a Depth TSV file or disable Show depth track.', 'DEPTH_INVALID'],
  ['Please upload at least one Depth TSV file or disable the depth track.', 'DEPTH_INVALID'],
  ['Please upload at least one BLAST outfmt 6/7 file for Pairwise Comparisons.', 'COMPARISON_INPUT'],
  ['Please upload at least one comparison FASTA file for Pairwise Comparisons.', 'COMPARISON_INPUT'],
  ['No records found', 'NO_RECORDS'], ['No records found for circular conservation reference.', 'NO_RECORDS'],
  ['A resolved Linear comparison plan is required.', 'COMPARISON_INPUT'],
  ['The diagram engine returned an invalid Result list.', 'RESULT_INVALID'],
  ['The generated artifact has no selected preview Result.', 'RESULT_INVALID']
]) NATIVE_VALIDATIONS.set(message, { code, stage: code === 'RESULT_INVALID' ? 'result-admission' : 'request-validation' });
// Fixed labels/options emitted by the existing current-option validators.
for (const [label, field, reason] of [
  ['Protein BLASTP Pairwise max hits', 'protein_blastp_max_hits', 'POSITIVE_INTEGER'],
  ['Protein BLASTP member hits per protein', 'orthogroup_member_max_hits', 'POSITIVE_INTEGER'],
  ['Collinear minimum anchors', 'collinear_min_anchors', 'POSITIVE_INTEGER'],
  ['Collinear maximum unit gap', 'collinear_max_unit_gap', 'NONNEGATIVE_INTEGER'],
  ['Collinear maximum diagonal drift', 'collinear_max_diagonal_drift', 'NONNEGATIVE_INTEGER'],
  ['Collinear maximum conflicts in merge gap', 'collinear_max_conflicts', 'NONNEGATIVE_INTEGER'],
  ['Collinear maximum paralog links per orthogroup', 'collinear_max_paralog_links_per_orthogroup', 'POSITIVE_INTEGER']
]) NATIVE_VALIDATIONS.set(`${label} must be ${reason === 'NONNEGATIVE_INTEGER' ? 'a non-negative integer' : 'a positive integer'}.`,
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field, reason } });
NATIVE_VALIDATIONS.set('Protein BLASTP Max target seqs must be a positive integer or None.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'protein_blastp_candidate_limit', reason: 'POSITIVE_OR_AUTO' } });
NATIVE_VALIDATIONS.set('Collinear orthogroup inference must be a boolean.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'collinear_infer_orthogroups', reason: 'BOOLEAN' } });
for (const [label, field, options] of [
  ['Circular multi-record size mode', 'circular_multi_record_size_mode', 'auto, linear, equal'],
  ['Linear track layout', 'linear_track_layout', 'above, middle, below'],
  ['Linear label placement', 'linear_label_placement', 'auto, above_feature'],
  ['Protein BLASTP mode', 'protein_blastp_mode', 'pairwise, orthogroup, collinear'],
  ['Collinear search scope', 'collinear_search_scope', 'adjacent, all'],
  ['Collinear unit mode', 'collinear_unit_mode', 'auto, cds, locus'],
  ['Collinear anchor mode', 'collinear_anchor_mode', 'all, one_to_one, rbh'],
  ['Collinear merge orientation', 'collinear_merge_orientation', 'strand, order, either'],
  ['Collinear color mode', 'collinear_color_mode', 'average_identity, orientation, orientation_identity'],
  ['Protein BLASTP orthogroup membership mode', 'orthogroup_membership_mode', 'anchor_core_v1']
]) {
  NATIVE_VALIDATIONS.set(`${label} must be one of: ${options}.`,
    { code: 'INPUT_INVALID', stage: 'request-validation', context: { field, reason: field.toUpperCase() } });
}
for (const message of [
  'No SVG result is available for export.', 'Interactive SVG export requires the committed feature catalog.',
  'The SVG export produced no content.', 'The current SVG could not be prepared for PDF export.',
  'The current SVG could not be parsed for PNG export.', 'The current SVG could not be parsed for PDF export.'
]) NATIVE_VALIDATIONS.set(message, { code: 'EXPORT_INPUT', stage: 'export-capture' });
for (const format of ['PNG', 'PDF']) NATIVE_VALIDATIONS.set(`The current SVG has no usable dimensions for ${format} export.`,
  { code: 'EXPORT_DIMENSIONS', operation: `export-${format.toLowerCase()}`, stage: 'export-capture' });
NATIVE_VALIDATIONS.set('The selected PNG DPI is invalid.', { code: 'PNG_DPI', operation: 'export-png', stage: 'export-conversion' });
for (const message of ['The browser could not initialize PNG conversion.', 'The browser could not load the SVG for PNG export.',
  'The browser produced no PNG export data.']) NATIVE_VALIDATIONS.set(message,
  { code: 'EXPORT_CONVERSION', operation: 'export-png', stage: 'export-conversion' });
for (const message of ['The vendored jsPDF library did not initialize.', 'The vendored svg2pdf library did not initialize.'])
  NATIVE_VALIDATIONS.set(message, { code: 'PDF_LIBRARY', operation: 'export-pdf', stage: 'initialization' });
for (const message of [
  'Select an exact reference feature before aligning this Similarity Group.',
  'An exact reference feature is required for alignment.',
  'The resolver did not resolve every record. Review the choices and retry.'
]) NATIVE_VALIDATIONS.set(message, { code: 'INPUT_INVALID', operation: 'align', stage: 'helper', context: { field: 'anchors', reason: 'TARGET_RECORD' } });
NATIVE_VALIDATIONS.set('Validated directions or reference placement changed. Review the updated preview and Apply again.',
  { code: 'INPUT_INVALID', operation: 'align', stage: 'helper', context: { field: 'direction', reason: 'SOURCE_VIEW_CONFLICT' } });
for (const message of ['Failed to load the vendored jsPDF library.', 'Failed to load the vendored svg2pdf library.'])
  NATIVE_VALIDATIONS.set(message, { code: 'PDF_LIBRARY', operation: 'export-pdf', stage: 'initialization' });
NATIVE_VALIDATIONS.set('A sequence file is required.', { code: 'INPUT_REQUIRED', stage: 'request-validation' });
NATIVE_VALIDATIONS.set('GFF3 and FASTA files are required.', { code: 'FASTA_REQUIRED', stage: 'request-validation' });
for (const message of ['Current session is missing its active form or advanced settings.',
  'Current session is missing its canonical configuration projection.']) NATIVE_VALIDATIONS.set(message,
    { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'config', reason: 'FIELDS' } });
// Native Session binding validation emits these finite messages. Resource IDs
// are discarded before the public model is constructed.
for (const message of [
  'Session webFiles.bindings must be an object.', 'Unsupported Web file binding schema.',
  'Web binding schema 2 requires c_gb.', 'Invalid Web file binding metadata.',
  'Web file bindings must be objects or null.',
  'A composite component must be an ordinary Web file binding.',
  'A Web file binding requires a canonical resourceId.', 'Unknown Web file binding field.',
  'Unsupported Web composite file binding.', 'Unknown or mixed composite binding fields.',
  'A composite binding requires at least two components.', 'Composite components require base64 resources.',
  'Unknown Web binding inventory field.'
]) NATIVE_VALIDATIONS.set(message, { code: 'INPUT_INVALID', stage: 'request-validation',
  context: { field: 'schema', reason: 'FIELDS' } });
for (const message of [
  'A canonical render request is required for adoptive ownership.',
  'Settings-only Session cannot contain biological sources.',
  'Settings-only Session cannot contain committed render artifacts.'
]) NATIVE_VALIDATIONS.set(message, { code: 'INPUT_INVALID', stage: 'request-validation',
  context: { field: 'schema', reason: 'FIELDS' } });
NATIVE_VALIDATIONS.set('A diagram generation request is already running.',
  { code: 'GENERATION_BUSY', stage: 'render' });
NATIVE_VALIDATIONS.set('Arrow head length ratio must be Auto or a positive finite number.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'arrow_head_length_ratio', reason: 'POSITIVE_OR_AUTO' } });
NATIVE_VALIDATIONS.set('Arrow shaft width ratio must be a finite number greater than 0 and at most 1.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'arrow_shaft_width_ratio', reason: 'POSITIVE_UNIT_INTERVAL' } });
NATIVE_VALIDATIONS.set('Normalize Record Lengths cannot be used when multiple records share the same Linear row. Turn Normalize off or assign each record to a separate row.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'config', reason: 'SEPARATE_LINEAR_ROWS' } });
NATIVE_VALIDATIONS.set('Current session active configuration config.form.keep_definition_left_aligned must be a boolean.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'keep_definition_left_aligned', reason: 'BOOLEAN' } });
for (const [message, reason] of [
  ['Circular region requires both Start and End coordinates.', 'BOTH_ENDPOINTS'],
  ['Circular region Start and End must be positive integers.', 'POSITIVE_INTEGER'],
  ['Circular region Start must not exceed End. Use Reverse complement to change display orientation.', 'ORDER']
]) NATIVE_VALIDATIONS.set(message, { code: 'REGION_INVALID', stage: 'request-validation',
  context: { field: 'region', reason } });
NATIVE_VALIDATIONS.set('Canonical renderRequest records are required.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'RECORDS_REQUIRED' } });
for (const message of ['The preserved comparison is missing a required resource.',
  'The saved comparison is missing a required resource.',
  'The saved Circular comparison is missing a required resource.']) NATIVE_VALIDATIONS.set(message,
  { code: 'COMPARISON_INPUT', stage: 'request-validation', context: { field: 'comparison', reason: 'REQUIRED' } });
const nativeValidation = (message) => {
  if (typeof message !== 'string') return null;
  if (/^Circular region End \([0-9]+\) exceeds the selected record length \([0-9]+\)\.$/.test(message)) return { code: 'REGION_INVALID', stage: 'request-validation', context: { field: 'region', reason: 'RECORD_BOUNDS' } };
  for (const [template, reason] of [
    [/^Session resource [\s\S]* has an invalid declared byte size\.$/, 'RESOURCE_SIZE'],
    [/^Custom Track Slots use an obsolete schema\. Recreate the slots with schema version [0-9]+\.$/, 'TRACK_SCHEMA'],
    [/^Session contains unclassified top-level field\(s\): [\s\S]*$/, 'SESSION_FIELDS']
  ]) if (template.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason } };
  if (/^Session version [0-9]+ is newer than this gbdraw supports \([0-9]+\)\.$/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'FIELDS' } };
  if (/^Invalid managed flag for (?:circular|linear)\.[a-z_]+\.$/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'config', reason: 'FIELDS' } };
  if (/^Missing canonical resource:/.test(message) || /^Session resource [\s\S]* has an unsupported encoded payload\.$/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'FIELDS' } };
  if (/^The SVG composition metadata is not valid JSON:/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'JSON_FORMAT' } };
  if (/^Unsupported session version: [0-9]+\.$/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'FIELDS' } };
  if (NATIVE_VALIDATIONS.has(message)) return NATIVE_VALIDATIONS.get(message);
  if (/^Could not read records from the (?:circular|Linear) input file\(s\)\.$/.test(message)) {
    return { code: 'INPUT_UNREADABLE', stage: 'request-validation' };
  }
  if (/^Invalid region (?:spec|coordinates) for LOSAT FASTA extraction: [\s\S]*$/.test(message)) {
    return { code: 'REGION_INVALID', stage: 'request-validation', context: { reason: 'REGION_FORMAT' } };
  }
  const pixel = /^(?:Circular|Linear) track slot '[\s\S]*' (height|spacing|inner_gap_px|outer_gap_px) must be (nonnegative|positive) finite number of pixels \(px optional\)\.$/.exec(message);
  if (pixel) return { code: 'TRACK_INVALID', stage: 'request-validation', context: {
    field: pixel[1], reason: pixel[2] === 'positive' ? 'PIXEL_POSITIVE' : 'PIXEL_NONNEGATIVE' } };
  const scalar = /^Circular track slot '[\s\S]*' (radius|width) must be a positive finite px or factor scalar\.$/.exec(message);
  if (scalar) return { code: 'TRACK_INVALID', stage: 'request-validation', context: { field: scalar[1], reason: 'POSITIVE_SCALAR' } };
  if (/^Circular track slot '[\s\S]*' uses obsolete field '(?:spacing|strict|compress|reserve)'\. Use inner_gap_px and outer_gap_px for physical gaps\.$/.test(message)) return { code: 'TRACK_INVALID', stage: 'request-validation', context: { reason: 'CIRCULAR_GAPS' } };
  const order = /^Start position \([0-9]+\) must be less than end position \([0-9]+\)\.$/.test(message);
  if (order) return { code: 'REGION_INVALID', stage: 'request-validation', context: { reason: 'STRICT_ORDER' } };
  const selector = /^Record selector #[0-9]+ is out of range \(loaded ([0-9]+) record\(s\)\)\.$/.exec(message);
  if (selector) return { code: 'RECORD_SELECTION', stage: 'request-validation', context: { reason: 'OUT_OF_RANGE', recordCount: Number(selector[1]) } };
  for (const [template, reason] of [
    [/^Record selector '[\s\S]*' did not match any record ID\.$/, 'NO_MATCH'],
    [/^Record selector '[\s\S]*' matched multiple records\. Use #index to disambiguate\.$/, 'AMBIGUOUS']
  ]) if (template.test(message)) return { code: 'RECORD_SELECTION', stage: 'request-validation', context: { reason } };
  const slot = /^(?:Circular|Linear) depth slot '[\s\S]*' references removed depth track index (unknown|[0-9]+)\. Select an existing Depth TSV or remove the slot\.$/.exec(message);
  if (slot) return { code: 'DEPTH_INVALID', stage: 'request-validation', context: {
    ...(slot[1] === 'unknown' ? {} : { seriesIndex: Number(slot[1]) }), reason: 'DEPTH_SERIES' } };
  for (const [template, reason] of [
    [/^Region annotation [\s\S]*: start and end must be positive integers \(1-based coordinates\)\.$/, 'POSITIVE_INTEGER'],
    [/^Choose a (?:valid )?target record for region annotation [\s\S]*\.$/, 'TARGET_RECORD'],
    [/^Region annotation [\s\S]*: clear the target record or enable Multi-record canvas\.$/, 'TARGET_MODE'],
    [/^Region annotation [\s\S]*: the selected target record is no longer available; choose it again\.$/, 'TARGET_RECORD'],
    [/^Region annotation [\s\S]*: record index #[0-9]+ is out of range\.$/, 'OUT_OF_RANGE'],
    [/^Region annotation [\s\S]*: record ID [\s\S]* is not available\.$/, 'NO_MATCH'],
    [/^Region annotation [\s\S]*: record ID [\s\S]* is duplicated; choose a #index entry\.$/, 'AMBIGUOUS']
  ]) if (template.test(message)) return { code: 'ANNOTATION_TARGET', stage: 'request-validation', context: { reason } };
  const emptyFasta = /^Pairwise comparison FASTA #([0-9]+) has no sequence data\.$/.exec(message);
  if (emptyFasta) return { code: 'COMPARISON_INPUT', stage: 'request-validation', context: { inputOrdinal: Number(emptyFasta[1]), reason: 'REQUIRED' } };
  const glyph = /^PDF fonts do not contain U\+([0-9A-F]{1,6})\. Use SVG to retain this text\.$/.exec(message);
  if (glyph) return { code: 'PDF_GLYPH', operation: 'export-pdf', stage: 'font-validation', context: { codepoint: parseInt(glyph[1], 16) } };
  const sequence = /^Sequence #([0-9]+): (Provide both Region start and end, or leave both empty\.|Region start\/end must be numbers\.|Region start\/end must be integers\.|Region start\/end must be >= 1\.|Missing GenBank file\.|GFF3 and FASTA are required\.)$/.exec(message);
  if (sequence) {
    const reasons = { 'Provide both Region start and end, or leave both empty.': 'BOTH_ENDPOINTS',
      'Region start/end must be numbers.': 'FINITE', 'Region start/end must be integers.': 'INTEGER',
      'Region start/end must be >= 1.': 'POSITIVE_INTEGER' };
    return { code: sequence[2] === 'Missing GenBank file.' ? 'INPUT_REQUIRED' : sequence[2] === 'GFF3 and FASTA are required.' ? 'FASTA_REQUIRED' : 'REGION_INVALID',
      stage: 'request-validation', context: { inputOrdinal: Number(sequence[1]), reason: reasons[sequence[2]] } };
  }
  const row = /^Invalid specific-color TSV at line ([0-9]+): (expected 4 or 5 columns\.|column ([1-4]) is required\.)$/.exec(message);
  if (row) return { code: 'TABLE_INVALID', stage: 'request-validation', context: { row: Number(row[1]),
    field: ['feature_type', 'qualifier', 'pattern', 'color'][Number(row[3]) - 1],
    reason: row[3] ? 'REQUIRED' : 'SPECIFIC_COLUMNS' } };
  const labelColumns = /^Invalid label TSV at line ([0-9]+): expected 5 columns, found ([0-9]+)\.$/.exec(message);
  if (labelColumns) return { code: 'TABLE_INVALID', stage: 'request-validation', context: { row: Number(labelColumns[1]), columnCount: 5, reason: 'FIELDS' } };
  const labelRequired = /^Invalid label TSV at line ([0-9]+): column ([1-4]) \((record_id|feature_type|qualifier|value)\) is required\.$/.exec(message);
  if (labelRequired) return { code: 'TABLE_INVALID', stage: 'request-validation', context: { row: Number(labelRequired[1]), field: labelRequired[3], reason: 'REQUIRED' } };
  const color = /^Invalid specific-color value at line ([0-9]+): [\s\S]*$/.exec(message);
  if (color) return { code: 'TABLE_INVALID', stage: 'request-validation', context: { row: Number(color[1]), field: 'color', reason: 'COLOR' } };
  const series = /^Depth series #([0-9]+) \(logical track index ([0-9]+)\) has no TSV source in any record\.(?: Add a TSV or remove the series\.)?$/.exec(message);
  if (series) return { code: 'DEPTH_INVALID', stage: 'request-validation', context: { seriesIndex: Number(series[2]), reason: 'REQUIRED' } };
  return null;
};
// The native track validator already identifies its failures. Project its
// first issue; never copy its message, slot ID, or arbitrary field path.
const TRACK_ISSUES = Object.freeze({
  slots_not_array: ['ARRAY'], mode_unsupported: ['TRACK_RENDERER'], depth_count_invalid: ['NONNEGATIVE_INTEGER', 'depth'],
  slot_not_object: ['OBJECT'], enabled_not_boolean: ['BOOLEAN'], linear_slots_empty: ['REQUIRED'],
  axis_out_of_range: ['TRACK_AXIS', 'axis'], id_required: ['REQUIRED'], id_duplicate: ['UNIQUE_IDS'],
  renderer_unsupported: ['TRACK_RENDERER', 'renderer'], side_unsupported: ['TRACK_SIDE', 'side'],
  params_not_object: ['OBJECT'], z_invalid: ['INTEGER', 'z'],
  generic_param: ['TRACK_PARAMS'], renderer_param_mismatch: ['TRACK_PARAMS'], feature_lane: ['TRACK_LANE'],
  feature_lane_side_conflict: ['TRACK_LANE', 'side'], overlay_renderer_unsupported: ['TRACK_SIDE', 'side'],
  conservation_overlay: ['TRACK_SIDE', 'side'], ticks_obsolete_param: ['TICK_PARAMETER'],
  depth_track_index: ['NONNEGATIVE_INTEGER', 'depth'], depth_source_missing: ['REQUIRED', 'depth'],
  depth_track_index_range: ['DEPTH_SERIES', 'depth'], annotation_set_required: ['ANNOTATION_SET', 'set_id'],
  annotation_set_unknown: ['ANNOTATION_SET', 'set_id'], annotation_marks: ['ANNOTATION_MARKS'],
  annotation_lane_gap: ['NONNEGATIVE', 'lane_gap_px'], annotation_padding: ['NONNEGATIVE', 'padding_px'],
  annotation_cover_anchor: ['BOOLEAN', 'cover_anchor'], annotation_overflow: ['OVERFLOW', 'overflow'],
  annotation_layer: ['LAYER', 'layer'], annotation_anchor_required: ['OVERLAY_ANCHOR', 'anchor_slot'],
  annotation_anchor_without_overlay: ['OVERLAY_SIDE', 'side'], conservation_unmanaged: ['CONSERVATION_SOURCE'],
  conservation_source_missing: ['CONSERVATION_SOURCE'], features_multiple: ['FEATURES_COUNT'],
  feature_underlay_features_count: ['FEATURES_COUNT'], annotation_anchor_unknown: ['OVERLAY_ANCHOR', 'anchor_slot'],
  annotation_anchor_ineligible: ['OVERLAY_ANCHOR', 'anchor_slot'], annotation_underlay_z: ['DRAW_ORDER', 'z'],
  annotation_foreground_z: ['DRAW_ORDER', 'z'], axis_side_conflict: ['TRACK_SIDE', 'side']
});
const typedTrackValidation = (error) => {
  const issue = Array.isArray(error.issues) ? error.issues[0] : null;
  if (issue?.code === 'geometry_invalid') {
    const geometry = nativeValidation(issue.message);
    return geometry?.code === 'TRACK_INVALID'
      ? { ...geometry, context: { ...geometry.context, slotIndex: issue.rowIndex } }
      : { code: 'VALIDATION_UNCLASSIFIED', stage: 'request-validation', context: {} };
  }
  const correction = issue && typeof issue.code === 'string' && issue.code.length <= 80 && Object.hasOwn(TRACK_ISSUES, issue.code) ? TRACK_ISSUES[issue.code] : null;
  return { code: correction ? 'TRACK_INVALID' : 'VALIDATION_UNCLASSIFIED', stage: 'request-validation',
    context: correction ? { reason: correction[0], field: correction[1], slotIndex: issue.rowIndex } : {} };
};
const identifier = (value, domain, fallback) => typeof value === 'string' && value.length <= 80 && domain.has(value) ? value : fallback;
const contextFor = (value) => {
  const context = {};
  if (!value || typeof value !== 'object' || Array.isArray(value)) return context;
  if (FIELDS.has(value.field)) context.field = value.field;
  if (typeof value.reason === 'string' && Object.hasOwn(REASONS, value.reason)) context.reason = value.reason;
  if (value.positionUnit === 'python-character') context.positionUnit = value.positionUnit;
  for (const key of ['position', 'row', 'inputOrdinal', 'recordIndex', 'seriesIndex', 'slotIndex', 'recordCount', 'columnCount', 'codepoint']) {
    if (Number.isSafeInteger(value[key]) && value[key] >= 0 && value[key] <= (key === 'codepoint' ? 0x10ffff : 10000000)) {
      if (key !== 'position' || context.positionUnit === 'python-character') context[key] = value[key];
    }
  }
  return context;
};

/** The sole public diagnostic model. Raw messages, causes, logs and stacks are never copied. */
export const normalizeUserFacingError = (value, {
  operation = 'unknown', stage = 'unknown', code = 'UNKNOWN', summaryLimit = 1000, detailLimit = 4000
} = {}) => {
  if (!value) return null;
  const object = value && typeof value === 'object' ? value : {};
  const source = object.code === 'CUSTOM_TRACK_PLAN_INVALID' ? typedTrackValidation(object) : CODES.has(object.code) ? object : nativeValidation(typeof value === 'string' ? value : object.message) || object;
  const result = {
    code: identifier(source.code, CODES, identifier(code, CODES, 'UNKNOWN')),
    operation: identifier(operation === 'unknown' ? source.operation : operation, OPERATIONS, 'unknown'),
    stage: identifier(source.stage === 'unknown' ? stage : source.stage, STAGES, identifier(stage, STAGES, 'unknown')),
    context: contextFor(source.context),
    secondary: []
  };
  for (const item of Array.isArray(source.secondary) ? source.secondary.slice(0, 2) : []) {
    if (item?.code !== 'CLEANUP_FAILED' || item?.stage !== 'cleanup') continue;
    result.secondary.push({ code: 'CLEANUP_FAILED', stage: 'cleanup' });
  }
  const [message, actions] = DEFINITIONS[result.code];
  const guidance = REASONS[result.context.reason] || '';
  const field = result.context.field ? ` Field: ${FIELD_LABELS[result.context.field] || result.context.field}.` : '';
  const position = result.context.position !== undefined ? ` Python character position ${result.context.position} (zero-based).` : '';
  const columns = result.context.columnCount !== undefined ? ` Required columns: ${result.context.columnCount}.` : '';
  const count = result.context.recordCount !== undefined ? ` Loaded records: ${result.context.recordCount}.` : '';
  const schemaGuidance = result.code === 'INPUT_INVALID' && result.context.field === 'schema'
    ? ' Load a supported Session file or recreate it with the current writer.' : '';
  const continuation = result.operation === 'align' && result.context.field === 'direction'
    && result.context.reason === 'SOURCE_VIEW_CONFLICT' ? ' Review the updated preview and Apply again.' : '';
  const decorationResult = result.code === 'DECORATION_CONTINUITY' && result.context.inputOrdinal !== undefined
    ? ` Result ${result.context.inputOrdinal}.` : '';
  result.summary = `${message}${decorationResult}${field}${guidance ? ` ${guidance}` : ''}${position}${columns}${count}${continuation}${schemaGuidance}`
    .slice(0, Number.isSafeInteger(summaryLimit) ? Math.max(0, Math.min(summaryLimit, 1000)) : 1000);
  const detail = [`Code: ${result.code}`, `Operation: ${result.operation}`, `Stage: ${result.stage}`,
    ...Object.entries(result.context).map(([key, item]) => `${key}: ${item}`),
    ...result.secondary.map((item) => `Secondary: ${item.code} / ${item.stage}`)].join('\n');
  result.details = [{ label: 'Diagnostics', text: detail.slice(0,
    Number.isSafeInteger(detailLimit) ? Math.max(0, Math.min(detailLimit, 4000)) : 4000) }].slice(0, 8);
  result.actions = [...actions];
  return result;
};

// Presentation of finite operation and transaction facts stays with public wording.
export const operationErrorTitle = (operation) => ({
  generate: 'Generation Error', align: 'Alignment error', evaluateRules: 'Rule error',
  'export-svg': 'SVG export error', 'export-png': 'PNG export error', 'export-pdf': 'PDF export error',
  listSequenceRecords: 'Input error', listGffFastaRecords: 'Input error',
  'feature-extraction': 'Feature preparation error'
}[operation] || 'Operation error');
export const generationRecoveryGuidance = (recovery) => ({
  'no-result': 'No successful Result is available yet. Correct the cause and retry.',
  preserved: 'The last successful Result and committed request are unchanged.',
  restored: 'The last successful Result and committed request were restored. Correct the cause and retry.',
  'restore-failed': 'Restoration could not complete. Save a Session if available, then reload before retrying.'
}[recovery] || '');
