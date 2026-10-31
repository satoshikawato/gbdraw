// @ts-check
// This module has no imports (tests/web/helpers/error-text-keys.mjs loads a
// copy of it alone), so it declares the identity of utils/i18n.js `msg`.
/** @param {string} text */
const msg = (text) => text;
const OPERATIONS = new Set(['unknown', 'generate', 'align', 'feature-extraction', 'export-svg', 'export-png', 'export-pdf', 'evaluateRules', 'embedFonts', 'pdfFonts', 'buildProteinLosatCacheKeys', 'computeLabelFitRow', 'convertMainSessionComparisonFrame', 'convertRingTablesToSearchFrame', 'convertStoredProteinResultFrames',
  'convertLosatpPairsToGenomicPayload', 'extractCdsProteinFasta', 'extractFirstFasta',
  'hydrateProteinLosatTsv', 'listGffFastaRecords', 'listSequenceRecords',
  'promoteLegacyLosatpCache', 'readAllSequenceRecords', 'readComparisonSequence', 'readFeatureOverrideTable', 'registerFonts', 'resolveLegacyProteinReferences',
  'resolveSimilarityAlignment', 'synthesizeColorRules', 'validateConfigOverrides', 'session-save', 'session-load']);
const STAGES = new Set(['unknown', 'initialization', 'resource-staging', 'request-validation',
  'helper', 'rule-validation', 'losat', 'render', 'result-admission', 'cleanup', 'export-capture', 'export-conversion', 'font-validation', 'transport', 'read', 'write', 'parse']);
const FIELDS = new Set(`legend title scale decorations pattern record_selector region start end sourceStart sourceEnd recordLength
recordIndex queryIndex subjectIndex depth min_depth max_depth window step tick font_size font_family plot_title_font_size height large_tick_interval small_tick_interval tick_font_size
inner_gap_px outer_gap_px radius width spacing arrow_head_length_ratio arrow_shaft_width_ratio keep_definition_left_aligned color action feature_type qualifier value record_id label_text
config configOverrides records anchors schema recordKey groupId direction sourceStrand role blast files
input comparison protein_blastp_max_hits orthogroup_member_max_hits bitscore evalue identity
losatp_max_hits losatp_max_target_seqs losatp_member_max_hits losatp_mode losatn_task record_gencodes
losat losat_search program pairs source losat_gencode comparison_sequence comparison_sequence_source comparison_fasta
conservation_blast_files conservation_dataframes conservation_reference conservation_sequence_files
conservation_labels conservation_colors conservation_losat_gencodes conservation_search_results
alignment_length wrap_length collinear_min_anchors collinear_max_gene_gap collinear_block_merge_gap
collinear_singleton_merge_gap collinear_max_diagonal_drift collinear_gap_penalty
collinear_nearby_duplicate_window collinear_constant_anchor_score collinear_infer_orthogroups
collinear_min_score collinear_min_block_span record_gap_px record_axis_height depth_window depth_step
depth_min depth_max min_gc max_gc gc_tick_interval gc_axis_font_size depth_tick_interval
depth_axis_font_size dinucleotide feature_width_circular depth_width_circular gc_content_width_circular
gc_content_radius_circular gc_skew_width_circular gc_skew_radius_circular conservation_ring_width
conservation_ring_gap center_reserved_radius multi_record_min_radius_ratio multi_record_column_gap_ratio
multi_record_row_gap_ratio protein_blastp_mode protein_blastp_candidate_limit collinear_search_scope collinear_unit_mode collinear_anchor_mode collinear_merge_orientation collinear_color_mode orthogroup_membership_mode collinear_max_unit_gap collinear_max_conflicts collinear_max_paralog_links_per_orthogroup circular_multi_record_size_mode linear_track_layout linear_label_placement set_id anchor_slot to_slot side renderer lane_gap_px padding_px cover_anchor overflow layer z axis match_height source fasta gff annotations featurePlacements output_prefix
placement strict compress reserve gapAfter gap_after innerRadius inner_radius outerRadius outer_radius`.split(/\s+/));
/** @type {Readonly<Record<string, string>>} */
const REASONS = Object.freeze({
  CROSS_ORIGIN_ISOLATION: 'This page is not cross-origin isolated.', SHARED_MEMORY: 'This browser does not provide SharedArrayBuffer.',
  WORKERS: 'This browser does not provide Web Workers.', THREADED_WASM: 'The threaded LOSAT runtime could not start.',
  DECORATION_IDENTITY: 'Source, region, mode, grouping or record identity changed, is unknown, or is ambiguous.',
  DECORATION_TARGET: 'The decoration target is missing or ambiguous.',
  DECORATION_METADATA: 'Composition metadata is invalid or unavailable.',
  ANNOTATION_SET: 'Select an existing annotation set.', ANNOTATION_MARKS: 'Use supported annotation mark filters.',
  OVERFLOW: 'Choose error, compress, or clip.', LAYER: 'Choose foreground or underlay.',
  OVERLAY_ANCHOR: 'Select an eligible enabled anchor for an annotation track whose placement is Overlay anchor (Circular) or On axis (Linear).',
  OVERLAY_SPAN: 'For Extend to, select an eligible enabled slot other than the anchor, or None.',
  OVERLAY_SIDE: 'Choose the placement Overlay anchor (Circular) or On axis (Linear) when an anchor or an Extend to slot is selected.',
  DRAW_ORDER: 'Place underlay annotations below the anchor and the Extend to slot, and foreground annotations above them.',
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
  COLLINEAR_ANCHOR_MODE: 'Choose all or rbh.',
  COLLINEAR_MERGE_ORIENTATION: 'Choose strand or order.',
  COLLINEAR_COLOR_MODE: 'Choose average_identity, orientation, or orientation_identity.',
  ORTHOGROUP_MEMBERSHIP_MODE: 'Choose anchor_core_v1.',

  UNKNOWN_CONFIG_PATH: 'Remove unknown configuration overrides or use a supported setting.',
  CIRCULAR_SETTING: 'It applies only to Circular diagrams. Reset it under Preserved session settings, or switch to Circular.',
  WRAP: 'A wrapped figure draws no comparisons. Remove the comparison inputs or set the wrap length to none.',
  LINEAR_SETTING: 'It applies only to Linear diagrams. Reset it under Preserved session settings, or switch to Linear.',
  TRACK_SCHEMA: 'Recreate Custom Track Slots with the current schema.',
  SESSION_FIELDS: 'Remove unsupported top-level Session fields.',
  RECORDS_REQUIRED: 'Session render requests require at least one record.',
  POSITIVE_UNIT_INTERVAL: 'Use a finite number greater than zero and at most one.',
  PIXEL_POSITIVE: 'Use a finite number of pixels greater than zero (px optional).',
  PIXEL_NONNEGATIVE: 'Use a finite number of pixels of zero or greater (px optional).',
  POSITIVE_SCALAR: 'Use a number greater than 0, in px or ×R.',
  CIRCULAR_GAPS: 'Use inner_gap_px and outer_gap_px for physical gaps.',
  OBSOLETE_TRACK_FIELD: 'Custom Track Slots no longer read this field. Use slot-level radius, width, inner_gap_px, outer_gap_px, side, and z fields.',
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
  SEARCH_FRAME: 'Use comparison coordinates within the selected and cropped record, counted on its source strand.',
  LOSAT_OPTION_PROGRAM: 'Use only the options of the selected LOSAT program.',
  LOSAT_PLAN: 'Choose record pairs that the LOSAT search can run: two or more records, LOSAT rows only with a LOSAT program, and LOSATP Pairwise for selected pairs.',
  // A Selected plan with a LOSATP mode that compares every adjacent pair.
  SELECTED_PAIRS_NEED_PAIRWISE: 'Similarity groups and Collinear blocks compare every adjacent pair with LOSAT, so they do not run on pairs set one by one in Selected pairs. Choose Run LOSAT for all adjacent pairs, or choose Pairwise matches.',
  LOSAT_TASK: 'Choose a LOSATN task that the selected LOSAT runtime supports.',
  // The LOSATP converter named another derived identity than Generate looked up (run-analysis.js).
  DERIVED_KEY_MISMATCH: 'The converted LOSATP comparison does not match the settings it was looked up with. If Retry fails the same way, save the Session and report the problem.',
  RING_LOSAT_PROGRAM: 'Choose LOSATN or TLOSATX for comparison rings.',
  RING_LOSAT_INPUT: 'Give one comparison sequence file per ring, without BLAST tables, and use the reference file as the subject.',
  SEQUENCE_MISSING: 'Use a comparison file that contains the sequence (FASTA, or GenBank or DDBJ with an ORIGIN sequence).',
  UNAVAILABLE: 'Install LOSAT or NCBI BLAST+, or select a runtime executable.',
  FAILED: 'The search runtime exited with an error; check its message and inputs.',
  OUTPUT: 'The search output was not BLAST outfmt 6 rows for the searched records.',
  JSON_FORMAT: 'Use valid JSON.', VISIBILITY_ACTION: 'Use show, off, or exclude_matching; accepted aliases are on, hide, false, and 0.',
  TARGET_RECORD: 'Choose an available target record.',
  DEPTH_SERIES: 'Select an existing Depth TSV or remove the slot.',
  DEPTH_SERIES_SOURCE: 'Attach a Depth TSV to this series, or remove the series.',
  BOTH_ENDPOINTS: 'Supply both region endpoints or leave both empty.', SPECIFIC_COLUMNS: 'Supply four or five tab-separated columns.',
  COLOR: 'Use none, a supported named color, or a hex color with 3 or 6 digits. Feature type rows of a Default colors table and configuration colors other than conservation also accept transparent, rgb(), and hsl(); comparison gradient colors take #RGB or #RRGGBB.',
  // Output Prefix (gbdraw/api/requests.py RenderOutputRequest).
  FILENAME: 'Use one file name without a folder: no / \\ : * ? " < > | or control characters, no trailing dot, and not a Windows device name such as CON or NUL.',
  FILENAME_LENGTH: 'Use at most 200 bytes in UTF-8, for example 200 ASCII characters.',
  FONT_NOT_LOADED: 'This font is no longer loaded in this page. Load it again in Font, or choose another font.',
  BOOLEAN: 'Use true or false.', INTEGER: 'Use an integer.', NONNEGATIVE: 'Use a finite value of zero or greater.',
  FINITE: 'Use a finite number.', POSITIVE_INTEGER: 'Use an integer greater than zero.',
  POSITIVE: 'Use a finite value greater than zero.', POSITIVE_OR_AUTO: 'Use Auto or a finite value greater than zero.',
  NONNEGATIVE_INTEGER: 'Use an integer of zero or greater.', PERCENT: 'Use a finite value between 0 and 100.',
  ARRAY: 'Use a list.', OBJECT: 'Use an object.', FIELDS: 'Check the required fields.', REQUIRED: 'Supply the required value.',
  ORDER: 'The start or minimum must not exceed the end or maximum.', RECORD_BOUNDS: 'Keep the region within the record length.',
  STRAND: 'Use -1, 1, or no strand.', DISPLAY_START_BOUNDS: 'Use a display start between 1 and the record length.', CROP_START_CONFLICT: 'Choose a crop or an explicit display start.',
  REFERENCE_REQUIRED: 'Supply the depth reference column.', REFERENCE_MISMATCH: 'Match depth references to the selected record.',
  WORKSPACE: 'Retry the operation.', OUT_OF_RANGE: 'Choose a record within the loaded range.',
  NO_MATCH: 'Choose an available record.', AMBIGUOUS: 'Use #<number> to distinguish records with the same ID.',
  SELECTOR_FORMAT: 'Use #<number> or a record ID.', SELECT_ONE: 'Select exactly one record.',
  NONE_DRAWN: 'Turn on at least one record to draw.',
  REGION_FORMAT: 'Use record_id:start-end[:rc] or #<number>:start-end[:rc].',
  CANNOT_FIT: 'Move the track, reduce widths, disable conflicting labels, or place it on the other side of the Axis.',
  CENTER_RESERVED: 'The Center Reserved Radius limits the inside tracks. Set a smaller Center Reserved Radius or place tracks outside.',
  SPLIT_LANES: 'Set that feature\'s Feature placement to Auto or Main, or use split feature lanes (Track Preset Middle).',
  OVERLAY_LANES: 'Set that feature\'s Feature placement to Auto or Main, or use Track Layout Features on axis with Separate Strands off.',
  // OV-1262: gbdraw/annotations/layout.py, an explicit lane under Overflow "Error on overflow".
  ANNOTATION_LANE_CONFLICT: 'Two region annotations overlap in the same lane of the annotation table, counting their labels and padding. Give one of them another lane in Region Annotations, or set the Annotations row\'s Overflow behavior to Compress or Clip in Custom Track Slots (Use custom stack).',
  THREE_COLUMNS: 'Supply at least three tab-separated columns.',
  POSITIVE_INTEGER_OR_AUTO: 'Use Auto or an integer greater than zero.',
  DINUCLEOTIDE: 'Use two letters from A, C, G, T, and U, for example GC.',
  DEPTH_VALUES: 'Use integer positions and numeric depth values.',
  SELECT_RECORD_FOR_REGION: 'Choose a Record before setting a region on a multi-record file.',
  DISCOVERY_PENDING: 'Wait for the record list to finish loading, then retry.',
  SESSION_FORMAT: 'Choose a gbdraw Session file.',
  SESSION_TOO_NEW: 'This Session was saved by a newer gbdraw.',
  SESSION_VERSION_UNSUPPORTED: 'This gbdraw does not read this Session version.',
  SESSION_TOO_OLD: 'This Session was saved by gbdraw 0.11.0 or older, or by gbdraw.app before 0.12.0, and this version cannot read it. Generate the diagram again from its input files.',
  FORMAT: 'Use GenBank or the required GFF3 and FASTA inputs.', NO_PROTEINS: 'Choose input containing CDS proteins.',
  // gbdraw/render/groups/linear/pairwise_match.py: a link end without a usable feature ID or index.
  EMPTY_ENDPOINT: 'Check the comparison inputs: a comparison link does not identify the feature at one of its ends.',
  INDEX_ALIGNMENT: 'Check the comparison inputs: a comparison link does not identify the feature at one of its ends.',
  SOURCE_INDEX: 'Check the comparison inputs: a comparison link does not identify the feature at one of its ends.',
  SOURCE_VIEW_CONFLICT: 'Check the comparison inputs and display transforms.',
  RECORD_ID: 'A table row names another displayed record in this column. Swap the query and subject columns of the table, or assign the table to the pair it describes.',
  BLAST_TSV_REQUIRED: 'Choose a BLAST TSV for this pair, or set the pair to No comparison or Run LOSAT.',
  PAIR_TOPOLOGY: 'Each selected pair must join two different available records on adjacent rows, and appear once.',
  OUTFMT7_FIELDS: 'List all 12 standard BLAST outfmt 6 fields in the "# Fields:" line, or remove the line so that the first 12 columns are read in the standard order.',
  EMPTY_FILE: 'The file is empty.', NOT_GENBANK: 'The file is not a GenBank/DDBJ flat file: it has no record header line.',
  FASTA_IN_GENBANK: 'This file looks like FASTA. Use GFF3 + FASTA input, or a GenBank/DDBJ flat file.',
  // OV-1258, OV-1257: what a GenBank/DDBJ, GFF3, or FASTA slot file lacks.
  GZIP_INPUT: 'The file is gzip-compressed. Decompress it first (for example with gunzip); GenBank/DDBJ input must be an uncompressed flat file.',
  GENBANK_PARSE: 'The GenBank/DDBJ file could not be parsed: it may be cut off or damaged. Download or export it again.',
  SEQUENCE_UNDEFINED: 'Sequence content is undefined: a record has no sequence after ORIGIN. Use a GenBank/DDBJ file that includes the sequence.',
  GFF_NO_FEATURES: 'The GFF3 file has no feature lines.',
  FASTA_NO_ENTRIES: 'The FASTA file has no sequence entry: no line starts with ">".',
  FONT_FILE: 'Choose a readable TrueType or OpenType font file (.ttf, .otf, .ttc, .otc) that names its font family.',
  FONT_EMBEDDING: 'Turn off Embed font in exported SVG in Font to export without the font, or choose another font.',
  FONT_RASTER: 'The PNG cannot draw this font. Choose another font.',
  FORCED_LABEL: 'Open the feature\'s popup and set Label visibility to Default, or change the setting that prevents its label.',
  // LABEL_LAYOUT of labels not placed for lack of room (CAPACITY-STOP): no fix to choose; Python's counts follow.
  LABEL_CAPACITY: 'Label fewer features, use a smaller Label Font Size, or draw a Linear map.'
});
// What a decline names as changed under an evaluation of rule inputs (OV-447,
// OV-449): the label table import and the Specific rules import read it here.
// A key of its own (OV-1415): the label table's notice shows it with `tm`.
export const RULE_INPUT_CHANGE_WORDS = msg("colors, color rules, the Legend, or the diagram's inputs");
/**
 * Per code: the message, the actions offered, and the correction guidance when the context gives no reason.
 * @type {Readonly<Record<string, readonly [string, readonly string[], string?]>>}
 */
const DEFINITIONS = Object.freeze({
  SESSION_SIZE_LIMIT: ['The Session exceeds the browser size limit. Use a smaller Session file.', ['select-input', 'retry']],
  SESSION_BROWSER_UNSUPPORTED: ['This browser cannot complete this Session operation. Use a browser with the required gzip and Worker support.', ['retry']],
  SESSION_IMPORT_UNAVAILABLE: ['This browser cannot load Sessions because Session import Workers are unavailable. Use a browser with Worker support.', ['retry']],
  DECORATION_CONTINUITY: ['Could not preserve decoration placement. Reset the affected position or use Reset Layout on the previous Result, or restore matching settings, then Generate again.', ['edit-input', 'retry', 'save-session']],
  GENERATION_BUSY: ['A diagram generation request is already running. Wait for it to finish before retrying.', ['retry', 'save-session']],
  UNKNOWN: ['The operation failed without recognized diagnostic information. Retry; if it continues, save a Session for investigation.', ['retry', 'save-session']],
  VALIDATION_UNCLASSIFIED: ['Input validation failed. Review the inputs before retrying.', ['review-input', 'retry']],
  // The engine output depends only on the request, so Retry would fail the same way.
  RENDER_FAILED: ['The diagram engine failed while drawing this diagram. Generating again with the same inputs and settings fails the same way. Change the inputs or settings, or save a Session for investigation.', ['save-session']],
  // Python does not draw every forced label (an underlay feature, Embedded Only);
  // Generate fails the same way until the feature's popup choice changes.
  LABEL_NOT_DRAWN: ['Label visibility is On for a feature, but this diagram does not draw its label.', ['edit-input']],
  INPUT_INVALID: ['An input value is invalid.', ['edit-input', 'retry']],
  // OV-644: a rule commit whose Legend row would take the key of a listed row
  // of another color; committing again fails the same way.
  LEGEND_ROW_CONFLICT: ['A color rule would draw a Legend item that already exists in another color. Change the rule\'s Legend name or color, or rename that Legend item.', ['edit-input']],
  // OV-1206: the Specific rules owner (`runSpecificRuleWrite`) retried while
  // the rule inputs kept changing and gave up. Its result names no cause, so
  // one phrase covers a redraw too: the label table decline words a redraw
  // apart (OV-449), which here would need a wider port from the owner.
  RULE_INPUTS_CHANGED: [`The Specific rules file was not applied: ${RULE_INPUT_CHANGE_WORDS} changed while it was applied. Retry the table import.`, ['retry']],
  INPUT_REQUIRED: ['Supply GenBank input or matching GFF3 and FASTA inputs.', ['select-input', 'retry']],
  FASTA_REQUIRED: ['Supply a matching FASTA input for each GFF3 input.', ['select-input', 'retry']],
  INPUT_UNREADABLE: ['An input could not be read. Replace or reselect it.', ['select-input', 'retry']],
  SESSION_SAVE_REQUIRES_GENERATE: ['Save Session needs current feature metadata for the loaded Result. Generate once, then Save Session.', ['generate']],
  LIVE_EDIT_REQUIRES_GENERATE: ['Generate the diagram to update its label placement.', ['generate']],
  NO_RECORDS: ['No records were found. Choose input containing records.', ['select-input', 'retry']],
  RECORD_SELECTION: ['The record selection is invalid.', ['select-record', 'retry']],
  REGION_INVALID: ['The region is invalid.', ['edit-region', 'retry']],
  DEPTH_INVALID: ['The depth input or settings are invalid.', ['edit-depth', 'disable-track', 'retry']],
  TABLE_INVALID: ['The table is invalid.', ['edit-table', 'retry']],
  // The third element is the guidance when the producer names no reason.
  COMPARISON_INPUT: ['The comparison input is invalid.', ['edit-comparison', 'retry'], 'Supply a comparison sequence file (FASTA, GenBank, or DDBJ) or BLAST outfmt 6/7 as required.'],
  LOSAT_RUNTIME: ['The comparison search could not run or returned unusable output. Check the LOSAT or NCBI BLAST+ runtime, then Generate again.', ['edit-comparison', 'retry']],
  LOSAT_THREADING_UNAVAILABLE: ['Threaded LOSAT execution is unavailable in this browser environment. Select Serial or Auto execution, then Generate again.', ['edit-comparison', 'retry']],
  COMPARISON_IDENTITY: ['The features at the ends of comparison links do not match the displayed features. Review the comparison inputs and display transforms; save a Session if it continues.', ['edit-comparison', 'retry', 'save-session']],
  // A read-only saved comparison fails the same way until Inherit, Replace, or Clear is chosen.
  COMPARISON_CHOICE_REQUIRED: ['Choose how to handle the saved comparison in the Comparison panel, then Generate again.', ['edit-comparison']],
  // Replace with an empty draft would clear the saved comparison, which Clear comparison does explicitly.
  COMPARISON_REPLACEMENT_EMPTY: ['Set up a comparison in the Comparison panel, or choose Clear comparison, then Generate again.', ['edit-comparison']],
  ANNOTATION_TARGET: ['The region annotation target is invalid.', ['edit-annotation', 'retry']],
  TRACK_INVALID: ['The track settings are invalid.', ['edit-track', 'retry']],
  TRACK_LAYOUT: ['A circular track does not fit.', ['edit-track', 'retry']],
  // Labels of any kind overlap, or radial labels lack room; nothing new was drawn (LP-D43,
  // SCOPE-ALL-LABELS). Generate opens the fixes Python found.
  LABEL_LAYOUT: ['Labels do not fit.', ['edit-input', 'retry'], 'Generate again to choose a fix, or change the label settings.'],
  FEATURE_PLACEMENT: ['A Feature placement uses a lane that the current feature track does not have.', ['edit-track', 'retry']],
  FEATURE_IDENTITY: ['A feature edit does not identify a record of the current inputs. Generate again; if it continues, save a Session for investigation.', ['retry', 'save-session']],
  // The same request fails the same way until the setting or the mode changes.
  MODE_SETTING: ['A preserved session setting does not apply to this diagram mode.', ['edit-input']],
  // The Session's Reset alignment evidence is bound to the records and plan of the latest Align.
  ALIGNMENT_RESET_EVIDENCE: ['The Reset alignment evidence of the latest Align does not match the current records or alignment plan. Undo the changes made after Align, or select Align… again.', ['edit-input']],
  REGEX_SYNTAX: ['Not a valid regular expression.', ['edit-pattern', 'retry']],
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

/** @type {Readonly<Record<string, string>>} */
const FIELD_LABELS = Object.freeze({ protein_blastp_max_hits: 'Protein BLASTP Pairwise max hits',
  losatp_max_hits: 'Protein BLASTP Pairwise max hits', losatp_max_target_seqs: 'Protein BLASTP Max target seqs',
  losatp_member_max_hits: 'Protein BLASTP max hits per protein for groups',
  output_prefix: 'Output Prefix',
  arrow_head_length_ratio: 'Arrow head length ratio', arrow_shaft_width_ratio: 'Arrow shaft width ratio',
  keep_definition_left_aligned: 'Lock Definition Column', window: 'Window', step: 'Step',
  depth_window: 'Depth Window', depth_step: 'Depth Step', dinucleotide: 'Dinucleotide',
  evalue: 'E-value', bitscore: 'Bitscore', identity: 'Identity', alignment_length: 'Alignment length', wrap_length: 'Wrap length',
  plot_title_font_size: 'Plot Title Font Size', font_family: 'Font', feature_width_circular: 'Feature Width',
  depth_width_circular: 'Depth Width', gc_content_width_circular: 'GC Content Width',
  gc_content_radius_circular: 'GC Content Radius', gc_skew_width_circular: 'GC Skew Width',
  gc_skew_radius_circular: 'GC Skew Radius', match_height: 'Pairwise Match Height',
  conservation_ring_width: 'Comparison ring width', conservation_ring_gap: 'Comparison ring gap',
  center_reserved_radius: 'Center Reserved Radius', multi_record_min_radius_ratio: 'Minimum radius ratio',
  multi_record_column_gap_ratio: 'Column gap ratio', multi_record_row_gap_ratio: 'Row gap ratio' });
const CODES = new Set(Object.keys(DEFINITIONS));
// Python exception classes the adapter may name for RENDER_FAILED (class names, never messages).
const EXCEPTION_TYPES = new Set(`ValidationError ParseError ConfigError ExportError GbdrawError ZeroDivisionError
FloatingPointError OverflowError ArithmeticError AssertionError AttributeError IndexError KeyError LookupError
NotImplementedError RecursionError RuntimeError TypeError UnboundLocalError NameError ValueError Exception`.split(/\s+/));
// Locators shown in the summary; indexes are zero-based, ordinals one-based.
/** @type {Readonly<Record<string, string>>} */
const ORDINAL_LABELS = Object.freeze({ DECORATION_CONTINUITY: 'Result {n}.', COMPARISON_INPUT: 'Comparison sequence {n}.' });
// The tables a Session stores as text resources, as the documentation names them.
/** @type {Readonly<Record<string, string>>} */
export const SESSION_TABLE_LABELS = Object.freeze({ 'default-colors': 'Default colors', 'specific-colors': 'Specific colors',
  'label-whitelist': 'Label whitelist', 'qualifier-priority': 'Qualifier priority', 'feature-visibility': 'Feature visibility',
  'label-overrides': 'Label overrides' });
const CONFIG_PATH = /^[a-z][a-z0-9_]*(?:\.[a-z0-9_]+)+$/;
// OV-1274, OV-1276: the Circular track rows with an off control that a failed
// placement names (`keptInsideTracks`, `slotTrack`; gbdraw/diagrams/circular/radial_layout.py
// `SWITCHED_TRACK_IDS`), and the remedy: one sentence that names the control that turns it
// off and the setting that narrows it, with the number of choices it offers.
/** @typedef {{ name: string, remedy: string, options: number }} SwitchedTrack */
/** @type {Readonly<Record<string, SwitchedTrack>>} */
const SWITCHED_TRACKS = Object.freeze({
  gc_content: { name: 'GC content', remedy: 'Select Hide GC Content, or reduce GC Content Width (with Custom Track Slots on, the row width).', options: 2 },
  gc_skew: { name: 'GC skew', remedy: 'Select Hide GC Skew, or reduce GC Skew Width (with Custom Track Slots on, the row width).', options: 2 },
  ticks: { name: 'the scale', remedy: 'Clear Show Coordinate Scale (with Custom Track Slots on, remove the Ticks row).', options: 1 }
});
// The sentences of a summary other than the table rows above: each is one part
// (`ErrorTextPart`) and a catalog key; `{name}` takes the part's value.
const PHRASES = Object.freeze({
  diagramLinear: 'Diagram: Linear.', diagramCircular: 'Diagram: Circular.', sessionTable: 'Session table: {table}.',
  sequence: 'Sequence {n}.', pair: 'Pair: #{query} to #{subject}.', featureCaption: 'Feature: {caption}.',
  line: 'Line {row}.', column: 'Column {column}.', trackRow: 'Track row {row}.', depthSeries: 'Depth series {series}.',
  setting: 'Setting: {path}.', feature: 'Feature: {feature}.', featuresAffected: 'Features affected: {count}.',
  field: 'Field: {field}.', band: 'Available band: {inner}–{outer} px.',
  sessionTooNew: 'Session version {version} is newer than this gbdraw supports ({latest}).',
  sessionUnsupported: 'Unsupported session version: {version}. This gbdraw reads Session versions {versions}.',
  keptBelow: 'Below it, {px} px are kept for the tracks inside it.',
  keptBelowAmong: 'Below it, {px} px are kept for the tracks inside it, among them {tracks}.',
  // OV-1389 review item 6: each remedy is a whole sentence of `SWITCHED_TRACKS`, after one of these.
  makeRoom: 'To make room, do the following.', makeRoomAny: 'To make room, do one of the following.',
  // OV-2179: the room the labels not placed need (CAPACITY-STOP), the largest of each side, as Python sends it.
  outerLabelRoom: 'Outer labels need about {ratio} times the room around the map.',
  innerLabelRoom: 'Inner labels need about {ratio} times the room around the map.',
  and: '{first} and {second}',
  position: 'At character {position}.', requiredColumns: 'Required columns: {count}.', loadedRecords: 'Loaded records: {count}.',
  exception: 'Python exception: {type}.', applyAgain: 'Review the updated preview and Apply again.',
  schema: 'Load a supported Session file or recreate it with the current writer.',
  notASession: 'This file is not a gbdraw Session: it is not valid JSON.',
  // OV-2240: Save refuses a Session that Load could not open (the cap of services/session-file.js).
  sessionNotSaved: 'The Session was not saved: it is {size} MB expanded, and the browser loads at most {limit} MB.',
  sessionFileNotSaved: 'The Session was not saved: its file would be {size} MB, and the browser loads files of at most {limit} MB.',
  largestInputs: 'Largest inputs: {inputs}.',
  sessionSaveRemedy: 'Remove or replace the largest inputs with smaller ones, then save again.',
  liveEditFailed: 'Live edit failed: direct edits already applied are kept; geometry may still need updating.',
  liveEditRetry: 'Retry the live edit or use Generate.',
  liveEditChange: 'Change the edit, or change the settings and use Generate.',
  operationError: 'Operation error'
});
/**
 * One sentence or phrase of a user-facing error text. `text` is the English
 * template and its catalog key; each `{name}` takes `params[name]`: a part is
 * shown through its own key, a string or a number as written (identifiers,
 * file names, counts). A `verbatim` text is no key: it is shown as written
 * (Python words it, D-09).
 * @typedef {object} ErrorTextPart
 * @property {string} text
 * @property {Readonly<Record<string, string | number | ErrorTextPart>>} [params]
 * @property {true} [verbatim]
 */
const PLACEHOLDER = /\{([A-Za-z_][A-Za-z0-9_]*)\}/g;
/**
 * @param {ErrorTextPart} part
 * @param {(text: string) => string} translate
 * @returns {string}
 */
const renderPart = (part, translate) => {
  const text = part.verbatim ? part.text : translate(part.text);
  const { params } = part;
  return params ? text.replace(PLACEHOLDER, (whole, /** @type {string} */ name) => {
    if (!Object.hasOwn(params, name)) return whole;
    const value = params[name];
    return typeof value === 'object' ? renderPart(value, translate) : String(value);
  }) : text;
};
/**
 * The one renderer of error text parts: the English text with no `translate`,
 * the shown text with the UI translation (`errorText` in utils/error-text.js).
 * @param {readonly ErrorTextPart[]} parts
 * @param {(text: string) => string} [translate] the shown text of a key
 * @param {string} [separator] between two parts
 * @returns {string}
 */
export const renderErrorText = (parts, translate = (text) => text, separator = ' ') =>
  parts.map((part) => renderPart(part, translate)).join(separator);
/**
 * The leading whole parts whose English text fits in `limit` characters, so the
 * summary and the shown text are one cut; a first part longer than `limit` is
 * cut as English text, shown as written.
 * @param {readonly ErrorTextPart[]} parts
 * @param {number} limit
 * @returns {ErrorTextPart[]}
 */
const capParts = (parts, limit) => {
  /** @type {ErrorTextPart[]} */
  const kept = [];
  let length = -1;
  for (const part of parts) {
    length += 1 + renderPart(part, (text) => text).length;
    if (length > limit) break;
    kept.push(part);
  }
  if (kept.length || !parts.length) return kept;
  const cut = renderErrorText(parts).slice(0, limit);
  return cut ? [{ text: cut, verbatim: true }] : [];
};
/**
 * Joins parts with `and` as English does: "A and B and C".
 * @param {ErrorTextPart[]} items at least one
 * @returns {ErrorTextPart}
 */
const andPart = (items) => items.reduce((first, second) => ({ text: PHRASES.and, params: { first, second } }));
// A rendered feature ID (gbdraw/svg/ids.py) and an INSDC or GFF3 feature type.
const FEATURE_ID = /^[A-Za-z_][A-Za-z0-9_.-]{0,159}$/;
const FEATURE_TYPE = /^[A-Za-z0-9_'-]{1,40}$/;

/**
 * The context the normalizer keeps: only these keys survive `contextFor`.
 * @typedef {object} UserFacingContext
 * @property {string} [field]
 * @property {string} [reason]
 * @property {'python-character'} [positionUnit]
 * @property {string} [configPath]
 * @property {string} [exceptionType]
 * @property {string} [featureId]
 * @property {string} [featureType]
 * @property {string} [featureCaption]
 * @property {string} [sessionTable]
 * @property {'circular' | 'linear'} [diagramMode]
 * @property {number} [featureStart]
 * @property {number} [featureEnd]
 * @property {number} [featureCount]
 * @property {number} [position]
 * @property {number} [row]
 * @property {number} [column]
 * @property {number} [inputOrdinal]
 * @property {number} [recordIndex]
 * @property {number} [seriesIndex]
 * @property {number} [slotIndex]
 * @property {number} [recordCount]
 * @property {number} [columnCount]
 * @property {number} [codepoint]
 * @property {number} [innerPx]
 * @property {number} [outerPx]
 * @property {number} [keptInsidePx] The room kept below the band for the inside rows placed after it.
 * @property {string} [keptInsideTracks] The rows with an off control it is kept for (`SWITCHED_TRACKS` keys, comma-separated).
 * @property {string} [slotTrack] The failed row, when it has an off control (a `SWITCHED_TRACKS` key).
 * @property {number} [placementIndex]
 * @property {number} [queryRecordIndex]
 * @property {number} [subjectRecordIndex]
 * @property {string} [labelFitCounts] Python's sentence of the label fit counts, shown as sent.
 * @property {number} [labelCapacityOuter] The room the outer labels not placed need, as a multiple of the room around
 *   the map (CAPACITY-STOP), to one decimal place as Python rounds it.
 * @property {number} [labelCapacityInner] The same for the inner labels.
 * @property {number} [sessionVersion] The version of a Session file this app does not read.
 * @property {number} [latestSessionVersion] The newest Session version this app reads.
 * @property {string} [readSessionVersions] The Session versions this app reads, as words ("27–33, 39–42, and 44").
 * @property {number} [sessionBytes] The expanded size in bytes of a Session that Save refused.
 * @property {number} [sessionByteLimit] The size in bytes that Load opens at most.
 * @property {'file'} [sessionSizeKind] Set when the two sizes are those of the gzip file, not of its expanded JSON.
 * @property {string} [largestInputs] The file names of the largest inputs of a Session that Save refused.
 */

/**
 * The public diagnostic model: the wording and the facts of one failure.
 * @typedef {object} UserFacingError
 * @property {string} code
 * @property {string} operation
 * @property {string} stage
 * @property {UserFacingContext} context
 * @property {{ code: string, stage: string }[]} secondary
 * @property {string} summary The English text of `summaryParts`.
 * @property {ErrorTextPart[]} summaryParts The sentences of the summary, shown through `errorText`.
 * @property {{ label: string, text: string }[]} details
 * @property {string[]} actions
 */

/**
 * @typedef {Record<string, any>} DiagnosticContext Bounded facts about the failure; `contextFor` keeps only the known keys.
 * @typedef {{ code: string, stage: string, operation?: string, context?: DiagnosticContext }} NativeValidation
 * @typedef {Error & { code: string, stage: string, operation?: string, context: DiagnosticContext }} DiagnosticError
 */

/**
 * The single JS producer contract for a user-correctable failure: a code and
 * bounded context; wording stays here. The message is a fixed identifier, so an
 * uncaught error never carries document values.
 * @param {string} code
 * @param {DiagnosticContext} [context]
 * @param {{ stage?: string, operation?: string }} [options]
 * @returns {DiagnosticError}
 */
export const diagnosticError = (code, context = {}, { stage = 'request-validation', operation } = {}) =>
  Object.assign(new Error(`${code}/${context.reason || ''}`), { code, stage, ...(operation && { operation }), context });
// Exact native JS validation messages with no document interpolation.
/** @type {Map<string, NativeValidation>} */
const NATIVE_VALIDATIONS = new Map([
  ['The diagram engine returned incompatible feature metadata. Reload the page and Generate again.', { code: 'FEATURE_METADATA', stage: 'result-admission' }],
  ['A File-like object with arrayBuffer() or text() is required.', { code: 'RESOURCE_INVALID', stage: 'resource-staging' }],
  ['A File-like object with arrayBuffer() is required.', { code: 'RESOURCE_INVALID', stage: 'resource-staging' }],
]);
for (const [message, field, reason] of [
  ['Depth minimum must be less than or equal to depth maximum.', 'min_depth', 'ORDER'],
  ['Depth large tick interval must be greater than 0.', 'large_tick_interval', 'POSITIVE'],
  ['Depth small tick interval must be greater than 0.', 'small_tick_interval', 'POSITIVE'],
  ['Depth tick font size must be greater than 0.', 'tick_font_size', 'POSITIVE'],
  ['GC content minimum percent must be a finite number.', 'min_gc', 'FINITE'],
  ['GC content maximum percent must be a finite number.', 'max_gc', 'FINITE'],
  ['GC content minimum percent must be less than or equal to maximum percent.', 'min_gc', 'ORDER'],
  ['GC content large tick interval must be greater than 0.', 'large_tick_interval', 'POSITIVE'],
  ['GC content small tick interval must be greater than 0.', 'small_tick_interval', 'POSITIVE'],
  ['GC content tick font size must be greater than 0.', 'tick_font_size', 'POSITIVE']
]) NATIVE_VALIDATIONS.set(message, { code: 'INPUT_INVALID', stage: 'request-validation', context: { field, reason } });
for (const [message, code] of [
  ['Please upload a Depth TSV file or disable Show depth track.', 'DEPTH_INVALID'],
  ['Please upload at least one Depth TSV file or disable the depth track.', 'DEPTH_INVALID'],
  ['Please upload at least one BLAST outfmt 6/7 file for Pairwise Comparisons.', 'COMPARISON_INPUT'],
  ['Please upload at least one comparison sequence file for Pairwise Comparisons.', 'COMPARISON_INPUT'],
  ['No records found', 'NO_RECORDS'], ['No records found for circular conservation reference.', 'NO_RECORDS'],
  ['A resolved Linear comparison plan is required.', 'COMPARISON_INPUT'],
  ['The diagram engine returned an invalid Result list.', 'RESULT_INVALID'],
  ['The generated artifact has no selected preview Result.', 'RESULT_INVALID']
]) NATIVE_VALIDATIONS.set(message, { code, stage: code === 'RESULT_INVALID' ? 'result-admission' : 'request-validation' });
// Fixed labels/options emitted by the existing current-option validators.
for (const [label, field, reason] of [
  ['Protein BLASTP Pairwise max hits', 'protein_blastp_max_hits', 'POSITIVE_INTEGER'],
  ['Protein BLASTP max hits per protein for groups', 'orthogroup_member_max_hits', 'POSITIVE_INTEGER'],
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
  ['Protein BLASTP mode', 'protein_blastp_mode', 'pairwise, orthogroup, collinear']
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
NATIVE_VALIDATIONS.set('Canonical renderRequest records are required.',
  { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'RECORDS_REQUIRED' } });
for (const message of ['The preserved comparison is missing a required resource.',
  'The saved comparison is missing a required resource.',
  'The saved Circular comparison is missing a required resource.']) NATIVE_VALIDATIONS.set(message,
  { code: 'COMPARISON_INPUT', stage: 'request-validation', context: { field: 'comparison', reason: 'REQUIRED' } });
/**
 * @param {unknown} message
 * @returns {NativeValidation | null}
 */
const nativeValidation = (message) => {
  if (typeof message !== 'string') return null;
  for (const [template, reason] of /** @type {[RegExp, string][]} */ ([
    [/^Custom Track Slots use an obsolete schema\. Recreate the slots with schema version [0-9]+\.$/, 'TRACK_SCHEMA'],
    [/^Session contains unclassified top-level field\(s\): [\s\S]*$/, 'SESSION_FIELDS']
  ])) if (template.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason } };
  if (/^Invalid managed flag for (?:circular|linear)\.[a-z_]+\.$/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'config', reason: 'FIELDS' } };
  if (/^Missing canonical resource:/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'FIELDS' } };
  if (/^The SVG composition metadata is not valid JSON:/.test(message)) return { code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'JSON_FORMAT' } };
  // has() holds here, and every stored entry is a NativeValidation.
  if (NATIVE_VALIDATIONS.has(message)) return /** @type {NativeValidation} */ (NATIVE_VALIDATIONS.get(message));
  if (/^Invalid region (?:spec|coordinates) for LOSAT FASTA extraction: [\s\S]*$/.test(message)) {
    return { code: 'REGION_INVALID', stage: 'request-validation', context: { reason: 'REGION_FORMAT' } };
  }
  const pixel = /^(?:Circular|Linear) track slot '[\s\S]*' (height|spacing|inner_gap_px|outer_gap_px) must be (nonnegative|positive) finite number of pixels \(px optional\)\.$/.exec(message);
  if (pixel) return { code: 'TRACK_INVALID', stage: 'request-validation', context: {
    field: pixel[1], reason: pixel[2] === 'positive' ? 'PIXEL_POSITIVE' : 'PIXEL_NONNEGATIVE' } };
  const scalar = /^Circular track slot '[\s\S]*' (radius|width) must be a number greater than 0, in px or ×R\.$/.exec(message);
  if (scalar) return { code: 'TRACK_INVALID', stage: 'request-validation', context: { field: scalar[1], reason: 'POSITIVE_SCALAR' } };
  if (/^Circular track slot '[\s\S]*' uses obsolete field '(?:spacing|strict|compress|reserve)'\. Use inner_gap_px and outer_gap_px for physical gaps\.$/.test(message)) return { code: 'TRACK_INVALID', stage: 'request-validation', context: { reason: 'CIRCULAR_GAPS' } };
  const selector = /^Record selector #[0-9]+ is out of range \(loaded ([0-9]+) record\(s\)\)\.$/.exec(message);
  if (selector) return { code: 'RECORD_SELECTION', stage: 'request-validation', context: { reason: 'OUT_OF_RANGE', recordCount: Number(selector[1]) } };
  for (const [template, reason] of /** @type {[RegExp, string][]} */ ([
    [/^Record selector '[\s\S]*' did not match any record ID\.$/, 'NO_MATCH'],
    [/^Record selector '[\s\S]*' matched multiple records\. Use #index to disambiguate\.$/, 'AMBIGUOUS']
  ])) if (template.test(message)) return { code: 'RECORD_SELECTION', stage: 'request-validation', context: { reason } };
  const glyph = /^PDF fonts do not contain U\+([0-9A-F]{1,6})\. Use SVG to retain this text\.$/.exec(message);
  if (glyph) return { code: 'PDF_GLYPH', operation: 'export-pdf', stage: 'font-validation', context: { codepoint: parseInt(glyph[1], 16) } };
  const labelColumns = /^Invalid label TSV at line ([0-9]+): expected 5 columns, found ([0-9]+)\.$/.exec(message);
  if (labelColumns) return { code: 'TABLE_INVALID', stage: 'request-validation', context: { row: Number(labelColumns[1]), columnCount: 5, reason: 'FIELDS' } };
  const labelRequired = /^Invalid label TSV at line ([0-9]+): column ([1-4]) \((record_id|feature_type|qualifier|value)\) is required\.$/.exec(message);
  if (labelRequired) return { code: 'TABLE_INVALID', stage: 'request-validation', context: { row: Number(labelRequired[1]), field: labelRequired[3], reason: 'REQUIRED' } };
  const series = /^Depth series #([0-9]+) \(logical track index ([0-9]+)\) has no TSV source in any record\.(?: Add a TSV or remove the series\.)?$/.exec(message);
  if (series) return { code: 'DEPTH_INVALID', stage: 'request-validation', context: { seriesIndex: Number(series[2]), reason: 'DEPTH_SERIES_SOURCE' } };
  return null;
};
// The native track validator already identifies its failures. Project its
// first issue; never copy its message, slot ID, or arbitrary field path.
/** @type {Readonly<Record<string, readonly string[]>>} */
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
  annotation_anchor_unknown: ['OVERLAY_ANCHOR', 'anchor_slot'],
  annotation_anchor_ineligible: ['OVERLAY_ANCHOR', 'anchor_slot'], annotation_underlay_z: ['DRAW_ORDER', 'z'],
  annotation_to_slot_unknown: ['OVERLAY_SPAN', 'to_slot'], annotation_to_slot_ineligible: ['OVERLAY_SPAN', 'to_slot'],
  annotation_to_slot_same: ['OVERLAY_SPAN', 'to_slot'], annotation_to_slot_without_overlay: ['OVERLAY_SIDE', 'side'],
  annotation_foreground_z: ['DRAW_ORDER', 'z'], axis_side_conflict: ['TRACK_SIDE', 'side']
});
/**
 * @param {Record<string, any>} error
 * @returns {NativeValidation}
 */
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
/**
 * @template Fallback
 * @param {unknown} value
 * @param {ReadonlySet<string>} domain
 * @param {Fallback} fallback
 * @returns {string | Fallback}
 */
const identifier = (value, domain, fallback) => typeof value === 'string' && value.length <= 80 && domain.has(value) ? value : fallback;
/**
 * @param {DiagnosticContext | null | undefined} value The producer's context, unvalidated.
 * @returns {UserFacingContext}
 */
const contextFor = (value) => {
  /** @type {UserFacingContext} */
  const context = {};
  if (!value || typeof value !== 'object' || Array.isArray(value)) return context;
  if (FIELDS.has(value.field)) context.field = value.field;
  if (typeof value.reason === 'string' && Object.hasOwn(REASONS, value.reason)) context.reason = value.reason;
  if (value.positionUnit === 'python-character') context.positionUnit = value.positionUnit;
  if (typeof value.configPath === 'string' && value.configPath.length <= 80 && CONFIG_PATH.test(value.configPath)) context.configPath = value.configPath;
  if (typeof value.exceptionType === 'string' && EXCEPTION_TYPES.has(value.exceptionType)) context.exceptionType = value.exceptionType;
  if (typeof value.featureId === 'string' && FEATURE_ID.test(value.featureId)) context.featureId = value.featureId;
  if (typeof value.featureType === 'string' && FEATURE_TYPE.test(value.featureType)) context.featureType = value.featureType;
  // One-based feature span; genomes may exceed the other numeric bounds.
  for (const key of /** @type {Array<'featureStart' | 'featureEnd' | 'featureCount'>} */ (['featureStart', 'featureEnd', 'featureCount'])) {
    if (Number.isSafeInteger(value[key]) && value[key] >= 1 && value[key] <= 1e12) context[key] = value[key];
  }
  // The Web names the feature of a reported placement row (nameFeaturePlacementFailure).
  const caption = typeof value.featureCaption === 'string' ? value.featureCaption.replace(/[\s\u0000-\u001f\u007f]+/g, ' ').trim().slice(0, 80).trimEnd() : '';
  if (caption) context.featureCaption = caption;
  if (typeof value.sessionTable === 'string' && Object.hasOwn(SESSION_TABLE_LABELS, value.sessionTable)) context.sessionTable = value.sessionTable;
  if (Number.isSafeInteger(value.sessionVersion)) context.sessionVersion = value.sessionVersion;
  if (Number.isSafeInteger(value.latestSessionVersion)) context.latestSessionVersion = value.latestSessionVersion;
  if (typeof value.readSessionVersions === 'string' && /^[0-9\u2013, and]{1,200}$/.test(value.readSessionVersions)) {
    context.readSessionVersions = value.readSessionVersions;
  }
  // A refused Save: both sizes together, or neither; the names are one line of at most 200 characters.
  if (Number.isSafeInteger(value.sessionBytes) && Number.isSafeInteger(value.sessionByteLimit)
    && value.sessionBytes > 0 && value.sessionByteLimit > 0) {
    context.sessionBytes = value.sessionBytes;
    context.sessionByteLimit = value.sessionByteLimit;
  }
  if (context.sessionBytes !== undefined && value.sessionSizeKind === 'file') context.sessionSizeKind = 'file';
  const inputs = typeof value.largestInputs === 'string' ? value.largestInputs.replace(/[\s\u0000-\u001f\u007f]+/g, ' ').trim().slice(0, 200).trimEnd() : '';
  if (inputs) context.largestInputs = inputs;
  // The diagram mode whose Result the failure concerns, when it is not the one shown (E1).
  if (value.diagramMode === 'circular' || value.diagramMode === 'linear') context.diagramMode = value.diagramMode;
  // Python words the label fit counts (gbdraw/labels/fit.py): one line of printable text, shown as sent.
  if (typeof value.labelFitCounts === 'string' && value.labelFitCounts.length <= 80 && /^[\x20-\x7e]+$/.test(value.labelFitCounts)) {
    context.labelFitCounts = value.labelFitCounts;
  }
  for (const key of /** @type {const} */ (['labelCapacityOuter', 'labelCapacityInner'])) {
    const room = value[key];
    if (typeof room === 'number' && Number.isFinite(room) && room > 0 && room <= 1000000) context[key] = room;
  }
  if (typeof value.keptInsideTracks === 'string' && value.keptInsideTracks.length <= 80) {
    const tracks = [...new Set(value.keptInsideTracks.split(','))].filter((track) => Object.hasOwn(SWITCHED_TRACKS, track));
    if (tracks.length) context.keptInsideTracks = tracks.join(',');
  }
  if (typeof value.slotTrack === 'string' && Object.hasOwn(SWITCHED_TRACKS, value.slotTrack)) context.slotTrack = value.slotTrack;
  for (const key of /** @type {Array<'position' | 'row' | 'column' | 'inputOrdinal' | 'recordIndex' | 'seriesIndex' | 'slotIndex' | 'recordCount' | 'columnCount' | 'codepoint' | 'innerPx' | 'outerPx' | 'keptInsidePx' | 'placementIndex' | 'queryRecordIndex' | 'subjectRecordIndex'>} */ ([
    'position', 'row', 'column', 'inputOrdinal', 'recordIndex', 'seriesIndex', 'slotIndex', 'recordCount', 'columnCount', 'codepoint', 'innerPx', 'outerPx', 'keptInsidePx', 'placementIndex', 'queryRecordIndex', 'subjectRecordIndex'
  ])) {
    if (Number.isSafeInteger(value[key]) && value[key] >= 0 && value[key] <= (key === 'codepoint' ? 0x10ffff : 10000000)) {
      if (key !== 'position' || context.positionUnit === 'python-character') context[key] = value[key];
    }
  }
  return context;
};

// OV-1253: the CLI's words for a Session version this app does not read
// (gbdraw/session_io.py `validate_session`). The Session reader
// (services/session-authority.js) names the versions it reads.
/**
 * @param {UserFacingContext} context
 * @returns {ErrorTextPart | null}
 */
const sessionVersionGuidance = ({ reason, sessionVersion, latestSessionVersion, readSessionVersions }) => {
  if (sessionVersion === undefined) return null;
  if (reason === 'SESSION_TOO_NEW' && latestSessionVersion !== undefined) {
    return { text: PHRASES.sessionTooNew, params: { version: sessionVersion, latest: latestSessionVersion } };
  }
  if (reason === 'SESSION_VERSION_UNSUPPORTED' && readSessionVersions) {
    return { text: PHRASES.sessionUnsupported, params: { version: sessionVersion, versions: readSessionVersions } };
  }
  return null;
};

// OV-1274, OV-1276: the room a Circular placement keeps below its band for the
// inside rows after it (other rows may share it), and the off controls and widths
// of the failed row and the kept rows that have them.
/**
 * @param {UserFacingContext} context
 * @returns {ErrorTextPart[]}
 */
const trackRoomParts = ({ keptInsidePx, keptInsideTracks, slotTrack }) => {
  const kept = keptInsidePx !== undefined && keptInsideTracks ? keptInsideTracks.split(',') : [];
  /** @type {ErrorTextPart[]} */
  const parts = keptInsidePx === undefined ? [] : [kept.length
    ? { text: PHRASES.keptBelowAmong, params: { px: keptInsidePx, tracks: andPart(kept.map((track) => ({ text: SWITCHED_TRACKS[track].name }))) } }
    : { text: PHRASES.keptBelow, params: { px: keptInsidePx } }];
  const tracks = [...new Set([...(slotTrack ? [slotTrack] : []), ...kept])].map((track) => SWITCHED_TRACKS[track]);
  if (!tracks.length) return parts;
  const options = tracks.reduce((count, track) => count + track.options, 0);
  parts.push({ text: options > 1 ? PHRASES.makeRoomAny : PHRASES.makeRoom }, ...tracks.map((track) => ({ text: track.remedy })));
  return parts;
};

/** @typedef {{ operation?: string | null, stage?: string, code?: string, summaryLimit?: number, detailLimit?: number }} UserFacingOptions A null operation is read as no known operation. */

/**
 * The sole public diagnostic model. Raw messages, causes, logs and stacks are never copied.
 * @param {string | Record<string, any> | null | undefined} value
 * @param {UserFacingOptions} [options]
 * @returns {UserFacingError | null} Null only for a falsy value.
 */
export const normalizeUserFacingError = (value, {
  operation = 'unknown', stage = 'unknown', code = 'UNKNOWN', summaryLimit = 1000, detailLimit = 4000
} = {}) => {
  if (!value) return null;
  const object = value && typeof value === 'object' ? value : {};
  const source = /** @type {Record<string, any>} */ (object.code === 'CUSTOM_TRACK_PLAN_INVALID' ? typedTrackValidation(object) : CODES.has(object.code) ? object : nativeValidation(typeof value === 'string' ? value : object.message) || object);
  // `summary`, `details`, and `actions` are assigned below, before the model is returned.
  const result = /** @type {UserFacingError} */ ({
    code: identifier(source.code, CODES, identifier(code, CODES, 'UNKNOWN')),
    operation: identifier(operation === 'unknown' ? source.operation : operation, OPERATIONS, 'unknown'),
    stage: identifier(source.stage === 'unknown' ? stage : source.stage, STAGES, identifier(stage, STAGES, 'unknown')),
    context: contextFor(source.context),
    secondary: /** @type {{ code: string, stage: string }[]} */ ([])
  });
  for (const item of Array.isArray(source.secondary) ? source.secondary.slice(0, 2) : []) {
    if (item?.code !== 'CLEANUP_FAILED' || item?.stage !== 'cleanup') continue;
    result.secondary.push({ code: 'CLEANUP_FAILED', stage: 'cleanup' });
  }
  const [message, actions] = DEFINITIONS[result.code];
  const { inputOrdinal, row, column, slotIndex, seriesIndex, innerPx, outerPx, configPath, featureCaption, sessionTable, diagramMode } = result.context;
  const { featureId, featureType, featureStart, featureEnd, featureCount, field, reason, position, columnCount, recordCount } = result.context;
  const featureSpan = featureStart !== undefined && featureEnd !== undefined ? `${featureStart}..${featureEnd}` : '';
  const featureName = [[featureType, featureSpan].filter(Boolean).join(' '), featureId !== undefined ? `ID ${featureId}` : '']
    .filter(Boolean).join(', ');
  const { queryRecordIndex, subjectRecordIndex, labelFitCounts, labelCapacityOuter, labelCapacityInner, exceptionType } = result.context;
  const { sessionBytes, sessionByteLimit, sessionSizeKind, largestInputs } = result.context;
  // A Save over the cap: the file would not load again, and saving again fails the same way.
  const saveOverCap = result.code === 'SESSION_SIZE_LIMIT' && result.operation === 'session-save'
    && sessionBytes !== undefined && sessionByteLimit !== undefined;
  // One part per sentence, in the order the summary shows them.
  /** @type {(ErrorTextPart | null | false | '' | undefined)[]} */
  const pieces = [
    saveOverCap
      ? { text: sessionSizeKind === 'file' ? PHRASES.sessionFileNotSaved : PHRASES.sessionNotSaved, params: { size: (/** @type {number} */ (sessionBytes) / 1048576).toFixed(1), limit: Math.round(/** @type {number} */ (sessionByteLimit) / 1048576) } }
      : { text: message },
    saveOverCap && largestInputs !== undefined && { text: PHRASES.largestInputs, params: { inputs: largestInputs } },
    saveOverCap && { text: PHRASES.sessionSaveRemedy },
    diagramMode !== undefined && { text: diagramMode === 'linear' ? PHRASES.diagramLinear : PHRASES.diagramCircular },
    sessionTable !== undefined && { text: PHRASES.sessionTable, params: { table: { text: SESSION_TABLE_LABELS[sessionTable] } } },
    inputOrdinal !== undefined && { text: ORDINAL_LABELS[result.code] || PHRASES.sequence, params: { n: inputOrdinal } },
    // A comparison pair as the comparison panel numbers it (#1 to #2).
    queryRecordIndex !== undefined && subjectRecordIndex !== undefined
      && { text: PHRASES.pair, params: { query: queryRecordIndex + 1, subject: subjectRecordIndex + 1 } },
    featureCaption !== undefined && { text: PHRASES.featureCaption, params: { caption: featureCaption } },
    row !== undefined && { text: PHRASES.line, params: { row } },
    column !== undefined && { text: PHRASES.column, params: { column } },
    slotIndex !== undefined && { text: PHRASES.trackRow, params: { row: slotIndex + 1 } },
    seriesIndex !== undefined && { text: PHRASES.depthSeries, params: { series: seriesIndex + 1 } },
    configPath !== undefined && { text: PHRASES.setting, params: { path: configPath } },
    featureName && { text: PHRASES.feature, params: { feature: featureName } },
    // `undefined > 1` is false: with no count there is no line.
    /** @type {number} */ (featureCount) > 1 && { text: PHRASES.featuresAffected, params: { count: /** @type {number} */ (featureCount) } },
    field && { text: PHRASES.field, params: { field: FIELD_LABELS[field] ? { text: FIELD_LABELS[field] } : field } },
    sessionVersionGuidance(result.context)
      || (REASONS[/** @type {string} */ (reason)] || DEFINITIONS[result.code][2] ? { text: REASONS[/** @type {string} */ (reason)] || /** @type {string} */ (DEFINITIONS[result.code][2]) } : null),
    innerPx !== undefined && outerPx !== undefined && { text: PHRASES.band, params: { inner: innerPx, outer: outerPx } },
    ...trackRoomParts(result.context),
    labelFitCounts && { text: labelFitCounts, verbatim: /** @type {const} */ (true) },
    labelCapacityOuter !== undefined && { text: PHRASES.outerLabelRoom, params: { ratio: labelCapacityOuter } },
    labelCapacityInner !== undefined && { text: PHRASES.innerLabelRoom, params: { ratio: labelCapacityInner } },
    // The context keeps Python's zero-based index; people count from one.
    position !== undefined && { text: PHRASES.position, params: { position: position + 1 } },
    columnCount !== undefined && { text: PHRASES.requiredColumns, params: { count: columnCount } },
    recordCount !== undefined && { text: PHRASES.loadedRecords, params: { count: recordCount } },
    exceptionType && { text: PHRASES.exception, params: { type: exceptionType } },
    result.operation === 'align' && field === 'direction' && reason === 'SOURCE_VIEW_CONFLICT' && { text: PHRASES.applyAgain },
    result.code === 'INPUT_INVALID' && field === 'schema' && { text: PHRASES.schema }
  ];
  // A file chosen in Load Session that is not JSON is not a Session;
  // every other Session reason keeps its correction.
  const notASession = result.operation === 'session-load' && result.code === 'INPUT_INVALID'
    && field === 'schema' && reason === 'JSON_FORMAT';
  result.summaryParts = capParts(notASession ? [{ text: PHRASES.notASession }]
    : pieces.filter(/** @returns {piece is ErrorTextPart} */ (piece) => Boolean(piece)),
  Number.isSafeInteger(summaryLimit) ? Math.max(0, Math.min(summaryLimit, 1000)) : 1000);
  result.summary = renderErrorText(result.summaryParts);
  const detail = [`Code: ${result.code}`, `Operation: ${result.operation}`, `Stage: ${result.stage}`,
    ...Object.entries(result.context).map(([key, item]) => `${key}: ${item}`),
    ...result.secondary.map((item) => `Secondary: ${item.code} / ${item.stage}`)].join('\n');
  result.details = [{ label: 'Diagnostics', text: detail.slice(0,
    Number.isSafeInteger(detailLimit) ? Math.max(0, Math.min(detailLimit, 4000)) : 4000) }].slice(0, 8);
  // A table stored in a Session is corrected in the Session file, not in this page.
  result.actions = sessionTable !== undefined || saveOverCap ? ['select-input'] : [...actions];
  return result;
};


// A value that a `catch` block receives is a failure even when it is falsy
// (`throw undefined`), so it normalizes to the unknown failure instead of to
// "no error" (OV-70).
/**
 * @param {unknown} value
 * @param {Parameters<typeof normalizeUserFacingError>[1]} [options]
 * @returns {UserFacingError}
 */
export const normalizeCaughtError = (value, options) => (
  // A truthy value always normalizes to a model.
  /** @type {UserFacingError} */ (normalizeUserFacingError(value || {}, options))
);

// A failure whose own recovery changes no input (B11, R6): the runtime failed,
// not the request, so the same request can succeed on Retry.
const PLAIN_RETRY_ACTIONS = new Set(['retry', 'reload', 'save-session']);
/** @param {UserFacingError} error */
export const retryCanSucceed = (error) =>
  error.actions.includes('retry') && error.actions.every((action) => PLAIN_RETRY_ACTIONS.has(action));

// A live edit rerenders the committed Session with the current editor tables, so
// repeating the edit sends the same request. Retry is offered only when it can
// succeed; otherwise the edit or the settings must change, and settings apply
// on Generate.
/** @param {UserFacingError | null | undefined} error */
export const liveEditFailure = (error) => {
  if (!error) return null;
  const retry = retryCanSucceed(error);
  const next = error.actions.includes('generate') ? [] : [{ text: retry ? PHRASES.liveEditRetry : PHRASES.liveEditChange }];
  const noteParts = [{ text: PHRASES.liveEditFailed }, ...error.summaryParts, ...next];
  return { ...error, actions: retry ? ['retry', 'generate'] : ['generate'], noteParts, note: renderErrorText(noteParts) };
};

// Presentation of finite operation and transaction facts stays with public wording.
/** @type {Readonly<Record<string, string>>} */
const OPERATION_TITLES = Object.freeze({
  generate: 'Generation Error', align: 'Alignment error', evaluateRules: 'Rule error', 'session-load': 'Session not loaded',
  'export-svg': 'SVG export error', 'export-png': 'PNG export error', 'export-pdf': 'PDF export error',
  listSequenceRecords: 'Input error', listGffFastaRecords: 'Input error', readComparisonSequence: 'Input error',
  readFeatureOverrideTable: 'Feature edits TSV error',
  'feature-extraction': 'Feature preparation error'
});
/** @type {Readonly<Record<string, string>>} */
const RECOVERY_GUIDANCE = Object.freeze({
  'no-result': 'No successful Result is available yet. Correct the cause and retry.',
  // The one notice that the shown Result is the last successful one (PD-OI-016).
  preserved: 'Showing the last successful Result; your changes were not applied.',
  restored: 'Showing the last successful Result again; your changes were not applied. Correct the cause and retry.',
  'restore-failed': 'Restoration could not complete. Save a Session if available, then reload before retrying.'
});
/** @param {string | null | undefined} operation */
export const operationErrorTitle = (operation) => (
  Object.hasOwn(OPERATION_TITLES, String(operation)) ? OPERATION_TITLES[String(operation)] : PHRASES.operationError);
/** @param {string | null | undefined} recovery */
export const generationRecoveryGuidance = (recovery) => (
  Object.hasOwn(RECOVERY_GUIDANCE, String(recovery)) ? RECOVERY_GUIDANCE[String(recovery)] : '');
