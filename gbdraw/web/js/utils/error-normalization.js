// @ts-check
// This module imports only leaves without imports of their own (the mode
// profiles, the persisted contracts and utils/i18n.js), so a test that copies
// it needs only those beside it.
import { MODE_PROFILE_DATA } from '../mode-profiles.generated.js';
import {
  INPUT_NAME_CHARACTERS, SESSION_CONTAINER_DAMAGED_MESSAGE, SESSION_CONTAINER_NO_DOCUMENT_MESSAGE
} from '../contracts/schema.generated.js';
import { msg, numberText } from './i18n.js';
const OPERATIONS = new Set(['unknown', 'generate', 'align', 'feature-extraction', 'export-svg', 'export-png', 'export-pdf', 'evaluateRules', 'embedFonts', 'pdfFonts', 'buildProteinLosatCacheKeys', 'checkExportLimit', 'computeLabelFitRow', 'convertRingTablesToSearchFrame', 'convertStoredProteinResultFrames',
  'addRegionAnnotation', 'convertLosatpPairsToGenomicPayload', 'decideAutoOrientation', 'deriveRegionDrawing', 'extractCdsProteinFasta', 'extractFirstFasta', 'findInOtherFiles',
  'hydrateProteinLosatTsv', 'listGffFastaRecords', 'listSequenceRecords', 'planCardWindows', 'planGenomeScales', 'planLinearUnits',
  'planRecordPlaces', 'planRingScales',
  'promoteLegacyLosatpCache', 'readAllSequenceRecords', 'readComparisonSequence', 'readEmbedSequences', 'readFeatureOverrideTable', 'readFeatureSequences', 'readMatchSequences', 'readRecordSpanSequences', 'readSearchQueries', 'registerFonts', 'resolveFeatureAnchor', 'resolveLegacyProteinReferences',
  'resolveSimilarityAlignment', 'searchHitAnnotationSet', 'synthesisQualifiers', 'synthesizeColorRules', 'upgradeLegacySession', 'validateConfigOverrides', 'session-save', 'session-load']);
const STAGES = new Set(['unknown', 'initialization', 'resource-staging', 'request-validation',
  'helper', 'rule-validation', 'losat', 'render', 'result-admission', 'cleanup', 'export-capture', 'export-conversion', 'font-validation', 'transport', 'read', 'write', 'parse', 'drawing-list']);
const FIELDS = new Set(`legend title scale decorations pattern record_selector region start end sourceStart sourceEnd recordLength
recordIndex queryIndex subjectIndex depth min_depth max_depth window step tick font_size font_family plot_title_font_size height large_tick_interval small_tick_interval tick_font_size
inner_gap_px outer_gap_px radius width spacing arrow_head_length_ratio arrow_shaft_width_ratio keep_definition_left_aligned color action feature_type qualifier value record_id label_text
config configOverrides records anchors schema recordKey groupId direction sourceStrand role blast files
input comparison protein_blastp_max_hits orthogroup_member_max_hits bitscore evalue identity
losatp_max_hits losatp_max_target_seqs losatp_member_max_hits losatp_mode losatn_task record_gencodes
losat losat_search program pairs source losat_gencode comparison_sequence comparison_sequence_source comparison_fasta
conservation_blast_files conservation_dataframes conservation_reference conservation_sequence_files
conservation_labels conservation_colors conservation_losat_gencodes conservation_search_results reference_gencode
alignment_length wrap_length anchor_feature collinear_min_anchors collinear_max_gene_gap collinear_block_merge_gap
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
  // A threaded run whose host could not run LOSAT's threads (services/losat-runtime.js `losatThreadsFailure`, OV-3316).
  THREADS_FAILED: 'The threaded LOSAT engine could not start or keep its threads running in this browser.',
  DECORATION_IDENTITY: 'Source, region, mode, grouping or record identity changed, is unknown, or is ambiguous.',
  DECORATION_TARGET: 'The decoration target is missing or ambiguous.',
  DECORATION_METADATA: 'Composition metadata is invalid or unavailable.',
  ANNOTATION_SET: 'Select an existing annotation set.', ANNOTATION_MARKS: 'Use supported annotation mark filters.',
  OVERFLOW: 'Choose Error on overflow, Compress, or Clip.', LAYER: 'Choose Foreground or Underlay.',
  OVERLAY_ANCHOR: 'Select an eligible enabled anchor for an annotation track whose placement is Overlay anchor (Circular) or On axis (Linear).',
  OVERLAY_SPAN: 'For Extend to, select an eligible enabled slot other than the anchor, or None.',
  OVERLAY_SIDE: 'Choose the placement Overlay anchor (Circular) or On axis (Linear) when an anchor or an Extend to slot is selected.',
  DRAW_ORDER: 'Place underlay annotations below the anchor and the Extend to slot, and foreground annotations above them.',
  TRACK_RENDERER: 'Choose a renderer supported by the current diagram mode.', TRACK_SIDE: 'Choose a side supported by the current diagram mode and renderer.',
  TRACK_PARAMS: 'Use the parameters supported by the selected renderer.', TRACK_AXIS: 'Select one eligible enabled feature track as the axis.',
  TRACK_LANE: 'Use a supported feature lane and matching side.', CONSERVATION_SOURCE: 'Use the supported conservation workflow and select a comparison source.',
  TICK_PARAMETER: 'Use large_tick_interval instead of tick_interval.', FEATURES_COUNT: 'Use the supported number of feature tracks for this layout.',
  UNDERLAY_FEATURE_SLOTS: 'Feature types drawn as Underlay (repeat_region by default) need exactly one enabled Features row. Keep one Features row enabled, or set those types to Arrow or Rectangle in Feature Types.',

  CIRCULAR_MULTI_RECORD_SIZE_MODE: 'Choose Auto (Default), Linear, or Equal.',
  // One circle draws neither of these yet; the words are the command's (`gbdraw circular --one_circle`).
  ONE_CIRCLE_DEPTH: 'Depth tracks are not drawn when several records share one circle.',
  ONE_CIRCLE_CONSERVATION: 'Conservation rings are not drawn when several records share one circle.',
  ONE_CIRCLE_REGIONS: 'Region tracks are not drawn when several records share one circle.',
  LINEAR_TRACK_LAYOUT: 'Choose Above, Features on axis, or Below.',
  LINEAR_LABEL_PLACEMENT: 'Choose Auto or Above Feature.',
  PROTEIN_BLASTP_MODE: 'Choose Pairwise matches, Similarity groups, or Collinear blocks.',
  COLLINEAR_SEARCH_SCOPE: 'Choose Adjacent pairs or All records.',
  COLLINEAR_UNIT_MODE: 'Choose Auto, CDS, or Locus.',
  COLLINEAR_ANCHOR_MODE: 'Choose All hits or Reciprocal best hits.',
  COLLINEAR_MERGE_ORIENTATION: 'Choose Strand or Gene order.',
  COLLINEAR_COLOR_MODE: 'Choose Average identity, Orientation, or Orientation + identity.',
  ORTHOGROUP_MEMBERSHIP_MODE: 'Choose anchor_core_v1.',

  UNKNOWN_CONFIG_PATH: 'Remove unknown configuration overrides or use a supported setting.',
  CIRCULAR_SETTING: 'It applies only to Circular diagrams. Reset it under Preserved session settings, or switch to Circular.',
  // A record cut by Wrap Rows (gbdraw/layout/linear_wrap.py): what rows cannot draw, and the ways out. The Wrap Rows
  // field shows the same text while it is on and cannot be changed (`reasonText`), except WRAP_TRACK_LAYOUT, whose
  // Track Layout step is spent then (`WRAPPED_BLOCKED_TEXT` in app/linear-wrap-field.js).
  WRAP: 'A wrapped figure draws no comparisons. Remove the comparison inputs (BLAST, comparison pairs, LOSAT, Similarity groups), or turn Wrap Rows off.',
  WRAP_TRACK_LAYOUT: 'Wrapped rows draw their tracks above or below the axis. Set Track Layout to Above or Below, move the feature track off the axis in Custom Track Slots, or turn Wrap Rows off.',
  WRAP_SCALE_STYLE: 'Each wrapped row draws a ruler on its axis. Set Scale Style to Ruler (Ticks), or turn Wrap Rows off.',
  WRAP_NORMALIZE_LENGTH: 'Wrapped rows draw at one scale. Turn Normalize Record Lengths off, or turn Wrap Rows off.',
  WRAP_GRID_ROW: 'A wrapped record needs a Record Layout row of its own. Give each record its own row, or turn Wrap Rows off.',
  // Align at feature's Resolve and gbdraw linear --anchor_feature (gbdraw/api/record_planning.py, gbdraw/layout/linear_wrap.py).
  // Python raises it for a cut record the anchor moves sideways and for records lined up with a cut reference the
  // anchor does not rotate (`refuse_wrapped_anchor_moves`); the one text is true for both.
  WRAP_ANCHOR_FEATURE: 'Align at feature cannot line up these records while Wrap Rows cuts a record into rows: it would move a cut record sideways, or line the others up with the feature of a cut record that it does not rotate. Only a complete circular record is aligned by rotation and can be wrapped. Turn Wrap Rows off, then use Align at feature.',
  // A cut record whose stored translation moves it sideways, whatever wrote it, or a Similarity alignment without
  // its comparison inputs (gbdraw/layout/linear_wrap.py). Web text names Web-visible sources only (OV-2472).
  WRAP_RECORD_POSITION: 'A wrapped record cannot keep a sideways move (from Align at feature, Similarity groups, or a loaded Session). Click Reset record positions under Wrap Rows, or turn Wrap Rows off.',
  // A PNG or PDF over its format's size limit (gbdraw/render/export_limits.py). The limit and the advice are the
  // EXPORT_LIMIT phrases: a row length or DPI is advised only when Python found one that fits (D-46).
  PNG_SIZE: 'Use SVG, which has no size limit.',
  PDF_SIZE: 'Use SVG, which has no size limit.',
  // A requested figure size the figure cannot fit (gbdraw/layout/figure_size.py). With Python's lengths in the
  // context, the FIGURE_SIZE phrases word the parts and the smallest size; this is the text without them.
  FIGURE_TOO_NARROW: 'The figure cannot be this narrow. Set a larger Figure width, or leave it empty (Auto).',
  FIGURE_TOO_SHORT: 'The figure cannot be this short. Set a larger Figure height, or leave it empty (Auto).',
  SELECTOR_NO_MATCH: 'No record has a feature that this selector names. Check the key and the value, for example gene=COX1.',
  SELECTOR_AMBIGUOUS: 'In one record, this selector names features that start at different bases. Use a selector that names one feature, such as locus_tag=<tag>.',
  // Align at feature reads the drawn request, so new Record Layout rows count after a Generate (OV-2501).
  ANCHOR_SHARED_ROW: 'Records that share a row cannot be aligned at a feature. Give each record its own row in Record Layout, click Generate Diagram, then use Align at feature again.',
  ANCHOR_STALE: 'The diagram changed after the preview. Align at feature again.',
  // Align at feature writes a rotated record's start to its File card's record (OV-2474). Generate splits a card that
  // draws every record of its file into one card per record (D-98).
  ANCHOR_RECORD_CARDS: 'Align at feature writes each record\'s start on its own File card, and these records share one File card (a file with several records). Click Generate Diagram, which gives each record its own File card, then use Align at feature again.',
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
  ADJACENT_ALL: 'Choose Adjacent pairs or All records.', BLASTP_MODE: 'Choose Pairwise matches, Similarity groups, or Collinear blocks.',
  UNIQUE_IDS: 'Use distinct record or protein identifiers.', MATCH_IDS: 'Use matching protein FASTA and metadata identifiers.',
  SEARCH_FRAME: 'Use comparison coordinates within the selected and cropped record, counted on its source strand.',
  LOSAT_OPTION_PROGRAM: 'Use only the options of the selected LOSAT program.',
  // A TLOSATX genetic code that is not an NCBI table LOSAT reads (gbdraw/api/options.py, services/linear-sources.js
  // `requireLosatGeneticCode`); `geneticCodes` lists the tables (PHRASES.geneticCodes).
  GENETIC_CODE: 'Use the number of an NCBI genetic code that LOSAT reads.',
  LOSAT_PLAN: 'Choose record pairs that the LOSAT search can run: two or more records, LOSAT rows only with a LOSAT program, and LOSATP Pairwise for selected pairs.',
  // A Selected plan with a LOSATP mode that compares every adjacent pair.
  SELECTED_PAIRS_NEED_PAIRWISE: 'Similarity groups and Collinear blocks compare every adjacent pair with LOSAT, so they do not run on pairs set one by one in Selected pairs. Choose Run LOSAT for all adjacent pairs, or choose Pairwise matches.',
  LOSAT_TASK: 'Choose a LOSATN task that the selected LOSAT runtime supports.',
  // The serial LOSAT engine has no task dc-megablast (services/losat.js, OV-2376).
  LOSAT_TASK_THREADED: 'LOSATN task dc-megablast runs only with threaded LOSAT. Select Threaded or Auto under Execution in a browser where Threaded is available, or choose task megablast or blastn.',
  // The LOSATP converter named another derived identity than Generate looked up (run-analysis.js).
  DERIVED_KEY_MISMATCH: 'The converted LOSATP comparison does not match the settings it was looked up with. If Retry fails the same way, save the Session and report the problem.',
  RING_LOSAT_PROGRAM: 'Choose LOSATN or TLOSATX for comparison rings.',
  RING_LOSAT_INPUT: 'Give each ring one or more comparison sequence files, without BLAST tables, and use the reference file as the subject.',
  SEQUENCE_MISSING: 'Use a comparison file that contains the sequence (FASTA, or GenBank or DDBJ with an ORIGIN sequence).',
  UNAVAILABLE: 'Install LOSAT or NCBI BLAST+, or select a runtime executable.',
  FAILED: 'The search runtime exited with an error; check its message and inputs.',
  // A trap of LOSAT's code (services/losat-runtime.js `losatTrapFailure`, OV-3327), in Sequence search and Find in
  // other files; Generate says it in its own first sentence (`losatTrapped`).
  TRAP: 'LOSAT stopped unexpectedly.',
  // A trap left after the threaded pool's retries (services/losat-runtime.js `losatThreadedTrapFailure`), in the same
  // two searches; Generate says it in its own first sentence (`losatThreadedTrapped`).
  THREADED_TRAP: 'LOSAT stopped unexpectedly in threaded execution; Serial or Auto execution may run the search.',
  OUTPUT: 'The search output was not BLAST outfmt 6 rows for the searched records.',
  JSON_FORMAT: 'Use valid JSON.', VISIBILITY_ACTION: 'Use show, off, or exclude_matching; accepted aliases are on, hide, false, and 0.',
  TARGET_RECORD: 'Choose an available target record.',
  // A drawing's name (app/drawing-list.js; Python's _identity_problem).
  BLANK_NAME: 'Enter a name for the drawing.',
  REPEATED_NAME: 'Another drawing has this name. Names that differ only in letter case are the same name.',
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
  ORDER: 'The start or minimum must not exceed the end or maximum.', INCREASING: 'Use a minimum below the maximum.', RECORD_BOUNDS: 'Keep the region within the record length.',
  START_AFTER_END: 'Choose a start before the end.',
  // A Linear card window with Across the origin ticked (PD-OI-124, services/session-request.js linearWindowsPayload).
  ORIGIN_ORDER: 'Across the origin is ticked: give a Start greater than the End, or clear the box.',
  DISPLAY_START_BOUNDS: 'Use a display start between 1 and the record length.', CROP_START_CONFLICT: 'Choose a crop or an explicit display start.',
  REFERENCE_REQUIRED: 'Fill in the reference name on every row of the Depth TSV.', REFERENCE_MISMATCH: 'Match depth references to the selected record.',
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
  // A Session file container (gbdraw/session_container.py, services/session-container.js) that Load refuses
  // whole: the sentence is the refusal (NOT_A_SESSION).
  CONTAINER_DAMAGED: SESSION_CONTAINER_DAMAGED_MESSAGE,
  CONTAINER_NO_DOCUMENT: SESSION_CONTAINER_NO_DOCUMENT_MESSAGE,
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
  // A BLAST TSV gives the coordinates of one drawn record (gbdraw/api/record_planning.py refuse_file_comparison_of_windows).
  // The browser search reads a crop as one span of its record (gbdraw/web/js/services/linear-units.js, OV-2739).
  WINDOWED_FILE_COMPARISON: 'A record of this pair is drawn in several windows, and a BLAST TSV gives the coordinates of one drawn record. Set the pair to Run LOSAT.',
  // Generate's Inherit of a saved comparison (services/imported-comparison-intent.js): its indexes count the units
  // the saved request drew, and this drawing draws other records or windows, or holds no saved comparison.
  INHERIT_UNITS_CHANGED: 'The drawn records or windows changed since the comparison was saved. Choose Replace or Clear.',
  INHERIT_SOURCE_MISSING: 'This drawing holds no saved comparison to inherit. Choose Replace or Clear.',
  PAIR_TOPOLOGY: 'Each selected pair must join two different available records on adjacent rows, and appear once.',
  OUTFMT7_FIELDS: 'List all 12 standard BLAST outfmt 6 fields in the "# Fields:" line, or remove the line so that the first 12 columns are read in the standard order.',
  EMPTY_FILE: 'The file is empty.', NOT_GENBANK: 'The file is not a GenBank/DDBJ flat file: it has no record header line.',
  FASTA_IN_GENBANK: 'This file looks like FASTA. Choose FASTA (or GFF3 + FASTA when it has a GFF3 file) as the input type, or use a GenBank/DDBJ flat file.',
  // OV-1258, OV-1257: what a GenBank/DDBJ, GFF3, or FASTA slot file lacks.
  GENBANK_PARSE: 'The GenBank/DDBJ file could not be parsed: it may be cut off or damaged. Download or export it again.',
  // A gzip input that expands past its cap (gbdraw/io/input_text.py MAX_EXPANDED_INPUT_BYTES); the cap is PHRASES.gzipLimit.
  GZIP_EXPANDED_SIZE: 'Decompress the gzip file and choose the plain file.',
  // A staging of the inputs for the diagram Worker that found no memory for their bytes (a RangeError in the page,
  // a failed read of a readable input in the Worker; OV-2749). It states the whole failure, so the code's message,
  // which asks to reselect the input, is not shown.
  INPUT_MEMORY: 'There was not enough memory to prepare the inputs for drawing. Close other tabs and retry, or save a Session and reload the page before retrying.',
  // A drawn diagram whose Worker reply the page cannot receive (a `messageerror` of the run in
  // services/diagram-generation.js, PD-OI-118). It states the whole failure, so the code's message is not shown.
  REPLY_TOO_LARGE: 'The diagram was drawn, but it is too large for this browser page to receive. Draw a smaller region, or turn off tracks or labels you do not need, or draw it with the gbdraw command line (Save Session, then `gbdraw circular --session` or `gbdraw linear --session`).',
  SEQUENCE_UNDEFINED: 'Sequence content is undefined: a record has no sequence after ORIGIN. Use a GenBank/DDBJ file that includes the sequence.',
  GFF_NO_FEATURES: 'The GFF3 file has no feature lines.',
  FASTA_NO_ENTRIES: 'The FASTA file has no sequence entry: no line starts with ">".',
  // gbdraw/io/genome.py require_nucleotide_records: a FASTA given without GFF3 holds nucleotide records only.
  NOT_NUCLEOTIDE: 'This FASTA file holds a sequence that is not nucleotide, such as a protein; gbdraw draws nucleotide records. Use a FASTA file of DNA or RNA sequences (IUPAC codes A C G T U R Y S W K M B D H V N).',
  // gbdraw/io/genome.py load_query_fasta: a Sequence search query is a nucleotide or a protein sequence.
  NOT_SEQUENCE: 'A query holds a character that is neither a nucleotide nor an amino acid code. Give DNA, RNA or protein sequences in FASTA: letters, and * for a stop.',
  EMPTY_SEQUENCE: 'A record of this FASTA file has no sequence. Give each FASTA record its sequence.',
  FONT_FILE: 'Choose a readable TrueType or OpenType font file (.ttf, .otf, .ttc, .otc) that names its font family.',
  FONT_EMBEDDING: 'Turn off Embed font in exported SVG in Font to export without the font, or choose another font.',
  FONT_RASTER: 'The PNG cannot draw this font. Choose another font.',
  FORCED_LABEL: 'Open the feature\'s popup and set Label visibility to Default, or change the setting that prevents its label.',
  // LABEL_LAYOUT of labels not placed for lack of room (CAPACITY-STOP): no fix to choose; Python's counts follow.
  LABEL_CAPACITY: 'Label fewer features, use a smaller Label Font Size, or draw a Linear map.'
});
// The actions of a reason in place of its code's, where Retry would fail the same way: a refused gzip file is
// replaced by its plain file, and a reply too large for the page is drawn from a saved Session (PD-OI-118).
/** @type {Readonly<Record<string, readonly string[]>>} */
const REASON_ACTIONS = Object.freeze({ GZIP_EXPANDED_SIZE: ['select-input'], REPLY_TOO_LARGE: ['save-session'],
  INPUT_MEMORY: ['retry', 'save-session', 'reload'] });
// The reasons that state the whole failure: the code's message is not shown with them.
const WHOLE_FAILURE_REASONS = new Set(['REPLY_TOO_LARGE', 'INPUT_MEMORY']);
// What a decline names as changed under an evaluation of rule inputs (OV-447,
// OV-449): the label table import and the Specific rules import read it here.
// A key of its own (OV-1415): the label table's notice shows it with `tm`.
export const RULE_INPUT_CHANGE_WORDS = msg("colors, color rules, the Legend, or the diagram's inputs");
// The code of an engine that ran out of memory (its DEFINITIONS entry below): the service that owns the Worker
// replaces that engine when a reply carries it.
export const OUT_OF_MEMORY_CODE = 'OUT_OF_MEMORY';
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
  // OUT_OF_MEMORY_CODE: Python's MEMORY_REFUSAL_CODE (gbdraw/web_support/error_adapter.py), a MemoryError by its
  // class. The engine's memory never shrinks, so services/diagram-generation.js replaces that engine and Retry runs
  // on a new one (OV-3027, PD-OI-125).
  [OUT_OF_MEMORY_CODE]: ['The diagram engine ran out of memory and was replaced with a new one. Generate again. If it runs out of memory again, draw a smaller region, turn off tracks or labels you do not need, or draw it with the gbdraw command line (Save Session, then `gbdraw circular --session` or `gbdraw linear --session`).', ['retry', 'save-session']],
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
  ALIGNMENT_RESET_EVIDENCE: ['What Reset alignment would restore no longer matches the current records or alignment plan. Undo the changes made after Align, or select Align… again.', ['edit-input']],
  REGEX_SYNTAX: ['Not a valid regular expression.', ['edit-pattern', 'retry']],
  RESOURCE_INVALID: ['Input resource preparation failed. Reselect the input and retry.', ['select-input', 'retry']],
  HELPER_PROTOCOL: ['The helper request is invalid. Retry the operation.', ['retry', 'save-session']],
  RUNTIME_INCOMPATIBLE: ['The diagram engine is incompatible with this Web app. Reload the page and retry; contact the site administrator if it continues.', ['reload', 'retry']],
  WORKER_INIT: ['The diagram runtime could not start. Save a Session and reload before retrying.', ['save-session', 'reload', 'retry']],
  FEATURE_METADATA: ['The diagram engine returned incompatible feature metadata. Reload the page and Generate again.', ['reload', 'retry']],
  RESULT_INVALID: ['The generated result could not be accepted. Retry the operation.', ['retry', 'save-session']],
  // A shown Result without the request it was made from (a Result-only Session 27-39): Python reads its sequences
  // from the request the next Generate commits (D5).
  SEQUENCE_REQUIRES_GENERATE: ['This Result was saved by an older gbdraw, so its sequences are read after the next Generate.', ['generate']],
  EXPORT_INPUT: ['The current diagram could not be prepared for export. Generate a diagram and retry.', ['generate', 'retry']],
  EXPORT_DIMENSIONS: ['The diagram has no usable export dimensions. Generate again or use SVG.', ['generate', 'use-svg']],
  PNG_DPI: ['The PNG DPI is invalid. Use a positive finite DPI and retry.', ['edit-dpi', 'retry']],
  EXPORT_CONVERSION: ['The browser could not convert the diagram for export. Retry or use SVG.', ['retry', 'use-svg']],
  PDF_LIBRARY: ['The PDF libraries could not start. Reload and retry or use SVG.', ['reload', 'retry', 'use-svg']],
  PDF_GLYPH: ['PDF fonts cannot represent this text. Use SVG to retain this text.', ['use-svg']],
  // Python checks the size before anything is drawn or downloaded (gbdraw/render/export_limits.py).
  EXPORT_LIMIT: ['The figure is larger than this file format allows, so no file was made.', ['use-svg']],
  CLEANUP_FAILED: ['Temporary resource cleanup failed. Save a Session and reload before retrying.', ['save-session', 'reload', 'retry']]
});

// The sequence reads of a shown Result (D5) give their own meaning to a code that names another operation's failure.
/** @type {Readonly<Record<string, Readonly<Record<string, readonly [string, readonly string[], string?]>>>>} */
const OPERATION_DEFINITIONS = Object.freeze({
  readFeatureSequences: Object.freeze(/** @type {Readonly<Record<string, readonly [string, readonly string[], string?]>>} */ ({
    FEATURE_IDENTITY: ['The request of the shown Result does not draw this feature, so its sequence cannot be read. Generate again; if it continues, save a Session for investigation.', ['retry', 'save-session']]
  }))
});

/** @type {Readonly<Record<string, string>>} */
const FIELD_LABELS = Object.freeze({ protein_blastp_max_hits: 'Pairwise display max hits',
  losatp_max_hits: 'Pairwise display max hits', losatp_max_target_seqs: 'LOSATP Max target seqs',
  losatp_member_max_hits: 'Max hits per protein for groups',
  output_prefix: 'Output Prefix',
  losat_gencode: 'LOSAT Gencode', conservation_losat_gencodes: 'Comparison gencode', reference_gencode: 'Reference gencode',
  arrow_head_length_ratio: 'Arrow head length ratio', arrow_shaft_width_ratio: 'Arrow shaft width ratio',
  keep_definition_left_aligned: 'Lock Definition Column', window: 'Window', step: 'Step',
  depth_window: 'Depth Window', depth_step: 'Depth Step', dinucleotide: 'Dinucleotide',
  evalue: 'E-value', bitscore: 'Bitscore', identity: 'Identity', alignment_length: 'Alignment length', wrap_length: 'Wrap Rows', anchor_feature: 'Feature selector',
  plot_title_font_size: 'Plot Title Font Size', font_family: 'Font', feature_width_circular: 'Feature Track Width',
  depth_width_circular: 'Depth Width', gc_content_width_circular: 'GC Content Width',
  gc_content_radius_circular: 'GC Content Radius', gc_skew_width_circular: 'GC Skew Width',
  gc_skew_radius_circular: 'GC Skew Radius', match_height: 'Match height',
  conservation_ring_width: 'Ring Width', conservation_ring_gap: 'Ring Gap',
  center_reserved_radius: 'Center Reserved Radius', multi_record_min_radius_ratio: 'Min Radius Ratio',
  multi_record_column_gap_ratio: 'Column Gap Ratio', multi_record_row_gap_ratio: 'Row Gap Ratio',
  orthogroup_member_max_hits: 'Max hits per protein for groups', protein_blastp_candidate_limit: 'LOSATP Max target seqs',
  protein_blastp_mode: 'LOSATP mode', collinear_min_anchors: 'Min block genes', collinear_max_unit_gap: 'Max unit gap',
  collinear_max_diagonal_drift: 'Diagonal drift', collinear_max_conflicts: 'Merge conflicts',
  collinear_max_paralog_links_per_orthogroup: 'Paralog links per group',
  collinear_infer_orthogroups: 'Build similarity groups with within-record searches', collinear_search_scope: 'Search between',
  collinear_unit_mode: 'Unit mode', collinear_anchor_mode: 'Anchor mode', collinear_merge_orientation: 'Merge orientation',
  collinear_color_mode: 'Color mode', circular_multi_record_size_mode: 'Record Size Mode', linear_track_layout: 'Track Layout',
  linear_label_placement: 'Label Placement', overflow: 'Overflow behavior', layer: 'Annotation layer', set_id: 'Annotation set',
  anchor_slot: 'Enabled drawable anchor', to_slot: 'Extend to', renderer: 'Renderer', side: 'Track placement', spacing: 'Spacing',
  depth_min: 'Depth Min', depth_max: 'Depth Max', min_depth: 'Depth Min', max_depth: 'Depth Max' });
const CODES = new Set(Object.keys(DEFINITIONS));
// Python exception classes the adapter may name for RENDER_FAILED (class names, never messages).
const EXCEPTION_TYPES = new Set(`ValidationError ParseError ConfigError ExportError GbdrawError ZeroDivisionError
FloatingPointError OverflowError ArithmeticError AssertionError AttributeError IndexError KeyError LookupError
NotImplementedError RecursionError RuntimeError TypeError UnboundLocalError NameError ValueError Exception`.split(/\s+/));
// Locators shown in the summary; indexes are zero-based, ordinals one-based.
/** @type {Readonly<Record<string, string>>} */
const ORDINAL_LABELS = Object.freeze({ DECORATION_CONTINUITY: 'Result {n}.', COMPARISON_INPUT: 'Comparison sequence {n}.' });
// Fields whose ordinal counts the Linear cards (PHRASES.sequence), not comparison sequences: a card's LOSAT Gencode.
const CARD_ORDINAL_FIELDS = new Set(['losat_gencode']);
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
/** @typedef {{ remedy: string, options: number }} SwitchedTrack */
/** @type {Readonly<Record<string, SwitchedTrack>>} */
const SWITCHED_TRACKS = Object.freeze({
  gc_content: { remedy: 'Select Hide GC Content, or reduce GC Content Width (with Custom Track Slots on, the row width).', options: 2 },
  gc_skew: { remedy: 'Select Hide GC Skew, or reduce GC Skew Width (with Custom Track Slots on, the row width).', options: 2 },
  ticks: { remedy: 'Clear Show Coordinate Scale (with Custom Track Slots on, remove the Ticks row).', options: 1 }
});
// The room kept for the tracks inside a failed row, one whole sentence per set of switched tracks,
// by their `SWITCHED_TRACKS` keys in that table's order (a sentence is never built from track names).
/** @type {Readonly<Record<string, string>>} */
const KEPT_INSIDE_AMONG = Object.freeze({
  gc_content: 'Room kept toward the center for the inside tracks (among them GC content): {px} px.',
  gc_skew: 'Room kept toward the center for the inside tracks (among them GC skew): {px} px.',
  ticks: 'Room kept toward the center for the inside tracks (among them the scale): {px} px.',
  'gc_content,gc_skew': 'Room kept toward the center for the inside tracks (among them GC content and GC skew): {px} px.',
  'gc_content,ticks': 'Room kept toward the center for the inside tracks (among them GC content and the scale): {px} px.',
  'gc_skew,ticks': 'Room kept toward the center for the inside tracks (among them GC skew and the scale): {px} px.',
  'gc_content,gc_skew,ticks': 'Room kept toward the center for the inside tracks (among them GC content, GC skew, and the scale): {px} px.'
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
  keptBelow: 'Room kept toward the center for the inside tracks: {px} px.',
  // OV-1389 review item 6: each remedy is a whole sentence of `SWITCHED_TRACKS`, after one of these.
  makeRoom: 'To make room, do the following.', makeRoomAny: 'To make room, do one of the following.',
  // OV-2179: the room the labels not placed need (CAPACITY-STOP), the largest of each side, as Python sends it.
  outerLabelRoom: 'Outer labels need about {ratio} times the room around the map.',
  innerLabelRoom: 'Inner labels need about {ratio} times the room around the map.',
  record: 'Record: {record}.', recordPosition: 'Record #{position}.', letter: 'First letter that is not a nucleotide code: {letter}.',
  queryLetter: 'First character that is neither a nucleotide nor an amino acid code: {letter}.',
  position: 'At character {position}.', requiredColumns: 'Required columns: {count}.', loadedRecords: 'Loaded records: {count}.',
  exception: 'Python exception: {type}.', applyAgain: 'Review the updated preview and Apply again.',
  schema: 'Load a Session file that this version of gbdraw supports.',
  notASession: 'This file is not a gbdraw Session: it is not valid JSON.',
  // Receipt session.file-container: the settings and results (`gbdraw-session.json`) over the browser's cap
  // (services/session-file.js), and browser storage too short for a Session's inputs; sizes in MB to one decimal.
  sessionNotSaved: 'The Session was not saved: its settings and results are {size} MB, and the browser loads at most {limit} MB.',
  sessionSaveRemedy: 'Draw a smaller region, or turn off tracks or labels you do not need, then save again.',
  sessionNotLoaded: 'The Session was not loaded: its settings and results are {size} MB, and the browser loads at most {limit} MB.',
  sessionCommandLine: 'The gbdraw command line can draw this Session.',
  storageNotLoaded: 'The Session was not loaded: its inputs need {need} MB of browser storage and {free} MB is free.',
  storageLoadRemedy: 'Free up disk space, then load it again.',
  storageNotSaved: 'The Session was not saved: it needs up to {need} MB of browser storage and {free} MB is free.',
  storageNotSavedNoEstimate: 'The Session was not saved: the browser has no room for it (up to {need} MB).',
  storageSaveRemedy: 'Free up disk space, or close private windows, then save again.',
  largestInputs: 'Largest inputs: {inputs}.',
  // INPUT_UNREADABLE of a gzip input over its expanded cap: the file the browser read, and the cap in MiB.
  inputFile: 'File: {name}.', gzipLimit: 'gbdraw reads gzip files that expand to at most {limit} MiB.',
  // LOSAT_RUNTIME of a LOSAT run that exited with an error (services/losat-runtime.js `losatRunFailure`): LOSAT's
  // own first error line, shown as written (`shownLosatMessage`).
  losatReported: 'LOSAT reported: {message}',
  // GENETIC_CODE: the NCBI tables LOSAT reads, as Python names them (`LOSAT_GENETIC_CODE_TEXT`).
  geneticCodes: 'Genetic codes LOSAT reads: {codes}.',
  // The first sentence of a LOSAT failure that repeats (`losatFailureRepeats`) on Generate, in place of the code's
  // message and the reason's sentence: the one action leads to the comparison settings, never to Retry (OV-3317).
  losatRefused: 'LOSAT stopped the comparison search. Generating again with the same settings and inputs stops the same way; change what LOSAT names below, then Generate again.',
  losatTrapped: 'LOSAT stopped unexpectedly during the comparison search. Generating again with the same settings and inputs stops the same way; change the comparison settings or inputs, then Generate again.',
  // The first sentence of a threaded trap on Generate, in place of the code's message and the reason's sentence: the
  // serial engine may run the search, so it names Serial and Auto and keeps Retry.
  losatThreadedTrapped: 'LOSAT stopped unexpectedly in threaded execution. Select Serial or Auto execution, or change the comparison settings or inputs, then Generate again.',
  // LOSAT_RUNTIME in Sequence search (OV-3288) and Find in other files: the code's message names Generate, which
  // neither search uses; each runs from its own Search button.
  searchFailed: 'The search could not run or returned unusable output. Check the search settings, then click Search again.',
  // LOSAT_THREADING_UNAVAILABLE in the same two searches.
  searchThreadsUnavailable: 'Threaded LOSAT execution is unavailable in this browser environment. Select Serial or Auto execution, then click Search again.',
  // EXPORT_LIMIT (gbdraw/render/export_limits.py): the limit and size in the format's unit, and what would fit.
  exportSizePx: 'PNG is limited to {limit} px per side, and this figure is {width} x {height} px.',
  exportSizePt: 'PDF is limited to {limit} pt per side, and this figure is {width} x {height} pt.',
  exportRowLength: 'A row length of {length} bp should fit (estimated from the drawn height).',
  exportDpi: 'A DPI of {dpi} or lower would fit.',
  // FIGURE_TOO_NARROW, FIGURE_TOO_SHORT: Python's lengths as its message writes them (`180 mm`).
  figureTooNarrow: 'The figure cannot be {requested} wide: the parts that do not shrink need {minimum} ({parts}).',
  figureTooShort: 'The figure cannot be {requested} tall at {atWidth} wide: it needs {minimum} ({parts}).',
  figureWidthRemedy: 'Set Figure width to at least {minimum}, or leave it empty (Auto).',
  figureWidthLegendRemedy: 'Set Figure width to at least {minimum}, place the Legend at the top or bottom, or leave it empty (Auto).',
  figureHeightRemedy: 'Set Figure height to at least {minimum}, or leave it empty (Auto).',
  figureParts2: '{first}, {second}', figureParts3: '{first}, {second}, {third}',
  figureParts4: '{first}, {second}, {third}, {fourth}',
  liveEditFailed: 'Live edit failed: direct edits already applied are kept; geometry may still need updating.',
  liveEditRetry: 'Retry the live edit or use Generate.',
  liveEditChange: 'Change the edit, or change the settings and use Generate.',
  operationError: 'Operation error'
});
// A code's message for one operation, where the code's own message names a step that operation does not have.
/** @type {Readonly<Record<string, Readonly<Record<string, string>>>>} */
const OPERATION_MESSAGES = Object.freeze({
  LOSAT_RUNTIME: Object.freeze({ searchHitAnnotationSet: PHRASES.searchFailed, findInOtherFiles: PHRASES.searchFailed }),
  LOSAT_THREADING_UNAVAILABLE: Object.freeze({
    searchHitAnnotationSet: PHRASES.searchThreadsUnavailable, findInOtherFiles: PHRASES.searchThreadsUnavailable
  })
});
// The reasons whose Load refusal says the file is not a Session, and that sentence alone.
/** @type {Readonly<Record<string, string>>} */
const NOT_A_SESSION = Object.freeze({
  JSON_FORMAT: PHRASES.notASession, CONTAINER_DAMAGED: REASONS.CONTAINER_DAMAGED, CONTAINER_NO_DOCUMENT: REASONS.CONTAINER_NO_DOCUMENT
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
 * LOSAT's words as the summary shows them: the hyphens that start an option name (`-searchsp`) are non-breaking, so a
 * line never ends between them and the name. The Diagnostics keep the words as LOSAT wrote them.
 * @param {string} message
 */
const shownLosatMessage = (message) => message.replace(/(^|[\s'"(<[])(-+)/g, (_whole, before, hyphens) => before + '\u2011'.repeat(hyphens.length));

/**
 * The correction text of a diagnostic reason: the catalog key a hint translates, the same text a refusal
 * with that reason shows. `null` for a reason this owner does not word.
 * @param {string} reason
 * @returns {string | null}
 */
export const reasonText = (reason) => (Object.hasOwn(REASONS, reason) ? REASONS[reason] : null);
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
 * @property {string} [recordName] the record a FASTA-only refusal names, as Python sends it (Python decides it can be shown)
 * @property {number} [recordPosition] the 1-based position of that record in its file, when Python cannot show its ID
 * @property {string} [sessionTable]
 * @property {'circular' | 'linear'} [diagramMode]
 * @property {string} [drawingId] the drawing whose Result the failure concerns, when it is not the one shown
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
 * @property {string} [geneticCodes] The NCBI genetic codes LOSAT reads, as Python names them ("1–6, 9–16, 21–33").
 * @property {string} [cardUid] The Linear card whose LOSAT Gencode was refused, whose input "Show in Comparison" focuses.
 * @property {number} [sessionBytes] The size in bytes of the settings and results (`gbdraw-session.json`) of a Session
 *   that Save or Load refused.
 * @property {number} [sessionByteLimit] The size in bytes of settings and results that the browser loads at most.
 * @property {number} [storageNeedBytes] The browser storage in bytes that a refused Load or Save needs.
 * @property {number} [storageFreeBytes] The free browser storage in bytes the browser estimated; absent with no estimate.
 * @property {string} [largestInputs] The file names of the largest inputs of a Session that Save refused.
 * @property {string} [inputName] The name of the input file that could not be read (a gzip file over its expanded cap).
 * @property {string} [losatMessage] The first line LOSAT wrote about a run that exited with an error, as written.
 * @property {number} [limit] An export format's size limit per side, in its unit (EXPORT_LIMIT), or the expanded size
 *   cap of a gzip input in MiB (GZIP_EXPANDED_SIZE).
 * @property {number} [width] The figure's width in the export format's unit (EXPORT_LIMIT).
 * @property {number} [height] The figure's height in the export format's unit (EXPORT_LIMIT).
 * @property {number} [rowLength] The longer Wrap Rows length in bp that should fit (EXPORT_LIMIT).
 * @property {number} [dpi] The highest PNG DPI that fits (EXPORT_LIMIT).
 * @property {string} [requested] A refused figure size (FIGURE_TOO_NARROW, FIGURE_TOO_SHORT): the requested side, as
 *   Python's message writes a length (`180 mm`); `minimum`, `atWidth` and the parts likewise.
 * @property {string} [minimum]
 * @property {string} [atWidth]
 * @property {string} [primary]
 * @property {string} [legend]
 * @property {string} [title]
 * @property {string} [padding]
 * @property {FigurePrimaryItem} [primaryItem]
 * @property {true} [legendRemedy]
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
  // Python alone decides whether a record ID can be shown (gbdraw/web_support/error_adapter.py); else it sends the position.
  if (typeof value.recordName === 'string' && value.recordName !== '') context.recordName = value.recordName;
  if (typeof value.sessionTable === 'string' && Object.hasOwn(SESSION_TABLE_LABELS, value.sessionTable)) context.sessionTable = value.sessionTable;
  if (Number.isSafeInteger(value.sessionVersion)) context.sessionVersion = value.sessionVersion;
  if (Number.isSafeInteger(value.latestSessionVersion)) context.latestSessionVersion = value.latestSessionVersion;
  if (typeof value.geneticCodes === 'string' && /^[0-9\u2013, ]{1,80}$/.test(value.geneticCodes)) context.geneticCodes = value.geneticCodes;
  if (typeof value.cardUid === 'string' && /^[A-Za-z0-9_.:~-]{1,128}$/.test(value.cardUid)) context.cardUid = value.cardUid;
  if (typeof value.readSessionVersions === 'string' && /^[0-9\u2013, and]{1,200}$/.test(value.readSessionVersions)) {
    context.readSessionVersions = value.readSessionVersions;
  }
  // A refused Save: both sizes together, or neither; the names are one line of at most 200 characters.
  if (Number.isSafeInteger(value.sessionBytes) && Number.isSafeInteger(value.sessionByteLimit)
    && value.sessionBytes > 0 && value.sessionByteLimit > 0) {
    context.sessionBytes = value.sessionBytes;
    context.sessionByteLimit = value.sessionByteLimit;
  }
  // A refused Load or Save for browser storage: what it needs, and what is free when the browser estimated it.
  if (Number.isSafeInteger(value.storageNeedBytes) && value.storageNeedBytes >= 0) {
    context.storageNeedBytes = value.storageNeedBytes;
    if (Number.isSafeInteger(value.storageFreeBytes) && value.storageFreeBytes >= 0) context.storageFreeBytes = value.storageFreeBytes;
  }
  const inputs = typeof value.largestInputs === 'string' ? value.largestInputs.replace(/[\s\u0000-\u001f\u007f]+/g, ' ').trim().slice(0, 200).trimEnd() : '';
  if (inputs) context.largestInputs = inputs;
  // Python's Generate refusal cuts the name at the same bound (gbdraw/web_support/error_adapter.py).
  const inputName = typeof value.inputName === 'string'
    ? value.inputName.replace(/[\s\u0000-\u001f\u007f]+/g, ' ').trim().slice(0, INPUT_NAME_CHARACTERS).trimEnd() : '';
  if (inputName) context.inputName = inputName;
  // LOSAT's own words: one line of at most 200 characters, shown as text.
  const losatMessage = typeof value.losatMessage === 'string' ? value.losatMessage.replace(/[\s\u0000-\u001f\u007f]+/g, ' ').trim().slice(0, 200).trimEnd() : '';
  if (losatMessage) context.losatMessage = losatMessage;
  // The diagram mode whose Result the failure concerns, when it is not the one shown (E1).
  if (value.diagramMode === 'circular' || value.diagramMode === 'linear') context.diagramMode = value.diagramMode;
  if (typeof value.drawingId === 'string' && /^[A-Za-z0-9_-]{1,64}$/.test(value.drawingId)) context.drawingId = value.drawingId;
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
  // A refused figure size: Python's lengths (`180 mm`), each a bounded number and unit, and its primary item.
  if (value.reason === 'FIGURE_TOO_NARROW' || value.reason === 'FIGURE_TOO_SHORT') {
    for (const key of /** @type {const} */ (['requested', 'minimum', 'atWidth', 'primary', 'legend', 'title', 'padding'])) {
      if (isFigureLengthText(value[key])) context[key] = value[key];
    }
    if (typeof value.primaryItem === 'string' && Object.hasOwn(FIGURE_PRIMARY_ITEMS, value.primaryItem)) {
      context.primaryItem = /** @type {FigurePrimaryItem} */ (value.primaryItem);
    }
    if (value.legendRemedy === true) context.legendRemedy = true;
  }
  for (const key of /** @type {Array<'position' | 'row' | 'column' | 'inputOrdinal' | 'recordIndex' | 'seriesIndex' | 'slotIndex' | 'recordCount' | 'columnCount' | 'codepoint' | 'innerPx' | 'outerPx' | 'keptInsidePx' | 'placementIndex' | 'queryRecordIndex' | 'subjectRecordIndex' | 'limit' | 'width' | 'height' | 'rowLength' | 'dpi' | 'recordPosition'>} */ ([
    'position', 'row', 'column', 'inputOrdinal', 'recordIndex', 'seriesIndex', 'slotIndex', 'recordCount', 'columnCount', 'codepoint', 'innerPx', 'outerPx', 'keptInsidePx', 'placementIndex', 'queryRecordIndex', 'subjectRecordIndex',
    'limit', 'width', 'height', 'rowLength', 'dpi', 'recordPosition'
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

/**
 * The sentences of a Session that Save or Load refused for the browser's
 * limits (receipt session.file-container), or null for any other failure:
 * browser storage too short (`storageNeedBytes`), or settings and results over
 * the cap (`sessionBytes`). A JSON Session's refusal without sizes keeps the
 * code's text.
 * @param {Pick<UserFacingError, 'code' | 'operation' | 'context'>} result
 * @returns {ErrorTextPart[] | null}
 */
const sessionLimitParts = ({ code, operation, context }) => {
  const { sessionBytes, sessionByteLimit, storageNeedBytes, storageFreeBytes, largestInputs } = context;
  const save = operation === 'session-save';
  if (code !== 'SESSION_SIZE_LIMIT' || (!save && operation !== 'session-load')) return null;
  /** @param {number} bytes */
  const megabytes = (bytes) => (bytes / 1048576).toFixed(1);
  if (storageNeedBytes !== undefined) {
    const need = megabytes(storageNeedBytes);
    const free = storageFreeBytes === undefined ? undefined : megabytes(storageFreeBytes);
    if (!save) {
      return free === undefined ? null : [{ text: PHRASES.storageNotLoaded, params: { need, free } },
        { text: PHRASES.storageLoadRemedy }, { text: PHRASES.sessionCommandLine }];
    }
    return [free === undefined ? { text: PHRASES.storageNotSavedNoEstimate, params: { need } }
      : { text: PHRASES.storageNotSaved, params: { need, free } },
    ...(largestInputs === undefined ? [] : [{ text: PHRASES.largestInputs, params: { inputs: largestInputs } }]),
    { text: PHRASES.storageSaveRemedy }];
  }
  if (sessionBytes === undefined || sessionByteLimit === undefined) return null;
  const params = { size: megabytes(sessionBytes), limit: Math.round(sessionByteLimit / 1048576) };
  return save ? [{ text: PHRASES.sessionNotSaved, params }, { text: PHRASES.sessionSaveRemedy }]
    : [{ text: PHRASES.sessionNotLoaded, params }, { text: PHRASES.sessionCommandLine }];
};

// The reasons only a Session load gives. `FIELDS` also comes from other work,
// so it counts as a load only in the Load Session operation.
const SESSION_LOAD_REASONS = new Set(['SESSION_FORMAT', 'SESSION_FIELDS', 'SESSION_TOO_NEW', 'SESSION_VERSION_UNSUPPORTED', 'SESSION_TOO_OLD']);
/**
 * @param {{ code: string, operation: string }} result
 * @param {string | undefined} field
 * @param {string | undefined} reason
 */
const sessionLoadRejected = ({ code, operation }, field, reason) => code === 'INPUT_INVALID' && field === 'schema'
  && (SESSION_LOAD_REASONS.has(String(reason)) || (reason === 'FIELDS' && operation === 'session-load'));

// A refused figure size (gbdraw/layout/figure_size.py): its primary item by Python's id, each part with its length.
/** @typedef {'record_ends' | 'records' | 'map' | 'maps'} FigurePrimaryItem */
/** @type {Readonly<Record<FigurePrimaryItem, string>>} */
const FIGURE_PRIMARY_ITEMS = Object.freeze({
  record_ends: 'record names and labels past the sequence ends {length}', records: 'records {length}',
  map: 'map at its smallest {length}', maps: 'maps at their smallest {length}'
});
const FIGURE_PARTS = Object.freeze({ legend: 'legend {length}', title: 'title {length}', padding: 'padding {length}' });
// A length as Python's message writes it (`180 mm`): a bounded number and one of Python's units
// (gbdraw/core/figure_length.py, generated into the mode profiles). The one check of it in the Web: a refused size
// here, and the larger-figure estimate in services/label-fit.js.
const FIGURE_LENGTH_TEXT = new RegExp(`^[0-9]{1,30}(\\.[0-9]{1,6})? (${MODE_PROFILE_DATA.figureLength.units.join('|')})$`);
/**
 * Whether `value` is a length as Python's messages write it (`180 mm`).
 * @param {unknown} value
 * @returns {value is string}
 */
export const isFigureLengthText = (value) => typeof value === 'string' && FIGURE_LENGTH_TEXT.test(value);
/**
 * The sentences of a refused figure size, in Python's numbers: what does not fit and the smallest size, then what to
 * do. None without Python's lengths (the reason's own text then stands).
 * @param {UserFacingContext} context
 * @returns {ErrorTextPart[]}
 */
const figureSizeParts = ({ reason, requested, minimum, atWidth, primary, primaryItem, legend, title, padding, legendRemedy }) => {
  const narrow = reason === 'FIGURE_TOO_NARROW';
  if ((!narrow && reason !== 'FIGURE_TOO_SHORT') || requested === undefined || minimum === undefined) return [];
  if (!narrow && atWidth === undefined) return [];
  /** @type {ErrorTextPart[]} */
  const parts = [
    ...(primary !== undefined && primaryItem !== undefined ? [{ text: FIGURE_PRIMARY_ITEMS[primaryItem], params: { length: primary } }] : []),
    ...(/** @type {const} */ (['legend', 'title', 'padding'])).flatMap((key) => {
      const length = { legend, title, padding }[key];
      return length === undefined ? [] : [{ text: FIGURE_PARTS[key], params: { length } }];
    })
  ];
  const joiner = [null, null, PHRASES.figureParts2, PHRASES.figureParts3, PHRASES.figureParts4][parts.length];
  if (parts.length === 0 || (parts.length > 1 && !joiner)) return [];
  const names = /** @type {const} */ (['first', 'second', 'third', 'fourth']);
  const listed = parts.length === 1 ? parts[0]
    : { text: /** @type {string} */ (joiner), params: Object.fromEntries(parts.map((part, index) => [names[index], part])) };
  return [
    narrow ? { text: PHRASES.figureTooNarrow, params: { requested, minimum, parts: listed } }
      : { text: PHRASES.figureTooShort, params: { requested, atWidth: /** @type {string} */ (atWidth), minimum, parts: listed } },
    { text: narrow ? (legendRemedy ? PHRASES.figureWidthLegendRemedy : PHRASES.figureWidthRemedy) : PHRASES.figureHeightRemedy, params: { minimum } }
  ];
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
  const parts = keptInsidePx === undefined ? [] : [{
    text: kept.length ? KEPT_INSIDE_AMONG[Object.keys(SWITCHED_TRACKS).filter((track) => kept.includes(track)).join(',')] : PHRASES.keptBelow,
    params: { px: keptInsidePx }
  }];
  const tracks = [...new Set([...(slotTrack ? [slotTrack] : []), ...kept])].map((track) => SWITCHED_TRACKS[track]);
  if (!tracks.length) return parts;
  const options = tracks.reduce((count, track) => count + track.options, 0);
  parts.push({ text: options > 1 ? PHRASES.makeRoomAny : PHRASES.makeRoom }, ...tracks.map((track) => ({ text: track.remedy })));
  return parts;
};

/**
 * Whether the same LOSAT search would stop the same way: LOSAT exited with an error and said why (`losatMessage`,
 * services/losat-runtime.js `losatRunFailure`, OV-3317), or LOSAT's code trapped (`losatTrapFailure`, OV-3327). The
 * failure offers the comparison settings and never Retry.
 * @param {{ code: string, context: UserFacingContext }} error
 * @returns {boolean}
 */
export const losatFailureRepeats = (error) => error.code === 'LOSAT_RUNTIME'
  && (error.context.losatMessage !== undefined || error.context.reason === 'TRAP');

/**
 * The sentence that lists the genetic codes LOSAT reads, as a refused code's Generation Error says it: a Linear card
 * marks a typed code LOSAT does not read with it before Generate (OV-3326).
 * @param {string} codes The tables as Python names them (`LOSAT_GENETIC_CODE_TEXT`).
 * @returns {ErrorTextPart[]}
 */
export const losatGeneticCodeHint = (codes) => [{ text: PHRASES.geneticCodes, params: { codes } }];

/**
 * Whether the same Generate would stop on the same comparison setting again: a LOSAT failure that repeats
 * (`losatFailureRepeats`), or a TLOSATX genetic code LOSAT does not read (OV-3326). The failure offers the comparison
 * settings and never Retry, in both modes.
 * @param {{ code: string, context: UserFacingContext }} error
 * @returns {boolean}
 */
export const comparisonCauseRepeats = (error) => losatFailureRepeats(error)
  || (error.code === 'COMPARISON_INPUT' && error.context.reason === 'GENETIC_CODE');

/**
 * The message and actions of a failure's code: the operation's own definition where it has one, else the code's.
 * @param {{ code: string, operation: string }} error
 * @returns {readonly [string, readonly string[], string?]}
 */
const definitionOf = (error) => (Object.hasOwn(OPERATION_DEFINITIONS, error.operation)
  && Object.hasOwn(OPERATION_DEFINITIONS[error.operation], error.code)
  ? OPERATION_DEFINITIONS[error.operation] : DEFINITIONS)[error.code];

/**
 * The first sentence of a failure: the code's message for its operation where the code's own names a step that
 * operation does not have (`OPERATION_MESSAGES`), else for a threaded trap or a LOSAT failure that repeats the sentence
 * that says so, which also stands for the reason's sentence (`coversReason`), else the code's message.
 * @param {{ code: string, operation: string, context: UserFacingContext }} error
 * @returns {{ text: string, coversReason: boolean }}
 */
const leadSentence = (error) => {
  const operationMessages = Object.hasOwn(OPERATION_MESSAGES, error.code) ? OPERATION_MESSAGES[error.code] : {};
  if (Object.hasOwn(operationMessages, error.operation)) return { text: operationMessages[error.operation], coversReason: false };
  if (error.code === 'LOSAT_RUNTIME' && error.context.reason === 'THREADED_TRAP') return { text: PHRASES.losatThreadedTrapped, coversReason: true };
  if (!losatFailureRepeats(error)) return { text: definitionOf(error)[0], coversReason: false };
  return { text: error.context.reason === 'TRAP' ? PHRASES.losatTrapped : PHRASES.losatRefused, coversReason: true };
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
  const [, actions] = definitionOf(result);
  const lead = leadSentence(result);
  const { inputOrdinal, row, column, slotIndex, seriesIndex, innerPx, outerPx, configPath, featureCaption, sessionTable, diagramMode } = result.context;
  const { featureId, featureType, featureStart, featureEnd, featureCount, field, reason, position, columnCount, recordCount } = result.context;
  const { recordName, recordPosition, codepoint } = result.context;
  // A letter Python names by its code point: itself when printable, else its U+ form.
  const letter = codepoint === undefined ? '' : codepoint > 0x20 && codepoint < 0x7f
    ? `"${String.fromCodePoint(codepoint)}"` : `U+${codepoint.toString(16).toUpperCase().padStart(4, '0')}`;
  const featureSpan = featureStart !== undefined && featureEnd !== undefined ? `${featureStart}..${featureEnd}` : '';
  const featureName = [[featureType, featureSpan].filter(Boolean).join(' '), featureId !== undefined ? `ID ${featureId}` : '']
    .filter(Boolean).join(', ');
  const { queryRecordIndex, subjectRecordIndex, labelFitCounts, labelCapacityOuter, labelCapacityInner, exceptionType } = result.context;
  const { inputName, losatMessage } = result.context;
  const { limit, width, height, rowLength, dpi, geneticCodes } = result.context;
  const exportSize = (reason === 'PNG_SIZE' || reason === 'PDF_SIZE') && limit !== undefined && width !== undefined
    && height !== undefined && { text: reason === 'PNG_SIZE' ? PHRASES.exportSizePx : PHRASES.exportSizePt,
    params: { limit: numberText(limit), width: numberText(width), height: numberText(height) } };
  // A Save or Load over the browser's limits: the refusal's own sentences (receipt session.file-container).
  const sessionLimit = sessionLimitParts(result);
  // A refused figure size names its control (Figure width, Figure height) in its own sentence, so it shows no
  // Setting line with the internal path.
  const figureSize = reason === 'FIGURE_TOO_NARROW' || reason === 'FIGURE_TOO_SHORT';
  const figureSizeSentences = figureSizeParts(result.context);
  // One part per sentence, in the order the summary shows them.
  /** @type {(ErrorTextPart | null | false | '' | undefined)[]} */
  const pieces = [
    !WHOLE_FAILURE_REASONS.has(/** @type {string} */ (reason)) && { text: lead.text },
    diagramMode !== undefined && { text: diagramMode === 'linear' ? PHRASES.diagramLinear : PHRASES.diagramCircular },
    sessionTable !== undefined && { text: PHRASES.sessionTable, params: { table: { text: SESSION_TABLE_LABELS[sessionTable] } } },
    inputOrdinal !== undefined && { text: (!CARD_ORDINAL_FIELDS.has(field ?? '') && ORDINAL_LABELS[result.code]) || PHRASES.sequence, params: { n: inputOrdinal } },
    // A comparison pair as the comparison panel numbers it (#1 to #2).
    queryRecordIndex !== undefined && subjectRecordIndex !== undefined
      && { text: PHRASES.pair, params: { query: queryRecordIndex + 1, subject: subjectRecordIndex + 1 } },
    featureCaption !== undefined && { text: PHRASES.featureCaption, params: { caption: featureCaption } },
    inputName !== undefined && { text: PHRASES.inputFile, params: { name: inputName } },
    reason === 'GZIP_EXPANDED_SIZE' && limit !== undefined && { text: PHRASES.gzipLimit, params: { limit } },
    recordName !== undefined && { text: PHRASES.record, params: { record: recordName } },
    recordName === undefined && recordPosition !== undefined && { text: PHRASES.recordPosition, params: { position: recordPosition } },
    row !== undefined && { text: PHRASES.line, params: { row } },
    column !== undefined && { text: PHRASES.column, params: { column } },
    slotIndex !== undefined && { text: PHRASES.trackRow, params: { row: slotIndex + 1 } },
    seriesIndex !== undefined && { text: PHRASES.depthSeries, params: { series: seriesIndex + 1 } },
    configPath !== undefined && !figureSize && { text: PHRASES.setting, params: { path: configPath } },
    featureName && { text: PHRASES.feature, params: { feature: featureName } },
    // `undefined > 1` is false: with no count there is no line.
    /** @type {number} */ (featureCount) > 1 && { text: PHRASES.featuresAffected, params: { count: /** @type {number} */ (featureCount) } },
    field && field !== 'schema' && { text: PHRASES.field, params: { field: FIELD_LABELS[field] ? { text: FIELD_LABELS[field] } : field } },
    exportSize,
    exportSize && rowLength !== undefined && { text: PHRASES.exportRowLength, params: { length: numberText(rowLength) } },
    exportSize && dpi !== undefined && { text: PHRASES.exportDpi, params: { dpi } },
    ...figureSizeSentences,
    !lead.coversReason && (sessionVersionGuidance(result.context) || (figureSizeSentences.length > 0 ? null
      : (REASONS[/** @type {string} */ (reason)] || DEFINITIONS[result.code][2] ? { text: REASONS[/** @type {string} */ (reason)] || /** @type {string} */ (DEFINITIONS[result.code][2]) } : null))),
    losatMessage !== undefined && { text: PHRASES.losatReported, params: { message: shownLosatMessage(losatMessage) } },
    reason === 'GENETIC_CODE' && geneticCodes !== undefined && { text: PHRASES.geneticCodes, params: { codes: geneticCodes } },
    reason === 'NOT_NUCLEOTIDE' && letter && { text: PHRASES.letter, params: { letter } },
    reason === 'NOT_SEQUENCE' && letter && { text: PHRASES.queryLetter, params: { letter } },
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
    sessionLoadRejected(result, field, reason) && { text: PHRASES.schema }
  ];
  // A file chosen in Load Session that is not JSON, or a damaged container or one without
  // gbdraw-session.json, is not a Session; every other Session reason keeps its correction.
  const notASession = result.operation === 'session-load' && result.code === 'INPUT_INVALID' && field === 'schema'
    && reason !== undefined && Object.hasOwn(NOT_A_SESSION, reason) ? NOT_A_SESSION[reason] : null;
  result.summaryParts = capParts(notASession ? [{ text: notASession }] : sessionLimit
    ?? pieces.filter(/** @returns {piece is ErrorTextPart} */ (piece) => Boolean(piece)),
  Number.isSafeInteger(summaryLimit) ? Math.max(0, Math.min(summaryLimit, 1000)) : 1000);
  result.summary = renderErrorText(result.summaryParts);
  const detail = [`Code: ${result.code}`, `Operation: ${result.operation}`, `Stage: ${result.stage}`,
    ...Object.entries(result.context).map(([key, item]) => `${key}: ${item}`),
    ...result.secondary.map((item) => `Secondary: ${item.code} / ${item.stage}`)].join('\n');
  result.details = [{ label: 'Diagnostics', text: detail.slice(0,
    Number.isSafeInteger(detailLimit) ? Math.max(0, Math.min(detailLimit, 4000)) : 4000) }].slice(0, 8);
  // A table stored in a Session is corrected in the Session file, not in this page.
  result.actions = sessionTable !== undefined || sessionLimit ? ['select-input']
    : comparisonCauseRepeats(result) ? ['edit-comparison']
    : [...(reason !== undefined && Object.hasOwn(REASON_ACTIONS, reason) ? REASON_ACTIONS[reason] : actions)];
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
