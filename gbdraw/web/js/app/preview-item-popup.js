// @ts-check
/** @import { DrawingState, LinearSeq, WebState } from '../state.js' */
/** @import { LegendRow } from '../services/drawing/paint.js' */
/** @import { DrawnTrack } from '../services/track-legend-rows.js' */
/** @import { SettingApplicationId } from '../services/setting-application.js' */
/** @import { RegionSelection } from '../services/region-axis.js' */
// The legend-entry and record-definition popups of the preview. A click (or
// Enter) on a Legend row or a record definition opens the match popup frame
// with a `previewItem`; this owner derives what the popup shows from the
// displayed Result and the current settings. It holds no edit of its own: the
// template calls the Legend panel's operations, and a record edit writes the
// setting the sidebar field writes, as one History step.
import { legendEntryKey } from '../services/drawing/paint.js';
import { previewItemPopupOf } from '../services/match-frame-popup.js';
import { plainTextLinearRecordLabel } from '../services/linear-comparisons.js';
import { specificRuleHasQualifier } from '../services/specific-color-rules.js';
import {
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle
} from '../services/linear-sources.js';
import { regionSelectionHeading } from '../services/region-axis.js';
import { readLinearSeqWindows } from '../services/drawing/inputs.js';
import {
  canonicalRecordReverseComplement, cropFormReverseComplement, linearCardCropForm, orientedCropForm
} from '../services/record-display-model.js';
import { orientLinearCard } from '../services/linear-card-windows.js';
import { appliedFeatureColorOf, colorSwatchView } from '../utils/color-utils.js';
import { msg, numberText, t, tn, tm, tx } from '../utils/i18n.js';

/**
 * @typedef {{ label: string, value: string }} PreviewItemFact
 * @typedef {{ feat?: unknown, qual?: unknown, val?: unknown }} ColorRuleLike
 * @typedef {{ start: number, end: number, reverseComplement: boolean }} RequestRecordWindow
 * @typedef {{
 *   recordKey?: unknown,
 *   cardinality?: unknown,
 *   presentation: { label?: unknown, subtitle?: unknown, reverseComplement: boolean },
 *   display?: { startCoordinate?: unknown } | null,
 *   windows: RequestRecordWindow[]
 * }} RequestRecordLike
 * @typedef {{ recordKey?: unknown, status?: unknown, anchorLabel?: unknown, rationaleLabel?: unknown }} AlignmentDecisionLike
 * @typedef {{ scope?: unknown, sourceUid?: unknown, recordId?: unknown, canonicalRecordKey?: unknown }} CommittedRecordRowLike
 * @typedef {'linear' | 'circular'} RecordScope
 * @typedef {'name' | 'subtitle' | 'reverse'} RecordField
 * @typedef {{
 *   name: string, subtitle: string, reverse: boolean, reverseAvailable: boolean,
 *   nameLabel: string, subtitleLabel: string, sequenceUid: string,
 *   fields: RecordPopupFields, settingsSelector: string
 * }} RecordDraft
 * @typedef {{ readonly key: string, readonly setting: SettingApplicationId }} RecordPopupField
 *   The draft key a popup field writes and the setting row (`SETTING_APPLICATION`) of its timing note.
 * @typedef {RecordPopupField & { readonly windows: string }} RecordPopupReverseField
 *   The Reverse complement field and the key of the windows its form lists, which hold the orientation while
 *   listed (`orientedCropForm`).
 * @typedef {{ readonly name: RecordPopupField, readonly subtitle: RecordPopupField, readonly reverse: RecordPopupReverseField }} RecordPopupFields
 * @typedef {ReturnType<typeof buildLegendEntryView>} LegendEntryView
 * @typedef {ReturnType<typeof buildRecordDefinitionView>} RecordDefinitionView
 */

/**
 * The members of the Web state this owner reads.
 * @typedef {object} PreviewItemState
 * @property {WebState['clickedPairwiseMatch']} clickedPairwiseMatch The open match-frame popup.
 * @property {{ value: Element | null }} svgContainer The preview's mounted Result.
 * @property {{ value: string }} mode
 * @property {WebState['extractedFeatures']} extractedFeatures
 * @property {() => DrawingState} activeDrawing The shown drawing, whose Result and record rows the popup shows.
 * @property {{ value: Record<string, Record<string, string>> | null }} paletteDefinitions
 * @property {{ value: string }} appliedPaletteName The palette of the displayed Result.
 */

/**
 * The root's way to bring a setting of the settings pane into view and focus it.
 * @typedef {(setting: HTMLElement) => unknown} RevealSetting
 */

/**
 * @typedef {object} PreviewItemPopupOptions
 * @property {PreviewItemState} state
 * @property {<T>(getter: () => T) => { readonly value: T }} computed Vue `computed`
 * @property {RevealSetting} revealSetting
 * @property {(entry: LegendRow) => string} readLegendEntryColor
 *   The color the Legend panel shows for a row.
 * @property {(caption: string) => string} [readShownLegendName]
 *   The name the Legend draws for the row that shows the key `caption` (never the key).
 * @property {(entry: LegendRow) => readonly ColorRuleLike[]} readLegendRowRules
 *   The color rules that draw a row of the displayed Result.
 * @property {(caption: string) => number} countLegendCaptionFeatures
 *   The number of features whose legend entry is `caption`.
 * @property {(entry: LegendRow) => DrawnTrack['kind'] | null} readLegendTrackKind
 *   The kind of the track whose row this is in the displayed Result (`trackLegendRows`); null: a feature row.
 * @property {(entry: LegendRow) => string | null} readLegendRowDefaultColorType
 *   The type whose default color draws the row and its Color sets (OV-1281), once the matches are ready.
 * @property {(type: string) => string | null} readUserDefaultColor
 *   The user default color of a feature type, or null when the user set none.
 * @property {() => { records?: readonly RequestRecordLike[] } | null} getCommittedRequest
 *   The request of the displayed Result.
 * @property {() => { records?: readonly AlignmentDecisionLike[] } | null} readAlignmentInspector
 *   The active similarity alignment as the drawer inspects it, or null.
 * @property {() => ReadonlyArray<CommittedRecordRowLike>} readCommittedRecordRows
 *   The drawn records of the displayed Result's request, with their File cards.
 * @property {() => CircularOwnRecordBlock | ''} readCircularOwnRecordBlock
 *   Why a Circular record's own fields (Record label, Subtitle, Reverse complement) are unavailable, or ''.
 * @property {(label: string, apply: () => void) => unknown} runUndoable History's one undoable step.
 * @property {() => unknown} sessionBusy A reason settings cannot change now, else a falsy value.
 * @property {(index: number, caption: string) => unknown} renameLegendRow
 *   The Legend panel's rename of row `index`, as one History step.
 * @property {() => Promise<unknown>} legendRenameSettled
 *   Resolves once the rename's choice dialog, if one is open, is answered.
 */

/** @param {unknown} value */
const text = (value) => String(value ?? '').trim();
/** @param {unknown} value */
const plain = (value) => (text(value) ? plainTextLinearRecordLabel(text(value)) : '');
/** @param {unknown} value */
const positive = (value) => {
  const number = Number(value);
  return value != null && Number.isFinite(number) && number > 0 ? number : null;
};
/** A Legend row's key, trimmed: the identity the popup finds the row by. @param {LegendRow} entry @returns {string} */
const rowKey = (entry) => text(legendEntryKey(entry));

/**
 * OV-1466: a rule without a feature type is a whole sentence of its own, so a
 * locale writes "any feature" where and how its sentence needs it. OV-2435: a
 * rule needs a qualifier (`specificRuleHasQualifier`); one without a qualifier
 * or a pattern is named by its feature type alone, so no sentence has an empty part.
 * @param {ColorRuleLike} rule
 */
const colorRuleText = (rule) => {
  const feature = text(rule?.feat);
  const qualifier = text(rule?.qual);
  const pattern = text(rule?.val);
  if (!specificRuleHasQualifier(rule) || !pattern) return feature ? t('{feature} color rule', { feature }) : t('Any feature color rule');
  return feature
    ? t('{feature} where {qualifier} matches “{pattern}”', { feature, qualifier, pattern })
    : t('Any feature where {qualifier} matches “{pattern}”', { qualifier, pattern });
};

// OV-1272: a track row names its plot; it lists no features, and its Color
// edits only the swatch (the plot's color is a Settings value).
/** @type {Readonly<Record<DrawnTrack['kind'], string>>} */
const TRACK_ROW_SUBTITLE = Object.freeze({ content: msg('Content plot'), skew: msg('Skew plot'), depth: msg('Depth plot') });
const TRACK_ROW_COLOR_SOURCE = msg('Changes the legend swatch only.');

/**
 * OIC-029: what the Legend item choice of the scope dialog and the Legend popup's Color do to a type's default
 * color. The one owner of that line; `DefaultColorScopeLine` draws it.
 * @typedef {{ kind: 'sets' | 'replaces', type: string, color: string | null }} DefaultColorScopeLine
 * @param {string | null | undefined} type
 * @param {string | null | undefined} replacedColor the user default color of `type`, if the user set one
 * @returns {DefaultColorScopeLine | null}
 */
export const defaultColorScopeLine = (type, replacedColor) => (
  type ? { kind: replacedColor ? 'replaces' : 'sets', type, color: replacedColor || null } : null
);

/**
 * The legend-entry popup view.
 * @param {{ entry: LegendRow, name: string, index: number, color: string, rules: readonly ColorRuleLike[],
 *   featureCount: number | null, featureTypes: ReadonlySet<string>, trackKind: DrawnTrack['kind'] | null,
 *   defaultColorType: string | null, replacedDefaultColor: string | null, pickerFallback: string }} input
 *   `name`: the name the Legend draws. `trackKind`: the kind of the track whose row it is; null: a feature row.
 *   `pickerFallback`: where the Color input starts on a row of color none.
 */
const buildLegendEntryView = ({
  entry, name, index, color, rules, featureCount, featureTypes, trackKind, defaultColorType, replacedDefaultColor, pickerFallback
}) => {
  const key = rowKey(entry);
  const ruleLines = trackKind ? [] : rules.map(colorRuleText);
  const defaultColorLine = trackKind || ruleLines.length ? null : defaultColorScopeLine(defaultColorType, replacedDefaultColor);
  // OV-1336: the OIC-029 line, once shown, says where the color comes from by itself; the caption gives way.
  // One rule line says where the color comes from by itself too; several get their count.
  const colorSource = defaultColorLine
    ? null
    : trackKind
      ? tm(TRACK_ROW_COLOR_SOURCE)
      : ruleLines.length
        ? (ruleLines.length > 1 ? tn(ruleLines.length, '{count} color rule', '{count} color rules') : null)
        : featureTypes.has(key)
          ? t('Default color of {type} features', { type: key })
          : t('Color set on individual features');
  return {
    kind: /** @type {const} */ ('legend-entry'),
    identity: `legend:${key}`,
    dialogLabel: t('Legend item: {name}', { name }),
    index,
    name,
    color,
    swatch: colorSwatchView(color, pickerFallback),
    title: name,
    subtitle: trackKind
      ? tm(TRACK_ROW_SUBTITLE[trackKind])
      : featureCount == null
        ? t('Counting features…')
        : tn(featureCount, '{count} feature', '{count} features', { count: numberText(featureCount) }),
    colorSource,
    defaultColorLine,
    ruleLines
  };
};

/** @param {readonly RequestRecordWindow[]} windows the request record's windows (request schema 10), none for the whole record */
const regionText = (windows) => (windows.length
  ? windows.map(({ start, end }) => `${numberText(start)}..${numberText(end)}`).join(', ')
  : t('Whole record'));

/** @param {AlignmentDecisionLike | null | undefined} decision */
const alignmentText = (decision) => {
  if (!decision) return t('Not aligned');
  if (decision.status === 'reference') return t('Reference ({anchor})', { anchor: text(decision.anchorLabel) });
  if (decision.status === 'aligned') return t('Aligned on {anchor}', { anchor: text(decision.anchorLabel) });
  const reason = text(decision.rationaleLabel);
  return reason ? t('Not moved: {reason}', { reason }) : t('Not moved');
};

/**
 * The record-definition popup view. `drawn` holds the definition lines the
 * Result draws, `record` the Result's request record (the facts are what the
 * Result shows), and `draft` the settings the sidebar edits, which apply on
 * Generate.
 * @param {{ scope: RecordScope, recordIndex: number, recordId: string,
 *   drawn: Record<string, string>, record: RequestRecordLike | null, draft: RecordDraft | null,
 *   editHint: string, alignment: AlignmentDecisionLike | null, alignmentShown: boolean }} input
 */
const buildRecordDefinitionView = ({
  scope, recordIndex, recordId, drawn, record, draft, editHint, alignment, alignmentShown
}) => {
  const accession = text(drawn.accession) || text(recordId);
  // The heading is the first name line the Result draws; a Result without a
  // definition is named by the request's label.
  const name = text(drawn.name)
    || (Object.values(drawn).some(Boolean) ? '' : plain(record?.presentation?.label));
  const named = Boolean(name) && name !== accession;
  const reversed = record !== null && canonicalRecordReverseComplement(record);
  const start = positive(record?.display?.startCoordinate);
  /** @type {PreviewItemFact[]} */
  const facts = [
    { label: t('Accession'), value: accession || '—' },
    // OV-502: only a Linear definition marks which line is the record's name.
    // A Circular definition's name lines may be an organism, a strain, a
    // replicon or a Record label, wrapped or not: they are the heading only.
    ...(scope === 'linear' && named ? [{ label: t('Organism'), value: name }] : []),
    { label: t('Length'), value: text(drawn.length) || '—' },
    { label: t('Reverse complement'), value: reversed ? t('Yes') : t('No') },
    { label: t('Display start'), value: start == null ? t('Record start') : numberText(start) },
    { label: t('Region'), value: regionText(record?.windows ?? []) }
  ];
  // The record's place in a similarity alignment, not a sequence alignment
  // (the match popup's Alignment section): a locale may name the two apart.
  if (alignmentShown) facts.push({ label: tx('similarity alignment', 'Alignment'), value: alignmentText(alignment) });
  return {
    kind: /** @type {const} */ ('record-definition'),
    identity: `record:${scope}:${recordIndex}:${recordId}`,
    dialogLabel: t('Record {accession}', { accession }),
    scope,
    recordIndex,
    title: name || accession,
    subtitle: [text(drawn.subtitle), named ? accession : ''].filter(Boolean).join(' · '),
    facts,
    draft,
    editHint
  };
};

/**
 * The texts of a definition group's lines by kind (`data-definition-line-kind`).
 * @param {Element | null | undefined} group @param {string} recordId
 */
const readDrawnDefinitionLines = (group, recordId) => {
  /** @type {Record<string, string>} */
  const lines = {};
  group?.querySelectorAll('[data-definition-line-kind]').forEach((element) => {
    const kind = element.getAttribute('data-definition-line-kind') || '';
    if (kind && !lines[kind]) lines[kind] = text(element.textContent);
  });
  if (Object.keys(lines).length) return lines;
  // A Circular definition marks no line kinds: its lines are the name lines,
  // then the accession, the length, and the GC content.
  const texts = Array.from(group?.querySelectorAll('text') || [], (element) => text(element.textContent)).filter(Boolean);
  const isFact = (/** @type {string} */ line) => /\sbp$|% GC$/.test(line);
  const named = texts.filter((line) => !isFact(line) && line !== recordId);
  lines.length = texts.find((line) => /\sbp$/.test(line)) || '';
  lines.name = named[0] || '';
  lines.subtitle = named.slice(1).join(' ');
  if (texts.includes(recordId)) lines.accession = recordId;
  return lines;
};

// The settings a record popup edits: those the sidebar writes for the drawn
// definition, each with the setting row that says when it applies (R1): the
// one owner of field -> setting key -> row for the popup and its notes. A
// Circular record selected on its own draws its Record label and Subtitle,
// which win over the drawing's Species and Strain.
/** @type {Readonly<Record<'linear' | 'circularRecord' | 'circularDrawing', RecordPopupFields>>} */
export const RECORD_POPUP_FIELDS = Object.freeze({
  linear: Object.freeze({
    name: { key: 'definition', setting: 'record-names' }, subtitle: { key: 'record_subtitle', setting: 'inputs' },
    reverse: { key: 'region_reverse', windows: 'windows', setting: 'inputs' }
  }),
  circularRecord: Object.freeze({
    name: { key: 'circular_record_label', setting: 'record-names' },
    subtitle: { key: 'circular_record_subtitle', setting: 'record-names' },
    reverse: { key: 'circular_reverse', windows: 'circular_windows', setting: 'inputs' }
  }),
  circularDrawing: Object.freeze({
    name: { key: 'species', setting: 'species' }, subtitle: { key: 'strain', setting: 'species-strain' },
    reverse: { key: 'circular_reverse', windows: 'circular_windows', setting: 'inputs' }
  })
});
/**
 * The write of a record popup's text field, or null when the field holds `next` already.
 * @param {Record<string, unknown>} target @param {string} key @param {string} next
 * @returns {(() => void) | null}
 */
const recordTextWrite = (target, key, next) => (target[key] === next ? null : () => { target[key] = next; });
/**
 * The write that draws a record popup's record in the direction `next` where its crop form holds the orientation:
 * its listed windows, else its Reverse complement (`orientedCropForm`, as Generate reads them); a Linear card's
 * through its windows' owner (`orientLinearCard`). Null when it is drawn so already.
 * @param {Record<string, unknown>} target @param {RecordPopupReverseField} field @param {boolean} next
 * @returns {(() => void) | null}
 */
const recordReverseWrite = (target, field, next) => {
  if (field === RECORD_POPUP_FIELDS.linear.reverse) {
    const card = /** @type {LinearSeq} */ (/** @type {unknown} */ (target));
    return cropFormReverseComplement(linearCardCropForm(card)) === next ? null : () => { orientLinearCard(card, next); };
  }
  const form = { windows: readLinearSeqWindows(target[field.windows]), reverse: Boolean(target[field.key]) };
  if (cropFormReverseComplement(form) === next) return null;
  return () => {
    const oriented = orientedCropForm(form, next);
    target[field.windows] = oriented.windows;
    target[field.key] = oriented.reverse;
  };
};
const CIRCULAR_RECORD = {
  fields: RECORD_POPUP_FIELDS.circularRecord,
  nameLabel: msg('Record label'), subtitleLabel: msg('Subtitle'),
  settingsSelector: '[data-circular-record-presentation] input[type="text"]'
};
const CIRCULAR_DRAWING = {
  fields: RECORD_POPUP_FIELDS.circularDrawing,
  nameLabel: msg('Species'), subtitleLabel: msg('Strain'), settingsSelector: '#circular-species'
};
// What a popup without edit fields says instead. A Circular record's sentence
// names the cause; the second and third are the sidebar's for the same cause.
const LINEAR_EDIT_HINT = msg('Edit this record in its File card in the sidebar.');
const CIRCULAR_EDIT_HINTS = {
  'not-inspected': msg('Inspect source records in the sidebar to edit this record\'s label.'),
  'multi-record-canvas': msg('Multi-Record Canvas uses a grid. Turn it off to select one record.'),
  'one-circle': msg('One circle draws every record. Turn it off to select one record.'),
  'no-single-record': msg('Select one record to edit crop, orientation, and title lines.')
};
/** @typedef {keyof typeof CIRCULAR_EDIT_HINTS} CircularOwnRecordBlock */

/**
 * What a Legend row's popup says about its color and features.
 * @typedef {{ trackKind: DrawnTrack['kind'] | null, rules: readonly ColorRuleLike[], featureCount: number | null,
 *   defaultColorType: string | null, replacedDefaultColor: string | null }} LegendRowFacts
 */

/**
 * The rules, feature count and default color behind a Legend row; a track row has none to read.
 * @param {Pick<PreviewItemPopupOptions, 'readLegendTrackKind' | 'readLegendRowRules' | 'countLegendCaptionFeatures'
 *   | 'readLegendRowDefaultColorType' | 'readUserDefaultColor'>} ports
 * @param {LegendRow} entry @param {boolean} matchesReady
 * @returns {LegendRowFacts}
 */
const legendRowFacts = (ports, entry, matchesReady) => {
  const trackKind = ports.readLegendTrackKind(entry);
  /** @type {LegendRowFacts} */
  const facts = { trackKind, rules: [], featureCount: null, defaultColorType: null, replacedDefaultColor: null };
  if (trackKind) return facts;
  try {
    facts.rules = ports.readLegendRowRules(entry);
    if (matchesReady) [facts.featureCount, facts.defaultColorType] = [ports.countLegendCaptionFeatures(rowKey(entry)), ports.readLegendRowDefaultColorType(entry)];
    if (facts.defaultColorType) facts.replacedDefaultColor = ports.readUserDefaultColor(facts.defaultColorType);
  } catch {
    // The rule matches are still preparing: the popup shows the color alone.
  }
  return facts;
};

/** The settings a Linear record popup edits, from the record's File card. @param {LinearSeq | null} sequence @returns {RecordDraft | null} */
const linearRecordDraft = (sequence) => sequence && {
  name: text(resolveLinearRecordEffectiveDefinition(sequence)),
  subtitle: text(resolveLinearRecordEffectiveSubtitle(sequence)),
  reverse: cropFormReverseComplement({ windows: sequence.windows, reverse: sequence.region_reverse }),
  reverseAvailable: true,
  nameLabel: t('Organism / strain'),
  subtitleLabel: t('Subtitle'),
  sequenceUid: sequence.uid,
  fields: RECORD_POPUP_FIELDS.linear,
  // The record's own name fields, on its File card or in its Record options.
  settingsSelector: `[data-linear-record-names=${JSON.stringify(sequence.uid)}] input[type="text"]`
};

/**
 * The Circular settings a record popup reads: the name and subtitle fields its draft names, and the crop form.
 * @typedef {Readonly<Record<string, unknown>> & {
 *   readonly circular_windows: readonly { reverseComplement: boolean }[], readonly circular_reverse: boolean
 * }} CircularRecordForm
 */

/**
 * The settings a Circular record popup edits: the record's own fields, else the drawing's Species and Strain.
 * @param {CircularRecordForm} form the shown drawing's settings @param {boolean} ownRecord
 * @returns {RecordDraft | null}
 */
const circularRecordDraft = (form, ownRecord) => {
  const target = ownRecord ? CIRCULAR_RECORD : CIRCULAR_DRAWING;
  // A Record label drawn while its sidebar fields are unavailable (for
  // example, a loaded Session's records are not inspected yet) overrides
  // Species and Strain: the popup offers neither (OV-432).
  if (!ownRecord && (text(form.circular_record_label) || text(form.circular_record_subtitle))) return null;
  return {
    ...target,
    nameLabel: tm(target.nameLabel),
    subtitleLabel: tm(target.subtitleLabel),
    name: text(form[target.fields.name.key]),
    subtitle: text(form[target.fields.subtitle.key]),
    reverse: cropFormReverseComplement({ windows: form.circular_windows, reverse: form.circular_reverse }),
    reverseAvailable: ownRecord,
    sequenceUid: ''
  };
};

/**
 * The request record a drawn record definition shows, and the File card (by its uid) that draws it, or null.
 * @param {RecordScope} scope @param {{ recordIndex: number, recordId: string }} item
 * @param {readonly RequestRecordLike[]} records the request records of the displayed Result
 * @param {PreviewItemPopupOptions['readCommittedRecordRows']} readCommittedRecordRows
 * @returns {{ record: RequestRecordLike | null, sourceUid: string | null }}
 */
const drawnRequestRecord = (scope, { recordIndex, recordId }, records, readCommittedRecordRows) => {
  // A Linear file may draw several records (one request record, or one
  // `<uid>:<n>` record each), so the drawn record's id, not its index, names
  // its request record and File card; an ambiguous id names neither. A file's
  // record list is read on demand (a loaded Session has none yet): without
  // rows the index names the record only when every request record draws
  // exactly one, and that record's key is its File card's.
  const rows = scope === 'linear'
    ? readCommittedRecordRows().filter((row) => row.scope === 'linear' && text(row.recordId) === recordId)
    : [];
  const row = rows.length === 1 ? rows[0] : null;
  const byIndex = scope !== 'linear'
    || (!rows.length && records.every((entry) => entry.cardinality === 'exactly_one'));
  const record = (row
    ? records.find((entry) => entry.recordKey === row.canonicalRecordKey)
    : byIndex && records[recordIndex]) || null;
  return { record, sourceUid: record ? text(row ? row.sourceUid : record.recordKey) : null };
};

/**
 * The Region popup's heading; its body is the Select region owner's (`app/region-select.js`).
 * @param {RegionSelection} item
 * @returns {{ kind: 'region', identity: string, dialogLabel: string } & ReturnType<typeof regionSelectionHeading>}
 */
const regionView = (item) => ({
  kind: /** @type {const} */ ('region'),
  identity: `region:${item.resultIndex}:${item.recordKey}:${item.span.start}:${item.span.length}`,
  dialogLabel: t('Region'),
  ...regionSelectionHeading(item)
});

/**
 * The key an open Legend popup follows after a rename: the renamed row's when it kept it, else the one key the rename
 * added (a rule row's rename), else the key of the row that has the typed name (Merge); null when none is drawn.
 * @param {readonly LegendRow[]} entries the rows after the rename @param {string} renamedKey
 * @param {ReadonlySet<string>} known the keys before the rename @param {string} caption the typed name
 * @returns {string | null}
 */
const keyFollowedAfterRename = (entries, renamedKey, known, caption) => {
  const added = entries.filter((entry) => !known.has(rowKey(entry)));
  const entry = entries.find((row) => rowKey(row) === renamedKey)
    || (added.length === 1 ? added[0] : null)
    || entries.find((row) => text(row.caption) === text(caption));
  return entry ? rowKey(entry) : null;
};

// Enter in a popup text field commits it as leaving the field does and keeps
// focus in the field, so a choice the commit asks returns focus there.
/** @param {Event} event @returns {void} */
const commitField = (event) => { event.target?.dispatchEvent(new Event('change')); };

/** @param {PreviewItemPopupOptions} options */
export const createPreviewItemPopup = (options) => {
  const {
    state, computed, revealSetting, readLegendEntryColor, readShownLegendName = (caption) => caption, getCommittedRequest,
    readCommittedRecordRows, readAlignmentInspector, readCircularOwnRecordBlock, runUndoable, sessionBusy, renameLegendRow,
    legendRenameSettled
  } = options;
  /** @returns {readonly LegendRow[]} */
  const legendEntries = () => state.activeDrawing().legendEntries.value;

  /** @param {string} key @param {boolean} matchesReady */
  const legendView = (key, matchesReady) => {
    const entries = legendEntries();
    const index = entries.findIndex((entry) => rowKey(entry) === key);
    const entry = entries[index];
    if (!entry) return null;
    const facts = legendRowFacts(options, entry, matchesReady);
    const featureTypes = new Set(state.extractedFeatures.value.map((feature) => text(feature?.type)));
    return buildLegendEntryView({
      entry, name: readShownLegendName(rowKey(entry)), index, color: readLegendEntryColor(entry), ...facts, featureTypes,
      // OV-1386: a row of color none opens the Color input at its type's palette color, as the Palette pane does.
      pickerFallback: appliedFeatureColorOf(
        { paletteDefinitions: state.paletteDefinitions, appliedPaletteName: state.appliedPaletteName }, facts.defaultColorType ?? undefined
      ) || '#cccccc'
    });
  };

  /** @param {{ recordIndex: number, recordId: string }} item */
  const recordView = ({ recordIndex, recordId }) => {
    /** @type {RecordScope} */
    const scope = state.mode.value === 'circular' ? 'circular' : 'linear';
    const { record, sourceUid } = drawnRequestRecord(scope, { recordIndex, recordId }, getCommittedRequest()?.records || [], readCommittedRecordRows);
    const sequence = sourceUid !== null
      ? state.activeDrawing().linearSeqs.find((entry) => entry.uid === sourceUid) || null
      : null;
    const group = state.svgContainer.value?.querySelector(
      `[data-gbdraw-role^="record-definition"][data-gbdraw-record-index="${recordIndex}"]`
    );
    const circularBlock = scope === 'circular' ? readCircularOwnRecordBlock() : '';
    return buildRecordDefinitionView({
      scope,
      recordIndex,
      recordId,
      drawn: readDrawnDefinitionLines(group, recordId),
      record,
      draft: scope === 'linear' ? linearRecordDraft(sequence) : circularRecordDraft(state.activeDrawing().form, !circularBlock),
      editHint: scope === 'linear' ? tm(LINEAR_EDIT_HINT) : (circularBlock && tm(CIRCULAR_EDIT_HINTS[circularBlock])),
      alignment: scope === 'linear'
        ? readAlignmentInspector()?.records?.find((decision) => decision.recordKey === record?.recordKey) || null
        : null,
      alignmentShown: scope === 'linear'
    });
  };

  const view = computed(() => {
    const opened = previewItemPopupOf(state.clickedPairwiseMatch.value);
    if (!opened) return null;
    const item = opened.previewItem;
    if (item.kind === 'region') return regionView(item);
    return item.kind === 'legend-entry' ? legendView(item.key, opened.matchesReady) : recordView(item);
  });

  /** @param {RecordDefinitionView} recordView @returns {Record<string, unknown> | null} */
  const recordTarget = (recordView) => (recordView.scope === 'circular'
    ? state.activeDrawing().form
    : state.activeDrawing().linearSeqs.find((entry) => entry.uid === recordView.draft?.sequenceUid) || null);

  /**
   * Writes the setting the sidebar field writes; it applies as its setting does.
   * @param {RecordDefinitionView} recordView @param {RecordField} field @param {unknown} value
   */
  const setRecordField = (recordView, field, value) => {
    const target = recordTarget(recordView);
    const fields = recordView.draft?.fields;
    if (!target || !fields || sessionBusy()) return null;
    const write = field === 'reverse' ? recordReverseWrite(target, fields.reverse, Boolean(value))
      : recordTextWrite(target, fields[field].key, String(value ?? ''));
    return write && runUndoable(field === 'reverse' ? msg('Change Reverse complement') : msg('Edit record definition'), write);
  };

  // Shows the record's settings in the settings pane: the field the popup's
  // first text box edits, focused.
  /** @param {RecordDefinitionView} recordView */
  const showRecordSettings = (recordView) => {
    const selector = recordView.draft?.settingsSelector;
    /** @type {HTMLElement | null} */
    const field = selector ? document.querySelector(selector) : null;
    if (!field) return false;
    state.clickedPairwiseMatch.value = null;
    void revealSetting(field);
    return true;
  };

  // The label fit dialog's "Edit Species / Strain…" (PD-OI-097): the settings
  // a click on the drawn record definition with this record index would show,
  // found by that index's drawn definition group (two records may share an ID).
  /** @param {number} recordIndex */
  const showRecordTextSettings = (recordIndex) => {
    if (!Number.isSafeInteger(recordIndex) || recordIndex < 0) return false;
    const group = state.svgContainer.value?.querySelector(
      `[data-gbdraw-role="record-definition"][data-gbdraw-record-index="${recordIndex}"][data-gbdraw-record-id]`
    );
    if (!group) return false;
    return showRecordSettings(recordView({ recordIndex, recordId: String(group.getAttribute('data-gbdraw-record-id')) }));
  };

  // A rename may ask a choice first (a name another row has) and may change
  // the row's key. Once the choice is answered, the open popup follows the
  // row by its key: it stays on a row that kept its key (Cancel, a Legend-only
  // rename, or Keep both of such a row), else moves to the one key the rename
  // added (a rule row's rename), else to the row it joined, which has the
  // typed name (Merge).
  /** @param {number} index @param {string} caption */
  const renameLegendEntry = async (index, caption) => {
    const opened = previewItemPopupOf(state.clickedPairwiseMatch.value);
    const before = legendEntries();
    const renamedKey = before[index] ? rowKey(before[index]) : '';
    const known = new Set(before.map(rowKey));
    const result = await renameLegendRow(index, caption);
    await legendRenameSettled();
    if (!opened || state.clickedPairwiseMatch.value !== opened || opened.previewItem.kind !== 'legend-entry') return result;
    const entries = legendEntries();
    const key = keyFollowedAfterRename(entries, renamedKey, known, caption);
    if (key !== null) opened.previewItem = { kind: 'legend-entry', key };
    return result;
  };

  return { view, setRecordField, showRecordSettings, showRecordTextSettings, renameLegendRow: renameLegendEntry, commitField };
};
