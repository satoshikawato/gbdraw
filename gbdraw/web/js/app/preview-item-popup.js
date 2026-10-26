// @ts-check
/** @import { DrawingState, LinearSeq } from '../state.js' */
/** @import { LegendRow } from '../services/drawing/paint.js' */
/** @import { DrawnTrack } from '../services/track-legend-rows.js' */
// The legend-entry and record-definition popups of the preview. A click (or
// Enter) on a Legend row or a record definition opens the match popup frame
// with a `previewItem`; this owner derives what the popup shows from the
// displayed Result and the current settings. It holds no edit of its own: the
// template calls the Legend panel's operations, and a record edit writes the
// setting the sidebar field writes, as one History step.
import { COMPARISON_LEGEND_SELECTOR, pythonLegendKey } from '../services/legend-svg.js';
import { legendEntryKey } from '../services/drawing/paint.js';
import { plainTextLinearRecordLabel } from '../services/linear-comparisons.js';
import {
  resolveLinearRecordEffectiveDefinition,
  resolveLinearRecordEffectiveSubtitle
} from '../services/linear-sources.js';
import { appliedFeatureColorOf, colorSwatchView } from '../utils/color-utils.js';
import { msg, t, tn, tm } from '../utils/i18n.js';

/**
 * The item a preview popup is about: a Legend row by its key, or the
 * definition of the Result's record at `recordIndex`.
 * @typedef {{ kind: 'legend-entry', key: string }
 *   | { kind: 'record-definition', recordIndex: number, recordId: string }} PreviewItem
 * @typedef {{ previewItem?: PreviewItem, matchesReady?: boolean }} PreviewItemPayload
 * @typedef {{ label: string, value: string }} PreviewItemFact
 * @typedef {{ feat?: unknown, qual?: unknown, val?: unknown }} ColorRuleLike
 * @typedef {{ start?: unknown, end?: unknown }} RegionLike
 * @typedef {{
 *   recordKey?: unknown,
 *   cardinality?: unknown,
 *   presentation?: { label?: unknown, subtitle?: unknown, reverseComplement?: unknown } | null,
 *   display?: { startCoordinate?: unknown } | null,
 *   region?: RegionLike | null
 * }} RequestRecordLike
 * @typedef {{ recordKey?: unknown, status?: unknown, anchorLabel?: unknown, rationaleLabel?: unknown }} AlignmentDecisionLike
 * @typedef {{ scope?: unknown, sourceUid?: unknown, recordId?: unknown, canonicalRecordKey?: unknown }} CommittedRecordRowLike
 * @typedef {'linear' | 'circular'} RecordScope
 * @typedef {'name' | 'subtitle' | 'reverse'} RecordField
 * @typedef {{
 *   name: string, subtitle: string, reverse: boolean, reverseAvailable: boolean,
 *   nameLabel: string, subtitleLabel: string, sequenceUid: string,
 *   fields: Record<RecordField, string>, settingsSelector: string
 * }} RecordDraft
 * @typedef {ReturnType<typeof buildLegendEntryView>} LegendEntryView
 * @typedef {ReturnType<typeof buildRecordDefinitionView>} RecordDefinitionView
 */

/**
 * The members of the Web state this owner reads.
 * @typedef {object} PreviewItemState
 * @property {{ value: PreviewItemPayload | null }} clickedPairwiseMatch The open match-frame popup.
 * @property {{ value: Element | null }} svgContainer The preview's mounted Result.
 * @property {{ value: string }} mode
 * @property {LinearSeq[]} linearSeqs
 * @property {{ value: ReadonlyArray<{ type?: unknown }> }} extractedFeatures
 * @property {() => DrawingState} activeDrawing
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

const NUMBER = new Intl.NumberFormat('en-US');
/** @param {unknown} value */
const text = (value) => String(value ?? '').trim();
/** @param {unknown} value */
const plain = (value) => (text(value) ? plainTextLinearRecordLabel(text(value)) : '');
/** @param {unknown} value */
const positive = (value) => {
  const number = Number(value);
  return value != null && Number.isFinite(number) && number > 0 ? number : null;
};

export const PREVIEW_ITEM_SELECTOR =
  'g[data-legend-key], [data-gbdraw-role="record-definition"], [data-gbdraw-role="record-definition-row"]';

/**
 * The Legend row or record definition an event in the preview reached.
 * @param {EventTarget | null} target @param {Element} svg
 * @returns {{ element: Element, item: PreviewItem } | null}
 */
export const previewItemTarget = (target, svg) => {
  const element = target instanceof Element ? target.closest(PREVIEW_ITEM_SELECTOR) : null;
  if (!element || !svg.contains(element) || element.closest(COMPARISON_LEGEND_SELECTOR)) return null;
  if (element.hasAttribute('data-legend-key')) {
    const key = text(pythonLegendKey(element));
    return key ? { element, item: { kind: 'legend-entry', key } } : null;
  }
  const index = text(element.getAttribute('data-gbdraw-record-index'));
  const recordIndex = Number(index);
  if (!index || !Number.isInteger(recordIndex) || recordIndex < 0) return null;
  return { element, item: { kind: 'record-definition', recordIndex, recordId: text(element.getAttribute('data-gbdraw-record-id')) } };
};

/**
 * OV-1466: a rule without a feature type is a whole sentence of its own, so a
 * locale writes "any feature" where and how its sentence needs it.
 * @param {ColorRuleLike} rule
 */
const colorRuleText = (rule) => {
  const feature = text(rule?.feat);
  const qualifier = text(rule?.qual);
  const pattern = text(rule?.val);
  if (!qualifier && !pattern) return feature ? t('{feature} color rule', { feature }) : t('Any feature color rule');
  const params = { qualifier: qualifier || t('a qualifier'), pattern };
  return feature
    ? t('{feature} where {qualifier} matches “{pattern}”', { feature, ...params })
    : t('Any feature where {qualifier} matches “{pattern}”', params);
};

// OV-1272: a track row names its plot; it lists no features, and its Color
// edits only the swatch (the plot's color is a Settings value).
/** @type {Readonly<Record<DrawnTrack['kind'], string>>} */
const TRACK_ROW_SUBTITLE = Object.freeze({ content: msg('Content plot'), skew: msg('Skew plot'), depth: msg('Depth plot') });
const TRACK_ROW_COLOR_SOURCE = msg('Changes this legend swatch only. The plot keeps the color set in Settings.');

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
  const key = text(legendEntryKey(entry));
  const ruleLines = trackKind ? [] : rules.map(colorRuleText);
  const defaultColorLine = trackKind || ruleLines.length ? null : defaultColorScopeLine(defaultColorType, replacedDefaultColor);
  // OV-1336: the OIC-029 line, once shown, says where the color comes from by itself; the caption gives way.
  const colorSource = defaultColorLine
    ? null
    : trackKind
      ? tm(TRACK_ROW_COLOR_SOURCE)
      : ruleLines.length
        ? tn(ruleLines.length, 'From this color rule:', 'From these {count} color rules:')
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
        : tn(featureCount, '{count} feature', '{count} features', { count: NUMBER.format(featureCount) }),
    colorSource,
    defaultColorLine,
    ruleLines
  };
};

/** @param {RegionLike | null | undefined} region */
const regionText = (region) => {
  const start = positive(region?.start);
  const end = positive(region?.end);
  if (start == null && end == null) return t('Whole record');
  return `${start == null ? '1' : NUMBER.format(start)}..${end == null ? 'end' : NUMBER.format(end)}`;
};

/** @param {AlignmentDecisionLike | null | undefined} decision */
const alignmentText = (decision) => {
  if (!decision) return t('Not aligned');
  if (decision.status === 'reference') return t('Reference ({anchor})', { anchor: text(decision.anchorLabel) });
  if (decision.status === 'aligned') return t('Aligned on {anchor}', { anchor: text(decision.anchorLabel) });
  return t('Not moved ({reason})', { reason: text(decision.rationaleLabel).toLowerCase() || t('skipped') });
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
  const reversed = Boolean(record?.presentation?.reverseComplement);
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
    { label: t('Display start'), value: start == null ? t('Record start') : NUMBER.format(start) },
    { label: t('Region'), value: regionText(record?.region) }
  ];
  if (alignmentShown) facts.push({ label: t('Alignment'), value: alignmentText(alignment) });
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
// definition. A Circular record selected on its own draws its Record label and
// Subtitle, which win over the drawing's Species and Strain.
const LINEAR_FIELDS = { name: 'definition', subtitle: 'record_subtitle', reverse: 'region_reverse' };
const CIRCULAR_RECORD = {
  fields: { name: 'circular_record_label', subtitle: 'circular_record_subtitle', reverse: 'circular_reverse' },
  nameLabel: msg('Record label'), subtitleLabel: msg('Subtitle'),
  settingsSelector: '[data-circular-record-presentation] input[type="text"]'
};
const CIRCULAR_DRAWING = {
  fields: { name: 'species', subtitle: 'strain', reverse: 'circular_reverse' },
  nameLabel: msg('Species'), subtitleLabel: msg('Strain'), settingsSelector: '#circular-species'
};
// What a popup without edit fields says instead. A Circular record's sentence
// names the cause; the second and third are the sidebar's for the same cause.
const LINEAR_EDIT_HINT = msg('Edit this record in its File card in the sidebar.');
const CIRCULAR_EDIT_HINTS = {
  'not-inspected': msg('Inspect source records in the sidebar to edit this record\'s label.'),
  'multi-record-canvas': msg('Multi-Record Canvas uses a grid. Turn it off to select one record.'),
  'no-single-record': msg('Select one record to edit crop, orientation, and title lines.')
};
/** @typedef {keyof typeof CIRCULAR_EDIT_HINTS} CircularOwnRecordBlock */

/** @param {PreviewItemPopupOptions} options */
export const createPreviewItemPopup = ({
  state, computed, revealSetting, readLegendEntryColor, readShownLegendName = (caption) => caption, readLegendRowRules, countLegendCaptionFeatures,
  readLegendTrackKind, readLegendRowDefaultColorType, readUserDefaultColor, getCommittedRequest, readCommittedRecordRows, readAlignmentInspector, readCircularOwnRecordBlock, runUndoable, sessionBusy, renameLegendRow,
  legendRenameSettled
}) => {
  /** @returns {readonly LegendRow[]} */
  const legendEntries = () => state.activeDrawing().legendEntries.value || [];

  /** @param {string} key @param {boolean} matchesReady */
  const legendView = (key, matchesReady) => {
    const entries = legendEntries();
    const index = entries.findIndex((entry) => text(legendEntryKey(entry)) === key);
    const entry = entries[index];
    if (!entry) return null;
    const trackKind = readLegendTrackKind(entry);
    /** @type {readonly ColorRuleLike[]} */
    let rules = [];
    /** @type {number | null} */
    let featureCount = null;
    /** @type {string | null} */
    let defaultColorType = null;
    /** @type {string | null} */
    let replacedDefaultColor = null;
    // A track row has no rules and no features to read.
    if (!trackKind) {
      try {
        rules = readLegendRowRules(entry);
        if (matchesReady) [featureCount, defaultColorType] = [countLegendCaptionFeatures(text(legendEntryKey(entry))), readLegendRowDefaultColorType(entry)];
        if (defaultColorType) replacedDefaultColor = readUserDefaultColor(defaultColorType);
      } catch {
        // The rule matches are still preparing: the popup shows the color alone.
      }
    }
    const featureTypes = new Set(state.extractedFeatures.value.map((feature) => text(feature?.type)));
    return buildLegendEntryView({
      entry, name: readShownLegendName(text(legendEntryKey(entry))), index, color: readLegendEntryColor(entry), rules, featureCount, featureTypes,
      trackKind, defaultColorType, replacedDefaultColor,
      // OV-1386: a row of color none opens the Color input at its type's palette color, as the Palette pane does.
      pickerFallback: appliedFeatureColorOf(
        { paletteDefinitions: state.paletteDefinitions, appliedPaletteName: state.appliedPaletteName }, defaultColorType ?? undefined
      ) || '#cccccc'
    });
  };

  /** @param {LinearSeq | null} sequence @returns {RecordDraft | null} */
  const linearDraft = (sequence) => sequence && {
    name: text(resolveLinearRecordEffectiveDefinition(sequence)),
    subtitle: text(resolveLinearRecordEffectiveSubtitle(sequence)),
    reverse: Boolean(sequence.region_reverse),
    reverseAvailable: true,
    nameLabel: t('Organism / strain'),
    subtitleLabel: t('Subtitle'),
    sequenceUid: sequence.uid,
    fields: LINEAR_FIELDS,
    // The record's own name fields, on its File card or in its Record options.
    settingsSelector: `[data-linear-record-names=${JSON.stringify(sequence.uid)}] input[type="text"]`
  };

  /** @param {boolean} ownRecord @returns {RecordDraft | null} */
  const circularDraft = (ownRecord) => {
    const form = state.activeDrawing().form;
    const target = ownRecord ? CIRCULAR_RECORD : CIRCULAR_DRAWING;
    // A Record label drawn while its sidebar fields are unavailable (for
    // example, a loaded Session's records are not inspected yet) overrides
    // Species and Strain: the popup offers neither (OV-432).
    if (!ownRecord && (text(form.circular_record_label) || text(form.circular_record_subtitle))) return null;
    return {
      ...target,
      nameLabel: tm(target.nameLabel),
      subtitleLabel: tm(target.subtitleLabel),
      name: text(form[target.fields.name]),
      subtitle: text(form[target.fields.subtitle]),
      reverse: Boolean(form.circular_reverse),
      reverseAvailable: ownRecord,
      sequenceUid: ''
    };
  };

  /** @param {{ recordIndex: number, recordId: string }} item */
  const recordView = ({ recordIndex, recordId }) => {
    /** @type {RecordScope} */
    const scope = state.mode.value === 'circular' ? 'circular' : 'linear';
    const records = getCommittedRequest()?.records || [];
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
    const sequence = record
      ? state.linearSeqs.find((entry) => entry.uid === text(row ? row.sourceUid : record.recordKey)) || null
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
      draft: scope === 'linear' ? linearDraft(sequence) : circularDraft(!circularBlock),
      editHint: scope === 'linear' ? tm(LINEAR_EDIT_HINT) : (circularBlock && tm(CIRCULAR_EDIT_HINTS[circularBlock])),
      alignment: scope === 'linear'
        ? readAlignmentInspector()?.records?.find((decision) => decision.recordKey === record?.recordKey) || null
        : null,
      alignmentShown: scope === 'linear'
    });
  };

  const view = computed(() => {
    const opened = state.clickedPairwiseMatch.value;
    const item = opened?.previewItem;
    if (!item) return null;
    return item.kind === 'legend-entry' ? legendView(item.key, Boolean(opened?.matchesReady)) : recordView(item);
  });

  /** @param {RecordDefinitionView} recordView @returns {Record<string, unknown> | null} */
  const recordTarget = (recordView) => (recordView.scope === 'circular'
    ? state.activeDrawing().form
    : state.linearSeqs.find((entry) => entry.uid === recordView.draft?.sequenceUid) || null);

  /**
   * Writes the setting the sidebar field writes; it applies on Generate.
   * @param {RecordDefinitionView} recordView @param {RecordField} field @param {unknown} value
   */
  const setRecordField = (recordView, field, value) => {
    const target = recordTarget(recordView);
    if (!target || sessionBusy()) return null;
    const key = recordView.draft?.fields[field];
    if (!key) return null;
    const next = field === 'reverse' ? Boolean(value) : String(value ?? '');
    if (target[key] === next) return null;
    return runUndoable(field === 'reverse' ? msg('Change Reverse complement') : msg('Edit record definition'), () => {
      target[key] = next;
    });
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
    const opened = state.clickedPairwiseMatch.value;
    const before = legendEntries();
    const renamedKey = before[index] ? text(legendEntryKey(before[index])) : '';
    const known = new Set(before.map((entry) => text(legendEntryKey(entry))));
    const result = await renameLegendRow(index, caption);
    await legendRenameSettled();
    if (state.clickedPairwiseMatch.value !== opened || opened?.previewItem?.kind !== 'legend-entry') return result;
    const entries = legendEntries();
    const added = entries.filter((entry) => !known.has(text(legendEntryKey(entry))));
    const entry = entries.find((row) => text(legendEntryKey(row)) === renamedKey)
      || (added.length === 1 ? added[0] : null)
      || entries.find((row) => text(row.caption) === text(caption));
    if (entry) opened.previewItem = { kind: 'legend-entry', key: text(legendEntryKey(entry)) };
    return result;
  };

  // Enter in a popup text field commits it as leaving the field does and keeps
  // focus in the field, so a choice the commit asks returns focus there.
  /** @param {Event} event */
  const commitField = (event) => { event.target?.dispatchEvent(new Event('change')); };

  return { view, setRecordField, showRecordSettings, showRecordTextSettings, renameLegendRow: renameLegendEntry, commitField };
};
