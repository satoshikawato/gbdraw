// @ts-check
/** @import { DrawingState, LinearSeq } from '../state.js' */
/** @import { LegendRow } from '../services/drawing/paint.js' */
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

/** @param {ColorRuleLike} rule */
const colorRuleText = (rule) => {
  const feature = text(rule?.feat) || 'Any feature';
  const qualifier = text(rule?.qual);
  const pattern = text(rule?.val);
  if (!qualifier && !pattern) return `${feature} color rule`;
  return `${feature} where ${qualifier || 'a qualifier'} matches “${pattern}”`;
};

/**
 * The legend-entry popup view.
 * @param {{ entry: LegendRow, name: string, index: number, color: string, rules: readonly ColorRuleLike[],
 *   featureCount: number | null, featureTypes: ReadonlySet<string> }} input `name`: the name the Legend draws.
 */
const buildLegendEntryView = ({ entry, name, index, color, rules, featureCount, featureTypes }) => {
  const key = text(legendEntryKey(entry));
  const ruleLines = rules.map(colorRuleText);
  const colorSource = ruleLines.length
    ? (ruleLines.length === 1 ? 'From this color rule:' : `From these ${ruleLines.length} color rules:`)
    : featureTypes.has(key)
      ? `Default color of ${key} features`
      : 'Set on individual features';
  return {
    kind: /** @type {const} */ ('legend-entry'),
    identity: `legend:${key}`,
    dialogLabel: `Legend item: ${name}`,
    index,
    name,
    color,
    title: name,
    subtitle: featureCount == null
      ? 'Counting features…'
      : `${NUMBER.format(featureCount)} feature${featureCount === 1 ? '' : 's'}`,
    colorSource,
    ruleLines
  };
};

/** @param {RegionLike | null | undefined} region */
const regionText = (region) => {
  const start = positive(region?.start);
  const end = positive(region?.end);
  if (start == null && end == null) return 'Whole record';
  return `${start == null ? '1' : NUMBER.format(start)}..${end == null ? 'end' : NUMBER.format(end)}`;
};

/** @param {AlignmentDecisionLike | null | undefined} decision */
const alignmentText = (decision) => {
  if (!decision) return 'Not aligned';
  if (decision.status === 'reference') return `Reference (${text(decision.anchorLabel)})`;
  if (decision.status === 'aligned') return `Aligned on ${text(decision.anchorLabel)}`;
  return `Not moved (${text(decision.rationaleLabel).toLowerCase() || 'skipped'})`;
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
    { label: 'Accession', value: accession || '—' },
    // OV-502: only a Linear definition marks which line is the record's name.
    // A Circular definition's name lines may be an organism, a strain, a
    // replicon or a Record label, wrapped or not: they are the heading only.
    ...(scope === 'linear' && named ? [{ label: 'Organism', value: name }] : []),
    { label: 'Length', value: text(drawn.length) || '—' },
    { label: 'Reverse complement', value: reversed ? 'Yes' : 'No' },
    { label: 'Display start', value: start == null ? 'Record start' : NUMBER.format(start) },
    { label: 'Region', value: regionText(record?.region) }
  ];
  if (alignmentShown) facts.push({ label: 'Alignment', value: alignmentText(alignment) });
  return {
    kind: /** @type {const} */ ('record-definition'),
    identity: `record:${scope}:${recordIndex}:${recordId}`,
    dialogLabel: `Record ${accession}`,
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
  nameLabel: 'Record label', subtitleLabel: 'Subtitle',
  settingsSelector: '[data-circular-record-presentation] input[type="text"]'
};
const CIRCULAR_DRAWING = {
  fields: { name: 'species', subtitle: 'strain', reverse: 'circular_reverse' },
  nameLabel: 'Species', subtitleLabel: 'Strain', settingsSelector: '#circular-species'
};
// What a popup without edit fields says instead. A Circular record's sentence
// names the cause; the second and third are the sidebar's for the same cause.
const LINEAR_EDIT_HINT = 'Edit this record in its File card in the sidebar.';
const CIRCULAR_EDIT_HINTS = {
  'not-inspected': 'Inspect source records in the sidebar to edit this record\'s label.',
  'multi-record-canvas': 'Multi-Record Canvas uses a grid. Turn it off to select one record.',
  'no-single-record': 'Select one record to edit crop, orientation, and title lines.'
};
/** @typedef {keyof typeof CIRCULAR_EDIT_HINTS} CircularOwnRecordBlock */

/** @param {PreviewItemPopupOptions} options */
export const createPreviewItemPopup = ({
  state, computed, revealSetting, readLegendEntryColor, readShownLegendName = (caption) => caption, readLegendRowRules, countLegendCaptionFeatures,
  getCommittedRequest, readCommittedRecordRows, readAlignmentInspector, readCircularOwnRecordBlock, runUndoable, sessionBusy, renameLegendRow,
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
    /** @type {readonly ColorRuleLike[]} */
    let rules = [];
    /** @type {number | null} */
    let featureCount = null;
    try {
      rules = readLegendRowRules(entry);
      if (matchesReady) featureCount = countLegendCaptionFeatures(text(legendEntryKey(entry)));
    } catch {
      // The rule matches are still preparing: the popup shows the color alone.
    }
    const featureTypes = new Set(state.extractedFeatures.value.map((feature) => text(feature?.type)));
    return buildLegendEntryView({
      entry, name: readShownLegendName(text(legendEntryKey(entry))), index, color: readLegendEntryColor(entry), rules, featureCount, featureTypes
    });
  };

  /** @param {LinearSeq | null} sequence @returns {RecordDraft | null} */
  const linearDraft = (sequence) => sequence && {
    name: text(resolveLinearRecordEffectiveDefinition(sequence)),
    subtitle: text(resolveLinearRecordEffectiveSubtitle(sequence)),
    reverse: Boolean(sequence.region_reverse),
    reverseAvailable: true,
    nameLabel: 'Organism / strain',
    subtitleLabel: 'Subtitle',
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
      editHint: scope === 'linear' ? LINEAR_EDIT_HINT : (circularBlock && CIRCULAR_EDIT_HINTS[circularBlock]),
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
    return runUndoable(field === 'reverse' ? 'Change Reverse complement' : 'Edit record definition', () => {
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

  return { view, setRecordField, showRecordSettings, renameLegendRow: renameLegendEntry, commitField };
};
