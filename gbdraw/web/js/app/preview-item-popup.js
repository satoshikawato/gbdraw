// @ts-check
/** @import { DrawingState, LegendEntry, LinearSeq } from '../state.js' */
// The legend-entry and record-definition popups of the preview. A click (or
// Enter) on a Legend row or a record definition opens the match popup frame
// with a `previewItem`; this owner derives what the popup shows from the
// displayed Result and the current settings. It holds no edit of its own: the
// template calls the Legend panel's operations, and a record edit writes the
// setting the sidebar field writes, as one History step.
import { COMPARISON_LEGEND_SELECTOR, legendEntryKey } from '../services/legend-svg.js';
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
 *   presentation?: { label?: unknown, subtitle?: unknown, reverseComplement?: unknown } | null,
 *   display?: { startCoordinate?: unknown } | null,
 *   region?: RegionLike | null
 * }} RequestRecordLike
 * @typedef {{ recordKey?: unknown, status?: unknown, anchorLabel?: unknown, rationaleLabel?: unknown }} AlignmentDecisionLike
 * @typedef {'linear' | 'circular'} RecordScope
 * @typedef {'name' | 'subtitle' | 'reverse'} RecordField
 * @typedef {{
 *   name: string, subtitle: string, reverse: boolean, reverseAvailable: boolean,
 *   nameLabel: string, subtitleLabel: string, sequenceUid: string
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
 * @typedef {object} PreviewItemPopupOptions
 * @property {PreviewItemState} state
 * @property {<T>(getter: () => T) => { readonly value: T }} computed Vue `computed`
 * @property {(entry: LegendEntry) => string} readLegendEntryColor
 *   The color the Legend panel shows for a row.
 * @property {(entry: LegendEntry) => readonly ColorRuleLike[]} readLegendRowRules
 *   The color rules that draw a row of the displayed Result.
 * @property {(caption: string) => number} countLegendCaptionFeatures
 *   The number of features whose legend entry is `caption`.
 * @property {() => { records?: readonly RequestRecordLike[] } | null} getCommittedRequest
 *   The request of the displayed Result.
 * @property {() => { records?: readonly AlignmentDecisionLike[] } | null} readAlignmentInspector
 *   The active similarity alignment as the drawer inspects it, or null.
 * @property {() => boolean} readCircularReverseAvailable Whether the Circular reverse setting is editable.
 * @property {(label: string, apply: () => void) => unknown} runUndoable History's one undoable step.
 * @property {() => unknown} sessionBusy A reason settings cannot change now, else a falsy value.
 * @property {(index: number, caption: string) => unknown} renameLegendRow
 *   The Legend panel's rename of row `index`, as one History step.
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
    const key = text(element.getAttribute('data-legend-key'));
    return key ? { element, item: { kind: 'legend-entry', key } } : null;
  }
  const index = text(element.getAttribute('data-gbdraw-record-index'));
  const recordIndex = Number(index);
  if (!index || !Number.isInteger(recordIndex) || recordIndex < 0) return null;
  return { element, item: { kind: 'record-definition', recordIndex, recordId: text(element.getAttribute('data-gbdraw-record-id')) } };
};

/** @param {ColorRuleLike} rule */
export const colorRuleText = (rule) => {
  const feature = text(rule?.feat) || 'Any feature';
  const qualifier = text(rule?.qual);
  const pattern = text(rule?.val);
  if (!qualifier && !pattern) return `${feature} color rule`;
  return `${feature} where ${qualifier || 'a qualifier'} matches “${pattern}”`;
};

/**
 * The legend-entry popup view.
 * @param {{ entry: LegendEntry, index: number, color: string, rules: readonly ColorRuleLike[],
 *   featureCount: number | null, featureTypes: ReadonlySet<string> }} input
 */
export const buildLegendEntryView = ({ entry, index, color, rules, featureCount, featureTypes }) => {
  const caption = text(entry.caption);
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
    dialogLabel: `Legend entry: ${caption}`,
    index,
    caption,
    color,
    title: caption,
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
 *   alignment: AlignmentDecisionLike | null, alignmentShown: boolean }} input
 */
export const buildRecordDefinitionView = ({
  scope, recordIndex, recordId, drawn, record, draft, alignment, alignmentShown
}) => {
  const accession = text(drawn.accession) || text(recordId);
  const organism = text(drawn.name) || plain(record?.presentation?.label) || accession;
  const reversed = Boolean(record?.presentation?.reverseComplement);
  const start = positive(record?.display?.startCoordinate);
  /** @type {PreviewItemFact[]} */
  const facts = [
    { label: 'Accession', value: accession || '—' },
    { label: 'Organism', value: organism || '—' },
    { label: 'Length', value: text(drawn.length) || '—' },
    { label: 'Reversed', value: reversed ? 'Yes' : 'No' },
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
    title: organism,
    subtitle: [text(drawn.subtitle), accession].filter(Boolean).join(' · '),
    facts,
    draft
  };
};

/**
 * The texts of a definition group's lines by kind (`data-definition-line-kind`).
 * @param {Element | null | undefined} group @param {string} recordId
 */
export const readDrawnDefinitionLines = (group, recordId) => {
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

/** @type {Record<RecordScope, Record<RecordField, string>>} */
const RECORD_FIELDS = {
  linear: { name: 'definition', subtitle: 'record_subtitle', reverse: 'region_reverse' },
  circular: { name: 'species', subtitle: 'strain', reverse: 'circular_reverse' }
};

/** @param {PreviewItemPopupOptions} options */
export const createPreviewItemPopup = ({
  state, computed, readLegendEntryColor, readLegendRowRules, countLegendCaptionFeatures,
  getCommittedRequest, readAlignmentInspector, readCircularReverseAvailable, runUndoable, sessionBusy, renameLegendRow
}) => {
  /** @param {string} key @param {boolean} matchesReady */
  const legendView = (key, matchesReady) => {
    /** @type {LegendEntry[]} */
    const entries = state.activeDrawing().legendEntries.value || [];
    const index = entries.findIndex((entry) => text(legendEntryKey(entry)) === key || text(entry.caption) === key);
    const entry = entries[index];
    if (!entry) return null;
    /** @type {readonly ColorRuleLike[]} */
    let rules = [];
    /** @type {number | null} */
    let featureCount = null;
    try {
      rules = readLegendRowRules(entry);
      if (matchesReady) featureCount = countLegendCaptionFeatures(text(entry.caption));
    } catch {
      // The rule matches are still preparing: the popup shows the color alone.
    }
    const featureTypes = new Set(state.extractedFeatures.value.map((feature) => text(feature?.type)));
    return buildLegendEntryView({ entry, index, color: readLegendEntryColor(entry), rules, featureCount, featureTypes });
  };

  /** @param {LinearSeq | null} sequence @returns {RecordDraft | null} */
  const linearDraft = (sequence) => sequence && {
    name: text(resolveLinearRecordEffectiveDefinition(sequence)),
    subtitle: text(resolveLinearRecordEffectiveSubtitle(sequence)),
    reverse: Boolean(sequence.region_reverse),
    reverseAvailable: true,
    nameLabel: 'Organism / strain',
    subtitleLabel: 'Subtitle',
    sequenceUid: sequence.uid
  };

  /** @returns {RecordDraft} */
  const circularDraft = () => {
    const form = state.activeDrawing().form;
    return {
      name: text(form.species),
      subtitle: text(form.strain),
      reverse: Boolean(form.circular_reverse),
      reverseAvailable: readCircularReverseAvailable(),
      nameLabel: 'Species',
      subtitleLabel: 'Strain',
      sequenceUid: ''
    };
  };

  /** @param {{ recordIndex: number, recordId: string }} item */
  const recordView = ({ recordIndex, recordId }) => {
    /** @type {RecordScope} */
    const scope = state.mode.value === 'circular' ? 'circular' : 'linear';
    const record = getCommittedRequest()?.records?.[recordIndex] || null;
    const sequence = scope === 'linear'
      ? state.linearSeqs.find((entry) => entry.uid === text(record?.recordKey)) || null
      : null;
    const group = state.svgContainer.value?.querySelector(
      `[data-gbdraw-role^="record-definition"][data-gbdraw-record-index="${recordIndex}"]`
    );
    return buildRecordDefinitionView({
      scope,
      recordIndex,
      recordId,
      drawn: readDrawnDefinitionLines(group, recordId),
      record,
      draft: scope === 'linear' ? linearDraft(sequence) : circularDraft(),
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
    const key = RECORD_FIELDS[recordView.scope][field];
    const next = field === 'reverse' ? Boolean(value) : String(value ?? '');
    if (target[key] === next) return null;
    return runUndoable(field === 'reverse' ? 'Change record orientation' : 'Edit record definition', () => {
      target[key] = next;
    });
  };

  // Opens the record's settings in the sidebar and focuses its first field.
  /** @param {RecordDefinitionView} recordView */
  const showRecordSettings = (recordView) => {
    const uid = recordView.draft?.sequenceUid;
    /** @type {HTMLElement | null} */
    let field = null;
    if (recordView.scope === 'linear' && uid) {
      const card = document.querySelector(`[data-linear-record-card="${CSS.escape(uid)}"]`);
      card?.closest('details')?.setAttribute('open', '');
      const options = card?.querySelector('[data-linear-record-options]');
      options?.setAttribute('open', '');
      field = options?.querySelector('input[type="text"]') || null;
    } else if (recordView.scope === 'circular') {
      field = document.getElementById('circular-species');
    }
    if (!field) return false;
    state.clickedPairwiseMatch.value = null;
    field.scrollIntoView({ block: 'center' });
    field.focus({ preventScroll: true });
    return true;
  };

  // A rename can change the row's key; the open popup then follows the row by
  // its new caption (after a Merge or Suffix choice too, when one row has it).
  /** @param {number} index @param {string} caption */
  const renameLegendEntry = async (index, caption) => {
    const result = await renameLegendRow(index, caption);
    const open = state.clickedPairwiseMatch.value;
    if (open?.previewItem?.kind !== 'legend-entry') return result;
    /** @type {LegendEntry[]} */
    const entries = state.activeDrawing().legendEntries.value || [];
    const named = entries.filter((entry) => text(entry.caption) === text(caption));
    const entry = named.length === 1 ? named[0] : entries[index];
    if (entry) open.previewItem = { kind: 'legend-entry', key: text(legendEntryKey(entry)) };
    return result;
  };

  return { view, setRecordField, showRecordSettings, renameLegendRow: renameLegendEntry };
};
