// @ts-check
// The proposal of Find in other files (PD-OI-123) as the dialog shows it. Python made every value (the windows, which
// are ticked at first and why the others are not, the genes found, the window numbers of the drawing, the alignments
// removed); these functions only word them. A window is found by its identity in the plan: its File and its column.
import { msg, t, tm, tn } from '../../utils/i18n.js';
import { diagnosticError } from '../../utils/error-normalization.js';
import { LOCATE_OPERATION, circledNumber, windowGeneSpan, windowSpanText } from '../../services/find-in-other-files.js';

/**
 * @import { FilePlanJson, LocateDrawingReady, LocateFileCountsReply, LocatePlanJson, LocatePlanReply,
 *   LocateWindowChoiceInput, LocateWindowRefused, PlannedWindowJson } from '../../contracts/locate.generated.js'
 */

/**
 * A window of the search, as the user keeps or edits it: ticked or not, and its span and orientation.
 * @typedef {object} WindowDraft
 * @property {string} key `fileKey#column`, the window's identity in the plan.
 * @property {string} fileKey
 * @property {number} column
 * @property {number} start
 * @property {number} end
 * @property {boolean} reverseComplement
 * @property {boolean} ticked
 * @property {boolean} edited Whether the user changed the span or orientation the search gave.
 */

/**
 * One window row: display values only.
 * @typedef {object} ProposalWindowView
 * @property {string} key
 * @property {boolean} ticked
 * @property {string} number The circled number it has in the drawing, or '' when it is not ticked.
 * @property {string} recordId
 * @property {string} span With thousands separators, and "(−)" when reverse-complemented.
 * @property {string} genes Its number of genes with the first and last gene.
 * @property {string[]} ends "contig start" and "contig end" where it reaches the end of a linear record.
 * @property {string} reason Why the search did not tick it, or ''.
 */

/**
 * One File of the proposal: display values only.
 * @typedef {object} ProposalFileView
 * @property {string} fileKey
 * @property {string} name
 * @property {string} summary Its number of windows and "N of M genes".
 * @property {ProposalWindowView[]} windows The windows the search ticked, in the drawing's order.
 * @property {ProposalWindowView[]} notSelected
 * @property {string} found How many genes were found in full, in part and not at all.
 * @property {string} partial The genes found in part, or ''.
 * @property {string} absent The genes not found, or ''.
 * @property {string} empty Why the File has no window (no window, or its search failed), or ''.
 */

/** The key of a window in the plan. @param {string} fileKey @param {number | null} column */
export const windowKey = (fileKey, column) => `${fileKey}#${column ?? 0}`;

/**
 * The windows of a search as drafts: the windows it ticked, then the others, each as the search gave it.
 * @param {LocatePlanJson} plan
 * @returns {WindowDraft[]}
 */
export const windowDrafts = (plan) => plan.files.flatMap((file) => [
  ...file.windows.map((item) => [item, true]), ...file.candidates.map((item) => [item, false])
].map(([item, ticked]) => {
  const { window } = /** @type {PlannedWindowJson} */ (item);
  const column = window.column ?? 0;
  return {
    key: windowKey(file.file_key, column), fileKey: file.file_key, column, start: window.start, end: window.end,
    reverseComplement: window.reverse_complement, ticked: Boolean(ticked), edited: false
  };
}));

/** The choices of the ticked drafts, in the plan's order. @param {readonly WindowDraft[]} drafts @returns {LocateWindowChoiceInput[]} */
export const windowChoices = (drafts) => drafts.filter((draft) => draft.ticked).map(({ fileKey, column, start, end, reverseComplement }) => ({
  fileKey, column, start, end, reverseComplement
}));

/** @type {Readonly<Record<NonNullable<PlannedWindowJson['reason']>, string>>} */
const NOT_SELECTED_TEXT = Object.freeze({
  COVERED_BY_OTHER_WINDOW: msg('Not selected: adds no gene that the selected windows do not cover.'),
  TOO_FEW_GENES: msg('Not selected: adds too few genes.'),
  MAX_WINDOWS: msg('Not selected: over the limit of {count} windows of a file.')
});

/**
 * What a window edit or reference window breaks, by Python's reason, with the window's span.
 * @param {LocateWindowRefused['reason']} reason @param {number} start @param {number} end
 */
export const windowRefusalText = (reason, start, end) => {
  const span = windowSpanText({ start, end, reverse_complement: false });
  return {
    WINDOW_WITHOUT_CDS: () => t('The window {span} holds no CDS.', { span }),
    WINDOWS_OVERLAP: () => t('Two windows of the record overlap at {span}.', { span }),
    WINDOW_PAST_END: () => t('The window {span} reaches past the end of the record.', { span }),
    WINDOW_ACROSS_LINEAR_ORIGIN: () => t('The window {span} crosses the origin, which a linear record does not have.', { span })
  }[reason]();
};

/**
 * @param {PlannedWindowJson} item The window as the search gave it.
 * @param {PlannedWindowJson | undefined} chosen The same window in the plan of the kept windows, when it is ticked.
 * @param {WindowDraft} draft
 * @param {number} number Its number in the drawing, or 0.
 * @param {Readonly<Record<string, string>>} names
 * @param {number} maxWindows
 * @returns {ProposalWindowView}
 */
const windowView = (item, chosen, draft, number, names, maxWindows) => {
  // A kept window shows the plan's span and ends; one the user edited and did not keep shows the edit, without ends
  // (Python works them out when it is kept).
  const shown = chosen ?? item;
  const span = chosen ? chosen.window
    : { start: draft.start, end: draft.end, reverse_complement: draft.reverseComplement };
  const genes = windowGeneSpan(item.genes, names);
  const ends = chosen || !draft.edited
    ? [shown.touches_start ? t('contig start') : '', shown.touches_end ? t('contig end') : ''].filter(Boolean) : [];
  return {
    key: draft.key, ticked: draft.ticked, number: number ? circledNumber(number) : '', recordId: item.record_id,
    span: windowSpanText(span),
    genes: genes.first === genes.last ? t('{count} · {first}', { count: genes.count, first: genes.first })
      : t('{count} · {first} to {last}', genes),
    ends,
    reason: item.reason ? tm(NOT_SELECTED_TEXT[item.reason], { count: maxWindows }) : ''
  };
};

/** @param {readonly string[]} ids @param {Readonly<Record<string, string>>} names */
const geneList = (ids, names) => ids.map((id) => names[id] ?? id).join(', ');

/**
 * Why a File has no window: its search failed (with the reason), or no window was found (with how many alignments
 * the thresholds removed, Python's count).
 * @param {FilePlanJson} file
 * @param {LocateFileCountsReply} counts
 */
const emptyText = (file, counts) => {
  const failed = file.warnings.find((warning) => warning.code === 'FILE_SEARCH_FAILED');
  if (failed) return t('The search of this file failed: {reason}', { reason: String(failed.args.error ?? '') });
  if (file.windows.length || file.candidates.length) return '';
  return tn(counts.alignmentsRemoved, 'No window found. The thresholds removed {count} alignment.',
    'No window found. The thresholds removed {count} alignments.');
};

/**
 * One File of a plan by its File key; a plan without it is not the plan of this search.
 * @template {{ file_key: string } | { fileKey: string }} T
 * @param {readonly T[]} items @param {string} fileKey
 * @returns {T}
 */
const byFileKey = (items, fileKey) => {
  const found = items.find((item) => ('file_key' in item ? item.file_key : item.fileKey) === fileKey);
  if (!found) throw diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation: LOCATE_OPERATION });
  return found;
};

/**
 * The proposal: each File in the drawing's order, with the windows of the search and their state now. Python gives
 * every number: the window numbers of the drawing, the genes found and the alignments removed (`counts`).
 * @param {LocatePlanReply} search The search's plan, its gene names.
 * @param {Pick<LocateDrawingReady, 'plan' | 'counts'>} chosen The plan of the kept windows (Python's
 *   `choose_locate_windows`) and the counts of its Files (`locate_file_counts`).
 * @param {readonly WindowDraft[]} drafts
 * @param {Readonly<Record<string, string>>} fileNames By File key.
 * @returns {ProposalFileView[]}
 */
export const proposalView = (search, chosen, drafts, fileNames) => {
  const names = search.geneNames;
  const draftOf = new Map(drafts.map((draft) => [draft.key, draft]));
  return search.plan.files.map((file) => {
    const kept = byFileKey(chosen.plan.files, file.file_key);
    const counts = byFileKey(chosen.counts, file.file_key);
    const numbers = new Map(counts.windowNumbers.map(({ column, number }) => [column, number]));
    const keptByKey = new Map(kept.windows.map((item) => [windowKey(file.file_key, item.window.column), item]));
    /** @param {PlannedWindowJson} item */
    const row = (item) => {
      const key = windowKey(file.file_key, item.window.column);
      const keptItem = keptByKey.get(key);
      return windowView(item, keptItem, /** @type {WindowDraft} */ (draftOf.get(key)),
        keptItem ? numbers.get(item.window.column ?? 0) ?? 0 : 0, names, search.plan.options.max_windows);
    };
    const { complete, partial, absent } = kept.coverage;
    return {
      fileKey: file.file_key,
      name: fileNames[file.file_key] ?? file.file_key,
      summary: t('{windows} · {found} of {total} genes', {
        windows: tn(kept.windows.length, '{count} window', '{count} windows'), found: counts.genesFound, total: counts.genesTotal
      }),
      windows: file.windows.map(row),
      notSelected: file.candidates.map(row),
      found: t('{complete} found in full, {partial} in part, {absent} not at all.', {
        complete: complete.length, partial: partial.length, absent: absent.length
      }),
      partial: geneList(partial, names),
      absent: geneList(absent, names),
      empty: emptyText(file, counts)
    };
  });
};

/**
 * The footer's count: the ticked windows and the Files they are in.
 * @param {readonly WindowDraft[]} drafts
 */
export const proposalCountText = (drafts) => {
  const ticked = drafts.filter((draft) => draft.ticked);
  return t('{windows} in {files}', {
    windows: tn(ticked.length, '{count} window', '{count} windows'),
    files: tn(new Set(ticked.map((draft) => draft.fileKey)).size, '{count} file', '{count} files')
  });
};
