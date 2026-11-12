// @ts-check
// The pure middles of a Label TSV load and export (the label table owner,
// `app/feature-editor/label-table.js`): the label intent a load compares
// before and after it waits, the label edits a loaded table makes, the
// notice each outcome shows, and the rows Export Label TSV writes. State-free
// (R13).
import { escapeRegexLiteral } from './feature-selector.js';
import { featureIdentityKeyOf } from './feature-placement.js';
import { serializeLabelOverrideRows } from './label-override-table.js';
import { RULE_INPUT_CHANGE_WORDS } from '../utils/error-normalization.js';
import { joinSentences, t, tm, tn } from '../utils/i18n.js';
import { normalizeTsvCell } from '../utils/tsv-cell.js';
/** @import { ExtractedFeature } from './feature-metadata-extraction.js' */
/** @import { LabelOverrideTsvRow, StoredLabelOverrideRow } from './label-override-table.js' */

/**
 * A label a Label TSV row can select: the rendered feature ID it is bound to
 * ('' when none), that feature's row (null when none), and the label's source text.
 * @typedef {{ featureId: string, feature: ExtractedFeature | null, sourceText: string }} LabelImportTarget
 */
/**
 * A label edit a loaded table makes: the label, the row that won it, and
 * whether the label's feature has a feature key.
 * @typedef {{ entry: LabelImportTarget, matchedRow: LabelOverrideTsvRow, tracked: boolean }} LabelImportEdit
 */
/**
 * The label edits a loaded table makes, and the matched labels it skips
 * because the edit needs a feature key the label lacks.
 * @typedef {{ applicable: LabelImportEdit[], skippedNonTrackableCount: number }} LabelImportPlan
 */
/**
 * Why a load is declined before its rows are matched: what changed while the
 * table loaded (OV-447, OV-449), or that no diagram is shown.
 * @typedef {'label-intent' | 'diagram' | 'rule-inputs' | 'no-diagram'} LabelImportDecline
 */
/**
 * The outcome of a Label TSV load the notice reports: declined for a change
 * or a missing diagram (`declined`), declined for labels a `location` or
 * `record_location` row cannot decide (`undecided`), or the plan, which
 * applies when it edits a label.
 * @typedef {{ declined: LabelImportDecline } | { undecided: number } | { plan: LabelImportPlan }} LabelImportOutcome
 */
/** @typedef {{ text: string, items?: string[], tone?: 'warning' }} LabelImportNotice */
/**
 * The drawing members that hold the label intent.
 * @typedef {{ featureOverrides: unknown, labelTextBulkOverrides: unknown }} LabelIntentHolder
 */

/**
 * The label intent of a drawing, as a load compares it before and after it waits.
 * @param {LabelIntentHolder} drawing
 */
export const labelIntentSignature = (drawing) => JSON.stringify([drawing.featureOverrides, drawing.labelTextBulkOverrides]);

// The import replaces the label intent once; the displayed Result shows it
// now and every other Result when it is displayed (B6, R3).
/**
 * @param {readonly LabelImportTarget[]} labels
 * @param {readonly LabelOverrideTsvRow[]} rows
 * @param {readonly number[]} winners The index of the row that wins each label, by label.
 * @returns {LabelImportPlan}
 */
export const planLabelImport = (labels, rows, winners) => {
  let skippedNonTrackableCount = 0;
  const applicable = labels.flatMap((entry, index) => {
    const matchedRow = rows[winners[index]];
    if (!matchedRow) return [];
    const tracked = Boolean(featureIdentityKeyOf(entry.feature));
    // A global `label` row follows the source text; any other row its
    // feature, as that feature's label edit (design Q4).
    if (matchedRow.isGlobalLabelRule ? !entry.sourceText : !tracked) {
      skippedNonTrackableCount += 1;
      return [];
    }
    return [{ entry, matchedRow, tracked }];
  });
  return { applicable, skippedNonTrackableCount };
};

// OV-447, OV-449 (Owner decisions): the table was evaluated against a state
// that changed meanwhile, so the load is declined and names the change. An
// edit of the label intent stays; it is named first, since it also redraws
// the diagram and makes the rule preparation stale.
/** @param {LabelImportDecline} reason */
const declinedSentence = (reason) => ({
  'label-intent': () => t('Not applied: label text, label visibility, or feature visibility changed while the table loaded. The changes were kept. Load the table again.'),
  diagram: () => t('Not applied: the displayed diagram changed while the table loaded. Load the table again.'),
  'rule-inputs': () => t('Not applied: {inputs} changed while the table loaded. Load the table again.', { inputs: tm(RULE_INPUT_CHANGE_WORDS) }),
  'no-diagram': () => t('No diagram is shown, so nothing was applied.')
})[reason]();

/**
 * The notice of a Label TSV load. A table that applies to no label is
 * declined whole, so the label edits stay (OV-131).
 * @param {number} rowCount The rows the table holds.
 * @param {LabelImportOutcome} outcome
 * @returns {LabelImportNotice}
 */
export const labelImportNotice = (rowCount, outcome) => {
  const loadedRowsText = tn(rowCount, 'Loaded {count} row.', 'Loaded {count} rows.');
  if ('declined' in outcome) return { text: joinSentences(loadedRowsText, declinedSentence(outcome.declined)), tone: 'warning' };
  if ('undecided' in outcome) {
    return { text: joinSentences(loadedRowsText, t('Not applied.')), items: [
      tn(outcome.undecided,
        'location and record_location rows cannot be matched to {count} label on a cropped, reverse-complemented, or rotated record until the diagram records where it drew it. Generate the diagram again, then load the table.',
        'location and record_location rows cannot be matched to {count} labels on cropped, reverse-complemented, or rotated records until the diagram records where it drew them. Generate the diagram again, then load the table.')
    ] };
  }
  const { applicable, skippedNonTrackableCount } = outcome.plan;
  if (applicable.length === 0) {
    return { text: joinSentences(loadedRowsText, t('Not applied; the existing label edits were kept.')), items: [
      skippedNonTrackableCount > 0
        ? tn(skippedNonTrackableCount, 'The {count} matched label lacked a feature key.', 'The {count} matched labels lacked a feature key.')
        : t('No row matched a label of the diagram.')
    ] };
  }
  return {
    text: joinSentences(loadedRowsText, tn(applicable.length, 'Applied to {count} label.', 'Applied to {count} labels.')),
    items: skippedNonTrackableCount > 0
      ? [tn(skippedNonTrackableCount, '{count} match lacked a feature key and was not applied.', '{count} matches lacked a feature key and were not applied.')]
      : []
  };
};

// Export Label TSV writes the label rules the next Generate sends: a saved
// label table, or the bulk label edits as `* * label` rows. Per-feature label
// edits are identity rows; Export Feature Edits TSV writes them (design Q4 6.4).
/**
 * @param {readonly StoredLabelOverrideRow[] | null | undefined} savedRows The saved label table.
 * @param {Readonly<Record<string, string>>} bulkOverrides The bulk label edits, by source text.
 * @returns {string[]} The table's lines, without line ends.
 */
export const labelTableExportLines = (savedRows, bulkOverrides) => {
  const savedTable = serializeLabelOverrideRows(savedRows);
  return savedTable
    ? savedTable.trimEnd().split('\n')
    : Object.keys(bulkOverrides).sort((a, b) => a.localeCompare(b)).filter(Boolean)
      .map((sourceText) => `*\t*\tlabel\t^${escapeRegexLiteral(sourceText)}$\t${normalizeTsvCell(bulkOverrides[sourceText])}`);
};
