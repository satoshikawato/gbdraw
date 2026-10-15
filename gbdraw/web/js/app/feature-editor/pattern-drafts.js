// @ts-check
// Transient display state for existing Color-rule patterns. Canonical rows stay
// with rule-actions; signatures only recognize unchanged rows after bulk restore.
/** @import { SpecificColorRule } from '../../services/specific-color-rules.js' */

/** @param {SpecificColorRule} rule */
const specificRuleRevision = (rule) => JSON.stringify([
  rule.feat, rule.qual, rule.val, rule.color, rule.cap, Boolean(rule.fromFile)
]);

/**
 * @typedef {object} SpecificRulePatternDraftsOptions
 * @property {() => SpecificColorRule[]} rules The canonical Color-rule rows of the shown mode's drawing
 *   (their owner is rule-actions).
 * @property {(value: Map<SpecificColorRule, PatternDraft>) => { value: Map<SpecificColorRule, PatternDraft> }} ref Vue `ref`
 * @property {() => void} changed Tells the rule owner a draft changed, after the change, so the
 *   owner releases the pattern commits it cancelled.
 */

/**
 * The text a row's pattern field shows while it differs from the accepted row.
 * `revision` and `pending` are set once a commit of the text begins.
 * @typedef {object} PatternDraft
 * @property {string} text
 * @property {number} [revision]
 * @property {boolean} [pending]
 * @property {unknown} error
 */

/**
 * @typedef {object} PatternDraftsSnapshot
 * @property {number} documentRevision
 * @property {{ index: number, accepted: string, draft: PatternDraft, id: string }[]} entries
 */

/** @param {SpecificRulePatternDraftsOptions} options */
export const createSpecificRulePatternDrafts = ({ rules, ref, changed }) => {
  const drafts = ref(/** @type {Map<SpecificColorRule, PatternDraft>} */ (new Map()));
  /** @type {WeakMap<SpecificColorRule, string>} */
  const fieldIds = new WeakMap();
  let nextFieldId = 0;
  let revision = 0;
  let documentRevision = 0;
  /** @param {SpecificColorRule} row */
  const fieldId = (row) => {
    if (!fieldIds.has(row)) fieldIds.set(row, `color-rule-pattern-${++nextFieldId}`);
    // Set just above when it was missing.
    return /** @type {string} */ (fieldIds.get(row));
  };
  /** @param {SpecificColorRule} row */
  const get = (row) => drafts.value.get(row) || null;
  /** @param {SpecificColorRule} row */
  const text = (row) => get(row)?.text ?? row.val;
  /** @param {SpecificColorRule} row @param {unknown} value */
  const edit = (row, value) => {
    if (!rules().includes(row)) return null;
    const valueText = String(value ?? '');
    if (get(row)?.text === valueText) return get(row);
    if (valueText === row.val) drafts.value.delete(row);
    else drafts.value.set(row, { text: valueText, revision: ++revision, pending: false, error: null });
    changed();
    return get(row);
  };
  /** @param {SpecificColorRule} row @param {unknown} value */
  const begin = (row, value) => {
    edit(row, value);
    if (!get(row)) drafts.value.set(row, { text: String(value ?? ''), error: null });
    // `edit` or the `set` just above stored a draft for the row.
    const draft = /** @type {PatternDraft} */ (get(row));
    draft.revision = ++revision;
    draft.pending = true;
    draft.error = null;
    changed();
    return draft.revision;
  };
  /** @param {SpecificColorRule} row @param {number} token */
  const isCurrent = (row, token) => rules().includes(row) && get(row)?.revision === token;
  const suspend = () => {
    for (const draft of drafts.value.values()) {
      if (!draft.pending) continue;
      draft.revision = ++revision;
      draft.pending = false;
    }
    changed();
  };
  /** @param {SpecificColorRule} row */
  const revert = (row) => {
    drafts.value.delete(row);
    changed();
  };
  const clear = () => {
    suspend();
    drafts.value.clear();
    documentRevision += 1;
  };
  const reconcile = () => {
    for (const row of drafts.value.keys()) if (!rules().includes(row)) drafts.value.delete(row);
  };
  const capture = () => {
    suspend();
    return { documentRevision, entries: [...drafts.value].map(([row, draft]) => ({
      index: rules().indexOf(row), accepted: specificRuleRevision(row), draft: { ...draft }, id: fieldId(row)
    })) };
  };
  /** @param {PatternDraftsSnapshot} snapshot */
  const restore = (snapshot) => {
    if (snapshot.documentRevision !== documentRevision) return;
    suspend();
    drafts.value.clear();
    for (const entry of snapshot.entries) {
      const row = rules()[entry.index];
      if (!row || specificRuleRevision(row) !== entry.accepted) continue;
      fieldIds.set(row, entry.id);
      drafts.value.set(row, { ...entry.draft, revision: ++revision, pending: false });
    }
  };
  return { get, text, fieldId, edit, begin, isCurrent, suspend, revert, clear, reconcile, capture, restore };
};
