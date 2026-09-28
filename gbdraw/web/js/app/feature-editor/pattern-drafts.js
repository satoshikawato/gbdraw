// Transient display state for existing Color-rule patterns. Canonical rows stay
// with rule-actions; signatures only recognize unchanged rows after bulk restore.
const specificRuleRevision = (rule) => JSON.stringify([
  rule.feat, rule.qual, rule.val, rule.color, rule.cap, Boolean(rule.fromFile)
]);

export const createSpecificRulePatternDrafts = ({ rules, ref, invalidate }) => {
  const drafts = ref(new Map());
  const fieldIds = new WeakMap();
  let nextFieldId = 0;
  let revision = 0;
  let documentRevision = 0;
  const fieldId = (row) => {
    if (!fieldIds.has(row)) fieldIds.set(row, `color-rule-pattern-${++nextFieldId}`);
    return fieldIds.get(row);
  };
  const get = (row) => drafts.value.get(row) || null;
  const text = (row) => get(row)?.text ?? row.val;
  const edit = (row, value) => {
    if (!rules.includes(row)) return null;
    const valueText = String(value ?? '');
    if (get(row)?.text === valueText) return get(row);
    invalidate();
    if (valueText === row.val) {
      drafts.value.delete(row);
      return null;
    }
    drafts.value.set(row, { text: valueText, revision: ++revision, pending: false, error: null });
    return get(row);
  };
  const begin = (row, value) => {
    edit(row, value);
    if (!get(row)) drafts.value.set(row, { text: String(value ?? ''), error: null });
    const draft = get(row);
    draft.revision = ++revision;
    draft.pending = true;
    draft.error = null;
    return draft.revision;
  };
  const isCurrent = (row, token) => rules.includes(row) && get(row)?.revision === token;
  const suspend = () => {
    invalidate();
    for (const draft of drafts.value.values()) {
      if (!draft.pending) continue;
      draft.revision = ++revision;
      draft.pending = false;
    }
  };
  const revert = (row) => {
    invalidate();
    drafts.value.delete(row);
  };
  const clear = () => {
    suspend();
    drafts.value.clear();
    documentRevision += 1;
  };
  const reconcile = () => {
    for (const row of drafts.value.keys()) if (!rules.includes(row)) drafts.value.delete(row);
  };
  const capture = () => {
    suspend();
    return { documentRevision, entries: [...drafts.value].map(([row, draft]) => ({
      index: rules.indexOf(row), accepted: specificRuleRevision(row), draft: { ...draft }, id: fieldId(row)
    })) };
  };
  const restore = (snapshot) => {
    if (snapshot.documentRevision !== documentRevision) return;
    suspend();
    drafts.value.clear();
    for (const entry of snapshot.entries) {
      const row = rules[entry.index];
      if (!row || specificRuleRevision(row) !== entry.accepted) continue;
      fieldIds.set(row, entry.id);
      drafts.value.set(row, { ...entry.draft, revision: ++revision, pending: false });
    }
  };
  return { get, text, fieldId, edit, begin, isCurrent, suspend, revert, clear, reconcile, capture, restore };
};
