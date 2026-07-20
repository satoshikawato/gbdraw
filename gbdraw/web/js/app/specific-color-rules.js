import { resolveColorToHex } from './color-utils.js';
import { parseSpecificRules } from './file-imports.js';

export const SPECIFIC_COLOR_FILE_OWNER = 'specific-color-file';

const normalizeText = (value) => String(value ?? '').trim();
const normalizeColor = (value) => String(resolveColorToHex(normalizeText(value)) || '').toLowerCase();

export const normalizeSpecificRule = (rule, { fromFile = Boolean(rule?.fromFile) } = {}) => ({
  feat: normalizeText(rule?.feat),
  qual: normalizeText(rule?.qual),
  val: normalizeText(rule?.val),
  color: normalizeColor(rule?.color),
  cap: normalizeText(rule?.cap),
  ...(fromFile ? { fromFile: true } : {})
});

export const specificRuleIdentity = (rule) => {
  const normalized = normalizeSpecificRule(rule);
  return [
    normalized.feat,
    normalized.qual,
    normalized.val,
    normalized.color,
    normalized.cap
  ];
};

const identityKey = (rule) => JSON.stringify(specificRuleIdentity(rule));

export const applySpecificRuleProvenance = (canonicalRules, storedRules) => {
  const storedFileIdentities = new Set(
    (Array.isArray(storedRules) ? storedRules : [])
      .filter((rule) => rule?.fromFile)
      .map(identityKey)
  );
  return (Array.isArray(canonicalRules) ? canonicalRules : []).map((rule) => {
    const normalized = normalizeSpecificRule(rule, { fromFile: false });
    return storedFileIdentities.has(identityKey(normalized))
      ? { ...normalized, fromFile: true }
      : normalized;
  });
};

export const buildLegendIntents = (rules, { conflictPolicy = 'reject' } = {}) => {
  const intents = [];
  const byCaption = new Map();
  const conflicts = [];

  (Array.isArray(rules) ? rules : []).forEach((rule, index) => {
    const normalized = normalizeSpecificRule(rule);
    if (!normalized.cap) return;
    const previous = byCaption.get(normalized.cap);
    if (previous && previous.color !== normalized.color) {
      const conflict = {
        caption: normalized.cap,
        previousColor: previous.color,
        nextColor: normalized.color,
        ruleIndex: index
      };
      conflicts.push(conflict);
      if (conflictPolicy === 'reject') {
        throw new Error(
          `Specific-color caption "${normalized.cap}" uses multiple colors (${previous.color}, ${normalized.color}).`
        );
      }
      previous.color = normalized.color;
      return;
    }
    if (previous) return;
    const intent = { caption: normalized.cap, color: normalized.color };
    byCaption.set(intent.caption, intent);
    intents.push(intent);
  });

  return { intents, conflicts };
};

export const diffLegendIntents = (currentEntries, desiredIntents) => {
  const desired = new Map(
    (Array.isArray(desiredIntents) ? desiredIntents : [])
      .map((entry) => [normalizeText(entry?.caption), normalizeColor(entry?.color)])
      .filter(([caption]) => caption)
  );
  const current = new Map();
  (Array.isArray(currentEntries) ? currentEntries : []).forEach((entry) => {
    const caption = normalizeText(entry?.caption);
    if (!caption || current.has(caption)) return;
    current.set(caption, { ...entry, caption, color: normalizeColor(entry?.color) });
  });

  const diff = { add: [], update: [], remove: [], unchanged: [] };
  current.forEach((entry, caption) => {
    if (!desired.has(caption)) {
      diff.remove.push(entry);
    } else if (desired.get(caption) === entry.color) {
      diff.unchanged.push(entry);
    } else {
      diff.update.push({ ...entry, color: desired.get(caption) });
    }
    desired.delete(caption);
  });
  desired.forEach((color, caption) => diff.add.push({ caption, color }));
  return diff;
};

export const prepareSpecificColorImport = (text, currentRules = []) => {
  const { rules } = parseSpecificRules(text);
  const { intents } = buildLegendIntents(rules, { conflictPolicy: 'reject' });
  const retainedRules = (Array.isArray(currentRules) ? currentRules : [])
    .filter((rule) => !rule?.fromFile)
    .map((rule) => normalizeSpecificRule(rule));
  const nextRules = [...retainedRules, ...rules.map((rule) => normalizeSpecificRule(rule, { fromFile: true }))];
  return {
    nextRules,
    intents,
    fileLegendCaptions: intents.map((intent) => intent.caption),
    importedCount: rules.length
  };
};

