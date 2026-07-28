const isGroup = (element) => String(element?.tagName || '').toLowerCase() === 'g';
const isSvg = (element) => String(element?.tagName || '').toLowerCase() === 'svg';
const hasSemanticRecordId = (group) =>
  isGroup(group) && group.hasAttribute?.('data-gbdraw-record-id');
const hasLegacyRecordId = (group) =>
  isGroup(group) && String(group.id || '').startsWith('record_');
const DEFINITION_ROLES = new Set(['record-definition', 'record-definition-row']);
const isDefinitionGroup = (group) =>
  DEFINITION_ROLES.has(
    String(group?.getAttribute?.('data-gbdraw-role') || '').trim()
  ) ||
  /_definition(?:_record_\d+)?(?:_row)?$/.test(String(group?.id || ''));

const semanticRecordIdentity = (group) => {
  const recordIndex = String(
    group?.getAttribute?.('data-gbdraw-record-index') || ''
  ).trim();
  if (recordIndex) return `index:${recordIndex}`;

  const recordId = String(
    group?.getAttribute?.('data-gbdraw-record-id') || ''
  ).trim();
  if (recordId) return `id:${recordId}`;

  return `group:${String(group?.id || '')}`;
};

export const isRecordGroup = (group) =>
  hasSemanticRecordId(group) || hasLegacyRecordId(group);

export const directRecordGroups = (svg) => {
  const groups = Array.from(svg?.children || [])
    .filter(isGroup)
    .filter((group) => !isDefinitionGroup(group));
  const semanticGroups = groups.filter(hasSemanticRecordId);
  return semanticGroups.length > 0
    ? semanticGroups
    : groups.filter(hasLegacyRecordId);
};

export const isMultiRecordCanvasSvg = (svg) => {
  const groups = directRecordGroups(svg);
  const semanticGroups = groups.filter(hasSemanticRecordId);
  if (semanticGroups.length > 0) {
    return new Set(semanticGroups.map(semanticRecordIdentity)).size > 1;
  }
  return groups.length > 1;
};

export const closestRecordGroup = (element) => {
  const ancestors = [];
  for (let current = element; current; current = current.parentElement) {
    if (isRecordGroup(current) && !isDefinitionGroup(current)) {
      ancestors.push(current);
    }
  }

  // Current Circular multi-record SVGs expose a semantic top-level wrapper.
  // Prefer it so record dragging moves the axis, tracks, and definition together.
  const semanticTopLevelGroup = ancestors.find(
    (group) => hasSemanticRecordId(group) && isSvg(group.parentElement)
  );
  if (semanticTopLevelGroup) return semanticTopLevelGroup;

  // Persisted SVGs from before semantic wrappers used top-level record_N IDs.
  const legacyTopLevelGroup = ancestors.find(
    (group) => hasLegacyRecordId(group) && isSvg(group.parentElement)
  );
  if (legacyTopLevelGroup) return legacyTopLevelGroup;

  return (
    ancestors.find(hasSemanticRecordId) ||
    ancestors.find(hasLegacyRecordId) ||
    null
  );
};
