const RESOURCE_REFERENCE_FIELDS = new Set([
  'resourceId',
  'gffResourceId',
  'fastaResourceId'
]);

export const isCanonicalResourceReferenceField = (field) => (
  RESOURCE_REFERENCE_FIELDS.has(field)
);

export const collectCanonicalResourceIds = (value, target = new Set()) => {
  if (Array.isArray(value)) {
    value.forEach((item) => collectCanonicalResourceIds(item, target));
    return target;
  }
  if (!value || typeof value !== 'object') return target;
  Object.entries(value).forEach(([key, item]) => {
    if (
      isCanonicalResourceReferenceField(key)
      && typeof item === 'string'
      && item.trim()
    ) {
      target.add(item.trim());
      return;
    }
    collectCanonicalResourceIds(item, target);
  });
  return target;
};
