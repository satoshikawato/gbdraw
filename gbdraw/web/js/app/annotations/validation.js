import {
  annotationRecordBinding,
  annotationRecordSelectorFromTarget,
  resolveAnnotationRecord
} from './record-selector.js';

const cleanText = (value) => String(value ?? '').trim();

export const validateAnnotationCoordinates = ({ start, end }) => (
  [start, end].every((value) => cleanText(value) !== '' && Number.isSafeInteger(Number(value)) && Number(value) >= 1)
    ? '' : 'start and end must be positive integers (1-based coordinates).'
);

const annotationLabel = (set, annotation) => (
  `${cleanText(set?.id) || 'annotations'}/${cleanText(annotation?.id) || 'region'}`
);

const allAnnotations = (sets) => (
  (Array.isArray(sets) ? sets : []).flatMap((set) => (
    (Array.isArray(set?.annotations) ? set.annotations : []).map((annotation) => ({ set, annotation }))
  ))
);

const selectorError = (selector, records) => {
  if (selector.kind === 'recordIndex') {
    return selector.index < records.length ? '' : `record index #${selector.index + 1} is out of range`;
  }
  const matches = records.filter((record) => record.recordId === selector.value);
  if (matches.length === 0) return `record ID ${JSON.stringify(selector.value)} is not available`;
  if (matches.length > 1) return `record ID ${JSON.stringify(selector.value)} is duplicated; choose a #index entry`;
  return '';
};

export const validateAnnotationRecordTargets = (sets, catalog) => {
  const annotations = allAnnotations(sets);
  if (annotations.length === 0) return '';
  const catalogIssue = Array.isArray(catalog?.issues) ? catalog.issues[0] : '';
  if (catalogIssue) return `Region annotations: ${catalogIssue}`;

  const records = Array.isArray(catalog?.records) ? catalog.records : [];
  for (const { set, annotation } of annotations) {
    const label = annotationLabel(set, annotation);
    if (annotation.target?.kind === 'coordinateSpan') {
      const coordinateError = validateAnnotationCoordinates(annotation.target);
      if (coordinateError) return `Region annotation ${label}: ${coordinateError}`;
    }
    const parsed = annotationRecordSelectorFromTarget(annotation);
    if (parsed.error) return `Choose a valid target record for region annotation ${label}.`;
    if (!parsed.selector) {
      if (catalog?.requiresSelection) return `Choose a target record for region annotation ${label}.`;
      continue;
    }
    if (catalog?.allowExplicitSelectors === false) {
      return `Region annotation ${label}: clear the target record or enable Multi-record canvas.`;
    }
    if (annotationRecordBinding(annotation) && !resolveAnnotationRecord(catalog, annotation)) {
      return `Region annotation ${label}: the selected target record is no longer available; choose it again.`;
    }
    const error = selectorError(parsed.selector, records);
    if (error) return `Region annotation ${label}: ${error}.`;
  }
  return '';
};
