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
// A failed target is a producer diagnostic; the annotation ID stays private.
const targetIssue = (reason) => ({ code: 'ANNOTATION_TARGET', context: { reason } });

const allAnnotations = (sets) => (
  (Array.isArray(sets) ? sets : []).flatMap((set) => (
    (Array.isArray(set?.annotations) ? set.annotations : []).map((annotation) => ({ set, annotation }))
  ))
);

const selectorReason = (selector, records) => {
  if (selector.kind === 'recordIndex') return selector.index < records.length ? '' : 'OUT_OF_RANGE';
  const matches = records.filter((record) => record.recordId === selector.value);
  if (matches.length === 0) return 'NO_MATCH';
  return matches.length > 1 ? 'AMBIGUOUS' : '';
};

/** Returns null or the first failure as a { code, context } diagnostic. */
export const validateAnnotationRecordTargets = (sets, catalog) => {
  const annotations = allAnnotations(sets);
  if (annotations.length === 0) return null;
  const catalogIssue = Array.isArray(catalog?.issues) ? catalog.issues[0] : null;
  if (catalogIssue) return catalogIssue;

  const records = Array.isArray(catalog?.records) ? catalog.records : [];
  for (const { annotation } of annotations) {
    if (annotation.target?.kind === 'coordinateSpan' && validateAnnotationCoordinates(annotation.target)) {
      return targetIssue('POSITIVE_INTEGER');
    }
    const parsed = annotationRecordSelectorFromTarget(annotation);
    if (parsed.error) return targetIssue('TARGET_RECORD');
    if (!parsed.selector) {
      if (catalog?.requiresSelection) return targetIssue('TARGET_RECORD');
      continue;
    }
    if (annotationRecordBinding(annotation) && !resolveAnnotationRecord(catalog, annotation)) {
      return targetIssue('TARGET_RECORD');
    }
    const reason = selectorReason(parsed.selector, records);
    if (reason) return targetIssue(reason);
  }
  return null;
};
