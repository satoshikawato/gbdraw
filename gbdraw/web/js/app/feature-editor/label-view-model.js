const text = (value) => String(value ?? '');
const trimmed = (value) => text(value).trim();

const firstText = (...values) => {
  for (const value of values) {
    if (Array.isArray(value)) {
      const nested = firstText(...value);
      if (nested) return nested;
      continue;
    }
    const normalized = trimmed(value);
    if (normalized) return normalized;
  }
  return '';
};

const qualifierValues = (feature, wanted) => {
  const qualifiers = feature?.qualifiers;
  if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) return [];
  const key = Object.keys(qualifiers).find((candidate) => (
    trimmed(candidate).toLowerCase() === wanted
  ));
  const value = key ? qualifiers[key] : [];
  return Array.isArray(value) ? value : [value];
};

/** Keep authored label whitespace; only non-string values are rejected. */
export const normalizeFeatureLabelText = (value) => (
  typeof value === 'string' ? value : null
);

export const resolveFeatureSourceAnnotationLabel = (feature = {}) => {
  const explicit = firstText(
    feature.sourceAnnotationLabel,
    feature.source_annotation_label
  );
  if (explicit) {
    return Object.freeze({
      text: explicit,
      qualifier: firstText(
        feature.sourceAnnotationQualifier,
        feature.source_annotation_qualifier
      ) || 'label'
    });
  }
  for (const qualifier of ['product', 'gene', 'locus_tag', 'note']) {
    const value = firstText(feature[qualifier], ...qualifierValues(feature, qualifier));
    if (value) return Object.freeze({ text: value, qualifier });
  }
  const type = firstText(feature.type);
  const start = Number(feature.start);
  const end = Number(feature.end);
  return Object.freeze({
    text: type && Number.isFinite(start) && Number.isFinite(end)
      ? `${type} at ${start}..${end}`
      : type,
    qualifier: type ? 'type' : ''
  });
};

const firstOwnValue = (map, keys) => {
  for (const key of keys) {
    if (key && Object.prototype.hasOwnProperty.call(map || {}, key)) {
      return { key, value: text(map[key]) };
    }
  }
  return null;
};

/**
 * Resolve the editable label independently of its SVG presentation.
 * `presentationText` is accepted only as a current-rendered fallback; the
 * persisted feature/source maps remain the semantic owners.
 */
export const resolveFeatureLabelViewModel = ({
  feature = {},
  renderedSvgIds = [],
  featureOverrides = {},
  bulkOverrides = {},
  featureOverrideSources = {},
  presentationText = null,
  presentationSourceText = ''
} = {}) => {
  const svgIds = [...new Set([
    ...(Array.isArray(renderedSvgIds) ? renderedSvgIds : []),
    feature.renderedFeatureSvgId,
    feature.rendered_feature_svg_id,
    feature.svgId,
    feature.svg_id
  ].map(trimmed).filter(Boolean))];
  const annotation = resolveFeatureSourceAnnotationLabel(feature);
  const explicit = firstOwnValue(featureOverrides, svgIds);
  const explicitSource = firstOwnValue(featureOverrideSources, svgIds);
  const sourceText = trimmed(explicitSource?.value)
    || trimmed(presentationSourceText)
    || annotation.text;
  const grouped = sourceText && Object.prototype.hasOwnProperty.call(bulkOverrides || {}, sourceText)
    ? text(bulkOverrides[sourceText])
    : null;
  const authoredPresentation = normalizeFeatureLabelText(presentationText);
  const effectiveText = explicit
    ? explicit.value
    : (grouped !== null
        ? grouped
        : (authoredPresentation !== null ? authoredPresentation : sourceText));
  const origin = explicit
    ? 'feature-override'
    : (grouped !== null
        ? 'source-group'
        : (authoredPresentation !== null && authoredPresentation !== sourceText
            ? 'rendered-svg'
            : 'source-annotation'));
  const originLabel = {
    'feature-override': 'Feature-specific label',
    'source-group': `Inherited from source label: ${sourceText}`,
    'rendered-svg': 'Current rendered label',
    'source-annotation': annotation.qualifier
      ? `Inherited from source annotation: ${annotation.qualifier}`
      : 'Inherited from source annotation'
  }[origin];
  return Object.freeze({
    effectiveText,
    sourceText,
    sourceQualifier: annotation.qualifier,
    explicitValue: explicit ? explicit.value : null,
    origin,
    originLabel,
    canReset: Boolean(explicit),
    renderedSvgIds: Object.freeze(svgIds)
  });
};
