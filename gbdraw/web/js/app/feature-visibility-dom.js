import { recordStructuralMetric } from '../services/runtime-test-hooks.js';

export const FEATURE_VISIBILITY_PREVIEW_ATTRIBUTE =
  'data-gbdraw-feature-visibility-preview';
export const FEATURE_RENDERER_DISPLAY_ATTRIBUTE =
  'data-gbdraw-feature-renderer-display';

const normalizeVisibilityMode = (value) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  if (normalized === 'suppress') return 'exclude_matching';
  return ['on', 'off', 'exclude_matching'].includes(normalized)
    ? normalized
    : 'default';
};

const setAttribute = (mutation, element, name, value) => (
  mutation
    ? mutation.setAttribute(element, name, value)
    : (element.setAttribute(name, value), true)
);

const removeAttribute = (mutation, element, name) => (
  mutation
    ? mutation.removeAttribute(element, name)
    : (element.removeAttribute(name), true)
);

/**
 * Apply an editor-owned Feature visibility decision without losing renderer
 * visibility. The first explicit decision records the renderer's display
 * attribute; returning to default restores that exact baseline.
 */
export const applyFeatureVisibility = (element, modeRaw, {
  markPreview = false,
  mutation = null,
  reason = 'feature-visibility'
} = {}) => {
  if (!element) return false;
  const mode = normalizeVisibilityMode(modeRaw);
  // Matching exclusion is orthogonal to visibility. For DOM projection it is
  // the same as no visibility decision, while its canonical value remains
  // available to matching/selection logic.
  const visibilityMode = mode === 'exclude_matching' ? 'default' : mode;
  const previewMode = element.getAttribute?.(FEATURE_VISIBILITY_PREVIEW_ATTRIBUTE);
  let changed = false;

  if (visibilityMode !== 'default') {
    if (markPreview && previewMode === null) {
      const rendererDisplay = element.getAttribute?.('display');
      setAttribute(
        mutation,
        element,
        FEATURE_RENDERER_DISPLAY_ATTRIBUTE,
        rendererDisplay ?? ''
      );
      recordStructuralMetric('featureRendererBaselineCaptureCount', 1, {
        reason: String(reason || 'feature-visibility'),
        elementOwner: element
      });
      changed = true;
    }
    if (markPreview && previewMode !== visibilityMode) {
      setAttribute(mutation, element, FEATURE_VISIBILITY_PREVIEW_ATTRIBUTE, visibilityMode);
      changed = true;
    }
    if (visibilityMode === 'off') {
      if (element.getAttribute?.('display') !== 'none') {
        setAttribute(mutation, element, 'display', 'none');
        changed = true;
      }
    } else if (element.getAttribute?.('display') !== null) {
      removeAttribute(mutation, element, 'display');
      changed = true;
    }
    return changed;
  }

  if (previewMode === null) return false;
  const rendererDisplay = element.getAttribute?.(FEATURE_RENDERER_DISPLAY_ATTRIBUTE);
  if (rendererDisplay) setAttribute(mutation, element, 'display', rendererDisplay);
  else removeAttribute(mutation, element, 'display');
  removeAttribute(mutation, element, FEATURE_VISIBILITY_PREVIEW_ATTRIBUTE);
  removeAttribute(mutation, element, FEATURE_RENDERER_DISPLAY_ATTRIBUTE);
  return true;
};
