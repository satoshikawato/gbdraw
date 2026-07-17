import { createDefaultAnnotationStyle } from './state.js';

export const updateAnnotationStyle = (style, patch = {}) => {
  const normalized = createDefaultAnnotationStyle(style);
  Object.assign(normalized, patch && typeof patch === 'object' ? patch : {});
  return normalized;
};

export const createHatchStyle = ({ angle = 45, spacing = 5, color = '#666666', width = 1, cross = false } = {}) => ({
  angle: Number(angle),
  spacing: Math.max(0.1, Number(spacing) || 5),
  color: String(color || '#666666'),
  width: Math.max(0.1, Number(width) || 1),
  cross: Boolean(cross)
});

export const clearHatchStyle = (style) => updateAnnotationStyle(style, { hatch: null });
