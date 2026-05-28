export const getDepthTrackFallbackLabel = (index) => {
  const numericIndex = Math.max(0, Number(index) || 0);
  return numericIndex === 0 ? 'Depth' : `Depth ${numericIndex + 1}`;
};

export const getDepthTrackFileBaseName = (file) => {
  const rawName = String(file?.name || '').split(/[\\/]/).pop().trim();
  if (!rawName) return '';
  const withoutExtension = rawName.replace(/\.[^.]+$/, '').trim();
  return withoutExtension || rawName;
};

export const getDepthTrackLabelFromFile = (file, index) => {
  const baseName = getDepthTrackFileBaseName(file);
  return baseName || getDepthTrackFallbackLabel(index);
};

export const isDepthTrackAutoLabel = (label, index, file = null) => {
  const normalized = String(label ?? '').trim();
  if (!normalized) return true;
  if (normalized === getDepthTrackFallbackLabel(index)) return true;
  const fileBaseName = getDepthTrackFileBaseName(file);
  return Boolean(fileBaseName && normalized === fileBaseName);
};
