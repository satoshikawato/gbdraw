const DEPTH_TRACK_FALLBACK_COLORS = [
  '#4A90E2',
  '#E45756',
  '#2CA02C',
  '#F28E2B',
  '#9467BD',
  '#8C564B',
  '#17BECF',
  '#7F7F7F'
];

const normalizePositiveNumber = (value) => {
  if (value === null || value === undefined || value === '') return null;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

const normalizeOptionalText = (value) => {
  const text = String(value ?? '').trim();
  return text.length > 0 ? text : null;
};

const cloneParams = (params = {}) => {
  if (!params || typeof params !== 'object' || Array.isArray(params)) return {};
  return { ...params };
};

const fallbackLabel = (index) => {
  const numericIndex = Math.max(0, Number(index) || 0);
  return numericIndex === 0 ? 'Depth' : `Depth ${numericIndex + 1}`;
};

const fallbackColor = (index, defaults = {}) => {
  if (typeof defaults.colorForIndex === 'function') {
    return String(defaults.colorForIndex(index) || '');
  }
  if (index === 0 && defaults.depthColor) return String(defaults.depthColor);
  return DEPTH_TRACK_FALLBACK_COLORS[index % DEPTH_TRACK_FALLBACK_COLORS.length];
};

export const depthFileBaseName = (file) => {
  const rawName = String(file?.name || '').split(/[\\/]/).pop().trim();
  if (!rawName) return '';
  return rawName.replace(/\.[^.]+$/, '').trim() || rawName;
};

export const depthLabelMatchesFile = (label, file) => {
  const labelText = String(label || '').toLowerCase();
  const fileText = depthFileBaseName(file).toLowerCase();
  if (!labelText || !fileText) return false;
  if (labelText.includes(fileText) || fileText.includes(labelText)) return true;
  const accessionTokens = fileText.match(/[a-z]{2,5}\d{5,}/g) || [];
  return accessionTokens.some((token) => labelText.includes(token));
};

export const depthFileSlotsFromValue = (value) => {
  if (Array.isArray(value)) return value.slice();
  return value ? [value] : [];
};

export const compactDepthFileSlots = (slots) => {
  const next = depthFileSlotsFromValue(slots);
  while (next.length > 0 && !next[next.length - 1]) {
    next.pop();
  }
  return next;
};

export const uploadedDepthFileCount = (value) => (
  depthFileSlotsFromValue(value).filter(Boolean).length
);

export const normalizeDepthTrackConfig = (entry, index, defaults = {}) => {
  const source = entry && typeof entry === 'object' && !Array.isArray(entry) ? entry : {};
  const hasHeight = Object.prototype.hasOwnProperty.call(source, 'height');
  const labelForIndex = typeof defaults.labelForIndex === 'function'
    ? defaults.labelForIndex
    : fallbackLabel;
  return {
    label: String(source.label ?? labelForIndex(index)),
    color: String(source.color || fallbackColor(index, defaults)),
    height: normalizePositiveNumber(hasHeight ? source.height : defaults.depthHeight),
    large_tick_interval: normalizePositiveNumber(
      source.large_tick_interval ?? source.tick_interval ?? defaults.largeTickInterval
    ),
    small_tick_interval: normalizePositiveNumber(
      source.small_tick_interval ?? defaults.smallTickInterval
    ),
    tick_font_size: normalizePositiveNumber(
      source.tick_font_size ?? defaults.tickFontSize
    )
  };
};

export const ensureDepthTrackConfigCount = (tracks, count, defaults = {}) => {
  const rawTracks = Array.isArray(tracks) ? tracks : [];
  const targetCount = Math.max(1, Number(count) || 1);
  const normalized = rawTracks.map((entry, index) => (
    normalizeDepthTrackConfig(entry, index, defaults)
  ));
  while (normalized.length < targetCount) {
    normalized.push(normalizeDepthTrackConfig(null, normalized.length, defaults));
  }
  return normalized;
};

export const reconcileDepthTracksToFiles = ({
  files,
  depthTracks,
  targetCount,
  defaults = {}
} = {}) => {
  const fileSlots = Array.isArray(files) ? files : depthFileSlotsFromValue(files);
  const rawTracks = Array.isArray(depthTracks) ? depthTracks : [];
  const usedTrackIndexes = new Set();
  return Array.from({ length: Math.max(1, Number(targetCount) || 1) }, (_, index) => {
    const file = fileSlots[index] || null;
    let sourceIndex = index;
    if (
      file &&
      rawTracks[index] &&
      !depthLabelMatchesFile(rawTracks[index]?.label, file)
    ) {
      const matchingIndex = rawTracks.findIndex((entry, candidateIndex) => (
        candidateIndex > index &&
        !usedTrackIndexes.has(candidateIndex) &&
        depthLabelMatchesFile(entry?.label, file)
      ));
      if (matchingIndex >= 0) sourceIndex = matchingIndex;
    }
    usedTrackIndexes.add(sourceIndex);
    return normalizeDepthTrackConfig(rawTracks[sourceIndex], index, defaults);
  });
};

export const removeDepthTrackAt = ({
  files,
  depthTracks,
  index,
  minCount = 1,
  defaults = {}
} = {}) => {
  const idx = Number(index);
  const originalFiles = depthFileSlotsFromValue(files);
  const originalTracks = Array.isArray(depthTracks) ? depthTracks.slice() : [];
  if (!Number.isInteger(idx) || idx < 0) {
    return {
      files: compactDepthFileSlots(originalFiles),
      depthTracks: ensureDepthTrackConfigCount(originalTracks, minCount, defaults),
      removed: false,
      removedFile: null,
      previousFiles: originalFiles
    };
  }
  const visibleCount = Math.max(minCount, originalFiles.length, 1);
  const nextFiles = originalFiles.slice();
  const nextTracks = originalTracks.slice();
  const removedFile = nextFiles[idx] || null;
  if (visibleCount <= 1) {
    nextFiles[0] = null;
    nextTracks.splice(0, nextTracks.length);
  } else {
    nextFiles.splice(idx, 1);
    if (idx < nextTracks.length) nextTracks.splice(idx, 1);
  }
  const compactedFiles = compactDepthFileSlots(nextFiles);
  const normalizedTracks = ensureDepthTrackConfigCount(
    nextTracks,
    Math.max(minCount, compactedFiles.length),
    defaults
  );
  return {
    files: compactedFiles,
    depthTracks: normalizedTracks,
    removed: true,
    removedFile,
    previousFiles: originalFiles
  };
};

export const depthSlotTrackIndex = (slot, fallbackIndex = null) => {
  const rawTrackIndex = slot?.params?.track_index;
  if (rawTrackIndex !== null && rawTrackIndex !== undefined && rawTrackIndex !== '') {
    const parsed = Number(rawTrackIndex);
    if (Number.isInteger(parsed) && parsed >= 0) return parsed;
  }
  const idMatch = String(slot?.id || '').trim().match(/^depth_(\d+)$/);
  if (idMatch) return Math.max(0, Number(idMatch[1]) - 1);
  if (fallbackIndex !== null && fallbackIndex !== undefined) {
    const fallback = Number(fallbackIndex);
    if (Number.isInteger(fallback) && fallback >= 0) return fallback;
  }
  return null;
};

const slotHasGeometry = (slot) => (
  normalizeOptionalText(slot?.width) !== null ||
  normalizeOptionalText(slot?.radius) !== null ||
  normalizeOptionalText(slot?.spacing) !== null ||
  normalizeOptionalText(slot?.inner_gap_px) !== null ||
  normalizeOptionalText(slot?.outer_gap_px) !== null ||
  Number(slot?.z || 0) !== 0
);

const paramsOnlyAllowedKeys = (params, allowedKeys) => {
  const allowed = new Set(allowedKeys);
  return Object.entries(cloneParams(params)).every(([key, value]) => (
    allowed.has(String(key)) || normalizeOptionalText(value) === null
  ));
};

export const isDefaultManagedDepthSlot = (slot) => {
  if (!slot || typeof slot !== 'object' || Array.isArray(slot)) return false;
  if (String(slot.renderer || '') !== 'depth') return false;
  if (slotHasGeometry(slot)) return false;
  const id = String(slot.id || '').trim();
  if (!/^depth(?:_\d+)?$/.test(id)) return false;
  return paramsOnlyAllowedKeys(slot.params, ['track_index', 'legend_label', 'managed']);
};

const cloneSlot = (slot) => ({
  ...slot,
  params: cloneParams(slot?.params)
});

const disableInvalidManualDepthSlot = (slot) => {
  const next = cloneSlot(slot);
  next.enabled = false;
  if (/^depth(?:_\d+)?$/.test(String(next.id || '').trim())) {
    next.id = `disabled_${next.id}`;
  }
  delete next.params.track_index;
  delete next.params.legend_label;
  return next;
};

export const reindexDepthSlots = ({
  slots,
  removedIndex,
  activeCount,
  managedPredicate = isDefaultManagedDepthSlot
} = {}) => {
  const idx = Number(removedIndex);
  const count = Math.max(0, Number(activeCount) || 0);
  if (!Array.isArray(slots) || !Number.isInteger(idx) || idx < 0) {
    return Array.isArray(slots) ? slots.map((slot) => cloneSlot(slot)) : [];
  }
  const nextSlots = [];
  slots.forEach((slot, ordinal) => {
    if (!slot || String(slot.renderer || '') !== 'depth') {
      nextSlots.push(cloneSlot(slot));
      return;
    }
    const oldTrackIndex = depthSlotTrackIndex(slot, ordinal);
    const isManaged = managedPredicate(slot);
    if (oldTrackIndex === null || oldTrackIndex === idx) {
      if (!isManaged) nextSlots.push(disableInvalidManualDepthSlot(slot));
      return;
    }
    const newTrackIndex = oldTrackIndex > idx ? oldTrackIndex - 1 : oldTrackIndex;
    if (newTrackIndex < 0 || newTrackIndex >= count) {
      if (!isManaged) nextSlots.push(disableInvalidManualDepthSlot(slot));
      return;
    }
    const next = cloneSlot(slot);
    next.params.track_index = newTrackIndex;
    nextSlots.push(next);
  });
  return nextSlots;
};

export const dropInvalidManagedDepthSlots = ({
  slots,
  activeCount,
  managedPredicate = isDefaultManagedDepthSlot
} = {}) => {
  const count = Math.max(0, Number(activeCount) || 0);
  if (!Array.isArray(slots)) return [];
  return slots
    .map((slot, ordinal) => {
      if (!slot || String(slot.renderer || '') !== 'depth') return cloneSlot(slot);
      const trackIndex = depthSlotTrackIndex(slot, ordinal);
      if (trackIndex !== null && trackIndex >= 0 && trackIndex < count) {
        const next = cloneSlot(slot);
        next.params.track_index = trackIndex;
        return next;
      }
      return managedPredicate(slot) ? null : disableInvalidManualDepthSlot(slot);
    })
    .filter(Boolean);
};

export const syncDepthSlotLabels = ({ slots, depthTracks, activeCount = null } = {}) => {
  if (!Array.isArray(slots)) return;
  const count = activeCount === null || activeCount === undefined
    ? (Array.isArray(depthTracks) ? depthTracks.length : 0)
    : Math.max(0, Number(activeCount) || 0);
  slots.forEach((slot, ordinal) => {
    if (!slot || String(slot.renderer || '') !== 'depth') return;
    const trackIndex = depthSlotTrackIndex(slot, ordinal);
    if (trackIndex === null || trackIndex < 0 || trackIndex >= count) {
      return;
    }
    const label = String(depthTracks?.[trackIndex]?.label ?? '').trim();
    slot.params = cloneParams(slot.params);
    slot.params.track_index = trackIndex;
    if (label) {
      slot.params.legend_label = label;
    } else {
      delete slot.params.legend_label;
    }
  });
};
