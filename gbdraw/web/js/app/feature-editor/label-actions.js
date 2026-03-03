import { buildFeatureMetadataMap, buildLabelOverrideRows, parseLabelOverrideTsv } from './label-override-table.js';

const FEATURE_SELECTOR = 'path[id^="f"], polygon[id^="f"], rect[id^="f"]';
const EXCLUDED_GROUP_SELECTOR =
  '#legend, #feature_legend, #pairwise_legend, #horizontal_legend, #vertical_legend, #length_bar, #tick';
const EDITABLE_LABEL_SELECTOR = 'text[data-label-editable="true"]';

const toNumber = (value, fallback = 0) => {
  const parsed = Number.parseFloat(value);
  return Number.isFinite(parsed) ? parsed : fallback;
};

const normalizeKeyToken = (value) => String(value ?? '').trim().toLowerCase();
const normalizePositionToken = (value) => {
  const token = String(value ?? '').trim();
  if (!token || token === '*') return '*';
  const match = token.match(/^(\d+)\.\.(\d+):(.+)$/);
  if (!match) return token;
  const strandRaw = String(match[3] ?? '').trim().toLowerCase();
  const strand =
    strandRaw === '+' || strandRaw === 'positive' || strandRaw === 'forward' || strandRaw === '1'
      ? '+'
      : strandRaw === '-' || strandRaw === 'negative' || strandRaw === 'reverse' || strandRaw === '-1'
        ? '-'
        : 'undefined';
  return `${match[1]}..${match[2]}:${strand}`;
};
const wildcardOrExactMatch = (ruleValue, actualValue) =>
  String(ruleValue || '') === '*' || String(ruleValue || '') === String(actualValue || '');
const makeSafeFilename = (name, fallback = 'gbdraw') => {
  const cleaned = String(name || '')
    .replace(/[^\w.-]+/g, '_')
    .replace(/^_+|_+$/g, '');
  return cleaned || fallback;
};
const downloadTextFile = (filename, text) => {
  const blob = new Blob([text], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
};

const getSvgPoint = (svg, element, x, y) => {
  if (!svg || !element) return { x, y };
  const point = svg.createSVGPoint();
  point.x = x;
  point.y = y;
  const ctm = element.getCTM();
  if (!ctm) return { x, y };
  const transformed = point.matrixTransform(ctm);
  return { x: transformed.x, y: transformed.y };
};

const isFinitePoint = (point) =>
  Boolean(point) && Number.isFinite(point.x) && Number.isFinite(point.y);

const getTextPathHref = (textPathEl) =>
  String(textPathEl?.getAttribute('href') || textPathEl?.getAttribute('xlink:href') || '').trim();

const resolveEmbeddedLabelPathElement = (svg, textEl) => {
  const textPathEl = textEl?.querySelector?.('textPath');
  if (!textPathEl || !svg) return null;
  const href = getTextPathHref(textPathEl);
  if (href.startsWith('#')) {
    const linkedEl = svg.getElementById(href.slice(1));
    if (linkedEl?.tagName?.toLowerCase() === 'path') return linkedEl;
  }
  const prev = textEl.previousElementSibling;
  if (prev?.tagName?.toLowerCase() === 'path') return prev;
  return null;
};

const getEmbeddedLabelAnchor = (svg, textEl) => {
  const pathEl = resolveEmbeddedLabelPathElement(svg, textEl);
  if (!pathEl) return null;
  try {
    const totalLength = pathEl.getTotalLength();
    if (Number.isFinite(totalLength) && totalLength > 0) {
      const midpoint = pathEl.getPointAtLength(totalLength / 2);
      return getSvgPoint(svg, pathEl, midpoint.x, midpoint.y);
    }
  } catch {
    // Fall back to path bbox center when path length APIs are unavailable.
  }
  try {
    const bbox = pathEl.getBBox();
    return getSvgPoint(svg, pathEl, bbox.x + bbox.width / 2, bbox.y + bbox.height / 2);
  } catch {
    return null;
  }
};

const getElementCenter = (svg, element) => {
  try {
    const bbox = element.getBBox();
    return getSvgPoint(svg, element, bbox.x + bbox.width / 2, bbox.y + bbox.height / 2);
  } catch {
    if (element?.tagName?.toLowerCase() === 'text') {
      const anchor = getEmbeddedLabelAnchor(svg, element);
      if (isFinitePoint(anchor)) return anchor;
    }
    return { x: 0, y: 0 };
  }
};

const getLabelText = (textEl) => {
  const textPath = textEl.querySelector('textPath');
  if (textPath) return textPath.textContent || '';
  return textEl.textContent || '';
};

const setLabelText = (textEl, value) => {
  const nextText = String(value ?? '');
  const textPath = textEl.querySelector('textPath');
  if (textPath) {
    textPath.textContent = nextText;
    return;
  }
  textEl.textContent = nextText;
};

const hasExcludedAncestor = (textEl) => Boolean(textEl.closest(EXCLUDED_GROUP_SELECTOR));

const getCircularFeatureAnchor = (svg, textEl) => {
  const inLabelsGroup = Boolean(textEl.closest('g#labels'));
  if (!inLabelsGroup) return null;

  const prev = textEl.previousElementSibling;
  const prev2 = prev ? prev.previousElementSibling : null;
  const lines = [prev, prev2].filter(
    (candidate) => candidate && candidate.tagName && candidate.tagName.toLowerCase() === 'line'
  );
  if (lines.length === 0) return null;

  const featureLine = lines[0];
  const x2 = toNumber(featureLine.getAttribute('x2'), NaN);
  const y2 = toNumber(featureLine.getAttribute('y2'), NaN);
  if (!Number.isFinite(x2) || !Number.isFinite(y2)) return null;
  return getSvgPoint(svg, featureLine, x2, y2);
};

const getLabelReferencePoint = (svg, textEl, mode) => {
  if (mode === 'circular') {
    const circularAnchor = getCircularFeatureAnchor(svg, textEl);
    if (isFinitePoint(circularAnchor)) return circularAnchor;
  }
  const embeddedAnchor = getEmbeddedLabelAnchor(svg, textEl);
  if (isFinitePoint(embeddedAnchor)) return embeddedAnchor;
  const center = getElementCenter(svg, textEl);
  return isFinitePoint(center) ? center : null;
};

const collectEditableLabelElements = (svg, mode) => {
  const labels = new Set();

  if (mode === 'circular') {
    svg.querySelectorAll('g#labels text').forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
    svg.querySelectorAll('text > textPath').forEach((textPathEl) => {
      const parentText = textPathEl.parentElement;
      if (!parentText) return;
      if (hasExcludedAncestor(parentText)) return;
      labels.add(parentText);
    });
  } else {
    svg.querySelectorAll('text[dominant-baseline="central"]').forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
  }

  return Array.from(labels).sort((a, b) => {
    const aCenter = getLabelReferencePoint(svg, a, mode) || { x: 0, y: 0 };
    const bCenter = getLabelReferencePoint(svg, b, mode) || { x: 0, y: 0 };
    if (Math.abs(aCenter.y - bCenter.y) > 1) return aCenter.y - bCenter.y;
    return aCenter.x - bCenter.x;
  });
};

const collectFeatureGeometry = (svg) => {
  const grouped = new Map();
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((el) => {
    const id = el.getAttribute('id');
    if (!id) return;
    const center = getElementCenter(svg, el);
    const groupId = el.closest('g[id]')?.id || '';
    if (!grouped.has(id)) {
      grouped.set(id, { id, x: 0, y: 0, n: 0, groupId });
    }
    const item = grouped.get(id);
    item.x += center.x;
    item.y += center.y;
    item.n += 1;
  });

  const all = [];
  const byGroup = new Map();
  grouped.forEach((item) => {
    if (item.n <= 0) return;
    const centroid = {
      id: item.id,
      x: item.x / item.n,
      y: item.y / item.n,
      groupId: item.groupId
    };
    all.push(centroid);
    if (!byGroup.has(centroid.groupId)) {
      byGroup.set(centroid.groupId, []);
    }
    byGroup.get(centroid.groupId).push(centroid);
  });

  return { all, byGroup };
};

const getFeatureCandidatesForLabel = (featureGeometry, textEl, mode) => {
  if (!featureGeometry || featureGeometry.all.length === 0) return [];
  const labelGroupId = textEl.closest('g[id]')?.id || '';
  let candidates = featureGeometry.all;
  if (mode === 'linear' && labelGroupId) {
    const grouped = featureGeometry.byGroup.get(labelGroupId);
    if (grouped && grouped.length > 0) {
      candidates = grouped;
    }
  }
  return candidates;
};

const getDistanceThreshold = (mode, kind) => {
  if (mode === 'linear') return kind === 'embedded' ? 700 : 540;
  return kind === 'embedded' ? 520 : 420;
};

const computeCandidateDistance = (referencePoint, candidate, mode) => {
  const dx = candidate.x - referencePoint.x;
  const dy = candidate.y - referencePoint.y;
  return mode === 'linear' ? Math.abs(dx) + Math.abs(dy) * 0.6 : Math.hypot(dx, dy);
};

const assignFeatureIdsToLabels = (svg, labelElements, featureGeometry, mode) => {
  const assignments = new Map();
  if (!featureGeometry || featureGeometry.all.length === 0) return assignments;

  const featureIds = new Set(featureGeometry.all.map((feature) => feature.id));
  const usedFeatureIds = new Set();
  const labelMeta = labelElements
    .map((textEl) => {
      const referencePoint = getLabelReferencePoint(svg, textEl, mode);
      if (!isFinitePoint(referencePoint)) return null;
      const candidates = getFeatureCandidatesForLabel(featureGeometry, textEl, mode);
      return {
        textEl,
        referencePoint,
        candidates,
        candidateById: new Map(candidates.map((candidate) => [candidate.id, candidate])),
        kind: textEl.querySelector('textPath') ? 'embedded' : 'regular'
      };
    })
    .filter(Boolean);

  labelMeta.forEach((meta) => {
    const existingId = String(meta.textEl.getAttribute('data-label-feature-id') || '').trim();
    if (!existingId || !featureIds.has(existingId)) return;
    if (!meta.candidateById.has(existingId)) return;
    if (usedFeatureIds.has(existingId)) return;
    assignments.set(meta.textEl, existingId);
    usedFeatureIds.add(existingId);
  });

  const edges = [];
  labelMeta.forEach((meta) => {
    if (assignments.has(meta.textEl)) return;
    const threshold = getDistanceThreshold(mode, meta.kind);
    meta.candidates.forEach((candidate) => {
      if (usedFeatureIds.has(candidate.id)) return;
      const distance = computeCandidateDistance(meta.referencePoint, candidate, mode);
      if (!Number.isFinite(distance) || distance > threshold) return;
      edges.push({ textEl: meta.textEl, featureId: candidate.id, distance });
    });
  });
  edges.sort((a, b) => a.distance - b.distance);

  const assignedLabels = new Set(assignments.keys());
  edges.forEach((edge) => {
    if (assignedLabels.has(edge.textEl)) return;
    if (usedFeatureIds.has(edge.featureId)) return;
    assignments.set(edge.textEl, edge.featureId);
    assignedLabels.add(edge.textEl);
    usedFeatureIds.add(edge.featureId);
  });

  return assignments;
};

const buildContextKey = (svg, mode) => {
  const ids = Array.from(
    new Set(
      Array.from(svg.querySelectorAll(FEATURE_SELECTOR))
        .map((el) => el.getAttribute('id'))
        .filter((value) => value && value.trim() !== '')
    )
  ).sort();
  return `${mode}:${ids.join(',')}`;
};

export const createFeatureLabelActions = ({ state }) => {
  const {
    mode,
    form,
    filterMode,
    manualWhitelist,
    results,
    selectedResultIndex,
    svgContainer,
    skipCaptureBaseConfig,
    editableLabels,
    extractedFeatures,
    clickedFeature,
    labelTextScopeDialog,
    labelTextFeatureOverrides,
    labelTextBulkOverrides,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    labelOverrideContextKey,
    labelOverrideBuildWarning,
    globalLabelModeDialog,
    autoLabelReflowEnabled,
    labelReflowRequestSeq,
    labelReflowRequestReason,
    labelReflowForceRequestSeq,
    labelReflowForceRequestReason,
    labelReflowLastError
  } = state;

  const clearOverrides = () => {
    Object.keys(labelTextFeatureOverrides).forEach((key) => delete labelTextFeatureOverrides[key]);
    Object.keys(labelTextBulkOverrides).forEach((key) => delete labelTextBulkOverrides[key]);
    Object.keys(labelTextFeatureOverrideSources).forEach((key) => delete labelTextFeatureOverrideSources[key]);
    Object.keys(labelVisibilityOverrides).forEach((key) => delete labelVisibilityOverrides[key]);
    labelOverrideBuildWarning.value = '';
  };

  const serializeCurrentSvg = (svg) => {
    const index = selectedResultIndex.value;
    if (index < 0 || index >= results.value.length) return;
    const serializer = new XMLSerializer();
    const serialized = serializer.serializeToString(svg);
    if (results.value[index]?.content === serialized) return;
    skipCaptureBaseConfig.value = true;
    results.value[index] = { ...results.value[index], content: serialized };
  };

  const queueLabelReflow = (reason, force = false) => {
    labelReflowLastError.value = null;
    const normalizedReason = String(reason || 'label-edit');
    if (force) {
      labelReflowForceRequestReason.value = normalizedReason;
      labelReflowForceRequestSeq.value += 1;
      return;
    }
    if (!autoLabelReflowEnabled.value) return;
    labelReflowRequestReason.value = normalizedReason;
    labelReflowRequestSeq.value += 1;
  };

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return normalized === 'on' || normalized === 'off' ? normalized : 'default';
  };

  const isGlobalLabelsOff = () => {
    if (mode.value === 'circular') {
      const labelsMode = String(form.labels_mode || 'none').trim().toLowerCase();
      return labelsMode === 'none';
    }
    const linearLabels = String(form.show_labels_linear || 'none').trim().toLowerCase();
    return linearLabels === 'none';
  };

  const enableGlobalLabels = () => {
    if (mode.value === 'circular') {
      form.labels_mode = 'out';
      return;
    }
    form.show_labels_linear = 'all';
  };

  const closeGlobalLabelModeDialog = () => {
    globalLabelModeDialog.show = false;
    globalLabelModeDialog.featureId = '';
    globalLabelModeDialog.featureType = '';
    globalLabelModeDialog.resolve = null;
  };

  const requestGlobalLabelModeChoice = (featureId, featureType) =>
    new Promise((resolve) => {
      if (!featureId) {
        resolve('show_all');
        return;
      }
      globalLabelModeDialog.show = true;
      globalLabelModeDialog.featureId = String(featureId || '');
      globalLabelModeDialog.featureType = String(featureType || '');
      globalLabelModeDialog.resolve = resolve;
    });

  const handleGlobalLabelModeChoice = (choiceRaw) => {
    if (!globalLabelModeDialog.show) return;
    const resolver = globalLabelModeDialog.resolve;
    const normalizedChoice = choiceRaw === 'whitelist_only' ? 'whitelist_only' : 'show_all';
    closeGlobalLabelModeDialog();
    if (typeof resolver === 'function') {
      resolver(normalizedChoice);
    }
  };

  const ensureWhitelistHashRule = (featureTypeRaw, featureIdRaw) => {
    const featureType = String(featureTypeRaw || '').trim();
    const featureId = String(featureIdRaw || '').trim();
    if (!featureType || !featureId) return;
    const exists = manualWhitelist.some((rule) => {
      return (
        normalizeKeyToken(rule?.feat) === normalizeKeyToken(featureType) &&
        normalizeKeyToken(rule?.qual) === 'hash' &&
        normalizeKeyToken(rule?.key) === normalizeKeyToken(featureId)
      );
    });
    if (exists) return;
    manualWhitelist.push({ feat: featureType, qual: 'hash', key: featureId });
  };

  const applyGlobalLabelModeChoice = (choice, featureType, featureId) => {
    if (choice === 'whitelist_only') {
      filterMode.value = 'Whitelist';
      ensureWhitelistHashRule(featureType, featureId);
    } else {
      filterMode.value = 'None';
    }
    enableGlobalLabels();
  };

  const resetLabelsToSourceText = (svg) => {
    let changed = false;
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
      const sourceText = textEl.getAttribute('data-label-source-text');
      if (sourceText === null) return;
      if (getLabelText(textEl) === sourceText) return;
      setLabelText(textEl, sourceText);
      changed = true;
    });
    return changed;
  };

  const applyStoredOverridesToSvg = (svg) => {
    let changed = false;
    const labelElements = svg.querySelectorAll(EDITABLE_LABEL_SELECTOR);
    labelElements.forEach((textEl) => {
      const sourceText = textEl.getAttribute('data-label-source-text') || getLabelText(textEl);
      const featureId = textEl.getAttribute('data-label-feature-id');
      const currentText = getLabelText(textEl);
      const desiredText =
        (featureId ? labelTextFeatureOverrides[featureId] : undefined) ??
        labelTextBulkOverrides[sourceText];
      if (desiredText === undefined) return;
      if (currentText === desiredText) return;
      setLabelText(textEl, desiredText);
      changed = true;
    });
    return changed;
  };

  const refreshEditableList = (svg) => {
    const nextEntries = [];
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl, index) => {
      const key = textEl.getAttribute('data-label-key');
      if (!key) return;
      const text = getLabelText(textEl);
      const sourceText = textEl.getAttribute('data-label-source-text') || text;
      const featureId = textEl.getAttribute('data-label-feature-id') || '';
      const kind = textEl.querySelector('textPath') ? 'embedded' : 'regular';
      nextEntries.push({
        key,
        idx: index + 1,
        text,
        sourceText,
        featureId,
        kind,
        draftText: text
      });
    });
    editableLabels.value = nextEntries;
  };

  const getEditableLabelByFeatureId = (featureId) => {
    const target = normalizeKeyToken(featureId);
    if (!target) return null;
    return (
      editableLabels.value.find(
        (entry) => normalizeKeyToken(entry?.featureId) === target
      ) || null
    );
  };

  const syncClickedFeatureLabelState = () => {
    if (!clickedFeature.value) return;
    const featureId = String(clickedFeature.value.svg_id || clickedFeature.value.id || '').trim();
    const entry = getEditableLabelByFeatureId(featureId);
    const fallbackText =
      (featureId ? labelTextFeatureOverrides[featureId] : undefined) ||
      clickedFeature.value.labelText ||
      clickedFeature.value.label ||
      '';
    const fallbackSource =
      (featureId ? labelTextFeatureOverrideSources[featureId] : undefined) ||
      clickedFeature.value.labelSourceText ||
      clickedFeature.value.label ||
      '';
    const visibilityMode = featureId ? normalizeVisibilityMode(labelVisibilityOverrides[featureId]) : 'default';
    clickedFeature.value.labelKey = entry?.key || '';
    clickedFeature.value.labelText = entry?.text ?? fallbackText;
    clickedFeature.value.labelSourceText = entry?.sourceText ?? fallbackSource;
    clickedFeature.value.labelVisibility = visibilityMode;
    clickedFeature.value.hasEditableLabel = Boolean(entry || featureId);
    clickedFeature.value.labelUnavailableReason = entry || featureId
      ? ''
      : 'No editable feature label for this feature in current diagram.';
  };

  const syncLabelEditor = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const contextKey = buildContextKey(svg, mode.value);
    if (labelOverrideContextKey.value && labelOverrideContextKey.value !== contextKey) {
      clearOverrides();
      labelTextScopeDialog.show = false;
      closeGlobalLabelModeDialog();
    }
    labelOverrideContextKey.value = contextKey;

    const featureGeometry = collectFeatureGeometry(svg);
    const labelElements = collectEditableLabelElements(svg, mode.value);
    const featureAssignments = assignFeatureIdsToLabels(svg, labelElements, featureGeometry, mode.value);
    labelElements.forEach((textEl, index) => {
      textEl.style.cursor = 'text';
      textEl.setAttribute('data-label-editable', 'true');
      textEl.setAttribute('data-label-key', `label-${index + 1}`);
      const currentText = getLabelText(textEl);
      if (!textEl.hasAttribute('data-label-source-text')) {
        textEl.setAttribute('data-label-source-text', currentText);
      }
      const featureId = featureAssignments.get(textEl);
      if (featureId) {
        textEl.setAttribute('data-label-feature-id', featureId);
      } else {
        textEl.removeAttribute('data-label-feature-id');
      }
    });

    const changed = applyStoredOverridesToSvg(svg);
    refreshEditableList(svg);
    syncClickedFeatureLabelState();
    if (changed) {
      serializeCurrentSvg(svg);
    }
  };

  const requestLabelTextChangeByKey = (labelKey, nextTextRaw) => {
    if (!labelKey) return;
    const entry = editableLabels.value.find((candidate) => candidate.key === labelKey);
    if (!entry) return;
    const nextText = String(nextTextRaw ?? '');
    if (entry.text === nextText) return;

    labelTextScopeDialog.show = true;
    labelTextScopeDialog.labelKey = entry.key;
    labelTextScopeDialog.newText = nextText;
    labelTextScopeDialog.sourceText = entry.sourceText || entry.text;
    labelTextScopeDialog.featureId = entry.featureId || '';
    labelTextScopeDialog.matchingCount = editableLabels.value.filter(
      (candidate) => candidate.sourceText === (entry.sourceText || entry.text)
    ).length;
    if (clickedFeature.value?.labelKey === entry.key) {
      clickedFeature.value.labelText = nextText;
    }
  };

  const closeLabelTextScopeDialog = () => {
    labelTextScopeDialog.show = false;
    labelTextScopeDialog.labelKey = '';
    labelTextScopeDialog.newText = '';
    labelTextScopeDialog.sourceText = '';
    labelTextScopeDialog.featureId = '';
    labelTextScopeDialog.matchingCount = 0;
  };

  const requestLabelTextChangeByFeatureId = (featureId, nextTextRaw) => {
    const entry = getEditableLabelByFeatureId(featureId);
    if (!entry) return false;
    requestLabelTextChangeByKey(entry.key, nextTextRaw);
    return true;
  };

  const applyClickedFeatureVisibilityOverride = () => {
    if (!clickedFeature.value) return;
    const featureId = String(clickedFeature.value.svg_id || clickedFeature.value.id || '').trim();
    if (!featureId) return false;
    const nextMode = normalizeVisibilityMode(clickedFeature.value.labelVisibility);
    const previousMode = normalizeVisibilityMode(labelVisibilityOverrides[featureId]);
    if (nextMode === 'default') {
      if (!Object.prototype.hasOwnProperty.call(labelVisibilityOverrides, featureId)) return false;
      delete labelVisibilityOverrides[featureId];
      clickedFeature.value.labelVisibility = 'default';
      return previousMode !== 'default';
    }
    labelVisibilityOverrides[featureId] = nextMode;
    clickedFeature.value.labelVisibility = nextMode;
    return previousMode !== nextMode;
  };

  const applyDirectFeatureLabelOverride = (featureId, labelTextRaw, sourceTextRaw, baselineTextRaw) => {
    const featureIdKey = String(featureId || '').trim();
    if (!featureIdKey) return false;
    const nextText = String(labelTextRaw ?? '');
    const baselineText = String(baselineTextRaw ?? '');
    const hasExistingOverride = Object.prototype.hasOwnProperty.call(labelTextFeatureOverrides, featureIdKey);
    const prevText = Object.prototype.hasOwnProperty.call(labelTextFeatureOverrides, featureIdKey)
      ? String(labelTextFeatureOverrides[featureIdKey] ?? '')
      : undefined;
    const sourceText = String(sourceTextRaw ?? '').trim();
    let changed = false;

    if (!hasExistingOverride && nextText === baselineText) {
      return false;
    }

    if (hasExistingOverride && nextText === baselineText) {
      delete labelTextFeatureOverrides[featureIdKey];
      delete labelTextFeatureOverrideSources[featureIdKey];
      return true;
    }

    if (prevText !== nextText || !hasExistingOverride) {
      labelTextFeatureOverrides[featureIdKey] = nextText;
      changed = true;
    }
    if (sourceText) {
      const prevSource = String(labelTextFeatureOverrideSources[featureIdKey] ?? '');
      if (prevSource !== sourceText) changed = true;
      labelTextFeatureOverrideSources[featureIdKey] = sourceText;
    }
    return changed;
  };

  const applyDirectTextToCurrentSvg = (featureId, nextText) => {
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const entry = getEditableLabelByFeatureId(featureId);
    if (!entry?.key) return false;
    const targetEl = svg.querySelector(`text[data-label-key="${CSS.escape(entry.key)}"]`);
    if (!targetEl) return false;
    const currentText = getLabelText(targetEl);
    if (currentText === nextText) return false;
    setLabelText(targetEl, nextText);
    serializeCurrentSvg(svg);
    syncLabelEditor();
    return true;
  };

  const updateClickedFeatureLabelText = async () => {
    if (!clickedFeature.value) return;
    const featureId = String(clickedFeature.value.svg_id || clickedFeature.value.id || '').trim();
    if (!featureId) return;

    const featureType = String(clickedFeature.value.feat?.type || '').trim();
    const nextText = String(clickedFeature.value.labelText ?? '');
    const sourceText = String(clickedFeature.value.labelSourceText || clickedFeature.value.label || '');
    const baselineText = sourceText;
    const visibilityChanged = applyClickedFeatureVisibilityOverride();
    const textChanged = applyDirectFeatureLabelOverride(featureId, nextText, sourceText, baselineText);

    if (clickedFeature.value.hasEditableLabel) {
      applyDirectTextToCurrentSvg(featureId, nextText);
    }

    const requiresGlobalSelection = isGlobalLabelsOff() && (visibilityChanged || textChanged);
    if (requiresGlobalSelection) {
      const choice = await requestGlobalLabelModeChoice(featureId, featureType);
      applyGlobalLabelModeChoice(choice, featureType, featureId);
      queueLabelReflow('global-off-label-apply', true);
      return;
    }

    if (visibilityChanged || (!clickedFeature.value.hasEditableLabel && textChanged)) {
      queueLabelReflow('label-visibility-apply', true);
      return;
    }

    if (textChanged) {
      queueLabelReflow('apply');
    }
  };

  const handleLabelTextScopeChoice = (choice) => {
    if (choice === 'cancel' || !labelTextScopeDialog.show) {
      closeLabelTextScopeDialog();
      return;
    }
    if (!svgContainer.value) {
      closeLabelTextScopeDialog();
      return;
    }
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) {
      closeLabelTextScopeDialog();
      return;
    }

    const targetKey = String(labelTextScopeDialog.labelKey || '');
    const sourceText = String(labelTextScopeDialog.sourceText || '');
    const newText = String(labelTextScopeDialog.newText ?? '');
    const featureId = String(labelTextScopeDialog.featureId || '');

    if (choice === 'all') {
      let matchedCount = 0;
      let hasUntrackableMatch = false;
      svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
        const candidateSource = textEl.getAttribute('data-label-source-text') || '';
        if (candidateSource !== sourceText) return;
        matchedCount += 1;
        const candidateFeatureId = textEl.getAttribute('data-label-feature-id');
        if (candidateFeatureId) {
          labelTextFeatureOverrides[candidateFeatureId] = newText;
          labelTextFeatureOverrideSources[candidateFeatureId] = candidateSource;
        } else {
          hasUntrackableMatch = true;
        }
        setLabelText(textEl, newText);
      });
      if (hasUntrackableMatch || matchedCount === 0) {
        labelTextBulkOverrides[sourceText] = newText;
      } else if (Object.prototype.hasOwnProperty.call(labelTextBulkOverrides, sourceText)) {
        delete labelTextBulkOverrides[sourceText];
      }
    } else if (choice === 'single' && targetKey) {
      const targetEl = svg.querySelector(`text[data-label-key="${CSS.escape(targetKey)}"]`);
      if (targetEl) {
        setLabelText(targetEl, newText);
      }
      if (featureId) {
        labelTextFeatureOverrides[featureId] = newText;
        labelTextFeatureOverrideSources[featureId] = sourceText;
      }
    } else {
      closeLabelTextScopeDialog();
      return;
    }

    serializeCurrentSvg(svg);
    closeLabelTextScopeDialog();
    syncLabelEditor();
    queueLabelReflow('apply');
  };

  const resetAllLabelTextOverrides = () => {
    clearOverrides();
    if (!svgContainer.value) {
      labelOverrideContextKey.value = '';
      editableLabels.value = [];
      closeLabelTextScopeDialog();
      closeGlobalLabelModeDialog();
      return;
    }
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    resetLabelsToSourceText(svg);

    closeLabelTextScopeDialog();
    closeGlobalLabelModeDialog();
    serializeCurrentSvg(svg);
    syncLabelEditor();
    queueLabelReflow('reset');
  };

  const getEntryMeta = (entry, metadataByFeatureId) => {
    const featureIdKey = normalizeKeyToken(entry?.featureId);
    const metadata = metadataByFeatureId.get(featureIdKey) || null;
    const record = String(metadata?.record || '').trim();
    const location = String(metadata?.location || '').trim();
    const position = normalizePositionToken(metadata?.position || '');
    const featureType = String(metadata?.featureType || '').trim();
    const qualifiers =
      metadata && metadata.qualifiers && typeof metadata.qualifiers === 'object'
        ? metadata.qualifiers
        : {};
    const recordLocation = record && position && position !== '*' ? `${record}:${position}` : '';
    return {
      record,
      location,
      position,
      featureType,
      qualifiers,
      recordLocation
    };
  };

  const getQualifierValuesForEntry = (entry, entryMeta, qualifierKeyRaw) => {
    const qualifierKey = String(qualifierKeyRaw || '').trim().toLowerCase();
    if (qualifierKey === 'label') return [String(entry?.sourceText || '')];
    if (qualifierKey === 'hash') return [String(entry?.featureId || '')];
    if (qualifierKey === 'location') return entryMeta.location ? [entryMeta.location] : [];
    if (qualifierKey === 'record_location') return entryMeta.recordLocation ? [entryMeta.recordLocation] : [];
    const values = entryMeta.qualifiers[qualifierKey];
    if (!Array.isArray(values)) return [];
    return values.map((value) => String(value));
  };

  const rowMatchesEntry = (row, entry, entryMeta) => {
    if (!wildcardOrExactMatch(row.recordId, entryMeta.record)) return false;
    if (!wildcardOrExactMatch(row.featureType, entryMeta.featureType)) return false;
    const values = getQualifierValuesForEntry(entry, entryMeta, row.qualifier);
    if (!values.length) return false;
    return values.some((value) => row.qualifierValuePattern.test(String(value || '')));
  };

  const loadLabelOverrideTable = async (event) => {
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;

    try {
      const text = await file.text();
      const rows = parseLabelOverrideTsv(text);

      if (!svgContainer.value) {
        clearOverrides();
        labelOverrideContextKey.value = '';
        editableLabels.value = [];
        closeLabelTextScopeDialog();
        closeGlobalLabelModeDialog();
        window.alert(`Loaded ${rows.length} row(s). No diagram is currently displayed.`);
        return;
      }

      const svg = svgContainer.value.querySelector('svg');
      if (!svg) {
        clearOverrides();
        labelOverrideContextKey.value = '';
        editableLabels.value = [];
        closeLabelTextScopeDialog();
        closeGlobalLabelModeDialog();
        window.alert(`Loaded ${rows.length} row(s). No diagram is currently displayed.`);
        return;
      }

      syncLabelEditor();
      const metadataByFeatureId = buildFeatureMetadataMap(extractedFeatures.value);
      const operations = [];
      editableLabels.value.forEach((entry) => {
        const entryMeta = getEntryMeta(entry, metadataByFeatureId);
        const matchedRow = rows.find((row) => rowMatchesEntry(row, entry, entryMeta));
        if (!matchedRow) return;
        operations.push({
          key: String(entry.key || ''),
          featureId: String(entry.featureId || ''),
          sourceText: String(entry.sourceText || ''),
          nextText: String(matchedRow.labelText ?? ''),
          isGlobalLabelRule: Boolean(matchedRow.isGlobalLabelRule)
        });
      });

      clearOverrides();
      const resetChanged = resetLabelsToSourceText(svg);

      let appliedCount = 0;
      let skippedNonTrackableCount = 0;
      const nextFeatureOverrides = {};
      const nextBulkOverrides = {};

      operations.forEach((operation) => {
        if (!operation.key) return;
        const target = svg.querySelector(`text[data-label-key="${CSS.escape(operation.key)}"]`);
        if (!target) return;
        setLabelText(target, operation.nextText);
        appliedCount += 1;

        if (operation.isGlobalLabelRule) {
          if (operation.sourceText) {
            nextBulkOverrides[operation.sourceText] = operation.nextText;
          }
          return;
        }
        if (operation.featureId) {
          nextFeatureOverrides[operation.featureId] = operation.nextText;
        } else {
          skippedNonTrackableCount += 1;
        }
      });

      Object.entries(nextFeatureOverrides).forEach(([featureId, labelText]) => {
        labelTextFeatureOverrides[featureId] = labelText;
      });
      operations.forEach((operation) => {
        if (!operation.featureId || !operation.sourceText) return;
        labelTextFeatureOverrideSources[operation.featureId] = operation.sourceText;
      });
      Object.entries(nextBulkOverrides).forEach(([sourceText, labelText]) => {
        if (!sourceText) return;
        labelTextBulkOverrides[sourceText] = labelText;
      });

      closeLabelTextScopeDialog();
      closeGlobalLabelModeDialog();
      if (resetChanged || appliedCount > 0) {
        serializeCurrentSvg(svg);
      }
      syncLabelEditor();
      queueLabelReflow('load');

      let message = `Loaded ${rows.length} row(s). Applied to ${appliedCount} label(s).`;
      if (skippedNonTrackableCount > 0) {
        message += ` ${skippedNonTrackableCount} match(es) lacked a feature hash and were not tracked for re-export.`;
      }
      window.alert(message);
    } catch (error) {
      console.error('Failed to load label TSV:', error);
      window.alert(`Failed to load label TSV. ${error?.message || 'Please check the 5-column TSV format.'}`);
    } finally {
      if (input) input.value = '';
    }
  };

  const downloadLabelOverrideTable = () => {
    const { rows, skippedFeatureCount, skippedFeatureSourceCount, skippedMissingSourceCount } = buildLabelOverrideRows(
      labelTextFeatureOverrides,
      labelTextBulkOverrides,
      {
        editableLabels: editableLabels.value,
        extractedFeatures: extractedFeatures.value,
        featureOverrideSources: labelTextFeatureOverrideSources,
        visibilityOverrides: labelVisibilityOverrides
      }
    );
    if (rows.length === 0) {
      if (skippedFeatureCount > 0 || skippedFeatureSourceCount > 0 || skippedMissingSourceCount > 0) {
        window.alert(
          `No exportable label edits. ${skippedFeatureCount} invalid feature hash edit(s) and ` +
          `${skippedFeatureSourceCount} source label lookup failure(s) were skipped.` +
          (skippedMissingSourceCount > 0
            ? ` ${skippedMissingSourceCount} feature override row(s) had missing source label context.`
            : '')
        );
      } else {
        window.alert('No label edits to export.');
      }
      return;
    }

    const selectedIdx = selectedResultIndex.value;
    const resultName =
      selectedIdx >= 0 && selectedIdx < results.value.length
        ? String(results.value[selectedIdx]?.name || '')
        : '';
    const outputName = `${makeSafeFilename(resultName, 'gbdraw')}.label_table.tsv`;
    downloadTextFile(outputName, `${rows.join('\n')}\n`);
    if (skippedFeatureCount > 0 || skippedFeatureSourceCount > 0 || skippedMissingSourceCount > 0) {
      window.alert(
        `Exported ${rows.length} row(s). ${skippedFeatureCount} invalid feature hash edit(s) and ` +
        `${skippedFeatureSourceCount} source label lookup failure(s) were skipped.` +
        (skippedMissingSourceCount > 0
          ? ` ${skippedMissingSourceCount} feature override row(s) had missing source label context.`
          : '')
      );
    }
  };

  return {
    downloadLabelOverrideTable,
    loadLabelOverrideTable,
    getEditableLabelByFeatureId,
    handleGlobalLabelModeChoice,
    handleLabelTextScopeChoice,
    requestLabelTextChangeByFeatureId,
    requestLabelTextChangeByKey,
    resetAllLabelTextOverrides,
    syncClickedFeatureLabelState,
    syncLabelEditor,
    updateClickedFeatureLabelText
  };
};
