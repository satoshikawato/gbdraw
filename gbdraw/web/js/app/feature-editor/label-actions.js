import {
  buildFeatureMetadataMap,
  buildFeatureUniquenessIndex,
  buildLabelOverrideRows,
  parseLabelOverrideTsv,
  selectStableFeatureKey
} from './label-override-table.js';
import { FEATURE_SELECTOR, getFeatureIdentity } from './svg-actions.js';
import { serializeCleanSvg } from '../../services/svg-serialization.js';
import { downloadTextFile } from '../../services/text-download.js';
import { COMPARISON_LEGEND_SELECTOR } from '../legend/utils.js';

export const EXCLUDED_GROUP_SELECTOR = [
  '#legend',
  '#feature_legend',
  COMPARISON_LEGEND_SELECTOR,
  '#horizontal_legend',
  '#vertical_legend',
  '#length_bar',
  'g[data-gbdraw-slot-renderer="ticks"]',
  'g[id="tick"]',
  'g[id^="tick_"]'
].join(', ');
const EDITABLE_LABEL_SELECTOR = 'text[data-label-editable="true"]';
const LABEL_VISIBILITY_PREVIEW_ATTRIBUTE = 'data-gbdraw-label-visibility-preview';

const toNumber = (value, fallback = 0) => {
  const parsed = Number.parseFloat(value);
  return Number.isFinite(parsed) ? parsed : fallback;
};

const normalizeKeyToken = (value) => String(value ?? '').trim().toLowerCase();
const escapeRegexLiteral = (value) => String(value ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
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

const getPhasedCircularFeatureLine = (svg, textEl) => {
  const textGroup = textEl.closest('g[id="label_text"], g[id^="label_text_"]');
  if (!textGroup) return null;
  const textGroups = Array.from(svg.querySelectorAll('g[id="label_text"], g[id^="label_text_"]'));
  const leaderGroups = Array.from(svg.querySelectorAll('g[id="label_leaders"], g[id^="label_leaders_"]'));
  const groupIndex = textGroups.indexOf(textGroup);
  if (groupIndex < 0 || groupIndex >= leaderGroups.length) return null;
  const textIndex = Array.from(textGroup.querySelectorAll('text')).indexOf(textEl);
  if (textIndex < 0) return null;
  const leaderLines = Array.from(leaderGroups[groupIndex].querySelectorAll('line'));
  return leaderLines[(2 * textIndex) + 1] || null;
};

const getCircularFeatureAnchor = (svg, textEl) => {
  const inLabelsGroup = Boolean(textEl.closest('g[id="labels"], g[id^="labels_"]'));
  let featureLine = null;
  if (inLabelsGroup) {
    const prev = textEl.previousElementSibling;
    const prev2 = prev ? prev.previousElementSibling : null;
    const lines = [prev, prev2].filter(
      (candidate) => candidate && candidate.tagName && candidate.tagName.toLowerCase() === 'line'
    );
    featureLine = lines[0] || null;
  } else {
    featureLine = getPhasedCircularFeatureLine(svg, textEl);
  }
  if (!featureLine) return null;

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
    svg.querySelectorAll('g[id="labels"] text, g[id^="labels_"] text').forEach((textEl) => {
      if (hasExcludedAncestor(textEl)) return;
      labels.add(textEl);
    });
    svg.querySelectorAll('g[id="label_text"] text, g[id^="label_text_"] text').forEach((textEl) => {
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
    const id = getFeatureIdentity(el);
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
        .map((el) => getFeatureIdentity(el))
        .filter((value) => value && value.trim() !== '')
    )
  ).sort();
  return `${mode}:${ids.join(',')}`;
};

const replaceObjectContents = (target, source) => {
  Object.keys(target || {}).forEach((key) => delete target[key]);
  Object.entries(source || {}).forEach(([key, value]) => {
    target[key] = value;
  });
};

export const createFeatureLabelActions = ({
  state,
  previewRuntime = null,
  serializeSvg = serializeCleanSvg
}) => {
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
    if (previewRuntime?.markActiveResultDirty?.('feature-label')) {
      skipCaptureBaseConfig.value = true;
      previewRuntime.flushActiveResult?.();
      return;
    }
    const index = selectedResultIndex.value;
    if (index < 0 || index >= results.value.length) return;
    const serialized = serializeSvg(svg);
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
    const normalizedChoice = choiceRaw === 'cancel'
      ? 'cancel'
      : (choiceRaw === 'whitelist_only' ? 'whitelist_only' : 'show_all');
    closeGlobalLabelModeDialog();
    if (typeof resolver === 'function') {
      resolver(normalizedChoice);
    }
  };

  const ensureWhitelistRuleForFeature = (featureTypeRaw, featureIdRaw) => {
    const featureType = String(featureTypeRaw || '').trim();
    const featureId = String(featureIdRaw || '').trim();
    if (!featureType || !featureId) return;
    const metadataByFeatureId = buildFeatureMetadataMap(extractedFeatures.value);
    const metadata = metadataByFeatureId.get(normalizeKeyToken(featureId));
    const selector = selectStableFeatureKey(
      {
        featureId,
        record: metadata?.record || '',
        featureType: metadata?.featureType || featureType,
        position: metadata?.position || '',
        qualifiers: metadata?.qualifiers || {}
      },
      buildFeatureUniquenessIndex(extractedFeatures.value)
    );
    const ruleFeatureType = String(metadata?.featureType || featureType).trim();
    const ruleQualifier = String(selector?.qualifier || 'hash').trim().toLowerCase();
    const ruleKey = `^${escapeRegexLiteral(String(selector?.value || featureId).trim())}$`;
    if (!ruleFeatureType || !ruleQualifier || !ruleKey) return;
    const exists = manualWhitelist.some((rule) => {
      return (
        normalizeKeyToken(rule?.feat) === normalizeKeyToken(ruleFeatureType) &&
        normalizeKeyToken(rule?.qual) === normalizeKeyToken(ruleQualifier) &&
        normalizeKeyToken(rule?.key) === normalizeKeyToken(ruleKey)
      );
    });
    if (exists) return;
    manualWhitelist.push({ feat: ruleFeatureType, qual: ruleQualifier, key: ruleKey });
  };

  const applyGlobalLabelModeChoice = (choice, featureType, featureId) => {
    if (choice === 'whitelist_only') {
      filterMode.value = 'Whitelist';
      ensureWhitelistRuleForFeature(featureType, featureId);
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

  const applyLabelVisibilityPreview = (textEl, modeRaw) => {
    if (!textEl) return false;
    const visibilityMode = normalizeVisibilityMode(modeRaw);
    const previewHidden = textEl.hasAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
    if (visibilityMode === 'off') {
      let changed = false;
      if (!previewHidden) {
        textEl.setAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE, 'off');
        changed = true;
      }
      if (textEl.getAttribute('display') !== 'none') {
        textEl.setAttribute('display', 'none');
        changed = true;
      }
      return changed;
    }
    if (visibilityMode === 'on') {
      const changed = previewHidden || textEl.getAttribute('display') === 'none';
      textEl.removeAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
      textEl.removeAttribute('display');
      return changed;
    }
    if (!previewHidden) return false;
    textEl.removeAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
    textEl.removeAttribute('display');
    return true;
  };

  const applyStoredVisibilityOverridesToSvg = (svg) => {
    let changed = false;
    svg.querySelectorAll(EDITABLE_LABEL_SELECTOR).forEach((textEl) => {
      const featureId = String(
        textEl.getAttribute('data-label-feature-id') || ''
      ).trim();
      const visibilityMode = featureId
        ? labelVisibilityOverrides[featureId]
        : 'default';
      changed = applyLabelVisibilityPreview(textEl, visibilityMode) || changed;
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

    const textChanged = applyStoredOverridesToSvg(svg);
    const visibilityChanged = applyStoredVisibilityOverridesToSvg(svg);
    refreshEditableList(svg);
    syncClickedFeatureLabelState();
    if (textChanged || visibilityChanged) {
      serializeCurrentSvg(svg);
    }
  };

  const reconcileLabelOverrides = () => {
    const svg = svgContainer.value?.querySelector?.('svg');
    if (!svg) return false;
    const resetChanged = resetLabelsToSourceText(svg);
    const overrideChanged = applyStoredOverridesToSvg(svg);
    const visibilityChanged = applyStoredVisibilityOverridesToSvg(svg);
    refreshEditableList(svg);
    syncClickedFeatureLabelState();
    if (resetChanged || overrideChanged || visibilityChanged) serializeCurrentSvg(svg);
    return resetChanged || overrideChanged || visibilityChanged;
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

  const applyDirectVisibilityToCurrentSvg = (featureId, visibilityMode) => {
    if (!svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;
    const entry = getEditableLabelByFeatureId(featureId);
    if (!entry?.key) return false;
    const targetEl = svg.querySelector(`text[data-label-key="${CSS.escape(entry.key)}"]`);
    if (!applyLabelVisibilityPreview(targetEl, visibilityMode)) return false;
    serializeCurrentSvg(svg);
    return true;
  };

  const updateClickedFeatureLabelVisibility = async () => {
    if (!clickedFeature.value) return false;
    const featureId = String(clickedFeature.value.svg_id || clickedFeature.value.id || '').trim();
    if (!featureId) return false;
    const visibilityMode = normalizeVisibilityMode(clickedFeature.value.labelVisibility);
    if (!applyClickedFeatureVisibilityOverride()) return false;
    applyDirectVisibilityToCurrentSvg(featureId, visibilityMode);
    if (visibilityMode === 'on' && isGlobalLabelsOff()) {
      const featureType = String(clickedFeature.value?.feat?.type || '').trim();
      const choice = await requestGlobalLabelModeChoice(featureId, featureType);
      if (choice === 'cancel') return false;
      applyGlobalLabelModeChoice(choice, featureType, featureId);
      queueLabelReflow('label-visibility-apply', true);
    } else {
      queueLabelReflow('label-visibility-apply', visibilityMode === 'on');
    }
    return true;
  };

  const resetAllLabelTextOverrides = () => {
    clearOverrides();
    if (!svgContainer.value) {
      labelOverrideContextKey.value = '';
      editableLabels.value = [];
      closeGlobalLabelModeDialog();
      return;
    }
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;
    resetLabelsToSourceText(svg);
    applyStoredVisibilityOverridesToSvg(svg);

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

  const elementHasVisibilityPreview = (element) => (
    element.hasAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE)
  );

  const planLabelOverrideImport = (svg, rows) => {
    const featureGeometry = collectFeatureGeometry(svg);
    const labelElements = collectEditableLabelElements(svg, mode.value);
    const featureAssignments = assignFeatureIdsToLabels(
      svg,
      labelElements,
      featureGeometry,
      mode.value
    );
    const entries = labelElements.map((element, index) => {
      const text = getLabelText(element);
      const sourceText = element.getAttribute('data-label-source-text') ?? text;
      const featureId = featureAssignments.get(element)
        || String(element.getAttribute('data-label-feature-id') || '').trim();
      return {
        element,
        key: `label-${index + 1}`,
        text,
        sourceText,
        featureId,
        kind: element.querySelector('textPath') ? 'embedded' : 'regular'
      };
    });
    const metadataByFeatureId = buildFeatureMetadataMap(extractedFeatures.value);
    const nextFeatureOverrides = {};
    const nextBulkOverrides = {};
    const nextFeatureOverrideSources = {};
    let appliedCount = 0;
    let skippedNonTrackableCount = 0;

    const preparedEntries = entries.map((entry) => {
      const entryMeta = getEntryMeta(entry, metadataByFeatureId);
      const matchedRow = rows.find((row) => rowMatchesEntry(row, entry, entryMeta));
      let nextText = entry.sourceText;
      if (matchedRow) {
        nextText = String(matchedRow.labelText ?? '');
        appliedCount += 1;
        if (matchedRow.isGlobalLabelRule) {
          if (entry.sourceText) nextBulkOverrides[entry.sourceText] = nextText;
        } else if (entry.featureId) {
          nextFeatureOverrides[entry.featureId] = nextText;
        } else {
          skippedNonTrackableCount += 1;
        }
        if (entry.featureId && entry.sourceText) {
          nextFeatureOverrideSources[entry.featureId] = entry.sourceText;
        }
      }
      return {
        ...entry,
        nextText,
        clearVisibilityPreview: elementHasVisibilityPreview(entry.element)
      };
    });

    return {
      appliedCount,
      contextKey: buildContextKey(svg, mode.value),
      nextBulkOverrides,
      nextFeatureOverrides,
      nextFeatureOverrideSources,
      preparedEntries,
      skippedNonTrackableCount
    };
  };

  const captureLabelImportState = (preparedEntries) => {
    const touchedAttributes = [
      'data-label-editable',
      'data-label-key',
      'data-label-source-text',
      'data-label-feature-id',
      LABEL_VISIBILITY_PREVIEW_ATTRIBUTE,
      'display',
      'style'
    ];
    const labelSnapshots = preparedEntries.map(({ element }) => ({
      element,
      text: getLabelText(element),
      attributes: Object.fromEntries(
        touchedAttributes.map((name) => [name, element.getAttribute(name)])
      )
    }));
    const resultsReference = results.value;
    const resultsSnapshot = [...resultsReference];
    const editableLabelsSnapshot = editableLabels.value;
    const clickedFeatureReference = clickedFeature.value;
    const clickedFeatureSnapshot = clickedFeatureReference
      ? { ...clickedFeatureReference }
      : null;
    const runtime = previewRuntime?.getActiveRuntime?.() || null;
    const runtimeSnapshot = runtime
      ? {
          dirty: Boolean(runtime.dirty),
          dirtyReasons: runtime.dirtyReasons,
          dirtyReasonValues: [...(runtime.dirtyReasons || [])]
        }
      : null;
    const globalDialogSnapshot = {
      show: globalLabelModeDialog.show,
      featureId: globalLabelModeDialog.featureId,
      featureType: globalLabelModeDialog.featureType,
      resolve: globalLabelModeDialog.resolve
    };
    const featureOverridesSnapshot = { ...labelTextFeatureOverrides };
    const bulkOverridesSnapshot = { ...labelTextBulkOverrides };
    const featureOverrideSourcesSnapshot = { ...labelTextFeatureOverrideSources };
    const visibilityOverridesSnapshot = { ...labelVisibilityOverrides };
    const contextKeySnapshot = labelOverrideContextKey.value;
    const buildWarningSnapshot = labelOverrideBuildWarning.value;
    const skipCaptureSnapshot = skipCaptureBaseConfig.value;

    return () => {
      labelSnapshots.forEach(({ element, text, attributes }) => {
        setLabelText(element, text);
        Object.entries(attributes).forEach(([name, value]) => {
          if (value === null) element.removeAttribute(name);
          else element.setAttribute(name, value);
        });
      });
      replaceObjectContents(labelTextFeatureOverrides, featureOverridesSnapshot);
      replaceObjectContents(labelTextBulkOverrides, bulkOverridesSnapshot);
      replaceObjectContents(labelTextFeatureOverrideSources, featureOverrideSourcesSnapshot);
      replaceObjectContents(labelVisibilityOverrides, visibilityOverridesSnapshot);
      labelOverrideContextKey.value = contextKeySnapshot;
      labelOverrideBuildWarning.value = buildWarningSnapshot;
      editableLabels.value = editableLabelsSnapshot;
      if (clickedFeatureReference && clickedFeatureSnapshot) {
        Object.keys(clickedFeatureReference).forEach((key) => delete clickedFeatureReference[key]);
        Object.assign(clickedFeatureReference, clickedFeatureSnapshot);
      }
      clickedFeature.value = clickedFeatureReference;
      Object.assign(globalLabelModeDialog, globalDialogSnapshot);
      if (!runtime && previewRuntime?.getActiveRuntime?.()) {
        previewRuntime.clearActiveRuntime?.();
      } else if (runtime && runtimeSnapshot) {
        runtime.dirty = runtimeSnapshot.dirty;
        runtime.dirtyReasons = runtimeSnapshot.dirtyReasons;
        runtime.dirtyReasons?.clear?.();
        runtimeSnapshot.dirtyReasonValues.forEach((reason) => runtime.dirtyReasons?.add?.(reason));
      }
      if (results.value !== resultsReference) results.value = resultsReference;
      resultsReference.splice(0, resultsReference.length, ...resultsSnapshot);
      skipCaptureBaseConfig.value = skipCaptureSnapshot;
    };
  };

  const commitLabelOverrideImport = (svg, plan) => {
    const currentElements = collectEditableLabelElements(svg, mode.value);
    if (
      currentElements.length !== plan.preparedEntries.length
      || currentElements.some((element, index) => element !== plan.preparedEntries[index].element)
    ) {
      throw new Error('Editable labels changed before the TSV import started.');
    }
    const restore = captureLabelImportState(plan.preparedEntries);
    try {
      plan.preparedEntries.forEach((entry) => {
        entry.element.style.cursor = 'text';
        entry.element.setAttribute('data-label-editable', 'true');
        entry.element.setAttribute('data-label-key', entry.key);
        entry.element.setAttribute('data-label-source-text', entry.sourceText);
        if (entry.featureId) entry.element.setAttribute('data-label-feature-id', entry.featureId);
        else entry.element.removeAttribute('data-label-feature-id');
        if (entry.clearVisibilityPreview) {
          entry.element.removeAttribute(LABEL_VISIBILITY_PREVIEW_ATTRIBUTE);
          entry.element.removeAttribute('display');
        }
        setLabelText(entry.element, entry.nextText);
      });
      replaceObjectContents(labelTextFeatureOverrides, plan.nextFeatureOverrides);
      replaceObjectContents(labelTextBulkOverrides, plan.nextBulkOverrides);
      replaceObjectContents(labelTextFeatureOverrideSources, plan.nextFeatureOverrideSources);
      replaceObjectContents(labelVisibilityOverrides, {});
      labelOverrideBuildWarning.value = '';
      labelOverrideContextKey.value = plan.contextKey;
      closeGlobalLabelModeDialog();
      serializeCurrentSvg(svg);
      refreshEditableList(svg);
      syncClickedFeatureLabelState();
    } catch (error) {
      restore();
      throw error;
    }
  };

  const loadLabelOverrideTable = async (event) => {
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;

    let message = '';
    let failure = null;
    try {
      const text = await file.text();
      const rows = parseLabelOverrideTsv(text);

      if (!svgContainer.value) {
        clearOverrides();
        labelOverrideContextKey.value = '';
        editableLabels.value = [];
        closeGlobalLabelModeDialog();
        message = `Loaded ${rows.length} row(s). No diagram is currently displayed.`;
      } else {
        const svg = svgContainer.value.querySelector('svg');
        if (!svg) {
          clearOverrides();
          labelOverrideContextKey.value = '';
          editableLabels.value = [];
          closeGlobalLabelModeDialog();
          message = `Loaded ${rows.length} row(s). No diagram is currently displayed.`;
        } else {
          const plan = planLabelOverrideImport(svg, rows);
          commitLabelOverrideImport(svg, plan);
          queueLabelReflow('load');

          message = `Loaded ${rows.length} row(s). Applied to ${plan.appliedCount} label(s).`;
          if (plan.skippedNonTrackableCount > 0) {
            message += ` ${plan.skippedNonTrackableCount} match(es) lacked a feature key and were not tracked for re-export.`;
          }
        }
      }
    } catch (error) {
      console.error('Failed to load label TSV:', error);
      failure = `Failed to load label TSV. ${error?.message || 'Please check the 5-column TSV format.'}`;
    } finally {
      if (input) input.value = '';
    }
    window.alert(failure || message);
  };

  const downloadLabelOverrideTable = () => {
    const { rows, skippedFeatureCount, skippedFeatureSourceCount, skippedMissingSourceCount, fallbackHashCount } = buildLabelOverrideRows(
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
          `No exportable label edits. ${skippedFeatureCount} invalid feature-targeted edit(s) and ` +
          `${skippedFeatureSourceCount} source label lookup failure(s) were skipped.` +
          (skippedMissingSourceCount > 0
            ? ` ${skippedMissingSourceCount} feature override row(s) had missing source label context.`
            : '') +
          (fallbackHashCount > 0
            ? ` ${fallbackHashCount} row(s) fell back to hash because non-hash keys were not unique.`
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
    downloadTextFile(outputName, `${rows.join('\n')}\n`, 'text/tab-separated-values');
    if (skippedFeatureCount > 0 || skippedFeatureSourceCount > 0 || skippedMissingSourceCount > 0 || fallbackHashCount > 0) {
      window.alert(
        `Exported ${rows.length} row(s). ${skippedFeatureCount} invalid feature-targeted edit(s) and ` +
        `${skippedFeatureSourceCount} source label lookup failure(s) were skipped.` +
        (skippedMissingSourceCount > 0
          ? ` ${skippedMissingSourceCount} feature override row(s) had missing source label context.`
          : '') +
        (fallbackHashCount > 0
          ? ` ${fallbackHashCount} row(s) fell back to hash because non-hash keys were not unique.`
          : '')
      );
    }
  };

  return {
    downloadLabelOverrideTable,
    loadLabelOverrideTable,
    getEditableLabelByFeatureId,
    handleGlobalLabelModeChoice,
    reconcileLabelOverrides,
    resetAllLabelTextOverrides,
    syncClickedFeatureLabelState,
    syncLabelEditor,
    updateClickedFeatureLabelVisibility
  };
};
