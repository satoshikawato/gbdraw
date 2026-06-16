import { resolveColorToHex } from '../color-utils.js';
import { getFeatureCaption } from '../feature-utils.js';

export const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
export const FEATURE_SELECTOR = [
  `path[${FEATURE_ID_ATTRIBUTE}]`,
  `polygon[${FEATURE_ID_ATTRIBUTE}]`,
  `rect[${FEATURE_ID_ATTRIBUTE}]`,
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');

export const getFeatureIdentity = (element) =>
  String(
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    element?.id ||
    ''
  ).trim();

const featureElementIndexCache = new WeakMap();

export const buildFeatureElementIndex = (svg, { markCursor = false } = {}) => {
  const indexed = new Map();
  if (!svg) return indexed;

  Array.from(svg.querySelectorAll(FEATURE_SELECTOR)).forEach((element) => {
    const id = getFeatureIdentity(element);
    if (!id) return;
    if (!indexed.has(id)) indexed.set(id, []);
    indexed.get(id).push(element);
    if (markCursor && element?.style) element.style.cursor = 'pointer';
  });
  featureElementIndexCache.set(svg, indexed);
  return indexed;
};

export const getFeatureElementIndex = (svg, options = {}) => {
  if (!svg) return new Map();
  if (options.rebuild || !featureElementIndexCache.has(svg)) {
    return buildFeatureElementIndex(svg, options);
  }
  return featureElementIndexCache.get(svg) || new Map();
};

export const clearFeatureElementIndex = (svg) => {
  if (svg) featureElementIndexCache.delete(svg);
};

export const getFeatureElements = (svg, featureId, featureIndex = null) => {
  const normalizedId = String(featureId || '').trim();
  if (!svg || !normalizedId) return [];

  const indexed = featureIndex || featureElementIndexCache.get(svg);
  const indexedElements = indexed?.get?.(normalizedId);
  if (indexedElements?.length) return indexedElements;

  const byId = svg.getElementById?.(normalizedId) || svg.querySelector?.(`#${CSS.escape(normalizedId)}`);
  return byId ? [byId] : [];
};

export const createFeatureSvgActions = ({
  state,
  getFeatureColor,
  getEffectiveLegendCaption,
  onFeaturePopupOpened = null
}) => {
  const {
    results,
    selectedResultIndex,
    extractedFeatures,
    featuresBySvgId,
    featureColorOverrides,
    featureVisibilityOverrides,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    featurePopupSize,
    skipCaptureBaseConfig,
    adv
  } = state;
  const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
  const formatDuration = (ms) => `${ms.toFixed(1)}ms`;
  let delegatedFeatureHandlers = null;

  const getOrthogroupIds = (value) =>
    Array.from(new Set(
      String(value || '')
        .split(';')
        .map((entry) => entry.trim())
        .filter(Boolean)
    ));

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return normalized === 'on' || normalized === 'off' ? normalized : 'default';
  };

  const getPopupPosition = (eventLike, popupWidth = 720, popupHeight = 520) => {
    const margin = 12;
    const fallbackX = window.innerWidth / 2;
    const fallbackY = window.innerHeight / 2;
    const resolvedPopupWidth = Math.min(popupWidth, Math.max(0, window.innerWidth - (2 * margin)));
    const resolvedPopupHeight = Math.min(popupHeight, Math.max(0, window.innerHeight - (2 * margin)));
    const rawX = Number.isFinite(eventLike?.clientX) ? eventLike.clientX + 10 : fallbackX;
    const rawY = Number.isFinite(eventLike?.clientY) ? eventLike.clientY + 10 : fallbackY;
    const maxX = Math.max(margin, window.innerWidth - resolvedPopupWidth - margin);
    const maxY = Math.max(margin, window.innerHeight - resolvedPopupHeight - margin);
    return {
      x: Math.min(Math.max(rawX, margin), maxX),
      y: Math.min(Math.max(rawY, margin), maxY)
    };
  };

  const buildFeatureLocation = (feat) => {
    const startNumeric = Number(feat.start);
    const endNumeric = Number(feat.end);
    const startPos = Number.isFinite(startNumeric) ? startNumeric + 1 : feat.start;
    const endPos = Number.isFinite(endNumeric) ? endNumeric : feat.end;
    return `${startPos}..${endPos}${feat.strand ? ` (${feat.strand})` : ''}`;
  };

  const normalizeStringArray = (value) => {
    if (Array.isArray(value)) {
      return value
        .filter((item) => item !== null && item !== undefined)
        .map((item) => String(item));
    }
    if (value === null || value === undefined || value === '') return [];
    return [String(value)];
  };

  const normalizeQualifierRows = (qualifiers) => {
    if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) return [];
    return Object.entries(qualifiers)
      .map(([key, value]) => {
        const values = normalizeStringArray(value);
        return {
          key: String(key || ''),
          values,
          copyText: values.join('\n'),
          displayValue: values.join('\n')
        };
      })
      .filter((row) => row.key && row.values.length > 0)
      .sort((left, right) => left.key.localeCompare(right.key));
  };

  const displayProteinId = (feat) => String(
    feat?.sourceProteinId || feat?.source_protein_id || feat?.proteinId || feat?.protein_id || ''
  ).trim();

  const buildOrthogroupDetailRows = (feat) => {
    const member = feat?.orthogroupMember || feat?.orthogroup_member || null;
    const proteinId = displayProteinId(feat) || displayProteinId(member);
    const rows = [
      { key: 'orthogroup_id', label: 'Orthogroup ID', value: feat?.orthogroupId || feat?.orthogroup_id },
      { key: 'orthogroup_members', label: 'Members', value: feat?.orthogroupMemberCount || feat?.orthogroup_member_count },
      { key: 'orthogroup_coverage', label: 'Record coverage', value: feat?.orthogroupRecordCoverage || feat?.orthogroup_record_coverage },
      { key: 'protein_id', label: 'Protein ID', value: proteinId }
    ];
    return rows.filter((row) => String(row.value === null || row.value === undefined ? '' : row.value) !== '');
  };

  const buildDetailRows = ({ defaultLabel, feat, locationText }) => {
    const rows = [
      { key: 'label', label: 'Label', value: defaultLabel },
      { key: 'record_id', label: 'Record ID', value: feat.record_id },
      { key: 'type', label: 'Feature type', value: feat.type },
      { key: 'location', label: 'Location', value: locationText }
    ];
    rows.push(...buildOrthogroupDetailRows(feat));
    return rows
      .map((row) => ({ ...row, value: row.value === null || row.value === undefined ? '' : String(row.value) }))
      .filter((row) => row.value !== '');
  };

  const buildFeatureLookup = () => {
    if (featuresBySvgId?.value instanceof Map) return featuresBySvgId.value;
    const indexed = new Map();
    const features = Array.isArray(extractedFeatures.value) ? extractedFeatures.value : [];
    for (const feat of features) {
      const svgId = String(feat?.svg_id || '').trim();
      if (!svgId || indexed.has(svgId)) continue;
      indexed.set(svgId, feat);
    }
    return indexed;
  };

  const getFeatureTarget = (target, svg) => {
    if (!target || typeof target.closest !== 'function') return null;
    const featureEl = target.closest(FEATURE_SELECTOR);
    if (!featureEl || !svg.contains(featureEl)) return null;
    return featureEl;
  };

  const cleanupDelegatedFeatureHandlers = () => {
    if (!delegatedFeatureHandlers?.cleanup) return;
    delegatedFeatureHandlers.cleanup();
    delegatedFeatureHandlers = null;
  };

  const buildClickedFeaturePayload = (feat, featureElement = null) => {
    const defaultLabel = getFeatureCaption(feat);
    const existingOverride = featureColorOverrides[feat.id];
    const effectiveCaption = String(getEffectiveLegendCaption?.(feat) || existingOverride?.caption || defaultLabel || '').trim();
    const locationText = buildFeatureLocation(feat);
    const locationParts = Array.isArray(feat.location_parts) ? feat.location_parts : [];
    const qualifierRows = normalizeQualifierRows(feat.qualifiers);
    const sequenceWarnings = normalizeStringArray(feat.sequence_warnings);
    const nucleotideSequence = String(feat.nucleotide_sequence || '');
    const aminoAcidSequence = String(feat.amino_acid_sequence || '');

    const currentColor = resolveColorToHex(
      featureElement?.getAttribute('fill') || getFeatureColor(feat)
    );
    const currentStrokeColor = featureElement?.getAttribute('stroke') || '#000000';
    const currentStrokeWidth = parseFloat(featureElement?.getAttribute('stroke-width')) || 0.5;

    const visibilityMode = normalizeVisibilityMode(featureVisibilityOverrides[feat.svg_id]);

    return {
      id: feat.id,
      svg_id: feat.svg_id,
      label: defaultLabel,
      location: locationText,
      locationParts,
      color: currentColor,
      feat,
      activeTab: 'edit',
      recordId: String(feat.record_id || ''),
      recordIdx: Number.isInteger(Number(feat.record_idx)) ? Number(feat.record_idx) : null,
      featureType: String(feat.type || ''),
      start: Number.isFinite(Number(feat.start)) ? Number(feat.start) : null,
      end: Number.isFinite(Number(feat.end)) ? Number(feat.end) : null,
      strand: String(feat.strand || ''),
      qualifiers: feat.qualifiers && typeof feat.qualifiers === 'object' ? feat.qualifiers : {},
      qualifierRows,
      sequenceWarnings,
      nucleotideSequence,
      aminoAcidSequence,
      detailRows: buildDetailRows({ defaultLabel, feat, locationText }),
      legendName: effectiveCaption,
      appliedLegendName: effectiveCaption,
      strokeColor: currentStrokeColor,
      strokeWidth: currentStrokeWidth,
      originalStrokeColor: currentStrokeColor,
      originalStrokeWidth: currentStrokeWidth,
      labelKey: '',
      labelText: '',
      labelSourceText: '',
      labelVisibility: 'default',
      featureVisibility: visibilityMode,
      proteinId: feat.proteinId || '',
      sourceProteinId: feat.sourceProteinId || '',
      orthogroupId: feat.orthogroupId || '',
      orthogroupMemberCount: feat.orthogroupMemberCount || 0,
      orthogroupRecordCoverage: feat.orthogroupRecordCoverage || 0,
      orthogroupRepresentative: Boolean(feat.orthogroupRepresentative),
      orthogroupMember: feat.orthogroupMember || null,
      hasEditableLabel: false,
      labelUnavailableReason: 'No editable feature label for this feature in current diagram.'
    };
  };

  const openFeatureEditorForFeature = (feat, eventLike = null) => {
    if (!feat || !feat.svg_id) return null;
    if (!svgContainer.value) return null;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return null;

    const featureElements = getFeatureElements(svg, feat.svg_id);
    const featureElement = featureElements[0] || null;
    clickedFeature.value = buildClickedFeaturePayload(feat, featureElement);
    if (featurePopupSize) {
      featurePopupSize.width = 0;
      featurePopupSize.height = 0;
    }

    const popupPosition = getPopupPosition(eventLike, adv?.rich_feature_popup === false ? 440 : 720);
    clickedFeaturePos.x = popupPosition.x;
    clickedFeaturePos.y = popupPosition.y;
    if (typeof onFeaturePopupOpened === 'function') {
      onFeaturePopupOpened();
    }
    return clickedFeature.value;
  };

  const applyInstantPreview = (feat, color) => {
    const svgId = feat.svg_id;
    if (!svgId) {
      console.log('No svg_id for feature', feat);
      return;
    }

    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    try {
      const elements = getFeatureElements(svg, svgId);
      let updated = elements.length > 0;

      if (updated) {
        elements.forEach((el) => el.setAttribute('fill', color));
      }

      if (updated) {
        const serializer = new XMLSerializer();
        const newContent = serializer.serializeToString(svg);
        skipCaptureBaseConfig.value = true;
        const idx = selectedResultIndex.value;
        if (idx >= 0 && results.value.length > idx) {
          results.value[idx] = { ...results.value[idx], content: newContent };
        }
        console.log(`Instant preview: updated ${elements.length} element(s) for ${svgId} to ${color}`);
      } else {
        console.log(`Instant preview: element ${svgId} not found in SVG`);
      }
    } catch (e) {
      console.error('Instant preview error:', e);
    }
  };

  const applyVisibilityPreviewBySvgId = (svgId, modeRaw) => {
    const mode = normalizeVisibilityMode(modeRaw);
    if (!svgId || !svgContainer.value) return false;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return false;

    try {
      const elements = getFeatureElements(svg, svgId);
      if (!elements || elements.length === 0) {
        console.log(`Instant preview: element ${svgId} not found for visibility update`);
        return false;
      }
      elements.forEach((el) => {
        if (mode === 'off') {
          el.setAttribute('display', 'none');
        } else {
          el.removeAttribute('display');
        }
      });

      const serializer = new XMLSerializer();
      const newContent = serializer.serializeToString(svg);
      skipCaptureBaseConfig.value = true;
      const idx = selectedResultIndex.value;
      if (idx >= 0 && results.value.length > idx) {
        results.value[idx] = { ...results.value[idx], content: newContent };
      }
      return true;
    } catch (e) {
      console.error('Instant visibility preview error:', e);
      return false;
    }
  };

  const attachSvgFeatureHandlers = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    if (delegatedFeatureHandlers && delegatedFeatureHandlers.svg !== svg) {
      cleanupDelegatedFeatureHandlers();
    }

    const queryStartedAt = getNow();
    const pathsByIdMap = getFeatureElementIndex(svg, { rebuild: true, markCursor: true });
    const featurePathCount = Array.from(pathsByIdMap.values()).reduce((sum, elements) => sum + elements.length, 0);
    const queryDuration = getNow() - queryStartedAt;

    const indexStartedAt = getNow();
    const featureLookup = buildFeatureLookup();
    const featureIdsByOrthogroupId = new Map();
    featureLookup.forEach((feat, svgId) => {
      if (!svgId) return;
      getOrthogroupIds(feat?.orthogroupId).forEach((orthogroupId) => {
        if (!featureIdsByOrthogroupId.has(orthogroupId)) {
          featureIdsByOrthogroupId.set(orthogroupId, new Set());
        }
        featureIdsByOrthogroupId.get(orthogroupId).add(svgId);
      });
    });
    const comparisonElementsByOrthogroupId = new Map();
    svg.querySelectorAll('[data-orthogroup-id]').forEach((element) => {
      if (element.matches?.(FEATURE_SELECTOR)) return;
      getOrthogroupIds(element.getAttribute('data-orthogroup-id')).forEach((orthogroupId) => {
        if (!comparisonElementsByOrthogroupId.has(orthogroupId)) {
          comparisonElementsByOrthogroupId.set(orthogroupId, []);
        }
        comparisonElementsByOrthogroupId.get(orthogroupId).push(element);
      });
    });
    const indexDuration = getNow() - indexStartedAt;

    if (!delegatedFeatureHandlers) {
      const handlerState = {
        svg,
        pathsByIdMap,
        featureLookup,
        featureIdsByOrthogroupId,
        comparisonElementsByOrthogroupId,
        activeHoverSvgId: null,
        activeHoverKey: '',
        cleanup: null
      };

      const setHoverStyle = (element, highlight) => {
        if (!element?.style) return;
        if (highlight) {
          if (!element.hasAttribute('data-gbdraw-hover-opacity')) {
            element.setAttribute('data-gbdraw-hover-opacity', element.style.opacity || '');
            element.setAttribute('data-gbdraw-hover-filter', element.style.filter || '');
          }
          element.style.opacity = '0.7';
          element.style.filter = 'brightness(1.2)';
          return;
        }
        if (element.hasAttribute('data-gbdraw-hover-opacity')) {
          element.style.opacity = element.getAttribute('data-gbdraw-hover-opacity') || '';
          element.style.filter = element.getAttribute('data-gbdraw-hover-filter') || '';
          element.removeAttribute('data-gbdraw-hover-opacity');
          element.removeAttribute('data-gbdraw-hover-filter');
        }
      };

      const setFeatureHover = (svgId, highlight) => {
        (handlerState.pathsByIdMap.get(svgId) || []).forEach((element) => {
          setHoverStyle(element, highlight);
        });
      };

      const getFeatureHoverKey = (svgId) => {
        const feat = handlerState.featureLookup.get(svgId);
        const orthogroupId = String(feat?.orthogroupId || '').trim();
        return orthogroupId ? `orthogroup:${orthogroupId}` : `feature:${svgId}`;
      };

      const setOrthogroupHover = (orthogroupId, highlight) => {
        const id = String(orthogroupId || '').trim();
        if (!id) return;
        (handlerState.featureIdsByOrthogroupId.get(id) || new Set()).forEach((featureId) => {
          setFeatureHover(featureId, highlight);
        });
        (handlerState.comparisonElementsByOrthogroupId.get(id) || []).forEach((element) => {
          setHoverStyle(element, highlight);
        });
      };

      const setHoverHighlight = (svgId, highlight) => {
        const feat = handlerState.featureLookup.get(svgId);
        const orthogroupId = String(feat?.orthogroupId || '').trim();
        if (orthogroupId) {
          setOrthogroupHover(orthogroupId, highlight);
          return;
        }
        setFeatureHover(svgId, highlight);
      };

      const handleMouseOver = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        const svgId = getFeatureIdentity(featureEl);
        if (!svgId) return;
        const hoverKey = getFeatureHoverKey(svgId);
        if (handlerState.activeHoverKey === hoverKey) return;
        if (handlerState.activeHoverSvgId) {
          setHoverHighlight(handlerState.activeHoverSvgId, false);
        }
        handlerState.activeHoverSvgId = svgId;
        handlerState.activeHoverKey = hoverKey;
        setHoverHighlight(svgId, true);
      };

      const handleMouseOut = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        const svgId = getFeatureIdentity(featureEl);
        if (!svgId || handlerState.activeHoverSvgId !== svgId) return;
        const relatedFeature = getFeatureTarget(e.relatedTarget, svg);
        if (relatedFeature && getFeatureHoverKey(getFeatureIdentity(relatedFeature)) === handlerState.activeHoverKey) return;
        setHoverHighlight(svgId, false);
        handlerState.activeHoverSvgId = null;
        handlerState.activeHoverKey = '';
      };

      const handleClick = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        const svgId = getFeatureIdentity(featureEl);
        if (!svgId) return;
        e.stopPropagation();
        const feat = handlerState.featureLookup.get(svgId);
        if (feat) {
          openFeatureEditorForFeature(feat, e);
        } else {
          console.log(`No feature found for svg_id: ${svgId}`);
        }
      };

      svg.addEventListener('mouseover', handleMouseOver);
      svg.addEventListener('mouseout', handleMouseOut);
      svg.addEventListener('click', handleClick);
      handlerState.cleanup = () => {
        svg.removeEventListener('mouseover', handleMouseOver);
        svg.removeEventListener('mouseout', handleMouseOut);
        svg.removeEventListener('click', handleClick);
      };
      delegatedFeatureHandlers = handlerState;
    } else {
      delegatedFeatureHandlers.pathsByIdMap = pathsByIdMap;
      delegatedFeatureHandlers.featureLookup = featureLookup;
      delegatedFeatureHandlers.featureIdsByOrthogroupId = featureIdsByOrthogroupId;
      delegatedFeatureHandlers.comparisonElementsByOrthogroupId = comparisonElementsByOrthogroupId;
    }

    console.groupCollapsed('post-gbdraw timing');
    console.info(`feature handler index querySelectorAll: ${formatDuration(queryDuration)}`);
    console.info(`feature handler index/delegation setup: ${formatDuration(indexDuration)}`);
    console.groupEnd();
    console.log(
      `Delegated feature handlers for ${featurePathCount} feature paths (${pathsByIdMap.size} unique features)`
    );
  };

  return {
    applyInstantPreview,
    applyVisibilityPreviewBySvgId,
    attachSvgFeatureHandlers,
    getFeatureElements,
    openFeatureEditorForFeature
  };
};
