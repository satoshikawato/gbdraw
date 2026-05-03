import { resolveColorToHex } from '../color-utils.js';
import { getFeatureCaption } from '../feature-utils.js';

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
    skipCaptureBaseConfig
  } = state;
  const FEATURE_SELECTOR = 'path[id^="f"], polygon[id^="f"], rect[id^="f"]';
  const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
  const formatDuration = (ms) => `${ms.toFixed(1)}ms`;
  let delegatedFeatureHandlers = null;

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return normalized === 'on' || normalized === 'off' ? normalized : 'default';
  };

  const getPopupPosition = (eventLike, popupWidth = 360, popupHeight = 360) => {
    const margin = 12;
    const fallbackX = window.innerWidth / 2;
    const fallbackY = window.innerHeight / 2;
    const rawX = Number.isFinite(eventLike?.clientX) ? eventLike.clientX + 10 : fallbackX;
    const rawY = Number.isFinite(eventLike?.clientY) ? eventLike.clientY + 10 : fallbackY;
    const maxX = Math.max(margin, window.innerWidth - popupWidth - margin);
    const maxY = Math.max(margin, window.innerHeight - popupHeight - margin);
    return {
      x: Math.min(Math.max(rawX, margin), maxX),
      y: Math.min(Math.max(rawY, margin), maxY)
    };
  };

  const buildFeatureLocation = (feat) => {
    const startPos = Number.isFinite(feat.start) ? feat.start + 1 : feat.start;
    const endPos = Number.isFinite(feat.end) ? feat.end : feat.end;
    return `${startPos}..${endPos}${feat.strand ? ` (${feat.strand})` : ''}`;
  };

  const getFeatureElements = (svg, svgId) => {
    if (!svg || !svgId) return [];
    return Array.from(svg.querySelectorAll(`#${CSS.escape(svgId)}`));
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
      location: buildFeatureLocation(feat),
      color: currentColor,
      feat,
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
      orthogroupRepresentative: Boolean(feat.orthogroupRepresentative),
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

    const popupPosition = getPopupPosition(eventLike);
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
      const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
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
      const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
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
    const featurePaths = Array.from(svg.querySelectorAll(FEATURE_SELECTOR));
    const queryDuration = getNow() - queryStartedAt;

    const indexStartedAt = getNow();
    const pathsByIdMap = new Map();
    featurePaths.forEach((path) => {
      const id = path.getAttribute('id');
      if (!id) return;
      if (!pathsByIdMap.has(id)) pathsByIdMap.set(id, []);
      pathsByIdMap.get(id).push(path);
      path.style.cursor = 'pointer';
    });
    const featureLookup = buildFeatureLookup();
    const indexDuration = getNow() - indexStartedAt;

    if (!delegatedFeatureHandlers) {
      const handlerState = {
        svg,
        pathsByIdMap,
        featureLookup,
        activeHoverSvgId: null,
        cleanup: null
      };

      const highlightFeature = (svgId, highlight) => {
        const paths = handlerState.pathsByIdMap.get(svgId) || [];
        paths.forEach((p) => {
          p.style.opacity = highlight ? '0.7' : '1';
          p.style.filter = highlight ? 'brightness(1.2)' : 'none';
        });
      };

      const handleMouseOver = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        const svgId = featureEl?.getAttribute('id');
        if (!svgId || handlerState.activeHoverSvgId === svgId) return;
        if (handlerState.activeHoverSvgId) {
          highlightFeature(handlerState.activeHoverSvgId, false);
        }
        handlerState.activeHoverSvgId = svgId;
        highlightFeature(svgId, true);
      };

      const handleMouseOut = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        const svgId = featureEl?.getAttribute('id');
        if (!svgId || handlerState.activeHoverSvgId !== svgId) return;
        const relatedFeature = getFeatureTarget(e.relatedTarget, svg);
        if (relatedFeature?.getAttribute('id') === svgId) return;
        highlightFeature(svgId, false);
        handlerState.activeHoverSvgId = null;
      };

      const handleClick = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        const svgId = featureEl?.getAttribute('id');
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
    }

    console.groupCollapsed('post-gbdraw timing');
    console.info(`feature handler querySelectorAll: ${formatDuration(queryDuration)}`);
    console.info(`feature handler index/delegation setup: ${formatDuration(indexDuration)}`);
    console.groupEnd();
    console.log(
      `Delegated feature handlers for ${featurePaths.length} feature paths (${pathsByIdMap.size} unique features)`
    );
  };

  return {
    applyInstantPreview,
    applyVisibilityPreviewBySvgId,
    attachSvgFeatureHandlers,
    openFeatureEditorForFeature
  };
};
