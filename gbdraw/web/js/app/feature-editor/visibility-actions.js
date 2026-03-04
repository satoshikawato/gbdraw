export const createFeatureVisibilityActions = ({ state, featureSvgActions }) => {
  const {
    clickedFeature,
    featureVisibilityOverrides,
    labelReflowForceRequestSeq,
    labelReflowForceRequestReason
  } = state;

  const { applyVisibilityPreviewBySvgId } = featureSvgActions;

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return normalized === 'on' || normalized === 'off' ? normalized : 'default';
  };

  const queueFeatureVisibilityReflow = (reason = 'feature-visibility-apply') => {
    labelReflowForceRequestReason.value = String(reason || 'feature-visibility-apply');
    labelReflowForceRequestSeq.value += 1;
  };

  const setFeatureVisibility = (feat, modeRaw, options = {}) => {
    const svgId = String(feat?.svg_id || '').trim();
    if (!svgId) return false;

    const triggerReflow = options.triggerReflow !== false;
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = normalizeVisibilityMode(featureVisibilityOverrides[svgId]);

    if (nextMode === 'default') {
      if (Object.prototype.hasOwnProperty.call(featureVisibilityOverrides, svgId)) {
        delete featureVisibilityOverrides[svgId];
      }
    } else {
      featureVisibilityOverrides[svgId] = nextMode;
    }

    if (clickedFeature.value?.svg_id === svgId) {
      clickedFeature.value.featureVisibility = nextMode;
    }

    applyVisibilityPreviewBySvgId(svgId, nextMode);

    if (triggerReflow && previousMode !== nextMode) {
      queueFeatureVisibilityReflow();
    }

    return previousMode !== nextMode;
  };

  const updateClickedFeatureVisibility = (modeRaw) => {
    if (!clickedFeature.value?.feat) return false;
    return setFeatureVisibility(clickedFeature.value.feat, modeRaw, { triggerReflow: true });
  };

  const getFeatureVisibility = (feat) => {
    const svgId = String(feat?.svg_id || '').trim();
    if (!svgId) return 'default';
    return normalizeVisibilityMode(featureVisibilityOverrides[svgId]);
  };

  return {
    getFeatureVisibility,
    setFeatureVisibility,
    updateClickedFeatureVisibility
  };
};
