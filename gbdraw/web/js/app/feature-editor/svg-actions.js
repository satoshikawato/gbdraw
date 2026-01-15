import { getFeatureCaption } from '../feature-utils.js';

export const createFeatureSvgActions = ({ state, getFeatureColor }) => {
  const {
    results,
    selectedResultIndex,
    extractedFeatures,
    featureColorOverrides,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    skipCaptureBaseConfig
  } = state;

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

  const attachSvgFeatureHandlers = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const featurePaths = svg.querySelectorAll('path[id^="f"], polygon[id^="f"]');

    const pathsByIdMap = {};
    featurePaths.forEach((path) => {
      const id = path.getAttribute('id');
      if (!pathsByIdMap[id]) pathsByIdMap[id] = [];
      pathsByIdMap[id].push(path);
    });

    const highlightFeature = (svgId, highlight) => {
      const paths = pathsByIdMap[svgId] || [];
      paths.forEach((p) => {
        p.style.opacity = highlight ? '0.7' : '1';
        p.style.filter = highlight ? 'brightness(1.2)' : 'none';
      });
    };

    featurePaths.forEach((path) => {
      const svgId = path.getAttribute('id');
      path.style.cursor = 'pointer';

      path.addEventListener('mouseenter', () => highlightFeature(svgId, true));
      path.addEventListener('mouseleave', () => highlightFeature(svgId, false));

      path.addEventListener('click', (e) => {
        e.stopPropagation();
        const feat = extractedFeatures.value.find((f) => f.svg_id === svgId);
        if (feat) {
          const currentColor = path.getAttribute('fill') || getFeatureColor(feat);
          const defaultLabel = getFeatureCaption(feat);

          const existingOverride = featureColorOverrides[feat.id];
          const existingCaption = existingOverride?.caption || '';

          const currentStrokeColor = path.getAttribute('stroke') || '#000000';
          const currentStrokeWidth = parseFloat(path.getAttribute('stroke-width')) || 0.5;

          clickedFeature.value = {
            id: feat.id,
            svg_id: svgId,
            label: defaultLabel,
            color: currentColor,
            feat: feat,
            legendName: existingCaption,
            strokeColor: currentStrokeColor,
            strokeWidth: currentStrokeWidth,
            originalStrokeColor: currentStrokeColor,
            originalStrokeWidth: currentStrokeWidth
          };
          clickedFeaturePos.x = Math.min(e.clientX + 10, window.innerWidth - 280);
          clickedFeaturePos.y = Math.min(e.clientY + 10, window.innerHeight - 100);
        } else {
          console.log(`No feature found for svg_id: ${svgId}`);
        }
      });
    });

    console.log(
      `Attached handlers to ${featurePaths.length} feature paths (${Object.keys(pathsByIdMap).length} unique features)`
    );
  };

  return {
    applyInstantPreview,
    attachSvgFeatureHandlers
  };
};
