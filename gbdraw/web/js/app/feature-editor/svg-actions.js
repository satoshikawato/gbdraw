import { resolveColorToHex } from '../color-utils.js';
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

    const featurePaths = svg.querySelectorAll('path[id^="f"], polygon[id^="f"], rect[id^="f"]');

    const pathsByIdMap = {};
    featurePaths.forEach((path) => {
      const id = path.getAttribute('id');
      if (!id) return;
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
      if (!svgId) return;
      path.style.cursor = 'pointer';

      path.addEventListener('mouseenter', () => highlightFeature(svgId, true));
      path.addEventListener('mouseleave', () => highlightFeature(svgId, false));

      path.addEventListener('click', (e) => {
        e.stopPropagation();
        const feat = extractedFeatures.value.find((f) => f.svg_id === svgId);
        if (feat) {
          const currentColor = resolveColorToHex(path.getAttribute('fill') || getFeatureColor(feat));
          const defaultLabel = getFeatureCaption(feat);

          const existingOverride = featureColorOverrides[feat.id];
          const existingCaption = existingOverride?.caption || '';

          const currentStrokeColor = path.getAttribute('stroke') || '#000000';
          const currentStrokeWidth = parseFloat(path.getAttribute('stroke-width')) || 0.5;
          const startPos = Number.isFinite(feat.start) ? feat.start + 1 : feat.start;
          const endPos = Number.isFinite(feat.end) ? feat.end : feat.end;
          const location = `${startPos}..${endPos}${feat.strand ? ` (${feat.strand})` : ''}`;

          clickedFeature.value = {
            id: feat.id,
            svg_id: svgId,
            label: defaultLabel,
            location,
            color: currentColor,
            feat: feat,
            legendName: existingCaption,
            strokeColor: currentStrokeColor,
            strokeWidth: currentStrokeWidth,
            originalStrokeColor: currentStrokeColor,
            originalStrokeWidth: currentStrokeWidth
          };
          const popupWidth = 360;
          const popupHeight = 260;
          const margin = 12;
          const maxX = Math.max(margin, window.innerWidth - popupWidth - margin);
          const maxY = Math.max(margin, window.innerHeight - popupHeight - margin);
          clickedFeaturePos.x = Math.min(Math.max(e.clientX + 10, margin), maxX);
          clickedFeaturePos.y = Math.min(Math.max(e.clientY + 10, margin), maxY);
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
