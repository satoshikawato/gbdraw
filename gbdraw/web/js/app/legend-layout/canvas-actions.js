import { getElementsBounds } from './transform-utils.js';

export const createLegendCanvasActions = ({ state }) => {
  const {
    svgContainer,
    canvasPadding,
    originalSvgStroke,
    mode,
    linearBaseConfig,
    circularBaseConfig,
    diagramElementOriginalTransforms,
    diagramElementBaseTransforms,
    form
  } = state;

  const applyCanvasPadding = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    let viewBox = svg.getAttribute('viewBox');
    if (!viewBox) {
      const width = parseFloat(svg.getAttribute('width')) || 800;
      const height = parseFloat(svg.getAttribute('height')) || 600;
      viewBox = `0 0 ${width} ${height}`;
    }

    const parts = viewBox.split(/\s+/).map(parseFloat);
    if (parts.length !== 4) return;

    if (!svg.dataset.originalViewBox) {
      svg.dataset.originalViewBox = viewBox;
    }

    const [origX, origY, origW, origH] = svg.dataset.originalViewBox.split(/\s+/).map(parseFloat);

    const newX = origX - canvasPadding.left;
    const newY = origY - canvasPadding.top;
    const newW = origW + canvasPadding.left + canvasPadding.right;
    const newH = origH + canvasPadding.top + canvasPadding.bottom;

    svg.setAttribute('viewBox', `${newX} ${newY} ${newW} ${newH}`);

    const origWidth = parseFloat(svg.dataset.originalWidth || svg.getAttribute('width')) || origW;
    const origHeight = parseFloat(svg.dataset.originalHeight || svg.getAttribute('height')) || origH;
    if (!svg.dataset.originalWidth) {
      svg.dataset.originalWidth = origWidth;
      svg.dataset.originalHeight = origHeight;
    }

    const scaleX = newW / origW;
    const scaleY = newH / origH;
    svg.setAttribute('width', origWidth * scaleX);
    svg.setAttribute('height', origHeight * scaleY);
  };

  const resetCanvasPadding = () => {
    canvasPadding.top = 0;
    canvasPadding.right = 0;
    canvasPadding.bottom = 0;
    canvasPadding.left = 0;

    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    if (svg.dataset.originalViewBox) {
      svg.setAttribute('viewBox', svg.dataset.originalViewBox);
    }
    if (svg.dataset.originalWidth) {
      svg.setAttribute('width', svg.dataset.originalWidth);
      svg.setAttribute('height', svg.dataset.originalHeight);
    }
  };

  const captureOriginalStroke = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const firstFeaturePath = svg.querySelector('path[id^="f"]');
    if (firstFeaturePath) {
      const strokeColor = firstFeaturePath.getAttribute('stroke');
      const strokeWidthAttr = firstFeaturePath.getAttribute('stroke-width');
      const strokeWidth = strokeWidthAttr !== null ? parseFloat(strokeWidthAttr) : null;

      originalSvgStroke.value = { color: strokeColor, width: strokeWidth };
      console.log(`Captured original stroke: color=${strokeColor}, width=${strokeWidth}`);
    }
  };

  const captureBaseConfig = () => {
    if (!svgContainer.value) return;
    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const viewBox = svg.getAttribute('viewBox');
    if (!viewBox) return;
    const parts = viewBox.split(/\s+/).map(parseFloat);
    if (parts.length !== 4) return;
    const [vbX, vbY, vbW, vbH] = parts;

    const legendGroup = svg.getElementById('legend');
    let legendWidth = 120;
    let legendHeight = 150;
    if (legendGroup) {
      const bbox = legendGroup.getBBox();
      legendWidth = bbox.width || legendWidth;
      legendHeight = bbox.height || legendHeight;
    }

    const isLinear = mode.value === 'linear';

    if (isLinear) {
      const verticalVb = svg.getAttribute('data-vertical-viewbox');
      const horizontalVb = svg.getAttribute('data-horizontal-viewbox');

      if (verticalVb && horizontalVb) {
        const vParts = verticalVb.split(/\s+/).map(parseFloat);
        const hParts = horizontalVb.split(/\s+/).map(parseFloat);
        linearBaseConfig.value.verticalViewBox = { x: vParts[0], y: vParts[1], w: vParts[2], h: vParts[3] };
        linearBaseConfig.value.horizontalViewBox = { x: hParts[0], y: hParts[1], w: hParts[2], h: hParts[3] };
      } else {
        let vLegendW = legendWidth;
        let vLegendH = legendHeight;
        let hLegendW = legendWidth;
        let hLegendH = legendHeight;

        if (legendGroup) {
          const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
          const verticalLegend = legendGroup.querySelector('#legend_vertical');
          if (horizontalLegend && verticalLegend) {
            const hDisplay = horizontalLegend.getAttribute('display');
            const vDisplay = verticalLegend.getAttribute('display');
            horizontalLegend.removeAttribute('display');
            verticalLegend.removeAttribute('display');

            const hBbox = horizontalLegend.getBBox();
            const vBbox = verticalLegend.getBBox();
            hLegendW = hBbox.width || legendWidth;
            hLegendH = hBbox.height || legendHeight;
            vLegendW = vBbox.width || legendWidth;
            vLegendH = vBbox.height || legendHeight;

            if (hDisplay) horizontalLegend.setAttribute('display', hDisplay);
            else horizontalLegend.removeAttribute('display');
            if (vDisplay) verticalLegend.setAttribute('display', vDisplay);
            else verticalLegend.removeAttribute('display');
          }
        }

        const genPos = form.legend;
        const isGeneratedVertical = genPos === 'left' || genPos === 'right';

        if (isGeneratedVertical) {
          linearBaseConfig.value.verticalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
          linearBaseConfig.value.horizontalViewBox = {
            x: vbX,
            y: vbY,
            w: vbW - vLegendW,
            h: vbH + hLegendH
          };
        } else {
          linearBaseConfig.value.horizontalViewBox = { x: vbX, y: vbY, w: vbW, h: vbH };
          linearBaseConfig.value.verticalViewBox = {
            x: vbX,
            y: vbY,
            w: vbW + vLegendW,
            h: vbH - hLegendH
          };
        }
      }

      if (legendGroup) {
        const horizontalLegend = legendGroup.querySelector('#legend_horizontal');
        const verticalLegend = legendGroup.querySelector('#legend_vertical');

        if (horizontalLegend && verticalLegend) {
          const hDisplay = horizontalLegend.getAttribute('display');
          const vDisplay = verticalLegend.getAttribute('display');
          horizontalLegend.removeAttribute('display');
          verticalLegend.removeAttribute('display');

          const hBbox = horizontalLegend.getBBox();
          const vBbox = verticalLegend.getBBox();

          linearBaseConfig.value.horizontalLegendWidth = hBbox.width || legendWidth;
          linearBaseConfig.value.horizontalLegendHeight = hBbox.height || legendHeight;
          linearBaseConfig.value.verticalLegendWidth = vBbox.width || legendWidth;
          linearBaseConfig.value.verticalLegendHeight = vBbox.height || legendHeight;

          if (hDisplay) horizontalLegend.setAttribute('display', hDisplay);
          if (vDisplay) verticalLegend.setAttribute('display', vDisplay);
        } else {
          linearBaseConfig.value.verticalLegendWidth = legendWidth;
          linearBaseConfig.value.verticalLegendHeight = legendHeight;
          linearBaseConfig.value.horizontalLegendWidth = legendWidth;
          linearBaseConfig.value.horizontalLegendHeight = legendHeight;
        }
      }

      linearBaseConfig.value.generatedPosition = form.legend;
      linearBaseConfig.value.diagramBaseTransforms = new Map(diagramElementOriginalTransforms.value);

      const baseHorizontal = linearBaseConfig.value.horizontalViewBox;
      if (baseHorizontal && baseHorizontal.w > 0) {
        svg.setAttribute('data-horizontal-wrap-width', `${baseHorizontal.w}`);
      }
    } else {
      const legendPos = form.legend;
      let baseVbW = vbW;
      let baseVbH = vbH;
      const diagramBounds = getElementsBounds(diagramElements.value);

      if (diagramBounds) {
        const marginLeft = diagramBounds.x - vbX;
        const marginRight = vbX + vbW - (diagramBounds.x + diagramBounds.width);
        const marginTop = diagramBounds.y - vbY;
        const marginBottom = vbY + vbH - (diagramBounds.y + diagramBounds.height);
        const symmetricMarginX = Math.max(0, Math.min(marginLeft, marginRight));
        const symmetricMarginY = Math.max(0, Math.min(marginTop, marginBottom));
        const derivedBaseWidth = diagramBounds.width + symmetricMarginX * 2;
        const derivedBaseHeight = diagramBounds.height + symmetricMarginY * 2;
        if (Number.isFinite(derivedBaseWidth) && derivedBaseWidth > 0) {
          baseVbW = derivedBaseWidth;
        }
        if (Number.isFinite(derivedBaseHeight) && derivedBaseHeight > 0) {
          baseVbH = derivedBaseHeight;
        }
      }

      const diagramCenterX = vbX + baseVbW / 2;
      const diagramCenterY = vbY + baseVbH / 2;

      circularBaseConfig.value = {
        viewBoxWidth: baseVbW,
        viewBoxHeight: baseVbH,
        generatedViewBoxWidth: vbW,
        generatedViewBoxHeight: vbH,
        diagramCenterX: diagramCenterX,
        diagramCenterY: diagramCenterY,
        legendWidth: legendWidth,
        legendHeight: legendHeight,
        generatedPosition: legendPos
      };
    }

    diagramElementBaseTransforms.value = new Map(diagramElementOriginalTransforms.value);
  };

  return {
    applyCanvasPadding,
    captureBaseConfig,
    captureOriginalStroke,
    resetCanvasPadding
  };
};
