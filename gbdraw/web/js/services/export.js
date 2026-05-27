import { state } from '../state.js';
import { setDpiInPng } from '../utils/png.js';

const getDownloadName = (extension) => {
  const baseName =
    state.results.value?.[state.selectedResultIndex.value]?.name ||
    (extension ? `gbdraw.${extension}` : 'gbdraw.svg');
  if (!extension) return baseName;
  const normalized = baseName.replace(/\.svg$/i, `.${extension}`);
  if (normalized === baseName && !baseName.toLowerCase().endsWith(`.${extension}`)) {
    return `${baseName}.${extension}`;
  }
  return normalized;
};

const cloneCurrentSvg = () => {
  const liveSvg = state.svgContainer.value?.querySelector('svg');
  if (liveSvg) {
    const clone = liveSvg.cloneNode(true);
    if (!clone.getAttribute('xmlns')) {
      clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    }
    if (!clone.getAttribute('xmlns:xlink')) {
      clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
    }
    return clone;
  }
  return getSvgFromString(state.svgContent.value);
};

const getCurrentSvgString = ({ interactive = false } = {}) => {
  const clone = cloneCurrentSvg();
  if (!clone) return state.svgContent.value;
  if (interactive) {
    enrichSvgWithStandaloneFeaturePopup(clone);
  }
  return new XMLSerializer().serializeToString(clone);
};

const getSvgFromString = (svgString) => {
  if (!svgString) return null;
  const svgEl = document.createElement('div');
  svgEl.innerHTML = svgString;
  return svgEl.querySelector('svg');
};

const getSvgDimensions = (svg) => {
  if (!svg) return null;
  let w = parseFloat(svg.getAttribute('width'));
  let h = parseFloat(svg.getAttribute('height'));
  if (!w || !h) {
    const viewBox = svg.getAttribute('viewBox');
    if (!viewBox) return null;
    const parts = viewBox.trim().split(/[\s,]+/);
    if (parts.length < 4) return null;
    w = parseFloat(parts[2]);
    h = parseFloat(parts[3]);
  }
  if (!w || !h) return null;
  return { width: w, height: h };
};

const SVG_NS = 'http://www.w3.org/2000/svg';
const FEATURE_SELECTOR = 'path[id^="f"], polygon[id^="f"], rect[id^="f"]';
const INTERACTIVE_METADATA_ID = 'gbdraw-interactive-feature-metadata';
const INTERACTIVE_STYLE_ID = 'gbdraw-interactive-feature-style';
const INTERACTIVE_SCRIPT_ID = 'gbdraw-interactive-feature-script';

const STANDALONE_INTERACTIVE_STYLE = `
.gbdraw-interactive-feature {
  cursor: pointer;
}
.gbdraw-feature-popup {
  overflow: visible;
  pointer-events: auto;
}
.gbdraw-feature-popup * {
  box-sizing: border-box;
}
.gfi {
  position: relative;
  width: 100%;
  height: 100%;
  display: flex;
  flex-direction: column;
  overflow: hidden;
  border: 1px solid #cbd5e1;
  border-radius: 10px;
  background: #ffffff;
  color: #334155;
  box-shadow: 0 18px 45px rgba(15, 23, 42, 0.22);
  font-family: Arial, Helvetica, sans-serif;
  font-size: var(--gfi-font-size, 13px);
}
.gfi-header {
  display: grid;
  grid-template-columns: minmax(0, 1fr) 28px;
  gap: 8px;
  align-items: start;
  padding: 10px 12px;
  border-bottom: 1px solid #e2e8f0;
}
.gfi-title {
  min-width: 0;
  font-weight: 700;
  color: #1e293b;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.gfi-subtitle {
  margin-top: 2px;
  color: #64748b;
  font-size: var(--gfi-subtitle-font-size, 11px);
  overflow-wrap: anywhere;
}
.gfi-close,
.gfi-copy,
.gfi-tab {
  border: 0;
  border-radius: 7px;
  background: transparent;
  color: #475569;
  cursor: pointer;
  font: inherit;
}
.gfi-close {
  width: 28px;
  height: 28px;
  font-size: var(--gfi-close-font-size, 20px);
  line-height: 1;
}
.gfi-close:hover,
.gfi-copy:hover,
.gfi-tab:hover {
  background: #eff6ff;
  color: #2563eb;
}
.gfi-tabs {
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 4px;
  margin: 8px 10px 0;
  padding: 4px;
  border-radius: 9px;
  background: #f1f5f9;
}
.gfi-tab {
  min-width: 0;
  padding: 6px 8px;
  font-size: var(--gfi-tab-font-size, 11px);
  font-weight: 700;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.gfi-tab.is-active {
  background: #ffffff;
  color: #0f172a;
  box-shadow: 0 1px 2px rgba(15, 23, 42, 0.12);
}
.gfi-content {
  flex: 1;
  min-height: 0;
  overflow: auto;
  padding: 10px 18px 18px 10px;
}
.gfi-row {
  display: grid;
  grid-template-columns: 86px minmax(0, 1fr) 50px;
  gap: 8px;
  align-items: start;
  padding: 7px 0;
  border-bottom: 1px solid #f1f5f9;
}
.gfi-key {
  color: #64748b;
  font-size: var(--gfi-key-font-size, 10px);
  font-weight: 700;
  text-transform: uppercase;
}
.gfi-value {
  min-width: 0;
  color: #1e293b;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: var(--gfi-value-font-size, 11px);
  overflow-wrap: anywhere;
  white-space: pre-wrap;
}
.gfi-copy {
  width: 48px;
  padding: 5px 0;
  background: #f8fafc;
  font-size: var(--gfi-copy-font-size, 10px);
  font-weight: 700;
}
.gfi-block {
  margin-bottom: 10px;
  border: 1px solid #e2e8f0;
  border-radius: 8px;
  overflow: hidden;
}
.gfi-block-title {
  display: flex;
  gap: 8px;
  align-items: center;
  padding: 7px 9px;
  background: #f8fafc;
  color: #475569;
  font-size: var(--gfi-block-title-font-size, 10px);
  font-weight: 700;
  text-transform: uppercase;
}
.gfi-block-title .gfi-copy {
  margin-left: auto;
}
.gfi-pre {
  max-height: 120px;
  margin: 0;
  padding: 8px 9px;
  overflow: auto;
  color: #1e293b;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: var(--gfi-pre-font-size, 11px);
  line-height: 1.45;
  overflow-wrap: anywhere;
  white-space: pre-wrap;
}
.gfi-empty,
.gfi-warning {
  padding: 10px;
  border-radius: 8px;
  background: #f8fafc;
  color: #64748b;
}
.gfi-warning {
  margin-bottom: 8px;
  border: 1px solid #fcd34d;
  background: #fffbeb;
  color: #92400e;
}
.gfi-resize-handle {
  position: absolute;
  right: 6px;
  bottom: 6px;
  width: 18px;
  height: 18px;
  border: 0;
  border-radius: 5px;
  background: rgba(248, 250, 252, 0.9);
  color: #64748b;
  cursor: nwse-resize;
}
.gfi-resize-handle:before {
  content: "";
  position: absolute;
  right: 4px;
  bottom: 4px;
  width: 9px;
  height: 9px;
  border-right: 2px solid currentColor;
  border-bottom: 2px solid currentColor;
  box-shadow: 4px 4px 0 -2px currentColor;
}
.gfi-resize-handle:hover {
  background: #eff6ff;
  color: #2563eb;
}
`;

const STANDALONE_INTERACTIVE_SCRIPT = `
(function () {
  'use strict';

  var FEATURE_SELECTOR = 'path[id^="f"], polygon[id^="f"], rect[id^="f"]';
  var SVG_NS = 'http://www.w3.org/2000/svg';
  var XHTML_NS = 'http://www.w3.org/1999/xhtml';
  var svg = document.currentScript && document.currentScript.ownerSVGElement
    ? document.currentScript.ownerSVGElement
    : document.documentElement;
  var metadata = svg.querySelector('#gbdraw-interactive-feature-metadata');
  var payload = null;
  var popup = null;
  var activePopupResize = null;
  var baseDevicePixelRatio = Math.max(1, Number(window.devicePixelRatio) || 1);

  try {
    payload = JSON.parse(metadata ? metadata.textContent || '{}' : '{}');
  } catch (error) {
    payload = {};
  }

  var features = Array.isArray(payload.features) ? payload.features : [];
  if (!features.length) return;

  var featuresById = new Map();
  features.forEach(function (feature) {
    var svgId = String(feature && feature.svg_id || '').trim();
    if (svgId && !featuresById.has(svgId)) {
      featuresById.set(svgId, feature);
    }
  });

  function escapeHtml(value) {
    return String(value == null ? '' : value).replace(/[&<>"']/g, function (character) {
      return {
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#39;'
      }[character];
    });
  }

  function normalizeArray(value) {
    if (Array.isArray(value)) {
      return value.filter(function (item) {
        return item !== null && item !== undefined;
      }).map(function (item) {
        return String(item);
      });
    }
    if (value === null || value === undefined || value === '') return [];
    return [String(value)];
  }

  function locationText(feature) {
    if (feature.location) return String(feature.location);
    var start = Number(feature.start);
    var end = Number(feature.end);
    var text = (Number.isFinite(start) ? start + 1 : '') + '..' + (Number.isFinite(end) ? end : '');
    if (feature.strand) text += ' (' + feature.strand + ')';
    return text;
  }

  function detailRows(feature) {
    return [
      ['Label', feature.label || ''],
      ['SVG ID', feature.svg_id || ''],
      ['Record ID', feature.record_id || ''],
      ['Record index', feature.record_idx === null || feature.record_idx === undefined ? '' : feature.record_idx],
      ['Type', feature.type || ''],
      ['Location', locationText(feature)],
      ['Strand', feature.strand || '']
    ].filter(function (row) {
      return String(row[1] == null ? '' : row[1]) !== '';
    });
  }

  function qualifierRows(feature) {
    var qualifiers = feature.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
      ? feature.qualifiers
      : {};
    return Object.keys(qualifiers).sort().map(function (key) {
      var values = normalizeArray(qualifiers[key]);
      return {
        key: key,
        values: values,
        text: values.join('\\n')
      };
    }).filter(function (row) {
      return row.key && row.values.length;
    });
  }

  function copyButton(value) {
    var index = copyValues.push(String(value == null ? '' : value)) - 1;
    return '<button type="button" class="gfi-copy" data-copy-index="' + index + '">Copy</button>';
  }

  var copyValues = [];

  function renderRows(rows) {
    if (!rows.length) {
      return '<div class="gfi-empty">No details available.</div>';
    }
    return rows.map(function (row) {
      return '<div class="gfi-row">' +
        '<div class="gfi-key">' + escapeHtml(row[0]) + '</div>' +
        '<div class="gfi-value">' + escapeHtml(row[1]) + '</div>' +
        copyButton(row[1]) +
        '</div>';
    }).join('');
  }

  function renderLocationParts(feature) {
    var parts = Array.isArray(feature.location_parts) ? feature.location_parts : [];
    if (!parts.length) return '';
    var rows = parts.map(function (part, index) {
      var display = String(part && part.display || '').trim();
      var strand = String(part && part.strand || '').trim();
      return ['Part ' + String(index + 1), strand ? display + ' (' + strand + ')' : display];
    }).filter(function (row) {
      return row[1];
    });
    if (!rows.length) return '';
    return '<div class="gfi-block"><div class="gfi-block-title">Location parts</div>' + renderRows(rows) + '</div>';
  }

  function renderQualifiers(feature) {
    var rows = qualifierRows(feature);
    if (!rows.length) {
      return '<div class="gfi-empty">No qualifiers available.</div>';
    }
    return rows.map(function (row) {
      return '<div class="gfi-block">' +
        '<div class="gfi-block-title"><span>' + escapeHtml(row.key) + '</span><span>' + row.values.length + '</span>' + copyButton(row.text) + '</div>' +
        '<pre class="gfi-pre">' + escapeHtml(row.text) + '</pre>' +
        '</div>';
    }).join('');
  }

  function renderSequenceBlock(title, sequence) {
    var text = String(sequence || '');
    if (!text) {
      return '<div class="gfi-block"><div class="gfi-block-title">' + escapeHtml(title) + '</div><div class="gfi-empty">No sequence available.</div></div>';
    }
    return '<div class="gfi-block">' +
      '<div class="gfi-block-title"><span>' + escapeHtml(title) + '</span>' + copyButton(text) + '</div>' +
      '<pre class="gfi-pre">' + escapeHtml(text) + '</pre>' +
      '</div>';
  }

  function renderSequences(feature) {
    var warnings = normalizeArray(feature.sequence_warnings);
    var warningHtml = warnings.map(function (warning) {
      return '<div class="gfi-warning">' + escapeHtml(warning) + '</div>';
    }).join('');
    return warningHtml +
      renderSequenceBlock('Nucleotide', feature.nucleotide_sequence) +
      renderSequenceBlock('Amino acid', feature.amino_acid_sequence);
  }

  function renderPopup(feature, activeTab) {
    copyValues = [];
    var tab = activeTab || 'details';
    var panel = '';
    if (tab === 'qualifiers') {
      panel = renderQualifiers(feature);
    } else if (tab === 'sequence') {
      panel = renderSequences(feature);
    } else {
      panel = renderRows(detailRows(feature)) + renderLocationParts(feature);
    }
    function tabButton(id, label) {
      return '<button type="button" class="gfi-tab' + (tab === id ? ' is-active' : '') + '" data-tab="' + id + '">' + label + '</button>';
    }
    return '<div class="gfi">' +
      '<div class="gfi-header">' +
      '<div><div class="gfi-title">' + escapeHtml(feature.label || feature.svg_id || 'Feature') + '</div>' +
      '<div class="gfi-subtitle">' + escapeHtml(locationText(feature)) + '</div></div>' +
      '<button type="button" class="gfi-close" data-close="true">x</button>' +
      '</div>' +
      '<div class="gfi-tabs">' +
      tabButton('details', 'Details') +
      tabButton('qualifiers', 'Qualifiers') +
      tabButton('sequence', 'Sequence') +
      '</div>' +
      '<div class="gfi-content">' + panel + '</div>' +
      '<button type="button" class="gfi-resize-handle" data-resize="true" title="Drag to resize" aria-label="Resize popup"></button>' +
      '</div>';
  }

  function stopPopupResize() {
    if (!activePopupResize) return;
    document.removeEventListener('mousemove', activePopupResize.onMove);
    document.removeEventListener('mouseup', activePopupResize.onEnd);
    activePopupResize = null;
  }

  function closePopup() {
    stopPopupResize();
    if (popup && popup.parentNode) {
      popup.parentNode.removeChild(popup);
    }
    popup = null;
  }

  function getViewRect() {
    var viewBox = svg.viewBox && svg.viewBox.baseVal;
    if (viewBox && viewBox.width && viewBox.height) {
      return { x: viewBox.x, y: viewBox.y, width: viewBox.width, height: viewBox.height };
    }
    var width = parseFloat(svg.getAttribute('width')) || 900;
    var height = parseFloat(svg.getAttribute('height')) || 650;
    return { x: 0, y: 0, width: width, height: height };
  }

  function getScreenScale() {
    var ctm = typeof svg.getScreenCTM === 'function' ? svg.getScreenCTM() : null;
    if (!ctm) return { x: 1, y: 1 };
    return {
      x: Math.hypot(ctm.a, ctm.b) || 1,
      y: Math.hypot(ctm.c, ctm.d) || 1
    };
  }

  function getViewportClientRect() {
    var doc = document.documentElement || {};
    var width = window.innerWidth || doc.clientWidth || 900;
    var height = window.innerHeight || doc.clientHeight || 650;
    return { left: 0, top: 0, right: width, bottom: height, width: width, height: height };
  }

  function clientPoint(clientX, clientY) {
    var point = typeof svg.createSVGPoint === 'function' ? svg.createSVGPoint() : null;
    var matrix = typeof svg.getScreenCTM === 'function' ? svg.getScreenCTM() : null;
    if (point && matrix) {
      try {
        point.x = clientX;
        point.y = clientY;
        return point.matrixTransform(matrix.inverse());
      } catch (error) {
        return null;
      }
    }
    return null;
  }

  function eventPoint(event) {
    if (event && Number.isFinite(event.clientX) && Number.isFinite(event.clientY)) {
      var converted = clientPoint(event.clientX, event.clientY);
      if (converted) return converted;
    }
    var view = getViewRect();
    return { x: view.x + view.width / 2, y: view.y + view.height / 2 };
  }

  function getVisibleViewRect() {
    var viewport = getViewportClientRect();
    var points = [
      clientPoint(viewport.left, viewport.top),
      clientPoint(viewport.right, viewport.top),
      clientPoint(viewport.right, viewport.bottom),
      clientPoint(viewport.left, viewport.bottom)
    ].filter(function (point) {
      return point && Number.isFinite(point.x) && Number.isFinite(point.y);
    });
    if (points.length !== 4) return getViewRect();

    var minX = Math.min.apply(null, points.map(function (point) { return point.x; }));
    var maxX = Math.max.apply(null, points.map(function (point) { return point.x; }));
    var minY = Math.min.apply(null, points.map(function (point) { return point.y; }));
    var maxY = Math.max.apply(null, points.map(function (point) { return point.y; }));
    if (maxX <= minX || maxY <= minY) return getViewRect();

    var view = getViewRect();
    var x1 = Math.max(view.x, minX);
    var y1 = Math.max(view.y, minY);
    var x2 = Math.min(view.x + view.width, maxX);
    var y2 = Math.min(view.y + view.height, maxY);
    if (x2 <= x1 || y2 <= y1) {
      return { x: minX, y: minY, width: maxX - minX, height: maxY - minY };
    }
    return { x: x1, y: y1, width: x2 - x1, height: y2 - y1 };
  }

  function getBrowserZoomScale(viewport) {
    var scales = [1];
    var visualViewport = window.visualViewport || null;
    var visualScale = visualViewport ? Number(visualViewport.scale) : 1;
    if (Number.isFinite(visualScale) && visualScale > 1) {
      scales.push(visualScale);
    }
    var devicePixelRatio = Number(window.devicePixelRatio) || 1;
    if (baseDevicePixelRatio > 0 && devicePixelRatio > baseDevicePixelRatio) {
      scales.push(devicePixelRatio / baseDevicePixelRatio);
    }
    var outerWidth = Number(window.outerWidth);
    var viewportWidth = viewport ? Number(viewport.width) : 0;
    if (Number.isFinite(outerWidth) && outerWidth > 0 && Number.isFinite(viewportWidth) && viewportWidth > 0) {
      var widthScale = outerWidth / viewportWidth;
      if (widthScale > 1.15) {
        scales.push(widthScale);
      }
    }
    return Math.max.apply(null, scales);
  }

  function getPopupTextScale(viewport, popupCssWidth, popupCssHeight) {
    var fitScale = Math.min(1, popupCssWidth / 460, popupCssHeight / 540);
    var zoomScale = getBrowserZoomScale(viewport);
    var textScale = Math.min(fitScale, 1 / zoomScale);
    return Math.max(0.38, Math.min(1, textScale));
  }

  function clampValue(value, min, max) {
    var safeMin = Number.isFinite(min) ? min : 0;
    var safeMax = Number.isFinite(max) ? Math.max(safeMin, max) : safeMin;
    var numeric = Number(value);
    if (!Number.isFinite(numeric)) return safeMin;
    return Math.min(Math.max(numeric, safeMin), safeMax);
  }

  function setPopupTextScale(root, textScale) {
    var scale = Math.max(0.38, Math.min(1, Number(textScale) || 1));
    root.style.setProperty('--gfi-text-scale', String(scale));
    [
      ['--gfi-font-size', 13],
      ['--gfi-subtitle-font-size', 11],
      ['--gfi-close-font-size', 20],
      ['--gfi-tab-font-size', 11],
      ['--gfi-key-font-size', 10],
      ['--gfi-value-font-size', 11],
      ['--gfi-copy-font-size', 10],
      ['--gfi-block-title-font-size', 10],
      ['--gfi-pre-font-size', 11]
    ].forEach(function (entry) {
      root.style.setProperty(entry[0], (entry[1] * scale).toFixed(3) + 'px');
    });
  }

  function openPopup(feature, event) {
    closePopup();
    var viewport = getViewportClientRect();
    var view = getVisibleViewRect();
    var scale = getScreenScale();
    var safeScaleX = Math.max(scale.x, 0.001);
    var safeScaleY = Math.max(scale.y, 0.001);
    var marginCss = 12;
    var marginX = marginCss / safeScaleX;
    var marginY = marginCss / safeScaleY;
    var popupCssWidth = Math.max(1, Math.min(460, viewport.width - marginCss * 2));
    var popupCssHeight = Math.max(1, Math.min(540, viewport.height - marginCss * 2));
    var popupTextScale = getPopupTextScale(viewport, popupCssWidth, popupCssHeight);
    var width = popupCssWidth / safeScaleX;
    var height = popupCssHeight / safeScaleY;
    var point = eventPoint(event);
    var x = point.x + marginX;
    var y = point.y + marginY;
    var minX = view.x + marginX;
    var minY = view.y + marginY;
    var maxX = view.x + view.width - width - marginX;
    var maxY = view.y + view.height - height - marginY;
    if (maxX < minX) maxX = minX;
    if (maxY < minY) maxY = minY;
    x = Math.min(Math.max(x, minX), maxX);
    y = Math.min(Math.max(y, minY), maxY);

    var activeTab = 'details';
    var foreignObject = document.createElementNS(SVG_NS, 'foreignObject');
    foreignObject.setAttribute('id', 'gbdraw-feature-popup');
    foreignObject.setAttribute('class', 'gbdraw-feature-popup');
    foreignObject.setAttribute('x', x);
    foreignObject.setAttribute('y', y);
    foreignObject.setAttribute('width', width);
    foreignObject.setAttribute('height', height);

    var root = document.createElementNS(XHTML_NS, 'div');
    root.setAttribute('xmlns', XHTML_NS);
    root.style.width = popupCssWidth + 'px';
    root.style.height = popupCssHeight + 'px';
    root.style.transformOrigin = '0 0';
    root.style.transform = 'scale(' + (1 / safeScaleX) + ', ' + (1 / safeScaleY) + ')';
    setPopupTextScale(root, popupTextScale);

    function getResizeLimits() {
      viewport = getViewportClientRect();
      var rect = typeof foreignObject.getBoundingClientRect === 'function'
        ? foreignObject.getBoundingClientRect()
        : null;
      var left = rect && Number.isFinite(rect.left) ? rect.left : marginCss;
      var top = rect && Number.isFinite(rect.top) ? rect.top : marginCss;
      var minWidth = Math.min(300, Math.max(1, viewport.width - marginCss * 2));
      var minHeight = Math.min(220, Math.max(1, viewport.height - marginCss * 2));
      return {
        minWidth: minWidth,
        minHeight: minHeight,
        maxWidth: Math.max(minWidth, viewport.width - left - marginCss),
        maxHeight: Math.max(minHeight, viewport.height - top - marginCss)
      };
    }

    function applyPopupCssSize(nextCssWidth, nextCssHeight) {
      var limits = getResizeLimits();
      popupCssWidth = clampValue(nextCssWidth, limits.minWidth, limits.maxWidth);
      popupCssHeight = clampValue(nextCssHeight, limits.minHeight, limits.maxHeight);
      width = popupCssWidth / safeScaleX;
      height = popupCssHeight / safeScaleY;
      foreignObject.setAttribute('width', width);
      foreignObject.setAttribute('height', height);
      root.style.width = popupCssWidth + 'px';
      root.style.height = popupCssHeight + 'px';
      setPopupTextScale(root, getPopupTextScale(viewport, popupCssWidth, popupCssHeight));
    }

    function startPopupResize(event) {
      stopPopupResize();
      var startClientX = event.clientX;
      var startClientY = event.clientY;
      var startWidth = popupCssWidth;
      var startHeight = popupCssHeight;
      var onMove = function (moveEvent) {
        applyPopupCssSize(
          startWidth + (moveEvent.clientX - startClientX),
          startHeight + (moveEvent.clientY - startClientY)
        );
        moveEvent.preventDefault();
      };
      var onEnd = function () {
        stopPopupResize();
      };
      activePopupResize = {
        onMove: onMove,
        onEnd: onEnd
      };
      document.addEventListener('mousemove', onMove);
      document.addEventListener('mouseup', onEnd);
      event.preventDefault();
      event.stopPropagation();
    }

    function redraw() {
      root.innerHTML = renderPopup(feature, activeTab);
    }

    root.addEventListener('mousedown', function (rootEvent) {
      rootEvent.stopPropagation();
      var closest = rootEvent.target && rootEvent.target.closest
        ? rootEvent.target.closest.bind(rootEvent.target)
        : function () { return null; };
      var resizeHandle = closest('[data-resize]');
      if (!resizeHandle) return;
      startPopupResize(rootEvent);
    });

    root.addEventListener('click', function (rootEvent) {
      rootEvent.stopPropagation();
      var closest = rootEvent.target && rootEvent.target.closest
        ? rootEvent.target.closest.bind(rootEvent.target)
        : function () { return null; };
      var closeButton = closest('[data-close]');
      if (closeButton) {
        closePopup();
        return;
      }
      var tabButton = closest('[data-tab]');
      if (tabButton) {
        activeTab = tabButton.getAttribute('data-tab') || 'details';
        redraw();
        return;
      }
      var copyTarget = closest('[data-copy-index]');
      if (copyTarget) {
        var index = Number(copyTarget.getAttribute('data-copy-index'));
        var value = Number.isFinite(index) ? copyValues[index] || '' : '';
        Promise.resolve(copyText(value, copyTarget)).catch(function () {});
      }
    });

    redraw();
    foreignObject.appendChild(root);
    svg.appendChild(foreignObject);
    popup = foreignObject;
  }

  async function copyText(value, button) {
    var text = String(value == null ? '' : value);
    if (navigator.clipboard && navigator.clipboard.writeText) {
      try {
        await navigator.clipboard.writeText(text);
        if (button) button.textContent = 'Copied';
        return;
      } catch (error) {
        // Fall back to a prompt below when clipboard permission is unavailable.
      }
    }
    window.prompt('Copy text', text);
  }

  function closestFeature(node) {
    while (node && node !== svg) {
      if (node.matches && node.matches(FEATURE_SELECTOR)) return node;
      node = node.parentNode;
    }
    return null;
  }

  svg.addEventListener('click', function (event) {
    if (popup && popup.contains(event.target)) return;
    var featureElement = closestFeature(event.target);
    if (!featureElement) {
      closePopup();
      return;
    }
    var feature = featuresById.get(String(featureElement.id || ''));
    if (!feature) return;
    event.preventDefault();
    event.stopPropagation();
    openPopup(feature, event);
  });

  document.addEventListener('keydown', function (event) {
    if (event.key === 'Escape') {
      closePopup();
    }
  });
}());
`;

const normalizeStringArray = (value) => {
  if (Array.isArray(value)) {
    return value
      .filter((item) => item !== null && item !== undefined)
      .map((item) => String(item));
  }
  if (value === null || value === undefined || value === '') return [];
  return [String(value)];
};

const normalizeQualifierMap = (qualifiers) => {
  const normalized = {};
  if (!qualifiers || typeof qualifiers !== 'object' || Array.isArray(qualifiers)) {
    return normalized;
  }
  Object.entries(qualifiers).forEach(([key, value]) => {
    const normalizedKey = String(key || '').trim();
    const values = normalizeStringArray(value);
    if (!normalizedKey || values.length === 0) return;
    normalized[normalizedKey] = values;
  });
  return normalized;
};

const normalizeLocationParts = (parts) => {
  if (!Array.isArray(parts)) return [];
  return parts
    .map((part) => {
      const start = Number(part?.start);
      const end = Number(part?.end);
      const display = String(part?.display || '').trim() ||
        `${Number.isFinite(start) ? start + 1 : ''}..${Number.isFinite(end) ? end : ''}`;
      return {
        start: Number.isFinite(start) ? start : null,
        end: Number.isFinite(end) ? end : null,
        strand: String(part?.strand || '').trim(),
        display
      };
    })
    .filter((part) => part.display && part.display !== '..');
};

const buildStandaloneFeatureLocation = (feature) => {
  const start = Number(feature?.start);
  const end = Number(feature?.end);
  const startText = Number.isFinite(start) ? String(start + 1) : String(feature?.start ?? '');
  const endText = Number.isFinite(end) ? String(end) : String(feature?.end ?? '');
  const strand = String(feature?.strand || '').trim();
  const range = `${startText}..${endText}`;
  return strand ? `${range} (${strand})` : range;
};

const firstQualifierValue = (feature, key) => {
  const qualifiers = feature?.qualifiers && typeof feature.qualifiers === 'object'
    ? feature.qualifiers
    : {};
  const values = normalizeStringArray(qualifiers[key]);
  return values.find((value) => value.trim()) || '';
};

const getStandaloneFeatureLabel = (feature) => {
  const candidates = [
    feature?.label,
    feature?.gene,
    feature?.locus_tag,
    firstQualifierValue(feature, 'gene'),
    firstQualifierValue(feature, 'locus_tag'),
    firstQualifierValue(feature, 'product'),
    feature?.product,
    feature?.type,
    feature?.svg_id
  ];
  for (const candidate of candidates) {
    const label = String(candidate || '').trim();
    if (label) return label;
  }
  return 'Feature';
};

const collectRenderedFeatureIds = (svg) => {
  const ids = new Set();
  if (!svg) return ids;
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = String(element.getAttribute('id') || '').trim();
    if (id) ids.add(id);
  });
  return ids;
};

const buildStandaloneFeaturePayloads = (svg) => {
  const renderedIds = collectRenderedFeatureIds(svg);
  if (renderedIds.size === 0) return [];

  const features = Array.isArray(state.extractedFeatures.value) ? state.extractedFeatures.value : [];
  const payloads = [];
  const seenIds = new Set();
  features.forEach((feature) => {
    const svgId = String(feature?.svg_id || '').trim();
    if (!svgId || !renderedIds.has(svgId) || seenIds.has(svgId)) return;
    seenIds.add(svgId);
    payloads.push({
      svg_id: svgId,
      label: getStandaloneFeatureLabel(feature),
      record_id: String(feature?.record_id || ''),
      record_idx: Number.isFinite(Number(feature?.record_idx)) ? Number(feature.record_idx) : null,
      type: String(feature?.type || ''),
      start: Number.isFinite(Number(feature?.start)) ? Number(feature.start) : null,
      end: Number.isFinite(Number(feature?.end)) ? Number(feature.end) : null,
      strand: String(feature?.strand || ''),
      location: buildStandaloneFeatureLocation(feature),
      qualifiers: normalizeQualifierMap(feature?.qualifiers),
      location_parts: normalizeLocationParts(feature?.location_parts),
      nucleotide_sequence: String(feature?.nucleotide_sequence || ''),
      amino_acid_sequence: String(feature?.amino_acid_sequence || ''),
      sequence_warnings: normalizeStringArray(feature?.sequence_warnings)
    });
  });
  return payloads;
};

const removeExistingStandaloneFeaturePopupAssets = (svg) => {
  [INTERACTIVE_METADATA_ID, INTERACTIVE_STYLE_ID, INTERACTIVE_SCRIPT_ID, 'gbdraw-feature-popup'].forEach((id) => {
    const element = svg.querySelector(`#${CSS.escape(id)}`);
    if (element?.parentNode) {
      element.parentNode.removeChild(element);
    }
  });
};

const addClassToken = (element, token) => {
  const existing = String(element.getAttribute('class') || '').trim();
  const tokens = new Set(existing ? existing.split(/\s+/) : []);
  tokens.add(token);
  element.setAttribute('class', Array.from(tokens).join(' '));
};

const enrichSvgWithStandaloneFeaturePopup = (svg) => {
  if (!svg || state.adv.rich_feature_popup === false) return false;

  removeExistingStandaloneFeaturePopupAssets(svg);
  const features = buildStandaloneFeaturePayloads(svg);
  if (features.length === 0) return false;

  const featureIds = new Set(features.map((feature) => feature.svg_id));
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = String(element.getAttribute('id') || '').trim();
    if (!featureIds.has(id)) return;
    element.setAttribute('data-gbdraw-interactive-feature', 'true');
    addClassToken(element, 'gbdraw-interactive-feature');
  });

  const metadata = document.createElementNS(SVG_NS, 'metadata');
  metadata.setAttribute('id', INTERACTIVE_METADATA_ID);
  metadata.setAttribute('data-schema', 'gbdraw-interactive-feature-popup-v1');
  metadata.textContent = JSON.stringify({
    schema: 'gbdraw-interactive-feature-popup-v1',
    features
  });

  const style = document.createElementNS(SVG_NS, 'style');
  style.setAttribute('id', INTERACTIVE_STYLE_ID);
  style.setAttribute('type', 'text/css');
  style.textContent = STANDALONE_INTERACTIVE_STYLE;

  const script = document.createElementNS(SVG_NS, 'script');
  script.setAttribute('id', INTERACTIVE_SCRIPT_ID);
  script.setAttribute('type', 'application/ecmascript');
  script.textContent = STANDALONE_INTERACTIVE_SCRIPT;

  svg.appendChild(metadata);
  svg.appendChild(style);
  svg.appendChild(script);
  svg.setAttribute('data-gbdraw-interactive-svg', 'true');
  return true;
};

const ensureSvgDefs = (svg) => {
  let defs = svg.querySelector('defs');
  if (!defs) {
    defs = document.createElementNS(SVG_NS, 'defs');
    svg.insertBefore(defs, svg.firstChild);
  }
  return defs;
};

const moveGradientsToDefs = (svg) => {
  const defs = ensureSvgDefs(svg);
  svg.querySelectorAll('linearGradient, radialGradient').forEach((gradient) => {
    if (gradient.parentNode !== defs) {
      defs.appendChild(gradient);
    }
  });
};

const copyAttributes = (target, source, attributes) => {
  attributes.forEach((attr) => {
    const value = source.getAttribute(attr);
    if (value !== null) {
      target.setAttribute(attr, value);
    }
  });
};

const applyFontAndPaintStyles = (target, source) => {
  if (!source) return;
  const style = window.getComputedStyle(source);
  [
    'font-family',
    'font-size',
    'font-weight',
    'font-style',
    'font-variant',
    'fill',
    'fill-opacity',
    'stroke',
    'stroke-opacity',
    'stroke-width'
  ].forEach((prop) => {
    const value = style.getPropertyValue(prop);
    if (value) {
      target.setAttribute(prop, value.trim());
    }
  });
};

const flattenTextPathsForPdf = (svg) => {
  const textElements = Array.from(svg.querySelectorAll('text')).filter((el) =>
    el.querySelector('textPath')
  );
  if (!textElements.length) return;

  textElements.forEach((textEl) => {
    const textPath = textEl.querySelector('textPath');
    if (!textPath) return;
    const textContent = textPath.textContent || textEl.textContent || '';
    const charCount =
      typeof textEl.getNumberOfChars === 'function' ? textEl.getNumberOfChars() : 0;
    if (!textContent || !charCount) return;

    const chars = Array.from(textContent);
    const group = document.createElementNS(SVG_NS, 'g');
    copyAttributes(group, textEl, ['transform', 'opacity', 'display', 'visibility']);
    applyFontAndPaintStyles(group, textPath);

    for (let i = 0; i < charCount && i < chars.length; i += 1) {
      const char = chars[i];
      if (char === '\n' || char === '\r') continue;
      let pos;
      try {
        pos = textEl.getStartPositionOfChar(i);
      } catch (error) {
        pos = null;
      }
      if (!pos) continue;
      const rotation = textEl.getRotationOfChar(i) || 0;
      const charText = document.createElementNS(SVG_NS, 'text');
      charText.textContent = char;
      charText.setAttribute('x', pos.x);
      charText.setAttribute('y', pos.y);
      charText.setAttribute('text-anchor', 'start');
      charText.setAttribute('dominant-baseline', 'alphabetic');
      if (rotation) {
        charText.setAttribute('transform', `rotate(${rotation} ${pos.x} ${pos.y})`);
      }
      group.appendChild(charText);
    }

    if (textEl.parentNode) {
      textEl.parentNode.replaceChild(group, textEl);
    }
  });
};

const normalizeTextBaselinesForPdf = (svg) => {
  svg.querySelectorAll('text').forEach((textEl) => {
    const baseline = textEl.getAttribute('dominant-baseline');
    if (!baseline || baseline === 'alphabetic' || baseline === 'auto') return;

    let before;
    try {
      before = textEl.getBBox();
    } catch (error) {
      return;
    }

    textEl.setAttribute('dominant-baseline', 'alphabetic');

    let after;
    try {
      after = textEl.getBBox();
    } catch (error) {
      textEl.setAttribute('dominant-baseline', baseline);
      return;
    }

    const dy = before.y - after.y;
    if (!Number.isFinite(dy) || Math.abs(dy) < 0.01) {
      textEl.removeAttribute('dominant-baseline');
      return;
    }

    const yAttr = textEl.getAttribute('y');
    if (yAttr) {
      const yValues = yAttr.split(/[\s,]+/).map((value) => parseFloat(value));
      const adjusted = yValues.map((value) =>
        Number.isFinite(value) ? value + dy : value
      );
      textEl.setAttribute('y', adjusted.join(' '));
    } else {
      const dyAttr = textEl.getAttribute('dy');
      const dyValue = dyAttr ? parseFloat(dyAttr) : 0;
      const nextDy = Number.isFinite(dyValue) ? dyValue + dy : dy;
      textEl.setAttribute('dy', nextDy);
    }

    textEl.removeAttribute('dominant-baseline');
  });
};

const prepareSvgForPdf = (svg) => {
  if (!svg) return null;
  const clone = svg.cloneNode(true);
  if (!clone.getAttribute('xmlns')) {
    clone.setAttribute('xmlns', SVG_NS);
  }
  if (!clone.getAttribute('xmlns:xlink')) {
    clone.setAttribute('xmlns:xlink', 'http://www.w3.org/1999/xlink');
  }

  const staging = document.createElement('div');
  staging.style.position = 'absolute';
  staging.style.left = '-10000px';
  staging.style.top = '-10000px';
  staging.style.visibility = 'hidden';
  staging.appendChild(clone);
  document.body.appendChild(staging);

  moveGradientsToDefs(clone);
  flattenTextPathsForPdf(clone);
  normalizeTextBaselinesForPdf(clone);

  return {
    svg: clone,
    cleanup: () => staging.remove()
  };
};

export const downloadSVG = () => {
  const svgString = getCurrentSvgString({ interactive: true });
  if (!svgString) return;
  const blob = new Blob([svgString], { type: 'image/svg+xml' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = getDownloadName('svg');
  a.click();
  URL.revokeObjectURL(url);
};

export const downloadPNG = () => {
  const svgString = getCurrentSvgString();
  if (!svgString) return;
  const svg = getSvgFromString(svgString);
  if (!svg) return;
  const dims = getSvgDimensions(svg);
  if (!dims) return;
  const canvas = document.createElement('canvas');
  const dpi = parseInt(state.downloadDpi.value, 10);
  const scale = dpi / 96;
  canvas.width = dims.width * scale;
  canvas.height = dims.height * scale;
  const ctx = canvas.getContext('2d');
  const img = new Image();
  const blob = new Blob([svgString], { type: 'image/svg+xml;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  img.onload = () => {
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    canvas.toBlob(async (pngBlob) => {
      const fixedBlob = await setDpiInPng(pngBlob, dpi);
      const downloadUrl = URL.createObjectURL(fixedBlob);
      const link = document.createElement('a');
      link.download = getDownloadName('png');
      link.href = downloadUrl;
      link.click();
      URL.revokeObjectURL(downloadUrl);
      URL.revokeObjectURL(url);
    }, 'image/png');
  };
  img.onerror = () => {
    console.error('Failed to load SVG for PNG conversion');
    URL.revokeObjectURL(url);
  };
  img.src = url;
};

export const downloadPDF = async () => {
  const svgString = getCurrentSvgString();
  if (!svgString) return;

  // 1. Create temporary SVG element to read dimensions
  const svg = getSvgFromString(svgString);
  if (!svg) return;
  const prepared = prepareSvgForPdf(svg);
  if (!prepared) return;
  const { svg: pdfSvg, cleanup } = prepared;

  // 2. Determine Width/Height (Use attributes or viewBox)
  const dims = getSvgDimensions(pdfSvg);
  if (!dims) {
    cleanup();
    return;
  }

  // 3. Initialize jsPDF (Use 'pt' as unit for vector fidelity)
  const { jsPDF } = window.jspdf;
  const doc = new jsPDF({
    orientation: dims.width > dims.height ? 'l' : 'p',
    unit: 'pt',
    format: [dims.width, dims.height]
  });

  // 4. Convert SVG to PDF
  try {
    await doc.svg(pdfSvg, {
      x: 0,
      y: 0,
      width: dims.width,
      height: dims.height
    });
    // 5. Save File
    doc.save(getDownloadName('pdf'));
  } finally {
    cleanup();
  }
};
