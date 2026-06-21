import { buildFeatureSequenceFastas } from '../app/feature-sequence-fasta.js';

const SVG_NS = 'http://www.w3.org/2000/svg';
const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
const FEATURE_SELECTOR = [
  `path[${FEATURE_ID_ATTRIBUTE}]`,
  `polygon[${FEATURE_ID_ATTRIBUTE}]`,
  `rect[${FEATURE_ID_ATTRIBUTE}]`,
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');
const STANDALONE_MATCH_SELECTOR = [
  'path[data-gbdraw-pairwise-match-id]',
  'path[data-match-kind]',
  'path[data-pairwise-match-style]'
].join(', ');
const INTERACTIVE_METADATA_ID = 'gbdraw-interactive-feature-metadata';
const INTERACTIVE_STYLE_ID = 'gbdraw-interactive-feature-style';
const INTERACTIVE_SCRIPT_ID = 'gbdraw-interactive-feature-script';
const INTERACTIVE_GLOW_FILTER_ID = 'gbdraw-interactive-feature-glow';
const INTERACTIVE_MATCH_GLOW_FILTER_ID = 'gbdraw-interactive-feature-match-glow';
const FEATURE_PART_SUFFIX_RE = /__part\d+$/;

const normalizeFeatureElementId = (value) =>
  String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');

export const stripEditorOnlyCursorStyles = (svg) => {
  if (!svg) return;
  svg.querySelectorAll('[style]').forEach((element) => {
    const style = element.getAttribute('style');
    if (!style || !/\bcursor\s*:/i.test(style)) return;
    element.style.removeProperty('cursor');
    if (!element.getAttribute('style')?.trim()) {
      element.removeAttribute('style');
    }
  });
};

const getElementFeatureId = (element) =>
  normalizeFeatureElementId(
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    ''
  );

const STANDALONE_INTERACTIVE_STYLE = `
.gbdraw-interactive-feature {
  cursor: pointer;
  transition: opacity 120ms ease, filter 120ms ease, stroke 120ms ease, stroke-width 120ms ease, stroke-opacity 120ms ease;
}
.gbdraw-interactive-feature:hover,
.gbdraw-interactive-feature.gbdraw-interactive-feature--hover {
  opacity: 0.82;
  filter: url(#gbdraw-interactive-feature-glow);
}
.gbdraw-interactive-feature.gbdraw-interactive-feature--dimmed {
  opacity: 0.18;
}
.gbdraw-interactive-feature.gbdraw-interactive-feature--match {
  opacity: 1;
  filter: url(#gbdraw-interactive-feature-match-glow);
  stroke: #fbbf24;
  stroke-linejoin: round;
  stroke-opacity: 0.6;
  stroke-width: 1;
  paint-order: stroke fill markers;
}
.gbdraw-interactive-feature.gbdraw-interactive-feature--active-match {
  opacity: 1;
  filter: url(#gbdraw-interactive-feature-glow);
  stroke: #f59e0b;
  stroke-opacity: 1;
  stroke-width: 2;
}
.gbdraw-interactive-orthogroup-link {
  transition: opacity 120ms ease, filter 120ms ease;
}
.gbdraw-interactive-orthogroup-link.gbdraw-interactive-orthogroup-link--hover {
  opacity: 0.92;
  filter: url(#gbdraw-interactive-feature-glow);
}
.gbdraw-interactive-pairwise-match {
  cursor: pointer;
  transition: opacity 120ms ease, filter 120ms ease, stroke 120ms ease, stroke-width 120ms ease, stroke-opacity 120ms ease;
}
.gbdraw-interactive-pairwise-match.gbdraw-interactive-pairwise-match--hover {
  opacity: 0.92;
  filter: url(#gbdraw-interactive-feature-glow);
  stroke: #f59e0b;
  stroke-opacity: 0.9;
  stroke-width: 1.5;
  paint-order: stroke fill markers;
}
.gbdraw-feature-popup {
  overflow: visible;
  pointer-events: auto;
}
.gbdraw-feature-hover-popup {
  overflow: visible;
  pointer-events: none;
}
.gbdraw-feature-popup * {
  box-sizing: border-box;
}
.gbdraw-feature-hover-popup * {
  box-sizing: border-box;
}
.gfhs {
  width: 100%;
  height: 100%;
  overflow: hidden;
  border: 1px solid #cbd5e1;
  border-radius: 8px;
  background: rgba(255, 255, 255, 0.98);
  color: #334155;
  box-shadow: 0 18px 45px rgba(15, 23, 42, 0.18);
  font-family: Arial, Helvetica, sans-serif;
  font-size: 12px;
  line-height: 1.35;
  padding: 10px 12px;
}
.gfhs-title {
  display: flex;
  gap: 8px;
  align-items: flex-start;
  min-width: 0;
  margin-bottom: 7px;
  padding-bottom: 7px;
  border-bottom: 1px solid #e2e8f0;
}
.gfhs-swatch {
  flex: 0 0 auto;
  width: 12px;
  height: 12px;
  margin-top: 3px;
  border: 1px solid rgba(15, 23, 42, 0.2);
  border-radius: 3px;
  background: #94a3b8;
}
.gfhs-text {
  min-width: 0;
}
.gfhs-heading {
  min-width: 0;
  overflow: hidden;
  color: #0f172a;
  font-weight: 700;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.gfhs-subtitle {
  margin-top: 2px;
  overflow: hidden;
  color: #64748b;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 10px;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.gfhs-row {
  display: grid;
  grid-template-columns: 76px minmax(0, 1fr);
  gap: 8px;
  align-items: start;
  margin-top: 4px;
}
.gfhs-key {
  color: #64748b;
  font-size: 10px;
  font-weight: 700;
  text-transform: uppercase;
}
.gfhs-value {
  min-width: 0;
  overflow-wrap: anywhere;
}
.gfhs-value.is-clamped {
  display: -webkit-box;
  overflow: hidden;
  -webkit-box-orient: vertical;
  -webkit-line-clamp: 2;
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
  font-size: 13px;
}
.gfi--simple .gfi-content {
  padding-right: 12px;
}
.gfi--simple .gfi-row {
  grid-template-columns: 82px minmax(0, 1fr) 50px;
}
.gfi-header {
  display: grid;
  grid-template-columns: minmax(0, 1fr) 28px;
  gap: 8px;
  align-items: start;
  padding: 10px 12px;
  border-bottom: 1px solid #e2e8f0;
  cursor: grab;
  user-select: none;
}
.gfi-header:active {
  cursor: grabbing;
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
  font-size: 11px;
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
  font-size: 20px;
  line-height: 1;
  cursor: pointer;
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
  font-size: 11px;
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
  overscroll-behavior: contain;
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
  font-size: 10px;
  font-weight: 700;
  text-transform: uppercase;
}
.gfi-value {
  min-width: 0;
  color: #1e293b;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 11px;
  overflow-wrap: anywhere;
  white-space: pre-wrap;
}
.gfi-copy {
  min-width: 48px;
  padding: 5px 0;
  background: #f8fafc;
  font-size: 10px;
  font-weight: 700;
  white-space: nowrap;
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
  font-size: 10px;
  font-weight: 700;
  text-transform: uppercase;
}
.gfi-block-title .gfi-copy {
  margin-left: auto;
}
.gfi-block-title .gfi-copy + .gfi-copy {
  margin-left: 0;
}
.gfi-pre {
  max-height: 120px;
  margin: 0;
  padding: 8px 9px;
  overflow: auto;
  color: #1e293b;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 11px;
  line-height: 1.45;
  overflow-wrap: anywhere;
  white-space: pre-wrap;
}
.gfi-table-wrap {
  max-height: 180px;
  overflow: auto;
}
.gfi-table {
  width: 100%;
  min-width: 620px;
  border-collapse: collapse;
  font-size: 11px;
}
.gfi-table th,
.gfi-table td {
  padding: 6px 8px;
  border-bottom: 1px solid #f1f5f9;
  text-align: left;
  vertical-align: top;
}
.gfi-table th {
  position: sticky;
  top: 0;
  z-index: 1;
  background: #f8fafc;
  color: #64748b;
  font-size: 10px;
  font-weight: 700;
  text-transform: uppercase;
}
.gfi-table td {
  color: #1e293b;
  overflow-wrap: anywhere;
}
.gfi-table td.gfi-mono {
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  white-space: nowrap;
}
.gfi-block-og-row {
  cursor: pointer;
}
.gfi-block-og-row:hover,
.gfi-block-og-row.is-active {
  background: #ecfdf5;
}
.gfi-block--selected-og {
  margin: 8px;
}
.gfi-block-actions {
  display: flex;
  justify-content: flex-end;
  gap: 6px;
  padding: 6px 8px;
  background: #f8fafc;
}
.gfi-seq-actions {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 4px;
  min-width: 118px;
}
.gfi-seq-actions .gfi-copy {
  min-width: 0;
  padding: 4px 5px;
}
.gfi-match-feature-table {
  width: 100%;
  min-width: 560px;
  border-collapse: collapse;
  font-size: 11px;
}
.gfi-match-feature-table th,
.gfi-match-feature-table td {
  padding: 6px 8px;
  border-bottom: 1px solid #f1f5f9;
  text-align: left;
  vertical-align: top;
}
.gfi-match-feature-table th {
  position: sticky;
  top: 0;
  z-index: 1;
  background: #ffffff;
  color: #94a3b8;
  font-size: 10px;
  font-weight: 700;
  text-transform: uppercase;
}
.gfi-match-feature-row {
  cursor: pointer;
}
.gfi-match-feature-row:hover {
  background: #eff6ff;
}
.gfi-match-feature-row.is-disabled {
  cursor: default;
}
.gfi-match-feature-row.is-disabled:hover {
  background: transparent;
}
.gfi-match-feature-main {
  min-width: 0;
  color: #334155;
  font-weight: 700;
  overflow-wrap: anywhere;
}
.gfi-match-feature-sub {
  margin-top: 2px;
  color: #94a3b8;
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 10px;
  overflow-wrap: anywhere;
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
.gbdraw-viewport-controls {
  pointer-events: auto;
}
.gbdraw-sticky-legend {
  pointer-events: none;
}
.gbdraw-sticky-legend-background {
  pointer-events: none;
}
.gbdraw-viewport-button {
  cursor: pointer;
}
.gbdraw-viewport-button rect {
  fill: #ffffff;
  stroke: #cbd5e1;
  stroke-width: 1;
  rx: 7;
  ry: 7;
}
.gbdraw-viewport-button text {
  fill: #475569;
  font-family: Arial, Helvetica, sans-serif;
  font-size: 13px;
  font-weight: 700;
  text-anchor: middle;
  dominant-baseline: middle;
  pointer-events: none;
  user-select: none;
}
.gbdraw-viewport-button:hover rect {
  fill: #eff6ff;
  stroke: #93c5fd;
}
.gbdraw-viewport-button.is-active rect {
  fill: #dbeafe;
  stroke: #2563eb;
}
.gbdraw-viewport-button.is-active text {
  fill: #1d4ed8;
}
.gbdraw-interactive-pan-enabled {
  cursor: grab;
}
.gbdraw-interactive-panning {
  cursor: grabbing;
}
.gbdraw-interactive-pan-enabled .gbdraw-interactive-feature {
  cursor: pointer;
}
.gbdraw-feature-search-controls {
  overflow: visible;
  pointer-events: auto;
}
.gbdraw-feature-search-controls * {
  box-sizing: border-box;
}
.gfs {
  width: 100%;
  height: 100%;
  padding: 7px;
  border: 1px solid #cbd5e1;
  border-radius: 8px;
  background: rgba(255, 255, 255, 0.96);
  color: #334155;
  box-shadow: 0 12px 28px rgba(15, 23, 42, 0.16);
  font-family: Arial, Helvetica, sans-serif;
  font-size: 11px;
  cursor: grab;
}
.gfs.is-dragging {
  cursor: grabbing;
}
.gfs-row {
  display: flex;
  gap: 5px;
  align-items: center;
  min-width: 0;
}
.gfs-row + .gfs-row {
  margin-top: 6px;
}
.gfs-input,
.gfs-select {
  height: 26px;
  min-width: 0;
  border: 1px solid #cbd5e1;
  border-radius: 6px;
  background: #ffffff;
  color: #1e293b;
  font: inherit;
}
.gfs-input {
  padding: 0 7px;
}
.gfs-query {
  flex: 1 1 126px;
}
.gfs-select {
  flex: 0 0 104px;
  padding: 0 5px;
}
.gfs-qualifier {
  flex: 1 1 104px;
}
.gfs-qualifier:disabled {
  opacity: 0.48;
}
.gfs-toggle {
  display: inline-flex;
  gap: 4px;
  align-items: center;
  white-space: nowrap;
  user-select: none;
}
.gfs-button {
  height: 26px;
  min-width: 30px;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  border: 1px solid #cbd5e1;
  border-radius: 6px;
  background: #f8fafc;
  color: #334155;
  cursor: pointer;
  font: inherit;
  font-weight: 700;
}
.gfs-button--clear {
  flex: 0 0 48px;
}
.gfs-button--open {
  flex: 0 0 44px;
}
.gfs-button--search {
  flex: 0 0 62px;
  border-color: #1d4ed8;
  background: #2563eb;
  color: #ffffff;
  box-shadow: 0 8px 16px rgba(37, 99, 235, 0.25);
}
.gfs-button:hover {
  background: #eff6ff;
  border-color: #93c5fd;
  color: #1d4ed8;
}
.gfs-button--search:hover {
  background: #1d4ed8;
  border-color: #1e40af;
  color: #ffffff;
}
.gfs-button:disabled {
  cursor: not-allowed;
  opacity: 0.45;
  box-shadow: none;
}
.gfs-button:disabled:hover {
  background: #f8fafc;
  border-color: #cbd5e1;
  color: #334155;
}
.gfs-button--search:disabled:hover {
  background: #2563eb;
  border-color: #1d4ed8;
  color: #ffffff;
}
.gfs-count {
  margin-left: auto;
  color: #475569;
  font-weight: 700;
  white-space: nowrap;
}
.gfs-match-detail {
  display: block;
  min-width: 0;
  margin-top: 5px;
  overflow: hidden;
  color: #64748b;
  font-size: 10px;
  line-height: 1.2;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.gfs input,
.gfs select,
.gfs button,
.gfs label {
  cursor: auto;
}
.gfs button,
.gfs label {
  cursor: pointer;
}
.gfs-match-detail strong {
  color: #334155;
  font-weight: 700;
}
.gfs.is-invalid .gfs-query {
  border-color: #ef4444;
  color: #991b1b;
}
.gfs.is-invalid .gfs-count {
  color: #b91c1c;
}
`;

const STANDALONE_INTERACTIVE_SCRIPT = `
(function () {
  'use strict';

  var FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
  var FEATURE_SELECTOR =
    'path[' + FEATURE_ID_ATTRIBUTE + '], polygon[' + FEATURE_ID_ATTRIBUTE + '], rect[' + FEATURE_ID_ATTRIBUTE + '], ' +
    'path[id^="f"], polygon[id^="f"], rect[id^="f"]';
  var MATCH_SELECTOR =
    'path[data-gbdraw-pairwise-match-id], path[data-match-kind], path[data-pairwise-match-style]';
  var SVG_NS = 'http://www.w3.org/2000/svg';
  var XHTML_NS = 'http://www.w3.org/1999/xhtml';
  var VIEWPORT_CONTROLS_ID = 'gbdraw-viewport-controls';
  var SEARCH_CONTROLS_ID = 'gbdraw-feature-search-controls';
  var svg = document.currentScript && document.currentScript.ownerSVGElement
    ? document.currentScript.ownerSVGElement
    : document.documentElement;
  var metadata = svg.querySelector('#gbdraw-interactive-feature-metadata');
  var payload = null;
  var popup = null;
  var hoverPopup = null;
  var viewportControls = null;
  var searchControls = null;
  var stickyLegend = null;
  var stickyLegendState = null;
  var activeCanvasPan = null;
  var suppressNextCanvasClick = false;
  var activePopupResize = null;
  var activePopupDrag = null;
  var activeSearchControlsDrag = null;
  var searchControlsOffsetCss = { x: 0, y: 0 };
  var hoverPopupTimer = null;
  var hoverPopupFrame = null;
  var hoverPopupFeatureId = '';
  var hoverPopupLastEvent = null;
  var activeHoverSvgId = null;
  var activeHoverKey = '';
  var activeHoverMatchId = '';
  var activeHoverMatchKey = '';
  var maxZoom = 16;
  var resizeFrame = null;
  var overlayRefreshFrame = null;
  var overlayRefreshTimer = null;
  var updateActivePopupViewportMetrics = null;
  var searchState = {
    query: '',
    field: 'all',
    qualifierKey: '',
    useRegex: false,
    matches: [],
    matchDetails: {},
    activeIndex: -1,
    error: ''
  };
  var pendingSearchState = {
    query: '',
    field: 'all',
    qualifierKey: '',
    useRegex: false
  };
  var baseDevicePixelRatio = Math.max(1, Number(window.devicePixelRatio) || 1);
  var FEATURE_PART_SUFFIX_RE = /__part\\d+$/;

  function normalizeFeatureElementId(value) {
    return String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');
  }

  function getElementFeatureId(element) {
    return normalizeFeatureElementId(
      element && (
        element.getAttribute(FEATURE_ID_ATTRIBUTE) ||
        element.getAttribute('id') ||
        element.id
      ) || ''
    );
  }

  try {
    payload = JSON.parse(metadata ? metadata.textContent || '{}' : '{}');
  } catch (error) {
    payload = {};
  }

  var features = Array.isArray(payload.features) ? payload.features : [];
  var orthogroups = Array.isArray(payload.orthogroups) ? payload.orthogroups : [];
  var matches = Array.isArray(payload.matches) ? payload.matches : [];
  var popupMode = payload.popup_mode === 'simple' ? 'simple' : 'rich';
  var searchFieldIds = {
    all: true,
    label: true,
    type: true,
    'record-id': true,
    location: true,
    strand: true,
    orthogroup: true,
    'qualifier-key': true,
    'qualifier-value': true,
    nucleotide: true,
    'amino-acid': true
  };
  var richSearchFields = {
    'qualifier-key': true,
    'qualifier-value': true,
    'nucleotide': true,
    'amino-acid': true
  };

  function isRichSearchField(field) {
    return Boolean(richSearchFields[String(field || '')]);
  }

  function normalizeSearchField(field) {
    var selected = String(field || 'all').trim() || 'all';
    if (!searchFieldIds[selected]) {
      return 'all';
    }
    if (popupMode === 'simple' && isRichSearchField(selected)) {
      return 'all';
    }
    return selected;
  }

  var featuresById = new Map();
  features.forEach(function (feature) {
    var svgId = String(feature && feature.svg_id || '').trim();
    if (svgId && !featuresById.has(svgId)) {
      featuresById.set(svgId, feature);
    }
  });
  var orthogroupsById = new Map();
  orthogroups.forEach(function (group) {
    [
      group && group.id,
      group && group.orthogroupId,
      group && group.orthogroup_id
    ].forEach(function (candidate) {
      var id = String(candidate || '').trim();
      if (id && !orthogroupsById.has(id)) {
        orthogroupsById.set(id, group);
      }
    });
  });
  var matchesById = new Map();
  matches.forEach(function (match) {
    var id = String(match && match.id || '').trim();
    if (id && !matchesById.has(id)) {
      matchesById.set(id, match);
    }
  });

  function getOrthogroupIds(value) {
    var seen = {};
    return String(value || '').split(';').map(function (entry) {
      return String(entry || '').trim();
    }).filter(function (entry) {
      if (!entry || seen[entry]) return false;
      seen[entry] = true;
      return true;
    });
  }

  var featureElementsById = new Map();
  Array.prototype.slice.call(svg.querySelectorAll(FEATURE_SELECTOR)).forEach(function (element) {
    var svgId = getElementFeatureId(element);
    if (!svgId || !featuresById.has(svgId)) return;
    if (!featureElementsById.has(svgId)) {
      featureElementsById.set(svgId, []);
    }
    featureElementsById.get(svgId).push(element);
  });

  var featureIdsByOrthogroupId = new Map();
  features.forEach(function (feature) {
    var svgId = String(feature && feature.svg_id || '').trim();
    if (!svgId) return;
    getOrthogroupIds(feature && feature.orthogroup_id).forEach(function (orthogroupId) {
      if (!featureIdsByOrthogroupId.has(orthogroupId)) {
        featureIdsByOrthogroupId.set(orthogroupId, new Set());
      }
      featureIdsByOrthogroupId.get(orthogroupId).add(svgId);
    });
  });

  var comparisonElementsByOrthogroupId = new Map();
  Array.prototype.slice.call(svg.querySelectorAll('[data-orthogroup-id]')).forEach(function (element) {
    if (element.matches && element.matches(FEATURE_SELECTOR)) return;
    getOrthogroupIds(element.getAttribute('data-orthogroup-id')).forEach(function (orthogroupId) {
      if (!comparisonElementsByOrthogroupId.has(orthogroupId)) {
        comparisonElementsByOrthogroupId.set(orthogroupId, []);
      }
      comparisonElementsByOrthogroupId.get(orthogroupId).push(element);
    });
    setClassToken(element, 'gbdraw-interactive-orthogroup-link', true);
  });

  var matchElementsById = new Map();
  var comparisonElementsByCollinearityBlockId = new Map();
  Array.prototype.slice.call(svg.querySelectorAll(MATCH_SELECTOR)).forEach(function (element) {
    var matchId = String(element.getAttribute('data-gbdraw-pairwise-match-id') || '').trim();
    if (matchId) {
      if (!matchElementsById.has(matchId)) {
        matchElementsById.set(matchId, []);
      }
      matchElementsById.get(matchId).push(element);
    }
    var blockId = String(element.getAttribute('data-collinearity-block-id') || '').trim();
    if (blockId) {
      if (!comparisonElementsByCollinearityBlockId.has(blockId)) {
        comparisonElementsByCollinearityBlockId.set(blockId, []);
      }
      comparisonElementsByCollinearityBlockId.get(blockId).push(element);
    }
    setClassToken(element, 'gbdraw-interactive-pairwise-match', true);
  });

  function setClassToken(element, token, enabled) {
    if (!element) return;
    if (element.classList && element.classList.toggle) {
      element.classList.toggle(token, Boolean(enabled));
      return;
    }
    var existing = String(element.getAttribute('class') || '').trim();
    var tokens = existing ? existing.split(/\\s+/) : [];
    var nextTokens = tokens.filter(function (entry) {
      return entry && entry !== token;
    });
    if (enabled) nextTokens.push(token);
    element.setAttribute('class', nextTokens.join(' '));
  }

  featureElementsById.forEach(function (elements) {
    elements.forEach(function (element) {
      setClassToken(element, 'gbdraw-interactive-feature', true);
    });
  });

  function setFeatureHighlight(svgId, highlight) {
    var elements = featureElementsById.get(String(svgId || '')) || [];
    elements.forEach(function (element) {
      setClassToken(element, 'gbdraw-interactive-feature--hover', highlight);
    });
  }

  function getFeatureHoverKey(svgId) {
    var feature = featuresById.get(String(svgId || '').trim());
    var orthogroupId = String(feature && feature.orthogroup_id || '').trim();
    return orthogroupId ? 'orthogroup:' + orthogroupId : 'feature:' + String(svgId || '').trim();
  }

  function setOrthogroupHover(orthogroupId, highlight) {
    var id = String(orthogroupId || '').trim();
    if (!id) return;
    var featureIds = featureIdsByOrthogroupId.get(id) || new Set();
    featureIds.forEach(function (featureId) {
      setFeatureHighlight(featureId, highlight);
    });
    var comparisonElements = comparisonElementsByOrthogroupId.get(id) || [];
    comparisonElements.forEach(function (element) {
      setClassToken(element, 'gbdraw-interactive-orthogroup-link--hover', highlight);
    });
  }

  function setCollinearityBlockHover(blockId, highlight) {
    var id = String(blockId || '').trim();
    if (!id) return;
    var comparisonElements = comparisonElementsByCollinearityBlockId.get(id) || [];
    comparisonElements.forEach(function (element) {
      setClassToken(element, 'gbdraw-interactive-pairwise-match--hover', highlight);
    });
  }

  function setHoverHighlight(svgId, highlight) {
    var feature = featuresById.get(String(svgId || '').trim());
    var orthogroupId = String(feature && feature.orthogroup_id || '').trim();
    if (orthogroupId) {
      setOrthogroupHover(orthogroupId, highlight);
      return;
    }
    setFeatureHighlight(svgId, highlight);
  }

  function matchAttr(element, name) {
    return String(element && element.getAttribute && element.getAttribute(name) || '').trim();
  }

  function getElementMatchId(element) {
    return matchAttr(element, 'data-gbdraw-pairwise-match-id');
  }

  function getMatchHoverKey(element) {
    var blockId = matchAttr(element, 'data-collinearity-block-id');
    if (blockId) return 'collinearity:' + blockId;
    var orthogroupId = matchAttr(element, 'data-orthogroup-id');
    if (orthogroupId) return 'orthogroup:' + orthogroupId;
    return 'match:' + (getElementMatchId(element) || matchAttr(element, 'd'));
  }

  function setMatchHover(element, highlight) {
    if (!element) return;
    setClassToken(element, 'gbdraw-interactive-pairwise-match--hover', highlight);
    var blockId = matchAttr(element, 'data-collinearity-block-id');
    if (blockId) {
      setCollinearityBlockHover(blockId, highlight);
      return;
    }
    var orthogroupId = matchAttr(element, 'data-orthogroup-id');
    if (orthogroupId) {
      setOrthogroupHover(orthogroupId, highlight);
    }
  }

  function clearActiveFeatureHover() {
    if (!activeHoverSvgId) return;
    setHoverHighlight(activeHoverSvgId, false);
    activeHoverSvgId = null;
    activeHoverKey = '';
  }

  function clearActiveMatchHover() {
    if (!activeHoverMatchId) return;
    var elements = matchElementsById.get(activeHoverMatchId) || [];
    elements.forEach(function (element) {
      setMatchHover(element, false);
    });
    activeHoverMatchId = '';
    activeHoverMatchKey = '';
  }

  function supportsStandaloneControls() {
    if (!document.createElementNS) return false;
    if (typeof window.SVGForeignObjectElement === 'undefined') return false;
    try {
      return document.createElementNS(SVG_NS, 'foreignObject') instanceof window.SVGForeignObjectElement;
    } catch (error) {
      return false;
    }
  }

  function normalizeSearchText(value) {
    return String(value == null ? '' : value).toLowerCase();
  }

  var NUCLEOTIDE_IUPAC = {
    A: 'A',
    C: 'C',
    G: 'G',
    T: 'TU',
    U: 'TU',
    R: 'AG',
    Y: 'CTU',
    S: 'GC',
    W: 'ATU',
    K: 'GTU',
    M: 'AC',
    B: 'CGTU',
    D: 'AGTU',
    H: 'ACTU',
    V: 'ACG',
    N: 'ACGTU',
    '-': '-'
  };
  var AMINO_ACID_IUPAC = {
    A: 'A',
    C: 'C',
    D: 'D',
    E: 'E',
    F: 'F',
    G: 'G',
    H: 'H',
    I: 'I',
    K: 'K',
    L: 'L',
    M: 'M',
    N: 'N',
    P: 'P',
    Q: 'Q',
    R: 'R',
    S: 'S',
    T: 'T',
    V: 'V',
    W: 'W',
    Y: 'Y',
    U: 'U',
    O: 'O',
    B: 'DN',
    Z: 'EQ',
    J: 'IL',
    X: 'ACDEFGHIKLMNPQRSTVWYUO',
    '*': '*',
    '-': '-'
  };

  function getIupacAlphabet(alphabet) {
    if (alphabet === 'nucleotide') return NUCLEOTIDE_IUPAC;
    if (alphabet === 'amino-acid') return AMINO_ACID_IUPAC;
    return null;
  }

  function uniqueCharacters(value) {
    var seen = {};
    var chars = [];
    String(value || '').split('').forEach(function (char) {
      if (!char || seen[char]) return;
      seen[char] = true;
      chars.push(char);
    });
    return chars;
  }

  function charSetsIntersect(left, right) {
    var rightSet = {};
    uniqueCharacters(right).forEach(function (char) {
      rightSet[char] = true;
    });
    return uniqueCharacters(left).some(function (char) {
      return Boolean(rightSet[char]);
    });
  }

  function escapeRegExpText(value) {
    return String(value).replace(/[.*+?^\${}()|[\]\\]/g, '\\\\$&');
  }

  function escapeRegExpClassChar(value) {
    return String(value).replace(/[\\\]\[\^-]/g, '\\\\$&');
  }

  function buildIupacQueryPattern(query, alphabet) {
    var map = getIupacAlphabet(alphabet);
    var normalized = String(query || '').replace(/\s+/g, '').toUpperCase();
    if (!map || !normalized) return null;
    var targetChars = Object.keys(map);
    var parts = [];
    normalized.split('').forEach(function (queryChar) {
      var querySet = map[queryChar] || queryChar;
      var matchingTargets = targetChars.filter(function (targetChar) {
        return charSetsIntersect(querySet, map[targetChar] || targetChar);
      });
      if (!matchingTargets.length) {
        parts.push(escapeRegExpText(queryChar));
      } else if (matchingTargets.length === 1) {
        parts.push(escapeRegExpText(matchingTargets[0]));
      } else {
        parts.push('[' + matchingTargets.map(escapeRegExpClassChar).join('') + ']');
      }
    });
    return parts.join('');
  }

  function appendSearchValues(values, value) {
    if (Array.isArray(value)) {
      value.forEach(function (entry) {
        appendSearchValues(values, entry);
      });
      return;
    }
    if (value === null || value === undefined) return;
    var text = String(value).trim();
    if (text) values.push(text);
  }

  function appendSearchItems(items, label, value, options) {
    if (Array.isArray(value)) {
      value.forEach(function (entry) {
        appendSearchItems(items, label, entry, options);
      });
      return;
    }
    if (value === null || value === undefined) return;
    var text = String(value).trim();
    if (text) {
      items.push({
        label: label,
        value: text,
        alphabet: options && options.alphabet ? String(options.alphabet) : ''
      });
    }
  }

  function getFeatureQualifiers(feature) {
    return feature && feature.qualifiers && typeof feature.qualifiers === 'object' && !Array.isArray(feature.qualifiers)
      ? feature.qualifiers
      : {};
  }

  function getQualifierSearchItems(feature, qualifierKey) {
    var target = normalizeSearchText(qualifierKey).trim();
    var qualifiers = getFeatureQualifiers(feature);
    var items = [];
    Object.keys(qualifiers).sort().forEach(function (key) {
      if (target && normalizeSearchText(key) !== target) return;
      appendSearchItems(items, 'Qualifier ' + key, qualifiers[key]);
    });
    return items;
  }

  function firstQualifierValue(feature, key) {
    var qualifiers = getFeatureQualifiers(feature);
    var values = [];
    appendSearchValues(values, qualifiers[key]);
    for (var i = 0; i < values.length; i += 1) {
      if (String(values[i] || '').trim()) return values[i];
    }
    return '';
  }

  function getFeatureLabelCandidates(feature) {
    return [
      feature && feature.display_label,
      feature && feature.displayLabel,
      feature && feature.label,
      feature && feature.gene,
      feature && feature.locus_tag,
      feature && feature.locusTag,
      feature && feature.product,
      firstQualifierValue(feature, 'gene'),
      firstQualifierValue(feature, 'locus_tag'),
      firstQualifierValue(feature, 'product'),
      feature && feature.search_labels,
      feature && feature.searchLabels,
      feature && feature.svg_id
    ];
  }

  function getLabelSearchItems(feature) {
    var items = [];
    var seen = {};
    getFeatureLabelCandidates(feature).forEach(function (value) {
      var values = [];
      appendSearchValues(values, value);
      values.forEach(function (entry) {
        var key = normalizeSearchText(entry);
        if (!key || seen[key]) return;
        seen[key] = true;
        appendSearchItems(items, key === normalizeSearchText(feature && feature.svg_id) ? 'SVG ID' : 'Label', entry);
      });
    });
    return items;
  }

  function buildFeatureLocation(feature) {
    var direct = String(feature && feature.location || '').trim();
    if (direct) return direct;
    var parts = Array.isArray(feature && feature.location_parts) ? feature.location_parts : [];
    var partText = parts.map(function (part) {
      return String(part && part.display || '').trim();
    }).filter(Boolean).join(', ');
    if (partText) return partText;
    var start = Number(feature && feature.start);
    var end = Number(feature && feature.end);
    var startText = Number.isFinite(start) ? String(start + 1) : String(feature && feature.start != null ? feature.start : '');
    var endText = Number.isFinite(end) ? String(end) : String(feature && feature.end != null ? feature.end : '');
    var range = startText + '..' + endText;
    var strand = String(feature && feature.strand || '').trim();
    return strand ? range + ' (' + strand + ')' : range;
  }

  function getOrthogroupById(orthogroupId) {
    var id = String(orthogroupId || '').trim();
    if (!id) return null;
    return orthogroupsById.get(id) || null;
  }

  function getFeatureOrthogroup(feature) {
    return getOrthogroupById(feature && (feature.orthogroup_id || feature.orthogroupId));
  }

  function getFeatureOrthogroupMember(feature, group) {
    if (feature && feature.orthogroup_member && typeof feature.orthogroup_member === 'object') {
      return feature.orthogroup_member;
    }
    if (feature && feature.orthogroupMember && typeof feature.orthogroupMember === 'object') {
      return feature.orthogroupMember;
    }
    var members = group && Array.isArray(group.members) ? group.members : [];
    var svgId = String(feature && feature.svg_id || '').trim();
    var recordIndex = Number(feature && (feature.record_idx || feature.recordIndex || feature.fileIdx));
    for (var i = 0; i < members.length; i += 1) {
      var member = members[i] || {};
      if (String(member.featureSvgId || member.feature_svg_id || '').trim() !== svgId) continue;
      if (!Number.isInteger(recordIndex) || Number(member.recordIndex || member.record_index) === recordIndex) {
        return member;
      }
    }
    return null;
  }

  function firstDisplayText() {
    for (var i = 0; i < arguments.length; i += 1) {
      var value = arguments[i];
      if (Array.isArray(value)) {
        var nested = firstDisplayText.apply(null, value);
        if (nested) return nested;
        continue;
      }
      var text = String(value == null ? '' : value).trim();
      if (text) return text;
    }
    return '';
  }

  function displayProteinId(feature, member, fallback) {
    return firstDisplayText(
      feature && (feature.source_protein_id || feature.sourceProteinId),
      member && (member.sourceProteinId || member.source_protein_id),
      firstQualifierValue(feature, 'protein_id'),
      feature && (feature.locus_tag || feature.locusTag),
      firstQualifierValue(feature, 'locus_tag'),
      member && (member.locusTag || member.locus_tag),
      feature && (feature.gene_id || feature.geneId),
      firstQualifierValue(feature, 'gene_id'),
      member && (member.geneId || member.gene_id),
      feature && (feature.old_locus_tag || feature.oldLocusTag),
      firstQualifierValue(feature, 'old_locus_tag'),
      member && (member.oldLocusTag || member.old_locus_tag),
      feature && feature.ID,
      firstQualifierValue(feature, 'ID'),
      feature && feature.Name,
      firstQualifierValue(feature, 'Name'),
      feature && feature.Parent,
      firstQualifierValue(feature, 'Parent'),
      feature && feature.gene,
      firstQualifierValue(feature, 'gene'),
      member && member.gene,
      feature && (feature.proteinId || feature.protein_id),
      member && (member.proteinId || member.protein_id),
      fallback
    );
  }

  function internalProteinId(feature, member) {
    return String(
      feature && feature.protein_id ||
      feature && feature.proteinId ||
      member && member.proteinId ||
      member && member.protein_id ||
      ''
    ).trim();
  }

  function getOrthogroupSearchItems(feature) {
    var orthogroupId = String(feature && (feature.orthogroup_id || feature.orthogroupId) || '').trim();
    var group = getFeatureOrthogroup(feature);
    var member = getFeatureOrthogroupMember(feature, group);
    var proteinId = displayProteinId(feature, member);
    var internalId = internalProteinId(feature, member);
    var items = [];
    appendSearchItems(items, 'Orthogroup ID', orthogroupId);
    appendSearchItems(items, 'Orthogroup name', group && (group.display_name || group.displayName || group.name));
    appendSearchItems(items, 'Orthogroup description', group && group.description);
    appendSearchItems(items, 'Protein ID', proteinId);
    if (internalId && internalId !== proteinId) {
      appendSearchItems(items, 'Internal protein ID', internalId);
    }
    appendSearchItems(items, 'Orthogroup member', member && member.label);
    appendSearchItems(items, 'Orthogroup member gene', member && member.gene);
    appendSearchItems(items, 'Orthogroup member product', member && member.product);
    appendSearchItems(items, 'Orthogroup member note', member && member.note);
    appendSearchItems(items, 'Orthogroup member protein ID', displayProteinId(null, member));
    return items;
  }

  function featureSearchItems(feature, field, qualifierKey) {
    var items = [];
    var selectedField = normalizeSearchField(field);
    var qualifiers = getFeatureQualifiers(feature);
    if (selectedField === 'label') {
      return getLabelSearchItems(feature);
    }
    if (selectedField === 'type') {
      appendSearchItems(items, 'Feature type', feature && feature.type);
      return items;
    }
    if (selectedField === 'record-id') {
      appendSearchItems(items, 'Record ID', feature && feature.record_id);
      appendSearchItems(items, 'Record ID', feature && feature.recordId);
      appendSearchItems(items, 'Record ID', feature && feature.displayRecordId);
      return items;
    }
    if (selectedField === 'location') {
      appendSearchItems(items, 'Location', buildFeatureLocation(feature));
      appendSearchItems(items, 'Start', feature && feature.start);
      appendSearchItems(items, 'End', feature && feature.end);
      return items;
    }
    if (selectedField === 'strand') {
      appendSearchItems(items, 'Strand', feature && feature.strand);
      return items;
    }
    if (selectedField === 'orthogroup') {
      return getOrthogroupSearchItems(feature);
    }
    if (selectedField === 'qualifier-key') {
      appendSearchItems(items, 'Qualifier key', Object.keys(qualifiers));
      return items;
    }
    if (selectedField === 'qualifier-value') {
      return getQualifierSearchItems(feature, qualifierKey);
    }
    if (selectedField === 'nucleotide') {
      appendSearchItems(items, 'Nucleotide sequence', feature && (feature.nucleotide_sequence || feature.nucleotideSequence), { alphabet: 'nucleotide' });
      return items;
    }
    if (selectedField === 'amino-acid') {
      appendSearchItems(items, 'Amino acid sequence', feature && (feature.amino_acid_sequence || feature.aminoAcidSequence), { alphabet: 'amino-acid' });
      return items;
    }

    Array.prototype.push.apply(items, getLabelSearchItems(feature));
    appendSearchItems(items, 'Record ID', feature && feature.record_id);
    appendSearchItems(items, 'Record ID', feature && feature.recordId);
    appendSearchItems(items, 'Record ID', feature && feature.displayRecordId);
    appendSearchItems(items, 'Feature type', feature && feature.type);
    appendSearchItems(items, 'Location', buildFeatureLocation(feature));
    appendSearchItems(items, 'Strand', feature && feature.strand);
    Array.prototype.push.apply(items, getOrthogroupSearchItems(feature));
    if (popupMode !== 'simple') {
      appendSearchItems(items, 'Qualifier key', Object.keys(qualifiers));
      Object.keys(qualifiers).sort().forEach(function (key) {
        appendSearchItems(items, 'Qualifier ' + key, qualifiers[key]);
      });
      appendSearchItems(items, 'Nucleotide sequence', feature && (feature.nucleotide_sequence || feature.nucleotideSequence), { alphabet: 'nucleotide' });
      appendSearchItems(items, 'Amino acid sequence', feature && (feature.amino_acid_sequence || feature.aminoAcidSequence), { alphabet: 'amino-acid' });
    }
    return items;
  }

  function featureSearchValues(feature, field, qualifierKey) {
    return featureSearchItems(feature, field, qualifierKey).map(function (item) {
      return item.value;
    });
  }

  function compileSearchMatcher(query, useRegex) {
    var trimmedQuery = String(query || '').trim();
    if (!trimmedQuery) {
      return { active: false, error: '', test: function () { return false; } };
    }
    if (useRegex) {
      try {
        var regex = new RegExp(trimmedQuery, 'i');
        return {
          active: true,
          error: '',
          match: function (value) {
            regex.lastIndex = 0;
            var match = String(value == null ? '' : value).match(regex);
            return match ? String(match[0] || '') : '';
          },
          test: function (values) {
            return values.some(function (value) {
              regex.lastIndex = 0;
              return regex.test(String(value == null ? '' : value));
            });
          }
        };
      } catch (error) {
        return { active: true, error: 'Invalid regex', test: function () { return false; } };
      }
    }
    var needle = normalizeSearchText(trimmedQuery);
    var sequenceMatchers = {};
    var getSequenceRegex = function (alphabet) {
      if (!getIupacAlphabet(alphabet)) return null;
      if (!Object.prototype.hasOwnProperty.call(sequenceMatchers, alphabet)) {
        var pattern = buildIupacQueryPattern(trimmedQuery, alphabet);
        sequenceMatchers[alphabet] = pattern ? new RegExp(pattern, 'i') : null;
      }
      return sequenceMatchers[alphabet];
    };
    return {
      active: true,
      error: '',
      match: function (value, alphabet) {
        var text = String(value == null ? '' : value);
        var sequenceRegex = getSequenceRegex(alphabet);
        if (sequenceRegex) {
          sequenceRegex.lastIndex = 0;
          var sequenceMatch = text.match(sequenceRegex);
          return sequenceMatch ? String(sequenceMatch[0] || '') : '';
        }
        var index = normalizeSearchText(text).indexOf(needle);
        return index === -1 ? '' : text.slice(index, index + trimmedQuery.length);
      },
      test: function (values) {
        return values.some(function (value) {
          return normalizeSearchText(value).indexOf(needle) !== -1;
        });
      }
    };
  }

  function featureMatchesSearch(feature, matcher, field, qualifierKey) {
    if (!matcher || !matcher.active || matcher.error) return false;
    return featureSearchMatches(feature, matcher, field, qualifierKey).length > 0;
  }

  function featureSearchMatches(feature, matcher, field, qualifierKey) {
    if (!matcher || !matcher.active || matcher.error) return [];
    return featureSearchItems(feature, field, qualifierKey).map(function (item) {
      var matchedText = matcher.match ? matcher.match(item.value, item.alphabet) : '';
      if (!matchedText) return null;
      return {
        label: item.label,
        value: String(item.value),
        match: matchedText
      };
    }).filter(Boolean);
  }

  function collapseWhitespace(value) {
    return String(value == null ? '' : value).replace(/\\s+/g, ' ').trim();
  }

  function searchMatchSnippet(value, matchText) {
    var text = collapseWhitespace(value);
    var match = collapseWhitespace(matchText);
    if (!text) return '';
    if (!match) return text.length > 80 ? text.slice(0, 77) + '...' : text;
    var lowerText = text.toLowerCase();
    var index = lowerText.indexOf(match.toLowerCase());
    if (index === -1) return text.length > 80 ? text.slice(0, 77) + '...' : text;
    var start = Math.max(0, index - 24);
    var end = Math.min(text.length, index + match.length + 24);
    return (start > 0 ? '...' : '') + text.slice(start, end) + (end < text.length ? '...' : '');
  }

  function formatSearchMatchDetail(detail) {
    if (!detail) return '';
    var label = collapseWhitespace(detail.label);
    var snippet = searchMatchSnippet(detail.value, detail.match);
    if (!label && !snippet) return '';
    if (!snippet) return label;
    return label ? label + ': ' + snippet : snippet;
  }

  function syncSearchControls() {
    if (!searchControls) return;
    var root = searchControls.querySelector('.gfs');
    var queryInput = searchControls.querySelector('[data-search-query]');
    var fieldSelect = searchControls.querySelector('[data-search-field]');
    var qualifierInput = searchControls.querySelector('[data-search-qualifier]');
    var regexInput = searchControls.querySelector('[data-search-regex]');
    var countText = searchControls.querySelector('[data-search-count]');
    var matchDetailText = searchControls.querySelector('[data-search-match-detail]');
    var prevButton = searchControls.querySelector('[data-search-prev]');
    var nextButton = searchControls.querySelector('[data-search-next]');
    var openButton = searchControls.querySelector('[data-search-open]');
    var clearButton = searchControls.querySelector('[data-search-clear]');
    searchState.field = normalizeSearchField(searchState.field);
    pendingSearchState.field = normalizeSearchField(pendingSearchState.field);
    if (queryInput && queryInput.value !== pendingSearchState.query) queryInput.value = pendingSearchState.query;
    if (fieldSelect && fieldSelect.value !== pendingSearchState.field) fieldSelect.value = pendingSearchState.field;
    if (qualifierInput && qualifierInput.value !== pendingSearchState.qualifierKey) qualifierInput.value = pendingSearchState.qualifierKey;
    if (regexInput) regexInput.checked = Boolean(pendingSearchState.useRegex);
    if (qualifierInput) {
      qualifierInput.disabled = popupMode === 'simple' || pendingSearchState.field !== 'qualifier-value';
      qualifierInput.style.display = popupMode === 'simple' ? 'none' : '';
    }
    var hasMatches = searchState.matches.length > 0;
    if (prevButton) prevButton.disabled = !hasMatches;
    if (nextButton) nextButton.disabled = !hasMatches;
    if (openButton) openButton.disabled = !hasMatches;
    if (clearButton) {
      clearButton.disabled = !pendingSearchState.query &&
        !pendingSearchState.qualifierKey &&
        pendingSearchState.field === 'all' &&
        !pendingSearchState.useRegex &&
        !searchState.query &&
        !searchState.qualifierKey &&
        !searchState.error &&
        !searchState.matches.length;
    }
    setClassToken(root, 'is-invalid', Boolean(searchState.error));
    if (countText) {
      if (searchState.error) {
        countText.textContent = searchState.error;
      } else if (!String(searchState.query || '').trim()) {
        countText.textContent = '0 / ' + String(featuresById.size) + ' features';
      } else {
        var current = searchState.activeIndex >= 0 ? searchState.activeIndex + 1 : 0;
        countText.textContent = String(current) + ' / ' + String(searchState.matches.length) + ' features';
      }
    }
    if (matchDetailText) {
      var activeId = searchState.activeIndex >= 0 ? searchState.matches[searchState.activeIndex] : '';
      var details = activeId ? searchState.matchDetails[activeId] || [] : [];
      var detailText = details.length ? formatSearchMatchDetail(details[0]) : '';
      matchDetailText.textContent = detailText ? 'Matched ' + detailText : '';
      matchDetailText.setAttribute('title', detailText ? 'Matched ' + detailText : '');
    }
  }

  function applySearchResults() {
    var queryActive = Boolean(String(searchState.query || '').trim()) && !searchState.error;
    var matchedIds = new Set(searchState.matches);
    var activeId = searchState.activeIndex >= 0 ? searchState.matches[searchState.activeIndex] : '';
    featureElementsById.forEach(function (elements, svgId) {
      var isMatch = queryActive && matchedIds.has(svgId);
      var isActive = isMatch && svgId === activeId;
      elements.forEach(function (element) {
        setClassToken(element, 'gbdraw-interactive-feature--match', isMatch);
        setClassToken(element, 'gbdraw-interactive-feature--active-match', isActive);
        setClassToken(element, 'gbdraw-interactive-feature--dimmed', queryActive && !isMatch);
      });
    });
    syncSearchControls();
  }

  function setSearchState(nextState) {
    nextState = nextState || {};
    var previousActiveId = searchState.activeIndex >= 0 ? searchState.matches[searchState.activeIndex] : '';
    searchState.query = Object.prototype.hasOwnProperty.call(nextState, 'query') ? String(nextState.query || '') : searchState.query;
    searchState.field = Object.prototype.hasOwnProperty.call(nextState, 'field')
      ? normalizeSearchField(nextState.field)
      : normalizeSearchField(searchState.field);
    searchState.qualifierKey = Object.prototype.hasOwnProperty.call(nextState, 'qualifierKey') ? String(nextState.qualifierKey || '') : searchState.qualifierKey;
    searchState.useRegex = Object.prototype.hasOwnProperty.call(nextState, 'useRegex') ? Boolean(nextState.useRegex) : searchState.useRegex;
    pendingSearchState.query = searchState.query;
    pendingSearchState.field = normalizeSearchField(searchState.field);
    pendingSearchState.qualifierKey = searchState.qualifierKey;
    pendingSearchState.useRegex = Boolean(searchState.useRegex);

    var matcher = compileSearchMatcher(searchState.query, searchState.useRegex);
    searchState.error = matcher.error;
    if (!matcher.active || matcher.error) {
      searchState.matches = [];
      searchState.matchDetails = {};
      searchState.activeIndex = -1;
      applySearchResults();
      return;
    }

    var nextMatchDetails = {};
    searchState.matches = features.filter(function (feature) {
      var svgId = String(feature && feature.svg_id || '').trim();
      if (!svgId || !featureElementsById.has(svgId)) return false;
      var details = featureSearchMatches(feature, matcher, searchState.field, searchState.qualifierKey);
      if (!details.length) return false;
      nextMatchDetails[svgId] = details;
      return true;
    }).map(function (feature) {
      return String(feature.svg_id);
    });
    searchState.matchDetails = nextMatchDetails;

    if (!searchState.matches.length) {
      searchState.activeIndex = -1;
    } else if (previousActiveId && searchState.matches.indexOf(previousActiveId) !== -1) {
      searchState.activeIndex = searchState.matches.indexOf(previousActiveId);
    } else {
      searchState.activeIndex = 0;
    }
    applySearchResults();
  }

  function clearAppliedSearchResults() {
    var hadAppliedSearch = Boolean(String(searchState.query || '').trim()) ||
      Boolean(searchState.qualifierKey) ||
      Boolean(searchState.error) ||
      searchState.matches.length > 0;
    searchState.query = '';
    searchState.field = 'all';
    searchState.qualifierKey = '';
    searchState.useRegex = false;
    searchState.matches = [];
    searchState.matchDetails = {};
    searchState.activeIndex = -1;
    searchState.error = '';
    if (hadAppliedSearch) {
      applySearchResults();
    } else {
      syncSearchControls();
    }
  }

  function setPendingSearchState(nextState) {
    nextState = nextState || {};
    pendingSearchState.query = Object.prototype.hasOwnProperty.call(nextState, 'query') ? String(nextState.query || '') : pendingSearchState.query;
    pendingSearchState.field = Object.prototype.hasOwnProperty.call(nextState, 'field')
      ? normalizeSearchField(nextState.field)
      : normalizeSearchField(pendingSearchState.field);
    pendingSearchState.qualifierKey = Object.prototype.hasOwnProperty.call(nextState, 'qualifierKey') ? String(nextState.qualifierKey || '') : pendingSearchState.qualifierKey;
    pendingSearchState.useRegex = Object.prototype.hasOwnProperty.call(nextState, 'useRegex') ? Boolean(nextState.useRegex) : pendingSearchState.useRegex;
    clearAppliedSearchResults();
  }

  function clearSearch() {
    pendingSearchState.query = '';
    pendingSearchState.field = 'all';
    pendingSearchState.qualifierKey = '';
    pendingSearchState.useRegex = false;
    searchState.query = '';
    searchState.field = 'all';
    searchState.qualifierKey = '';
    searchState.useRegex = false;
    searchState.matches = [];
    searchState.matchDetails = {};
    searchState.activeIndex = -1;
    searchState.error = '';
    applySearchResults();
  }

  function getFeatureBounds(svgId) {
    var elements = featureElementsById.get(String(svgId || '')) || [];
    var minX = Infinity;
    var minY = Infinity;
    var maxX = -Infinity;
    var maxY = -Infinity;

    function addPoint(point) {
      if (!point || !Number.isFinite(point.x) || !Number.isFinite(point.y)) return;
      minX = Math.min(minX, point.x);
      minY = Math.min(minY, point.y);
      maxX = Math.max(maxX, point.x);
      maxY = Math.max(maxY, point.y);
    }

    elements.forEach(function (element) {
      var rect = typeof element.getBoundingClientRect === 'function' ? element.getBoundingClientRect() : null;
      if (rect && Number.isFinite(rect.left) && Number.isFinite(rect.right) && rect.right > rect.left && rect.bottom > rect.top) {
        addPoint(clientPoint(rect.left, rect.top));
        addPoint(clientPoint(rect.right, rect.top));
        addPoint(clientPoint(rect.right, rect.bottom));
        addPoint(clientPoint(rect.left, rect.bottom));
        return;
      }
      if (typeof element.getBBox !== 'function') return;
      try {
        var bbox = element.getBBox();
        addPoint({ x: bbox.x, y: bbox.y });
        addPoint({ x: bbox.x + bbox.width, y: bbox.y + bbox.height });
      } catch (error) {
        // Ignore elements that cannot report a box in this viewer.
      }
    });

    if (!Number.isFinite(minX) || !Number.isFinite(minY) || maxX <= minX || maxY <= minY) {
      return null;
    }
    return { x: minX, y: minY, width: maxX - minX, height: maxY - minY };
  }

  function centerFeatureInView(svgId) {
    var bounds = getFeatureBounds(svgId);
    if (!bounds) return;
    var view = getViewRect();
    var targetWidth = Math.min(view.width, Math.max(bounds.width * 6, originalViewRect.width / 8));
    var targetHeight = Math.min(view.height, Math.max(bounds.height * 6, originalViewRect.height / 8));
    targetWidth = Math.max(targetWidth, bounds.width * 1.4);
    targetHeight = Math.max(targetHeight, bounds.height * 1.4);
    var targetRect = fitRectToAspect({
      x: bounds.x + bounds.width / 2 - targetWidth / 2,
      y: bounds.y + bounds.height / 2 - targetHeight / 2,
      width: targetWidth,
      height: targetHeight
    }, homeViewRect.width / homeViewRect.height);
    setSvgViewRect(targetRect);
  }

  function setActiveMatch(index, options) {
    if (!searchState.matches.length) {
      searchState.activeIndex = -1;
      applySearchResults();
      return;
    }
    var count = searchState.matches.length;
    var nextIndex = ((Number(index) || 0) % count + count) % count;
    searchState.activeIndex = nextIndex;
    applySearchResults();
    if (options && options.center) {
      centerFeatureInView(searchState.matches[nextIndex]);
    }
    if (options && options.openPopup && popup) {
      openActiveMatchPopup({ center: false });
    }
  }

  function openActiveMatchPopup(options) {
    if (!searchState.matches.length) return;
    if (searchState.activeIndex < 0) {
      setActiveMatch(0, { center: true });
    }
    var svgId = searchState.matches[searchState.activeIndex];
    var feature = featuresById.get(String(svgId || ''));
    if (!feature) return;
    if (!options || options.center !== false) {
      centerFeatureInView(svgId);
    }
    openPopup(feature, null);
  }

  function createXhtmlNode(tagName, attrs) {
    var node = document.createElementNS(XHTML_NS, tagName);
    Object.keys(attrs || {}).forEach(function (key) {
      if (key === 'text') {
        node.textContent = attrs[key];
      } else if (key === 'className') {
        node.setAttribute('class', attrs[key]);
      } else {
        node.setAttribute(key, attrs[key]);
      }
    });
    return node;
  }

  function updateSearchControlsPosition() {
    if (!searchControls) return;
    var visibleView = getVisibleViewRect();
    var view = getViewRect();
    var screenScale = getScreenScale();
    var fallbackSize = getSvgClientSize();
    var unit = Math.max(
      1 / Math.max(screenScale.x, 0.001),
      1 / Math.max(screenScale.y, 0.001),
      view.width / fallbackSize.width,
      view.height / fallbackSize.height
    );
    var margin = 12 * unit;
    var controlWidth = Number(searchControls.getAttribute('width')) || 430;
    var x = Math.max(
      visibleView.x + margin,
      visibleView.x + visibleView.width - (controlWidth * unit) - margin
    );
    searchControls.setAttribute(
      'transform',
      'translate(' +
        formatSvgNumber(x + (searchControlsOffsetCss.x * unit)) +
        ', ' +
        formatSvgNumber(visibleView.y + margin + (searchControlsOffsetCss.y * unit)) +
        ') scale(' + formatSvgNumber(unit) + ')'
    );
  }

  function stopSearchControlsDrag() {
    if (!activeSearchControlsDrag) return;
    document.removeEventListener('mousemove', activeSearchControlsDrag.onMove, true);
    document.removeEventListener('mouseup', activeSearchControlsDrag.onEnd, true);
    window.removeEventListener('mouseup', activeSearchControlsDrag.onEnd, true);
    window.removeEventListener('blur', activeSearchControlsDrag.onEnd);
    if (activeSearchControlsDrag.root && activeSearchControlsDrag.root.classList) {
      activeSearchControlsDrag.root.classList.remove('is-dragging');
    }
    activeSearchControlsDrag = null;
  }

  function startSearchControlsDrag(event, root) {
    if (event.button !== 0) return;
    stopSearchControlsDrag();
    var startClientX = Number(event.clientX) || 0;
    var startClientY = Number(event.clientY) || 0;
    var startOffsetX = Number(searchControlsOffsetCss.x) || 0;
    var startOffsetY = Number(searchControlsOffsetCss.y) || 0;
    var onMove = function (moveEvent) {
      searchControlsOffsetCss.x = startOffsetX + ((Number(moveEvent.clientX) || 0) - startClientX);
      searchControlsOffsetCss.y = startOffsetY + ((Number(moveEvent.clientY) || 0) - startClientY);
      updateSearchControlsPosition();
      moveEvent.preventDefault();
    };
    var onEnd = function () {
      stopSearchControlsDrag();
    };
    activeSearchControlsDrag = {
      root: root,
      onMove: onMove,
      onEnd: onEnd
    };
    if (root && root.classList) root.classList.add('is-dragging');
    document.addEventListener('mousemove', onMove, true);
    document.addEventListener('mouseup', onEnd, true);
    window.addEventListener('mouseup', onEnd, true);
    window.addEventListener('blur', onEnd);
    event.preventDefault();
    event.stopPropagation();
  }

  function setupSearchControls() {
    var existing = svg.querySelector('#' + SEARCH_CONTROLS_ID);
    if (existing && existing.parentNode) {
      existing.parentNode.removeChild(existing);
    }
    if (!supportsStandaloneControls()) return;

    searchControls = document.createElementNS(SVG_NS, 'foreignObject');
    searchControls.setAttribute('id', SEARCH_CONTROLS_ID);
    searchControls.setAttribute('class', 'gbdraw-feature-search-controls');
    searchControls.setAttribute('data-popup-mode', popupMode);
    searchControls.setAttribute('width', '430');
    searchControls.setAttribute('height', '86');

    var root = createXhtmlNode('div', {
      xmlns: XHTML_NS,
      className: 'gfs'
    });
    var firstRow = createXhtmlNode('div', { className: 'gfs-row' });
    var secondRow = createXhtmlNode('div', { className: 'gfs-row' });
    var queryInput = createXhtmlNode('input', {
      className: 'gfs-input gfs-query',
      type: 'search',
      placeholder: 'Search features',
      'aria-label': 'Search features',
      'data-search-query': 'true'
    });
    var fieldSelect = createXhtmlNode('select', {
      className: 'gfs-select',
      'aria-label': 'Search field',
      'data-search-field': 'true'
    });
    var searchFieldOptions = [
      ['all', 'All'],
      ['label', 'Label'],
      ['type', 'Feature type'],
      ['record-id', 'Record ID'],
      ['location', 'Location'],
      ['strand', 'Strand'],
      ['orthogroup', 'Orthogroup'],
      ['qualifier-key', 'Qualifier key'],
      ['qualifier-value', 'Qualifier value'],
      ['nucleotide', 'Nucleotide'],
      ['amino-acid', 'Amino acid']
    ];
    if (popupMode === 'simple') {
      searchFieldOptions = searchFieldOptions.filter(function (entry) {
        return !isRichSearchField(entry[0]);
      });
    }
    searchFieldOptions.forEach(function (entry) {
      var option = createXhtmlNode('option', { value: entry[0], text: entry[1] });
      fieldSelect.appendChild(option);
    });
    var searchButton = createXhtmlNode('button', {
      type: 'button',
      className: 'gfs-button gfs-button--search',
      title: 'Search features',
      'aria-label': 'Search features',
      'data-search-apply': 'true',
      text: 'Search'
    });
    var qualifierInput = createXhtmlNode('input', {
      className: 'gfs-input gfs-qualifier',
      type: 'text',
      placeholder: 'Qualifier key',
      'aria-label': 'Qualifier key for qualifier value search',
      'data-search-qualifier': 'true'
    });
    var regexLabel = createXhtmlNode('label', { className: 'gfs-toggle' });
    var regexInput = createXhtmlNode('input', {
      type: 'checkbox',
      'data-search-regex': 'true'
    });
    regexLabel.appendChild(regexInput);
    regexLabel.appendChild(document.createTextNode('Regex'));
    var prevButton = createXhtmlNode('button', {
      type: 'button',
      className: 'gfs-button',
      title: 'Previous match',
      'aria-label': 'Previous match',
      'data-search-prev': 'true',
      text: '<'
    });
    var nextButton = createXhtmlNode('button', {
      type: 'button',
      className: 'gfs-button',
      title: 'Next match',
      'aria-label': 'Next match',
      'data-search-next': 'true',
      text: '>'
    });
    var openButton = createXhtmlNode('button', {
      type: 'button',
      className: 'gfs-button gfs-button--open',
      title: 'Open active feature',
      'aria-label': 'Open active feature',
      'data-search-open': 'true',
      text: 'Open'
    });
    var clearButton = createXhtmlNode('button', {
      type: 'button',
      className: 'gfs-button gfs-button--clear',
      title: 'Clear search',
      'aria-label': 'Clear search',
      'data-search-clear': 'true',
      text: 'Clear'
    });
    var countText = createXhtmlNode('span', {
      className: 'gfs-count',
      'data-search-count': 'true',
      text: '0 / ' + String(featuresById.size) + ' features'
    });
    var matchDetailText = createXhtmlNode('div', {
      className: 'gfs-match-detail',
      'data-search-match-detail': 'true'
    });

    firstRow.appendChild(queryInput);
    firstRow.appendChild(fieldSelect);
    firstRow.appendChild(searchButton);
    if (popupMode !== 'simple') {
      secondRow.appendChild(qualifierInput);
    }
    secondRow.appendChild(regexLabel);
    secondRow.appendChild(prevButton);
    secondRow.appendChild(nextButton);
    secondRow.appendChild(openButton);
    secondRow.appendChild(clearButton);
    secondRow.appendChild(countText);
    root.appendChild(firstRow);
    root.appendChild(secondRow);
    root.appendChild(matchDetailText);
    searchControls.appendChild(root);

    root.addEventListener('mousedown', function (event) {
      var target = event.target && event.target.closest
        ? event.target.closest('input, textarea, select, button, label, a')
        : null;
      if (target && root.contains(target)) return;
      startSearchControlsDrag(event, root);
    });

    ['mousedown', 'mouseup', 'click', 'dblclick', 'keydown', 'keyup', 'keypress', 'wheel', 'touchstart', 'touchmove'].forEach(function (eventName) {
      root.addEventListener(eventName, function (event) {
        event.stopPropagation();
      }, { passive: eventName === 'touchstart' || eventName === 'touchmove' });
    });
    queryInput.addEventListener('input', function () {
      setPendingSearchState({ query: queryInput.value });
    });
    fieldSelect.addEventListener('change', function () {
      setPendingSearchState({ field: fieldSelect.value });
    });
    qualifierInput.addEventListener('input', function () {
      setPendingSearchState({ qualifierKey: qualifierInput.value });
    });
    regexInput.addEventListener('change', function () {
      setPendingSearchState({ useRegex: regexInput.checked });
    });
    searchButton.addEventListener('click', function () {
      setSearchState({
        query: pendingSearchState.query,
        field: pendingSearchState.field,
        qualifierKey: pendingSearchState.qualifierKey,
        useRegex: pendingSearchState.useRegex
      });
      queryInput.focus();
    });
    prevButton.addEventListener('click', function () {
      setActiveMatch(searchState.activeIndex - 1, { center: true, openPopup: Boolean(popup) });
    });
    nextButton.addEventListener('click', function () {
      setActiveMatch(searchState.activeIndex + 1, { center: true, openPopup: Boolean(popup) });
    });
    openButton.addEventListener('click', function () {
      openActiveMatchPopup();
    });
    clearButton.addEventListener('click', function () {
      clearSearch();
      queryInput.focus();
    });
    root.addEventListener('keydown', function (event) {
      if (event.key === 'Escape') {
        clearSearch();
        queryInput.focus();
        event.preventDefault();
      } else if (event.key === 'Enter' && searchState.matches.length) {
        openActiveMatchPopup();
        event.preventDefault();
      }
    });

    svg.appendChild(searchControls);
    syncStandaloneOverlayOrder();
    syncSearchControls();
    updateSearchControlsPosition();
  }

  function closestSearchControls(target) {
    return Boolean(searchControls && target && searchControls.contains(target));
  }

  function formatSvgNumber(value) {
    var numeric = Number(value);
    if (!Number.isFinite(numeric)) return '0';
    return String(Math.round(numeric * 1000) / 1000);
  }

  function copyViewRect(rect) {
    return {
      x: Number(rect && rect.x) || 0,
      y: Number(rect && rect.y) || 0,
      width: Math.max(1, Number(rect && rect.width) || 1),
      height: Math.max(1, Number(rect && rect.height) || 1)
    };
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

  function parseViewBoxRect(value) {
    var parts = String(value || '').trim().split(/[\\s,]+/).map(function (part) {
      return Number(part);
    });
    if (parts.length < 4 || parts.some(function (part) { return !Number.isFinite(part); })) {
      return null;
    }
    if (parts[2] <= 0 || parts[3] <= 0) return null;
    return { x: parts[0], y: parts[1], width: parts[2], height: parts[3] };
  }

  function parseOriginalViewRectFromSvg() {
    return parseViewBoxRect(svg.getAttribute('data-gbdraw-original-viewbox')) ||
      parseViewBoxRect(svg.getAttribute('viewBox')) ||
      copyViewRect(getViewRect());
  }

  function getViewportAspect() {
    var doc = document.documentElement || {};
    var width = window.innerWidth || doc.clientWidth || 0;
    var height = window.innerHeight || doc.clientHeight || 0;
    if (!Number.isFinite(width) || !Number.isFinite(height) || width <= 0 || height <= 0) {
      return 0;
    }
    return width / height;
  }

  function fitRectToAspect(rect, targetAspect) {
    var source = copyViewRect(rect);
    var sourceAspect = source.width / source.height;
    if (!Number.isFinite(targetAspect) || targetAspect <= 0 || !Number.isFinite(sourceAspect) || sourceAspect <= 0) {
      return source;
    }
    if (targetAspect > sourceAspect) {
      var width = source.height * targetAspect;
      return {
        x: source.x + source.width / 2 - width / 2,
        y: source.y,
        width: width,
        height: source.height
      };
    }
    var height = source.width / targetAspect;
    return {
      x: source.x,
      y: source.y + source.height / 2 - height / 2,
      width: source.width,
      height: height
    };
  }

  function rectsNearlyEqual(a, b) {
    return Math.abs(a.x - b.x) < 0.5 &&
      Math.abs(a.y - b.y) < 0.5 &&
      Math.abs(a.width - b.width) < 0.5 &&
      Math.abs(a.height - b.height) < 0.5;
  }

  var originalViewRect = parseOriginalViewRectFromSvg();
  var homeViewRect = fitRectToAspect(originalViewRect, getViewportAspect());

  function isElementHidden(element) {
    if (!element) return true;
    if (String(element.getAttribute('display') || '').trim().toLowerCase() === 'none') return true;
    if (element.style && String(element.style.display || '').trim().toLowerCase() === 'none') return true;
    if (typeof window.getComputedStyle === 'function') {
      try {
        var style = window.getComputedStyle(element);
        if (style && (style.display === 'none' || style.visibility === 'hidden')) return true;
      } catch (error) {
        // Ignore style lookup failures and let getBBox decide renderability.
      }
    }
    return false;
  }

  function getElementBBox(element) {
    if (!element || typeof element.getBBox !== 'function') return null;
    try {
      var bbox = element.getBBox();
      if (!bbox || !Number.isFinite(bbox.x) || !Number.isFinite(bbox.y)) return null;
      if (!Number.isFinite(bbox.width) || !Number.isFinite(bbox.height)) return null;
      if (bbox.width <= 0 || bbox.height <= 0) return null;
      return {
        x: bbox.x,
        y: bbox.y,
        width: bbox.width,
        height: bbox.height
      };
    } catch (error) {
      return null;
    }
  }

  function getElementBoundsInSvg(element, bbox) {
    if (!element || !bbox) return null;
    var point = typeof svg.createSVGPoint === 'function' ? svg.createSVGPoint() : null;
    var elementMatrix = typeof element.getScreenCTM === 'function' ? element.getScreenCTM() : null;
    var svgMatrix = typeof svg.getScreenCTM === 'function' ? svg.getScreenCTM() : null;
    if (!point || !elementMatrix || !svgMatrix) {
      return { x: bbox.x, y: bbox.y, width: bbox.width, height: bbox.height };
    }
    var svgInverse = null;
    try {
      svgInverse = svgMatrix.inverse();
    } catch (error) {
      return { x: bbox.x, y: bbox.y, width: bbox.width, height: bbox.height };
    }

    function convert(x, y) {
      point.x = x;
      point.y = y;
      var screenPoint = point.matrixTransform(elementMatrix);
      return screenPoint.matrixTransform(svgInverse);
    }

    var points = [
      convert(bbox.x, bbox.y),
      convert(bbox.x + bbox.width, bbox.y),
      convert(bbox.x + bbox.width, bbox.y + bbox.height),
      convert(bbox.x, bbox.y + bbox.height)
    ].filter(function (entry) {
      return entry && Number.isFinite(entry.x) && Number.isFinite(entry.y);
    });
    if (points.length !== 4) return null;
    var minX = Math.min.apply(null, points.map(function (entry) { return entry.x; }));
    var maxX = Math.max.apply(null, points.map(function (entry) { return entry.x; }));
    var minY = Math.min.apply(null, points.map(function (entry) { return entry.y; }));
    var maxY = Math.max.apply(null, points.map(function (entry) { return entry.y; }));
    if (maxX <= minX || maxY <= minY) return null;
    return { x: minX, y: minY, width: maxX - minX, height: maxY - minY };
  }

  function inferStickyLegendAlignment(center, start, size, low, high) {
    if (!Number.isFinite(size) || size <= 0) return 'center';
    var ratio = (center - start) / size;
    if (ratio < low) return 'left';
    if (ratio > high) return 'right';
    return 'center';
  }

  function inferStickyLegendVerticalAlignment(center, start, size, low, high) {
    if (!Number.isFinite(size) || size <= 0) return 'center';
    var ratio = (center - start) / size;
    if (ratio < low) return 'top';
    if (ratio > high) return 'bottom';
    return 'center';
  }

  function getStickyLegendFallbackMargin() {
    var shortSide = Math.min(Number(homeViewRect.width) || 0, Number(homeViewRect.height) || 0);
    return Math.max(8, Math.min(18, shortSide * 0.018 || 12));
  }

  function updateStickyLegendHomeMetrics() {
    if (!stickyLegendState) return;
    var bounds = stickyLegendState.initialBounds;
    var centerX = bounds.x + bounds.width / 2;
    var centerY = bounds.y + bounds.height / 2;
    var fallbackMargin = getStickyLegendFallbackMargin();
    var alignmentView = originalViewRect || homeViewRect;
    stickyLegendState.anchorX = inferStickyLegendAlignment(
      centerX,
      alignmentView.x,
      alignmentView.width,
      1 / 3,
      2 / 3
    );
    stickyLegendState.anchorY = inferStickyLegendVerticalAlignment(
      centerY,
      alignmentView.y,
      alignmentView.height,
      1 / 3,
      2 / 3
    );
    stickyLegendState.centerOffsetX = centerX - (homeViewRect.x + homeViewRect.width / 2);
    stickyLegendState.centerOffsetY = centerY - (homeViewRect.y + homeViewRect.height / 2);
    stickyLegendState.baseMargin = fallbackMargin;
    stickyLegendState.margins = {
      left: fallbackMargin,
      right: fallbackMargin,
      top: fallbackMargin,
      bottom: fallbackMargin
    };
  }

  function getStickyLegendScale(view) {
    var scale = Number(view && view.width) / Number(homeViewRect.width);
    if (!Number.isFinite(scale) || scale <= 0) {
      scale = Number(view && view.height) / Number(homeViewRect.height);
    }
    if (!Number.isFinite(scale) || scale <= 0) return 1;
    return Math.max(scale, 0.000001);
  }

  function getStickyLegendMargin(axis, scale) {
    if (!stickyLegendState) return 0;
    var margins = stickyLegendState.margins || {};
    var margin = Number(margins[axis]);
    var base = Number(stickyLegendState.baseMargin) || 8;
    return Math.max(base, Number.isFinite(margin) ? margin : base) * scale;
  }

  function clampStickyLegendCoordinate(value, viewStart, viewSize, itemSize, margin) {
    if (itemSize + margin * 2 > viewSize) {
      return viewStart + (viewSize - itemSize) / 2;
    }
    var min = viewStart + margin;
    var max = viewStart + viewSize - itemSize - margin;
    return clampValue(value, min, max);
  }

  function ensureStickyLegendBackground(legend, bbox) {
    if (!legend || !bbox) return;
    var existing = legend.querySelector('#gbdraw-sticky-legend-background');
    if (existing && existing.parentNode) {
      existing.parentNode.removeChild(existing);
    }
    var padding = Math.max(8, Math.min(18, Math.min(bbox.width, bbox.height) * 0.08 || 10));
    var rect = createSvgNode('rect', {
      id: 'gbdraw-sticky-legend-background',
      'class': 'gbdraw-sticky-legend-background',
      x: formatSvgNumber(bbox.x - padding),
      y: formatSvgNumber(bbox.y - padding),
      width: formatSvgNumber(bbox.width + padding * 2),
      height: formatSvgNumber(bbox.height + padding * 2),
      rx: formatSvgNumber(Math.min(12, padding)),
      ry: formatSvgNumber(Math.min(12, padding)),
      fill: '#ffffff',
      'fill-opacity': '0.84',
      stroke: '#cbd5e1',
      'stroke-opacity': '0.7',
      'stroke-width': '1'
    });
    legend.insertBefore(rect, legend.firstChild || null);
  }

  function updateStickyLegendPosition() {
    if (!stickyLegend || !stickyLegendState) return;
    var visibleView = getVisibleViewRect();
    var bbox = stickyLegendState.localBBox;
    var scale = getStickyLegendScale(getViewRect());
    var width = bbox.width * scale;
    var height = bbox.height * scale;
    var marginX = stickyLegendState.baseMargin * scale;
    var marginY = stickyLegendState.baseMargin * scale;
    var targetX = visibleView.x + visibleView.width / 2 + stickyLegendState.centerOffsetX * scale - width / 2;
    var targetY = visibleView.y + visibleView.height / 2 + stickyLegendState.centerOffsetY * scale - height / 2;

    if (stickyLegendState.anchorX === 'left') {
      marginX = getStickyLegendMargin('left', scale);
      targetX = visibleView.x + marginX;
    } else if (stickyLegendState.anchorX === 'right') {
      marginX = getStickyLegendMargin('right', scale);
      targetX = visibleView.x + visibleView.width - width - marginX;
    }

    if (stickyLegendState.anchorY === 'top') {
      marginY = getStickyLegendMargin('top', scale);
      targetY = visibleView.y + marginY;
    } else if (stickyLegendState.anchorY === 'bottom') {
      marginY = getStickyLegendMargin('bottom', scale);
      targetY = visibleView.y + visibleView.height - height - marginY;
    }

    targetX = clampStickyLegendCoordinate(targetX, visibleView.x, visibleView.width, width, marginX);
    targetY = clampStickyLegendCoordinate(targetY, visibleView.y, visibleView.height, height, marginY);
    stickyLegend.setAttribute(
      'transform',
      'translate(' +
        formatSvgNumber(targetX - bbox.x * scale) +
        ', ' +
        formatSvgNumber(targetY - bbox.y * scale) +
        ') scale(' +
        formatSvgNumber(scale) +
        ')'
    );
  }

  function syncStandaloneOverlayOrder() {
    if (stickyLegend && stickyLegend.parentNode) {
      stickyLegend.parentNode.appendChild(stickyLegend);
    }
    if (hoverPopup && hoverPopup.parentNode === svg) {
      svg.appendChild(hoverPopup);
    }
    if (popup && popup.parentNode === svg) {
      svg.appendChild(popup);
    }
    if (viewportControls && viewportControls.parentNode === svg) {
      svg.appendChild(viewportControls);
    }
    if (searchControls && searchControls.parentNode === svg) {
      svg.appendChild(searchControls);
    }
  }

  function setupStickyLegend() {
    var legend = svg.querySelector('#legend');
    if (!legend || isElementHidden(legend)) return;
    var contentBBox = getElementBBox(legend);
    if (!contentBBox) return;
    ensureStickyLegendBackground(legend, contentBBox);
    var bbox = getElementBBox(legend) || contentBBox;
    var initialBounds = getElementBoundsInSvg(legend, bbox);
    if (!initialBounds) return;
    stickyLegend = legend;
    stickyLegendState = {
      localBBox: bbox,
      initialBounds: initialBounds,
      originalDisplay: legend.getAttribute('display')
    };
    setClassToken(stickyLegend, 'gbdraw-sticky-legend', true);
    stickyLegend.setAttribute('data-gbdraw-sticky-legend', 'true');
    updateStickyLegendHomeMetrics();
    syncStandaloneOverlayOrder();
    updateStickyLegendPosition();
  }

  function getSvgClientSize() {
    var rect = typeof svg.getBoundingClientRect === 'function' ? svg.getBoundingClientRect() : null;
    var width = rect && Number.isFinite(rect.width) && rect.width > 0 ? rect.width : window.innerWidth || originalViewRect.width;
    var height = rect && Number.isFinite(rect.height) && rect.height > 0 ? rect.height : window.innerHeight || originalViewRect.height;
    return {
      width: Math.max(1, width),
      height: Math.max(1, height)
    };
  }

  function normalizeViewRect(rect) {
    var minWidth = homeViewRect.width / maxZoom;
    var minHeight = homeViewRect.height / maxZoom;
    var width = clampValue(Number(rect && rect.width), minWidth, homeViewRect.width);
    var height = clampValue(Number(rect && rect.height), minHeight, homeViewRect.height);
    var minX = homeViewRect.x;
    var minY = homeViewRect.y;
    var maxX = homeViewRect.x + homeViewRect.width - width;
    var maxY = homeViewRect.y + homeViewRect.height - height;
    var x = maxX <= minX ? minX : clampValue(Number(rect && rect.x), minX, maxX);
    var y = maxY <= minY ? minY : clampValue(Number(rect && rect.y), minY, maxY);
    return { x: x, y: y, width: width, height: height };
  }

  function refreshActivePopupForViewport() {
    if (typeof updateActivePopupViewportMetrics === 'function') {
      updateActivePopupViewportMetrics();
      return;
    }
    keepPopupWithinVisibleView();
  }

  function refreshViewportOverlays() {
    updateViewportControlsPosition();
    refreshActivePopupForViewport();
    scheduleHoverPopupPosition(hoverPopupLastEvent);
  }

  function scheduleViewportOverlayRefresh() {
    var requestFrame = window.requestAnimationFrame || function (callback) {
      return window.setTimeout(callback, 16);
    };
    if (overlayRefreshFrame === null) {
      overlayRefreshFrame = requestFrame(function () {
        overlayRefreshFrame = null;
        refreshViewportOverlays();
        requestFrame(refreshViewportOverlays);
      });
    }
    if (overlayRefreshTimer === null) {
      overlayRefreshTimer = window.setTimeout(function () {
        overlayRefreshTimer = null;
        refreshViewportOverlays();
      }, 80);
    }
  }

  function setSvgViewRect(rect) {
    var next = normalizeViewRect(rect);
    svg.setAttribute(
      'viewBox',
      [next.x, next.y, next.width, next.height].map(formatSvgNumber).join(' ')
    );
    refreshViewportOverlays();
    scheduleViewportOverlayRefresh();
    return next;
  }

  function zoomViewBy(factor, anchorPoint) {
    var view = getViewRect();
    var safeFactor = Number(factor);
    if (!Number.isFinite(safeFactor) || safeFactor <= 0) return;
    var currentScale = view.width / homeViewRect.width;
    var nextScale = clampValue(currentScale * safeFactor, 1 / maxZoom, 1);
    var width = homeViewRect.width * nextScale;
    var height = homeViewRect.height * nextScale;
    var anchor = anchorPoint || { x: view.x + view.width / 2, y: view.y + view.height / 2 };
    var ratioX = width / view.width;
    var ratioY = height / view.height;
    closeHoverPopup();
    setSvgViewRect({
      x: anchor.x - (anchor.x - view.x) * ratioX,
      y: anchor.y - (anchor.y - view.y) * ratioY,
      width: width,
      height: height
    });
  }

  function resetViewport() {
    closeHoverPopup();
    closePopup();
    stopCanvasPan();
    setSvgViewRect(homeViewRect);
  }

  function refitViewportToWindow() {
    var previousHomeViewRect = copyViewRect(homeViewRect);
    var currentView = getViewRect();
    var wasAtHome = rectsNearlyEqual(currentView, previousHomeViewRect);
    var centerX = currentView.x + currentView.width / 2;
    var centerY = currentView.y + currentView.height / 2;
    var currentScale = clampValue(currentView.width / previousHomeViewRect.width, 1 / maxZoom, 1);
    homeViewRect = fitRectToAspect(originalViewRect, getViewportAspect());
    updateStickyLegendHomeMetrics();
    if (wasAtHome) {
      setSvgViewRect(homeViewRect);
      return;
    }
    var width = homeViewRect.width * currentScale;
    var height = homeViewRect.height * currentScale;
    setSvgViewRect({
      x: centerX - width / 2,
      y: centerY - height / 2,
      width: width,
      height: height
    });
  }

  function scheduleViewportRefit() {
    if (resizeFrame !== null) return;
    var requestFrame = window.requestAnimationFrame || function (callback) {
      return window.setTimeout(callback, 16);
    };
    resizeFrame = requestFrame(function () {
      resizeFrame = null;
      refitViewportToWindow();
    });
  }

  function scheduleInitialViewportRefresh() {
    var initialView = copyViewRect(getViewRect());
    var requestFrame = window.requestAnimationFrame || function (callback) {
      return window.setTimeout(callback, 16);
    };

    function refresh(refitIfStillInitial) {
      if (refitIfStillInitial && rectsNearlyEqual(getViewRect(), initialView)) {
        refitViewportToWindow();
        return;
      }
      updateViewportControlsPosition();
    }

    requestFrame(function () {
      refresh(false);
      requestFrame(function () {
        refresh(true);
      });
    });
    window.setTimeout(function () {
      refresh(true);
    }, 80);
  }

  function createSvgNode(tagName, attrs) {
    var node = document.createElementNS(SVG_NS, tagName);
    Object.keys(attrs || {}).forEach(function (key) {
      node.setAttribute(key, attrs[key]);
    });
    return node;
  }

  function createViewportButton(button, x) {
    var group = createSvgNode('g', {
      'class': 'gbdraw-viewport-button',
      'data-action': button.action,
      'role': 'button',
      'aria-label': button.title,
      'tabindex': '0',
      'transform': 'translate(' + x + ', 0)'
    });
    var title = createSvgNode('title', {});
    title.textContent = button.title;
    var rect = createSvgNode('rect', {
      x: 0,
      y: 0,
      width: button.width,
      height: 30
    });
    var text = createSvgNode('text', {
      x: button.width / 2,
      y: 15
    });
    text.textContent = button.label;
    group.appendChild(title);
    group.appendChild(rect);
    group.appendChild(text);
    return group;
  }

  function closestViewportButton(target) {
    var node = target;
    while (node && node !== svg) {
      if (node.getAttribute && node.getAttribute('data-action')) return node;
      node = node.parentNode;
    }
    return null;
  }

  function updateViewportControlsPosition() {
    updateStickyLegendPosition();
    if (!viewportControls) return;
    var visibleView = getVisibleViewRect();
    var view = getViewRect();
    var screenScale = getScreenScale();
    var fallbackSize = getSvgClientSize();
    var unit = Math.max(
      1 / Math.max(screenScale.x, 0.001),
      1 / Math.max(screenScale.y, 0.001),
      view.width / fallbackSize.width,
      view.height / fallbackSize.height
    );
    var margin = 12 * unit;
    viewportControls.setAttribute(
      'transform',
      'translate(' + formatSvgNumber(visibleView.x + margin) + ', ' + formatSvgNumber(visibleView.y + margin) + ') scale(' + formatSvgNumber(unit) + ')'
    );
    updateSearchControlsPosition();
  }

  function setupViewportControls() {
    var existing = svg.querySelector('#' + VIEWPORT_CONTROLS_ID);
    if (existing && existing.parentNode) {
      existing.parentNode.removeChild(existing);
    }
    viewportControls = createSvgNode('g', {
      id: VIEWPORT_CONTROLS_ID,
      'class': 'gbdraw-viewport-controls',
      'data-gbdraw-viewport-controls': 'true'
    });
    var buttons = [
      { action: 'zoom-in', label: '+', title: 'Zoom in', width: 32 },
      { action: 'zoom-out', label: '-', title: 'Zoom out', width: 32 },
      { action: 'reset', label: 'Original', title: 'Return to original view', width: 62 }
    ];
    var x = 0;
    buttons.forEach(function (button) {
      viewportControls.appendChild(createViewportButton(button, x));
      x += button.width + 6;
    });
    viewportControls.addEventListener('mousedown', function (event) {
      event.preventDefault();
      event.stopPropagation();
    });
    viewportControls.addEventListener('click', function (event) {
      var button = closestViewportButton(event.target);
      if (!button) return;
      var action = button.getAttribute('data-action');
      if (action === 'zoom-in') {
        zoomViewBy(0.8);
      } else if (action === 'zoom-out') {
        zoomViewBy(1.25);
      } else if (action === 'reset') {
        resetViewport();
      }
      event.preventDefault();
      event.stopPropagation();
    });
    viewportControls.addEventListener('keydown', function (event) {
      if (event.key !== 'Enter' && event.key !== ' ') return;
      var button = closestViewportButton(event.target);
      if (!button) return;
      button.dispatchEvent(new MouseEvent('click', { bubbles: true, cancelable: true }));
      event.preventDefault();
    });
    svg.appendChild(viewportControls);
    syncStandaloneOverlayOrder();
    updateViewportControlsPosition();
  }

  function stopCanvasPan() {
    if (!activeCanvasPan) return;
    document.removeEventListener('mousemove', activeCanvasPan.onMove);
    document.removeEventListener('mouseup', activeCanvasPan.onEnd);
    activeCanvasPan = null;
    setClassToken(svg, 'gbdraw-interactive-panning', false);
  }

  function startCanvasPan(event) {
    if (event.button !== 0) return;
    if (popup && popup.contains(event.target)) return;
    if (closestSearchControls(event.target)) return;
    if (closestViewportButton(event.target)) return;
    stopCanvasPan();
    var startView = getViewRect();
    var size = getSvgClientSize();
    var unitX = startView.width / size.width;
    var unitY = startView.height / size.height;
    var startClientX = event.clientX;
    var startClientY = event.clientY;
    var didPan = false;
    var panThresholdSq = 16;
    var onMove = function (moveEvent) {
      var dx = moveEvent.clientX - startClientX;
      var dy = moveEvent.clientY - startClientY;
      if (!didPan && ((dx * dx) + (dy * dy)) >= panThresholdSq) {
        didPan = true;
        if (activeCanvasPan) activeCanvasPan.didPan = true;
        closeHoverPopup();
        closePopup();
        setClassToken(svg, 'gbdraw-interactive-panning', true);
      }
      if (!didPan) return;
      setSvgViewRect({
        x: startView.x - dx * unitX,
        y: startView.y - dy * unitY,
        width: startView.width,
        height: startView.height
      });
      moveEvent.preventDefault();
    };
    var onEnd = function () {
      var wasPanned = Boolean(activeCanvasPan && activeCanvasPan.didPan);
      stopCanvasPan();
      if (wasPanned) {
        suppressNextCanvasClick = true;
        window.setTimeout(function () {
          suppressNextCanvasClick = false;
        }, 0);
      }
    };
    activeCanvasPan = {
      onMove: onMove,
      onEnd: onEnd,
      didPan: false
    };
    document.addEventListener('mousemove', onMove);
    document.addEventListener('mouseup', onEnd);
  }

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

  function sequenceKindLabel(sequenceKind) {
    return sequenceKind === 'aa' ? 'aa' : 'nt';
  }

  function sequenceExtension(sequenceKind) {
    return sequenceKindLabel(sequenceKind) === 'aa' ? 'faa' : 'fna';
  }

  function makeSafeFilename(value, fallback) {
    var cleaned = String(value || '').trim().replace(/[^\\w.-]+/g, '_').replace(/^_+|_+$/g, '');
    return cleaned || fallback || 'sequence';
  }

  function normalizeSequence(value) {
    return String(value || '').replace(/\\s+/g, '').trim();
  }

  function normalizeFastaHeaderText(value) {
    return String(value || '').replace(/\\s+/g, ' ').trim();
  }

  function normalizeFastaId(value) {
    return normalizeFastaHeaderText(value).replace(/^>+/, '').replace(/\\s+/g, '_') || 'sequence';
  }

  function wrapFastaSequence(sequence) {
    var text = normalizeSequence(sequence);
    var lines = [];
    var width = 60;
    for (var i = 0; i < text.length; i += width) {
      lines.push(text.slice(i, i + width));
    }
    return lines.join('\\n');
  }

  function formatFastaEntry(id, description, sequence) {
    var wrapped = wrapFastaSequence(sequence);
    if (!wrapped) return '';
    var headerDescription = normalizeFastaHeaderText(description);
    var header = '>' + normalizeFastaId(id);
    if (headerDescription) header += ' ' + headerDescription;
    return header + '\\n' + wrapped;
  }

  function featureDescription(feature) {
    return firstDisplayText(
      feature && feature.product,
      firstQualifierValue(feature, 'product'),
      feature && feature.gene,
      firstQualifierValue(feature, 'gene'),
      feature && (feature.locus_tag || feature.locusTag),
      firstQualifierValue(feature, 'locus_tag'),
      feature && (feature.display_label || feature.displayLabel || feature.label),
      feature && feature.type
    );
  }

  function memberFeatureSvgId(memberOrRow) {
    return String(memberOrRow && (
      memberOrRow.featureSvgId ||
      memberOrRow.feature_svg_id ||
      memberOrRow.svgId ||
      memberOrRow.svg_id
    ) || '').trim();
  }

  function featureForMember(memberOrRow) {
    var svgId = memberFeatureSvgId(memberOrRow);
    return svgId ? featuresById.get(svgId) || null : null;
  }

  function featureSequenceFilename(feature, sequenceKind) {
    var label = sequenceKindLabel(sequenceKind);
    var id = firstDisplayText(
      label === 'aa' ? displayProteinId(feature, null, '') : '',
      feature && (feature.display_label || feature.displayLabel || feature.label),
      feature && feature.svg_id,
      'feature'
    );
    return makeSafeFilename(id + '_' + label, 'feature_' + label) + '.' + sequenceExtension(label);
  }

  function featureFasta(feature, sequenceKind) {
    if (!feature) return '';
    var label = sequenceKindLabel(sequenceKind);
    var existing = label === 'aa'
      ? firstDisplayText(feature.amino_acid_fasta, feature.aminoAcidFasta)
      : firstDisplayText(feature.nucleotide_fasta, feature.nucleotideFasta);
    if (existing) return existing;
    var sequence = label === 'aa'
      ? firstDisplayText(feature.amino_acid_sequence, feature.aminoAcidSequence)
      : firstDisplayText(feature.nucleotide_sequence, feature.nucleotideSequence);
    if (!normalizeSequence(sequence)) return '';
    var id = label === 'aa'
      ? displayProteinId(feature, null, feature.svg_id || 'protein')
      : firstDisplayText(feature.record_id, feature.recordId, feature.svg_id || 'record') + ':' + locationText(feature);
    return formatFastaEntry(id, featureDescription(feature), sequence);
  }

  function memberFasta(memberOrRow, sequenceKind) {
    return featureFasta(featureForMember(memberOrRow), sequenceKind);
  }

  function memberSequenceFilename(memberOrRow, sequenceKind, orthogroupId) {
    var label = sequenceKindLabel(sequenceKind);
    var feature = featureForMember(memberOrRow);
    var memberId = firstDisplayText(
      memberOrRow && (memberOrRow.sourceProteinId || memberOrRow.source_protein_id),
      memberOrRow && (memberOrRow.proteinId || memberOrRow.protein_id),
      feature && displayProteinId(feature, null, ''),
      memberFeatureSvgId(memberOrRow),
      'member'
    );
    var groupId = firstDisplayText(orthogroupId, memberOrRow && (memberOrRow.orthogroupId || memberOrRow.orthogroup_id), 'orthogroup');
    return makeSafeFilename(groupId + '_' + memberId + '_' + label, 'orthogroup_member_' + label) + '.' + sequenceExtension(label);
  }

  function groupFasta(memberRows, sequenceKind) {
    return (Array.isArray(memberRows) ? memberRows : []).map(function (row) {
      return String(memberFasta(row, sequenceKind) || '').trim();
    }).filter(function (text) {
      return text;
    }).join('\\n');
  }

  function groupFastaCount(memberRows, sequenceKind) {
    return (Array.isArray(memberRows) ? memberRows : []).filter(function (row) {
      return String(memberFasta(row, sequenceKind) || '').trim();
    }).length;
  }

  function groupSequenceFilename(orthogroupId, displayName, sequenceKind) {
    var label = sequenceKindLabel(sequenceKind);
    var id = firstDisplayText(orthogroupId, 'orthogroup');
    var name = makeSafeFilename(firstDisplayText(displayName, id), id);
    return makeSafeFilename(id + '_' + name + '_' + label, 'orthogroup_' + label) + '.' + sequenceExtension(label);
  }

  function downloadText(filename, text, mimeType) {
    var safeFilename = filename || 'download.txt';
    var blob = new Blob([String(text == null ? '' : text)], { type: mimeType || 'text/plain;charset=utf-8' });
    var url = URL.createObjectURL(blob);
    var link = document.createElementNS(XHTML_NS, 'a');
    link.setAttribute('href', url);
    link.setAttribute('download', safeFilename);
    link.href = url;
    link.download = safeFilename;
    link.style.display = 'none';
    var parent = document.body || document.documentElement;
    parent.appendChild(link);
    link.click();
    if (link.parentNode) {
      link.parentNode.removeChild(link);
    }
    window.setTimeout(function () {
      URL.revokeObjectURL(url);
    }, 0);
  }

  function locationText(feature) {
    if (feature.location) return String(feature.location);
    var start = Number(feature.start);
    var end = Number(feature.end);
    var text = (Number.isFinite(start) ? start + 1 : '') + '..' + (Number.isFinite(end) ? end : '');
    if (feature.strand) text += ' (' + feature.strand + ')';
    return text;
  }

  function memberLocationText(member) {
    if (!member || typeof member !== 'object') return '';
    var start = Number(member.start);
    var end = Number(member.end);
    var text = (Number.isFinite(start) ? start + 1 : '') + '..' + (Number.isFinite(end) ? end : '');
    if (member.strand) text += ' (' + member.strand + ')';
    return text === '..' ? '' : text;
  }

  function getFeatureSearchDetail(feature) {
    var svgId = String(feature && feature.svg_id || '').trim();
    if (!svgId || !searchState.query) return '';
    var details = searchState.matchDetails[svgId] || [];
    return details.length ? formatSearchMatchDetail(details[0]) : '';
  }

  function orthogroupRows(feature) {
    var group = getFeatureOrthogroup(feature);
    var member = getFeatureOrthogroupMember(feature, group);
    var memberCount = Number(feature && feature.orthogroup_member_count);
    var recordCoverage = Number(feature && feature.orthogroup_record_coverage);
    var proteinId = displayProteinId(feature, member);
    var rows = [
      ['Orthogroup ID', feature && feature.orthogroup_id || ''],
      ['Orthogroup name', group && (group.display_name || group.name) || ''],
      ['Members', Number.isFinite(memberCount) && memberCount > 0 ? String(memberCount) : (group && group.member_count ? String(group.member_count) : '')],
      ['Record coverage', Number.isFinite(recordCoverage) && recordCoverage > 0 ? String(recordCoverage) : (group && group.record_coverage_count ? String(group.record_coverage_count) : '')],
      ['Protein ID', proteinId]
    ];
    return rows.filter(function (row) {
      return String(row[1] == null ? '' : row[1]) !== '';
    });
  }

  function detailRows(feature) {
    var searchDetail = getFeatureSearchDetail(feature);
    var rows = [
      ['Search match', searchDetail],
      ['Label', feature.display_label || feature.label || ''],
      ['Record ID', feature.record_id || ''],
      ['Type', feature.type || ''],
      ['Location', locationText(feature)]
    ];
    Array.prototype.push.apply(rows, orthogroupRows(feature));
    return rows.filter(function (row) {
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

  function copyButton(value, label) {
    var index = copyValues.push(String(value == null ? '' : value)) - 1;
    return '<button type="button" class="gfi-copy" data-copy-index="' + index + '">' + escapeHtml(label || 'Copy') + '</button>';
  }

  function downloadButton(filename, text, label) {
    var index = downloadValues.push({
      filename: filename,
      text: String(text == null ? '' : text),
      type: 'text/plain;charset=utf-8'
    }) - 1;
    return '<button type="button" class="gfi-copy" data-download-index="' + index + '">' + escapeHtml(label || 'DL') + '</button>';
  }

  function renderGroupSequenceActions(memberRows, options) {
    var rows = Array.isArray(memberRows) ? memberRows : [];
    var orthogroupId = firstDisplayText(options && options.orthogroupId, rows[0] && (rows[0].orthogroupId || rows[0].orthogroup_id));
    var displayName = firstDisplayText(options && options.displayName, rows[0] && (rows[0].displayName || rows[0].display_name), orthogroupId);
    var html = [];
    ['nt', 'aa'].forEach(function (sequenceKind) {
      var text = groupFasta(rows, sequenceKind);
      if (!text) return;
      var count = groupFastaCount(rows, sequenceKind);
      var suffix = count > 1 ? ' (' + count + ')' : '';
      html.push(copyButton(text, 'Copy ' + sequenceKind + suffix));
      html.push(downloadButton(groupSequenceFilename(orthogroupId, displayName, sequenceKind), text, 'DL ' + sequenceKind + suffix));
    });
    return html.join('');
  }

  function renderMemberSequenceActions(memberOrRow, orthogroupId) {
    var html = [];
    ['nt', 'aa'].forEach(function (sequenceKind) {
      var text = String(memberFasta(memberOrRow, sequenceKind) || '').trim();
      if (!text) return;
      html.push(copyButton(text, 'Copy ' + sequenceKind));
      html.push(downloadButton(memberSequenceFilename(memberOrRow, sequenceKind, orthogroupId), text, 'DL ' + sequenceKind));
    });
    return html.length ? '<div class="gfi-seq-actions">' + html.join('') + '</div>' : '';
  }

  var copyValues = [];
  var downloadValues = [];

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

  function memberProductOrNote(member) {
    return String(member && (member.product || member.note) || '').trim();
  }

  function orthogroupMemberTableRows(members, group) {
    var orthogroupId = String(group && (group.id || group.orthogroupId || group.orthogroup_id) || '').trim();
    var displayName = String(group && (group.display_name || group.displayName || group.name) || '').trim();
    return members.map(function (member) {
      var feature = featureForMember(member);
      return {
        featureSvgId: memberFeatureSvgId(member),
        orthogroupId: orthogroupId,
        displayName: displayName,
        record: firstDisplayText(member && (member.recordId || member.record_id), feature && (feature.record_id || feature.recordId)),
        coordinates: firstDisplayText(memberLocationText(member), feature ? locationText(feature) : ''),
        proteinId: displayProteinId(feature, member),
        productOrNote: firstDisplayText(memberProductOrNote(member), featureDescription(feature))
      };
    }).filter(function (row) {
      return row.record || row.coordinates || row.proteinId || row.productOrNote || row.featureSvgId;
    });
  }

  function renderOrthogroupMembers(feature) {
    var group = getFeatureOrthogroup(feature);
    if (!group) return '';
    var members = Array.isArray(group.members) ? group.members : [];
    if (!members.length) return '';
    var rows = orthogroupMemberTableRows(members, group);
    if (!rows.length) return '';
    var orthogroupId = String(group && (group.id || group.orthogroupId || group.orthogroup_id) || feature && (feature.orthogroup_id || feature.orthogroupId) || '').trim();
    var displayName = String(group && (group.display_name || group.displayName || group.name) || '').trim();
    var text = ['Record\\tCoordinates (+/-)\\tProtein ID\\tProduct / note'].concat(
      rows.map(function (row) {
        return [row.record, row.coordinates, row.proteinId, row.productOrNote].join('\\t');
      })
    ).join('\\n');
    var sequenceActions = renderGroupSequenceActions(rows, {
      orthogroupId: orthogroupId,
      displayName: displayName
    });
    var body = rows.map(function (row) {
      var sequenceActions = renderMemberSequenceActions(row, orthogroupId);
      return '<tr>' +
        '<td class="gfi-mono">' + escapeHtml(row.record) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(row.coordinates) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(row.proteinId) + '</td>' +
        '<td>' + escapeHtml(row.productOrNote) + '</td>' +
        '<td>' + sequenceActions + '</td>' +
        '</tr>';
    }).join('');
    return '<div class="gfi-block">' +
      '<div class="gfi-block-title"><span>Orthogroup members</span><span>' + members.length + '</span>' + copyButton(text) + '</div>' +
      '<div class="gfi-table-wrap"><table class="gfi-table gfi-og-members-table">' +
      '<thead><tr><th>Record</th><th>Coordinates (+/-)</th><th>Protein ID</th><th>Product / note</th><th>Seq</th></tr></thead>' +
      '<tbody>' + body + '</tbody></table></div>' +
      (sequenceActions ? '<div class="gfi-block-actions">' + sequenceActions + '</div>' : '') +
      '</div>';
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

  function renderSequenceBlock(title, sequence, fasta, sequenceKind, feature) {
    var text = String(sequenceKind ? featureFasta(feature, sequenceKind) : '');
    if (!text) text = String(fasta || sequence || '');
    if (!text) {
      return '<div class="gfi-block"><div class="gfi-block-title">' + escapeHtml(title) + '</div><div class="gfi-empty">No sequence available.</div></div>';
    }
    var label = sequenceKindLabel(sequenceKind);
    var download = sequenceKind
      ? downloadButton(featureSequenceFilename(feature, label), text, 'DL ' + label)
      : '';
    return '<div class="gfi-block">' +
      '<div class="gfi-block-title"><span>' + escapeHtml(title) + '</span>' + copyButton(text) + download + '</div>' +
      '<pre class="gfi-pre">' + escapeHtml(text) + '</pre>' +
      '</div>';
  }

  function renderSequences(feature) {
    var warnings = normalizeArray(feature.sequence_warnings);
    var warningHtml = warnings.map(function (warning) {
      return '<div class="gfi-warning">' + escapeHtml(warning) + '</div>';
    }).join('');
    return warningHtml +
      renderSequenceBlock('Nucleotide', feature.nucleotide_sequence, feature.nucleotide_fasta || feature.nucleotideFasta, 'nt', feature) +
      renderSequenceBlock('Amino acid', feature.amino_acid_sequence, feature.amino_acid_fasta || feature.aminoAcidFasta, 'aa', feature);
  }

  function renderSimplePopup(feature) {
    return '<div class="gfi gfi--simple">' +
      '<div class="gfi-header" data-drag-handle="true">' +
      '<div><div class="gfi-title">' + escapeHtml(feature.display_label || feature.label || feature.svg_id || 'Feature') + '</div>' +
      '<div class="gfi-subtitle">' + escapeHtml(locationText(feature)) + '</div></div>' +
      '<button type="button" class="gfi-close" data-close="true">x</button>' +
      '</div>' +
      '<div class="gfi-content">' + renderRows(detailRows(feature)) + '</div>' +
      '</div>';
  }

  function renderPopup(feature, activeTab) {
    copyValues = [];
    downloadValues = [];
    if (popupMode === 'simple') {
      return renderSimplePopup(feature);
    }
    var tab = activeTab || 'details';
    var panel = '';
    if (tab === 'qualifiers') {
      panel = renderQualifiers(feature);
    } else if (tab === 'sequence') {
      panel = renderSequences(feature);
    } else {
      panel = renderRows(detailRows(feature)) + renderOrthogroupMembers(feature);
    }
    function tabButton(id, label) {
      return '<button type="button" class="gfi-tab' + (tab === id ? ' is-active' : '') + '" data-tab="' + id + '">' + label + '</button>';
    }
    return '<div class="gfi">' +
      '<div class="gfi-header" data-drag-handle="true">' +
      '<div><div class="gfi-title">' + escapeHtml(feature.display_label || feature.label || feature.svg_id || 'Feature') + '</div>' +
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

  function normalizeMatchRows(rows) {
    return (Array.isArray(rows) ? rows : []).map(function (row) {
      if (Array.isArray(row)) return [row[0], row[1]];
      return [row && row.label, row && row.value];
    }).filter(function (row) {
      return String(row[0] == null ? '' : row[0]).trim() &&
        String(row[1] == null ? '' : row[1]).trim();
    });
  }

  function normalizeMatchMemberRows(rows) {
    return (Array.isArray(rows) ? rows : []).map(function (row) {
      return {
        featureSvgId: String(row && (row.featureSvgId || row.feature_svg_id || row.svgId || row.svg_id) || '').trim(),
        orthogroupId: String(row && (row.orthogroupId || row.orthogroup_id) || '').trim(),
        displayName: String(row && (row.displayName || row.display_name) || '').trim(),
        record: String(row && (row.record || row.record_id) || '').trim(),
        coordinates: String(row && row.coordinates || '').trim(),
        proteinId: String(row && (row.proteinId || row.protein_id) || '').trim(),
        productOrNote: String(row && (row.productOrNote || row.product_or_note) || '').trim()
      };
    }).filter(function (row) {
      return row.record || row.coordinates || row.proteinId || row.productOrNote || row.featureSvgId;
    });
  }

  function normalizeMatchBlockOrthogroups(groups) {
    return (Array.isArray(groups) ? groups : []).map(function (group) {
      var id = String(group && group.id || '').trim();
      var displayName = String(group && (group.displayName || group.display_name) || '').trim();
      var memberRows = normalizeMatchMemberRows(group && (group.memberRows || group.member_rows)).map(function (row) {
        if (!row.orthogroupId) row.orthogroupId = id;
        if (!row.displayName) row.displayName = displayName;
        return row;
      });
      return {
        id: id,
        displayName: displayName,
        memberCount: String(group && (group.memberCount || group.member_count) || '').trim(),
        recordCoverage: String(group && (group.recordCoverage || group.record_coverage) || '').trim(),
        queryMember: String(group && (group.queryMember || group.query_member) || '').trim(),
        subjectMember: String(group && (group.subjectMember || group.subject_member) || '').trim(),
        detailRows: normalizeMatchRows(group && (group.detailRows || group.detail_rows)),
        memberRows: memberRows,
        member_copy_text: String(group && (group.member_copy_text || group.memberCopyText) || '').trim()
      };
    }).filter(function (group) {
      return group.id;
    });
  }

  function normalizeMatchFeatureRows(rows) {
    return (Array.isArray(rows) ? rows : []).map(function (row, index) {
      var svgId = String(row && (row.svgId || row.svg_id) || '').trim();
      var canOpen = Boolean(row && (row.canOpen || row.can_open)) || Boolean(svgId && featuresById.has(svgId));
      return {
        key: String(row && row.key || svgId || index),
        svgId: svgId,
        canOpen: canOpen,
        label: String(row && row.label || '').trim(),
        record: String(row && row.record || '').trim(),
        location: String(row && row.location || '').trim(),
        proteinId: String(row && (row.proteinId || row.protein_id) || '').trim(),
        locusId: String(row && (row.locusId || row.locus_id) || '').trim(),
        displayName: String(row && (row.displayName || row.display_name) || '').trim(),
        product: String(row && row.product || '').trim(),
        copyText: String(row && (row.copyText || row.copy_text) || '').trim()
      };
    }).filter(function (row) {
      return row.svgId || row.label || row.record || row.location || row.product;
    });
  }

  function renderMatchBlockOrthogroupTable(groups, selectedId) {
    if (!groups.length) return '';
    var body = groups.map(function (group) {
      var isSelected = selectedId && group.id === selectedId;
      return '<tr class="gfi-block-og-row' + (isSelected ? ' is-active' : '') + '" data-block-og-id="' + escapeHtml(group.id) + '">' +
        '<td class="gfi-mono">' + escapeHtml(group.id) + '</td>' +
        '<td>' + escapeHtml(group.displayName) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(group.memberCount) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(group.recordCoverage) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(group.queryMember) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(group.subjectMember) + '</td>' +
        '</tr>';
    }).join('');
    return '<div class="gfi-table-wrap"><table class="gfi-table gfi-block-og-table">' +
      '<thead><tr><th>OG ID</th><th>Display name</th><th>Members</th><th>Record coverage</th><th>Query member</th><th>Subject member</th></tr></thead>' +
      '<tbody>' + body + '</tbody></table></div>';
  }

  function renderMatchFeatureTable(rows) {
    if (!rows.length) return '';
    var body = rows.map(function (row) {
      var sub = [row.locusId, row.displayName].filter(function (value) { return value; }).join(' / ');
      var copyText = row.copyText || [row.record, row.location, row.proteinId, row.locusId, row.displayName, row.product].join('\\t');
      var featureAttr = row.svgId ? ' data-match-feature-id="' + escapeHtml(row.svgId) + '"' : '';
      return '<tr class="gfi-match-feature-row' + (row.canOpen ? '' : ' is-disabled') + '"' + featureAttr + '>' +
        '<td><div class="gfi-match-feature-main">' + escapeHtml(row.label || row.proteinId || row.svgId || 'Feature') + '</div>' +
        (sub ? '<div class="gfi-match-feature-sub">' + escapeHtml(sub) + '</div>' : '') + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(row.record) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(row.location) + '</td>' +
        '<td>' + escapeHtml(row.product) + '</td>' +
        '<td>' + copyButton(copyText) + '</td>' +
        '</tr>';
    }).join('');
    return '<div class="gfi-table-wrap"><table class="gfi-match-feature-table">' +
      '<thead><tr><th>Feature</th><th>Record</th><th>Location</th><th>Product</th><th>Copy</th></tr></thead>' +
      '<tbody>' + body + '</tbody></table></div>';
  }

  function renderMatchMemberTable(section, rows) {
    if (!rows.length) return '';
    var orthogroupId = String(section && (section.id || section.orthogroupId || section.orthogroup_id) || rows[0] && (rows[0].orthogroupId || rows[0].orthogroup_id) || '').trim();
    var displayName = String(section && (section.displayName || section.display_name || section.name) || rows[0] && (rows[0].displayName || rows[0].display_name) || '').trim();
    var copyText = String(section && (section.member_copy_text || section.memberCopyText) || '').trim();
    if (!copyText) {
      copyText = ['Record\\tCoordinates (+/-)\\tProtein ID\\tProduct / note'].concat(
        rows.map(function (row) {
          return [row.record, row.coordinates, row.proteinId, row.productOrNote].join('\\t');
        })
      ).join('\\n');
    }
    var body = rows.map(function (row) {
      var sequenceActions = renderMemberSequenceActions(row, orthogroupId);
      return '<tr>' +
        '<td class="gfi-mono">' + escapeHtml(row.record) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(row.coordinates) + '</td>' +
        '<td class="gfi-mono">' + escapeHtml(row.proteinId) + '</td>' +
        '<td>' + escapeHtml(row.productOrNote) + '</td>' +
        '<td>' + sequenceActions + '</td>' +
        '</tr>';
    }).join('');
    var sequenceActions = renderGroupSequenceActions(rows, {
      orthogroupId: orthogroupId,
      displayName: displayName
    });
    return '<div class="gfi-table-wrap"><table class="gfi-table gfi-og-members-table">' +
      '<thead><tr><th>Record</th><th>Coordinates (+/-)</th><th>Protein ID</th><th>Product / note</th><th>Seq</th></tr></thead>' +
      '<tbody>' + body + '</tbody></table></div>' +
      '<div class="gfi-block-actions">' + copyButton(copyText) + sequenceActions + '</div>';
  }

  function renderMatchSections(match) {
    var sections = Array.isArray(match && match.sections) ? match.sections : [];
    if (!sections.length) {
      return '<div class="gfi-empty">No match details available.</div>';
    }
    var selectedBlockOrthogroupId = String(match && (match.selected_block_orthogroup_id || match.selectedBlockOrthogroupId) || '').trim();
    return sections.map(function (section) {
      var rows = normalizeMatchRows(section && section.rows);
      var memberRows = normalizeMatchMemberRows(section && (section.member_rows || section.memberRows));
      var blockOrthogroups = normalizeMatchBlockOrthogroups(section && (section.block_orthogroups || section.blockOrthogroups));
      var featureRows = normalizeMatchFeatureRows(section && (section.feature_rows || section.featureRows));
      if (!rows.length && !memberRows.length && !blockOrthogroups.length && !featureRows.length) return '';
      var selectedBlockOrthogroup = blockOrthogroups.find(function (group) {
        return group.id === selectedBlockOrthogroupId;
      });
      var selectedHtml = selectedBlockOrthogroup
        ? '<div class="gfi-block gfi-block--selected-og">' +
          '<div class="gfi-block-title">' + escapeHtml('Selected orthogroup') + '</div>' +
          renderRows(selectedBlockOrthogroup.detailRows) +
          renderMatchMemberTable(selectedBlockOrthogroup, selectedBlockOrthogroup.memberRows) +
          '</div>'
        : '';
      return '<div class="gfi-block">' +
        '<div class="gfi-block-title">' + escapeHtml(section && section.title || 'Details') + '</div>' +
        (rows.length && !featureRows.length ? renderRows(rows) : '') +
        renderMatchFeatureTable(featureRows) +
        renderMatchBlockOrthogroupTable(blockOrthogroups, selectedBlockOrthogroupId) +
        renderMatchMemberTable(section, memberRows) +
        selectedHtml +
        '</div>';
    }).join('') || '<div class="gfi-empty">No match details available.</div>';
  }

  function renderMatchPopup(match) {
    copyValues = [];
    downloadValues = [];
    return '<div class="gfi gfi--simple">' +
      '<div class="gfi-header" data-drag-handle="true">' +
      '<div><div class="gfi-title">' + escapeHtml(match && match.title || 'Pairwise match') + '</div>' +
      '<div class="gfi-subtitle">' + escapeHtml(match && (match.subtitle || match.id) || '') + '</div></div>' +
      '<button type="button" class="gfi-close" data-close="true">x</button>' +
      '</div>' +
      '<div class="gfi-content">' + renderMatchSections(match) + '</div>' +
      '<button type="button" class="gfi-resize-handle" data-resize="true" title="Drag to resize" aria-label="Resize popup"></button>' +
      '</div>';
  }

  function stopPopupResize() {
    if (!activePopupResize) return;
    document.removeEventListener('mousemove', activePopupResize.onMove, true);
    document.removeEventListener('mouseup', activePopupResize.onEnd, true);
    window.removeEventListener('mouseup', activePopupResize.onEnd, true);
    window.removeEventListener('blur', activePopupResize.onEnd);
    activePopupResize = null;
  }

  function stopPopupDrag() {
    if (!activePopupDrag) return;
    document.removeEventListener('mousemove', activePopupDrag.onMove, true);
    document.removeEventListener('mouseup', activePopupDrag.onEnd, true);
    window.removeEventListener('mouseup', activePopupDrag.onEnd, true);
    window.removeEventListener('blur', activePopupDrag.onEnd);
    activePopupDrag = null;
  }

  function closePopup() {
    stopPopupResize();
    stopPopupDrag();
    updateActivePopupViewportMetrics = null;
    if (popup && popup.parentNode) {
      popup.parentNode.removeChild(popup);
    }
    popup = null;
  }

  function closeHoverPopup() {
    if (hoverPopupTimer) {
      window.clearTimeout(hoverPopupTimer);
      hoverPopupTimer = null;
    }
    if (hoverPopupFrame) {
      var cancelFrame = window.cancelAnimationFrame || window.clearTimeout;
      cancelFrame(hoverPopupFrame);
      hoverPopupFrame = null;
    }
    if (hoverPopup && hoverPopup.parentNode) {
      hoverPopup.parentNode.removeChild(hoverPopup);
    }
    hoverPopup = null;
    hoverPopupFeatureId = '';
    hoverPopupLastEvent = null;
  }

  function firstValue(value) {
    var values = normalizeArray(value);
    for (var i = 0; i < values.length; i += 1) {
      var text = String(values[i] || '').trim();
      if (text) return text;
    }
    return '';
  }

  function firstQualifierValue(feature, key) {
    var qualifiers = getFeatureQualifiers(feature);
    var normalizedKey = String(key || '').trim().toLowerCase();
    if (!normalizedKey) return '';
    var direct = feature && feature[normalizedKey];
    if (direct !== null && direct !== undefined && direct !== '') return firstValue(direct);
    var exact = qualifiers[normalizedKey];
    if (exact !== null && exact !== undefined && exact !== '') return firstValue(exact);
    var keys = Object.keys(qualifiers);
    for (var i = 0; i < keys.length; i += 1) {
      if (String(keys[i]).toLowerCase() === normalizedKey) {
        return firstValue(qualifiers[keys[i]]);
      }
    }
    return '';
  }

  function featureLengthText(feature) {
    var start = Number(feature && feature.start);
    var end = Number(feature && feature.end);
    if (!Number.isFinite(start) || !Number.isFinite(end) || end < start) return '';
    return String(Math.round(end - start).toLocaleString()) + ' bp';
  }

  function getFeatureDisplayColor(feature, svgId) {
    var direct = String(feature && (feature.fill_color || feature.color) || '').trim();
    if (direct) return direct;
    var elements = featureElementsById.get(String(svgId || feature && feature.svg_id || '').trim()) || [];
    for (var i = 0; i < elements.length; i += 1) {
      var fill = String(elements[i].getAttribute('fill') || elements[i].style && elements[i].style.fill || '').trim();
      if (fill && fill.toLowerCase() !== 'none') return fill;
    }
    return '#94a3b8';
  }

  function hoverTitle(feature) {
    var primary = firstQualifierValue(feature, 'gene') ||
      firstQualifierValue(feature, 'locus_tag') ||
      firstQualifierValue(feature, 'product') ||
      String(feature && (feature.display_label || feature.label || feature.svg_id) || '').trim();
    var type = String(feature && feature.type || 'Feature').trim() || 'Feature';
    return primary && primary !== type ? type + ': ' + primary : type;
  }

  function hoverRows(feature) {
    var primary = String(feature && (feature.display_label || feature.label || '') || '').trim();
    var product = firstQualifierValue(feature, 'product') || String(feature && feature.product || '').trim();
    var gene = firstQualifierValue(feature, 'gene') || String(feature && feature.gene || '').trim();
    var locus = firstQualifierValue(feature, 'locus_tag') || String(feature && feature.locus_tag || '').trim();
    var note = firstQualifierValue(feature, 'note') || String(feature && feature.note || '').trim();
    var rows = [];
    if (gene && gene !== primary) rows.push(['Gene', gene]);
    if (locus && locus !== primary) rows.push(['Locus', locus]);
    if (product && product !== primary) rows.push(['Product', product, true]);
    if (note && note !== primary && note !== product) rows.push(['Note', note, true]);
    rows.push(['Length', featureLengthText(feature)]);
    rows.push(['Location', locationText(feature)]);
    rows.push(['Record', feature && feature.record_id || '']);
    rows.push(['Orthogroup', feature && feature.orthogroup_id || '']);
    return rows.filter(function (row) {
      return String(row[1] == null ? '' : row[1]).trim() !== '';
    });
  }

  function renderHoverPopupHtml(feature, svgId) {
    var rows = hoverRows(feature);
    var color = getFeatureDisplayColor(feature, svgId);
    var rowHtml = rows.slice(0, 7).map(function (row) {
      return '<div class="gfhs-row">' +
        '<div class="gfhs-key">' + escapeHtml(row[0]) + '</div>' +
        '<div class="gfhs-value' + (row[2] ? ' is-clamped' : '') + '">' + escapeHtml(row[1]) + '</div>' +
        '</div>';
    }).join('');
    return '<div class="gfhs">' +
      '<div class="gfhs-title">' +
      '<div class="gfhs-swatch" style="background:' + escapeHtml(color) + '"></div>' +
      '<div class="gfhs-text"><div class="gfhs-heading">' + escapeHtml(hoverTitle(feature)) + '</div>' +
      '<div class="gfhs-subtitle">' + escapeHtml(locationText(feature) || String(svgId || '')) + '</div></div>' +
      '</div>' +
      rowHtml +
      '</div>';
  }

  function renderMatchHoverPopupHtml(match) {
    var rows = normalizeMatchRows(match && match.hover_rows).slice(0, 6);
    var rowHtml = rows.map(function (row) {
      return '<div class="gfhs-row">' +
        '<div class="gfhs-key">' + escapeHtml(row[0]) + '</div>' +
        '<div class="gfhs-value">' + escapeHtml(row[1]) + '</div>' +
        '</div>';
    }).join('');
    var color = String(match && match.fill || '#94a3b8');
    return '<div class="gfhs">' +
      '<div class="gfhs-title">' +
      '<div class="gfhs-swatch" style="background:' + escapeHtml(color) + '"></div>' +
      '<div class="gfhs-text"><div class="gfhs-heading">' + escapeHtml(match && match.title || 'Pairwise match') + '</div>' +
      '<div class="gfhs-subtitle">' + escapeHtml(match && (match.subtitle || match.id) || '') + '</div></div>' +
      '</div>' +
      rowHtml +
      '</div>';
  }

  function getHoverPopupCssMetrics(viewport, rowCount) {
    var zoomScale = getBrowserZoomScale(viewport);
    var margin = 12;
    var width = Math.min(340, Math.max(1, viewport.width * zoomScale - margin * 2));
    var height = Math.min(250, Math.max(118, 68 + Math.min(Math.max(Number(rowCount) || 0, 1), 7) * 24));
    height = Math.min(height, Math.max(1, viewport.height * zoomScale - margin * 2));
    return {
      zoomScale: zoomScale,
      margin: margin,
      width: width,
      height: height
    };
  }

  function positionHoverPopup(event) {
    if (!hoverPopup) return;
    var viewport = getViewportClientRect();
    var rowCount = Number(hoverPopup.getAttribute('data-row-count')) || 1;
    var metrics = getHoverPopupCssMetrics(viewport, rowCount);
    var scale = getScreenScale();
    var effectiveScaleX = Math.max(scale.x, 0.001) * metrics.zoomScale;
    var effectiveScaleY = Math.max(scale.y, 0.001) * metrics.zoomScale;
    var width = metrics.width / effectiveScaleX;
    var height = metrics.height / effectiveScaleY;
    var offsetX = 14 / effectiveScaleX;
    var offsetY = 14 / effectiveScaleY;
    var marginX = metrics.margin / effectiveScaleX;
    var marginY = metrics.margin / effectiveScaleY;
    var point = eventPoint(event);
    var view = getVisibleViewRect();
    var x = point.x + offsetX;
    var y = point.y + offsetY;
    if (x + width + marginX > view.x + view.width) x = point.x - width - offsetX;
    if (y + height + marginY > view.y + view.height) y = point.y - height - offsetY;
    x = clampValue(x, view.x + marginX, view.x + view.width - width - marginX);
    y = clampValue(y, view.y + marginY, view.y + view.height - height - marginY);
    hoverPopup.setAttribute('x', x);
    hoverPopup.setAttribute('y', y);
    hoverPopup.setAttribute('width', width);
    hoverPopup.setAttribute('height', height);
    var root = hoverPopup.firstElementChild;
    if (root && root.style) {
      root.style.width = metrics.width + 'px';
      root.style.height = metrics.height + 'px';
      root.style.transformOrigin = '0 0';
      root.style.transform = 'scale(' + (1 / effectiveScaleX) + ', ' + (1 / effectiveScaleY) + ')';
    }
  }

  function scheduleHoverPopupPosition(event) {
    hoverPopupLastEvent = event || hoverPopupLastEvent;
    if (!hoverPopup || hoverPopupFrame) return;
    var requestFrame = window.requestAnimationFrame || function (callback) {
      return window.setTimeout(callback, 16);
    };
    hoverPopupFrame = requestFrame(function () {
      hoverPopupFrame = null;
      positionHoverPopup(hoverPopupLastEvent);
    });
  }

  function hoverPopupAllowed() {
    if (popup || activeCanvasPan) return false;
    if (window.matchMedia && !window.matchMedia('(hover: hover) and (pointer: fine)').matches) {
      return false;
    }
    return true;
  }

  function showHoverPopup(feature, svgId, event) {
    if (!feature || !hoverPopupAllowed()) {
      closeHoverPopup();
      return;
    }
    closeHoverPopup();
    var rows = hoverRows(feature);
    var foreignObject = document.createElementNS(SVG_NS, 'foreignObject');
    foreignObject.setAttribute('id', 'gbdraw-feature-hover-popup');
    foreignObject.setAttribute('class', 'gbdraw-feature-hover-popup');
    foreignObject.setAttribute('data-row-count', String(rows.length || 1));
    var root = document.createElementNS(XHTML_NS, 'div');
    root.setAttribute('xmlns', XHTML_NS);
    root.innerHTML = renderHoverPopupHtml(feature, svgId);
    foreignObject.appendChild(root);
    svg.appendChild(foreignObject);
    hoverPopup = foreignObject;
    hoverPopupFeatureId = String(svgId || '').trim();
    hoverPopupLastEvent = event;
    positionHoverPopup(event);
    syncStandaloneOverlayOrder();
  }

  function scheduleHoverPopup(feature, svgId, event) {
    if (!feature || !hoverPopupAllowed()) {
      closeHoverPopup();
      return;
    }
    hoverPopupLastEvent = event || hoverPopupLastEvent;
    if (hoverPopup && hoverPopupFeatureId === String(svgId || '').trim()) {
      scheduleHoverPopupPosition(event);
      return;
    }
    if (hoverPopupTimer) {
      window.clearTimeout(hoverPopupTimer);
      hoverPopupTimer = null;
    }
    hoverPopupFeatureId = String(svgId || '').trim();
    hoverPopupTimer = window.setTimeout(function () {
      hoverPopupTimer = null;
      showHoverPopup(feature, svgId, hoverPopupLastEvent);
    }, 180);
  }

  function showMatchHoverPopup(match, event) {
    if (!match || !hoverPopupAllowed()) {
      closeHoverPopup();
      return;
    }
    closeHoverPopup();
    var rows = normalizeMatchRows(match.hover_rows);
    var foreignObject = document.createElementNS(SVG_NS, 'foreignObject');
    foreignObject.setAttribute('id', 'gbdraw-feature-hover-popup');
    foreignObject.setAttribute('class', 'gbdraw-feature-hover-popup');
    foreignObject.setAttribute('data-row-count', String(rows.length || 1));
    var root = document.createElementNS(XHTML_NS, 'div');
    root.setAttribute('xmlns', XHTML_NS);
    root.innerHTML = renderMatchHoverPopupHtml(match);
    foreignObject.appendChild(root);
    svg.appendChild(foreignObject);
    hoverPopup = foreignObject;
    hoverPopupFeatureId = 'match:' + String(match.id || '').trim();
    hoverPopupLastEvent = event;
    positionHoverPopup(event);
    syncStandaloneOverlayOrder();
  }

  function scheduleMatchHoverPopup(match, event) {
    if (!match || !hoverPopupAllowed()) {
      closeHoverPopup();
      return;
    }
    hoverPopupLastEvent = event || hoverPopupLastEvent;
    var hoverId = 'match:' + String(match.id || '').trim();
    if (hoverPopup && hoverPopupFeatureId === hoverId) {
      scheduleHoverPopupPosition(event);
      return;
    }
    if (hoverPopupTimer) {
      window.clearTimeout(hoverPopupTimer);
      hoverPopupTimer = null;
    }
    hoverPopupFeatureId = hoverId;
    hoverPopupTimer = window.setTimeout(function () {
      hoverPopupTimer = null;
      showMatchHoverPopup(match, hoverPopupLastEvent);
    }, 180);
  }

  function keepPopupWithinVisibleView() {
    if (!popup) return;
    var viewport = getViewportClientRect();
    var view = getVisibleViewRect();
    var metrics = getPopupCssMetrics(viewport);
    var scale = getScreenScale();
    var effectiveScaleX = Math.max(scale.x, 0.001) * metrics.zoomScale;
    var effectiveScaleY = Math.max(scale.y, 0.001) * metrics.zoomScale;
    var marginX = metrics.margin / effectiveScaleX;
    var marginY = metrics.margin / effectiveScaleY;
    var width = parseFloat(popup.getAttribute('width')) || 1;
    var height = parseFloat(popup.getAttribute('height')) || 1;
    var minX = view.x + marginX;
    var minY = view.y + marginY;
    var maxX = view.x + view.width - width - marginX;
    var maxY = view.y + view.height - height - marginY;
    if (maxX < minX) maxX = minX;
    if (maxY < minY) maxY = minY;
    popup.setAttribute('x', clampValue(parseFloat(popup.getAttribute('x')), minX, maxX));
    popup.setAttribute('y', clampValue(parseFloat(popup.getAttribute('y')), minY, maxY));
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

  function getPopupCssMetrics(viewport) {
    var zoomScale = getBrowserZoomScale(viewport);
    var margin = 12;
    var maxWidth = popupMode === 'simple' ? 360 : 460;
    var maxHeight = popupMode === 'simple' ? 260 : 540;
    var minWidth = popupMode === 'simple' ? 240 : 300;
    var minHeight = popupMode === 'simple' ? 160 : 220;
    var availableVisualWidth = Math.max(1, viewport.width * zoomScale - margin * 2);
    var availableVisualHeight = Math.max(1, viewport.height * zoomScale - margin * 2);
    return {
      zoomScale: zoomScale,
      margin: margin,
      width: Math.max(1, Math.min(maxWidth, availableVisualWidth)),
      height: Math.max(1, Math.min(maxHeight, availableVisualHeight)),
      minWidth: Math.max(1, Math.min(minWidth, availableVisualWidth)),
      minHeight: Math.max(1, Math.min(minHeight, availableVisualHeight))
    };
  }

  function clampValue(value, min, max) {
    var safeMin = Number.isFinite(min) ? min : 0;
    var safeMax = Number.isFinite(max) ? Math.max(safeMin, max) : safeMin;
    var numeric = Number(value);
    if (!Number.isFinite(numeric)) return safeMin;
    return Math.min(Math.max(numeric, safeMin), safeMax);
  }

  function openPopup(feature, event, popupKind) {
    if (!supportsStandaloneControls()) return;
    closeHoverPopup();
    closePopup();
    var kind = popupKind === 'match' ? 'match' : 'feature';
    var viewport = getViewportClientRect();
    var view = getVisibleViewRect();
    var scale = getScreenScale();
    var safeScaleX = Math.max(scale.x, 0.001);
    var safeScaleY = Math.max(scale.y, 0.001);
    var metrics = getPopupCssMetrics(viewport);
    var marginCss = metrics.margin;
    var effectiveScaleX = safeScaleX * metrics.zoomScale;
    var effectiveScaleY = safeScaleY * metrics.zoomScale;
    var marginX = marginCss / effectiveScaleX;
    var marginY = marginCss / effectiveScaleY;
    var popupCssWidth = metrics.width;
    var popupCssHeight = metrics.height;
    var width = popupCssWidth / effectiveScaleX;
    var height = popupCssHeight / effectiveScaleY;
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
    foreignObject.setAttribute('data-gbdraw-popup-css-width', String(popupCssWidth));
    foreignObject.setAttribute('data-gbdraw-popup-css-height', String(popupCssHeight));

    var root = document.createElementNS(XHTML_NS, 'div');
    root.setAttribute('xmlns', XHTML_NS);
    root.style.width = popupCssWidth + 'px';
    root.style.height = popupCssHeight + 'px';
    root.style.transformOrigin = '0 0';
    root.style.transform = 'scale(' + (1 / effectiveScaleX) + ', ' + (1 / effectiveScaleY) + ')';

    function getResizeLimits() {
      viewport = getViewportClientRect();
      var resizeMetrics = getPopupCssMetrics(viewport);
      var rect = typeof foreignObject.getBoundingClientRect === 'function'
        ? foreignObject.getBoundingClientRect()
        : null;
      var layoutMargin = resizeMetrics.margin / resizeMetrics.zoomScale;
      var left = rect && Number.isFinite(rect.left) ? rect.left : layoutMargin;
      var top = rect && Number.isFinite(rect.top) ? rect.top : layoutMargin;
      var minWidth = resizeMetrics.minWidth;
      var minHeight = resizeMetrics.minHeight;
      return {
        minWidth: minWidth,
        minHeight: minHeight,
        maxWidth: Math.max(minWidth, (viewport.width - left) * resizeMetrics.zoomScale - resizeMetrics.margin),
        maxHeight: Math.max(minHeight, (viewport.height - top) * resizeMetrics.zoomScale - resizeMetrics.margin)
      };
    }

    function getPopupPositionLimits() {
      var currentViewport = getViewportClientRect();
      var currentView = getVisibleViewRect();
      var currentMetrics = getPopupCssMetrics(currentViewport);
      var currentScale = getScreenScale();
      var currentEffectiveScaleX = Math.max(currentScale.x, 0.001) * currentMetrics.zoomScale;
      var currentEffectiveScaleY = Math.max(currentScale.y, 0.001) * currentMetrics.zoomScale;
      var marginX = currentMetrics.margin / currentEffectiveScaleX;
      var marginY = currentMetrics.margin / currentEffectiveScaleY;
      var currentWidth = parseFloat(foreignObject.getAttribute('width')) || width;
      var currentHeight = parseFloat(foreignObject.getAttribute('height')) || height;
      var minX = currentView.x + marginX;
      var minY = currentView.y + marginY;
      var maxX = currentView.x + currentView.width - currentWidth - marginX;
      var maxY = currentView.y + currentView.height - currentHeight - marginY;
      if (maxX < minX) maxX = minX;
      if (maxY < minY) maxY = minY;
      return {
        minX: minX,
        minY: minY,
        maxX: maxX,
        maxY: maxY
      };
    }

    function setPopupPosition(nextX, nextY) {
      var limits = getPopupPositionLimits();
      foreignObject.setAttribute('x', clampValue(nextX, limits.minX, limits.maxX));
      foreignObject.setAttribute('y', clampValue(nextY, limits.minY, limits.maxY));
    }

    function applyPopupCssSize(nextCssWidth, nextCssHeight) {
      var limits = getResizeLimits();
      popupCssWidth = clampValue(nextCssWidth, limits.minWidth, limits.maxWidth);
      popupCssHeight = clampValue(nextCssHeight, limits.minHeight, limits.maxHeight);
      var currentScale = getScreenScale();
      safeScaleX = Math.max(currentScale.x, 0.001);
      safeScaleY = Math.max(currentScale.y, 0.001);
      metrics = getPopupCssMetrics(viewport);
      effectiveScaleX = safeScaleX * metrics.zoomScale;
      effectiveScaleY = safeScaleY * metrics.zoomScale;
      width = popupCssWidth / effectiveScaleX;
      height = popupCssHeight / effectiveScaleY;
      foreignObject.setAttribute('width', width);
      foreignObject.setAttribute('height', height);
      foreignObject.setAttribute('data-gbdraw-popup-css-width', String(popupCssWidth));
      foreignObject.setAttribute('data-gbdraw-popup-css-height', String(popupCssHeight));
      root.style.width = popupCssWidth + 'px';
      root.style.height = popupCssHeight + 'px';
      root.style.transform = 'scale(' + (1 / effectiveScaleX) + ', ' + (1 / effectiveScaleY) + ')';
    }

    updateActivePopupViewportMetrics = function () {
      if (!popup || popup !== foreignObject) return;
      applyPopupCssSize(popupCssWidth, popupCssHeight);
      keepPopupWithinVisibleView();
    };

    function startPopupResize(event) {
      if (event.button !== 0) return;
      stopPopupDrag();
      stopPopupResize();
      var startClientX = event.clientX;
      var startClientY = event.clientY;
      var startWidth = popupCssWidth;
      var startHeight = popupCssHeight;
      var onMove = function (moveEvent) {
        if (typeof moveEvent.buttons === 'number' && (moveEvent.buttons & 1) !== 1) {
          stopPopupResize();
          return;
        }
        var dragZoomScale = getBrowserZoomScale(getViewportClientRect());
        applyPopupCssSize(
          startWidth + (moveEvent.clientX - startClientX) * dragZoomScale,
          startHeight + (moveEvent.clientY - startClientY) * dragZoomScale
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
      document.addEventListener('mousemove', onMove, true);
      document.addEventListener('mouseup', onEnd, true);
      window.addEventListener('mouseup', onEnd, true);
      window.addEventListener('blur', onEnd);
      event.preventDefault();
      event.stopPropagation();
    }

    function startPopupDrag(event) {
      if (event.button !== 0) return;
      stopPopupDrag();
      stopPopupResize();
      var startPoint = eventPoint(event);
      var startX = parseFloat(foreignObject.getAttribute('x')) || 0;
      var startY = parseFloat(foreignObject.getAttribute('y')) || 0;
      var onMove = function (moveEvent) {
        if (typeof moveEvent.buttons === 'number' && (moveEvent.buttons & 1) !== 1) {
          stopPopupDrag();
          return;
        }
        var currentPoint = eventPoint(moveEvent);
        setPopupPosition(
          startX + (currentPoint.x - startPoint.x),
          startY + (currentPoint.y - startPoint.y)
        );
        moveEvent.preventDefault();
      };
      var onEnd = function () {
        stopPopupDrag();
      };
      activePopupDrag = {
        onMove: onMove,
        onEnd: onEnd
      };
      document.addEventListener('mousemove', onMove, true);
      document.addEventListener('mouseup', onEnd, true);
      window.addEventListener('mouseup', onEnd, true);
      window.addEventListener('blur', onEnd);
      event.preventDefault();
      event.stopPropagation();
    }

    function redraw() {
      root.innerHTML = kind === 'match' ? renderMatchPopup(feature) : renderPopup(feature, activeTab);
    }

    function closestFromTarget(target, selector) {
      var node = target && target.closest ? target : target && target.parentElement;
      return node && node.closest ? node.closest(selector) : null;
    }

    root.addEventListener('mousedown', function (rootEvent) {
      rootEvent.stopPropagation();
      var resizeHandle = closestFromTarget(rootEvent.target, '[data-resize]');
      if (resizeHandle) {
        startPopupResize(rootEvent);
        return;
      }
      var interactiveTarget = closestFromTarget(rootEvent.target, 'button, input, textarea, select, a, [data-close], [data-copy-index], [data-download-index], [data-tab], [data-block-og-id], [data-match-feature-id]');
      if (interactiveTarget) return;
      var dragHandle = closestFromTarget(rootEvent.target, '[data-drag-handle]');
      if (dragHandle) {
        startPopupDrag(rootEvent);
      }
    });

    root.addEventListener('click', function (rootEvent) {
      rootEvent.stopPropagation();
      var closeButton = closestFromTarget(rootEvent.target, '[data-close]');
      if (closeButton) {
        closePopup();
        return;
      }
      var tabButton = closestFromTarget(rootEvent.target, '[data-tab]');
      if (tabButton) {
        activeTab = tabButton.getAttribute('data-tab') || 'details';
        redraw();
        return;
      }
      var blockOrthogroupRow = closestFromTarget(rootEvent.target, '[data-block-og-id]');
      if (blockOrthogroupRow && kind === 'match') {
        feature.selected_block_orthogroup_id = blockOrthogroupRow.getAttribute('data-block-og-id') || '';
        redraw();
        return;
      }
      var copyTarget = closestFromTarget(rootEvent.target, '[data-copy-index]');
      if (copyTarget) {
        var index = Number(copyTarget.getAttribute('data-copy-index'));
        var value = Number.isFinite(index) ? copyValues[index] || '' : '';
        Promise.resolve(copyText(value, copyTarget)).catch(function () {});
        return;
      }
      var downloadTarget = closestFromTarget(rootEvent.target, '[data-download-index]');
      if (downloadTarget) {
        var downloadIndex = Number(downloadTarget.getAttribute('data-download-index'));
        var payload = Number.isFinite(downloadIndex) ? downloadValues[downloadIndex] : null;
        if (payload && payload.text) {
          downloadText(payload.filename, payload.text, payload.type);
        }
        return;
      }
      var matchFeatureRow = closestFromTarget(rootEvent.target, '[data-match-feature-id]');
      if (matchFeatureRow && kind === 'match') {
        var featureId = String(matchFeatureRow.getAttribute('data-match-feature-id') || '').trim();
        var matchedFeature = featureId ? featuresById.get(featureId) : null;
        if (matchedFeature) {
          openPopup(matchedFeature, rootEvent, 'feature');
        }
        return;
      }
    });

    root.addEventListener('wheel', function (rootEvent) {
      rootEvent.stopPropagation();
    }, { passive: true });

    ['mouseup', 'dblclick', 'touchstart', 'touchmove'].forEach(function (eventName) {
      root.addEventListener(eventName, function (rootEvent) {
        rootEvent.stopPropagation();
      }, { passive: eventName === 'touchstart' || eventName === 'touchmove' });
    });

    redraw();
    foreignObject.appendChild(root);
    svg.appendChild(foreignObject);
    popup = foreignObject;
    if (viewportControls) {
      svg.appendChild(viewportControls);
    }
    if (searchControls) {
      svg.appendChild(searchControls);
    }
    syncStandaloneOverlayOrder();
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

  function closestMatch(node) {
    while (node && node !== svg) {
      if (node.matches && node.matches(MATCH_SELECTOR)) return node;
      node = node.parentNode;
    }
    return null;
  }

  setSvgViewRect(homeViewRect);
  setClassToken(svg, 'gbdraw-interactive-pan-enabled', true);
  setupStickyLegend();
  setupViewportControls();
  setupSearchControls();
  window.addEventListener('scroll', updateViewportControlsPosition, { passive: true });
  window.addEventListener('resize', scheduleViewportRefit);
  if (window.visualViewport) {
    window.visualViewport.addEventListener('scroll', updateViewportControlsPosition, { passive: true });
    window.visualViewport.addEventListener('resize', updateViewportControlsPosition, { passive: true });
  }
  scheduleInitialViewportRefresh();

  svg.addEventListener('mousedown', startCanvasPan);

  svg.addEventListener('wheel', function (event) {
    var factor = event.deltaY > 0 ? 1.15 : 0.87;
    zoomViewBy(factor, eventPoint(event));
    event.preventDefault();
  }, { passive: false });

  svg.addEventListener('mouseover', function (event) {
    if (popup && popup.contains(event.target)) return;
    if (closestSearchControls(event.target)) return;
    var featureElement = closestFeature(event.target);
    if (featureElement) {
      clearActiveMatchHover();
      var svgId = getElementFeatureId(featureElement);
      if (!svgId || !featureElementsById.has(svgId)) return;
      var feature = featuresById.get(svgId);
      var hoverKey = getFeatureHoverKey(svgId);
      if (activeHoverKey === hoverKey) {
        scheduleHoverPopup(feature, svgId, event);
        return;
      }
      if (activeHoverSvgId) {
        setHoverHighlight(activeHoverSvgId, false);
      }
      activeHoverSvgId = svgId;
      activeHoverKey = hoverKey;
      setHoverHighlight(svgId, true);
      scheduleHoverPopup(feature, svgId, event);
      return;
    }
    var matchElement = closestMatch(event.target);
    if (!matchElement) return;
    clearActiveFeatureHover();
    var matchId = getElementMatchId(matchElement);
    if (!matchId || !matchesById.has(matchId)) return;
    var match = matchesById.get(matchId);
    var matchHoverKey = getMatchHoverKey(matchElement);
    if (activeHoverMatchKey === matchHoverKey) {
      scheduleMatchHoverPopup(match, event);
      return;
    }
    clearActiveMatchHover();
    activeHoverMatchId = matchId;
    activeHoverMatchKey = matchHoverKey;
    setMatchHover(matchElement, true);
    scheduleMatchHoverPopup(match, event);
  });

  svg.addEventListener('mousemove', function (event) {
    if (popup && popup.contains(event.target)) {
      closeHoverPopup();
      return;
    }
    if (closestSearchControls(event.target) || closestViewportButton(event.target)) {
      closeHoverPopup();
      return;
    }
    var featureElement = closestFeature(event.target);
    if (featureElement) {
      var svgId = getElementFeatureId(featureElement);
      if (!svgId || !featureElementsById.has(svgId)) {
        closeHoverPopup();
        return;
      }
      scheduleHoverPopup(featuresById.get(svgId), svgId, event);
      return;
    }
    var matchElement = closestMatch(event.target);
    if (matchElement) {
      var matchId = getElementMatchId(matchElement);
      var match = matchesById.get(matchId);
      if (match) {
        scheduleMatchHoverPopup(match, event);
        return;
      }
    }
    closeHoverPopup();
  });

  svg.addEventListener('mouseout', function (event) {
    var featureElement = closestFeature(event.target);
    if (featureElement) {
      var svgId = getElementFeatureId(featureElement);
      if (!svgId || activeHoverSvgId !== svgId) return;
      var relatedFeature = closestFeature(event.relatedTarget);
      if (relatedFeature && getFeatureHoverKey(getElementFeatureId(relatedFeature)) === activeHoverKey) return;
      clearActiveFeatureHover();
      closeHoverPopup();
      return;
    }
    var matchElement = closestMatch(event.target);
    if (!matchElement) return;
    var matchId = getElementMatchId(matchElement);
    if (!matchId || activeHoverMatchId !== matchId) return;
    var relatedMatch = closestMatch(event.relatedTarget);
    if (relatedMatch && getMatchHoverKey(relatedMatch) === activeHoverMatchKey) return;
    clearActiveMatchHover();
    closeHoverPopup();
  });

  svg.addEventListener('click', function (event) {
    if (popup && popup.contains(event.target)) return;
    if (closestSearchControls(event.target)) return;
    if (suppressNextCanvasClick) {
      suppressNextCanvasClick = false;
      event.preventDefault();
      event.stopPropagation();
      return;
    }
    var featureElement = closestFeature(event.target);
    if (!featureElement) {
      var matchElement = closestMatch(event.target);
      if (!matchElement) {
        closeHoverPopup();
        closePopup();
        return;
      }
      var matchId = getElementMatchId(matchElement);
      var match = matchesById.get(matchId);
      if (!match) return;
      event.preventDefault();
      event.stopPropagation();
      closeHoverPopup();
      openPopup(match, event, 'match');
      return;
    }
    var svgId = getElementFeatureId(featureElement);
    var feature = featuresById.get(svgId);
    if (!feature) return;
    var matchIndex = searchState.matches.indexOf(svgId);
    if (matchIndex !== -1) {
      setActiveMatch(matchIndex, { center: false });
    }
    event.preventDefault();
    event.stopPropagation();
    closeHoverPopup();
    openPopup(feature, event, 'feature');
  });

  document.addEventListener('keydown', function (event) {
    if (event.key === 'Escape') {
      closeHoverPopup();
      closePopup();
      stopCanvasPan();
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

const normalizeFeatureIdKey = (value) => String(value || '').trim().toLowerCase();

const normalizeStandaloneContext = (options = {}) => ({
  popupMode: normalizeStandalonePopupMode(options.popupMode),
  features: Array.isArray(options.features) ? options.features : [],
  editableLabels: Array.isArray(options.editableLabels) ? options.editableLabels : [],
  labelTextFeatureOverrides:
    options.labelTextFeatureOverrides && typeof options.labelTextFeatureOverrides === 'object'
      ? options.labelTextFeatureOverrides
      : {},
  labelTextBulkOverrides:
    options.labelTextBulkOverrides && typeof options.labelTextBulkOverrides === 'object'
      ? options.labelTextBulkOverrides
      : {},
  featureOrthogroupIndex: options.featureOrthogroupIndex instanceof Map
    ? options.featureOrthogroupIndex
    : new Map(),
  orthogroups: Array.isArray(options.orthogroups) ? options.orthogroups : [],
  orthogroupNameOverrides:
    options.orthogroupNameOverrides && typeof options.orthogroupNameOverrides === 'object'
      ? options.orthogroupNameOverrides
      : {},
  orthogroupDescriptionOverrides:
    options.orthogroupDescriptionOverrides && typeof options.orthogroupDescriptionOverrides === 'object'
      ? options.orthogroupDescriptionOverrides
      : {},
  legendEntries: Array.isArray(options.legendEntries) ? options.legendEntries : [],
  currentColors:
    options.currentColors && typeof options.currentColors === 'object'
      ? options.currentColors
      : {}
});

const getEditableLabelEntryForStandaloneFeature = (feature, context) => {
  const entries = Array.isArray(context?.editableLabels) ? context.editableLabels : [];
  const candidates = [
    normalizeFeatureIdKey(feature?.svg_id),
    normalizeFeatureIdKey(feature?.id)
  ].filter(Boolean);
  if (candidates.length === 0) return null;
  return entries.find((entry) => candidates.includes(normalizeFeatureIdKey(entry?.featureId))) || null;
};

const getLabelTextFeatureOverride = (feature, context) => {
  const overrides = context?.labelTextFeatureOverrides || {};
  const candidates = [
    String(feature?.svg_id || '').trim(),
    String(feature?.id || '').trim()
  ].filter(Boolean);
  for (const candidate of candidates) {
    if (Object.prototype.hasOwnProperty.call(overrides, candidate)) {
      const text = String(overrides[candidate] ?? '').trim();
      if (text) return text;
    }
  }
  return '';
};

const getStandaloneDisplayLabel = (feature, fallbackLabel, context) => {
  const editableEntry = getEditableLabelEntryForStandaloneFeature(feature, context);
  const editableText = String(editableEntry?.text || '').trim();
  if (editableText) return editableText;

  const directOverride = getLabelTextFeatureOverride(feature, context);
  if (directOverride) return directOverride;

  const sourceCandidates = [
    editableEntry?.sourceText,
    fallbackLabel,
    feature?.product,
    feature?.gene,
    feature?.locus_tag
  ].map((value) => String(value || '').trim()).filter(Boolean);
  for (const sourceText of sourceCandidates) {
    const bulkText = String(context?.labelTextBulkOverrides?.[sourceText] ?? '').trim();
    if (bulkText) return bulkText;
  }

  return fallbackLabel;
};

const getStandaloneSearchLabels = (feature, fallbackLabel, displayLabel) => {
  const labels = [
    displayLabel,
    fallbackLabel,
    feature?.label,
    feature?.gene,
    feature?.locus_tag,
    feature?.product,
    firstQualifierValue(feature, 'gene'),
    firstQualifierValue(feature, 'locus_tag'),
    firstQualifierValue(feature, 'product'),
    feature?.svg_id
  ];
  const seen = new Set();
  return labels
    .map((value) => String(value || '').trim())
    .filter((value) => {
      const key = value.toLowerCase();
      if (!value || seen.has(key)) return false;
      seen.add(key);
      return true;
    });
};

const buildStandaloneOrthogroupIndexKey = (recordIndex, svgId) => `${Number(recordIndex)}:${String(svgId || '').trim()}`;

const getStandaloneFeatureOrthogroupEntry = (feature, context) => {
  const directOrthogroupId = String(feature?.orthogroupId || feature?.orthogroup_id || '').trim();
  if (directOrthogroupId) {
    return {
      orthogroupId: directOrthogroupId,
      orthogroupMemberCount: Number(feature?.orthogroupMemberCount || feature?.orthogroup_member_count || 0),
      orthogroupRecordCoverage: Number(feature?.orthogroupRecordCoverage || feature?.orthogroup_record_coverage || 0),
      proteinId: String(feature?.proteinId || feature?.protein_id || '').trim(),
      sourceProteinId: String(feature?.sourceProteinId || feature?.source_protein_id || '').trim(),
      orthogroupRepresentative: Boolean(feature?.orthogroupRepresentative || feature?.orthogroup_representative),
      orthogroupMember: feature?.orthogroupMember || feature?.orthogroup_member || null
    };
  }

  const index = context?.featureOrthogroupIndex instanceof Map ? context.featureOrthogroupIndex : new Map();
  const svgId = String(feature?.svg_id || '').trim();
  if (!svgId || index.size === 0) return null;

  const candidateRecordIndexes = [
    feature?.fileIdx,
    feature?.file_idx,
    feature?.record_idx,
    feature?.recordIndex,
    feature?.record_index
  ];
  for (const recordIndexRaw of candidateRecordIndexes) {
    const recordIndex = Number(recordIndexRaw);
    if (!Number.isInteger(recordIndex)) continue;
    const entry = index.get(buildStandaloneOrthogroupIndexKey(recordIndex, svgId));
    if (entry) return entry;
  }
  return index.get(svgId) || null;
};

const normalizeStandaloneOrthogroupMember = (member) => {
  if (!member || typeof member !== 'object' || Array.isArray(member)) return null;
  const start = Number(member.start);
  const end = Number(member.end);
  const recordIndex = Number(member.recordIndex ?? member.record_index);
  return {
    orthogroupId: String(member.orthogroupId || member.orthogroup_id || ''),
    proteinId: String(member.proteinId || member.protein_id || ''),
    sourceProteinId: String(member.sourceProteinId || member.source_protein_id || ''),
    recordIndex: Number.isInteger(recordIndex) ? recordIndex : null,
    recordId: String(member.recordId || member.record_id || ''),
    featureIndex: Number.isFinite(Number(member.featureIndex ?? member.feature_index))
      ? Number(member.featureIndex ?? member.feature_index)
      : null,
    label: String(member.label || ''),
    featureSvgId: String(member.featureSvgId || member.feature_svg_id || ''),
    start: Number.isFinite(start) ? start : null,
    end: Number.isFinite(end) ? end : null,
    strand: String(member.strand || ''),
    representative: Boolean(member.representative),
    gene: String(member.gene || ''),
    locusTag: String(member.locusTag || member.locus_tag || ''),
    geneId: String(member.geneId || member.gene_id || ''),
    oldLocusTag: String(member.oldLocusTag || member.old_locus_tag || ''),
    product: String(member.product || ''),
    note: String(member.note || '')
  };
};

const normalizeStandaloneOrthogroupCandidate = (candidate) => {
  if (!candidate || typeof candidate !== 'object' || Array.isArray(candidate)) return null;
  return {
    text: String(candidate.text || ''),
    source: String(candidate.source || ''),
    memberCount: Number.isFinite(Number(candidate.memberCount)) ? Number(candidate.memberCount) : 0,
    recordCoverageCount: Number.isFinite(Number(candidate.recordCoverageCount)) ? Number(candidate.recordCoverageCount) : 0,
    representativeCount: Number.isFinite(Number(candidate.representativeCount)) ? Number(candidate.representativeCount) : 0,
    score: Number.isFinite(Number(candidate.score)) ? Number(candidate.score) : 0
  };
};

const getStandaloneOrthogroupDisplayName = (group, context) => {
  const id = String(group?.id || '').trim();
  const override = id ? String(context?.orthogroupNameOverrides?.[id] || '').trim() : '';
  return override || String(group?.display_name || group?.displayName || group?.name || id).trim();
};

const getStandaloneOrthogroupDescription = (group, context) => {
  const id = String(group?.id || '').trim();
  const override = id ? String(context?.orthogroupDescriptionOverrides?.[id] || '').trim() : '';
  return override || String(group?.description || '').trim();
};

const buildStandaloneOrthogroupPayloads = (features, context) => {
  const neededIds = new Set(
    (Array.isArray(features) ? features : [])
      .map((feature) => String(feature?.orthogroup_id || feature?.orthogroupId || '').trim())
      .filter(Boolean)
  );
  if (neededIds.size === 0) return [];

  const groups = Array.isArray(context?.orthogroups) ? context.orthogroups : [];
  return groups
    .filter((group) => neededIds.has(String(group?.id || '').trim()))
    .map((group) => {
      const members = Array.isArray(group?.members)
        ? group.members.map(normalizeStandaloneOrthogroupMember).filter(Boolean)
        : [];
      const memberCount = Number(group?.member_count || group?.memberCount || members.length || 0);
      const recordCoverageFallback = new Set(
        members
          .map((member) => Number(member.recordIndex))
          .filter((recordIndex) => Number.isInteger(recordIndex))
      ).size;
      const recordCoverageCount = Number(group?.record_coverage_count || group?.recordCoverageCount || recordCoverageFallback || 0);
      return {
        id: String(group?.id || ''),
        name: String(group?.name || ''),
        display_name: getStandaloneOrthogroupDisplayName(group, context),
        description: getStandaloneOrthogroupDescription(group, context),
        nameConfidence: String(group?.nameConfidence || ''),
        nameCandidates: Array.isArray(group?.nameCandidates)
          ? group.nameCandidates.map(normalizeStandaloneOrthogroupCandidate).filter(Boolean)
          : [],
        member_count: Number.isFinite(memberCount) ? memberCount : members.length,
        record_coverage_count: Number.isFinite(recordCoverageCount) ? recordCoverageCount : 0,
        members
      };
    });
};

const collectRenderedFeatureIds = (svg) => {
  const ids = new Set();
  if (!svg) return ids;
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = getElementFeatureId(element);
    if (id) ids.add(id);
  });
  return ids;
};

const collectRenderedFeatureEntries = (svg) => {
  const entries = new Map();
  if (!svg) return entries;
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = getElementFeatureId(element);
    if (!id) return;
    if (!entries.has(id)) {
      entries.set(id, { id, elements: [] });
    }
    entries.get(id).elements.push(element);
  });
  return entries;
};

const normalizeColorKey = (value) => String(value || '').trim().toLowerCase();

const getRenderedElementFill = (element) => {
  const fill = String(element?.getAttribute?.('fill') || element?.style?.fill || '').trim();
  if (fill && fill.toLowerCase() !== 'none') return fill;
  return '';
};

const getRenderedFeatureFill = (entry) => {
  const elements = Array.isArray(entry?.elements) ? entry.elements : [];
  for (const element of elements) {
    const fill = getRenderedElementFill(element);
    if (fill) return fill;
  }
  return '';
};

const collectLegendCaptionsByColor = (svg, context) => {
  const captionsByColor = new Map();
  const addCaption = (color, caption) => {
    const key = normalizeColorKey(color);
    const text = String(caption || '').trim();
    if (key && text && !captionsByColor.has(key)) {
      captionsByColor.set(key, text);
    }
  };

  const legendEntries = Array.isArray(context?.legendEntries) ? context.legendEntries : [];
  legendEntries.forEach((entry) => addCaption(entry?.color, entry?.caption || entry?.originalCaption));

  const paletteColors = context?.currentColors && typeof context.currentColors === 'object'
    ? context.currentColors
    : {};
  Object.entries(paletteColors).forEach(([caption, color]) => addCaption(color, caption));

  if (!svg) return captionsByColor;
  const legendRoots = Array.from(svg.querySelectorAll('#legend, g[id*="legend"], [data-gbdraw-sticky-legend]'));
  legendRoots.forEach((root) => {
    root.querySelectorAll('rect, path, polygon').forEach((shape) => {
      const color = getRenderedElementFill(shape);
      if (!color) return;
      const group = shape.closest?.('g') || root;
      const groupText = group?.querySelector?.('text') || null;
      const siblingText = shape.nextElementSibling?.matches?.('text') ? shape.nextElementSibling : null;
      const textElement = groupText || siblingText;
      const caption = String(textElement?.textContent || '').trim();
      addCaption(color, caption);
    });
  });
  return captionsByColor;
};

const buildFallbackStandaloneFeaturePayload = (svgId, entry, captionsByColor) => {
  const fillColor = getRenderedFeatureFill(entry);
  const caption = captionsByColor.get(normalizeColorKey(fillColor)) || 'Feature';
  const label = caption === 'Feature' ? String(svgId) : caption;
  const searchLabels = Array.from(new Set([label, caption, svgId].map((value) => String(value || '').trim()).filter(Boolean)));
  return {
    svg_id: String(svgId || ''),
    label,
    display_label: label,
    search_labels: searchLabels,
    record_id: '',
    record_idx: null,
    type: caption,
    start: null,
    end: null,
    strand: '',
    location: '',
    locus_tag: '',
    gene_id: '',
    old_locus_tag: '',
    fill_color: fillColor,
    orthogroup_id: '',
    orthogroup_member_count: 0,
    orthogroup_record_coverage: 0,
    protein_id: '',
    source_protein_id: '',
    orthogroup_representative: false,
    qualifiers: {},
    location_parts: [],
    nucleotide_sequence: '',
    amino_acid_sequence: '',
    sequence_warnings: [],
    nucleotide_fasta: '',
    amino_acid_fasta: '',
    orthogroup_member: null
  };
};

const normalizeStandalonePopupMode = (popupMode) => (
  popupMode === 'simple' ? 'simple' : 'rich'
);

const buildStandaloneFeaturePayloads = (svg, options = {}) => {
  const renderedEntries = collectRenderedFeatureEntries(svg);
  const renderedIds = new Set(renderedEntries.keys());
  if (renderedIds.size === 0) return [];

  const context = normalizeStandaloneContext(options);
  const normalizedPopupMode = normalizeStandalonePopupMode(context.popupMode);
  const features = context.features;
  const payloads = [];
  const seenIds = new Set();
  features.forEach((feature) => {
    const svgId = String(feature?.svg_id || '').trim();
    if (!svgId || !renderedIds.has(svgId) || seenIds.has(svgId)) return;
    seenIds.add(svgId);
    const fallbackLabel = getStandaloneFeatureLabel(feature);
    const displayLabel = getStandaloneDisplayLabel(feature, fallbackLabel, context);
    const orthogroupEntry = getStandaloneFeatureOrthogroupEntry(feature, context);
    const orthogroupMember = normalizeStandaloneOrthogroupMember(orthogroupEntry?.orthogroupMember);
    const sequenceFastas = buildFeatureSequenceFastas(feature);
    const payload = {
      svg_id: svgId,
      label: fallbackLabel,
      display_label: displayLabel,
      search_labels: getStandaloneSearchLabels(feature, fallbackLabel, displayLabel),
      record_id: String(feature?.record_id || ''),
      record_idx: Number.isFinite(Number(feature?.record_idx)) ? Number(feature.record_idx) : null,
      type: String(feature?.type || ''),
      start: Number.isFinite(Number(feature?.start)) ? Number(feature.start) : null,
      end: Number.isFinite(Number(feature?.end)) ? Number(feature.end) : null,
      strand: String(feature?.strand || ''),
      location: buildStandaloneFeatureLocation(feature),
      locus_tag: String(feature?.locus_tag || feature?.locusTag || ''),
      gene_id: String(feature?.gene_id || feature?.geneId || ''),
      old_locus_tag: String(feature?.old_locus_tag || feature?.oldLocusTag || ''),
      orthogroup_id: String(orthogroupEntry?.orthogroupId || ''),
      orthogroup_member_count: Number.isFinite(Number(orthogroupEntry?.orthogroupMemberCount))
        ? Number(orthogroupEntry.orthogroupMemberCount)
        : 0,
      orthogroup_record_coverage: Number.isFinite(Number(orthogroupEntry?.orthogroupRecordCoverage))
        ? Number(orthogroupEntry.orthogroupRecordCoverage)
        : 0,
      protein_id: String(orthogroupEntry?.proteinId || feature?.proteinId || feature?.protein_id || ''),
      source_protein_id: String(
        orthogroupEntry?.sourceProteinId ||
        feature?.sourceProteinId ||
        feature?.source_protein_id ||
        firstQualifierValue(feature, 'protein_id') ||
        ''
      ),
      orthogroup_representative: Boolean(orthogroupEntry?.orthogroupRepresentative)
    };
    if (normalizedPopupMode === 'rich') {
      Object.assign(payload, {
        qualifiers: normalizeQualifierMap(feature?.qualifiers),
        location_parts: normalizeLocationParts(feature?.location_parts),
        nucleotide_sequence: String(feature?.nucleotide_sequence || ''),
        amino_acid_sequence: String(feature?.amino_acid_sequence || ''),
        sequence_warnings: normalizeStringArray(feature?.sequence_warnings),
        nucleotide_fasta: sequenceFastas.nucleotideFasta,
        amino_acid_fasta: sequenceFastas.aminoAcidFasta,
        orthogroup_member: orthogroupMember
      });
    }
    payloads.push(payload);
  });
  if (payloads.length < renderedIds.size) {
    const captionsByColor = collectLegendCaptionsByColor(svg, context);
    renderedEntries.forEach((entry, svgId) => {
      if (seenIds.has(svgId)) return;
      seenIds.add(svgId);
      payloads.push(buildFallbackStandaloneFeaturePayload(svgId, entry, captionsByColor));
    });
  }
  return payloads;
};

const standaloneAttr = (element, name) =>
  String(element?.getAttribute?.(name) || '').trim();

const firstStandaloneText = (...values) => {
  for (const value of values) {
    const text = String(value === null || value === undefined ? '' : value).trim();
    if (text) return text;
  }
  return '';
};

const splitStandaloneMetadataValues = (value) => firstStandaloneText(value)
  .split(';')
  .map((entry) => firstStandaloneText(entry))
  .filter(Boolean);

const uniqueStandaloneMetadataValues = (value) => {
  const seen = new Set();
  const values = [];
  splitStandaloneMetadataValues(value).forEach((entry) => {
    if (seen.has(entry)) return;
    seen.add(entry);
    values.push(entry);
  });
  return values;
};

const standaloneGeneratedProteinIdPattern = /^(?:gbd_r\d+_cds\d+|p_.+_\d+_\d+_-?\d+_[0-9a-f]{12}(?:_\d+)?)$/i;
const standaloneGeneratedUnitIdPattern = /^gbd_r\d+_unit\d+$/i;

const isStandaloneInternalDisplayId = (value) => {
  const text = firstStandaloneText(value);
  return Boolean(text && (
    standaloneGeneratedProteinIdPattern.test(text) ||
    standaloneGeneratedUnitIdPattern.test(text)
  ));
};

const addUniqueStandaloneDisplayText = (values, value) => {
  const text = firstStandaloneText(value);
  if (!text || isStandaloneInternalDisplayId(text) || values.includes(text)) return;
  values.push(text);
};

const resolveStandaloneDisplayProteinId = (feature, member = null, fallback = '') => firstStandaloneText(
  feature?.sourceProteinId,
  feature?.source_protein_id,
  member?.sourceProteinId,
  member?.source_protein_id,
  firstQualifierValue(feature, 'protein_id'),
  feature?.locusTag,
  feature?.locus_tag,
  firstQualifierValue(feature, 'locus_tag'),
  member?.locusTag,
  member?.locus_tag,
  feature?.geneId,
  feature?.gene_id,
  firstQualifierValue(feature, 'gene_id'),
  member?.geneId,
  member?.gene_id,
  feature?.oldLocusTag,
  feature?.old_locus_tag,
  firstQualifierValue(feature, 'old_locus_tag'),
  member?.oldLocusTag,
  member?.old_locus_tag,
  feature?.ID,
  firstQualifierValue(feature, 'ID'),
  feature?.Name,
  firstQualifierValue(feature, 'Name'),
  feature?.Parent,
  firstQualifierValue(feature, 'Parent'),
  feature?.gene,
  firstQualifierValue(feature, 'gene'),
  member?.gene,
  feature?.proteinId,
  feature?.protein_id,
  member?.proteinId,
  member?.protein_id,
  fallback
);

const addStandaloneMatchRow = (rows, label, value) => {
  const text = String(value === null || value === undefined ? '' : value).trim();
  if (!text) return;
  rows.push([label, text]);
};

const standaloneIntervalText = (start, end) => {
  const startText = String(start || '').trim();
  const endText = String(end || '').trim();
  if (startText && endText) return `${startText}..${endText}`;
  return startText || endText;
};

const standaloneMemberLocationText = (member) => {
  if (!member || typeof member !== 'object') return '';
  const start = Number(member.start);
  const end = Number(member.end);
  const startText = Number.isFinite(start) ? String(start + 1) : String(member.start ?? '').trim();
  const endText = Number.isFinite(end) ? String(end) : String(member.end ?? '').trim();
  const range = startText && endText ? `${startText}..${endText}` : startText || endText;
  const strand = String(member.strand || '').trim();
  return range && strand ? `${range} (${strand})` : range;
};

const buildStandaloneMatchMemberRows = (orthogroup) => {
  const members = Array.isArray(orthogroup?.members) ? orthogroup.members : [];
  const orthogroupId = firstStandaloneText(orthogroup?.id, orthogroup?.orthogroupId, orthogroup?.orthogroup_id);
  const displayName = firstStandaloneText(orthogroup?.display_name, orthogroup?.displayName, orthogroup?.name);
  return members
    .map((member) => ({
      featureSvgId: standaloneMemberFeatureSvgId(member),
      feature_svg_id: standaloneMemberFeatureSvgId(member),
      orthogroupId,
      orthogroup_id: orthogroupId,
      displayName,
      display_name: displayName,
      record: firstStandaloneText(member?.recordId, member?.record_id),
      coordinates: standaloneMemberLocationText(member),
      proteinId: resolveStandaloneDisplayProteinId(null, member),
      productOrNote: firstStandaloneText(member?.product, member?.note)
    }))
    .filter((row) => row.record || row.coordinates || row.proteinId || row.productOrNote || row.featureSvgId);
};

const standaloneMemberCopyText = (memberRows) => {
  if (!Array.isArray(memberRows) || memberRows.length === 0) return '';
  return [
    'Record\tCoordinates (+/-)\tProtein ID\tProduct / note',
    ...memberRows.map((row) => [row.record, row.coordinates, row.proteinId, row.productOrNote].join('\t'))
  ].join('\n');
};

const standaloneMatchKind = (element) => {
  const explicit = standaloneAttr(element, 'data-match-kind').toLowerCase();
  if (explicit === 'pairwise' || explicit === 'orthogroup' || explicit === 'collinear') return explicit;
  if (standaloneAttr(element, 'data-collinearity-block-id')) return 'collinear';
  if (standaloneAttr(element, 'data-orthogroup-id')) return 'orthogroup';
  return 'pairwise';
};

const STANDALONE_MATCH_TITLES = {
  pairwise: 'Pairwise match',
  orthogroup: 'Orthogroup match',
  collinear: 'Collinearity block'
};

const standaloneOrthogroupTitle = (orthogroupId, displayName) => {
  const id = String(orthogroupId || '').trim();
  const name = String(displayName || '').trim();
  if (id && name) return `${id}:${name}`;
  return id || name || STANDALONE_MATCH_TITLES.orthogroup;
};

const buildStandaloneMatchSection = (title, rows) => ({
  title,
  rows: rows.filter((row) => String(row?.[1] || '').trim())
});

const standaloneMemberFeatureSvgId = (member) => firstStandaloneText(member?.featureSvgId, member?.feature_svg_id);

const getStandaloneOrthogroupById = (orthogroups, orthogroupId) => {
  const id = firstStandaloneText(orthogroupId);
  if (!id) return null;
  return (Array.isArray(orthogroups) ? orthogroups : [])
    .find((entry) => firstStandaloneText(entry?.id, entry?.orthogroupId, entry?.orthogroup_id) === id) || null;
};

const standaloneGroupHasFeatureSvgId = (group, featureSvgId) => {
  const ids = splitStandaloneMetadataValues(featureSvgId);
  if (ids.length === 0) return false;
  return (Array.isArray(group?.members) ? group.members : [])
    .some((member) => ids.includes(standaloneMemberFeatureSvgId(member)));
};

const getStandaloneOrthogroupPayload = (orthogroups, orthogroupId, queryFeatureSvgId = '', subjectFeatureSvgId = '') => {
  const id = String(orthogroupId || '').trim();
  const groups = Array.isArray(orthogroups) ? orthogroups : [];
  const direct = getStandaloneOrthogroupById(groups, id);
  if (direct) return direct;
  return groups.find((entry) => (
    standaloneGroupHasFeatureSvgId(entry, queryFeatureSvgId) ||
    standaloneGroupHasFeatureSvgId(entry, subjectFeatureSvgId)
  )) || null;
};

const standaloneFeatureOrthogroupId = (feature) => firstStandaloneText(
  feature?.orthogroup_id,
  feature?.orthogroupId,
  feature?.orthogroup_member?.orthogroupId,
  feature?.orthogroup_member?.orthogroup_id,
  feature?.orthogroupMember?.orthogroupId,
  feature?.orthogroupMember?.orthogroup_id
);

const standaloneFeatureProduct = (feature) => firstStandaloneText(
  feature?.product,
  feature?.orthogroup_member?.product,
  feature?.orthogroupMember?.product,
  feature?.note,
  feature?.label,
  feature?.display_label,
  feature?.displayLabel
);

const standaloneFeatureToMember = (feature) => {
  if (!feature) return null;
  const member = feature?.orthogroup_member || feature?.orthogroupMember || {};
  return {
    recordId: firstStandaloneText(member?.recordId, member?.record_id, feature?.record_id),
    start: member?.start ?? feature?.start,
    end: member?.end ?? feature?.end,
    strand: firstStandaloneText(member?.strand, feature?.strand),
    proteinId: firstStandaloneText(member?.proteinId, member?.protein_id, feature?.protein_id, feature?.proteinId),
    sourceProteinId: firstStandaloneText(member?.sourceProteinId, member?.source_protein_id, feature?.source_protein_id, feature?.sourceProteinId),
    product: firstStandaloneText(member?.product, standaloneFeatureProduct(feature)),
    note: firstStandaloneText(member?.note, feature?.note),
    featureSvgId: firstStandaloneText(member?.featureSvgId, member?.feature_svg_id, feature?.svg_id)
  };
};

const buildStandaloneFallbackOrthogroup = ({ orthogroupId, queryFeature, subjectFeature, features }) => {
  const id = String(orthogroupId || '').trim();
  if (!id) return null;
  const matchingFeatures = (Array.isArray(features) ? features : [])
    .filter((feature) => standaloneFeatureOrthogroupId(feature) === id);
  const fallbackFeatures = matchingFeatures.length
    ? matchingFeatures
    : [queryFeature, subjectFeature].filter(Boolean);
  if (!fallbackFeatures.length) return null;
  const members = fallbackFeatures.map(standaloneFeatureToMember).filter(Boolean);
  const firstFeature = fallbackFeatures[0] || {};
  const recordCoverageFallback = new Set(
    fallbackFeatures
      .map((feature) => firstStandaloneText(feature?.record_id, feature?.recordId, feature?.record_idx, feature?.recordIndex))
      .filter(Boolean)
  ).size;
  return {
    id,
    name: standaloneFeatureProduct(firstFeature),
    member_count: firstStandaloneText(firstFeature?.orthogroup_member_count, firstFeature?.orthogroupMemberCount, members.length),
    record_coverage_count: firstStandaloneText(firstFeature?.orthogroup_record_coverage, firstFeature?.orthogroupRecordCoverage, recordCoverageFallback),
    members
  };
};

const getStandaloneGroupMemberForFeatureSvgId = (group, featureSvgId) => {
  const id = firstStandaloneText(featureSvgId);
  if (!id) return null;
  const members = Array.isArray(group?.members) ? group.members : [];
  return members.find((member) => standaloneMemberFeatureSvgId(member) === id) || null;
};

const resolveStandaloneFeatureSectionProteinIds = ({
  feature,
  featureSvgIds,
  featuresById,
  group,
  fallbackProteinIds,
  locusId,
  displayName
}) => {
  const values = [];
  const addFeatureProteinId = (candidateFeature, member = null) => {
    const text = resolveStandaloneDisplayProteinId(candidateFeature, member, '');
    addUniqueStandaloneDisplayText(values, text);
  };

  const ids = splitStandaloneMetadataValues(featureSvgIds);
  ids.forEach((featureSvgId) => {
    const candidateFeature = featuresById?.get?.(featureSvgId) || null;
    const member = getStandaloneGroupMemberForFeatureSvgId(group, featureSvgId);
    addFeatureProteinId(candidateFeature, member);
  });
  if (feature) {
    addFeatureProteinId(feature, ids.length === 1 ? getStandaloneGroupMemberForFeatureSvgId(group, ids[0]) : null);
  }
  if (values.length === 0) {
    splitStandaloneMetadataValues(locusId).forEach((value) => addUniqueStandaloneDisplayText(values, value));
  }
  if (values.length === 0) {
    splitStandaloneMetadataValues(displayName).forEach((value) => addUniqueStandaloneDisplayText(values, value));
  }
  if (values.length === 0) {
    splitStandaloneMetadataValues(fallbackProteinIds).forEach((value) => addUniqueStandaloneDisplayText(values, value));
  }
  return values.join('; ');
};

const firstStandaloneNonInternalDisplayText = (...values) => {
  for (const value of values) {
    const text = firstStandaloneText(value);
    if (text && !isStandaloneInternalDisplayId(text)) return text;
  }
  return '';
};

const standaloneFeatureLocationText = (feature) => {
  const direct = firstStandaloneText(feature?.location);
  if (direct && direct !== '..') return direct;
  const built = buildStandaloneFeatureLocation(feature);
  return built && built !== '..' ? built : '';
};

const standaloneFeatureRowLocusId = (feature, fallbackLocusId) => firstStandaloneText(
  feature?.locusTag,
  feature?.locus_tag,
  firstQualifierValue(feature, 'locus_tag'),
  feature?.geneId,
  feature?.gene_id,
  firstQualifierValue(feature, 'gene_id'),
  fallbackLocusId
);

const standaloneFeatureRowDisplayName = (feature, fallbackDisplayName) => firstStandaloneText(
  fallbackDisplayName,
  feature?.displayLabel,
  feature?.display_label,
  feature?.label,
  feature?.gene,
  firstQualifierValue(feature, 'gene'),
  feature?.locus_tag,
  feature?.locusTag,
  firstQualifierValue(feature, 'locus_tag'),
  standaloneFeatureProduct(feature)
);

const buildStandaloneFeatureListRows = ({
  featureSvgIds,
  featuresById,
  group,
  recordId,
  interval,
  proteinId,
  locusId,
  displayName
}) => {
  const featureIds = uniqueStandaloneMetadataValues(featureSvgIds);
  const proteinIds = splitStandaloneMetadataValues(proteinId);
  const locusIds = splitStandaloneMetadataValues(locusId);
  const displayNames = splitStandaloneMetadataValues(displayName);
  const count = Math.max(featureIds.length, proteinIds.length, locusIds.length, displayNames.length);
  if (count === 0) return [];

  return Array.from({ length: count }, (_unused, index) => {
    const svgId = featureIds[index] || '';
    const feature = svgId ? featuresById?.get?.(svgId) || null : null;
    const member = svgId ? getStandaloneGroupMemberForFeatureSvgId(group, svgId) : null;
    const fallbackProteinId = firstStandaloneNonInternalDisplayText(
      locusIds[index],
      displayNames[index],
      proteinIds[index]
    );
    const resolvedProteinId = resolveStandaloneDisplayProteinId(feature, member, '');
    const displayProteinId = firstStandaloneText(
      isStandaloneInternalDisplayId(resolvedProteinId) ? '' : resolvedProteinId,
      fallbackProteinId,
      resolvedProteinId,
      proteinIds[index]
    );
    const rowRecord = firstStandaloneText(feature?.record_id, feature?.recordId, member?.recordId, member?.record_id, recordId);
    const rowLocation = firstStandaloneText(
      standaloneFeatureLocationText(feature),
      standaloneMemberLocationText(member),
      count === 1 ? interval : ''
    );
    const rowLocusId = standaloneFeatureRowLocusId(feature, locusIds[index]);
    const rowDisplayName = standaloneFeatureRowDisplayName(feature, displayNames[index]);
    const product = standaloneFeatureProduct(feature);
    const label = firstStandaloneText(displayProteinId, rowDisplayName, rowLocusId, product, svgId, `Feature ${index + 1}`);
    const copyText = [
      rowRecord,
      rowLocation,
      displayProteinId,
      rowLocusId,
      rowDisplayName,
      product
    ].join('\t');
    return {
      key: `${svgId || 'feature'}-${index}`,
      svg_id: svgId,
      svgId,
      can_open: Boolean(feature?.svg_id),
      canOpen: Boolean(feature?.svg_id),
      label,
      record: rowRecord,
      location: rowLocation,
      protein_id: displayProteinId,
      proteinId: displayProteinId,
      locus_id: rowLocusId,
      locusId: rowLocusId,
      display_name: rowDisplayName,
      displayName: rowDisplayName,
      product,
      type: firstStandaloneText(feature?.type),
      copy_text: copyText,
      copyText
    };
  }).filter((row) => (
    row.svg_id ||
    row.record ||
    row.location ||
    row.protein_id ||
    row.locus_id ||
    row.display_name ||
    row.product
  ));
};

const resolveStandaloneBlockMemberLabels = ({
  group,
  featureSvgIds,
  featuresById
}) => {
  if (!group) return '';
  const values = [];
  uniqueStandaloneMetadataValues(featureSvgIds).forEach((featureSvgId) => {
    const feature = featuresById?.get?.(featureSvgId) || null;
    const member = group ? getStandaloneGroupMemberForFeatureSvgId(group, featureSvgId) : null;
    if (!member) return;
    addUniqueStandaloneDisplayText(values, resolveStandaloneDisplayProteinId(feature, member, ''));
  });
  return values.join('; ');
};

const buildStandaloneOrthogroupDetailRows = ({
  orthogroupId,
  displayName,
  description,
  memberCount,
  recordCoverage
}) => {
  const rows = [];
  addStandaloneMatchRow(rows, 'Orthogroup ID', orthogroupId);
  addStandaloneMatchRow(rows, 'Display name', displayName);
  addStandaloneMatchRow(rows, 'Description', description);
  addStandaloneMatchRow(rows, 'Members', memberCount);
  addStandaloneMatchRow(rows, 'Record coverage', recordCoverage);
  return rows;
};

const buildStandaloneBlockOrthogroups = ({
  orthogroupIds,
  orthogroups,
  featuresById,
  queryFeatureSvgId,
  subjectFeatureSvgId
}) => orthogroupIds.map((orthogroupId) => {
  const group = getStandaloneOrthogroupById(orthogroups, orthogroupId);
  const displayName = firstStandaloneText(group?.display_name, group?.displayName, group?.name);
  const description = firstStandaloneText(group?.description);
  const memberCount = firstStandaloneText(
    group?.member_count,
    group?.memberCount,
    Array.isArray(group?.members) ? group.members.length : ''
  );
  const recordCoverage = firstStandaloneText(group?.record_coverage_count, group?.recordCoverage);
  const memberRows = buildStandaloneMatchMemberRows(group);
  return {
    id: orthogroupId,
    display_name: displayName,
    displayName,
    description,
    member_count: memberCount,
    memberCount,
    record_coverage: recordCoverage,
    recordCoverage,
    query_member: resolveStandaloneBlockMemberLabels({ group, featureSvgIds: queryFeatureSvgId, featuresById }),
    subject_member: resolveStandaloneBlockMemberLabels({ group, featureSvgIds: subjectFeatureSvgId, featuresById }),
    detail_rows: buildStandaloneOrthogroupDetailRows({
      orthogroupId,
      displayName,
      description,
      memberCount,
      recordCoverage
    }),
    member_rows: memberRows,
    member_copy_text: standaloneMemberCopyText(memberRows)
  };
});

const buildStandaloneMatchFeatureSection = ({
  title,
  feature,
  recordId,
  interval,
  proteinId,
  locusId,
  displayName,
  featureSvgIds,
  featuresById,
  group
}) => {
  const rows = [];
  const featureRows = buildStandaloneFeatureListRows({
    featureSvgIds,
    featuresById,
    group,
    recordId,
    interval,
    proteinId,
    locusId,
    displayName
  });
  const displayProteinIds = resolveStandaloneFeatureSectionProteinIds({
    feature,
    featureSvgIds,
    featuresById,
    group,
    fallbackProteinIds: proteinId,
    locusId,
    displayName
  });
  addStandaloneMatchRow(rows, 'Record', firstStandaloneText(feature?.record_id, recordId));
  addStandaloneMatchRow(rows, 'Location', firstStandaloneText(feature?.location, interval));
  addStandaloneMatchRow(rows, displayProteinIds.includes(';') ? 'Protein IDs' : 'Protein ID', displayProteinIds);
  addStandaloneMatchRow(rows, 'Type', feature?.type);
  addStandaloneMatchRow(rows, 'Locus ID', locusId);
  addStandaloneMatchRow(rows, 'Display name', displayName);
  return {
    ...buildStandaloneMatchSection(title, rows),
    feature_rows: featureRows,
    featureRows
  };
};

const buildStandaloneMatchPayloads = (svg, { features = [], orthogroups = [] } = {}) => {
  if (!svg) return [];
  const featuresById = new Map();
  (Array.isArray(features) ? features : []).forEach((feature) => {
    const svgId = String(feature?.svg_id || '').trim();
    if (svgId && !featuresById.has(svgId)) featuresById.set(svgId, feature);
  });
  return Array.from(svg.querySelectorAll(STANDALONE_MATCH_SELECTOR)).map((element, index) => {
    let id = standaloneAttr(element, 'data-gbdraw-pairwise-match-id');
    if (!id) {
      id = `pairwise_match_${index + 1}`;
      element.setAttribute('data-gbdraw-pairwise-match-id', id);
    }
    const matchKind = standaloneMatchKind(element);
    const orthogroupId = standaloneAttr(element, 'data-orthogroup-id');
    const orthogroupIds = uniqueStandaloneMetadataValues(orthogroupId);
    const blockId = standaloneAttr(element, 'data-collinearity-block-id');
    const queryFeatureSvgId = standaloneAttr(element, 'data-query-feature-svg-id');
    const subjectFeatureSvgId = standaloneAttr(element, 'data-subject-feature-svg-id');
    const queryFeature = featuresById.get(queryFeatureSvgId) || null;
    const subjectFeature = featuresById.get(subjectFeatureSvgId) || null;
    const orthogroup = matchKind === 'collinear'
      ? null
      : getStandaloneOrthogroupPayload(orthogroups, orthogroupId, queryFeatureSvgId, subjectFeatureSvgId) ||
        buildStandaloneFallbackOrthogroup({
          orthogroupId,
          queryFeature,
          subjectFeature,
          features
        });
    const blockOrthogroups = matchKind === 'collinear'
      ? buildStandaloneBlockOrthogroups({
        orthogroupIds,
        orthogroups,
        featuresById,
        queryFeatureSvgId,
        subjectFeatureSvgId
      })
      : [];
    const qInterval = standaloneIntervalText(standaloneAttr(element, 'data-qstart'), standaloneAttr(element, 'data-qend'));
    const sInterval = standaloneIntervalText(standaloneAttr(element, 'data-sstart'), standaloneAttr(element, 'data-send'));
    const summaryRows = [];
    addStandaloneMatchRow(summaryRows, 'Query record', firstStandaloneText(standaloneAttr(element, 'data-query-record-id'), standaloneAttr(element, 'data-query')));
    addStandaloneMatchRow(summaryRows, 'Subject record', firstStandaloneText(standaloneAttr(element, 'data-subject-record-id'), standaloneAttr(element, 'data-subject')));
    addStandaloneMatchRow(summaryRows, 'Query interval', qInterval);
    addStandaloneMatchRow(summaryRows, 'Subject interval', sInterval);
    addStandaloneMatchRow(summaryRows, 'Orientation', firstStandaloneText(standaloneAttr(element, 'data-collinearity-orientation'), standaloneAttr(element, 'data-orientation')));

    const alignmentRows = [];
    addStandaloneMatchRow(alignmentRows, 'Identity', standaloneAttr(element, 'data-identity'));
    addStandaloneMatchRow(alignmentRows, 'Alignment length', standaloneAttr(element, 'data-alignment-length'));
    addStandaloneMatchRow(alignmentRows, 'E-value', standaloneAttr(element, 'data-evalue'));
    addStandaloneMatchRow(alignmentRows, 'Bit score', standaloneAttr(element, 'data-bitscore'));
    addStandaloneMatchRow(alignmentRows, 'Mismatches', standaloneAttr(element, 'data-mismatches'));
    addStandaloneMatchRow(alignmentRows, 'Gap opens', standaloneAttr(element, 'data-gap-opens'));

    const orthogroupDisplayName = firstStandaloneText(
      orthogroup?.display_name,
      orthogroup?.displayName,
      orthogroup?.name,
      standaloneFeatureProduct(queryFeature),
      standaloneFeatureProduct(subjectFeature)
    );
    const orthogroupRows = [];
    if (matchKind !== 'collinear') {
      addStandaloneMatchRow(orthogroupRows, 'Orthogroup ID', orthogroupId);
      addStandaloneMatchRow(orthogroupRows, 'Display name', orthogroupDisplayName);
      addStandaloneMatchRow(orthogroupRows, 'Description', orthogroup?.description);
      addStandaloneMatchRow(orthogroupRows, 'Members', firstStandaloneText(orthogroup?.member_count, orthogroup?.memberCount));
      addStandaloneMatchRow(orthogroupRows, 'Record coverage', firstStandaloneText(orthogroup?.record_coverage_count, orthogroup?.recordCoverage));
    }
    const orthogroupMemberRows = buildStandaloneMatchMemberRows(orthogroup);

    const blockRows = [];
    addStandaloneMatchRow(blockRows, 'Block ID', blockId);
    addStandaloneMatchRow(blockRows, 'Kind', standaloneAttr(element, 'data-collinearity-block-kind'));
    addStandaloneMatchRow(blockRows, 'Orientation', standaloneAttr(element, 'data-collinearity-orientation'));
    addStandaloneMatchRow(blockRows, 'Color mode', standaloneAttr(element, 'data-collinearity-color-mode'));
    if (matchKind === 'collinear') {
      addStandaloneMatchRow(blockRows, 'Average identity', standaloneAttr(element, 'data-identity'));
      addStandaloneMatchRow(blockRows, 'Aligned length', standaloneAttr(element, 'data-alignment-length'));
    }
    addStandaloneMatchRow(blockRows, 'Block score', standaloneAttr(element, 'data-collinearity-block-score'));
    addStandaloneMatchRow(blockRows, 'Block e-value', standaloneAttr(element, 'data-collinearity-block-evalue'));
    addStandaloneMatchRow(blockRows, 'Anchor', [
      standaloneAttr(element, 'data-collinearity-anchor-index'),
      standaloneAttr(element, 'data-collinearity-anchor-count')
    ].filter(Boolean).join(' / '));

    if (matchKind === 'orthogroup') {
      const summarySection = {
        ...buildStandaloneMatchSection('Summary', orthogroupRows),
        member_rows: orthogroupMemberRows,
        member_copy_text: standaloneMemberCopyText(orthogroupMemberRows)
      };
      const hoverRows = [];
      addStandaloneMatchRow(hoverRows, 'Kind', matchKind);
      addStandaloneMatchRow(hoverRows, 'Orthogroup', orthogroupId);
      addStandaloneMatchRow(hoverRows, 'Display name', orthogroupDisplayName);
      addStandaloneMatchRow(hoverRows, 'Members', firstStandaloneText(orthogroup?.member_count, orthogroup?.memberCount));

      return {
        id,
        title: standaloneOrthogroupTitle(orthogroupId, orthogroupDisplayName),
        subtitle: '',
        match_kind: matchKind,
        orthogroup_id: orthogroupId,
        collinearity_block_id: blockId,
        fill: firstStandaloneText(element.getAttribute('fill'), '#94a3b8'),
        sections: (summarySection.rows.length || summarySection.member_rows.length) ? [summarySection] : [],
        hover_rows: hoverRows
      };
    }

    const sections = [buildStandaloneMatchSection('Summary', summaryRows)];
    if (matchKind !== 'collinear') {
      sections.push(buildStandaloneMatchSection('Alignment', alignmentRows));
    }
    if (orthogroupRows.length || matchKind === 'orthogroup') {
      sections.push({
        ...buildStandaloneMatchSection('Orthogroup', orthogroupRows),
        member_rows: orthogroupMemberRows,
        member_copy_text: standaloneMemberCopyText(orthogroupMemberRows)
      });
    }
    if (matchKind === 'collinear') {
      const blockOrthogroupRows = [];
      addStandaloneMatchRow(blockOrthogroupRows, 'Number of orthogroups covered', String(orthogroupIds.length));
      sections.push({
        ...buildStandaloneMatchSection('Orthogroups covered', blockOrthogroupRows),
        block_orthogroups: blockOrthogroups
      });
    }
    if (blockRows.length || matchKind === 'collinear') {
      sections.push(buildStandaloneMatchSection('Collinearity', blockRows));
    }
    sections.push(buildStandaloneMatchFeatureSection({
      title: 'Query',
      feature: queryFeature,
      recordId: standaloneAttr(element, 'data-query-record-id'),
      interval: qInterval,
      proteinId: standaloneAttr(element, 'data-query-protein-id'),
      locusId: standaloneAttr(element, 'data-query-locus-id'),
      displayName: standaloneAttr(element, 'data-query-display-name'),
      featureSvgIds: queryFeatureSvgId,
      featuresById,
      group: orthogroup
    }));
    sections.push(buildStandaloneMatchFeatureSection({
      title: 'Subject',
      feature: subjectFeature,
      recordId: standaloneAttr(element, 'data-subject-record-id'),
      interval: sInterval,
      proteinId: standaloneAttr(element, 'data-subject-protein-id'),
      locusId: standaloneAttr(element, 'data-subject-locus-id'),
      displayName: standaloneAttr(element, 'data-subject-display-name'),
      featureSvgIds: subjectFeatureSvgId,
      featuresById,
      group: orthogroup
    }));

    const findRow = (sectionTitle, rowLabel) => {
      const section = sections.find((entry) => entry.title === sectionTitle);
      const row = section?.rows.find((entry) => entry[0] === rowLabel);
      return row?.[1] || '';
    };
    const hoverRows = [];
    addStandaloneMatchRow(hoverRows, 'Kind', matchKind);
    addStandaloneMatchRow(hoverRows, 'Identity', findRow('Alignment', 'Identity') || findRow('Collinearity', 'Average identity'));
    addStandaloneMatchRow(hoverRows, 'Query', findRow('Summary', 'Query interval'));
    addStandaloneMatchRow(hoverRows, 'Subject', findRow('Summary', 'Subject interval'));
    if (matchKind === 'collinear') {
      addStandaloneMatchRow(hoverRows, 'Orthogroups', String(blockOrthogroups.length || orthogroupIds.length));
    } else {
      addStandaloneMatchRow(hoverRows, 'Orthogroup', orthogroupId);
    }
    addStandaloneMatchRow(hoverRows, 'Block', blockId);

    return {
      id,
      title: STANDALONE_MATCH_TITLES[matchKind] || STANDALONE_MATCH_TITLES.pairwise,
      subtitle: firstStandaloneText(blockId, orthogroupId, id),
      match_kind: matchKind,
      orthogroup_id: orthogroupId,
      collinearity_block_id: blockId,
      block_orthogroup_count: blockOrthogroups.length || (matchKind === 'collinear' ? orthogroupIds.length : 0),
      block_orthogroups: blockOrthogroups,
      fill: firstStandaloneText(element.getAttribute('fill'), '#94a3b8'),
      sections: sections.filter((section) => (
        section.rows.length > 0 ||
        (Array.isArray(section.feature_rows) && section.feature_rows.length > 0)
      )),
      hover_rows: hoverRows
    };
  });
};

const ensureSvgDefs = (svg) => {
  let defs = svg.querySelector('defs');
  if (!defs) {
    defs = document.createElementNS(SVG_NS, 'defs');
    svg.insertBefore(defs, svg.firstChild);
  }
  return defs;
};

const appendStandaloneFeatureGlowFilter = (
  defs,
  {
    id,
    color,
    opacity,
    blurStdDeviation,
    slope,
    extent
  }
) => {
  const existing = defs.querySelector(`#${CSS.escape(id)}`);
  if (existing?.parentNode) {
    existing.parentNode.removeChild(existing);
  }

  const filter = document.createElementNS(SVG_NS, 'filter');
  filter.setAttribute('id', id);
  filter.setAttribute('x', `-${extent}%`);
  filter.setAttribute('y', `-${extent}%`);
  filter.setAttribute('width', `${100 + (extent * 2)}%`);
  filter.setAttribute('height', `${100 + (extent * 2)}%`);
  filter.setAttribute('color-interpolation-filters', 'sRGB');

  const componentTransfer = document.createElementNS(SVG_NS, 'feComponentTransfer');
  componentTransfer.setAttribute('in', 'SourceGraphic');
  componentTransfer.setAttribute('result', 'gbdrawBrightenedFeature');
  ['R', 'G', 'B'].forEach((channel) => {
    const func = document.createElementNS(SVG_NS, `feFunc${channel}`);
    func.setAttribute('type', 'linear');
    func.setAttribute('slope', String(slope));
    componentTransfer.appendChild(func);
  });

  const blur = document.createElementNS(SVG_NS, 'feGaussianBlur');
  blur.setAttribute('in', 'SourceAlpha');
  blur.setAttribute('stdDeviation', String(blurStdDeviation));
  blur.setAttribute('result', 'gbdrawFeatureGlowBlur');

  const flood = document.createElementNS(SVG_NS, 'feFlood');
  flood.setAttribute('flood-color', color);
  flood.setAttribute('flood-opacity', String(opacity));
  flood.setAttribute('result', 'gbdrawFeatureGlowColor');

  const composite = document.createElementNS(SVG_NS, 'feComposite');
  composite.setAttribute('in', 'gbdrawFeatureGlowColor');
  composite.setAttribute('in2', 'gbdrawFeatureGlowBlur');
  composite.setAttribute('operator', 'in');
  composite.setAttribute('result', 'gbdrawFeatureGlow');

  const merge = document.createElementNS(SVG_NS, 'feMerge');
  ['gbdrawFeatureGlow', 'gbdrawBrightenedFeature'].forEach((resultName) => {
    const mergeNode = document.createElementNS(SVG_NS, 'feMergeNode');
    mergeNode.setAttribute('in', resultName);
    merge.appendChild(mergeNode);
  });

  filter.appendChild(componentTransfer);
  filter.appendChild(blur);
  filter.appendChild(flood);
  filter.appendChild(composite);
  filter.appendChild(merge);
  defs.appendChild(filter);
};

const ensureStandaloneFeatureGlowFilter = (svg) => {
  const defs = ensureSvgDefs(svg);
  appendStandaloneFeatureGlowFilter(defs, {
    id: INTERACTIVE_GLOW_FILTER_ID,
    color: '#2563eb',
    opacity: 0.85,
    blurStdDeviation: 3,
    slope: 1.2,
    extent: 35
  });
  appendStandaloneFeatureGlowFilter(defs, {
    id: INTERACTIVE_MATCH_GLOW_FILTER_ID,
    color: '#fbbf24',
    opacity: 0.32,
    blurStdDeviation: 1.5,
    slope: 1.04,
    extent: 25
  });
};

const removeExistingStandaloneInteractivityAssets = (svg) => {
  [
    INTERACTIVE_METADATA_ID,
    INTERACTIVE_STYLE_ID,
    INTERACTIVE_SCRIPT_ID,
    INTERACTIVE_GLOW_FILTER_ID,
    INTERACTIVE_MATCH_GLOW_FILTER_ID,
    'gbdraw-viewport-controls',
    'gbdraw-feature-search-controls',
    'gbdraw-feature-popup',
    'gbdraw-feature-hover-popup'
  ].forEach((id) => {
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

const parseSvgViewBox = (value) => {
  const parts = String(value || '')
    .trim()
    .split(/[\s,]+/)
    .map((part) => Number(part));
  if (parts.length < 4 || parts.some((part) => !Number.isFinite(part))) return null;
  if (parts[2] <= 0 || parts[3] <= 0) return null;
  return {
    x: parts[0],
    y: parts[1],
    width: parts[2],
    height: parts[3]
  };
};

const formatSvgViewBox = (rect) =>
  [rect.x, rect.y, rect.width, rect.height].map((value) => String(value)).join(' ');

const resolveStandaloneOriginalGeometry = (svg) => {
  const dataViewBox = parseSvgViewBox(svg.getAttribute('data-gbdraw-original-viewbox'));
  const currentViewBox = parseSvgViewBox(svg.getAttribute('viewBox'));
  const width = parseFloat(svg.getAttribute('width'));
  const height = parseFloat(svg.getAttribute('height'));
  const fallbackViewBox =
    Number.isFinite(width) && Number.isFinite(height) && width > 0 && height > 0
      ? { x: 0, y: 0, width, height }
      : null;
  const viewBox = dataViewBox || currentViewBox || fallbackViewBox || { x: 0, y: 0, width: 900, height: 650 };
  const originalWidth = svg.getAttribute('data-gbdraw-original-width') ||
    svg.getAttribute('width') ||
    `${viewBox.width}px`;
  const originalHeight = svg.getAttribute('data-gbdraw-original-height') ||
    svg.getAttribute('height') ||
    `${viewBox.height}px`;
  return {
    viewBox,
    originalWidth,
    originalHeight
  };
};

const applyStandaloneViewportRoot = (svg) => {
  svg.setAttribute('width', '100vw');
  svg.setAttribute('height', '100vh');
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  svg.style.setProperty('width', '100vw');
  svg.style.setProperty('height', '100vh');
  svg.style.setProperty('display', 'block');
  svg.style.setProperty('background', '#ffffff');
};

export const enrichSvgWithStandaloneInteractivity = (svg, options = {}) => {
  if (!svg) return false;

  const context = normalizeStandaloneContext(options);
  const originalGeometry = resolveStandaloneOriginalGeometry(svg);
  const normalizedPopupMode = normalizeStandalonePopupMode(context.popupMode);
  removeExistingStandaloneInteractivityAssets(svg);
  const features = buildStandaloneFeaturePayloads(svg, {
    ...context,
    popupMode: normalizedPopupMode
  });
  const orthogroups = buildStandaloneOrthogroupPayloads(features, context);
  const matches = buildStandaloneMatchPayloads(svg, { features, orthogroups });

  const featureIds = new Set(features.map((feature) => feature.svg_id));
  svg.querySelectorAll(FEATURE_SELECTOR).forEach((element) => {
    const id = getElementFeatureId(element);
    if (!featureIds.has(id)) return;
    element.setAttribute('data-gbdraw-interactive-feature', 'true');
    addClassToken(element, 'gbdraw-interactive-feature');
  });
  const matchIds = new Set(matches.map((match) => match.id));
  svg.querySelectorAll(STANDALONE_MATCH_SELECTOR).forEach((element) => {
    const id = standaloneAttr(element, 'data-gbdraw-pairwise-match-id');
    if (!matchIds.has(id)) return;
    element.setAttribute('data-gbdraw-interactive-match', 'true');
    addClassToken(element, 'gbdraw-interactive-pairwise-match');
  });

  const metadata = document.createElementNS(SVG_NS, 'metadata');
  metadata.setAttribute('id', INTERACTIVE_METADATA_ID);
  metadata.setAttribute('data-schema', 'gbdraw-interactive-feature-popup-v1');
  metadata.setAttribute('data-popup-mode', normalizedPopupMode);
  metadata.textContent = JSON.stringify({
    schema: 'gbdraw-interactive-feature-popup-v1',
    popup_mode: normalizedPopupMode,
    features,
    orthogroups,
    matches
  });

  const style = document.createElementNS(SVG_NS, 'style');
  style.setAttribute('id', INTERACTIVE_STYLE_ID);
  style.setAttribute('type', 'text/css');
  style.textContent = STANDALONE_INTERACTIVE_STYLE;

  const script = document.createElementNS(SVG_NS, 'script');
  script.setAttribute('id', INTERACTIVE_SCRIPT_ID);
  script.setAttribute('type', 'application/ecmascript');
  script.textContent = STANDALONE_INTERACTIVE_SCRIPT;

  ensureStandaloneFeatureGlowFilter(svg);
  svg.appendChild(metadata);
  svg.appendChild(style);
  svg.appendChild(script);
  svg.setAttribute('data-gbdraw-original-viewbox', formatSvgViewBox(originalGeometry.viewBox));
  svg.setAttribute('data-gbdraw-original-width', originalGeometry.originalWidth);
  svg.setAttribute('data-gbdraw-original-height', originalGeometry.originalHeight);
  svg.setAttribute('data-gbdraw-interactive-svg', 'true');
  applyStandaloneViewportRoot(svg);
  return true;
};
