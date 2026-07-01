import { resolveColorToHex } from '../color-utils.js';
import { getFeatureCaption, resolveDisplayProteinId } from '../feature-utils.js';
import {
  PAIRWISE_MATCH_SELECTOR,
  buildPairwiseMatchHoverRows,
  buildPairwiseMatchPayload
} from '../pairwise-match-popup.js';
import { buildFeatureSequenceFastas } from '../feature-sequence-fasta.js';

export const FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id';
export const FEATURE_SELECTOR = [
  `path[${FEATURE_ID_ATTRIBUTE}]`,
  `polygon[${FEATURE_ID_ATTRIBUTE}]`,
  `rect[${FEATURE_ID_ATTRIBUTE}]`,
  'path[id^="f"]',
  'polygon[id^="f"]',
  'rect[id^="f"]'
].join(', ');

const FEATURE_PART_SUFFIX_RE = /__part\d+$/;

export const normalizeFeatureIdentity = (value) =>
  String(value || '').trim().replace(FEATURE_PART_SUFFIX_RE, '');

export const getFeatureIdentity = (element) =>
  normalizeFeatureIdentity(
    element?.getAttribute?.(FEATURE_ID_ATTRIBUTE) ||
    element?.getAttribute?.('id') ||
    element?.id ||
    ''
  );

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
    isPanning,
    orthogroups,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides,
    extractedFeatures,
    featuresBySvgId,
    featureColorOverrides,
    featureVisibilityOverrides,
    svgContainer,
    clickedFeature,
    clickedFeaturePos,
    clickedPairwiseMatch,
    clickedPairwiseMatchPos,
    featurePopupSize,
    skipCaptureBaseConfig,
    adv
  } = state;
  const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
  const formatDuration = (ms) => `${ms.toFixed(1)}ms`;
  let delegatedFeatureHandlers = null;
  const hoverSummaryState = {
    element: null,
    timer: null,
    frame: null,
    visible: false,
    activeSvgId: '',
    lastEvent: null
  };

  const getOrthogroupIds = (value) =>
    Array.from(new Set(
      String(value || '')
        .split(';')
        .map((entry) => entry.trim())
        .filter(Boolean)
    ));

  const normalizeVisibilityMode = (value) => {
    const normalized = String(value || '').trim().toLowerCase();
    return ['on', 'off', 'suppress'].includes(normalized) ? normalized : 'default';
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

  const getQualifierFirstValue = (feat, key) => {
    const normalizedKey = String(key || '').trim().toLowerCase();
    if (!normalizedKey) return '';
    const directValue = feat?.[normalizedKey];
    const qualifierValue = feat?.qualifiers && typeof feat.qualifiers === 'object'
      ? feat.qualifiers[normalizedKey]
      : null;
    const values = normalizeStringArray(directValue || qualifierValue)
      .map((value) => String(value || '').trim())
      .filter(Boolean);
    return values[0] || '';
  };

  const getHoverSummaryPrimaryLabel = (feat) => (
    getQualifierFirstValue(feat, 'gene') ||
    getQualifierFirstValue(feat, 'locus_tag') ||
    getQualifierFirstValue(feat, 'product') ||
    getFeatureCaption(feat) ||
    ''
  );

  const formatFeatureLength = (feat) => {
    const startNumeric = Number(feat?.start);
    const endNumeric = Number(feat?.end);
    if (!Number.isFinite(startNumeric) || !Number.isFinite(endNumeric)) return '';
    const length = Math.max(0, Math.round(endNumeric - startNumeric));
    return `${length.toLocaleString()} bp`;
  };

  const createHoverSummaryElement = (tagName, className = '', text = '') => {
    const element = document.createElement(tagName);
    if (className) element.className = className;
    if (text !== '') element.textContent = text;
    return element;
  };

  const addHoverSummaryRow = (container, label, value, { clamp = false } = {}) => {
    const normalizedValue = String(value === null || value === undefined ? '' : value).trim();
    if (!normalizedValue) return;
    const row = createHoverSummaryElement('div', 'feature-hover-summary-row');
    row.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-label', label));
    row.appendChild(createHoverSummaryElement(
      'div',
      `feature-hover-summary-value${clamp ? ' is-clamped' : ''}`,
      normalizedValue
    ));
    container.appendChild(row);
  };

  const buildHoverSummaryRows = (feat, primaryLabel) => {
    const product = getQualifierFirstValue(feat, 'product');
    const gene = getQualifierFirstValue(feat, 'gene');
    const locusTag = getQualifierFirstValue(feat, 'locus_tag');
    const note = getQualifierFirstValue(feat, 'note');
    const locationText = buildFeatureLocation(feat);
    const effectiveCaption = String(getEffectiveLegendCaption?.(feat) || '').trim();
    const rows = [];

    if (gene && gene !== primaryLabel) rows.push(['Gene', gene]);
    if (locusTag && locusTag !== primaryLabel) rows.push(['Locus', locusTag]);
    if (product && product !== primaryLabel) rows.push(['Product', product, true]);
    if (note && note !== primaryLabel && note !== product) rows.push(['Note', note, true]);
    rows.push(['Length', formatFeatureLength(feat)]);
    rows.push(['Location', locationText]);
    rows.push(['Record', feat?.record_id || '']);
    if (feat?.orthogroupId) rows.push(['Orthogroup', feat.orthogroupId]);
    if (effectiveCaption && effectiveCaption !== primaryLabel) rows.push(['Legend', effectiveCaption]);
    return rows;
  };

  const buildOrthogroupDetailRows = (feat) => {
    const member = feat?.orthogroupMember || feat?.orthogroup_member || null;
    const proteinId = resolveDisplayProteinId(feat, member);
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
    const matchEl = target.closest(PAIRWISE_MATCH_SELECTOR);
    if (matchEl && svg.contains(matchEl)) return null;
    const featureEl = target.closest(FEATURE_SELECTOR);
    if (!featureEl || !svg.contains(featureEl)) return null;
    return featureEl;
  };

  const getPairwiseMatchTarget = (target, svg) => {
    if (!target || typeof target.closest !== 'function') return null;
    const matchEl = target.closest(PAIRWISE_MATCH_SELECTOR);
    if (!matchEl || !svg.contains(matchEl)) return null;
    return matchEl;
  };

  const getTopmostSvgTarget = (eventLike, svg, selector) => {
    if (!svg || !selector || !Number.isFinite(eventLike?.clientX) || !Number.isFinite(eventLike?.clientY)) {
      return null;
    }
    const stack = typeof document.elementsFromPoint === 'function'
      ? document.elementsFromPoint(eventLike.clientX, eventLike.clientY)
      : [];
    for (const element of stack) {
      if (!element || !svg.contains(element)) continue;
      const target = element.matches?.(selector)
        ? element
        : element.closest?.(selector);
      if (target && svg.contains(target)) return target;
    }
    return null;
  };

  const getFeatureClickTarget = (eventLike, svg) =>
    getTopmostSvgTarget(eventLike, svg, FEATURE_SELECTOR) || getFeatureTarget(eventLike?.target, svg);

  const getPairwiseMatchClickTarget = (eventLike, svg) =>
    getTopmostSvgTarget(eventLike, svg, PAIRWISE_MATCH_SELECTOR) || getPairwiseMatchTarget(eventLike?.target, svg);

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
    const { nucleotideFasta, aminoAcidFasta } = buildFeatureSequenceFastas(feat, {
      nucleotideSequence,
      aminoAcidSequence
    });

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
      nucleotideFasta,
      aminoAcidFasta,
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
      proteinId: feat.proteinId || feat.protein_id || '',
      sourceProteinId: feat.sourceProteinId || feat.source_protein_id || '',
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

    hideHoverSummary();
    if (clickedPairwiseMatch) clickedPairwiseMatch.value = null;
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

  const hoverSummaryIsAllowed = () => {
    if (clickedFeature.value) return false;
    if (clickedPairwiseMatch?.value) return false;
    if (isPanning?.value) return false;
    if (window.matchMedia && !window.matchMedia('(hover: hover) and (pointer: fine)').matches) {
      return false;
    }
    return true;
  };

  const ensureHoverSummaryElement = () => {
    if (hoverSummaryState.element?.isConnected) return hoverSummaryState.element;
    const element = createHoverSummaryElement('div', 'feature-hover-summary');
    element.hidden = true;
    element.setAttribute('role', 'tooltip');
    document.body.appendChild(element);
    hoverSummaryState.element = element;
    return element;
  };

  const positionHoverSummary = (eventLike = hoverSummaryState.lastEvent) => {
    const element = hoverSummaryState.element;
    if (!element || element.hidden || !eventLike) return;
    const margin = 12;
    const offset = 14;
    const clientX = Number.isFinite(eventLike.clientX) ? eventLike.clientX : window.innerWidth / 2;
    const clientY = Number.isFinite(eventLike.clientY) ? eventLike.clientY : window.innerHeight / 2;
    const rect = element.getBoundingClientRect();
    let x = clientX + offset;
    let y = clientY + offset;

    if (x + rect.width + margin > window.innerWidth) x = clientX - rect.width - offset;
    if (y + rect.height + margin > window.innerHeight) y = clientY - rect.height - offset;
    x = Math.min(Math.max(x, margin), Math.max(margin, window.innerWidth - rect.width - margin));
    y = Math.min(Math.max(y, margin), Math.max(margin, window.innerHeight - rect.height - margin));
    element.style.left = `${Math.round(x)}px`;
    element.style.top = `${Math.round(y)}px`;
  };

  const scheduleHoverSummaryPosition = (eventLike) => {
    hoverSummaryState.lastEvent = {
      clientX: Number.isFinite(eventLike?.clientX) ? eventLike.clientX : hoverSummaryState.lastEvent?.clientX,
      clientY: Number.isFinite(eventLike?.clientY) ? eventLike.clientY : hoverSummaryState.lastEvent?.clientY
    };
    if (!hoverSummaryState.visible || hoverSummaryState.frame) return;
    hoverSummaryState.frame = window.requestAnimationFrame(() => {
      hoverSummaryState.frame = null;
      positionHoverSummary();
    });
  };

  const renderHoverSummary = (feat, featureElement, eventLike) => {
    if (!feat || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    const element = ensureHoverSummaryElement();
    const primaryLabel = getHoverSummaryPrimaryLabel(feat);
    const featureType = String(feat?.type || 'Feature').trim() || 'Feature';
    const titleText = primaryLabel ? `${featureType}: ${primaryLabel}` : featureType;
    const locationText = buildFeatureLocation(feat);
    const color = resolveColorToHex(
      featureElement?.getAttribute?.('fill') || getFeatureColor(feat) || '#94a3b8'
    ) || '#94a3b8';

    element.replaceChildren();
    const title = createHoverSummaryElement('div', 'feature-hover-summary-title');
    const swatch = createHoverSummaryElement('div', 'feature-hover-summary-swatch');
    swatch.style.backgroundColor = color;
    const titleTextWrap = createHoverSummaryElement('div', 'feature-hover-summary-text');
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-heading', titleText));
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-subtitle', locationText));
    title.appendChild(swatch);
    title.appendChild(titleTextWrap);
    element.appendChild(title);

    buildHoverSummaryRows(feat, primaryLabel).forEach(([label, value, clamp]) => {
      addHoverSummaryRow(element, label, value, { clamp: Boolean(clamp) });
    });

    element.hidden = false;
    hoverSummaryState.visible = true;
    hoverSummaryState.activeSvgId = String(feat.svg_id || '').trim();
    scheduleHoverSummaryPosition(eventLike);
    positionHoverSummary(eventLike);
  };

  const scheduleHoverSummary = (feat, featureElement, eventLike) => {
    if (hoverSummaryState.timer) {
      window.clearTimeout(hoverSummaryState.timer);
      hoverSummaryState.timer = null;
    }
    if (!feat || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    scheduleHoverSummaryPosition(eventLike);
    const show = () => {
      hoverSummaryState.timer = null;
      renderHoverSummary(feat, featureElement, eventLike);
    };
    if (hoverSummaryState.visible) {
      show();
      return;
    }
    hoverSummaryState.timer = window.setTimeout(show, 180);
  };

  const renderMatchHoverSummary = (payload, eventLike) => {
    if (!payload || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    const element = ensureHoverSummaryElement();
    const color = resolveColorToHex(payload.fill || '#94a3b8') || '#94a3b8';

    element.replaceChildren();
    const title = createHoverSummaryElement('div', 'feature-hover-summary-title');
    const swatch = createHoverSummaryElement('div', 'feature-hover-summary-swatch');
    swatch.style.backgroundColor = color;
    const titleTextWrap = createHoverSummaryElement('div', 'feature-hover-summary-text');
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-heading', payload.title));
    titleTextWrap.appendChild(createHoverSummaryElement('div', 'feature-hover-summary-subtitle', payload.subtitle || payload.id || ''));
    title.appendChild(swatch);
    title.appendChild(titleTextWrap);
    element.appendChild(title);

    buildPairwiseMatchHoverRows(payload).forEach((row) => {
      addHoverSummaryRow(element, row.label, row.value, { clamp: row.label === 'Query' || row.label === 'Subject' });
    });

    element.hidden = false;
    hoverSummaryState.visible = true;
    hoverSummaryState.activeSvgId = String(payload.id || '').trim();
    scheduleHoverSummaryPosition(eventLike);
    positionHoverSummary(eventLike);
  };

  const scheduleMatchHoverSummary = (payload, eventLike) => {
    if (hoverSummaryState.timer) {
      window.clearTimeout(hoverSummaryState.timer);
      hoverSummaryState.timer = null;
    }
    if (!payload || !hoverSummaryIsAllowed()) {
      hideHoverSummary();
      return;
    }
    scheduleHoverSummaryPosition(eventLike);
    const show = () => {
      hoverSummaryState.timer = null;
      renderMatchHoverSummary(payload, eventLike);
    };
    if (hoverSummaryState.visible) {
      show();
      return;
    }
    hoverSummaryState.timer = window.setTimeout(show, 180);
  };

  function hideHoverSummary() {
    if (hoverSummaryState.timer) {
      window.clearTimeout(hoverSummaryState.timer);
      hoverSummaryState.timer = null;
    }
    if (hoverSummaryState.frame) {
      window.cancelAnimationFrame(hoverSummaryState.frame);
      hoverSummaryState.frame = null;
    }
    if (hoverSummaryState.element) {
      hoverSummaryState.element.hidden = true;
    }
    hoverSummaryState.visible = false;
    hoverSummaryState.activeSvgId = '';
    hoverSummaryState.lastEvent = null;
  }

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

  const buildMatchPayload = (matchElement, featureLookup) => buildPairwiseMatchPayload(matchElement, {
    featureLookup,
    orthogroups: orthogroups?.value,
    orthogroupNameOverrides,
    orthogroupDescriptionOverrides
  });

  const openPairwiseMatchPopup = (matchElement, eventLike, featureLookup) => {
    if (!matchElement || !clickedPairwiseMatch || !clickedPairwiseMatchPos) return null;
    const payload = buildMatchPayload(matchElement, featureLookup);
    if (!payload) return null;
    hideHoverSummary();
    clickedFeature.value = null;
    clickedPairwiseMatch.value = payload;
    const popupPosition = getPopupPosition(eventLike, 460, 520);
    clickedPairwiseMatchPos.x = popupPosition.x;
    clickedPairwiseMatchPos.y = popupPosition.y;
    return payload;
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
        if (mode === 'off' || mode === 'suppress') {
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
    const pairwiseMatchElements = Array.from(svg.querySelectorAll(PAIRWISE_MATCH_SELECTOR));
    const comparisonElementsByCollinearityBlockId = new Map();
    pairwiseMatchElements.forEach((element) => {
      if (element?.style) element.style.cursor = 'pointer';
      const blockId = String(element.getAttribute('data-collinearity-block-id') || '').trim();
      if (!blockId) return;
      if (!comparisonElementsByCollinearityBlockId.has(blockId)) {
        comparisonElementsByCollinearityBlockId.set(blockId, []);
      }
      comparisonElementsByCollinearityBlockId.get(blockId).push(element);
    });
    const indexDuration = getNow() - indexStartedAt;

    if (!delegatedFeatureHandlers) {
      const handlerState = {
        svg,
        pathsByIdMap,
        featureLookup,
        featureIdsByOrthogroupId,
        comparisonElementsByOrthogroupId,
        pairwiseMatchElements,
        comparisonElementsByCollinearityBlockId,
        activeHoverSvgId: null,
        activeHoverKey: '',
        activeMatchHoverElement: null,
        activeMatchHoverKey: '',
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

      const setCollinearityBlockHover = (blockId, highlight) => {
        const id = String(blockId || '').trim();
        if (!id) return;
        (handlerState.comparisonElementsByCollinearityBlockId.get(id) || []).forEach((element) => {
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

      const matchAttr = (element, name) => String(element?.getAttribute?.(name) || '').trim();

      const getMatchHoverKey = (matchElement) => {
        const blockId = matchAttr(matchElement, 'data-collinearity-block-id');
        if (blockId) return `collinearity:${blockId}`;
        const orthogroupId = matchAttr(matchElement, 'data-orthogroup-id');
        if (orthogroupId) return `orthogroup:${orthogroupId}`;
        return `match:${matchAttr(matchElement, 'data-gbdraw-pairwise-match-id') || matchAttr(matchElement, 'd')}`;
      };

      const setMatchHover = (matchElement, highlight) => {
        if (!matchElement) return;
        const blockId = matchAttr(matchElement, 'data-collinearity-block-id');
        const orthogroupId = matchAttr(matchElement, 'data-orthogroup-id');
        setHoverStyle(matchElement, highlight);
        if (blockId) {
          setCollinearityBlockHover(blockId, highlight);
          return;
        }
        if (orthogroupId) {
          setOrthogroupHover(orthogroupId, highlight);
        }
      };

      const clearActiveFeatureHover = () => {
        if (!handlerState.activeHoverSvgId) return;
        setHoverHighlight(handlerState.activeHoverSvgId, false);
        handlerState.activeHoverSvgId = null;
        handlerState.activeHoverKey = '';
      };

      const clearActiveMatchHover = () => {
        if (!handlerState.activeMatchHoverElement) return;
        setMatchHover(handlerState.activeMatchHoverElement, false);
        handlerState.activeMatchHoverElement = null;
        handlerState.activeMatchHoverKey = '';
      };

      const handleMouseOver = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        if (featureEl) {
          clearActiveMatchHover();
          const svgId = getFeatureIdentity(featureEl);
          if (!svgId) return;
          const hoverKey = getFeatureHoverKey(svgId);
          if (handlerState.activeHoverKey !== hoverKey && handlerState.activeHoverSvgId) {
            setHoverHighlight(handlerState.activeHoverSvgId, false);
          }
          if (handlerState.activeHoverKey !== hoverKey) {
            setHoverHighlight(svgId, true);
          }
          handlerState.activeHoverSvgId = svgId;
          handlerState.activeHoverKey = hoverKey;
          scheduleHoverSummary(handlerState.featureLookup.get(svgId), featureEl, e);
          return;
        }
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (!matchEl) return;
        clearActiveFeatureHover();
        const matchKey = getMatchHoverKey(matchEl);
        if (handlerState.activeMatchHoverKey !== matchKey && handlerState.activeMatchHoverElement) {
          setMatchHover(handlerState.activeMatchHoverElement, false);
        }
        if (handlerState.activeMatchHoverKey !== matchKey) {
          setMatchHover(matchEl, true);
        }
        handlerState.activeMatchHoverElement = matchEl;
        handlerState.activeMatchHoverKey = matchKey;
        scheduleMatchHoverSummary(buildMatchPayload(matchEl, handlerState.featureLookup), e);
      };

      const handleMouseMove = (e) => {
        if (clickedFeature.value || clickedPairwiseMatch?.value || isPanning?.value) {
          hideHoverSummary();
          return;
        }
        if (hoverSummaryState.visible || hoverSummaryState.timer) {
          scheduleHoverSummaryPosition(e);
        }
      };

      const handleMouseOut = (e) => {
        const featureEl = getFeatureTarget(e.target, svg);
        if (featureEl) {
          const svgId = getFeatureIdentity(featureEl);
          if (!svgId || handlerState.activeHoverSvgId !== svgId) return;
          const relatedFeature = getFeatureTarget(e.relatedTarget, svg);
          if (relatedFeature && getFeatureHoverKey(getFeatureIdentity(relatedFeature)) === handlerState.activeHoverKey) return;
          clearActiveFeatureHover();
          hideHoverSummary();
          return;
        }
        const matchEl = getPairwiseMatchTarget(e.target, svg);
        if (!matchEl || handlerState.activeMatchHoverElement !== matchEl) return;
        const relatedMatch = getPairwiseMatchTarget(e.relatedTarget, svg);
        if (relatedMatch && getMatchHoverKey(relatedMatch) === handlerState.activeMatchHoverKey) return;
        clearActiveMatchHover();
        hideHoverSummary();
      };

      const handleClick = (e) => {
        const featureEl = getFeatureClickTarget(e, svg);
        if (featureEl) {
          const svgId = getFeatureIdentity(featureEl);
          if (!svgId) return;
          e.stopPropagation();
          hideHoverSummary();
          const feat = handlerState.featureLookup.get(svgId);
          if (feat) {
            openFeatureEditorForFeature(feat, e);
          } else {
            console.log(`No feature found for svg_id: ${svgId}`);
          }
          return;
        }
        const matchEl = getPairwiseMatchClickTarget(e, svg);
        if (matchEl) {
          e.stopPropagation();
          e.preventDefault();
          openPairwiseMatchPopup(matchEl, e, handlerState.featureLookup);
        }
      };

      svg.addEventListener('mouseover', handleMouseOver);
      svg.addEventListener('mousemove', handleMouseMove);
      svg.addEventListener('mouseout', handleMouseOut);
      svg.addEventListener('click', handleClick);
      handlerState.cleanup = () => {
        svg.removeEventListener('mouseover', handleMouseOver);
        svg.removeEventListener('mousemove', handleMouseMove);
        svg.removeEventListener('mouseout', handleMouseOut);
        svg.removeEventListener('click', handleClick);
        clearActiveFeatureHover();
        clearActiveMatchHover();
        hideHoverSummary();
      };
      delegatedFeatureHandlers = handlerState;
    } else {
      delegatedFeatureHandlers.pathsByIdMap = pathsByIdMap;
      delegatedFeatureHandlers.featureLookup = featureLookup;
      delegatedFeatureHandlers.featureIdsByOrthogroupId = featureIdsByOrthogroupId;
      delegatedFeatureHandlers.comparisonElementsByOrthogroupId = comparisonElementsByOrthogroupId;
      delegatedFeatureHandlers.pairwiseMatchElements = pairwiseMatchElements;
      delegatedFeatureHandlers.comparisonElementsByCollinearityBlockId = comparisonElementsByCollinearityBlockId;
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
