const EMPTY_LAYOUT = Object.freeze({
  definitionMaxRadius: null,
  orderedVisibleTrackIds: [],
  tracksById: {}
});

const parseNumericAttr = (element, name) => {
  if (!element) return null;
  const raw = String(element.getAttribute(name) || '').trim();
  if (!raw) return null;
  const value = Number(raw);
  return Number.isFinite(value) ? value : null;
};

const normalizeTrackMetrics = (group) => {
  if (!group) return null;
  const trackId = String(group.getAttribute('data-track-id') || '').trim();
  if (!trackId) return null;

  const innerRadius = parseNumericAttr(group, 'data-track-inner-radius');
  const outerRadius = parseNumericAttr(group, 'data-track-outer-radius');
  const centerRadius = parseNumericAttr(group, 'data-track-center-radius');
  const centerFactor = parseNumericAttr(group, 'data-track-center-factor');
  const widthAttr = parseNumericAttr(group, 'data-track-width');

  if (innerRadius === null || outerRadius === null || centerRadius === null) {
    return null;
  }

  return {
    id: trackId,
    kind: String(group.getAttribute('data-track-kind') || '').trim(),
    metric: String(group.getAttribute('data-track-metric') || '').trim() || null,
    innerRadius,
    outerRadius,
    centerRadius,
    centerFactor,
    width: widthAttr !== null ? widthAttr : Math.max(0, outerRadius - innerRadius),
    gapPx: null,
    gapKind: null,
    gapTargetId: null
  };
};

export const extractCircularTrackLayoutMetricsFromSvg = (svgMarkup, { tracks = [] } = {}) => {
  if (!svgMarkup || typeof svgMarkup !== 'string') {
    return EMPTY_LAYOUT;
  }

  try {
    const parser = new DOMParser();
    const doc = parser.parseFromString(svgMarkup, 'image/svg+xml');
    const svg = doc.querySelector('svg');
    if (!svg) return EMPTY_LAYOUT;

    const tracksById = {};
    svg.querySelectorAll('g[data-track-id]').forEach((group) => {
      const metrics = normalizeTrackMetrics(group);
      if (!metrics || tracksById[metrics.id]) return;
      tracksById[metrics.id] = metrics;
    });

    const definitionGroup = Array.from(svg.querySelectorAll('g[data-definition-max-radius]')).find((group) => {
      const groupId = String(group.getAttribute('id') || '').trim();
      return groupId !== 'plot_title';
    }) || null;
    const definitionMaxRadius = parseNumericAttr(definitionGroup, 'data-definition-max-radius');

    const visibleTrackIds = Array.isArray(tracks)
      ? tracks
          .filter((track) => track?.show === true)
          .map((track) => String(track?.id || '').trim())
          .filter((trackId) => trackId && tracksById[trackId])
      : Object.keys(tracksById);

    const orderedVisibleTrackIds = [...visibleTrackIds].sort((leftId, rightId) => {
      const leftOuter = Number(tracksById[leftId]?.outerRadius ?? -Infinity);
      const rightOuter = Number(tracksById[rightId]?.outerRadius ?? -Infinity);
      return rightOuter - leftOuter;
    });

    orderedVisibleTrackIds.forEach((trackId, index) => {
      const currentTrack = tracksById[trackId];
      if (!currentTrack) return;

      const nextTrackId = orderedVisibleTrackIds[index + 1] || null;
      if (nextTrackId) {
        const nextTrack = tracksById[nextTrackId];
        if (!nextTrack) return;
        currentTrack.gapPx = currentTrack.innerRadius - nextTrack.outerRadius;
        currentTrack.gapKind = 'track';
        currentTrack.gapTargetId = nextTrackId;
        return;
      }

      if (definitionMaxRadius !== null) {
        currentTrack.gapPx = currentTrack.innerRadius - definitionMaxRadius;
        currentTrack.gapKind = 'definition';
      }
    });

    return {
      definitionMaxRadius,
      orderedVisibleTrackIds,
      tracksById
    };
  } catch (_error) {
    return EMPTY_LAYOUT;
  }
};
