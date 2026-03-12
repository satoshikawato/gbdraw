const DEFAULT_FEATURE_TYPES = ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'repeat_region'];
const BUILTIN_TRACK_IDS = ['features', 'gc_content', 'gc_skew'];
const ANALYSIS_METRICS = ['content', 'skew'];

const normalizeFeatureTypes = (value, fallback = []) => {
  const source = Array.isArray(value) && value.length > 0 ? value : fallback;
  const out = [];
  const seen = new Set();
  source.forEach((item) => {
    const featureType = String(item || '').trim();
    if (!featureType || seen.has(featureType)) return;
    seen.add(featureType);
    out.push(featureType);
  });
  return out;
};

const normalizePlacement = (value) => {
  if (!value || typeof value !== 'object' || Array.isArray(value)) return null;
  const placement = {};
  ['radius', 'inner_radius', 'outer_radius', 'width'].forEach((key) => {
    const raw = String(value[key] ?? '').trim();
    if (!raw) return;
    placement[key] = raw;
  });
  if (value.z !== undefined && value.z !== null && value.z !== '') {
    const numericZ = Number(value.z);
    if (Number.isInteger(numericZ)) {
      placement.z = numericZ;
    }
  }
  return Object.keys(placement).length > 0 ? placement : null;
};

const normalizeAnalysisMetric = (value, fallback = 'content') => {
  const metric = String(value || '').trim().toLowerCase();
  return ANALYSIS_METRICS.includes(metric) ? metric : fallback;
};

const normalizeDinucleotide = (value, fallback = 'AT') => {
  const dinucleotide = String(value || '')
    .trim()
    .toUpperCase()
    .replace(/\s+/g, '');
  if (/^[A-Z]{2}$/.test(dinucleotide)) return dinucleotide;
  return String(fallback || 'AT').trim().toUpperCase();
};

const toPxString = (value) => {
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0 ? `${numeric}px` : null;
};

const toRatioString = (value) => {
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0 ? String(numeric) : null;
};

const fromPxString = (value) => {
  const raw = String(value ?? '').trim();
  if (!raw) return null;
  if (!raw.endsWith('px')) return null;
  const numeric = Number(raw.slice(0, -2));
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

const fromRatioString = (value) => {
  const raw = String(value ?? '').trim();
  if (!raw) return null;
  const numeric = Number(raw.endsWith('%') ? raw.slice(0, -1) / 100 : raw);
  return Number.isFinite(numeric) && numeric > 0 ? numeric : null;
};

const buildFeaturesPlacementFromLegacy = (adv) => {
  const width = toPxString(adv?.feature_width_circular);
  return width ? { width } : null;
};

const buildGcPlacementFromLegacy = (widthValue, radiusValue) => {
  const placement = {};
  const width = toPxString(widthValue);
  const radius = toRatioString(radiusValue);
  if (radius) placement.radius = radius;
  if (width) placement.width = width;
  return Object.keys(placement).length > 0 ? placement : null;
};

export const createDefaultCircularTracksFromLegacy = ({ adv = {}, form = {} } = {}) => {
  const featureTypes = normalizeFeatureTypes(adv?.features, DEFAULT_FEATURE_TYPES);
  return [
    {
      id: 'features',
      kind: 'features',
      show: true,
      placement: buildFeaturesPlacementFromLegacy(adv),
      params: {
        feature_types: featureTypes
      }
    },
    {
      id: 'gc_content',
      kind: 'gc_content',
      show: form?.suppress_gc !== true,
      placement: buildGcPlacementFromLegacy(adv?.gc_content_width_circular, adv?.gc_content_radius_circular),
      params: null
    },
    {
      id: 'gc_skew',
      kind: 'gc_skew',
      show: form?.suppress_skew !== true,
      placement: buildGcPlacementFromLegacy(adv?.gc_skew_width_circular, adv?.gc_skew_radius_circular),
      params: null
    }
  ];
};

export const createDefaultCustomCircularTrack = (id, featureType = 'CDS') => ({
  id,
  kind: 'custom',
  show: true,
  placement: null,
  params: {
    caption: 'Custom Track',
    feature_types: normalizeFeatureTypes([featureType], ['CDS']),
    strand_mode: 'all',
    rules: [],
    match_all: true
  }
});

export const createDefaultCircularAnalysisTrack = (id, metric = 'content', dinucleotide = 'AT') => {
  const normalizedMetric = normalizeAnalysisMetric(metric, 'content');
  const normalizedDinucleotide = normalizeDinucleotide(dinucleotide, 'AT');
  return {
    id,
    kind: 'analysis',
    show: true,
    placement: null,
    params: {
      caption: '',
      metric: normalizedMetric,
      dinucleotide: normalizedDinucleotide
    }
  };
};

export const normalizeCircularTracks = (value, { adv = {}, form = {} } = {}) => {
  const defaults = createDefaultCircularTracksFromLegacy({ adv, form });
  const fallbackFeatureTypes = defaults[0].params.feature_types;
  if (!Array.isArray(value)) {
    return defaults;
  }

  const normalizedEntries = [];
  const seenIds = new Set();
  const seenBuiltins = new Set();

  value.forEach((entry) => {
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return;
    const kind = String(entry.kind || '').trim().toLowerCase();
    if (!['features', 'gc_content', 'gc_skew', 'custom', 'analysis'].includes(kind)) return;

    let id = String(entry.id || '').trim();
    if (BUILTIN_TRACK_IDS.includes(kind)) {
      id = kind;
      if (seenBuiltins.has(kind)) return;
      seenBuiltins.add(kind);
    } else {
      if (!id) return;
    }
    if (!id || seenIds.has(id)) return;
    seenIds.add(id);

    const show = entry.show === undefined ? true : entry.show === true;
    const placement = normalizePlacement(entry.placement);

    if (kind === 'features') {
      normalizedEntries.push({
        id,
        kind,
        show,
        placement,
        params: {
          feature_types: normalizeFeatureTypes(entry.params?.feature_types, fallbackFeatureTypes)
        }
      });
      return;
    }

    if (kind === 'custom') {
      const featureTypes = normalizeFeatureTypes(entry.params?.feature_types, []);
      const caption = String(entry.params?.caption || '').trim();
      if (!caption || featureTypes.length === 0) return;
      const strandMode = ['all', 'positive', 'negative'].includes(String(entry.params?.strand_mode || '').trim().toLowerCase())
        ? String(entry.params?.strand_mode || '').trim().toLowerCase()
        : 'all';
      const rules = Array.isArray(entry.params?.rules)
        ? entry.params.rules
            .map((rule) => ({
              qualifier: String(rule?.qualifier || '').trim(),
              pattern: String(rule?.pattern || '').trim()
            }))
            .filter((rule) => rule.qualifier !== '')
        : [];
      normalizedEntries.push({
        id,
        kind,
        show,
        placement,
        params: {
          caption,
          feature_types: featureTypes,
          strand_mode: strandMode,
          rules,
          match_all: true
        }
      });
      return;
    }

    if (kind === 'analysis') {
      const metric = normalizeAnalysisMetric(entry.params?.metric, 'content');
      const dinucleotide = normalizeDinucleotide(entry.params?.dinucleotide, 'AT');
      const caption = String(entry.params?.caption ?? '').trim();
      normalizedEntries.push({
        id,
        kind,
        show,
        placement,
        params: {
          caption,
          metric,
          dinucleotide
        }
      });
      return;
    }

    normalizedEntries.push({
      id,
      kind,
      show,
      placement,
      params: null
    });
  });

  const explicitKinds = new Set(normalizedEntries.filter((entry) => entry.kind !== 'custom').map((entry) => entry.kind));
  const ordered = [];
  if (!explicitKinds.has('features')) {
    ordered.push(defaults[0]);
  }
  ordered.push(...normalizedEntries);
  if (!explicitKinds.has('gc_content')) {
    ordered.push(defaults[1]);
  }
  if (!explicitKinds.has('gc_skew')) {
    ordered.push(defaults[2]);
  }
  return ordered;
};

export const syncLegacyCircularTrackControls = (tracks, { adv = {}, form = {} } = {}) => {
  const normalizedTracks = normalizeCircularTracks(tracks, { adv, form });
  const featuresTrack = getCircularTrackById(normalizedTracks, 'features');
  const gcContentTrack = getCircularTrackById(normalizedTracks, 'gc_content');
  const gcSkewTrack = getCircularTrackById(normalizedTracks, 'gc_skew');

  if (adv && typeof adv === 'object') {
    adv.feature_width_circular = fromPxString(featuresTrack?.placement?.width);
    adv.gc_content_width_circular = fromPxString(gcContentTrack?.placement?.width);
    adv.gc_content_radius_circular = fromRatioString(gcContentTrack?.placement?.radius);
    adv.gc_skew_width_circular = fromPxString(gcSkewTrack?.placement?.width);
    adv.gc_skew_radius_circular = fromRatioString(gcSkewTrack?.placement?.radius);
  }

  if (form && typeof form === 'object') {
    form.suppress_gc = gcContentTrack?.show === false;
    form.suppress_skew = gcSkewTrack?.show === false;
  }

  return normalizedTracks;
};

export const getCircularTrackById = (tracks, id) => {
  if (!Array.isArray(tracks)) return null;
  return tracks.find((track) => String(track?.id || '') === String(id || '')) || null;
};

export const getDeclaredCircularFeatureTypes = (tracks) => {
  const seen = new Set();
  const out = [];
  (Array.isArray(tracks) ? tracks : []).forEach((track) => {
    const featureTypes = Array.isArray(track?.params?.feature_types) ? track.params.feature_types : [];
    featureTypes.forEach((featureType) => {
      const normalized = String(featureType || '').trim();
      if (!normalized || seen.has(normalized)) return;
      seen.add(normalized);
      out.push(normalized);
    });
  });
  return out;
};

export const getEnabledCircularFeatureTypes = (tracks) => {
  const seen = new Set();
  const out = [];
  let sawFeaturesTrack = false;

  (Array.isArray(tracks) ? tracks : []).forEach((track) => {
    const kind = String(track?.kind || '').trim().toLowerCase();
    if (kind === 'features') {
      sawFeaturesTrack = true;
    }
    if (track?.show !== true) return;
    if (kind !== 'features' && kind !== 'custom') return;
    const featureTypes = Array.isArray(track?.params?.feature_types) ? track.params.feature_types : [];
    featureTypes.forEach((featureType) => {
      const normalized = String(featureType || '').trim();
      if (!normalized || seen.has(normalized)) return;
      seen.add(normalized);
      out.push(normalized);
    });
  });

  if (!sawFeaturesTrack) {
    DEFAULT_FEATURE_TYPES.forEach((featureType) => {
      if (seen.has(featureType)) return;
      seen.add(featureType);
      out.push(featureType);
    });
  }
  return out;
};

export const nextCircularCustomTrackId = (tracks) => {
  const used = new Set((Array.isArray(tracks) ? tracks : []).map((track) => String(track?.id || '')));
  let index = 1;
  while (used.has(`custom_${index}`)) {
    index += 1;
  }
  return `custom_${index}`;
};

export const nextCircularAnalysisTrackId = (tracks) => {
  const used = new Set((Array.isArray(tracks) ? tracks : []).map((track) => String(track?.id || '')));
  let index = 1;
  while (used.has(`analysis_${index}`)) {
    index += 1;
  }
  return `analysis_${index}`;
};

export const BUILTIN_CIRCULAR_TRACK_IDS = BUILTIN_TRACK_IDS;
export const DEFAULT_CIRCULAR_FEATURE_TYPES = DEFAULT_FEATURE_TYPES;
