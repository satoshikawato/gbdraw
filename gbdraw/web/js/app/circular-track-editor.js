import {
  BUILTIN_CIRCULAR_TRACK_IDS,
  createDefaultCircularAnalysisTrack,
  createDefaultCustomCircularTrack,
  getCircularTrackById,
  getDeclaredCircularFeatureTypes,
  nextCircularAnalysisTrackId,
  nextCircularCustomTrackId,
  normalizeCircularTracks,
  syncLegacyCircularTrackControls
} from '../utils/circular-tracks.js';

const { computed, reactive, ref } = window.Vue;

const cloneTracks = (tracks) => JSON.parse(JSON.stringify(Array.isArray(tracks) ? tracks : []));

const normalizeFeatureType = (value) => String(value || '').trim();

const normalizeStrandMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['all', 'positive', 'negative'].includes(normalized) ? normalized : 'all';
};

const normalizeAnalysisMetric = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['content', 'skew'].includes(normalized) ? normalized : 'content';
};

const normalizeNewCircularTrackKind = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return normalized === 'analysis' ? 'analysis' : 'custom';
};

const normalizeDinucleotide = (value, fallback = 'AT') => {
  const normalized = String(value || '')
    .trim()
    .toUpperCase()
    .replace(/\s+/g, '');
  if (/^[A-Z]{2}$/.test(normalized)) return normalized;
  return String(fallback || 'AT').trim().toUpperCase();
};

const sanitizeDinucleotideDraft = (value) =>
  String(value || '')
    .toUpperCase()
    .replace(/[^A-Z]/g, '')
    .slice(0, 2);

const toPlacementInputValue = (placement, field) => {
  const raw = String(placement?.[field] ?? '').trim();
  if (!raw) return '';
  return raw.endsWith('px') ? raw.slice(0, -2) : raw;
};

export const createCircularTrackEditor = ({ state }) => {
  const { circularTracks, featureKeys, form, adv } = state;

  const pendingFeatureTypeByTrackId = reactive({});
  const pendingRuleQualifierByTrackId = reactive({});
  const pendingRulePatternByTrackId = reactive({});
  const pendingDinucleotideByTrackId = reactive({});
  const newCircularTrackKind = ref('custom');

  const commitTracks = (nextTracks) => {
    const normalizedTracks = normalizeCircularTracks(nextTracks, { adv, form });
    circularTracks.value = normalizedTracks;
    syncLegacyCircularTrackControls(normalizedTracks, { adv, form });
    return normalizedTracks;
  };

  const getTrackIndex = (trackId) =>
    circularTracks.value.findIndex((track) => String(track?.id || '') === String(trackId || ''));

  const updateTrack = (trackId, updater) => {
    const nextTracks = cloneTracks(circularTracks.value);
    const trackIndex = nextTracks.findIndex((track) => String(track?.id || '') === String(trackId || ''));
    if (trackIndex < 0) return;
    const currentTrack = nextTracks[trackIndex];
    const updatedTrack = updater(currentTrack);
    if (!updatedTrack) return;
    nextTracks[trackIndex] = updatedTrack;
    commitTracks(nextTracks);
  };

  const circularFeatureShapeTypes = computed(() => getDeclaredCircularFeatureTypes(circularTracks.value));

  const getCircularTrackKindLabel = (track) => {
    const kind = String(track?.kind || '').trim().toLowerCase();
    if (kind === 'features') return 'Features';
    if (kind === 'gc_content') return 'GC Content';
    if (kind === 'gc_skew') return 'GC Skew';
    if (kind === 'analysis') return 'Analysis';
    if (kind === 'custom') return 'Custom';
    return kind || 'Track';
  };

  const isBuiltinCircularTrack = (track) => BUILTIN_CIRCULAR_TRACK_IDS.includes(String(track?.id || ''));

  const canMoveCircularTrackUp = (index) => Number(index) > 0;

  const canMoveCircularTrackDown = (index) => {
    const normalizedIndex = Number(index);
    return Number.isInteger(normalizedIndex) && normalizedIndex >= 0 && normalizedIndex < circularTracks.value.length - 1;
  };

  const moveCircularTrack = (index, direction) => {
    const currentIndex = Number(index);
    const nextIndex = currentIndex + Number(direction);
    if (!Number.isInteger(currentIndex) || !Number.isInteger(nextIndex)) return;
    if (currentIndex < 0 || currentIndex >= circularTracks.value.length) return;
    if (nextIndex < 0 || nextIndex >= circularTracks.value.length) return;
    const nextTracks = cloneTracks(circularTracks.value);
    const [moved] = nextTracks.splice(currentIndex, 1);
    nextTracks.splice(nextIndex, 0, moved);
    commitTracks(nextTracks);
  };

  const addCircularCustomTrack = () => {
    const featureType = circularFeatureShapeTypes.value[0] || 'CDS';
    const trackId = nextCircularCustomTrackId(circularTracks.value);
    const nextTracks = cloneTracks(circularTracks.value);
    nextTracks.push(createDefaultCustomCircularTrack(trackId, featureType));
    commitTracks(nextTracks);
    pendingFeatureTypeByTrackId[trackId] = featureType;
    pendingRuleQualifierByTrackId[trackId] = 'product';
    pendingRulePatternByTrackId[trackId] = '';
  };

  const addCircularAnalysisTrack = () => {
    const trackId = nextCircularAnalysisTrackId(circularTracks.value);
    const nextTracks = cloneTracks(circularTracks.value);
    nextTracks.push(createDefaultCircularAnalysisTrack(trackId, 'content', 'AT'));
    commitTracks(nextTracks);
  };

  const addCircularTrack = () => {
    if (normalizeNewCircularTrackKind(newCircularTrackKind.value) === 'analysis') {
      addCircularAnalysisTrack();
      return;
    }
    addCircularCustomTrack();
  };

  const deleteCircularTrack = (trackId) => {
    const track = getCircularTrackById(circularTracks.value, trackId);
    if (!track || isBuiltinCircularTrack(track)) return;
    const nextTracks = cloneTracks(circularTracks.value).filter((entry) => String(entry?.id || '') !== String(trackId || ''));
    commitTracks(nextTracks);
    delete pendingFeatureTypeByTrackId[trackId];
    delete pendingRuleQualifierByTrackId[trackId];
    delete pendingRulePatternByTrackId[trackId];
    delete pendingDinucleotideByTrackId[trackId];
  };

  const setCircularTrackShow = (trackId, show) => {
    updateTrack(trackId, (track) => ({ ...track, show: show === true }));
  };

  const getCircularTrackPlacementValue = (trackId, field) => {
    const track = getCircularTrackById(circularTracks.value, trackId);
    return toPlacementInputValue(track?.placement, field);
  };

  const setCircularTrackPlacementValue = (trackId, field, rawValue) => {
    const normalizedField = String(field || '').trim();
    if (!['width', 'radius'].includes(normalizedField)) return;
    const trimmed = String(rawValue ?? '').trim();

    updateTrack(trackId, (track) => {
      const placement = track?.placement && typeof track.placement === 'object'
        ? { ...track.placement }
        : {};
      if (!trimmed) {
        delete placement[normalizedField];
      } else {
        const numeric = Number(trimmed);
        if (!Number.isFinite(numeric) || numeric <= 0) {
          delete placement[normalizedField];
        } else if (normalizedField === 'width') {
          placement.width = `${numeric}px`;
        } else {
          placement.radius = String(numeric);
        }
      }
      return {
        ...track,
        placement: Object.keys(placement).length > 0 ? placement : null
      };
    });
  };

  const setCircularTrackCaption = (trackId, value) => {
    const nextCaption = String(value ?? '');
    updateTrack(trackId, (track) => {
      const kind = String(track?.kind || '').trim().toLowerCase();
      if (kind !== 'custom' && kind !== 'analysis') {
        return track;
      }
      if (kind === 'custom' && !nextCaption.trim()) {
        return track;
      }
      return {
        ...track,
        params: {
          ...(track.params || {}),
          caption: nextCaption.trimStart()
        }
      };
    });
  };

  const setCircularTrackAnalysisMetric = (trackId, value) => {
    updateTrack(trackId, (track) => {
      if (String(track?.kind || '').trim().toLowerCase() !== 'analysis') {
        return track;
      }
      const nextMetric = normalizeAnalysisMetric(value);
      return {
        ...track,
        params: {
          ...(track.params || {}),
          metric: nextMetric,
          dinucleotide: normalizeDinucleotide(track?.params?.dinucleotide, 'AT')
        }
      };
    });
  };

  const getCircularTrackDinucleotideInputValue = (trackId) => {
    if (typeof pendingDinucleotideByTrackId[trackId] === 'string') {
      return pendingDinucleotideByTrackId[trackId];
    }
    const track = getCircularTrackById(circularTracks.value, trackId);
    return String(track?.params?.dinucleotide || 'AT');
  };

  const setCircularTrackDinucleotide = (trackId, value) => {
    const draft = sanitizeDinucleotideDraft(value);
    pendingDinucleotideByTrackId[trackId] = draft;
    if (draft.length !== 2) {
      return;
    }
    updateTrack(trackId, (track) => {
      if (String(track?.kind || '').trim().toLowerCase() !== 'analysis') {
        return track;
      }
      const metric = normalizeAnalysisMetric(track?.params?.metric);
      return {
        ...track,
        params: {
          ...(track.params || {}),
          metric,
          dinucleotide: draft
        }
      };
    });
  };

  const finalizeCircularTrackDinucleotide = (trackId) => {
    const draft = pendingDinucleotideByTrackId[trackId];
    if (typeof draft !== 'string') {
      return;
    }
    if (draft.length === 2) {
      delete pendingDinucleotideByTrackId[trackId];
      return;
    }
    delete pendingDinucleotideByTrackId[trackId];
  };

  const setCircularTrackStrandMode = (trackId, value) => {
    updateTrack(trackId, (track) => {
      if (String(track?.kind || '').trim().toLowerCase() !== 'custom') {
        return track;
      }
      return {
        ...track,
        params: {
          ...(track.params || {}),
          strand_mode: normalizeStrandMode(value)
        }
      };
    });
  };

  const getPendingCircularTrackFeatureType = (trackId) => {
    const track = getCircularTrackById(circularTracks.value, trackId);
    const featureTypes = Array.isArray(track?.params?.feature_types) ? track.params.feature_types : [];
    const fallback = featureTypes[0] || featureKeys[0] || 'CDS';
    return String(pendingFeatureTypeByTrackId[trackId] || fallback);
  };

  const setPendingCircularTrackFeatureType = (trackId, value) => {
    pendingFeatureTypeByTrackId[trackId] = normalizeFeatureType(value);
  };

  const addCircularTrackFeatureType = (trackId) => {
    const nextFeatureType = normalizeFeatureType(getPendingCircularTrackFeatureType(trackId));
    if (!nextFeatureType) return;
    updateTrack(trackId, (track) => {
      const featureTypes = Array.isArray(track?.params?.feature_types) ? [...track.params.feature_types] : [];
      if (!featureTypes.includes(nextFeatureType)) {
        featureTypes.push(nextFeatureType);
      }
      return {
        ...track,
        params: {
          ...(track.params || {}),
          feature_types: featureTypes
        }
      };
    });
  };

  const removeCircularTrackFeatureType = (trackId, featureType) => {
    const normalizedFeatureType = normalizeFeatureType(featureType);
    if (!normalizedFeatureType) return;
    updateTrack(trackId, (track) => {
      const featureTypes = Array.isArray(track?.params?.feature_types)
        ? track.params.feature_types.filter((item) => normalizeFeatureType(item) !== normalizedFeatureType)
        : [];
      if (featureTypes.length === 0) {
        return track;
      }
      return {
        ...track,
        params: {
          ...(track.params || {}),
          feature_types: featureTypes
        }
      };
    });
  };

  const getPendingCircularTrackRuleQualifier = (trackId) =>
    String(pendingRuleQualifierByTrackId[trackId] || 'product');

  const setPendingCircularTrackRuleQualifier = (trackId, value) => {
    pendingRuleQualifierByTrackId[trackId] = String(value || '').trim();
  };

  const getPendingCircularTrackRulePattern = (trackId) => String(pendingRulePatternByTrackId[trackId] || '');

  const setPendingCircularTrackRulePattern = (trackId, value) => {
    pendingRulePatternByTrackId[trackId] = String(value ?? '');
  };

  const addCircularTrackRule = (trackId) => {
    const qualifier = String(getPendingCircularTrackRuleQualifier(trackId) || '').trim();
    if (!qualifier) return;
    const pattern = String(getPendingCircularTrackRulePattern(trackId) || '');
    updateTrack(trackId, (track) => {
      if (String(track?.kind || '').trim().toLowerCase() !== 'custom') {
        return track;
      }
      const rules = Array.isArray(track?.params?.rules) ? [...track.params.rules] : [];
      rules.push({ qualifier, pattern });
      return {
        ...track,
        params: {
          ...(track.params || {}),
          rules
        }
      };
    });
    pendingRuleQualifierByTrackId[trackId] = qualifier;
    pendingRulePatternByTrackId[trackId] = '';
  };

  const updateCircularTrackRule = (trackId, ruleIndex, field, value) => {
    const normalizedField = String(field || '').trim();
    if (!['qualifier', 'pattern'].includes(normalizedField)) return;
    if (normalizedField === 'qualifier' && !String(value || '').trim()) return;
    updateTrack(trackId, (track) => {
      if (String(track?.kind || '').trim().toLowerCase() !== 'custom') {
        return track;
      }
      const rules = Array.isArray(track?.params?.rules) ? [...track.params.rules] : [];
      if (ruleIndex < 0 || ruleIndex >= rules.length) return track;
      rules[ruleIndex] = {
        ...rules[ruleIndex],
        [normalizedField]: normalizedField === 'qualifier' ? String(value || '').trim() : String(value ?? '')
      };
      return {
        ...track,
        params: {
          ...(track.params || {}),
          rules
        }
      };
    });
  };

  const removeCircularTrackRule = (trackId, ruleIndex) => {
    updateTrack(trackId, (track) => {
      if (String(track?.kind || '').trim().toLowerCase() !== 'custom') {
        return track;
      }
      const rules = Array.isArray(track?.params?.rules) ? [...track.params.rules] : [];
      if (ruleIndex < 0 || ruleIndex >= rules.length) return track;
      rules.splice(ruleIndex, 1);
      return {
        ...track,
        params: {
          ...(track.params || {}),
          rules
        }
      };
    });
  };

  commitTracks(circularTracks.value);

  return {
    newCircularTrackKind,
    circularFeatureShapeTypes,
    getCircularTrackKindLabel,
    isBuiltinCircularTrack,
    canMoveCircularTrackUp,
    canMoveCircularTrackDown,
    moveCircularTrack,
    addCircularTrack,
    addCircularCustomTrack,
    addCircularAnalysisTrack,
    deleteCircularTrack,
    setCircularTrackShow,
    getCircularTrackPlacementValue,
    setCircularTrackPlacementValue,
    setCircularTrackCaption,
    setCircularTrackAnalysisMetric,
    getCircularTrackDinucleotideInputValue,
    setCircularTrackDinucleotide,
    finalizeCircularTrackDinucleotide,
    setCircularTrackStrandMode,
    getPendingCircularTrackFeatureType,
    setPendingCircularTrackFeatureType,
    addCircularTrackFeatureType,
    removeCircularTrackFeatureType,
    getPendingCircularTrackRuleQualifier,
    setPendingCircularTrackRuleQualifier,
    getPendingCircularTrackRulePattern,
    setPendingCircularTrackRulePattern,
    addCircularTrackRule,
    updateCircularTrackRule,
    removeCircularTrackRule
  };
};
