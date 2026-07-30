import { MODE_PROFILE_DATA } from './mode-profiles.generated.js';

const MODE_NAMES = Object.freeze(['circular', 'linear']);
const COMPARISON_STATE_FIELDS = Object.freeze([
  'min_bitscore',
  'evalue',
  'identity',
  'alignment_length'
]);
const PROFILE_MANAGED_ADV_FIELDS = Object.freeze([
  ...COMPARISON_STATE_FIELDS,
  'axis_stroke_color'
]);
const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

const normalizeMode = (mode) => {
  const normalized = String(mode || '').trim().toLowerCase();
  if (!MODE_NAMES.includes(normalized)) {
    throw new TypeError(`Unsupported diagram mode: ${String(mode)}`);
  }
  return normalized;
};

const formatEvalue = (value) => {
  const numeric = Number(value);
  if (!Number.isFinite(numeric)) return String(value);
  if (numeric === 0) return '0';
  return numeric
    .toExponential()
    .replace(/\.?0+e/, 'e')
    .replace(/e\+?(-?)0+(\d+)/, 'e$1$2');
};

const valuesEquivalent = (left, right) => {
  const leftBlank = left === null || left === undefined ||
    (typeof left === 'string' && left.trim() === '');
  const rightBlank = right === null || right === undefined ||
    (typeof right === 'string' && right.trim() === '');
  if (leftBlank || rightBlank) return leftBlank && rightBlank;

  const leftNumber = Number(left);
  const rightNumber = Number(right);
  if (Number.isFinite(leftNumber) && Number.isFinite(rightNumber)) {
    const tolerance = Number.EPSILON *
      Math.max(1, Math.abs(leftNumber), Math.abs(rightNumber)) * 8;
    return Math.abs(leftNumber - rightNumber) <= tolerance;
  }
  return String(left ?? '').trim() === String(right ?? '').trim();
};

export const MODE_PROFILE_VERSION = MODE_PROFILE_DATA.version;
export const MODE_PROFILE_STATE_SCHEMA = 1;
export const MODE_DEFAULT_FEATURE_TYPES = Object.freeze([...MODE_PROFILE_DATA.featureTypes]);

export const modeProfile = (mode) => MODE_PROFILE_DATA.modes[normalizeMode(mode)];

export const trackDefaultsForMode = (mode) => {
  const tracks = modeProfile(mode).tracks;
  return {
    gc: Boolean(tracks.gc),
    skew: Boolean(tracks.skew)
  };
};

export const comparisonFiltersForMode = (mode) => {
  const comparison = modeProfile(mode).comparison;
  return {
    bitscore: comparison.bitscore,
    evalue: formatEvalue(comparison.evalue),
    identity: comparison.identity,
    alignment_length: comparison.alignmentLength
  };
};

export const comparisonStateForMode = (mode) => {
  const filters = comparisonFiltersForMode(mode);
  return {
    min_bitscore: filters.bitscore,
    evalue: filters.evalue,
    identity: filters.identity,
    alignment_length: filters.alignment_length
  };
};

export const comparisonProfileDefault = (mode, field) => {
  const defaults = comparisonStateForMode(mode);
  if (!Object.prototype.hasOwnProperty.call(defaults, field)) {
    throw new TypeError(`Unsupported comparison field: ${String(field)}`);
  }
  return defaults[field];
};

export const managedAdvStateForMode = (mode) => ({
  ...comparisonStateForMode(mode),
  axis_stroke_color: modeProfile(mode).linearAxisColor
});

export const effectiveLinearAxisColor = ({
  axisColor = null,
  rulerOnAxis = false,
  managed = false
} = {}) => {
  const profile = modeProfile('linear');
  const normalized = axisColor === null || axisColor === undefined
    ? ''
    : String(axisColor).trim();
  if (!managed && normalized) return normalized;
  return (
    rulerOnAxis
      ? profile.linearRulerAxisColor
      : profile.linearAxisColor
  ) || normalized || null;
};

const readManagedAdvState = (source) => Object.fromEntries(
  PROFILE_MANAGED_ADV_FIELDS.map((field) => [field, source?.[field] ?? null])
);

const writeManagedAdvState = (target, values) => {
  PROFILE_MANAGED_ADV_FIELDS.forEach((field) => {
    target[field] = values[field];
  });
};

const managedFlagsFor = (values, defaults) => Object.fromEntries(
  PROFILE_MANAGED_ADV_FIELDS.map((field) => [
    field,
    valuesEquivalent(values[field], defaults[field])
  ])
);

export const createModeProfileStateManager = (initialMode, initialState) => {
  const snapshots = new Map();
  let activeMode = normalizeMode(initialMode);
  let activeBaseline = null;
  let activeManaged = null;
  let stateSource = initialState;

  const install = (mode, snapshot, target = null) => {
    activeMode = normalizeMode(mode);
    activeBaseline = { ...snapshot.values };
    activeManaged = { ...snapshot.managed };
    snapshots.set(activeMode, {
      values: { ...snapshot.values },
      managed: { ...snapshot.managed }
    });
    if (target) writeManagedAdvState(target, snapshot.values);
  };

  const defaultSnapshot = (mode) => {
    const values = managedAdvStateForMode(mode);
    return {
      values,
      managed: Object.fromEntries(
        PROFILE_MANAGED_ADV_FIELDS.map((field) => [field, true])
      )
    };
  };

  const detectEdits = (source) => {
    if (!activeBaseline || !activeManaged) return;
    PROFILE_MANAGED_ADV_FIELDS.forEach((field) => {
      if (
        activeManaged[field] &&
        !valuesEquivalent(source?.[field], activeBaseline[field])
      ) {
        activeManaged[field] = false;
      }
    });
  };

  const reset = (mode, source = null) => {
    const normalizedMode = normalizeMode(mode);
    const defaults = managedAdvStateForMode(normalizedMode);
    const values = source ? readManagedAdvState(source) : defaults;
    if (source) stateSource = source;
    snapshots.clear();
    install(normalizedMode, {
      values,
      managed: managedFlagsFor(values, defaults)
    });
  };

  const invalidate = (mode) => {
    snapshots.clear();
    activeMode = normalizeMode(mode);
    activeBaseline = null;
    activeManaged = null;
  };

  const transition = (source, previousMode, nextMode) => {
    const previous = normalizeMode(previousMode);
    const next = normalizeMode(nextMode);
    if (previous === next) return [];
    stateSource = source;

    if (activeMode === previous && activeBaseline && activeManaged) {
      detectEdits(source);
    } else {
      activeMode = previous;
      activeBaseline = readManagedAdvState(source);
      activeManaged = Object.fromEntries(
        PROFILE_MANAGED_ADV_FIELDS.map((field) => [field, false])
      );
    }
    snapshots.set(previous, {
      values: readManagedAdvState(source),
      managed: { ...activeManaged }
    });

    const targetSnapshot = snapshots.get(next) || defaultSnapshot(next);
    const before = readManagedAdvState(source);
    install(next, targetSnapshot, source);
    return PROFILE_MANAGED_ADV_FIELDS.filter(
      (field) => !valuesEquivalent(before[field], targetSnapshot.values[field])
    );
  };

  const isManaged = (source, field) => {
    if (!PROFILE_MANAGED_ADV_FIELDS.includes(field)) return false;
    if (!activeBaseline || !activeManaged) return false;
    stateSource = source;
    detectEdits(source);
    return activeManaged[field] === true;
  };

  const snapshotCurrentMode = () => {
    if (!activeBaseline || !activeManaged) {
      const values = readManagedAdvState(stateSource);
      const defaults = managedAdvStateForMode(activeMode);
      activeBaseline = { ...values };
      activeManaged = managedFlagsFor(values, defaults);
    } else {
      detectEdits(stateSource);
    }
    snapshots.set(activeMode, {
      values: readManagedAdvState(stateSource),
      managed: { ...activeManaged }
    });
  };

  const exportState = () => {
    snapshotCurrentMode();
    return {
      schema: MODE_PROFILE_STATE_SCHEMA,
      activeMode,
      profiles: Object.fromEntries(
        MODE_NAMES.map((mode) => {
          const snapshot = snapshots.get(mode) || defaultSnapshot(mode);
          return [
            mode,
            {
              values: { ...snapshot.values },
              managed: { ...snapshot.managed }
            }
          ];
        })
      )
    };
  };

  const normalizeImportedSnapshot = (mode, candidate) => {
    const defaults = defaultSnapshot(mode);
    if (candidate === undefined || candidate === null) return defaults;
    if (!isPlainObject(candidate)) {
      throw new TypeError(`Invalid ${mode} mode profile snapshot.`);
    }
    const valuesSource = candidate.values;
    const managedSource = candidate.managed;
    if (valuesSource !== undefined && !isPlainObject(valuesSource)) {
      throw new TypeError(`Invalid ${mode} mode profile values.`);
    }
    if (managedSource !== undefined && !isPlainObject(managedSource)) {
      throw new TypeError(`Invalid ${mode} mode profile managed flags.`);
    }
    const values = { ...defaults.values };
    const managed = { ...defaults.managed };
    PROFILE_MANAGED_ADV_FIELDS.forEach((field) => {
      if (Object.prototype.hasOwnProperty.call(valuesSource || {}, field)) {
        values[field] = valuesSource[field] ?? null;
      }
      if (Object.prototype.hasOwnProperty.call(managedSource || {}, field)) {
        if (typeof managedSource[field] !== 'boolean') {
          throw new TypeError(`Invalid managed flag for ${mode}.${field}.`);
        }
        managed[field] = managedSource[field];
      }
    });
    return { values, managed };
  };

  const importState = (payload, mode = null, target = stateSource) => {
    const nextMode = normalizeMode(mode || payload?.activeMode || activeMode);
    stateSource = target || stateSource;
    snapshots.clear();

    if (payload === null || payload === undefined) {
      MODE_NAMES.forEach((profileMode) => {
        snapshots.set(profileMode, defaultSnapshot(profileMode));
      });
      const defaults = managedAdvStateForMode(nextMode);
      const values = readManagedAdvState(stateSource);
      snapshots.set(nextMode, {
        values,
        managed: managedFlagsFor(values, defaults)
      });
    } else {
      if (!isPlainObject(payload)) {
        throw new TypeError('Invalid mode profile state.');
      }
      if (payload.schema !== MODE_PROFILE_STATE_SCHEMA) {
        throw new TypeError(`Unsupported mode profile state schema: ${String(payload.schema)}.`);
      }
      if (!isPlainObject(payload.profiles)) {
        throw new TypeError('Invalid mode profile state profiles.');
      }
      MODE_NAMES.forEach((profileMode) => {
        snapshots.set(
          profileMode,
          normalizeImportedSnapshot(profileMode, payload.profiles[profileMode])
        );
      });
    }

    install(nextMode, snapshots.get(nextMode) || defaultSnapshot(nextMode), stateSource);
  };

  reset(activeMode, initialState);
  return { exportState, importState, invalidate, isManaged, reset, transition };
};

export { MODE_PROFILE_DATA };
