import { normalizeSpecificRule } from '../app/specific-color-rules.js';

const isPlainObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const stableValue = (value) => {
  if (Array.isArray(value)) return value.map(stableValue);
  if (!isPlainObject(value)) return value;
  return Object.fromEntries(
    Object.keys(value).sort().map((key) => [key, stableValue(value[key])])
  );
};

const fnv1a64 = (text) => {
  let hash = 0xcbf29ce484222325n;
  for (let index = 0; index < text.length; index += 1) {
    const code = text.charCodeAt(index);
    hash ^= BigInt(code & 0xff);
    hash = BigInt.asUintN(64, hash * 0x100000001b3n);
    hash ^= BigInt(code >>> 8);
    hash = BigInt.asUintN(64, hash * 0x100000001b3n);
  }
  return hash.toString(16).padStart(16, '0');
};

const normalizedRules = (rules) => (Array.isArray(rules) ? rules : [])
  .map((rule) => {
    const normalized = normalizeSpecificRule(rule, { fromFile: false });
    return {
      feat: normalized.feat,
      qual: normalized.qual,
      val: normalized.val,
      color: normalized.color,
      cap: normalized.cap,
      ...(normalized.match ? { match: normalized.match } : {})
    };
  });

const normalizedColors = (colors) => Object.fromEntries(
  Object.entries(isPlainObject(colors) ? colors : {})
    .map(([key, value]) => [String(key), String(value ?? '').trim().toLowerCase()])
    .filter(([key, value]) => key && value)
    .sort(([left], [right]) => left.localeCompare(right))
);

export const normalizeStyleSnapshot = ({
  rules = [],
  appliedPaletteName = 'default',
  appliedPaletteColors = {}
} = {}) => Object.freeze({
  rules: Object.freeze(normalizedRules(rules).map(Object.freeze)),
  appliedPaletteName: String(appliedPaletteName || 'default'),
  appliedPaletteColors: Object.freeze(normalizedColors(appliedPaletteColors))
});

export const styleFingerprint = (snapshot) => (
  `sf1_${fnv1a64(JSON.stringify(stableValue(normalizeStyleSnapshot(snapshot))))}`
);

export const captureStyleSnapshot = (state) => {
  const semantic = normalizeStyleSnapshot({
    rules: state?.manualSpecificRules,
    appliedPaletteName: state?.appliedPaletteName?.value,
    appliedPaletteColors: state?.appliedPaletteColors?.value
  });
  const fingerprint = styleFingerprint(semantic);
  return Object.freeze({
    documentEpoch: Number(state?.documentEpoch?.value || 0),
    resultGenerationKey: Number(state?.resultGenerationKey?.value || 0),
    revision: Number(state?.semanticStyleRevision?.value || 0),
    fingerprint,
    semantic
  });
};

export const styleSnapshotIsCurrent = (state, snapshot) => Boolean(
  snapshot
  && Number(state?.documentEpoch?.value || 0) === snapshot.documentEpoch
  && Number(state?.resultGenerationKey?.value || 0) === snapshot.resultGenerationKey
  && Number(state?.semanticStyleRevision?.value || 0) === snapshot.revision
  && styleFingerprint({
    rules: state?.manualSpecificRules,
    appliedPaletteName: state?.appliedPaletteName?.value,
    appliedPaletteColors: state?.appliedPaletteColors?.value
  }) === snapshot.fingerprint
);

export const advanceStyleRevision = (
  state,
  { expected = null, resultFingerprints = null } = {}
) => {
  if (expected && !styleSnapshotIsCurrent(state, expected)) {
    throw new Error('The document style changed while this operation was being prepared.');
  }
  const fingerprint = styleFingerprint({
    rules: state?.manualSpecificRules,
    appliedPaletteName: state?.appliedPaletteName?.value,
    appliedPaletteColors: state?.appliedPaletteColors?.value
  });
  state.semanticStyleRevision.value += 1;
  state.semanticStyleFingerprint.value = fingerprint;
  if (resultFingerprints !== null && state.validatedStyleFingerprintByResultKey) {
    state.validatedStyleFingerprintByResultKey.value = Object.freeze({
      ...(resultFingerprints || {})
    });
  }
  return Object.freeze({
    revision: state.semanticStyleRevision.value,
    fingerprint
  });
};

export const advanceDocumentEpoch = (state) => {
  state.documentEpoch.value += 1;
  state.semanticStyleRevision.value = 0;
  const fingerprint = styleFingerprint({
    rules: state?.manualSpecificRules,
    appliedPaletteName: state?.appliedPaletteName?.value,
    appliedPaletteColors: state?.appliedPaletteColors?.value
  });
  state.semanticStyleFingerprint.value = fingerprint;
  state.validatedStyleFingerprintByResultKey.value = Object.freeze({});
  if (state.pendingFeatureFillPlan) state.pendingFeatureFillPlan.value = null;
  if (state.featureFillPlanStatus) state.featureFillPlanStatus.value = 'idle';
  if (state.featureFillPlanProgress) state.featureFillPlanProgress.value = null;
  if (state.pendingFeatureStrokePlan) state.pendingFeatureStrokePlan.value = null;
  if (state.featureStrokePlanStatus) state.featureStrokePlanStatus.value = 'idle';
  if (state.featureStrokePlanProgress) state.featureStrokePlanProgress.value = null;
  if (state.pendingFeatureLabelPlan) state.pendingFeatureLabelPlan.value = null;
  if (state.featureLabelPlanStatus) state.featureLabelPlanStatus.value = 'idle';
  if (state.featureLabelPlanProgress) state.featureLabelPlanProgress.value = null;
  return state.documentEpoch.value;
};
