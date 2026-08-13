const INSTANCE_HASH_PREFIX = 'fi1_';
const INSTANCE_HASH_DOMAIN = 'gbdraw-feature-instance-v1\0';
const BASE32_ALPHABET = 'abcdefghijklmnopqrstuvwxyz234567';

export const FEATURE_INSTANCE_HASH_QUALIFIER = '__gbdraw_instance_hash__';
export const FEATURE_INSTANCE_HASH_PATTERN = /^fi1_[a-z2-7]{26}$/;

const text = (value) => String(value ?? '').trim();

const appendUint32 = (target, value) => {
  target.push(
    (value >>> 24) & 0xff,
    (value >>> 16) & 0xff,
    (value >>> 8) & 0xff,
    value & 0xff
  );
};

const appendLengthFramed = (target, value, encoder) => {
  const bytes = encoder.encode(value);
  appendUint32(target, bytes.byteLength);
  target.push(...bytes);
};

const base32WithoutPadding = (bytes) => {
  let bits = 0;
  let value = 0;
  let encoded = '';
  bytes.forEach((byte) => {
    value = (value << 8) | byte;
    bits += 8;
    while (bits >= 5) {
      encoded += BASE32_ALPHABET[(value >>> (bits - 5)) & 31];
      bits -= 5;
      value &= (1 << bits) - 1;
    }
  });
  if (bits > 0) encoded += BASE32_ALPHABET[(value << (5 - bits)) & 31];
  return encoded;
};

export const isFeatureInstanceHash = (value) => (
  FEATURE_INSTANCE_HASH_PATTERN.test(text(value))
);

export const featureInstanceIdentityBytes = (recordKey, biologicalFeatureId) => {
  const normalizedRecordKey = text(recordKey);
  const normalizedFeatureId = text(biologicalFeatureId);
  if (!normalizedRecordKey || !normalizedFeatureId) return null;

  const encoder = new TextEncoder();
  const bytes = [...encoder.encode(INSTANCE_HASH_DOMAIN)];
  appendLengthFramed(bytes, normalizedRecordKey, encoder);
  appendLengthFramed(bytes, normalizedFeatureId, encoder);
  return new Uint8Array(bytes);
};

export const deriveFeatureInstanceHash = async (
  recordKey,
  biologicalFeatureId,
  { crypto = globalThis.crypto } = {}
) => {
  const bytes = featureInstanceIdentityBytes(recordKey, biologicalFeatureId);
  if (!bytes) return null;
  if (!crypto?.subtle?.digest) {
    throw new Error('Feature instance identity requires Web Crypto SHA-256.');
  }
  const digest = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes));
  return `${INSTANCE_HASH_PREFIX}${base32WithoutPadding(digest.slice(0, 16))}`;
};

export const assessFeatureInstanceHashCapability = (catalog) => {
  if (catalog?.schema !== 3 || !Array.isArray(catalog?.items)) {
    return {
      exactScopeAvailable: false,
      canUpgradeFromCanonicalPairs: false,
      featureCount: 0,
      missingFeatureKeys: [],
      invalidFeatureKeys: ['(invalid feature catalog)'],
      diagnostic: 'Generate to enable exact feature edits'
    };
  }
  const featureEntries = Array.isArray(catalog?.items)
    ? catalog.items.flatMap((item) => {
      const recordKeys = new Set(
        Array.isArray(item?.recordKeys) ? item.recordKeys.map(text).filter(Boolean) : []
      );
      return Array.isArray(item?.biologicalFeatures)
        ? item.biologicalFeatures.map((feature) => ({ feature, recordKeys }))
        : [];
    })
    : [];
  const features = featureEntries.map(({ feature }) => feature);
  const missingFeatureKeys = [];
  const invalidFeatureKeys = [];
  const seenCanonicalPairs = new Set();
  const seenHashes = new Set();
  let ambiguous = false;

  featureEntries.forEach(({ feature, recordKeys }) => {
    const recordKey = text(feature?.recordKey);
    const featureId = text(feature?.biologicalFeatureId);
    const canonicalKey = recordKey && featureId
      ? JSON.stringify([recordKey, featureId])
      : '';
    const diagnosticKey = canonicalKey || '(missing canonical feature identity)';
    if (
      !canonicalKey
      || !recordKeys.has(recordKey)
      || seenCanonicalPairs.has(canonicalKey)
    ) {
      ambiguous = true;
      invalidFeatureKeys.push(diagnosticKey);
      return;
    }
    seenCanonicalPairs.add(canonicalKey);

    const instanceHash = text(feature?.instanceHash);
    if (!instanceHash) {
      missingFeatureKeys.push(diagnosticKey);
    } else if (!isFeatureInstanceHash(instanceHash) || seenHashes.has(instanceHash)) {
      ambiguous = true;
      invalidFeatureKeys.push(diagnosticKey);
    } else {
      seenHashes.add(instanceHash);
    }
  });

  return {
    exactScopeAvailable: features.length > 0
      && !ambiguous
      && missingFeatureKeys.length === 0
      && invalidFeatureKeys.length === 0,
    canUpgradeFromCanonicalPairs: features.length > 0
      && !ambiguous
      && invalidFeatureKeys.length === 0,
    featureCount: features.length,
    missingFeatureKeys,
    invalidFeatureKeys,
    diagnostic: features.length === 0
      ? 'Generate to enable exact feature edits'
      : (
        ambiguous || invalidFeatureKeys.length > 0
          ? 'Generate to enable exact feature edits'
          : ''
      )
  };
};

export const upgradeFeatureCatalogInstanceHashes = async (catalog, options = {}) => {
  const capability = assessFeatureInstanceHashCapability(catalog);
  const upgradedCatalog = cloneJson(catalog);
  if (upgradedCatalog?.schema !== 3 || !Array.isArray(upgradedCatalog?.items)) {
    return { catalog: upgradedCatalog, ...capability, upgraded: false };
  }

  let addedInstanceHash = false;
  for (const item of upgradedCatalog.items) {
    for (const feature of Array.isArray(item?.biologicalFeatures)
      ? item.biologicalFeatures
      : []) {
      const instanceHash = text(feature?.instanceHash);
      if (instanceHash && !isFeatureInstanceHash(instanceHash)) continue;
      if (!instanceHash && !capability.canUpgradeFromCanonicalPairs) continue;
      const expectedHash = await deriveFeatureInstanceHash(
        feature.recordKey,
        feature.biologicalFeatureId,
        options
      );
      if (instanceHash) {
        if (expectedHash && instanceHash !== expectedHash) {
          throw new Error(
            'Feature instance hash does not match its canonical record and Feature identity.'
          );
        }
        continue;
      }
      if (capability.canUpgradeFromCanonicalPairs && expectedHash) {
        feature.instanceHash = expectedHash;
        addedInstanceHash = true;
      }
    }
  }
  const upgradedCapability = assessFeatureInstanceHashCapability(upgradedCatalog);
  return {
    catalog: upgradedCatalog,
    ...upgradedCapability,
    upgraded: addedInstanceHash && upgradedCapability.exactScopeAvailable
  };
};
import { cloneJsonData as cloneJson } from './json-clone.js';
