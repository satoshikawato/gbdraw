const normalizedValue = (value, fallback) => {
  const normalized = String(value ?? '').trim().toLowerCase();
  return normalized || fallback;
};

const isPlainObject = (value) => (
  Boolean(value) && typeof value === 'object' && !Array.isArray(value)
);

const hasOwn = (value, key) => (
  isPlainObject(value) && Object.prototype.hasOwnProperty.call(value, key)
);

const migrateOwnField = (source, obsoleteName, currentName) => {
  if (!hasOwn(source, obsoleteName)) return source;
  const migrated = { ...source };
  if (!hasOwn(migrated, currentName)) {
    migrated[currentName] = migrated[obsoleteName];
  }
  delete migrated[obsoleteName];
  return migrated;
};

const requireCurrentValue = (value, fallback, supported, label) => {
  const normalized = normalizedValue(value, fallback);
  if (!supported.includes(normalized)) {
    throw new Error(`${label} must be one of: ${supported.join(', ')}.`);
  }
  return normalized;
};

export const requireCurrentCircularMultiRecordSizeMode = (value) => (
  requireCurrentValue(
    value,
    'auto',
    ['auto', 'linear', 'equal'],
    'Circular multi-record size mode'
  )
);

export const requireCurrentLinearTrackLayout = (value) => (
  requireCurrentValue(
    value,
    'middle',
    ['above', 'middle', 'below'],
    'Linear track layout'
  )
);

export const requireCurrentLinearLabelPlacement = (value) => (
  requireCurrentValue(
    value,
    'auto',
    ['auto', 'above_feature'],
    'Linear label placement'
  )
);

export const migratePersistedCircularMultiRecordSizeMode = (value) => {
  const normalized = normalizedValue(value, 'auto');
  return requireCurrentCircularMultiRecordSizeMode(
    normalized === 'sqrt' ? 'auto' : normalized
  );
};

export const migratePersistedLinearTrackLayout = (value) => {
  const normalized = normalizedValue(value, 'middle');
  if (normalized === 'spreadout') return 'above';
  if (normalized === 'tuckin') return 'below';
  return requireCurrentLinearTrackLayout(normalized);
};

export const migratePersistedLinearLabelPlacement = (value) => {
  const normalized = normalizedValue(value, 'auto');
  return requireCurrentLinearLabelPlacement(
    normalized === 'on_feature' ? 'above_feature' : normalized
  );
};

export const migratePersistedWebStateFieldNames = (config) => {
  if (!isPlainObject(config)) return config;

  let adv = migrateOwnField(
    config.adv,
    'depth_tick_interval',
    'depth_large_tick_interval'
  );
  if (isPlainObject(adv) && Array.isArray(adv.depth_tracks)) {
    adv = {
      ...adv,
      depth_tracks: adv.depth_tracks.map((track) => (
        migrateOwnField(track, 'tick_interval', 'large_tick_interval')
      ))
    };
  }

  let losat = config.losat;
  if (isPlainObject(losat) && isPlainObject(losat.blastp)) {
    const blastp = migrateOwnField(
      losat.blastp,
      'collinearMaxGeneGap',
      'collinearMaxUnitGap'
    );
    if (blastp !== losat.blastp) {
      losat = { ...losat, blastp };
    }
  }

  return {
    ...config,
    ...(adv === undefined ? {} : { adv }),
    ...(losat === undefined ? {} : { losat })
  };
};

export const requireCurrentWebStateFieldNames = (config) => {
  const adv = isPlainObject(config) ? config.adv : null;
  if (hasOwn(adv, 'depth_tick_interval')) {
    throw new Error(
      'Web state field adv.depth_tick_interval is obsolete; use adv.depth_large_tick_interval.'
    );
  }
  if (Array.isArray(adv?.depth_tracks)) {
    adv.depth_tracks.forEach((track, index) => {
      if (hasOwn(track, 'tick_interval')) {
        throw new Error(
          `Web state field adv.depth_tracks[${index}].tick_interval is obsolete; ` +
          'use large_tick_interval.'
        );
      }
    });
  }
  const blastp = isPlainObject(config?.losat) ? config.losat.blastp : null;
  if (hasOwn(blastp, 'collinearMaxGeneGap')) {
    throw new Error(
      'Web state field losat.blastp.collinearMaxGeneGap is obsolete; ' +
      'use losat.blastp.collinearMaxUnitGap.'
    );
  }
  return config;
};
