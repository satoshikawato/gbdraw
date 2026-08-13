import {
  normalizePaletteColors,
  resolveColorToHex
} from '../app/color-utils.js';
import {
  INSTANCE_HASH_QUALIFIER,
  parseColorTable,
  parseSpecificRules,
  serializeSpecificRules
} from '../app/file-imports.js';
import { resolveOrderedSpecificRule } from '../app/feature-editor/fill-view-model.js';
import { FEATURE_SEMANTIC_SCOPE_QUALIFIER } from '../app/feature-editor/semantic-fill-selectors.js';
import {
  CANONICAL_REQUEST_SCHEMA,
  decodeCanonicalResourceText
} from './session-request.js';

const V40_SESSION_VERSION = 40;
const V40_CANONICAL_REQUEST_SCHEMA = 5;
const COLOR_PATTERN = /^(?:none|#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})|[a-z]+)$/i;

const isObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const text = (value) => String(value ?? '').trim();

const comparableColor = (value) => text(resolveColorToHex(text(value))).toLowerCase();

const comparableColorMap = (value) => {
  if (!isObject(value)) return null;
  const normalized = normalizePaletteColors(value);
  const entries = Object.entries(normalized).map(([key, color]) => {
    const normalizedKey = text(key);
    const normalizedColor = comparableColor(color);
    if (!normalizedKey || !normalizedColor || !COLOR_PATTERN.test(normalizedColor)) {
      throw new Error('Legacy color recovery found an invalid applied color value.');
    }
    return [normalizedKey, normalizedColor];
  });
  entries.sort(([left], [right]) => left.localeCompare(right));
  return Object.fromEntries(entries);
};

const sameJson = (left, right) => JSON.stringify(left) === JSON.stringify(right);

const canonicalColorRef = (colors, primaryField, fallbackField, label) => {
  const primary = colors?.[primaryField];
  const fallback = colors?.[fallbackField];
  const ids = [primary, fallback]
    .filter((entry) => isObject(entry) && text(entry.resourceId))
    .map((entry) => text(entry.resourceId));
  if (new Set(ids).size > 1) {
    throw new Error(`Legacy color recovery found ambiguous canonical ${label} resources.`);
  }
  return primary?.resourceId ? primary : (fallback?.resourceId ? fallback : null);
};

const strictColorResourceMap = (resources, resourceId) => {
  const resourceText = decodeCanonicalResourceText(resources, resourceId);
  const seen = new Set();
  for (const [index, rawLine] of resourceText.split(/\r?\n/).entries()) {
    const line = rawLine.trim();
    if (!line || line.startsWith('#') || line.startsWith('[')) continue;
    const cells = rawLine.split('\t').map((cell) => cell.trim());
    if (
      cells.length === 2
      && cells[0].toLowerCase() === 'feature_type'
      && cells[1].toLowerCase() === 'color'
    ) continue;
    if (
      cells.length !== 2
      || !cells[0]
      || !cells[1]
      || !COLOR_PATTERN.test(cells[1])
      || seen.has(cells[0])
    ) {
      throw new Error(
        `Legacy default-color recovery found an invalid TSV row at line ${index + 1}.`
      );
    }
    seen.add(cells[0]);
  }
  const parsed = parseColorTable(resourceText);
  if (parsed.count !== seen.size) {
    throw new Error('Legacy default-color recovery could not parse the complete saved table.');
  }
  return comparableColorMap(parsed.colors);
};

const comparableRule = (rule) => ({
  feat: text(rule?.feat),
  qual: text(rule?.qual),
  val: text(rule?.val),
  color: comparableColor(rule?.color),
  cap: text(rule?.cap),
  ...([INSTANCE_HASH_QUALIFIER, FEATURE_SEMANTIC_SCOPE_QUALIFIER].includes(text(rule?.qual))
    ? { match: rule?.match === 'literal' ? 'literal' : 'regex' }
    : {})
});

const comparableStoredRules = (rules, schema) => {
  if (!Array.isArray(rules)) return null;
  const serialized = serializeSpecificRules(rules, { schema });
  const parsed = parseSpecificRules(serialized, { schema }).rules.map(comparableRule);
  if (parsed.length !== rules.length) {
    throw new Error('Legacy specific-color recovery found incomplete or duplicate editor rules.');
  }
  return parsed;
};

const comparableResourceRules = (resources, resourceId, schema) => (
  parseSpecificRules(
    decodeCanonicalResourceText(resources, resourceId),
    { schema }
  ).rules.map(comparableRule)
);

const validateCaptionColors = (rules, legendEntries) => {
  const captionColors = new Map();
  rules.forEach((rule) => {
    if (!rule.cap) return;
    const previous = captionColors.get(rule.cap);
    if (previous && previous !== rule.color) {
      throw new Error(
        `Legacy color recovery found conflicting colors for caption ${JSON.stringify(rule.cap)}.`
      );
    }
    captionColors.set(rule.cap, rule.color);
  });

  const savedLegendColors = new Map();
  (Array.isArray(legendEntries) ? legendEntries : []).forEach((entry) => {
    const caption = text(entry?.caption);
    if (!caption) return;
    const color = comparableColor(entry?.color);
    const previous = savedLegendColors.get(caption);
    if ((previous && previous !== color) || (captionColors.has(caption) && captionColors.get(caption) !== color)) {
      throw new Error(
        `Legacy color recovery found conflicting saved legend state for ${JSON.stringify(caption)}.`
      );
    }
    savedLegendColors.set(caption, color);
  });
};

const referencedResourceIds = (value, target = new Set()) => {
  if (Array.isArray(value)) {
    value.forEach((entry) => referencedResourceIds(entry, target));
    return target;
  }
  if (!isObject(value)) return target;
  if (text(value.resourceId)) target.add(text(value.resourceId));
  Object.values(value).forEach((entry) => referencedResourceIds(entry, target));
  return target;
};

const matchingUnreferencedResources = ({
  resources,
  renderRequest,
  kindPattern,
  parse,
  expected
}) => {
  const referenced = referencedResourceIds(renderRequest);
  const matches = [];
  Object.entries(resources || {}).forEach(([resourceId, resource]) => {
    if (referenced.has(resourceId) || !isObject(resource)) return;
    const identity = `${text(resource.kind)} ${text(resource.name)}`;
    if (!kindPattern.test(identity)) return;
    try {
      if (sameJson(parse(resourceId), expected)) matches.push(resourceId);
    } catch {
      // An unrelated or malformed resource is not recovery evidence.
    }
  });
  return matches.sort();
};

const biologicalFeatureKey = (feature) => {
  const recordKey = text(feature?.recordKey);
  const featureId = text(feature?.biologicalFeatureId);
  return recordKey && featureId ? `${recordKey}\u0000${featureId}` : '';
};

const escapeRegExp = (value) => value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

const savedSvgFill = (resultContent, svgId) => {
  if (!svgId || /["<>]/.test(svgId)) {
    throw new Error('Legacy color recovery found an invalid rendered Feature ID.');
  }
  const expression = new RegExp(`<[^>]*\\sid="${escapeRegExp(svgId)}"[^>]*>`, 'g');
  const matches = String(resultContent || '').match(expression) || [];
  if (matches.length !== 1) {
    throw new Error(`Legacy color recovery could not resolve saved Feature ${JSON.stringify(svgId)}.`);
  }
  const fill = matches[0].match(/\bfill="([^"]+)"/i)?.[1] || '';
  if (!fill) {
    throw new Error(`Legacy color recovery found a saved Feature without fill: ${JSON.stringify(svgId)}.`);
  }
  return comparableColor(fill);
};

const replaceSavedSvgFill = (resultContent, svgId, color) => {
  const expression = new RegExp(`<[^>]*\\sid="${escapeRegExp(svgId)}"[^>]*>`, 'g');
  const matches = String(resultContent || '').match(expression) || [];
  if (matches.length !== 1 || !/\bfill="[^"]+"/i.test(matches[0])) {
    throw new Error(`Legacy color recovery could not reproject saved Feature ${JSON.stringify(svgId)}.`);
  }
  const replacement = matches[0].replace(/\bfill="[^"]+"/i, `fill="${color}"`);
  return String(resultContent).replace(matches[0], replacement);
};

const featureDefaultColor = (colors, feature) => (
  colors?.[text(feature?.type)] || colors?.default || ''
);

const validateSavedColorEvidence = ({
  session,
  expectedRules,
  canonicalDefaults,
  recoveredDefaults
}) => {
  const catalog = session.editorState?.featureCatalog;
  const results = session.results;
  if (
    !isObject(catalog)
    || catalog.schema !== 3
    || !Array.isArray(catalog.items)
    || !Array.isArray(results)
    || catalog.items.length !== results.length
  ) {
    throw new Error('Legacy color recovery requires a complete schema-3 Feature catalogue.');
  }

  const globalRecordKeys = [];
  catalog.items.forEach((item) => {
    if (!Array.isArray(item?.recordKeys)) {
      throw new Error('Legacy color recovery found incomplete catalogue record identity.');
    }
    item.recordKeys.forEach((recordKeyRaw) => {
      const recordKey = text(recordKeyRaw);
      if (!recordKey) {
        throw new Error('Legacy color recovery found incomplete catalogue record identity.');
      }
      if (!globalRecordKeys.includes(recordKey)) globalRecordKeys.push(recordKey);
    });
  });

  const featureByAlias = new Map();
  const expectedByFeatureKey = new Map();
  const matchedFeatureKeys = new Set();
  const globalBiologicalFeatureKeys = new Set();
  const observedDefaultSources = new Map();
  const fillProjectionByResult = new Map();
  const addAlias = (alias, featureKey) => {
    const normalized = text(alias);
    if (!normalized) return;
    if (!featureByAlias.has(normalized)) featureByAlias.set(normalized, new Set());
    featureByAlias.get(normalized).add(featureKey);
  };

  catalog.items.forEach((item, resultIndex) => {
    if (
      item?.resultIndex !== resultIndex
      || text(item?.resultName) !== text(results[resultIndex]?.name)
      || !Array.isArray(item?.biologicalFeatures)
      || !Array.isArray(item?.features)
    ) {
      throw new Error('Legacy color recovery found ambiguous catalogue Result ownership.');
    }
    const biologicalByKey = new Map();
    item.biologicalFeatures.forEach((feature) => {
      const featureKey = biologicalFeatureKey(feature);
      if (
        !featureKey
        || biologicalByKey.has(featureKey)
        || globalBiologicalFeatureKeys.has(featureKey)
      ) {
        throw new Error('Legacy color recovery found ambiguous biological Feature identity.');
      }
      biologicalByKey.set(featureKey, feature);
      globalBiologicalFeatureKeys.add(featureKey);
    });
    const renderedFeatureKeys = new Set();
    item.features.forEach((rendered) => {
      const featureKey = biologicalFeatureKey(rendered);
      const biological = biologicalByKey.get(featureKey);
      const svgId = text(rendered?.svgId);
      if (!biological || !svgId) {
        throw new Error('Legacy color recovery found an unresolved rendered Feature.');
      }
      renderedFeatureKeys.add(featureKey);
      const resolved = resolveOrderedSpecificRule(biological, expectedRules);
      const expectedRule = resolved?.rule || null;
      const catalogFill = comparableColor(rendered.fillColor);
      const resultFill = savedSvgFill(results[resultIndex].content, svgId);
      if (!catalogFill || catalogFill !== resultFill) {
        throw new Error('Legacy color recovery found conflicting catalogue and Result Feature fills.');
      }
      if (expectedRule) {
        const expectedColor = comparableColor(expectedRule.color);
        if (catalogFill !== expectedColor) {
          throw new Error('Legacy color recovery found a saved Result that conflicts with editor rules.');
        }
        matchedFeatureKeys.add(featureKey);
        expectedByFeatureKey.set(featureKey, {
          color: expectedColor,
          caption: text(expectedRule.cap)
        });
        if (!fillProjectionByResult.has(resultIndex)) fillProjectionByResult.set(resultIndex, new Map());
        fillProjectionByResult.get(resultIndex).set(svgId, expectedColor);
      } else {
        const canonicalColor = comparableColor(featureDefaultColor(canonicalDefaults, biological));
        const recoveredColor = comparableColor(featureDefaultColor(recoveredDefaults, biological));
        const source = catalogFill === canonicalColor && catalogFill === recoveredColor
          ? 'both'
          : (catalogFill === canonicalColor
              ? 'canonical'
              : (catalogFill === recoveredColor ? 'recovered' : ''));
        if (!source) {
          throw new Error('Legacy color recovery found an unexplained saved default Feature fill.');
        }
        const featureType = text(biological.type);
        if (!observedDefaultSources.has(featureType)) observedDefaultSources.set(featureType, new Set());
        observedDefaultSources.get(featureType).add(source);
        if (!fillProjectionByResult.has(resultIndex)) fillProjectionByResult.set(resultIndex, new Map());
        fillProjectionByResult.get(resultIndex).set(svgId, recoveredColor || catalogFill);
      }

      addAlias(featureKey, featureKey);
      addAlias(svgId, featureKey);
      addAlias(biological.stableFeatureId, featureKey);
      const recordIndex = globalRecordKeys.indexOf(text(biological.recordKey));
      const sourceFeatureIndex = Number(biological.sourceFeatureIndex);
      if (recordIndex >= 0 && Number.isSafeInteger(sourceFeatureIndex) && sourceFeatureIndex >= 0) {
        addAlias(`file${recordIndex}_f${sourceFeatureIndex}`, featureKey);
      }
    });
    biologicalByKey.forEach((_feature, featureKey) => {
      if (!renderedFeatureKeys.has(featureKey)) return;
      addAlias(featureKey, featureKey);
    });
  });

  observedDefaultSources.forEach((sources) => {
    if (sources.has('canonical') && sources.has('recovered')) {
      throw new Error('Legacy color recovery found mixed default-color revisions in saved Results.');
    }
  });

  const projectArtifacts = (featureColorOverrides) => {
    const projectedCatalog = typeof structuredClone === 'function'
      ? structuredClone(catalog)
      : JSON.parse(JSON.stringify(catalog));
    const projectedResults = results.map((result) => ({ ...result }));
    projectedCatalog.items.forEach((item, resultIndex) => {
      const fills = fillProjectionByResult.get(resultIndex) || new Map();
      item.features.forEach((rendered) => {
        const svgId = text(rendered.svgId);
        if (fills.has(svgId)) rendered.fillColor = fills.get(svgId);
      });
      fills.forEach((color, svgId) => {
        projectedResults[resultIndex].content = replaceSavedSvgFill(
          projectedResults[resultIndex].content,
          svgId,
          color
        );
      });
    });
    return {
      featureColorOverrides,
      catalog: projectedCatalog,
      results: projectedResults
    };
  };

  const overrides = session.features?.featureColorOverrides;
  if (!isObject(overrides)) {
    if (matchedFeatureKeys.size > 0) {
      throw new Error('Legacy color recovery is missing its derived Feature override evidence.');
    }
    return projectArtifacts({});
  }
  const overrideFeatureKeys = new Set();
  Object.entries(overrides).forEach(([alias, override]) => {
    const targets = featureByAlias.get(alias);
    if (!targets || targets.size !== 1 || !isObject(override)) {
      throw new Error('Legacy color recovery found an ambiguous derived Feature override.');
    }
    const [featureKey] = targets;
    const expected = expectedByFeatureKey.get(featureKey);
    if (
      !expected
      || expected.color !== comparableColor(override.color)
      || expected.caption !== text(override.caption)
      || overrideFeatureKeys.has(featureKey)
    ) {
      throw new Error('Legacy color recovery found a conflicting derived Feature override.');
    }
    overrideFeatureKeys.add(featureKey);
  });
  if (
    overrideFeatureKeys.size !== matchedFeatureKeys.size
    || [...matchedFeatureKeys].some((featureKey) => !overrideFeatureKeys.has(featureKey))
  ) {
    throw new Error('Legacy color recovery found incomplete derived Feature override evidence.');
  }
  const featureColorOverrides = Object.fromEntries(
    [...expectedByFeatureKey.entries()].map(([featureKey, override]) => [
      featureKey,
      { ...override }
    ])
  );
  return projectArtifacts(featureColorOverrides);
};

const recoveredResourceRef = (resourceId, resource) => ({
  resourceId,
  representation: /canonical/i.test(text(resource?.kind)) ? 'canonicalTsv' : 'file'
});

/**
 * Recover the branch-backed v40 editor color authority without mutating the
 * envelope, request, catalogue, Results, or resource entries supplied by the caller.
 */
export const recoverV40ColorAuthority = (session) => {
  if (!isObject(session) || Number(session.version) !== V40_SESSION_VERSION) {
    return {
      session,
      recovered: false,
      specificStatus: 'not-applicable',
      defaultStatus: 'not-applicable'
    };
  }
  if (
    !isObject(session.renderRequest)
    || Number(session.renderRequest.schema) !== V40_CANONICAL_REQUEST_SCHEMA
    || !isObject(session.resources)
  ) {
    throw new Error('Session v40 color recovery requires a canonical schema-5 request.');
  }

  const request = session.renderRequest;
  const colors = request.diagramOptions?.colors;
  if (!isObject(colors)) {
    throw new Error('Session v40 color recovery requires canonical color options.');
  }
  const specificRef = canonicalColorRef(colors, 'colorTableFile', 'colorTable', 'specific-color');
  const defaultRef = canonicalColorRef(colors, 'defaultColorsFile', 'defaultColors', 'default-color');
  if (Array.isArray(session.config?.rules)) {
    const reservedQualifier = session.config.rules.find(
      (rule) => [INSTANCE_HASH_QUALIFIER, FEATURE_SEMANTIC_SCOPE_QUALIFIER]
        .includes(text(rule?.qual))
    )?.qual;
    if (reservedQualifier) {
      const label = reservedQualifier === INSTANCE_HASH_QUALIFIER
        ? 'instance'
        : 'semantic';
      throw new Error(
        `Session v40 cannot safely promote a schema-5 reserved ${label} selector; Generate again.`
      );
    }
  }
  const canonicalRules = specificRef
    ? comparableResourceRules(session.resources, specificRef.resourceId, V40_CANONICAL_REQUEST_SCHEMA)
    : [];
  const storedRules = Array.isArray(session.config?.rules) && session.config.rules.length > 0
    ? comparableStoredRules(session.config.rules, V40_CANONICAL_REQUEST_SCHEMA)
    : null;
  const expectedRules = storedRules || canonicalRules;
  const reservedRule = [...canonicalRules, ...expectedRules].find(
    (rule) => [INSTANCE_HASH_QUALIFIER, FEATURE_SEMANTIC_SCOPE_QUALIFIER]
      .includes(rule.qual)
  );
  if (reservedRule) {
    const label = reservedRule.qual === INSTANCE_HASH_QUALIFIER
      ? 'instance'
      : 'semantic';
    throw new Error(
      `Session v40 cannot safely promote a schema-5 reserved ${label} selector; Generate again.`
    );
  }
  validateCaptionColors(expectedRules, session.editorState?.legend?.entries);

  const canonicalDefaults = defaultRef
    ? strictColorResourceMap(session.resources, defaultRef.resourceId)
    : comparableColorMap({});
  const configColors = comparableColorMap(session.config?.colors);
  const appliedColors = comparableColorMap(session.ui?.appliedPaletteColors);
  const supportedAppliedColors = configColors && appliedColors && sameJson(configColors, appliedColors)
    ? appliedColors
    : null;
  const configPaletteName = text(session.config?.palette);
  const appliedPaletteName = text(session.ui?.appliedPaletteName);
  const canonicalPaletteName = text(colors.defaultColorsPalette) || 'default';
  const supportedPaletteName = (
    configPaletteName
    && appliedPaletteName
    && configPaletteName === appliedPaletteName
  ) ? appliedPaletteName : '';

  const specificMismatch = storedRules !== null && !sameJson(canonicalRules, storedRules);
  const defaultMismatch = supportedAppliedColors !== null && (
    !sameJson(canonicalDefaults, supportedAppliedColors)
    || (supportedPaletteName && canonicalPaletteName !== supportedPaletteName)
  );
  if (
    !supportedAppliedColors
    && configColors
    && appliedColors
    && !sameJson(configColors, appliedColors)
    && (
      !sameJson(canonicalDefaults, configColors)
      || !sameJson(canonicalDefaults, appliedColors)
    )
  ) {
    throw new Error('Session v40 color recovery found conflicting applied color authorities.');
  }

  let recoveredSpecificId = '';
  if (specificMismatch) {
    const matches = matchingUnreferencedResources({
      resources: session.resources,
      renderRequest: request,
      kindPattern: /(?:specific.*color|color.*table)/i,
      parse: (resourceId) => comparableResourceRules(
        session.resources,
        resourceId,
        V40_CANONICAL_REQUEST_SCHEMA
      ),
      expected: storedRules
    });
    if (matches.length !== 1) {
      throw new Error(
        matches.length === 0
          ? 'Session v40 specific-color mismatch has no deterministic saved resource support.'
          : 'Session v40 specific-color mismatch has ambiguous saved resource support.'
      );
    }
    [recoveredSpecificId] = matches;
  }

  let recoveredDefaultId = '';
  if (defaultMismatch) {
    if (!supportedPaletteName) {
      throw new Error('Session v40 default-color recovery found conflicting palette names.');
    }
    const matches = matchingUnreferencedResources({
      resources: session.resources,
      renderRequest: request,
      kindPattern: /(?:default.*color|color.*default)/i,
      parse: (resourceId) => strictColorResourceMap(session.resources, resourceId),
      expected: supportedAppliedColors
    });
    if (matches.length !== 1) {
      throw new Error(
        matches.length === 0
          ? 'Session v40 default-color mismatch has no deterministic saved resource support.'
          : 'Session v40 default-color mismatch has ambiguous saved resource support.'
      );
    }
    [recoveredDefaultId] = matches;
  }

  const recoveredProjection = specificMismatch || defaultMismatch
    ? validateSavedColorEvidence({
      session,
      expectedRules,
      canonicalDefaults,
      recoveredDefaults: supportedAppliedColors || canonicalDefaults
    })
    : null;

  const nextColors = specificMismatch || defaultMismatch
    ? {
        ...colors,
        ...(specificMismatch
          ? {
              colorTable: null,
              colorTableFile: recoveredResourceRef(
                recoveredSpecificId,
                session.resources[recoveredSpecificId]
              )
            }
          : {}),
        ...(defaultMismatch
          ? {
              defaultColors: null,
              defaultColorsFile: recoveredResourceRef(
                recoveredDefaultId,
                session.resources[recoveredDefaultId]
              ),
              defaultColorsPalette: supportedPaletteName
            }
          : {})
      }
    : colors;
  const nextRequest = {
    ...request,
    schema: CANONICAL_REQUEST_SCHEMA,
    ...(nextColors === colors
      ? {}
      : {
          diagramOptions: {
            ...request.diagramOptions,
            colors: nextColors
          }
        })
  };
  return {
    session: {
      ...session,
      renderRequest: nextRequest,
      ...(recoveredProjection === null
        ? {}
        : {
            features: {
              ...(isObject(session.features) ? session.features : {}),
              featureColorOverrides: recoveredProjection.featureColorOverrides
            },
            editorState: {
              ...(isObject(session.editorState) ? session.editorState : {}),
              featureCatalog: recoveredProjection.catalog
            },
            results: recoveredProjection.results
          })
    },
    recovered: specificMismatch || defaultMismatch,
    specificStatus: specificMismatch ? 'recovered' : 'consistent',
    defaultStatus: defaultMismatch ? 'recovered' : 'consistent',
    supportingResourceIds: Object.freeze({
      specific: recoveredSpecificId || null,
      defaults: recoveredDefaultId || null
    })
  };
};
