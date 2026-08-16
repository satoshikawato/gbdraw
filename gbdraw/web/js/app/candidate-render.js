import { createEditorSvgProjection } from './editor-svg-projection.js';
import { featureStateFromCatalog } from '../services/feature-catalog.js';
import {
  migrateLegacyFeatureOverrides
} from '../services/feature-override-identity.js';
import { cloneJsonValue } from '../services/json-clone.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from '../services/runtime-test-hooks.js';
import { ingestSvgResults } from '../services/svg-result-ingestion.js';

export const prepareReflowResultCommit = ({
  results,
  suppressPairwiseIdentityLegend = false,
  features = [],
  featureColorOverrides = {},
  featureStrokeOverrides = {},
  featureVisibilityOverrides = {},
  featureVisibilityManualRules = [],
  labelTextFeatureOverrides = {},
  labelTextBulkOverrides = {},
  labelTextFeatureOverrideSources = {},
  labelVisibilityOverrides = {},
  legendColorOverrides = {},
  legendStrokeOverrides = {},
  legendEntries = [],
  deletedLegendEntries = [],
  originalLegendOrder = [],
  addedLegendCaptions = [],
  manualSpecificRules = [],
  sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify
}) => {
  const projection = createEditorSvgProjection({
    features,
    featureColorOverrides,
    featureStrokeOverrides,
    featureVisibilityOverrides,
    featureVisibilityManualRules,
    labelTextFeatureOverrides,
    labelTextBulkOverrides,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    legendColorOverrides,
    legendStrokeOverrides,
    legendEntries,
    deletedLegendEntries,
    originalLegendOrder,
    addedLegendCaptions,
    manualSpecificRules,
    suppressPairwiseIdentityLegend
  });
  return {
    results: ingestSvgResults(results, {
      sanitizer,
      transformSvg: (svg) => projection.project(svg).changed
    })
  };
};

export const prepareCandidateRenderCommit = ({
  results,
  catalog,
  mode = '',
  featureColorOverrides,
  featureStrokeOverrides,
  featureVisibilityOverrides = {},
  featureVisibilityManualRules = [],
  labelTextFeatureOverrides = {},
  labelTextBulkOverrides = {},
  labelTextFeatureOverrideSources = {},
  labelVisibilityOverrides = {},
  legendColorOverrides = {},
  legendStrokeOverrides = {},
  legendEntries = [],
  deletedLegendEntries = [],
  originalLegendOrder = [],
  addedLegendCaptions = [],
  manualSpecificRules = [],
  legacyFeatures = [],
  preparedFeatureState = null,
  suppressPairwiseIdentityLegend = false,
  transformSvg = null,
  sanitizer = globalThis.DOMPurify || globalThis.window?.DOMPurify
}) => {
  recordSessionLifecycleEvent('feature-catalog-adoption-start');
  const featureState = preparedFeatureState || featureStateFromCatalog(catalog, { mode });
  recordSessionLifecycleEvent('feature-catalog-adoption-end');
  const fillOverrides = cloneJsonValue(featureColorOverrides, {});
  const strokeOverrides = cloneJsonValue(featureStrokeOverrides, {});
  const migrationOptions = (overrideKind) => ({
    legacyFeatures,
    onDiagnostic: (diagnostic) => {
      const detail = { overrideKind };
      recordStructuralMetric('legacyFeatureOverrideMigrationCallCount', 1, detail);
      recordStructuralMetric(
        'legacyFeatureOverrideCurrentDescriptorCount',
        diagnostic.currentDescriptorCount,
        detail
      );
      recordStructuralMetric(
        'legacyFeatureOverrideLegacyFeatureCount',
        diagnostic.legacyFeatureCount,
        detail
      );
      recordStructuralMetric(
        'legacyFeatureOverrideLegacyFeaturesVisited',
        diagnostic.legacyFeaturesVisited,
        detail
      );
      recordStructuralMetric(
        'legacyFeatureOverrideKeysNeedingMigration',
        diagnostic.legacyKeysNeedingMigration,
        detail
      );
      recordStructuralMetric(
        'legacyFeatureOverrideFullDescriptorComparisonCount',
        diagnostic.fullDescriptorComparisons,
        detail
      );
      recordStructuralMetric(
        'legacyFeatureOverrideIndexedDescriptorComparisonCount',
        diagnostic.indexedDescriptorComparisons,
        detail
      );
      recordStructuralMetric(
        'legacyFeatureOverrideScanSkipCount',
        diagnostic.skippedLegacyFeatureScan ? 1 : 0,
        detail
      );
    }
  });

  recordSessionLifecycleEvent('legacy-feature-override-migration-start');
  migrateLegacyFeatureOverrides(
    fillOverrides,
    featureState.extractedFeatures,
    migrationOptions('fill')
  );
  migrateLegacyFeatureOverrides(
    strokeOverrides,
    featureState.extractedFeatures,
    migrationOptions('stroke')
  );
  recordSessionLifecycleEvent('legacy-feature-override-migration-end');
  recordSessionLifecycleEvent('rule-derived-fill-override-start');
  const projection = createEditorSvgProjection({
    features: featureState.extractedFeatures,
    featureColorOverrides: fillOverrides,
    featureStrokeOverrides: strokeOverrides,
    featureVisibilityOverrides,
    featureVisibilityManualRules,
    labelTextFeatureOverrides,
    labelTextBulkOverrides,
    labelTextFeatureOverrideSources,
    labelVisibilityOverrides,
    legendColorOverrides,
    legendStrokeOverrides,
    legendEntries,
    deletedLegendEntries,
    originalLegendOrder,
    addedLegendCaptions,
    manualSpecificRules,
    suppressPairwiseIdentityLegend
  });
  recordSessionLifecycleEvent('rule-derived-fill-override-end');

  const itemsByIndex = new Map(
    catalog.items.map((item) => [item.resultIndex, item])
  );
  recordSessionLifecycleEvent('svg-admission-start');
  const processedResults = ingestSvgResults(results, {
    sanitizer,
    transformSvg: (svg, { resultIndex }) => {
      const item = itemsByIndex.get(resultIndex);
      if (!item) {
        throw new Error('The diagram engine returned incomplete feature metadata.');
      }
      const projectionResult = projection.project(svg, {
        item,
        requireFeatureBindings: true
      });
      const callerChanged = typeof transformSvg === 'function'
        ? Boolean(transformSvg(svg, { resultIndex }))
        : false;
      return projectionResult.changed || callerChanged;
    }
  });
  recordSessionLifecycleEvent('svg-admission-end');

  return {
    results: processedResults,
    featureState,
    featureColorOverrides: projection.featureColorOverrides,
    featureStrokeOverrides: projection.featureStrokeOverrides
  };
};
