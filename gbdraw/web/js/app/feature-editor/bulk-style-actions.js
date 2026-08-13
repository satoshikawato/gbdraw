import { catalogResultKey } from '../../services/feature-catalog.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { parseCommittedSvgResultRoot } from '../../services/svg-result-ingestion.js';
import { styleFingerprint } from '../../services/style-revision.js';
import {
  prepareFeatureFillLegendProjections,
  preserveResultLocalNonFeatureLegendEntries
} from '../legend/feature-fill-projection.js';
import {
  buildFeatureBulkStyleCommand,
  normalizeFeatureBulkStyleSnapshot,
  replaceFeatureBulkStyleSnapshot
} from './bulk-style-command.js';
import { prepareFeatureFillResultProjection } from './fill-result-projection.js';
import { extractResultLocalLegendEntries } from './fill-actions.js';
import {
  buildNonFeaturePaletteProjectionContext,
  changedNonFeaturePaletteColorKeys,
  prepareNonFeaturePaletteProjection
} from '../svg-styles.js';

const text = (value) => String(value ?? '').trim();

const resultEntries = (entriesByResult, resultKey) => (
  entriesByResult instanceof Map
    ? entriesByResult.get(resultKey) || []
    : entriesByResult?.[resultKey] || []
);

const projectionForResult = (side, resultKey) => (
  (Array.isArray(side?.legendProjections) ? side.legendProjections : [])
    .find((projection) => text(projection?.resultKey) === resultKey) || null
);

const replaceSvgRootContents = (target, source) => {
  if (!target || !source) throw new Error('The mounted bulk style SVG projection is unavailable.');
  for (const attribute of [...target.attributes]) target.removeAttribute(attribute.name);
  for (const attribute of [...source.attributes]) target.setAttribute(attribute.name, attribute.value);
  target.replaceChildren(...[...source.childNodes].map((node) => node.cloneNode(true)));
};

const currentStyle = (state) => normalizeFeatureBulkStyleSnapshot({
  rules: state.manualSpecificRules,
  appliedPaletteName: state.appliedPaletteName?.value,
  appliedPaletteColors: state.appliedPaletteColors?.value
});

/**
 * Own the browser adapters for complete rule/palette snapshot transactions.
 * Presentation callers receive only requestFeatureBulkStyleChange().
 */
export const createFeatureBulkStyleActions = ({
  state,
  history,
  previewRuntime,
  featureSvgActions,
  legendActions,
  prepareGeometryProjection = null,
  onError = null
} = {}) => {
  if (!state || !history) throw new Error('Bulk Feature style actions require state and History.');
  let requestSequence = 0;

  const currentCatalog = () => state.featureCatalog?.value || { schema: 3, items: [] };

  const nonFeatureProjectionContext = () => {
    const mode = text(state.mode?.value).toLowerCase();
    const trackSlots = mode === 'linear'
      ? state.adv?.linear_track_slots
      : state.adv?.circular_track_slots;
    return buildNonFeaturePaletteProjectionContext({
      trackSlots: Array.isArray(trackSlots) ? trackSlots : [],
      defaultDinucleotide: 'GC'
    });
  };

  const mountedContext = () => {
    const catalog = currentCatalog();
    const resultIndex = Number(state.selectedResultIndex?.value ?? 0);
    if (
      !Number.isInteger(resultIndex)
      || resultIndex < 0
      || resultIndex >= (catalog.items?.length || 0)
    ) {
      return null;
    }
    return {
      resultIndex,
      resultKey: catalogResultKey(catalog.items?.[resultIndex]),
      svg: state.svgContainer?.value?.querySelector?.('svg') || null
    };
  };

  const captureRuntimeState = () => {
    const mounted = mountedContext();
    const runtime = previewRuntime?.getActiveRuntime?.() || null;
    return {
      resultIndex: mounted?.resultIndex ?? null,
      resultKey: mounted?.resultKey || '',
      mountedRoot: mounted?.svg || null,
      svg: mounted?.svg?.cloneNode?.(true) || null,
      dirty: Boolean(runtime?.dirty),
      dirtyReasons: [...(runtime?.dirtyReasons || [])]
    };
  };

  const resultLocalLegendPresentation = () => {
    const catalog = currentCatalog();
    const results = state.results.value;
    const mounted = mountedContext();
    const entriesByResult = {};
    (Array.isArray(catalog?.items) ? catalog.items : []).forEach((item, resultIndex) => {
      const resultKey = catalogResultKey(item);
      if (!resultKey || !results?.[resultIndex]) {
        throw new Error(`Bulk Feature style legend Result is unavailable at index ${resultIndex}.`);
      }
      const source = (
        mounted
        && mounted.resultIndex === resultIndex
        && mounted.resultKey === resultKey
        && mounted.svg
      )
        ? mounted.svg
        : parseCommittedSvgResultRoot(results[resultIndex]);
      entriesByResult[resultKey] = extractResultLocalLegendEntries(source);
    });
    return entriesByResult;
  };

  const nonFeatureAffectedResults = ({ before, after }) => {
    if (changedNonFeaturePaletteColorKeys(
      before?.appliedPaletteColors,
      after?.appliedPaletteColors
    ).length === 0) return [];
    const catalog = currentCatalog();
    const results = state.results.value;
    const mounted = mountedContext();
    const context = nonFeatureProjectionContext();
    return (Array.isArray(catalog?.items) ? catalog.items : []).flatMap((item, resultIndex) => {
      const resultKey = catalogResultKey(item);
      if (!resultKey || !results?.[resultIndex]) {
        throw new Error(`Non-Feature palette Result is unavailable at index ${resultIndex}.`);
      }
      const source = (
        mounted
        && mounted.resultIndex === resultIndex
        && mounted.resultKey === resultKey
        && mounted.svg
      )
        ? mounted.svg
        : parseCommittedSvgResultRoot(results[resultIndex]);
      const prepared = prepareNonFeaturePaletteProjection({
        svg: source,
        beforeColors: before.appliedPaletteColors,
        afterColors: after.appliedPaletteColors,
        context,
        strict: true
      });
      return prepared.changed ? [resultKey] : [];
    });
  };

  const restoreRuntimeState = (snapshot) => {
    if (!snapshot) return true;
    const mounted = mountedContext();
    if (!mounted) return snapshot.mountedRoot === null;
    const sameOwner = (
      mounted.resultIndex === snapshot.resultIndex
      && mounted.resultKey === snapshot.resultKey
      && mounted.svg === snapshot.mountedRoot
    );
    const restored = sameOwner
      ? snapshot.svg
      : parseCommittedSvgResultRoot(state.results.value?.[mounted.resultIndex]);
    if (!restored || !mounted.svg) return mounted.svg === null;
    replaceSvgRootContents(mounted.svg, restored);
    previewRuntime?.clearActiveRuntime?.();
    const runtime = previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
    previewRuntime?.invalidatePreviewIndexes?.('feature-bulk-style-rollback');
    if (runtime) {
      runtime.dirty = sameOwner ? snapshot.dirty : false;
      runtime.dirtyReasons = new Set(sameOwner ? snapshot.dirtyReasons : []);
    }
    featureSvgActions?.attachSvgFeatureHandlers?.();
    legendActions?.setupLegendDrag?.();
    return true;
  };

  const prepareLegendProjection = async ({
    from,
    to,
    catalog,
    results,
    affectedResultKeys,
    projections,
    mounted,
    existingEntriesByResult,
    change
  }) => {
    const sourcesByResult = new Map();
    for (const resultKey of affectedResultKeys) {
      const resultIndex = catalog.items.findIndex((item) => catalogResultKey(item) === resultKey);
      if (resultIndex < 0 || !results[resultIndex]) {
        throw new Error(`Bulk Feature style legend Result is unavailable: ${resultKey}`);
      }
      const source = resultIndex === mounted.resultIndex && mounted.svg
        ? mounted.svg
        : parseCommittedSvgResultRoot(results[resultIndex]);
      sourcesByResult.set(resultKey, source);
    }
    const featureAffected = new Set(change?.featureAffectedResultKeys || affectedResultKeys);
    const nonFeatureAffected = new Set(change?.nonFeatureAffectedResultKeys || []);
    const featureProjections = projections
      .filter((projection) => featureAffected.has(projection.resultKey))
      .map((projection) => preserveResultLocalNonFeatureLegendEntries({
        afterProjection: projection,
        beforeProjection: projectionForResult(from, projection.resultKey),
        existingEntries: resultEntries(existingEntriesByResult, projection.resultKey)
      }));
    const featurePrepared = featureProjections.length > 0
      ? await prepareFeatureFillLegendProjections({
          sourcesByResult,
          projections: featureProjections,
          allowAbsentLegend: true
        })
      : { staged: new Map(), counters: {} };
    const staged = new Map(featurePrepared.staged);
    affectedResultKeys.forEach((resultKey) => {
      if (staged.has(resultKey)) return;
      const source = sourcesByResult.get(resultKey);
      staged.set(resultKey, {
        svg: source.cloneNode(true),
        entries: cloneJson(resultEntries(existingEntriesByResult, resultKey)),
        counters: {}
      });
    });
    const nonFeatureCounters = {
      gcPaths: 0,
      skewPaths: 0,
      comparisonPaths: 0,
      legendSwatches: 0,
      gradientStops: 0
    };
    const context = nonFeatureProjectionContext();
    nonFeatureAffected.forEach((resultKey) => {
      const current = staged.get(resultKey);
      const projected = prepareNonFeaturePaletteProjection({
        svg: current?.svg,
        beforeColors: from.appliedPaletteColors,
        afterColors: to.appliedPaletteColors,
        context,
        strict: true
      });
      if (!projected.changed) {
        throw new Error(`Non-Feature palette projection became stale for Result ${resultKey}.`);
      }
      Object.keys(nonFeatureCounters).forEach((key) => {
        nonFeatureCounters[key] += Number(projected.counters?.[key] || 0);
      });
      staged.set(resultKey, {
        ...current,
        svg: projected.svg,
        entries: extractResultLocalLegendEntries(projected.svg)
      });
    });
    return {
      staged,
      selectedEntries: cloneJson(
        staged.get(mounted.resultKey)?.entries
        || state.legendEntries?.value
        || []
      ),
      counters: {
        ...(featurePrepared.counters || {}),
        nonFeature: nonFeatureCounters
      },
      retainedForHistory: {
        entriesByResult: Object.fromEntries([...staged].map(([resultKey, entry]) => (
          [resultKey, entry.entries]
        )))
      }
    };
  };

  const prepareResultProjection = ({
    results,
    catalog,
    fillsByTargetKey,
    affectedResultKeys,
    mounted,
    preparedSvgByResultKey,
    targetFeatureKeysByResult,
    legend
  }) => {
    const projection = prepareFeatureFillResultProjection({
      results,
      catalog,
      fillsByTargetKey,
      affectedResultKeys,
      mounted,
      preparedSvgByResultKey,
      targetFeatureKeysByResult
    });
    return {
      ...projection,
      mountedSvg: projection.preparedMountedSvg,
      selectedLegendEntries: legend?.selectedEntries || [],
      retainedForHistory: {
        counters: projection.counters,
        legend: legend?.retainedForHistory || null
      }
    };
  };

  const commitMountedProjection = ({ prepared }) => {
    const mounted = mountedContext();
    const staged = prepared?.preparedMountedSvg || prepared?.mountedSvg || null;
    if (!staged) return true;
    if (!mounted.svg) return false;
    replaceSvgRootContents(mounted.svg, staged);
    previewRuntime?.clearActiveRuntime?.();
    return true;
  };

  const reconcile = ({ prepared }) => {
    const mounted = mountedContext();
    if (prepared?.projection?.preparedMountedSvg || prepared?.projection?.mountedSvg) {
      previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
      previewRuntime?.invalidatePreviewIndexes?.('feature-bulk-style-commit');
      featureSvgActions?.attachSvgFeatureHandlers?.();
      legendActions?.setupLegendDrag?.();
    }
    return true;
  };

  const requestFeatureBulkStyleChange = async ({
    after = null,
    replacement = null,
    presentationPatch = null,
    writerKind = 'bulk-style',
    label = 'Change feature styles'
  } = {}) => {
    const sequence = ++requestSequence;
    const before = currentStyle(state);
    if (
      !text(state.semanticStyleFingerprint?.value)
      && Number(state.semanticStyleRevision?.value || 0) === 0
      && (state.results?.value?.length || 0) === 0
    ) {
      state.semanticStyleFingerprint.value = styleFingerprint(before);
      state.validatedStyleFingerprintByResultKey.value = Object.freeze({});
    }
    const requestedAfter = after
      ? normalizeFeatureBulkStyleSnapshot(after)
      : replaceFeatureBulkStyleSnapshot(before, replacement || {});
    try {
      const committed = await history.runUndoableCommand(label, async () => {
        const existingEntriesByResult = resultLocalLegendPresentation();
        const nonFeatureAffectedResultKeys = nonFeatureAffectedResults({
          before,
          after: requestedAfter
        });
        const command = await buildFeatureBulkStyleCommand({
          state,
          catalog: currentCatalog(),
          before,
          after: requestedAfter,
          writerKind,
          label,
          prepareLegendProjection,
          prepareResultProjection,
          prepareGeometryProjection: typeof prepareGeometryProjection === 'function'
            ? (context) => prepareGeometryProjection({ ...context, state })
            : null,
          getMountedContext: mountedContext,
          commitMountedProjection,
          restoreMountedProjection: () => {
            previewRuntime?.clearActiveRuntime?.();
            return true;
          },
          captureRuntimeState,
          restoreRuntimeState,
          reconcile,
          refreshPresentation: () => true,
          existingEntriesByResult,
          nonFeatureAffectedResultKeys,
          presentationPatch,
          selectedFeatureTypes: Array.from(state.adv?.features || [])
        });
        return sequence === requestSequence ? command : { noop: true };
      });
      return sequence === requestSequence && committed === true;
    } catch (error) {
      if (sequence !== requestSequence) return false;
      if (typeof onError === 'function') onError(error, { writerKind, label });
      else console.error('Bulk Feature style command failed.', error);
      return false;
    }
  };

  return Object.freeze({ requestFeatureBulkStyleChange });
};
