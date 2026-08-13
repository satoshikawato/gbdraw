import {
  planFeatureFillChange,
  resolveFeatureFillPlan,
  stableFeatureTargetKey
} from './fill-scope-plan.js';
import { buildFeatureFillCommand } from './fill-command.js';
import { prepareFeatureFillResultProjection } from './fill-result-projection.js';
import {
  deriveUsedFeatureFillGroupsByResult,
  prepareFeatureFillLegendProjections,
  preserveResultLocalNonFeatureLegendEntries
} from '../legend/feature-fill-projection.js';
import {
  catalogResultKey
} from '../../services/feature-catalog.js';
import { parseCommittedSvgResultRoot } from '../../services/svg-result-ingestion.js';
import { captureStyleSnapshot } from '../../services/style-revision.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { getAllFeatureLegendGroups, parseTransformXY } from '../legend/utils.js';

const text = (value) => String(value ?? '').trim();

const replaceSvgRootContents = (target, source) => {
  if (!target || !source) throw new Error('The mounted SVG projection is unavailable.');
  for (const attribute of [...target.attributes]) target.removeAttribute(attribute.name);
  for (const attribute of [...source.attributes]) {
    target.setAttribute(attribute.name, attribute.value);
  }
  target.replaceChildren(...[...source.childNodes].map((node) => node.cloneNode(true)));
};

const selectedFeatureKeys = (features) => (Array.isArray(features) ? features : [])
  .map(stableFeatureTargetKey)
  .filter(Boolean);

const directLegendEntryGroups = (targetGroup) => {
  const direct = Array.from(targetGroup?.children || []).filter((child) => (
    text(child?.localName ?? child?.tagName).toLowerCase() === 'g'
    && child.hasAttribute?.('data-legend-key')
  ));
  return direct.length > 0
    ? direct
    : Array.from(targetGroup?.querySelectorAll?.('g[data-legend-key]') || []);
};

const elementIsDisplayed = (element, root) => {
  let current = element;
  while (current) {
    if (
      text(current?.getAttribute?.('display')).toLowerCase() === 'none'
      || text(current?.style?.display).toLowerCase() === 'none'
    ) return false;
    if (current === root) break;
    current = current.parentElement || null;
  }
  return true;
};

const legendTargetGroups = (svg) => {
  const groups = typeof svg?.getElementById === 'function'
    ? getAllFeatureLegendGroups(svg)
    : [];
  if (groups.length > 0) return groups;
  const entry = svg?.querySelector?.('g[data-legend-key]') || null;
  return entry?.parentElement ? [entry.parentElement] : [];
};

const legendEntrySwatch = (group) => Array.from(group?.querySelectorAll?.('path') || [])
  .find((path) => {
    const fill = text(path?.getAttribute?.('fill')).toLowerCase();
    return fill && fill !== 'none' && !fill.startsWith('url(');
  }) || null;

/** @internal Read Result-local presentation hints without promoting them to semantic state. */
export const extractResultLocalLegendEntries = (svg) => {
  const groups = legendTargetGroups(svg);
  const targetGroup = groups.find((group) => elementIsDisplayed(group, svg)) || groups[0] || null;
  if (!targetGroup) return [];
  return directLegendEntryGroups(targetGroup).map((group) => {
    const caption = text(group?.getAttribute?.('data-legend-key'));
    const swatch = legendEntrySwatch(group);
    const label = group?.querySelector?.('text') || null;
    const groupPosition = parseTransformXY(group?.getAttribute?.('transform'));
    const contentPosition = parseTransformXY(
      (label || swatch)?.getAttribute?.('transform')
    );
    const strokeColor = text(swatch?.getAttribute?.('stroke'));
    const rawStrokeWidth = swatch?.getAttribute?.('stroke-width');
    const strokeWidth = rawStrokeWidth === null || rawStrokeWidth === undefined || rawStrokeWidth === ''
      ? null
      : Number(rawStrokeWidth);
    return {
      caption,
      originalCaption: caption,
      color: text(swatch?.getAttribute?.('fill')) || '#cccccc',
      owner: text(group?.getAttribute?.('data-legend-owner')) || 'feature',
      xPos: groupPosition.x + contentPosition.x,
      yPos: groupPosition.y + contentPosition.y,
      showStroke: Boolean(
        strokeColor
        && strokeColor.toLowerCase() !== 'none'
        && (!Number.isFinite(strokeWidth) || strokeWidth > 0)
      ),
      ...(strokeColor ? { strokeColor } : {}),
      ...(Number.isFinite(strokeWidth) ? { strokeWidth } : {}),
      featureIds: []
    };
  }).filter((entry) => entry.caption);
};

const resultBinding = (catalog, index) => {
  const item = catalog?.items?.[index];
  return item ? { item, resultKey: catalogResultKey(item), resultIndex: index } : null;
};

const planDiagnostic = (plan, fallback) => text(plan?.diagnostics?.[0]?.message) || fallback;

export const createFeatureFillActions = ({
  state,
  history,
  previewRuntime,
  featureSvgActions,
  legendActions
} = {}) => {
  if (!state || !history) throw new Error('Feature fill actions require state and History.');
  let requestSequence = 0;

  const clearPending = ({ status = 'idle', progress = null } = {}) => {
    state.pendingFeatureFillPlan.value = null;
    state.featureFillPlanStatus.value = status;
    state.featureFillPlanProgress.value = progress;
  };

  const currentCatalog = () => state.featureCatalog?.value || null;

  const plannerSnapshot = () => {
    const style = captureStyleSnapshot(state);
    return {
      catalog: currentCatalog(),
      features: state.biologicalFeatures?.value,
      extractedFeatures: state.extractedFeatures?.value,
      editableLabels: state.editableLabels?.value,
      specificRules: state.manualSpecificRules,
      paletteColors: state.appliedPaletteColors?.value,
      selectedFeatureKeys: selectedFeatureKeys(state.selectedFeatures?.value),
      exactScopeAvailable: state.featureExactScopeAvailable?.value === true,
      resultGenerationKey: String(style.resultGenerationKey),
      documentEpoch: String(style.documentEpoch),
      styleFingerprint: style.fingerprint
    };
  };

  const mountedContext = () => {
    const resultIndex = Number(state.selectedResultIndex?.value || 0);
    const binding = resultBinding(currentCatalog(), resultIndex);
    return {
      resultIndex,
      resultKey: binding?.resultKey || '',
      svg: state.svgContainer?.value?.querySelector?.('svg') || null
    };
  };

  const captureRuntimeState = () => {
    const mounted = mountedContext();
    const runtime = previewRuntime?.getActiveRuntime?.() || null;
    return {
      resultIndex: mounted.resultIndex,
      resultKey: mounted.resultKey,
      mountedRoot: mounted.svg,
      svg: mounted.svg?.cloneNode?.(true) || null,
      dirty: Boolean(runtime?.dirty),
      dirtyReasons: [...(runtime?.dirtyReasons || [])]
    };
  };

  const restoreRuntimeState = (snapshot) => {
    if (!snapshot) return true;
    const mounted = mountedContext();
    const sameMountedOwner = (
      mounted.resultIndex === snapshot.resultIndex
      && mounted.resultKey === snapshot.resultKey
      && mounted.svg === snapshot.mountedRoot
    );
    const restoredSvg = sameMountedOwner
      ? snapshot.svg
      : parseCommittedSvgResultRoot(state.results.value?.[mounted.resultIndex]);
    if (!restoredSvg || !mounted.svg) return false;
    replaceSvgRootContents(mounted.svg, restoredSvg);
    previewRuntime?.clearActiveRuntime?.();
    const runtime = previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
    previewRuntime?.invalidatePreviewIndexes?.('feature-fill-rollback');
    if (runtime) {
      runtime.dirty = sameMountedOwner ? snapshot.dirty : false;
      runtime.dirtyReasons = new Set(sameMountedOwner ? snapshot.dirtyReasons || [] : []);
    }
    featureSvgActions?.attachSvgFeatureHandlers?.();
    legendActions?.setupLegendDrag?.();
    return true;
  };

  const prepareLegendProjection = async ({
    catalog,
    rules,
    paletteColors,
    affectedResultKeys,
    mounted
  }) => {
    const results = state.results.value;
    const sourcesByResult = new Map();
    const existingEntriesByResult = {};
    for (const resultKey of affectedResultKeys) {
      const resultIndex = catalog.items.findIndex((item) => catalogResultKey(item) === resultKey);
      if (resultIndex < 0 || !results[resultIndex]) {
        throw new Error(`Feature fill legend Result is unavailable: ${resultKey}`);
      }
      const source = resultIndex === mounted.resultIndex && mounted.svg
        ? mounted.svg
        : parseCommittedSvgResultRoot(results[resultIndex]);
      sourcesByResult.set(resultKey, source);
      existingEntriesByResult[resultKey] = extractResultLocalLegendEntries(source);
    }
    const manualLegendEntries = cloneJson(state.manualLegendEntries?.value || []);
    const beforeDerived = deriveUsedFeatureFillGroupsByResult({
      catalog,
      rules: state.manualSpecificRules,
      paletteColors: state.appliedPaletteColors?.value || {},
      manualLegendEntries,
      existingEntriesByResult
    });
    const derived = deriveUsedFeatureFillGroupsByResult({
      catalog,
      rules,
      paletteColors,
      manualLegendEntries,
      existingEntriesByResult
    });
    const projections = derived.projections
      .filter((projection) => affectedResultKeys.includes(projection.resultKey))
      .map((projection) => preserveResultLocalNonFeatureLegendEntries({
        beforeProjection: beforeDerived.projections.find(
          (before) => before.resultKey === projection.resultKey
        ),
        afterProjection: projection,
        existingEntries: existingEntriesByResult[projection.resultKey]
      }));
    const prepared = await prepareFeatureFillLegendProjections({
      sourcesByResult,
      projections,
      allowAbsentLegend: true
    });
    return {
      ...prepared,
      selectedEntries: cloneJson(
        prepared.staged.get(mounted.resultKey)?.entries
        || state.legendEntries?.value
        || []
      ),
      retainedForHistory: {
        entriesByResult: Object.fromEntries(
          projections.map((projection) => [projection.resultKey, projection.entries])
        )
      }
    };
  };

  const prepareResultProjection = ({
    results,
    catalog,
    fillsByTargetKey,
    affectedResultKeys,
    mounted,
    legend,
    preparedSvgByResultKey,
    targetFeatureKeysByResult
  }) => {
    const projection = prepareFeatureFillResultProjection({
      results,
      catalog,
      fillsByTargetKey,
      affectedResultKeys,
      mountedResultIndex: mounted.resultIndex,
      mountedSvg: mounted.svg,
      preparedSvgByResultKey,
      targetFeatureKeysByResult
    });
    return {
      ...projection,
      mountedSvg: projection.preparedMountedSvg
        || legend?.staged?.get?.(mounted.resultKey)?.svg
        || null,
      selectedLegendEntries: legend?.selectedEntries || [],
      retainedForHistory: {
        counters: projection.counters,
        legend: legend?.retainedForHistory || null
      }
    };
  };

  const commitMountedProjection = ({ prepared }) => {
    const mounted = mountedContext();
    if (!prepared?.mountedSvg || !mounted.svg) return true;
    replaceSvgRootContents(mounted.svg, prepared.mountedSvg);
    previewRuntime?.clearActiveRuntime?.();
    return true;
  };

  const reconcileMountedProjection = ({ prepared } = {}) => {
    if (!prepared?.projection?.mountedSvg) return true;
    const mounted = mountedContext();
    previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
    previewRuntime?.invalidatePreviewIndexes?.('feature-fill-commit');
    featureSvgActions?.attachSvgFeatureHandlers?.();
    legendActions?.setupLegendDrag?.();
    return true;
  };

  const executeResolvedPlan = async (plan) => {
    const sequence = ++requestSequence;
    state.featureFillPlanStatus.value = 'preparing';
    state.featureFillPlanProgress.value = {
      completed: 0,
      total: plan.affectedResultKeys.length,
      message: `Preparing ${plan.affectedResultKeys.length} affected Result(s)`
    };
    try {
      const committed = await history.runUndoableCommand('Change feature fill', async () => {
        const command = await buildFeatureFillCommand({
          plan,
          state,
          catalog: currentCatalog(),
          prepareResultProjection,
          prepareLegendProjection,
          getMountedContext: mountedContext,
          commitMountedProjection,
          restoreMountedProjection: () => {
            previewRuntime?.clearActiveRuntime?.();
            return true;
          },
          captureRuntimeState,
          restoreRuntimeState,
          reconcile: reconcileMountedProjection,
          refreshPresentation: () => {
            state.clickedFeature.value = null;
            state.selectedFeatures?.value?.forEach?.(() => {});
            return true;
          }
        });
        if (sequence !== requestSequence || state.pendingFeatureFillPlan.value?.token !== plan.token) {
          return { noop: true };
        }
        return command;
      });
      if (sequence !== requestSequence) return false;
      clearPending({
        status: committed ? 'committed' : 'stale',
        progress: committed
          ? { completed: plan.affectedResultKeys.length, total: plan.affectedResultKeys.length }
          : { completed: 0, total: plan.affectedResultKeys.length, message: 'The edit became stale.' }
      });
      return committed;
    } catch (error) {
      if (sequence !== requestSequence) return false;
      clearPending({
        status: 'failed',
        progress: { completed: 0, total: plan.affectedResultKeys.length, message: error?.message || 'Feature fill failed.' }
      });
      console.error('Feature fill command failed.', error);
      return false;
    }
  };

  const requestFeatureFillChange = async (feat, value, requestedCaption = '', options = {}) => {
    if (!feat) return false;
    const snapshot = plannerSnapshot();
    const targetFeatureKey = stableFeatureTargetKey(feat);
    const binding = resultBinding(currentCatalog(), Number(feat?.resultIndex));
    const token = `fill-request:${snapshot.documentEpoch}:${++requestSequence}`;
    const plan = planFeatureFillChange(snapshot, {
      targetFeatureKey,
      targetFeatureKeys: options.targetFeatureKeys,
      value,
      requestedCaption,
      source: options.source || 'popup',
      originResultIndex: Number.isInteger(Number(feat?.resultIndex))
        ? Number(feat.resultIndex)
        : state.selectedResultIndex.value,
      originResultKey: text(feat?.resultKey) || binding?.resultKey || '',
      resultGenerationKey: snapshot.resultGenerationKey,
      documentEpoch: snapshot.documentEpoch,
      styleFingerprint: snapshot.styleFingerprint,
      semanticScope: options.semanticScope,
      resultExtent: options.resultExtent,
      token
    });
    if (plan.status === 'invalid' || plan.status === 'noop') {
      clearPending({
        status: plan.status,
        progress: plan.status === 'invalid'
          ? { completed: 0, total: 0, message: planDiagnostic(plan, 'Feature fill is unavailable.') }
          : null
      });
      return false;
    }
    const selectedPlan = plan.status === 'conflict' && plan.selectedCandidateId
      ? resolveFeatureFillPlan(plan, plan.selectedCandidateId)
      : plan;
    state.pendingFeatureFillPlan.value = Object.freeze(cloneJson(selectedPlan));
    state.featureFillPlanStatus.value = selectedPlan.status;
    state.featureFillPlanProgress.value = null;
    if (options.closePopupOnDialog && selectedPlan.status !== 'ready') state.clickedFeature.value = null;
    if (selectedPlan.status !== 'ready') return false;
    return executeResolvedPlan(selectedPlan.semanticScope
      ? selectedPlan
      : resolveFeatureFillPlan(selectedPlan));
  };

  const resolveFeatureFillScope = (token, candidateId, conflictChoiceId = null) => {
    const pending = state.pendingFeatureFillPlan.value;
    if (!pending || text(pending.token) !== text(token)) return false;
    const resolved = resolveFeatureFillPlan(pending, candidateId, { conflictChoiceId });
    if (resolved?.status === 'conflict' && resolved?.selectedCandidateId && !conflictChoiceId) {
      state.pendingFeatureFillPlan.value = Object.freeze(cloneJson(resolved));
      return false;
    }
    if (resolved?.status !== 'ready') return false;
    state.pendingFeatureFillPlan.value = Object.freeze(cloneJson(resolved));
    return executeResolvedPlan(resolved);
  };

  const cancelFeatureFillScope = (token = null) => {
    const pending = state.pendingFeatureFillPlan.value;
    if (token && text(pending?.token) !== text(token)) return false;
    requestSequence += 1;
    clearPending();
    return true;
  };

  const requestSelectedFeatureFillChange = (features, value, requestedCaption) => {
    const targets = Array.isArray(features) ? features : [];
    if (targets.length === 0) return false;
    return requestFeatureFillChange(targets[0], value, requestedCaption, {
      source: 'selection-toolbar',
      semanticScope: 'selected-features',
      targetFeatureKeys: selectedFeatureKeys(targets)
    });
  };

  return {
    cancelFeatureFillScope,
    requestFeatureFillChange,
    requestSelectedFeatureFillChange,
    resolveFeatureFillScope
  };
};
