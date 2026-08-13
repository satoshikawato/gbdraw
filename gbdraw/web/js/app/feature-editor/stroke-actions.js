import { catalogResultKey } from '../../services/feature-catalog.js';
import { getFeatureOverride } from '../../services/feature-override-identity.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { parseCommittedSvgResultRoot } from '../../services/svg-result-ingestion.js';
import { captureStyleSnapshot } from '../../services/style-revision.js';
import { buildFeatureStrokeCommand } from './stroke-command.js';
import { prepareFeatureStrokeResultProjection } from './stroke-result-projection.js';
import {
  planFeatureStrokeChange,
  resolveFeatureStrokePlan,
  stableFeatureTargetKey
} from './stroke-scope-plan.js';
import { resolveFeatureStrokeViewModel } from './stroke-view-model.js';

const text = (value) => String(value ?? '').trim();

const requireRef = (value, label) => {
  if (!value || typeof value !== 'object' || !('value' in value)) {
    throw new Error(`Feature stroke actions require ${label}.`);
  }
  return value;
};

const replaceSvgRootContents = (target, source) => {
  if (!target || !source) throw new Error('The mounted SVG stroke projection is unavailable.');
  for (const attribute of [...target.attributes]) target.removeAttribute(attribute.name);
  for (const attribute of [...source.attributes]) {
    target.setAttribute(attribute.name, attribute.value);
  }
  target.replaceChildren(...[...source.childNodes].map((node) => node.cloneNode(true)));
};

const selectedFeatureKeys = (features) => [...new Set(
  (Array.isArray(features) ? features : []).map(stableFeatureTargetKey).filter(Boolean)
)];

const resultBinding = (catalog, index) => {
  const item = catalog?.items?.[index];
  return item
    ? { item, resultIndex: index, resultKey: catalogResultKey(item) }
    : null;
};

const fallbackLegendCaption = (feature) => text(
  feature?.effectiveLegendCaption
  ?? feature?.effective_legend_caption
  ?? feature?.legendCaption
  ?? feature?.legend_caption
  ?? feature?.caption
  ?? feature?.type
);

const captionOverride = (overrides, caption) => {
  const exactKey = Object.keys(overrides || {}).find(
    (key) => text(key).toLowerCase() === text(caption).toLowerCase()
  );
  return exactKey ? overrides[exactKey] : null;
};

const diagnosticMessage = (plan, fallback) => text(plan?.diagnostics?.[0]?.message) || fallback;

/**
 * Presentation-only Feature stroke facade. Every accepted mutation is built
 * after scope resolution and enters History through buildFeatureStrokeCommand.
 */
export const createFeatureStrokeActions = ({
  state,
  history,
  previewRuntime = null,
  featureSvgActions = null,
  legendActions = null,
  pendingFeatureStrokePlan,
  featureStrokePlanStatus,
  featureStrokePlanProgress,
  getEffectiveLegendCaption = null,
  commandBuilder = buildFeatureStrokeCommand,
  resultProjection = prepareFeatureStrokeResultProjection
} = {}) => {
  if (!state || !history?.runUndoableCommand) {
    throw new Error('Feature stroke actions require state and History.');
  }
  const pendingPlan = requireRef(pendingFeatureStrokePlan, 'a pending-plan ref');
  const planStatus = requireRef(featureStrokePlanStatus, 'a plan-status ref');
  const planProgress = requireRef(featureStrokePlanProgress, 'a plan-progress ref');
  if (typeof commandBuilder !== 'function' || typeof resultProjection !== 'function') {
    throw new Error('Feature stroke actions require command and Result projection owners.');
  }
  let requestSequence = 0;

  const currentCatalog = () => state.featureCatalog?.value || null;

  const effectiveCaption = (featureLike) => {
    const feature = featureLike?.feat || featureLike;
    const resolved = typeof getEffectiveLegendCaption === 'function'
      ? getEffectiveLegendCaption(feature)
      : '';
    return text(resolved) || fallbackLegendCaption(feature);
  };

  const getFeatureStrokeViewModel = (featureLike) => {
    const feature = featureLike?.feat || featureLike;
    if (!feature) return null;
    const caption = effectiveCaption(feature);
    const globalDefault = state.originalSvgStroke?.value || {};
    return resolveFeatureStrokeViewModel({
      exactOverride: getFeatureOverride(state.featureStrokeOverrides, feature) || null,
      captionOverride: captionOverride(state.legendStrokeOverrides, caption),
      legendCaption: caption,
      svgDefaultStroke: {
        color: featureLike?.originalStrokeColor ?? globalDefault.color,
        width: featureLike?.originalStrokeWidth ?? globalDefault.width
      }
    });
  };

  const plannerFeatures = () => (Array.isArray(state.extractedFeatures?.value)
    ? state.extractedFeatures.value
    : []).map((feature) => ({
      ...cloneJson(feature),
      legendCaption: effectiveCaption(feature)
    }));

  const plannerSnapshot = () => {
    const style = captureStyleSnapshot(state);
    return {
      catalog: currentCatalog(),
      features: plannerFeatures(),
      featureStrokeOverrides: cloneJson(state.featureStrokeOverrides || {}),
      legendStrokeOverrides: cloneJson(state.legendStrokeOverrides || {}),
      originalSvgStroke: cloneJson(state.originalSvgStroke?.value || {}),
      selectedFeatureKeys: selectedFeatureKeys(state.selectedFeatures?.value),
      exactScopeAvailable: state.featureExactScopeAvailable?.value === true,
      resultGenerationKey: style.resultGenerationKey,
      documentEpoch: style.documentEpoch,
      styleRevision: style.revision,
      styleFingerprint: style.fingerprint
    };
  };

  const mountedContext = () => {
    const resultIndex = Number(state.selectedResultIndex?.value ?? 0);
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
    if (!mounted.svg && !snapshot.mountedRoot) {
      previewRuntime?.clearActiveRuntime?.();
      return true;
    }
    if (!mounted.svg) return false;
    const sameMountedOwner = (
      mounted.resultIndex === snapshot.resultIndex
      && mounted.resultKey === snapshot.resultKey
      && mounted.svg === snapshot.mountedRoot
    );
    const restoredSvg = sameMountedOwner
      ? snapshot.svg
      : parseCommittedSvgResultRoot(state.results.value?.[mounted.resultIndex]);
    if (!restoredSvg) return false;
    replaceSvgRootContents(mounted.svg, restoredSvg);
    previewRuntime?.clearActiveRuntime?.();
    const runtime = previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
    previewRuntime?.invalidatePreviewIndexes?.('feature-stroke-rollback');
    if (runtime) {
      runtime.dirty = sameMountedOwner ? snapshot.dirty : false;
      runtime.dirtyReasons = new Set(sameMountedOwner ? snapshot.dirtyReasons : []);
    }
    featureSvgActions?.attachSvgFeatureHandlers?.();
    legendActions?.setupLegendDrag?.();
    return true;
  };

  const commitMountedProjection = ({ prepared }) => {
    const mounted = mountedContext();
    if (!mounted.svg || !prepared?.preparedMountedSvg) {
      throw new Error('The affected mounted Feature stroke projection is incomplete.');
    }
    replaceSvgRootContents(mounted.svg, prepared.preparedMountedSvg);
    previewRuntime?.clearActiveRuntime?.();
    return true;
  };

  const restoreMountedProjection = () => {
    previewRuntime?.clearActiveRuntime?.();
    return true;
  };

  const reconcileMountedProjection = ({ prepared }) => {
    if (!prepared?.projection?.preparedMountedSvg) return true;
    const mounted = mountedContext();
    previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
    previewRuntime?.invalidatePreviewIndexes?.('feature-stroke-commit');
    featureSvgActions?.attachSvgFeatureHandlers?.();
    legendActions?.setupLegendDrag?.();
    return true;
  };

  const refreshPresentation = () => {
    const clicked = state.clickedFeature?.value;
    if (!clicked) return true;
    const viewModel = getFeatureStrokeViewModel(clicked);
    if (!viewModel) return true;
    state.clickedFeature.value = {
      ...clicked,
      strokeColor: viewModel.effectiveColor,
      strokeWidth: viewModel.effectiveWidth
    };
    return true;
  };

  const clearPending = ({ status = 'idle', progress = null } = {}) => {
    pendingPlan.value = null;
    planStatus.value = status;
    planProgress.value = progress;
  };

  const executeResolvedPlan = async (plan) => {
    const sequence = ++requestSequence;
    planStatus.value = 'preparing';
    planProgress.value = {
      completed: 0,
      total: plan.affectedResultKeys.length,
      message: `Preparing ${plan.affectedResultKeys.length} affected Result(s)`
    };
    let commandWasNoop = false;
    try {
      const committed = await history.runUndoableCommand('Change feature stroke', async () => {
        const command = await commandBuilder({
          plan,
          state,
          catalog: currentCatalog(),
          prepareResultProjection: resultProjection,
          getMountedContext: mountedContext,
          commitMountedProjection,
          restoreMountedProjection,
          captureRuntimeState,
          restoreRuntimeState,
          reconcile: reconcileMountedProjection,
          refreshPresentation
        });
        commandWasNoop = Boolean(command?.noop);
        if (sequence !== requestSequence || pendingPlan.value?.token !== plan.token) {
          return { noop: true };
        }
        return command;
      });
      if (sequence !== requestSequence) return false;
      const status = committed ? 'committed' : (commandWasNoop ? 'noop' : 'stale');
      clearPending({
        status,
        progress: committed
          ? { completed: plan.affectedResultKeys.length, total: plan.affectedResultKeys.length }
          : (commandWasNoop
              ? null
              : {
                  completed: 0,
                  total: plan.affectedResultKeys.length,
                  message: 'The edit became stale.'
                })
      });
      return committed;
    } catch (error) {
      if (sequence !== requestSequence) return false;
      clearPending({
        status: 'failed',
        progress: {
          completed: 0,
          total: plan.affectedResultKeys.length,
          message: error?.message || 'Feature stroke failed.'
        }
      });
      console.error('Feature stroke command failed.', error);
      return false;
    }
  };

  const requestFeatureStrokeChange = async (featureLike, value, options = {}) => {
    const feature = featureLike?.feat || featureLike;
    if (!feature) return false;
    const snapshot = plannerSnapshot();
    const explicitResultIndex = Number(feature?.resultIndex ?? feature?.result_index);
    const originResultIndex = Number.isInteger(explicitResultIndex)
      ? explicitResultIndex
      : Number(state.selectedResultIndex?.value ?? 0);
    const binding = resultBinding(currentCatalog(), originResultIndex);
    const token = `stroke-request:${snapshot.documentEpoch}:${++requestSequence}`;
    const plan = planFeatureStrokeChange(snapshot, {
      targetFeatureKey: stableFeatureTargetKey(feature),
      targetFeatureKeys: options.targetFeatureKeys,
      value,
      source: options.source || 'popup',
      originResultIndex,
      originResultKey: text(feature?.resultKey ?? feature?.result_key) || binding?.resultKey || '',
      resultGenerationKey: snapshot.resultGenerationKey,
      documentEpoch: snapshot.documentEpoch,
      styleRevision: snapshot.styleRevision,
      styleFingerprint: snapshot.styleFingerprint,
      semanticScope: options.semanticScope,
      token
    });
    if (plan.status === 'invalid') {
      clearPending({
        status: 'invalid',
        progress: {
          completed: 0,
          total: 0,
          message: diagnosticMessage(plan, 'Feature stroke is unavailable.')
        }
      });
      return false;
    }
    pendingPlan.value = Object.freeze(cloneJson(plan));
    planStatus.value = plan.status;
    planProgress.value = null;
    if (plan.status !== 'ready') return false;
    return executeResolvedPlan(resolveFeatureStrokePlan(plan));
  };

  const resolveFeatureStrokeScope = (token, candidateId) => {
    const pending = pendingPlan.value;
    if (!pending || text(pending.token) !== text(token)) return false;
    const resolved = resolveFeatureStrokePlan(pending, candidateId);
    if (resolved?.status !== 'ready') return false;
    pendingPlan.value = Object.freeze(cloneJson(resolved));
    return executeResolvedPlan(resolved);
  };

  const cancelFeatureStrokeScope = (token = null) => {
    const pending = pendingPlan.value;
    if (token && text(pending?.token) !== text(token)) return false;
    requestSequence += 1;
    clearPending();
    return true;
  };

  const requestSelectedFeatureStrokeChange = (features, value) => {
    const targets = Array.isArray(features) ? features.filter(Boolean) : [];
    if (targets.length === 0) return false;
    return requestFeatureStrokeChange(targets[0], value, {
      source: 'selection-toolbar',
      semanticScope: targets.length > 1 ? 'selected-features' : 'single',
      targetFeatureKeys: selectedFeatureKeys(targets)
    });
  };

  return Object.freeze({
    cancelFeatureStrokeScope,
    getFeatureStrokeViewModel,
    requestFeatureStrokeChange,
    requestSelectedFeatureStrokeChange,
    resolveFeatureStrokeScope
  });
};
