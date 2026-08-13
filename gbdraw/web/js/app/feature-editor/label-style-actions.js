import { catalogResultKey } from '../../services/feature-catalog.js';
import { cloneJsonData as cloneJson } from '../../services/json-clone.js';
import { parseCommittedSvgResultRoot } from '../../services/svg-result-ingestion.js';
import { captureStyleSnapshot } from '../../services/style-revision.js';
import { buildFeatureLabelCommand } from './label-command.js';
import { prepareFeatureLabelResultProjection } from './label-result-projection.js';
import {
  planFeatureLabelChange,
  resolveFeatureLabelPlan,
  stableFeatureTargetKey
} from './label-scope-plan.js';
import { resolveFeatureLabelViewModel } from './label-view-model.js';

const text = (value) => String(value ?? '').trim();

const requireRef = (value, label) => {
  if (!value || typeof value !== 'object' || !('value' in value)) {
    throw new Error(`Feature label actions require ${label}.`);
  }
  return value;
};

const replaceSvgRootContents = (target, source) => {
  if (!target || !source) throw new Error('The mounted SVG label projection is unavailable.');
  for (const attribute of [...target.attributes]) target.removeAttribute(attribute.name);
  for (const attribute of [...source.attributes]) target.setAttribute(attribute.name, attribute.value);
  target.replaceChildren(...[...source.childNodes].map((node) => node.cloneNode(true)));
};

const resultBinding = (catalog, index) => {
  const item = catalog?.items?.[index];
  return item ? { item, resultIndex: index, resultKey: catalogResultKey(item) } : null;
};

const renderedSvgIds = (featureLike) => {
  const feature = featureLike?.feat || featureLike;
  return [...new Set([
    featureLike?.svg_id,
    featureLike?.svgId,
    featureLike?.renderedFeatureSvgId,
    featureLike?.rendered_feature_svg_id,
    feature?.svg_id,
    feature?.svgId,
    feature?.renderedFeatureSvgId,
    feature?.rendered_feature_svg_id
  ].map(text).filter(Boolean))];
};

const selectedFeatureKeys = (features) => [...new Set(
  (Array.isArray(features) ? features : []).map((entry) => (
    stableFeatureTargetKey(entry?.feat || entry)
  )).filter(Boolean)
)];

const diagnosticMessage = (plan, fallback) => text(plan?.diagnostics?.[0]?.message) || fallback;

/**
 * Presentation facade for label text. It exposes request/resolve/cancel only;
 * all mutation is delegated to the atomic command after scope resolution.
 */
export const createFeatureLabelStyleActions = ({
  state,
  history,
  previewRuntime = null,
  featureSvgActions = null,
  legendActions = null,
  pendingFeatureLabelPlan,
  featureLabelPlanStatus,
  featureLabelPlanProgress,
  prepareLabelReflow = null,
  refreshLabelPresentation = null,
  commandBuilder = buildFeatureLabelCommand,
  resultProjection = prepareFeatureLabelResultProjection
} = {}) => {
  if (!state || !history?.runUndoableCommand) {
    throw new Error('Feature label actions require state and History.');
  }
  const pendingPlan = requireRef(pendingFeatureLabelPlan, 'a pending-plan ref');
  const planStatus = requireRef(featureLabelPlanStatus, 'a plan-status ref');
  const planProgress = requireRef(featureLabelPlanProgress, 'a plan-progress ref');
  if (typeof commandBuilder !== 'function' || typeof resultProjection !== 'function') {
    throw new Error('Feature label actions require command and Result projection owners.');
  }
  let requestSequence = 0;

  const currentCatalog = () => state.featureCatalog?.value || null;

  const presentationLabelsBySvgId = () => Object.fromEntries(
    (Array.isArray(state.editableLabels?.value) ? state.editableLabels.value : [])
      .map((entry) => [text(entry?.featureId), {
        text: String(entry?.text ?? ''),
        sourceText: String(entry?.sourceText ?? ''),
        labelKey: text(entry?.key)
      }])
      .filter(([svgId]) => svgId)
  );

  const getFeatureLabelViewModel = (featureLike) => {
    const feature = featureLike?.feat || featureLike;
    if (!feature) return null;
    const ids = renderedSvgIds(featureLike);
    const presented = ids.map((id) => presentationLabelsBySvgId()[id]).find(Boolean) || null;
    return resolveFeatureLabelViewModel({
      feature,
      renderedSvgIds: ids,
      featureOverrides: state.labelTextFeatureOverrides,
      bulkOverrides: state.labelTextBulkOverrides,
      featureOverrideSources: state.labelTextFeatureOverrideSources,
      presentationText: presented ? presented.text : null,
      presentationSourceText: presented?.sourceText || ''
    });
  };

  const plannerSnapshot = () => {
    const style = captureStyleSnapshot(state);
    return {
      catalog: currentCatalog(),
      biologicalFeatures: state.biologicalFeatures?.value,
      extractedFeatures: state.extractedFeatures?.value,
      editableLabels: state.editableLabels?.value,
      presentationLabelsBySvgId: presentationLabelsBySvgId(),
      labelTextFeatureOverrides: cloneJson(state.labelTextFeatureOverrides || {}),
      labelTextBulkOverrides: cloneJson(state.labelTextBulkOverrides || {}),
      labelTextFeatureOverrideSources: cloneJson(state.labelTextFeatureOverrideSources || {}),
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
    const sameOwner = mounted.resultIndex === snapshot.resultIndex
      && mounted.resultKey === snapshot.resultKey && mounted.svg === snapshot.mountedRoot;
    const restored = sameOwner
      ? snapshot.svg
      : parseCommittedSvgResultRoot(state.results.value?.[mounted.resultIndex]);
    if (!restored) return false;
    replaceSvgRootContents(mounted.svg, restored);
    previewRuntime?.clearActiveRuntime?.();
    const runtime = previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
    previewRuntime?.invalidatePreviewIndexes?.('feature-label-rollback');
    if (runtime) {
      runtime.dirty = sameOwner ? snapshot.dirty : false;
      runtime.dirtyReasons = new Set(sameOwner ? snapshot.dirtyReasons : []);
    }
    featureSvgActions?.attachSvgFeatureHandlers?.();
    legendActions?.setupLegendDrag?.();
    return true;
  };

  const commitMountedProjection = ({ prepared }) => {
    const mounted = mountedContext();
    if (!mounted.svg || !prepared?.preparedMountedSvg) {
      throw new Error('The affected mounted Feature label projection is incomplete.');
    }
    replaceSvgRootContents(mounted.svg, prepared.preparedMountedSvg);
    previewRuntime?.clearActiveRuntime?.();
    return true;
  };

  const reconcileMountedProjection = ({ prepared }) => {
    if (prepared?.projection?.preparedMountedSvg) {
      const mounted = mountedContext();
      previewRuntime?.mountResultSvg?.(mounted.resultIndex, mounted.svg);
      previewRuntime?.invalidatePreviewIndexes?.('feature-label-commit');
      featureSvgActions?.attachSvgFeatureHandlers?.();
      legendActions?.setupLegendDrag?.();
    }
    if (typeof refreshLabelPresentation === 'function') {
      const refreshed = refreshLabelPresentation();
      if (refreshed && typeof refreshed.then === 'function') {
        throw new Error('Feature label presentation refresh must be synchronous.');
      }
      if (refreshed === false) return false;
    }
    return true;
  };

  const refreshPresentation = () => {
    const clicked = state.clickedFeature?.value;
    if (!clicked) return true;
    const viewModel = getFeatureLabelViewModel(clicked);
    if (!viewModel) return true;
    state.clickedFeature.value = {
      ...clicked,
      labelText: viewModel.effectiveText,
      labelSourceText: viewModel.sourceText
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
      const committed = await history.runUndoableCommand('Change feature label', async () => {
        const command = await commandBuilder({
          plan,
          state,
          catalog: currentCatalog(),
          prepareResultProjection: resultProjection,
          prepareLabelReflow,
          getMountedContext: mountedContext,
          commitMountedProjection,
          restoreMountedProjection: () => {
            previewRuntime?.clearActiveRuntime?.();
            return true;
          },
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
      clearPending({
        status: committed ? 'committed' : (commandWasNoop ? 'noop' : 'stale'),
        progress: committed
          ? { completed: plan.affectedResultKeys.length, total: plan.affectedResultKeys.length }
          : null
      });
      return committed;
    } catch (error) {
      if (sequence !== requestSequence) return false;
      clearPending({
        status: 'failed',
        progress: {
          completed: 0,
          total: plan.affectedResultKeys.length,
          message: error?.message || 'Feature label failed.'
        }
      });
      console.error('Feature label command failed.', error);
      return false;
    }
  };

  const requestFeatureLabelTextChange = async (featureLike, newText, options = {}) => {
    const feature = featureLike?.feat || featureLike;
    if (!feature) return false;
    const snapshot = plannerSnapshot();
    const explicitIndex = Number(feature?.resultIndex ?? feature?.result_index);
    const originResultIndex = Number.isInteger(explicitIndex)
      ? explicitIndex
      : Number(state.selectedResultIndex?.value ?? 0);
    const binding = resultBinding(currentCatalog(), originResultIndex);
    const token = `label-request:${snapshot.documentEpoch}:${++requestSequence}`;
    const plan = planFeatureLabelChange(snapshot, {
      targetFeatureKey: stableFeatureTargetKey(feature),
      targetFeatureKeys: options.targetFeatureKeys,
      newText,
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
        progress: { completed: 0, total: 0, message: diagnosticMessage(plan, 'Feature label is unavailable.') }
      });
      return false;
    }
    pendingPlan.value = Object.freeze(cloneJson(plan));
    planStatus.value = plan.status;
    planProgress.value = null;
    if (plan.status !== 'ready') return false;
    return executeResolvedPlan(resolveFeatureLabelPlan(plan));
  };

  const resolveFeatureLabelScope = (token, candidateId) => {
    const pending = pendingPlan.value;
    if (!pending || text(pending.token) !== text(token)) return false;
    const resolved = resolveFeatureLabelPlan(pending, candidateId);
    if (resolved?.status !== 'ready') return false;
    pendingPlan.value = Object.freeze(cloneJson(resolved));
    return executeResolvedPlan(resolved);
  };

  const cancelFeatureLabelScope = (token = null) => {
    const pending = pendingPlan.value;
    if (token && text(pending?.token) !== text(token)) return false;
    requestSequence += 1;
    clearPending();
    return true;
  };

  const requestSelectedFeatureLabelTextChange = (features, newText) => {
    const targets = Array.isArray(features) ? features.filter(Boolean) : [];
    if (targets.length === 0) return false;
    return requestFeatureLabelTextChange(targets[0], newText, {
      source: 'selection-toolbar',
      semanticScope: targets.length > 1 ? 'selected-features' : 'single',
      targetFeatureKeys: selectedFeatureKeys(targets)
    });
  };

  return Object.freeze({
    cancelFeatureLabelScope,
    getFeatureLabelViewModel,
    requestFeatureLabelTextChange,
    requestSelectedFeatureLabelTextChange,
    resolveFeatureLabelScope
  });
};
