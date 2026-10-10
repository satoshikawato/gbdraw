// @ts-check
/** @import { DrawingState } from '../../state.js' */
import {
  applyFeatureVisibilityOverrideChanges,
  buildFeatureVisibilityChanges,
  createDefaultFeatureVisibilityRule,
  featureDrawnContext,
  featureVisibilityQualifierSuggestions,
  getFeatureVisibilityOverride,
  normalizeFeatureVisibilityRule,
  normalizeVisibilityMode,
  removeEditorQualifierFeatureVisibilityRule,
  requestFeatureVisibilityRules,
  resolveFeatureDrawn,
  resultLegendSources,
  sameLegendSources,
  serializeFeatureVisibilityRules,
  setFeatureVisibilityOverride,
  upsertEditorQualifierFeatureVisibilityRule,
} from '../../services/feature-visibility.js';
import { isInternalProteinDisplayId } from '../../services/feature-utils.js';
import { downloadTextFile } from '../../services/text-download.js';
import { normalizeUserFacingError } from '../../utils/error-normalization.js';
import { resolveUniqueOrthogroupMemberForFeature } from '../../services/feature-identity.js';
import { featureIdentityKeyOf } from '../../services/feature-placement.js';
import { resultCatalogFeatures, stableFeatureOverrideKey } from '../../services/feature-catalog.js';

/**
 * The label owner's reaction the composition root registers once that owner exists.
 * @typedef {object} VisibilityActionsPorts
 * @property {(options?: { reflow?: boolean, rerender?: boolean, labels?: boolean }) => boolean} applyFeatureVisibilityToLabels
 *   Projects the stored visibility overrides (`labels`: the whole label intent) onto the labels and queues the label reflow.
 */

/**
 * @typedef {object} FeatureVisibilityActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(options: { domains: readonly string[] }) => boolean} showEditorIntent
 *   The root's port of the editor intent onto the displayed Result (R13), which
 *   compiles Feature visibility with `resolveFeatureDrawn`, as Generate draws it.
 * @property {VisibilityActionsPorts} ports
 * @property {(index: number) => any} selectResult The preview owner's selection of the Result a History command targets.
 * @property {import('../rule-matching.js').RulePreparation | null} [rulePreparation]
 * @property {() => ({ diagramOptions?: Record<string, any> } | null)} [getCommittedRequest]
 *   The committed canonical request (Python owns the option fields, R7).
 */

// The paint domain this owner shows through the port.
const VISIBILITY_DOMAINS = Object.freeze(['featureVisibility']);

/** @param {FeatureVisibilityActionsOptions} options */
export const createFeatureVisibilityActions = ({
  state,
  // R13: the root's port of the editor intent onto the displayed Result.
  showEditorIntent,
  // R13: the composition root registers `applyFeatureVisibilityToLabels` once
  // the label owner exists; this owner only calls it.
  ports,
  selectResult,
  rulePreparation = null,
  getCommittedRequest = () => null
}) => {
  // The preview owner selects the Result a History command targets (R13).
  if (typeof selectResult !== 'function') {
    throw new Error('createFeatureVisibilityActions requires PreviewRuntime Result selection.');
  }
  const {
    clickedFeature,
    extractedFeatures,
    orthogroups,
    featureVisibilityScopeDialog,
    resultGenerationKey,
    results,
    selectedResultIndex,
    svgContainer
  } = state;

  const ruleFields = new Set(['recordId', 'featureType', 'qualifier', 'value', 'action']);

  const normalizeText = (value) => String(value ?? '').trim();

  const firstText = (...values) => {
    for (const value of values) {
      if (Array.isArray(value)) {
        const found = firstText(...value);
        if (found) return found;
        continue;
      }
      const text = normalizeText(value);
      if (text) return text;
    }
    return '';
  };

  const getFeatureType = (feat) => normalizeText(feat?.type || feat?.featureType || feat?.feature_type);

  const getQualifierValue = (feat, key) => {
    const normalizedKey = normalizeText(key).toLowerCase();
    if (!normalizedKey) return '';
    const qualifiers = feat?.qualifiers && typeof feat.qualifiers === 'object' ? feat.qualifiers : {};
    return firstText(feat?.[normalizedKey], qualifiers[normalizedKey]);
  };

  const orthogroupIdFor = (source) => {
    const ids = new Set([
      source?.id,
      source?.orthogroupId,
      source?.orthogroup_id
    ].map(normalizeText).filter(Boolean));
    return ids.size === 1 ? ids.values().next().value : '';
  };

  const findOrthogroup = (orthogroupId) => {
    const id = normalizeText(orthogroupId);
    if (!id) return null;
    const matches = (Array.isArray(orthogroups?.value) ? orthogroups.value : [])
      .filter((group) => orthogroupIdFor(group) === id);
    return matches.length === 1 ? matches[0] : null;
  };

  const uniqueFeaturesBySvgId = (features) => {
    const seen = new Set();
    return features.filter((feat) => {
      const svgId = normalizeText(feat?.svg_id ?? feat?.svgId ?? feat?.id);
      if (!svgId || seen.has(svgId)) return false;
      seen.add(svgId);
      return true;
    });
  };

  const getOrthogroupMemberFeatures = (feat) => {
    const orthogroupId = orthogroupIdFor({
      orthogroupId: feat?.orthogroupId,
      orthogroup_id: feat?.orthogroup_id
    });
    if (!orthogroupId) return [];
    const group = findOrthogroup(orthogroupId);
    if (!group) return [];
    const members = Array.isArray(group.members) ? group.members : [];
    if (!resolveUniqueOrthogroupMemberForFeature(feat, members)) return [];

    const featuresByMember = new Map();
    (Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []).forEach((candidate) => {
      const candidateGroupId = orthogroupIdFor({
        orthogroupId: candidate?.orthogroupId,
        orthogroup_id: candidate?.orthogroup_id
      });
      if (candidateGroupId !== orthogroupId) return;
      const member = resolveUniqueOrthogroupMemberForFeature(candidate, members);
      if (!member) return;
      const matches = featuresByMember.get(member) || [];
      matches.push(candidate);
      featuresByMember.set(member, matches);
    });
    return members.flatMap((member) => {
      const matches = featuresByMember.get(member) || [];
      return matches.length === 1 ? matches : [];
    });
  };

  const buildVisibilityScopes = (feat) => {
    const featureType = getFeatureType(feat);
    const scopes = [{
      id: 'feature',
      label: 'This feature',
      description: 'One rule for only this feature.'
    }];
    const orthogroupMembers = getOrthogroupMemberFeatures(feat);
    if (orthogroupMembers.length > 1) {
      scopes.push({
        id: 'orthogroup',
        label: `Current similarity-group members (${orthogroupMembers.length})`,
        description: 'One rule per current similarity-group member.',
        features: orthogroupMembers
      });
    }
    const product = getQualifierValue(feat, 'product');
    if (featureType && product) {
      scopes.push({
        id: 'product',
        label: `Exact product: ${product}`,
        description: 'Exact product qualifier rule.',
        featureType,
        qualifier: 'product',
        value: product
      });
    }
    const proteinIdCandidate = firstText(
      feat?.sourceProteinId,
      feat?.source_protein_id,
      getQualifierValue(feat, 'protein_id'),
      feat?.proteinId,
      feat?.protein_id
    );
    const proteinId = isInternalProteinDisplayId(proteinIdCandidate) ? '' : proteinIdCandidate;
    if (featureType && proteinId) {
      scopes.push({
        id: 'protein_id',
        label: `Exact protein ID: ${proteinId}`,
        description: 'Exact protein_id qualifier rule.',
        featureType,
        qualifier: 'protein_id',
        value: proteinId
      });
    }
    return scopes;
  };

  // The one site that reaches the label owner (R13): it hides a hidden
  // feature's label with it, as Generate does, and queues the label reflow
  // unless the caller declines it (F-3); `rerender` asks for the automatic
  // rerender, which draws a feature the Result does not draw (R-5).
  const followLabels = (options) => ports.applyFeatureVisibilityToLabels(options);

  // The visibility table holds rules only; per-feature edits are identity rows.
  /** @param {DrawingState} drawing */
  const visibilityRuleRows = (drawing) => (
    Array.isArray(drawing.featureVisibilityRules?.value)
      ? drawing.featureVisibilityRules.value
      : drawing.featureVisibilityManualRules.map((rule) => normalizeFeatureVisibilityRule(rule))
  );

  const sameFeature = (left, right) => {
    const key = featureIdentityKeyOf(left);
    return Boolean(key) && key === featureIdentityKeyOf(right);
  };

  const isRuleScope = (scope) => scope?.id === 'product' || scope?.id === 'protein_id';

  // OV-42 (Owner decision 2026-10-06, option A): whether a Result's Legend
  // derives from other features than Python drew at the last render: the
  // Legend source of the features each Result's catalog draws against the
  // source of the features the resolver draws now, under the current rules.
  // Python redraws the Legend in the automatic rerender.
  const legendSourceChanged = (context) => !sameLegendSources(
    resultLegendSources(state, context, { asRendered: true }),
    resultLegendSources(state, context)
  );

  // The live projection of what the label rerender and Generate draw (R3):
  // the port shows the resolver's answer for every feature the displayed
  // Result draws, and an unknown answer keeps Python's. A feature the Result
  // does not draw and the resolver draws needs geometry the Result does not
  // have, so the action asks for the rerender, which draws what Generate draws
  // (R-5); so does an unknown answer for a feature the action names
  // (`targeted`), and, unless the caller declines it (`legend`), a changed
  // Legend source of any Result (OV-42). The resolver reads Python's rule
  // matches, so an action prepares them first (`prepareDrawn`).
  /** @param {DrawingState} drawing */
  const needsRerenderFor = (drawing, features, { targeted = false, legend = true } = {}) => {
    const context = featureDrawnContext(drawing, { diagramOptions: getCommittedRequest()?.diagramOptions });
    const catalogFeatures = resultCatalogFeatures(state);
    return features.some((feature) => {
      if (!catalogFeatures || catalogFeatures.renderedByIdentity.has(stableFeatureOverrideKey(feature))) return false;
      const drawn = resolveFeatureDrawn(feature, context);
      return drawn === true || (targeted && drawn === null);
    }) || (legend && legendSourceChanged(context));
  };
  // Every catalog feature of the displayed Result, as drawn where it is drawn.
  const displayedFeatures = () => {
    const catalogFeatures = resultCatalogFeatures(state);
    return catalogFeatures
      ? catalogFeatures.biological.map((feature) => (
        catalogFeatures.renderedByIdentity.get(stableFeatureOverrideKey(feature)) || feature
      ))
      : uniqueFeaturesBySvgId(Array.isArray(extractedFeatures.value) ? extractedFeatures.value : []);
  };
  // Shows the projection, unless the caller shows it in its own compile
  // (`show: false`); returns whether the Result changed (or is to change) and
  // whether the rerender must draw a feature or the Legend.
  /**
   * @param {DrawingState} drawing
   * @param {Record<string, any>[]} features
   * @param {{ targeted?: boolean, legend?: boolean, show?: boolean }} [projection]
   */
  const projectDrawn = (drawing, features, { show = true, ...projection } = {}) => ({
    updated: show ? showEditorIntent({ domains: VISIBILITY_DOMAINS }) : true,
    needsRerender: needsRerenderFor(drawing, features, projection)
  });
  // Also the reaction the popup and Label On ask this owner for: the matches
  // `resolveFeatureDrawn` reads are this owner's read model (R13).
  /**
   * @param {{ strict?: boolean }} [options]
   * @returns {boolean | Promise<boolean | { error: any }>}
   */
  const prepareDrawn = (options) => rulePreparation?.prepareDrawn?.(options) ?? true;

  /** @param {DrawingState} drawing */
  const updateClickedFeatureVisibilityFromRules = (drawing, features) => {
    const clicked = clickedFeature.value?.feat;
    if (!clicked || !features.some((feature) => sameFeature(feature, clicked))) return;
    clickedFeature.value.featureVisibility = getFeatureVisibilityOverride(drawing.featureOverrides, clicked);
  };

  const nextFrame = () => /** @type {Promise<void>} */ (new Promise((resolve) => {
    if (typeof window !== 'undefined' && typeof window.requestAnimationFrame === 'function') {
      window.requestAnimationFrame(() => resolve());
    } else {
      setTimeout(resolve, 0);
    }
  }));

  const ensureCommandTargetResult = async (resultIndex, generationKey) => {
    if (String(resultGenerationKey?.value ?? '') !== String(generationKey ?? '')) return false;
    const targetIndex = Number(resultIndex);
    if (!Number.isInteger(targetIndex) || targetIndex < 0) return false;
    const resultCount = Array.isArray(results?.value) ? results.value.length : 0;
    if (targetIndex >= resultCount) return false;
    if (Number(selectedResultIndex?.value || 0) !== targetIndex) {
      selectResult(targetIndex);
      await nextFrame();
    }
    return Boolean(svgContainer?.value?.querySelector?.('svg'));
  };

  const buildSelectedFeaturesVisibilityCommand = (features, modeRaw) => {
    const drawing = state.activeDrawing();
    const targetFeatures = uniqueFeaturesBySvgId(Array.isArray(features) ? features : [])
      .filter((feature) => featureIdentityKeyOf(feature));
    if (targetFeatures.length === 0) return null;

    const changes = buildFeatureVisibilityChanges(targetFeatures, modeRaw, drawing.featureOverrides);
    if (changes.length === 0) return null;
    const changedKeys = new Set(changes.map(featureIdentityKeyOf));
    const changedFeatures = targetFeatures.filter((feature) => changedKeys.has(featureIdentityKeyOf(feature)));
    const commandResultIndex = Number(selectedResultIndex?.value || 0);
    const commandGenerationKey = String(resultGenerationKey?.value ?? '');

    const applyChangeSet = async (direction) => {
      if (!(await ensureCommandTargetResult(commandResultIndex, commandGenerationKey))) return false;
      await prepareDrawn();
      const useAfter = direction === 'apply';
      const overrideChanged = changes.some((change) => change.before !== change.after);
      applyFeatureVisibilityOverrideChanges(
        drawing.featureOverrides,
        changes.map((change) => ({
          ...change,
          mode: useAfter ? change.after : change.before
        }))
      );
      const { updated, needsRerender } = projectDrawn(drawing, changedFeatures, { targeted: true });
      if (!updated && !overrideChanged) return false;
      updateClickedFeatureVisibilityFromRules(drawing, targetFeatures);
      followLabels({ rerender: needsRerender });
      return true;
    };

    return {
      label: 'Change selected feature visibility',
      resultIndex: commandResultIndex,
      resultGenerationKey: commandGenerationKey,
      changes,
      apply: () => applyChangeSet('apply'),
      revert: () => applyChangeSet('revert'),
      estimateBytes: () => JSON.stringify({ resultIndex: commandResultIndex, resultGenerationKey: commandGenerationKey, changes }).length * 2
    };
  };

  const showClickedFeatureVisibility = (targetFeatures, mode) => {
    const clicked = clickedFeature.value?.feat;
    if (clicked && targetFeatures.some((feature) => sameFeature(feature, clicked))) {
      clickedFeature.value.featureVisibility = mode;
    }
  };

  // A feature or similarity-group scope: the identity rows of its features.
  // Returns whether the rerender must draw one of them.
  /** @param {DrawingState} drawing */
  const applyFeatureVisibilityScope = (drawing, feat, modeRaw, scope) => {
    const nextMode = normalizeVisibilityMode(modeRaw);
    const targetFeatures = scope?.id === 'orthogroup'
      ? uniqueFeaturesBySvgId(scope.features || [])
      : [feat];
    targetFeatures.forEach((targetFeat) => {
      setFeatureVisibilityOverride(drawing.featureOverrides, targetFeat, nextMode);
    });
    const { needsRerender } = projectDrawn(drawing, targetFeatures, { targeted: true });
    showClickedFeatureVisibility(targetFeatures, nextMode);
    return needsRerender;
  };

  // The live projection of the visibility rules and edits onto every
  // displayed feature (R3), with the same resolver as the visibility actions,
  // so Undo and Redo of a scoped hide leave the preview as the action left it
  // (FE-04). It prepares Python's matches first. A rule table that Generate
  // rejects (an invalid regex) changes nothing, as a failed Generate keeps its
  // Result, and reports Generate's error until a projection reads a table that
  // Generate accepts. A later projection supersedes one that still waits for
  // its matches. Resolves to `projectDrawn`'s answer, or null when it projects
  // nothing.
  let projectionRun = 0;
  let ruleTableError = null;
  /** @param {DrawingState} drawing */
  const runProjection = async (drawing, { legend = true, show = true } = {}) => {
    const run = ++projectionRun;
    const prepared = /** @type {any} */ (await prepareDrawn());
    if (run !== projectionRun) return null;
    if (prepared?.error) {
      if (!state.errorLog) return null;
      state.errorLog.value = normalizeUserFacingError(prepared.error, { operation: 'evaluateRules', stage: 'rule-validation' });
      ruleTableError = state.errorLog.value; // As the error log holds it (a reactive copy).
      return null;
    }
    if (ruleTableError && state.errorLog?.value === ruleTableError) state.errorLog.value = null;
    ruleTableError = null;
    return projectDrawn(drawing, displayedFeatures(), { legend, show });
  };

  // The one projection of this domain (R3), which History apply, the display
  // of a Result, and Load Feature Edits TSV call through
  // `projectMountedEditorIntent`. The label of a feature it hides or shows
  // follows the feature (F-3, OV-35). A History step or a loaded table
  // (`rerender`) that draws a feature the Result does not draw, or changes a
  // Legend source, reruns the rerender, as the action did; a Result display
  // does not, so a feature Python does not draw cannot repeat it. A loaded
  // table (`reflow`) also places the labels, as a visibility edit does. A
  // caller that also projects the label intent (`labels`) has the label owner
  // project it in the same follow, so the edit queues one label request.
  // Returns whether the Result changed. A caller that shows the visibility in
  // its own compile passes `show: false`.
  const projectFeatureVisibility = async ({ rerender = false, reflow = false, show = true, labels = false } = {}) => {
    const drawing = state.activeDrawing();
    const projection = await runProjection(drawing, { legend: rerender, show });
    const drawsMissing = rerender && Boolean(projection?.needsRerender);
    if (labels || (projection && (projection.updated || reflow || drawsMissing))) {
      followLabels({ reflow: reflow && Boolean(projection), rerender: drawsMissing, labels });
    }
    return Boolean(projection?.updated);
  };

  // OV-19 (PD-OI-066, R10): the one transition of a visibility rule edit, from
  // the Features panel or the popup's product and protein ID scopes. `edit`
  // changes a copy of the rules; the transition writes it, then projects the
  // rules and the labels of the features they hide or show, as the label
  // rerender and Generate draw them; a rule that draws a feature the Result
  // does not draw asks for the rerender (R-5). The draft keeps a rule that
  // Generate rejects. An edit that leaves the rules Generate reads as they
  // were (a rule added, or a field set before its value) projects nothing.
  /** @param {DrawingState} drawing */
  const editFeatureVisibilityRules = async (drawing, edit) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const rules = [...drawing.featureVisibilityManualRules];
    if (edit(rules) === false) return false;
    const requested = JSON.stringify(requestFeatureVisibilityRules(drawing.featureVisibilityManualRules));
    drawing.featureVisibilityManualRules.splice(0, drawing.featureVisibilityManualRules.length, ...rules);
    if (JSON.stringify(requestFeatureVisibilityRules(rules)) === requested) return true;
    const projection = await runProjection(drawing);
    if (projection?.updated || projection?.needsRerender) followLabels({ rerender: projection.needsRerender });
    return true;
  };

  const clearFeatureVisibilityScopeDialog = ({ restorePrevious = false } = {}) => {
    if (restorePrevious && clickedFeature.value) {
      clickedFeature.value.featureVisibility = featureVisibilityScopeDialog.previousMode || 'default';
    }
    featureVisibilityScopeDialog.show = false;
    featureVisibilityScopeDialog.feat = null;
    featureVisibilityScopeDialog.mode = 'default';
    featureVisibilityScopeDialog.previousMode = 'default';
    featureVisibilityScopeDialog.scopes = [];
  };

  /** @param {DrawingState} drawing */
  const setFeatureVisibilityIn = (drawing, feat, modeRaw, options = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!featureIdentityKeyOf(feat)) return false;

    const triggerReflow = options.triggerReflow !== false;
    const scope = options.scope || { id: 'feature' };
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getFeatureVisibilityOverride(drawing.featureOverrides, feat);

    const needsRerender = applyFeatureVisibilityScope(drawing, feat, nextMode, scope);

    if (previousMode !== nextMode) {
      followLabels({ reflow: triggerReflow, rerender: needsRerender });
    }

    return previousMode !== nextMode;
  };

  /**
   * @param {Record<string, any>} feat
   * @param {string} modeRaw
   * @param {Record<string, any>} [options]
   */
  const setFeatureVisibility = (feat, modeRaw, options = {}) => {
    const drawing = state.activeDrawing();
    return setFeatureVisibilityIn(drawing, feat, modeRaw, options);
  };

  const setSelectedFeaturesVisibility = async (features, modeRaw) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const command = buildSelectedFeaturesVisibilityCommand(features, modeRaw);
    if (!command) return false;
    return command.apply();
  };

  const updateClickedFeatureVisibility = async (modeRaw) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!clickedFeature.value?.feat) return false;
    const feat = clickedFeature.value.feat;
    const scopes = buildVisibilityScopes(feat);
    const nextMode = normalizeVisibilityMode(modeRaw);
    const previousMode = getFeatureVisibilityOverride(drawing.featureOverrides, feat);
    if (scopes.length <= 1) {
      await prepareDrawn();
      return setFeatureVisibilityIn(drawing, feat, nextMode, { triggerReflow: true, scope: scopes[0] });
    }
    featureVisibilityScopeDialog.show = true;
    featureVisibilityScopeDialog.feat = feat;
    featureVisibilityScopeDialog.mode = nextMode;
    featureVisibilityScopeDialog.previousMode = previousMode;
    featureVisibilityScopeDialog.scopes = scopes;
    return false;
  };

  const handleFeatureVisibilityScopeChoice = async (scopeId) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (scopeId === 'cancel' || !featureVisibilityScopeDialog.show) {
      clearFeatureVisibilityScopeDialog({ restorePrevious: true });
      return false;
    }
    const feat = featureVisibilityScopeDialog.feat;
    const nextMode = featureVisibilityScopeDialog.mode;
    const scope = (featureVisibilityScopeDialog.scopes || []).find((entry) => entry.id === scopeId);
    if (!feat || !scope) {
      clearFeatureVisibilityScopeDialog({ restorePrevious: true });
      return false;
    }
    const previousMode = featureVisibilityScopeDialog.previousMode;
    clearFeatureVisibilityScopeDialog();
    if (isRuleScope(scope)) {
      const ruleInput = { featureType: scope.featureType, qualifier: scope.qualifier, value: scope.value, label: scope.label };
      await editFeatureVisibilityRules(drawing, (rules) => {
        if (nextMode === 'default') removeEditorQualifierFeatureVisibilityRule(rules, ruleInput);
        else upsertEditorQualifierFeatureVisibilityRule(rules, ruleInput, nextMode);
      });
      showClickedFeatureVisibility([feat], nextMode);
      return previousMode !== nextMode;
    }
    await prepareDrawn();
    const needsRerender = applyFeatureVisibilityScope(drawing, feat, nextMode, scope);
    if (previousMode !== nextMode) followLabels({ rerender: needsRerender });
    return previousMode !== nextMode;
  };

  const setFeatureVisibilityRuleField = (index, field, value) => {
    const drawing = state.activeDrawing();
    return editFeatureVisibilityRules(drawing, (rules) => {
      if (!ruleFields.has(field) || !rules[index]) return false;
      rules[index] = normalizeFeatureVisibilityRule({ ...rules[index], [field]: value });
    });
  };

  /** @param {DrawingState} drawing */
  const moveFeatureVisibilityRule = (drawing, index, offset) => editFeatureVisibilityRules(drawing, (rules) => {
    const target = index + offset;
    if (!rules[index] || target < 0 || target >= rules.length) return false;
    rules.splice(target, 0, ...rules.splice(index, 1));
  });

  const addFeatureVisibilityRule = () => {
    const drawing = state.activeDrawing();
    return editFeatureVisibilityRules(drawing, (rules) => {
      rules.push(createDefaultFeatureVisibilityRule());
    });
  };

  const removeFeatureVisibilityRule = (index) => {
    const drawing = state.activeDrawing();
    return editFeatureVisibilityRules(drawing, (rules) => {
      if (!rules[index]) return false;
      rules.splice(index, 1);
    });
  };

  const downloadFeatureVisibilityRulesTsv = () => {
    const drawing = state.activeDrawing();
    const text = serializeFeatureVisibilityRules(visibilityRuleRows(drawing));
    if (!text.trim()) {
      alert('No valid feature visibility rules to export.');
      return;
    }
    downloadTextFile('gbdraw_feature_visibility_table.tsv', text);
  };

  const featureVisibilityRuleDetail = (rule) => {
    const normalized = normalizeFeatureVisibilityRule(rule);
    if (normalized.source === 'editor' && normalized.featureId) {
      const qualifier = normalized.qualifier.toLowerCase() === 'hash'
        ? 'hash fallback'
        : normalized.qualifier.toLowerCase();
      return normalized.label ? `${normalized.label} (${qualifier})` : `${normalized.featureId} (${qualifier})`;
    }
    if (normalized.qualifier.toLowerCase() === 'hash') return normalized.value;
    return '';
  };

  return {
    addFeatureVisibilityRule,
    downloadFeatureVisibilityRulesTsv,
    featureVisibilityQualifierSuggestions,
    featureVisibilityRuleDetail,
    buildSelectedFeaturesVisibilityCommand,
    handleFeatureVisibilityScopeChoice,
    moveFeatureVisibilityRuleDown: (index) => {
      const drawing = state.activeDrawing();
      return moveFeatureVisibilityRule(drawing, index, 1);
    },
    moveFeatureVisibilityRuleUp: (index) => {
      const drawing = state.activeDrawing();
      return moveFeatureVisibilityRule(drawing, index, -1);
    },
    prepareDrawnFeatureMatches: prepareDrawn,
    projectFeatureVisibility,
    removeFeatureVisibilityRule,
    setFeatureVisibility,
    setSelectedFeaturesVisibility,
    setFeatureVisibilityRuleField,
    updateClickedFeatureVisibility
  };
};
