// @ts-check
import {
  applyCanvasPaddingToSvg,
  applyCompositionUserDeltas,
  bindCompositionMetadata,
  compositionUserDeltas,
  hasCanvasPadding,
  normalizeCanvasPadding
} from './composition-actions.js';
import { recordStructuralMetric } from '../../services/runtime-test-hooks.js';
import { LENGTH_BAR_ID } from '../../services/drawn-groups.js';
import { diagnosticError } from '../../utils/error-normalization.js';
/** @import { SvgResultTransform } from '../../services/svg-result-ingestion.js' */
/** @import { CompositionBinding } from './composition-actions.js' */

/** @typedef {'legend' | 'title' | 'scale'} DecorationRole */
/** @typedef {Record<DecorationRole, number[] | null>} DecorationDeltas The user offset of each decoration role. */

/**
 * The candidate Results' catalog, which a plan reads for the number of Results and each one's record keys.
 * @typedef {{ catalog: { items: Array<{ recordKeys?: unknown }> } }} DecorationAdmission
 */
/**
 * One transform per candidate Result (null where it keeps the offsets the generator drew).
 * @typedef {(candidate: Record<string, unknown>, admission: DecorationAdmission) => Array<SvgResultTransform | null>} DecorationTransformPlan
 */

/** @type {(options?: { field?: string, inputOrdinal?: number, diagnosticReason?: string }) => never} */
const fail = ({ field = 'decorations', inputOrdinal = 1, diagnosticReason = 'DECORATION_METADATA' } = {}) => {
  throw diagnosticError('DECORATION_CONTINUITY', { field, inputOrdinal, reason: diagnosticReason }, { stage: 'result-admission' });
};

/**
 * @param {unknown} left
 * @param {unknown} right
 * @returns {boolean}
 */
const sameValue = (left, right) => {
  if (left === right) return left !== null && left !== undefined;
  if (Array.isArray(left) && Array.isArray(right)) {
    return left.length === right.length && left.every((value, index) => sameValue(value, right[index]));
  }
  if (left && right && typeof left === 'object' && typeof right === 'object') {
    // Payload owners (File/Blob) are compared by reference only.
    const leftFields = /** @type {Record<string, unknown>} */ (left); // a non-null object, read by its own keys
    const rightFields = /** @type {Record<string, unknown>} */ (right);
    if (typeof leftFields.arrayBuffer === 'function' || typeof rightFields.arrayBuffer === 'function') return false;
    const keys = Object.keys(leftFields);
    return keys.length === Object.keys(rightFields).length
      && keys.every(key => Object.hasOwn(rightFields, key) && (
        leftFields[key] === rightFields[key] || sameValue(leftFields[key], rightFields[key])
      ));
  }
  return false;
};

/** @param {unknown} delta */
const nonzero = delta => Array.isArray(delta) && delta.some(value => value !== 0);
/**
 * @param {Element} svg
 * @returns {DecorationDeltas}
 */
const readDecorations = (svg) => {
  const binding = bindCompositionMetadata(svg);
  const deltas = compositionUserDeltas(svg);
  const scaleIndexes = binding.primary.targets.flatMap((target, index) => target.id === LENGTH_BAR_ID ? [index] : []);
  const decorations = { legend: deltas.legend, title: deltas.title,
    scale: scaleIndexes.length === 1 ? deltas.primary[scaleIndexes[0]] : null };
  if (scaleIndexes.length > 1) fail({ field: 'scale', diagnosticReason: 'DECORATION_TARGET' });
  Object.entries(decorations).forEach(([role, delta]) => {
    if (delta && (delta.length !== 2 || !delta.every(Number.isFinite))) fail({ field: role, diagnosticReason: 'FINITE' });
  });
  return decorations;
};

/**
 * @typedef {object} DecorationContinuityOptions
 * @property {Record<string, any> | null | undefined} canonical The committed canonical session (Python owns its fields, R7); null before one is committed.
 * @property {Array<{ content: string }>} [results] The displayed Results.
 * @property {{ items?: Array<{ recordKeys?: any }> } | null} [catalog] The feature catalog of the displayed Results.
 * @property {Element | null} [mountedSvg] The mounted root of the selected Result.
 * @property {number} [selectedResultIndex]
 * @property {typeof DOMParser} [parser]
 * @property {(source: Record<string, any> | null | undefined, recordKeys: any) => any} projectRecordIdentity
 *   The root's projection of a record's identity out of a canonical request or session.
 * @property {Record<string, any> | null} [canvasPadding]
 */

/**
 * Capture only small deltas, source/record bindings, and the one canvas
 * padding; never retain SVG roots. The padding reaches every candidate Result
 * (D-09, PD-OI-064).
 * @param {DecorationContinuityOptions} options
 */
export const captureDecorationContinuity = ({ canonical, results = [], catalog, mountedSvg = null,
  selectedResultIndex = 0, parser = globalThis.DOMParser, projectRecordIdentity, canvasPadding = null } = /** @type {DecorationContinuityOptions} */ ({})) => {
  const padding = hasCanvasPadding(canvasPadding) ? normalizeCanvasPadding(canvasPadding) : null;
  // filter(Boolean) drops the null that `map` returns for a Result with no decoration offset.
  const snapshots = /** @type {Array<{ deltas: Record<string, any>, identity: any, index: number }>} */ (results.map((result, index) => {
    /** @type {DecorationDeltas} */
    let deltas;
    try {
      let svg = index === selectedResultIndex ? mountedSvg : null;
      if (!svg) {
        recordStructuralMetric('applicationSvgParseCount', 1, { phase: 'decoration-snapshot', resultIndex: index });
        svg = new parser().parseFromString(result.content, 'image/svg+xml').documentElement;
      }
      deltas = readDecorations(svg);
    }
    catch (error) {
      if (error.code === 'DECORATION_CONTINUITY') {
        error.context.inputOrdinal = index + 1;
        throw error;
      }
      fail({ inputOrdinal: index + 1 });
    }
    if (!Object.values(deltas).some(nonzero)) return null;
    return { deltas, identity: projectRecordIdentity(canonical, catalog?.items?.[index]?.recordKeys), index };
  }).filter(Boolean));
  const pad = padding ? (/** @type {Element} */ svg) => { applyCanvasPaddingToSvg(svg, padding); } : null;
  if (!snapshots.length) {
    return pad ? /** @type {DecorationTransformPlan} */ ((_candidate, admission) => admission.catalog.items.map(() => pad)) : null;
  }
  const oldIdentities = (catalog?.items || []).map(item => projectRecordIdentity(canonical, item.recordKeys));
  return /** @type {DecorationTransformPlan} */ ((candidate, admission) => {
    const identities = admission.catalog.items.map(item => projectRecordIdentity(candidate, item.recordKeys));
    /** @type {Array<SvgResultTransform | null>} */
    const transforms = identities.map(() => null);
    for (const snapshot of snapshots) {
      // The keys of `snapshot.deltas` are exactly the decoration roles (`readDecorations`).
      const roles = /** @type {DecorationRole[]} */ (
        Object.entries(snapshot.deltas).filter(([, delta]) => nonzero(delta)).map(([role]) => role)
      );
      const matches = identities.flatMap((identity, index) => snapshot.identity && identity
        && sameValue(snapshot.identity, identity) ? [index] : []);
      const oldMatches = oldIdentities.filter(identity => snapshot.identity && identity
        && sameValue(snapshot.identity, identity));
      if (matches.length !== 1 || oldMatches.length !== 1 || transforms[matches[0]]) {
        fail({ inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_IDENTITY' });
      }
      transforms[matches[0]] = (svg) => {
        /** @type {CompositionBinding} */
        let binding;
        try { binding = bindCompositionMetadata(svg); }
        catch { fail({ inputOrdinal: snapshot.index + 1 }); }
        /** @type {{ primary?: Array<number[] | null>, legend?: number[] | null, title?: number[] | null }} */
        const deltas = {};
        for (const role of roles) {
          const automatic = role === 'scale' ? binding.metadata.primary.automaticTranslation
            : binding[role].metadata?.automaticTranslation;
          // `roles` holds only roles whose delta `nonzero` accepted: a non-empty array.
          const delta = /** @type {number[]} */ (snapshot.deltas[role]);
          if (automatic && !delta.every((value, index) => Number.isFinite(value + automatic[index]))) {
            fail({ field: role, inputOrdinal: snapshot.index + 1, diagnosticReason: 'FINITE' });
          }
          if (role === 'scale') {
            const indexes = binding.primary.targets.flatMap((element, index) => element.id === LENGTH_BAR_ID ? [index] : []);
            if (indexes.length !== 1) fail({ field: 'scale', inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_TARGET' });
            deltas.primary = [];
            deltas.primary[indexes[0]] = snapshot.deltas.scale;
          } else {
            if (binding[role].targets.length !== 1) fail({ field: role, inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_TARGET' });
            deltas[role] = delta;
          }
        }
        applyCompositionUserDeltas(svg, deltas);
      };
    }
    return pad ? transforms.map((transform) => (transform
      ? (svg, context) => { transform(svg, context); pad(svg); }
      : pad)) : transforms;
  });
};
