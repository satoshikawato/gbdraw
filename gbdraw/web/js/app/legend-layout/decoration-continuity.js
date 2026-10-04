import {
  applyCanvasPaddingToSvg,
  applyCompositionUserDeltas,
  bindCompositionMetadata,
  compositionUserDeltas,
  hasCanvasPadding,
  normalizeCanvasPadding
} from './composition-actions.js';
import { recordStructuralMetric } from '../../services/runtime-test-hooks.js';
import { diagnosticError } from '../../services/error-normalization.js';

const fail = ({ field = 'decorations', inputOrdinal = 1, diagnosticReason = 'DECORATION_METADATA' } = {}) => {
  throw diagnosticError('DECORATION_CONTINUITY', { field, inputOrdinal, reason: diagnosticReason }, { stage: 'result-admission' });
};

const sameValue = (left, right) => {
  if (left === right) return left !== null && left !== undefined;
  if (Array.isArray(left) && Array.isArray(right)) {
    return left.length === right.length && left.every((value, index) => sameValue(value, right[index]));
  }
  if (left && right && typeof left === 'object' && typeof right === 'object') {
    // Payload owners (File/Blob) are compared by reference only.
    if (typeof left.arrayBuffer === 'function' || typeof right.arrayBuffer === 'function') return false;
    const keys = Object.keys(left);
    return keys.length === Object.keys(right).length
      && keys.every(key => Object.hasOwn(right, key) && (
        left[key] === right[key] || sameValue(left[key], right[key])
      ));
  }
  return false;
};

const nonzero = delta => Array.isArray(delta) && delta.some(value => value !== 0);
const readDecorations = (svg) => {
  const binding = bindCompositionMetadata(svg);
  const deltas = compositionUserDeltas(svg);
  const scaleIndexes = binding.primary.targets.flatMap((target, index) => target.id === 'length_bar' ? [index] : []);
  const decorations = { legend: deltas.legend, title: deltas.title,
    scale: scaleIndexes.length === 1 ? deltas.primary[scaleIndexes[0]] : null };
  if (scaleIndexes.length > 1) fail({ field: 'scale', diagnosticReason: 'DECORATION_TARGET' });
  Object.entries(decorations).forEach(([role, delta]) => {
    if (delta && (delta.length !== 2 || !delta.every(Number.isFinite))) fail({ field: role, diagnosticReason: 'FINITE' });
  });
  return decorations;
};

/**
 * Capture only small deltas, source/record bindings, and the one canvas
 * padding; never retain SVG roots. The padding reaches every candidate Result
 * (D-09, PD-OI-064).
 */
export const captureDecorationContinuity = ({ canonical, results = [], catalog, mountedSvg = null,
  selectedResultIndex = 0, parser = globalThis.DOMParser, projectRecordIdentity, canvasPadding = null } = {}) => {
  const padding = hasCanvasPadding(canvasPadding) ? normalizeCanvasPadding(canvasPadding) : null;
  const snapshots = results.map((result, index) => {
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
  }).filter(Boolean);
  const pad = padding ? (svg) => { applyCanvasPaddingToSvg(svg, padding); } : null;
  if (!snapshots.length) {
    return pad ? (_candidate, admission) => admission.catalog.items.map(() => pad) : null;
  }
  const oldIdentities = (catalog?.items || []).map(item => projectRecordIdentity(canonical, item.recordKeys));
  return (candidate, admission) => {
    const identities = admission.catalog.items.map(item => projectRecordIdentity(candidate, item.recordKeys));
    const transforms = identities.map(() => null);
    for (const snapshot of snapshots) {
      const roles = Object.entries(snapshot.deltas).filter(([, delta]) => nonzero(delta)).map(([role]) => role);
      const matches = identities.flatMap((identity, index) => snapshot.identity && identity
        && sameValue(snapshot.identity, identity) ? [index] : []);
      const oldMatches = oldIdentities.filter(identity => snapshot.identity && identity
        && sameValue(snapshot.identity, identity));
      if (matches.length !== 1 || oldMatches.length !== 1 || transforms[matches[0]]) {
        fail({ inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_IDENTITY' });
      }
      transforms[matches[0]] = (svg) => {
        let binding;
        try { binding = bindCompositionMetadata(svg); }
        catch { fail({ inputOrdinal: snapshot.index + 1 }); }
        const deltas = {};
        for (const role of roles) {
          const automatic = role === 'scale' ? binding.metadata.primary.automaticTranslation
            : binding[role].metadata?.automaticTranslation;
          if (automatic && !snapshot.deltas[role].every((value, index) => Number.isFinite(value + automatic[index]))) {
            fail({ field: role, inputOrdinal: snapshot.index + 1, diagnosticReason: 'FINITE' });
          }
          if (role === 'scale') {
            const indexes = binding.primary.targets.flatMap((element, index) => element.id === 'length_bar' ? [index] : []);
            if (indexes.length !== 1) fail({ field: 'scale', inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_TARGET' });
            deltas.primary = [];
            deltas.primary[indexes[0]] = snapshot.deltas.scale;
          } else {
            if (binding[role].targets.length !== 1) fail({ field: role, inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_TARGET' });
            deltas[role] = snapshot.deltas[role];
          }
        }
        applyCompositionUserDeltas(svg, deltas);
      };
    }
    return pad ? transforms.map((transform) => (transform
      ? (svg, context) => { transform(svg, context); pad(svg); }
      : pad)) : transforms;
  };
};
