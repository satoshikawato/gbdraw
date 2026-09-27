import {
  applyCompositionUserDeltas,
  bindCompositionMetadata,
  compositionUserDeltas
} from './composition-actions.js';
import { recordStructuralMetric } from '../../services/runtime-test-hooks.js';

const fail = (target, reason, { field = 'decorations', inputOrdinal = 1, diagnosticReason = 'DECORATION_METADATA' } = {}) => {
  const error = new Error(`Cannot preserve ${target} placement: ${reason}. The previous Result is unchanged. Reset this item's position or use Reset Layout on the previous Result, or restore the matching settings, then Generate again.`);
  error.code = 'DECORATION_CONTINUITY';
  error.stage = 'result-admission';
  error.context = { field, inputOrdinal, reason: diagnosticReason };
  throw error;
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
  if (scaleIndexes.length > 1) fail('Linear scale', 'duplicate targets', { field: 'scale', diagnosticReason: 'DECORATION_TARGET' });
  Object.entries(decorations).forEach(([role, delta]) => {
    if (delta && (delta.length !== 2 || !delta.every(Number.isFinite))) fail(role, 'non-finite offset', { field: role, diagnosticReason: 'FINITE' });
  });
  return decorations;
};

/** Capture only small deltas and source/record bindings; never retain SVG roots. */
export const captureDecorationContinuity = ({ canonical, results = [], catalog, mountedSvg = null,
  selectedResultIndex = 0, parser = globalThis.DOMParser, projectRecordIdentity } = {}) => {
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
      fail(`Result ${index + 1} decorations`, error.message, { inputOrdinal: index + 1 });
    }
    if (!Object.values(deltas).some(nonzero)) return null;
    return { deltas, identity: projectRecordIdentity(canonical, catalog?.items?.[index]?.recordKeys), index };
  }).filter(Boolean);
  if (!snapshots.length) return null;
  const oldIdentities = (catalog?.items || []).map(item => projectRecordIdentity(canonical, item.recordKeys));
  return (candidate, admission) => {
    const identities = admission.catalog.items.map(item => projectRecordIdentity(candidate, item.recordKeys));
    const transforms = identities.map(() => null);
    for (const snapshot of snapshots) {
      const roles = Object.entries(snapshot.deltas).filter(([, delta]) => nonzero(delta)).map(([role]) => role);
      const target = `Result ${snapshot.index + 1} ${roles.join('/')}`;
      const matches = identities.flatMap((identity, index) => snapshot.identity && identity
        && sameValue(snapshot.identity, identity) ? [index] : []);
      const oldMatches = oldIdentities.filter(identity => snapshot.identity && identity
        && sameValue(snapshot.identity, identity));
      if (matches.length !== 1 || oldMatches.length !== 1 || transforms[matches[0]]) {
        fail(target, 'source, region, mode, grouping or record identity is changed, unknown or ambiguous', { inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_IDENTITY' });
      }
      transforms[matches[0]] = (svg) => {
        let binding;
        try { binding = bindCompositionMetadata(svg); }
        catch (error) { fail(target, error.message, { inputOrdinal: snapshot.index + 1 }); }
        const deltas = {};
        for (const role of roles) {
          const automatic = role === 'scale' ? binding.metadata.primary.automaticTranslation
            : binding[role].metadata?.automaticTranslation;
          if (automatic && !snapshot.deltas[role].every((value, index) => Number.isFinite(value + automatic[index]))) {
            fail(target, 'non-finite placement', { field: role, inputOrdinal: snapshot.index + 1, diagnosticReason: 'FINITE' });
          }
          if (role === 'scale') {
            const indexes = binding.primary.targets.flatMap((element, index) => element.id === 'length_bar' ? [index] : []);
            if (indexes.length !== 1) fail(target, 'Linear scale target is missing or ambiguous', { field: 'scale', inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_TARGET' });
            deltas.primary = [];
            deltas.primary[indexes[0]] = snapshot.deltas.scale;
          } else {
            if (binding[role].targets.length !== 1) fail(target, `${role} target is missing or ambiguous`, { field: role, inputOrdinal: snapshot.index + 1, diagnosticReason: 'DECORATION_TARGET' });
            deltas[role] = snapshot.deltas[role];
          }
        }
        applyCompositionUserDeltas(svg, deltas);
      };
    }
    return transforms;
  };
};
