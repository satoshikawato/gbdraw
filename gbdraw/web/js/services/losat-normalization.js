// @ts-check
import { ORTHOGROUP_INFERENCE_VERSION } from '../contracts/schema.generated.js';

/** @param {unknown} value */
export const normalizeCollinearAnchorMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  /** @type {Record<string, string>} */
  const aliases = {
    all_hits: 'all',
    raw: 'all',
    top_n: 'all',
    topn: 'all',
    one2one: 'one_to_one',
    one_to_ones: 'one_to_one',
    mutual_best: 'one_to_one',
    top1: 'one_to_one',
    top_1: 'one_to_one',
    reciprocal_best: 'rbh',
    strict_rbh: 'rbh'
  };
  const resolved = aliases[normalized] || normalized;
  return ['all', 'one_to_one', 'rbh'].includes(resolved) ? resolved : 'rbh';
};

/** @param {unknown} value */
export const normalizeOrthogroupMembershipMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  /** @type {Record<string, string>} */
  const aliases = {
    legacy: ORTHOGROUP_INFERENCE_VERSION,
    rbh: ORTHOGROUP_INFERENCE_VERSION,
    rbh_only: ORTHOGROUP_INFERENCE_VERSION,
    merge: ORTHOGROUP_INFERENCE_VERSION,
    family: ORTHOGROUP_INFERENCE_VERSION,
    family_merge: ORTHOGROUP_INFERENCE_VERSION,
    local_split: ORTHOGROUP_INFERENCE_VERSION,
    density_split: ORTHOGROUP_INFERENCE_VERSION,
    outparalog_split: ORTHOGROUP_INFERENCE_VERSION,
    distribution_split: ORTHOGROUP_INFERENCE_VERSION,
    orthogroups: ORTHOGROUP_INFERENCE_VERSION,
    anchor_core: ORTHOGROUP_INFERENCE_VERSION
  };
  const resolved = aliases[normalized] || normalized;
  return resolved === ORTHOGROUP_INFERENCE_VERSION ? resolved : ORTHOGROUP_INFERENCE_VERSION;
};

/** @param {unknown} value */
export const normalizeCollinearSearchScope = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  /** @type {Record<string, string>} */
  const aliases = {
    adjacent_pairs: 'adjacent',
    adjacent_pair: 'adjacent',
    neighbor_pairs: 'adjacent',
    neighboring_pairs: 'adjacent',
    all_records: 'all',
    all_pairs: 'all',
    all_record_pairs: 'all',
    global: 'all'
  };
  const resolved = aliases[normalized] || normalized;
  return ['adjacent', 'all'].includes(resolved) ? resolved : 'adjacent';
};

/** @param {unknown} value */
export const collinearGroupScopeForEvidenceScope = (value) =>
  normalizeCollinearSearchScope(value) === 'all' ? 'global_collinear' : 'adjacent_local';

/** @param {unknown} value */
export const normalizeGroupMetadataScope = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  if (normalized === 'global_collinear' || normalized === 'all_records') return 'global_collinear';
  if (normalized === 'adjacent_local' || normalized === 'local_collinear') return 'adjacent_local';
  if (normalized === 'record_local') return 'record_local';
  return 'global';
};

/** @param {unknown} value */
export const groupMetadataScopeLabel = (value) => {
  const scope = normalizeGroupMetadataScope(value);
  if (scope === 'global_collinear') return 'Collinearity-backed global evidence';
  if (scope === 'adjacent_local') return 'Local collinear group';
  if (scope === 'record_local') return 'Record-specific similarity group';
  return 'Cross-record similarity group';
};
