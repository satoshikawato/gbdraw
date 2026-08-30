export const normalizeCollinearAnchorMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
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

export const normalizeOrthogroupMembershipMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  const aliases = {
    legacy: 'anchor_core_v1',
    rbh: 'anchor_core_v1',
    rbh_only: 'anchor_core_v1',
    merge: 'anchor_core_v1',
    family: 'anchor_core_v1',
    family_merge: 'anchor_core_v1',
    local_split: 'anchor_core_v1',
    density_split: 'anchor_core_v1',
    outparalog_split: 'anchor_core_v1',
    distribution_split: 'anchor_core_v1',
    orthogroups: 'anchor_core_v1',
    anchor_core: 'anchor_core_v1'
  };
  const resolved = aliases[normalized] || normalized;
  return resolved === 'anchor_core_v1' ? resolved : 'anchor_core_v1';
};

export const normalizeCollinearSearchScope = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
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

export const collinearGroupScopeForEvidenceScope = (value) =>
  normalizeCollinearSearchScope(value) === 'all' ? 'global_collinear' : 'adjacent_local';

export const normalizeGroupMetadataScope = (value) => {
  const normalized = String(value || '').trim().toLowerCase().replace(/-/g, '_');
  if (normalized === 'global_collinear' || normalized === 'all_records') return 'global_collinear';
  if (normalized === 'adjacent_local' || normalized === 'local_collinear') return 'adjacent_local';
  if (normalized === 'record_local') return 'record_local';
  return 'global';
};

export const groupMetadataScopeLabel = (value) => {
  const scope = normalizeGroupMetadataScope(value);
  if (scope === 'global_collinear') return 'Collinearity-backed global evidence';
  if (scope === 'adjacent_local') return 'Local collinear group';
  if (scope === 'record_local') return 'Record-specific similarity group';
  return 'Cross-record similarity group';
};
