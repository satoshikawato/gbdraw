export const normalizeCollinearAnchorMode = (_value) => 'rbh';

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
  if (scope === 'record_local') return 'Species-specific orthogroup';
  return 'Cross-record orthogroup';
};
