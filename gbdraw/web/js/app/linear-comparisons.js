export const LINEAR_COMPARISON_MODES = Object.freeze({
  NONE: 'none',
  ADJACENT: 'adjacent',
  SELECTED: 'selected'
});

export const LINEAR_COMPARISON_SOURCES = Object.freeze({
  LOSAT: 'losat',
  UPLOAD: 'upload'
});

const VALID_MODES = new Set(Object.values(LINEAR_COMPARISON_MODES));
const VALID_SOURCES = new Set(Object.values(LINEAR_COMPARISON_SOURCES));
const VALID_LOSAT_PROGRAMS = new Set(['blastn', 'tblastx', 'blastp']);

let comparisonCounter = 0;

const cleanUid = (value) => String(value || '').trim();

const normalizeMode = (value) => {
  const mode = String(value || '').trim().toLowerCase();
  return VALID_MODES.has(mode) ? mode : LINEAR_COMPARISON_MODES.ADJACENT;
};

const normalizeSource = (value, fallback = LINEAR_COMPARISON_SOURCES.LOSAT) => {
  const source = String(value || '').trim().toLowerCase();
  return VALID_SOURCES.has(source) ? source : fallback;
};

const createComparisonId = () => {
  comparisonCounter += 1;
  if (globalThis.crypto?.randomUUID) {
    return `linear-comparison-${globalThis.crypto.randomUUID()}`;
  }
  return `linear-comparison-${Date.now()}-${comparisonCounter}`;
};

export const linearComparisonEdgeKey = (queryUid, subjectUid) => (
  `${cleanUid(queryUid)}->${cleanUid(subjectUid)}`
);

export const createDefaultLinearComparisonPlan = () => ({
  mode: LINEAR_COMPARISON_MODES.ADJACENT,
  defaultSource: LINEAR_COMPARISON_SOURCES.LOSAT,
  edges: []
});

export const createLinearComparisonEdge = ({
  id = '',
  queryUid = '',
  subjectUid = '',
  included = true,
  fileActive = false,
  losatFilenameActive = false,
  source = LINEAR_COMPARISON_SOURCES.UPLOAD,
  file = null,
  losatFilename = ''
} = {}) => ({
  id: cleanUid(id) || createComparisonId(),
  queryUid: cleanUid(queryUid),
  subjectUid: cleanUid(subjectUid),
  included: included === true,
  fileActive: fileActive === true,
  losatFilenameActive: losatFilenameActive === true,
  source: normalizeSource(source, LINEAR_COMPARISON_SOURCES.UPLOAD),
  file: file || null,
  losatFilename: String(losatFilename || '')
});

const normalizeEdge = (edge, index) => {
  const source = edge && typeof edge === 'object' ? edge : {};
  return createLinearComparisonEdge({
    ...source,
    id: cleanUid(source.id) || `linear-comparison-draft-${index + 1}`,
    included: source.included === true,
    fileActive: source.fileActive === true,
    losatFilenameActive: source.losatFilenameActive === true
  });
};

export const normalizeLinearComparisonPlan = (plan = {}) => {
  const source = plan && typeof plan === 'object' && !Array.isArray(plan) ? plan : {};
  const edges = (Array.isArray(source.edges) ? source.edges : []).map(normalizeEdge);
  const usedIds = new Set();
  edges.forEach((edge, index) => {
    const baseId = cleanUid(edge.id) || `linear-comparison-draft-${index + 1}`;
    let id = baseId;
    let suffix = 2;
    while (usedIds.has(id)) {
      id = `${baseId}-${suffix}`;
      suffix += 1;
    }
    edge.id = id;
    usedIds.add(id);
  });
  return {
    mode: normalizeMode(source.mode),
    defaultSource: normalizeSource(source.defaultSource),
    edges
  };
};

export const reconcileLinearComparisonPlan = (plan, sequences = []) => {
  const normalized = normalizeLinearComparisonPlan(plan);
  const validUids = new Set(
    (Array.isArray(sequences) ? sequences : []).map((sequence) => cleanUid(sequence?.uid)).filter(Boolean)
  );
  normalized.edges = normalized.edges.filter((edge) => (
    validUids.has(edge.queryUid) && validUids.has(edge.subjectUid)
  ));
  return normalized;
};

const positiveRow = (value, fallback) => {
  const row = Number(value);
  return Number.isInteger(row) && row > 0 ? row : fallback;
};

export const normalizeLinearComparisonRows = (sequences = [], layout = []) => {
  const storedRows = new Map(
    (Array.isArray(layout) ? layout : []).map((entry) => [
      cleanUid(entry?.uid),
      Number(entry?.row)
    ])
  );
  return (Array.isArray(sequences) ? sequences : []).map((sequence, index) => ({
    uid: cleanUid(sequence?.uid),
    row: positiveRow(storedRows.get(cleanUid(sequence?.uid)), index + 1)
  }));
};

export const adjacentRowPairs = (sequences = [], layout = [], allPairs = false) => {
  const rows = new Map(normalizeLinearComparisonRows(sequences, layout).map((entry) => [entry.uid, entry.row]));
  const grouped = new Map();
  (Array.isArray(sequences) ? sequences : []).forEach((sequence, index) => {
    const uid = cleanUid(sequence?.uid);
    const row = rows.get(uid) || index + 1;
    if (!grouped.has(row)) grouped.set(row, []);
    grouped.get(row).push(uid);
  });
  const orderedRows = [...grouped.keys()].sort((left, right) => left - right);
  const pairs = [];
  for (let index = 0; index < orderedRows.length - 1; index += 1) {
    const upper = grouped.get(orderedRows[index]);
    const lower = grouped.get(orderedRows[index + 1]);
    if (allPairs) {
      upper.forEach((queryUid) => {
        lower.forEach((subjectUid) => pairs.push([queryUid, subjectUid]));
      });
      continue;
    }
    const pairCount = Math.min(upper.length, lower.length);
    for (let pairIndex = 0; pairIndex < pairCount; pairIndex += 1) {
      pairs.push([upper[pairIndex], lower[pairIndex]]);
    }
  }
  return pairs;
};

const validationIssue = (code, message, edge = null) => ({
  code,
  message,
  edgeId: cleanUid(edge?.id),
  edgeKey: edge ? linearComparisonEdgeKey(edge.queryUid, edge.subjectUid) : ''
});

export const validateLinearComparisonEdges = ({
  edges = [],
  sequences = [],
  layout = []
} = {}) => {
  const sequenceList = Array.isArray(sequences) ? sequences : [];
  const validUids = new Set(sequenceList.map((sequence) => cleanUid(sequence?.uid)).filter(Boolean));
  const rows = new Map(normalizeLinearComparisonRows(sequenceList, layout).map((entry) => [entry.uid, entry.row]));
  const seen = new Set();
  const issues = [];

  (Array.isArray(edges) ? edges : []).forEach((edge) => {
    const edgeKey = linearComparisonEdgeKey(edge?.queryUid, edge?.subjectUid);
    if (seen.has(edgeKey)) {
      issues.push(validationIssue('duplicate', 'Selected comparisons cannot contain the same directed pair more than once.', edge));
    }
    seen.add(edgeKey);

    const queryUid = cleanUid(edge?.queryUid);
    const subjectUid = cleanUid(edge?.subjectUid);
    if (!validUids.has(queryUid) || !validUids.has(subjectUid)) {
      issues.push(validationIssue('missing-uid', 'A selected comparison references a record that is no longer available.', edge));
      return;
    }
    if (queryUid === subjectUid) {
      issues.push(validationIssue('self', 'A comparison cannot connect a record to itself.', edge));
      return;
    }
    const queryRow = rows.get(queryUid);
    const subjectRow = rows.get(subjectUid);
    if (queryRow === subjectRow) {
      issues.push(validationIssue('same-row', 'Comparisons must connect records in different rows.', edge));
    } else if (Math.abs(queryRow - subjectRow) !== 1) {
      issues.push(validationIssue('non-adjacent', 'Comparisons must connect adjacent rows.', edge));
    }
    if (
      normalizeSource(edge?.source, LINEAR_COMPARISON_SOURCES.UPLOAD) === LINEAR_COMPARISON_SOURCES.UPLOAD &&
      (edge?.fileActive !== true || !edge?.file)
    ) {
      issues.push(validationIssue(
        'missing-upload',
        'Every selected uploaded comparison requires an active BLAST TSV file.',
        edge
      ));
    }
  });
  return issues;
};

const resolvedEdge = ({ edge, queryUid, subjectUid, queryIndex, subjectIndex, source, ordinal }) => Object.freeze({
  id: cleanUid(edge?.id) || `linear-derived-${linearComparisonEdgeKey(queryUid, subjectUid)}`,
  queryUid,
  subjectUid,
  queryIndex,
  subjectIndex,
  edgeKey: linearComparisonEdgeKey(queryUid, subjectUid),
  ordinal,
  source,
  included: true,
  file: edge?.file || null,
  fileActive: edge?.fileActive === true,
  losatFilename: String(edge?.losatFilename || ''),
  losatFilenameActive: edge?.losatFilenameActive === true
});

export const resolveLinearComparisonPlan = ({
  plan = createDefaultLinearComparisonPlan(),
  sequences = [],
  layout = [],
  losatProgram = 'blastn',
  blastpMode = 'orthogroup'
} = {}) => {
  const normalized = normalizeLinearComparisonPlan(plan);
  const sequenceList = Array.isArray(sequences) ? sequences : [];
  const indexByUid = new Map(sequenceList.map((sequence, index) => [cleanUid(sequence?.uid), index]));
  const draftByEdgeKey = new Map();
  normalized.edges.forEach((edge) => {
    const edgeKey = linearComparisonEdgeKey(edge.queryUid, edge.subjectUid);
    if (!draftByEdgeKey.has(edgeKey)) draftByEdgeKey.set(edgeKey, edge);
  });

  let activeDrafts = [];
  if (normalized.mode === LINEAR_COMPARISON_MODES.SELECTED) {
    activeDrafts = normalized.edges.filter((edge) => edge.included === true);
  } else if (normalized.mode === LINEAR_COMPARISON_MODES.ADJACENT) {
    activeDrafts = adjacentRowPairs(sequenceList, layout).map(([queryUid, subjectUid]) => {
      const edgeKey = linearComparisonEdgeKey(queryUid, subjectUid);
      const draft = draftByEdgeKey.get(edgeKey);
      return {
        ...(draft || {
          id: `linear-derived-${edgeKey}`,
          file: null,
          fileActive: false,
          losatFilename: '',
          losatFilenameActive: false
        }),
        queryUid,
        subjectUid,
        included: true,
        source: normalized.defaultSource
      };
    });
    if (normalized.defaultSource === LINEAR_COMPARISON_SOURCES.UPLOAD) {
      activeDrafts = activeDrafts.filter((edge) => edge.fileActive === true && Boolean(edge.file));
    }
  }

  const issues = normalized.mode === LINEAR_COMPARISON_MODES.SELECTED
    ? validateLinearComparisonEdges({ edges: activeDrafts, sequences: sequenceList, layout })
    : [];
  const normalizedProgram = VALID_LOSAT_PROGRAMS.has(String(losatProgram || '').trim().toLowerCase())
    ? String(losatProgram).trim().toLowerCase()
    : 'blastn';
  const normalizedBlastpMode = String(blastpMode || 'orthogroup').trim().toLowerCase();
  const selectedHasLosat = normalized.mode === LINEAR_COMPARISON_MODES.SELECTED &&
    activeDrafts.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT);
  if (selectedHasLosat && normalizedProgram === 'blastp' && normalizedBlastpMode !== 'pairwise') {
    issues.push(validationIssue(
      'selected-losat-requires-pairwise',
      'Selected LOSAT comparisons require LOSATP Pairwise mode. Similarity groups and Collinear blocks use All adjacent LOSAT.',
      null
    ));
  }

  const edges = activeDrafts.map((edge, ordinal) => resolvedEdge({
    edge,
    queryUid: cleanUid(edge.queryUid),
    subjectUid: cleanUid(edge.subjectUid),
    queryIndex: indexByUid.get(cleanUid(edge.queryUid)),
    subjectIndex: indexByUid.get(cleanUid(edge.subjectUid)),
    source: normalizeSource(edge.source),
    ordinal
  }));
  const frozenIssues = Object.freeze(issues.map((issue) => Object.freeze(issue)));
  const frozenEdges = Object.freeze(edges);
  const hasLosatIntent = edges.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT);
  const hasUploadIntent = edges.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD);

  return Object.freeze({
    mode: normalized.mode,
    defaultSource: normalized.defaultSource,
    edges: frozenEdges,
    errors: frozenIssues,
    error: frozenIssues[0]?.message || '',
    valid: frozenIssues.length === 0,
    hasComparisonIntent: edges.length > 0,
    hasLosatIntent,
    hasUploadIntent
  });
};

export const hasLinearComparisonIntent = (resolution) => Boolean(resolution?.hasComparisonIntent);
export const hasActiveLinearLosatIntent = (resolution) => Boolean(resolution?.hasLosatIntent);
export const hasActiveLinearUploadIntent = (resolution) => Boolean(resolution?.hasUploadIntent);

const hasRetainedPayload = (edge) => Boolean(
  edge?.file || String(edge?.losatFilename || '').trim()
);

export const materializeResolvedEdgesAsSelectedPlan = (plan, resolution) => {
  const normalized = normalizeLinearComparisonPlan(plan);
  const existingByKey = new Map(
    normalized.edges.map((edge) => [linearComparisonEdgeKey(edge.queryUid, edge.subjectUid), edge])
  );
  const activeKeys = new Set();
  const selectedEdges = (Array.isArray(resolution?.edges) ? resolution.edges : []).map((edge) => {
    activeKeys.add(edge.edgeKey);
    const existing = existingByKey.get(edge.edgeKey);
    return createLinearComparisonEdge({
      ...(existing || edge),
      id: existing?.id || '',
      queryUid: edge.queryUid,
      subjectUid: edge.subjectUid,
      included: true,
      source: edge.source
    });
  });
  normalized.edges.forEach((edge) => {
    const edgeKey = linearComparisonEdgeKey(edge.queryUid, edge.subjectUid);
    if (activeKeys.has(edgeKey) || !hasRetainedPayload(edge)) return;
    selectedEdges.push({ ...edge, included: false });
  });
  return {
    mode: LINEAR_COMPARISON_MODES.SELECTED,
    defaultSource: normalized.defaultSource,
    edges: selectedEdges
  };
};

export const buildPairwiseLosatJobSpecs = ({
  resolution,
  program = 'blastn',
  blastpMode = 'pairwise'
} = {}) => {
  if (!resolution || !Array.isArray(resolution.edges)) {
    throw new Error('A resolved Linear comparison plan is required.');
  }
  if (resolution.error) throw new Error(resolution.error);
  const normalizedProgram = String(program || '').trim().toLowerCase();
  if (!VALID_LOSAT_PROGRAMS.has(normalizedProgram)) {
    throw new Error(`Unsupported LOSAT program: ${program}`);
  }
  if (normalizedProgram === 'blastp' && String(blastpMode || '').trim().toLowerCase() !== 'pairwise') {
    if (resolution.mode === LINEAR_COMPARISON_MODES.SELECTED && resolution.hasLosatIntent) {
      throw new Error('Selected LOSAT comparisons require LOSATP Pairwise mode.');
    }
    return Object.freeze([]);
  }
  return Object.freeze(resolution.edges
    .filter((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT)
    .map((edge) => Object.freeze({
      edgeKey: edge.edgeKey,
      ordinal: edge.ordinal,
      queryUid: edge.queryUid,
      subjectUid: edge.subjectUid,
      queryIndex: edge.queryIndex,
      subjectIndex: edge.subjectIndex,
      program: normalizedProgram
    })));
};
