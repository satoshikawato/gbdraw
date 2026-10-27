// @ts-check
import { drawnLinearSequences, omittedLinearUids } from './record-draw-selection.js';
import { msg } from '../utils/i18n.js';
/** @import { FieldLeaf } from '../utils/field-reading.js' */
/** @import { InputFileRef } from './history-files.js' */
/** @typedef {(leaf: FieldLeaf) => void} PlanLossCollector Told each stored value the plan does not use as stored. */
/**
 * @typedef {'none' | 'adjacent' | 'selected'} LinearComparisonMode
 * @typedef {'losat' | 'upload'} LinearComparisonSource
 */
/**
 * One draft pair of the Linear comparison plan. `file` names the uploaded
 * BLAST TSV by ref (`state.inputFiles`), or is null. A Load candidate's edge
 * holds its loaded file until `bindLoadedInputs` registers it at install.
 * @typedef {object} LinearComparisonEdge
 * @property {string} id
 * @property {string} queryUid
 * @property {string} subjectUid
 * @property {boolean} included
 * @property {boolean} fileActive
 * @property {boolean} losatFilenameActive
 * @property {LinearComparisonSource} source
 * @property {InputFileRef | null} file
 * @property {string} losatFilename
 */
/**
 * The options of `createLinearComparisonEdge`; every field may be absent or
 * unnormalized.
 * @typedef {object} LinearComparisonEdgeInit
 * @property {unknown} [id]
 * @property {unknown} [queryUid]
 * @property {unknown} [subjectUid]
 * @property {unknown} [included]
 * @property {unknown} [fileActive]
 * @property {unknown} [losatFilenameActive]
 * @property {unknown} [source]
 * @property {unknown} [file]
 * @property {unknown} [losatFilename]
 */
/**
 * @typedef {object} LinearComparisonPlan
 * @property {LinearComparisonMode} mode
 * @property {LinearComparisonSource} defaultSource
 * @property {LinearComparisonEdge[]} edges
 */
/**
 * One edge `resolveLinearComparisonPlan` derives from the plan.
 * @typedef {object} LinearComparisonResolvedEdge
 * @property {string} id
 * @property {string} queryUid
 * @property {string} subjectUid
 * @property {number | undefined} queryIndex `undefined` when no record card has the uid (an issue).
 * @property {number | undefined} subjectIndex `undefined` when no record card has the uid (an issue).
 * @property {string} edgeKey
 * @property {number} ordinal
 * @property {LinearComparisonSource} source
 * @property {boolean} included
 * @property {InputFileRef | null} file
 * @property {boolean} fileActive
 * @property {string} losatFilename
 * @property {boolean} losatFilenameActive
 */
/**
 * @typedef {object} LinearComparisonIssue
 * @property {string} code
 * @property {string} message
 * @property {string} edgeId
 * @property {string} edgeKey
 */
/**
 * @typedef {object} LinearComparisonResolution
 * @property {LinearComparisonMode} mode
 * @property {LinearComparisonSource} defaultSource
 * @property {readonly LinearComparisonResolvedEdge[]} edges
 * @property {readonly LinearComparisonIssue[]} errors
 * @property {string} error
 * @property {boolean} valid
 * @property {boolean} hasComparisonIntent
 * @property {boolean} hasLosatIntent
 * @property {boolean} hasUploadIntent
 */
/**
 * A record of a drawn Linear row.
 * @template S
 * @typedef {{ uid: string, index: number, sequence: S }} LinearRowRecord
 */
/**
 * A record card as the comparison reads it; the card's other fields are kept.
 * @typedef {{ uid?: unknown }} LinearComparisonSequence
 */
/**
 * A saved record row of the drawing.
 * @typedef {{ uid?: unknown, row?: unknown }} LinearComparisonLayoutRow
 */
/**
 * One record-pair search a Generate requests.
 * @typedef {object} LosatJobSpec
 * @property {string} edgeKey
 * @property {number} ordinal
 * @property {string} queryUid
 * @property {string} subjectUid
 * @property {number} queryIndex
 * @property {number} subjectIndex
 * @property {string} program
 */
export const LINEAR_COMPARISON_MODES = Object.freeze({
  NONE: 'none',
  ADJACENT: 'adjacent',
  SELECTED: 'selected'
});

export const LINEAR_COMPARISON_SOURCES = Object.freeze({
  LOSAT: 'losat',
  UPLOAD: 'upload'
});

/**
 * The fields of a Session 27-44 draft that held its comparisons before the
 * comparison plan: the migration to the plan reads them and removes them
 * (`migrateLegacyLinearComparisonDraft`, Python's `session_io`), so no
 * settings reader reads them (dotted paths in the draft `config`).
 * @type {readonly string[]}
 */
export const LEGACY_COMPARISON_SOURCE_FIELDS = Object.freeze(['blastSource', 'adv.blastSource', 'linearRecordLayout.comparisons']);

/** @type {Set<string>} */
const VALID_MODES = new Set(Object.values(LINEAR_COMPARISON_MODES));
/** @type {Set<string>} */
const VALID_SOURCES = new Set(Object.values(LINEAR_COMPARISON_SOURCES));
const VALID_LOSAT_PROGRAMS = new Set(['blastn', 'tblastx', 'blastp']);

let comparisonCounter = 0;

/** @param {unknown} value */
const cleanUid = (value) => String(value || '').trim();

/** @param {unknown} value */
export const plainTextLinearRecordLabel = (value) => {
  let label = String(value ?? '');
  /** @type {string | undefined} */
  let previous;
  do {
    previous = label;
    label = label.replace(/<[^<>]*>/g, '');
  } while (label !== previous);
  return label.replace(/[<>]/g, '').trim() || 'Record';
};

// An absent or unknown mode is the Web default (No comparison): a plan never
// starts a comparison its writer did not state (B15).
/**
 * @param {unknown} value
 * @returns {LinearComparisonMode}
 */
const normalizeMode = (value) => {
  const mode = String(value || '').trim().toLowerCase();
  return VALID_MODES.has(mode) ? /** @type {LinearComparisonMode} */ (mode) : LINEAR_COMPARISON_MODES.NONE;
};

/**
 * @param {unknown} value
 * @param {LinearComparisonSource} [fallback]
 * @returns {LinearComparisonSource}
 */
const normalizeSource = (value, fallback = LINEAR_COMPARISON_SOURCES.LOSAT) => {
  const source = String(value || '').trim().toLowerCase();
  return VALID_SOURCES.has(source) ? /** @type {LinearComparisonSource} */ (source) : fallback;
};

const createComparisonId = () => {
  comparisonCounter += 1;
  if (globalThis.crypto?.randomUUID) {
    return `linear-comparison-${globalThis.crypto.randomUUID()}`;
  }
  return `linear-comparison-${Date.now()}-${comparisonCounter}`;
};

/**
 * @param {unknown} queryUid
 * @param {unknown} subjectUid
 */
export const linearComparisonEdgeKey = (queryUid, subjectUid) => (
  `${cleanUid(queryUid)}->${cleanUid(subjectUid)}`
);

/** @returns {LinearComparisonPlan} */
export const createDefaultLinearComparisonPlan = () => ({
  mode: LINEAR_COMPARISON_MODES.NONE,
  defaultSource: LINEAR_COMPARISON_SOURCES.LOSAT,
  edges: []
});

/**
 * @param {LinearComparisonEdgeInit} [init]
 * @returns {LinearComparisonEdge}
 */
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
  // A binding's ref, or the loaded file of a Load candidate (see the typedef).
  file: /** @type {InputFileRef | null} */ (file || null),
  losatFilename: String(losatFilename || '')
});

/**
 * An edge's uploaded BLAST TSV, the one place that sets it: a file is the
 * edge's active upload; none leaves the edge with no file and none active.
 * The pair card's file input and its Remove set it, and so does Load where it
 * drops the edge's file binding (OV-899).
 * @param {LinearComparisonEdge} edge Changed in place.
 * @param {InputFileRef | null} file The ref of the file (`state.inputFiles.registerFile`).
 */
export const setComparisonEdgeFile = (edge, file) => {
  edge.file = file || null;
  edge.fileActive = Boolean(file);
};

/** @param {unknown} value */
const isStored = (value) => value !== undefined && value !== null;
/**
 * @param {Readonly<Record<string, unknown>>} source
 * @param {number} index
 * @param {PlanLossCollector} collect
 */
const normalizeEdge = (source, index, collect) => {
  const normalized = createLinearComparisonEdge({
    ...source,
    id: cleanUid(source.id) || `linear-comparison-draft-${index + 1}`,
    included: source.included === true,
    fileActive: source.fileActive === true,
    losatFilenameActive: source.losatFilenameActive === true
  });
  /** @type {const} */ (['included', 'fileActive', 'losatFilenameActive', 'source', 'losatFilename']).forEach((field) => {
    if (isStored(source[field]) && source[field] !== normalized[field]) {
      collect({ at: `.edges[${index}].${field}`, stored: source[field], used: normalized[field] });
    }
  });
  // A field an edge does not have is not read (E-OD03, E-OD12).
  Object.keys(source).filter((field) => !Object.hasOwn(normalized, field)).forEach((field) => {
    collect({ at: `.edges[${index}].${field}`, stored: source[field], used: undefined });
  });
  return normalized;
};

/**
 * @param {unknown} [plan] unvalidated plan (a session, a draft, or a live plan)
 * @param {PlanLossCollector} [collect] Load lists what it is told (OV-764).
 * @returns {LinearComparisonPlan}
 */
export const normalizeLinearComparisonPlan = (plan = {}, collect = () => {}) => {
  const source = /** @type {{ edges?: unknown, mode?: unknown, defaultSource?: unknown }} */ (
    plan && typeof plan === 'object' && !Array.isArray(plan) ? plan : {}
  );
  if (isStored(source.edges) && !Array.isArray(source.edges)) collect({ at: '.edges', stored: source.edges, used: [] });
  // An edge that is no object pairs no records: it is dropped, and the
  // generated ids keep the stored positions.
  const edges = (Array.isArray(source.edges) ? source.edges : []).flatMap((edge, index) => {
    if (edge && typeof edge === 'object' && !Array.isArray(edge)) return [normalizeEdge(edge, index, collect)];
    collect({ at: `.edges[${index}]`, stored: edge, used: undefined });
    return [];
  });
  /** @type {Set<string>} */
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
  const normalized = {
    mode: normalizeMode(source.mode),
    defaultSource: normalizeSource(source.defaultSource),
    edges
  };
  /** @type {const} */ (['mode', 'defaultSource']).forEach((field) => {
    if (isStored(source[field]) && source[field] !== normalized[field]) collect({ at: `.${field}`, stored: source[field], used: normalized[field] });
  });
  return normalized;
};

/**
 * @param {unknown} plan
 * @param {readonly unknown[]} [_sequences]
 */
export const reconcileLinearComparisonPlan = (plan, _sequences = []) => (
  normalizeLinearComparisonPlan(plan)
);

/**
 * @param {unknown} value
 * @param {number} fallback
 */
const positiveRow = (value, fallback) => {
  const row = Number(value);
  return Number.isInteger(row) && row > 0 ? row : fallback;
};

/**
 * @param {readonly LinearComparisonSequence[] | null} [sequences]
 * @param {readonly LinearComparisonLayoutRow[] | null} [layout]
 */
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

/**
 * @template {LinearComparisonSequence} S
 * @param {readonly S[] | null} [sequences]
 * @param {readonly LinearComparisonLayoutRow[] | null} [layout]
 */
const linearComparisonRowTopology = (sequences = [], layout = []) => {
  const sequenceList = Array.isArray(sequences) ? sequences : [];
  const normalizedRows = normalizeLinearComparisonRows(sequenceList, layout);
  const rowByUid = new Map(normalizedRows.map((entry) => [entry.uid, entry.row]));
  /** @type {Map<number, LinearRowRecord<S>[]>} */
  const recordsByRow = new Map();
  sequenceList.forEach((sequence, index) => {
    const uid = cleanUid(sequence?.uid);
    const row = rowByUid.get(uid) || index + 1;
    if (!recordsByRow.has(row)) recordsByRow.set(row, []);
    // The line above makes sure the row has its list.
    /** @type {LinearRowRecord<S>[]} */ (recordsByRow.get(row)).push({ uid, index, sequence });
  });
  const orderedRows = [...recordsByRow.keys()].sort((left, right) => left - right);
  return {
    rowByUid,
    recordsByRow,
    orderedRows,
    rowPositionByNumber: new Map(orderedRows.map((row, index) => [row, index]))
  };
};

/**
 * @param {readonly LinearComparisonSequence[] | null} [sequences]
 * @param {readonly LinearComparisonLayoutRow[] | null} [layout]
 * @param {boolean} [allPairs]
 */
export const adjacentRowPairs = (sequences = [], layout = [], allPairs = false) => {
  const { recordsByRow, orderedRows } = linearComparisonRowTopology(sequences, layout);
  /** @type {[string, string][]} */
  const pairs = [];
  for (let index = 0; index < orderedRows.length - 1; index += 1) {
    // Every ordered row is a key of `recordsByRow` (linearComparisonRowTopology).
    const upper = /** @type {LinearRowRecord<LinearComparisonSequence>[]} */ (recordsByRow.get(orderedRows[index]));
    const lower = /** @type {LinearRowRecord<LinearComparisonSequence>[]} */ (recordsByRow.get(orderedRows[index + 1]));
    if (allPairs) {
      upper.forEach(({ uid: queryUid }) => {
        lower.forEach(({ uid: subjectUid }) => pairs.push([queryUid, subjectUid]));
      });
      continue;
    }
    const pairCount = Math.min(upper.length, lower.length);
    for (let pairIndex = 0; pairIndex < pairCount; pairIndex += 1) {
      pairs.push([upper[pairIndex].uid, lower[pairIndex].uid]);
    }
  }
  return pairs;
};

/**
 * The display pairs of a collinear comparison laid out by a CLI Session's
 * complete row grid: every pair across adjacent occupied rows of the drawn
 * records, as Python closes up a row that records turned OFF leave empty.
 * Null when the drawn records do not all hold a canonical grid cell; the
 * comparison plan's edges decide then.
 * @param {readonly { uid: string }[]} sequences The drawn records, in request order.
 * @param {readonly { uid: string, row: number, canonicalRow?: number, canonicalColumn?: number,
 *   canonicalCardinality?: string }[]} rows The drawing's record rows.
 * @returns {((queryIndex: number, subjectIndex: number) => boolean) | null}
 */
export const canonicalGridDisplayPairs = (sequences, rows) => {
  const cells = sequences.map((seq) => rows.find((entry) => entry.uid === seq.uid));
  if (cells.length < 2 || !cells.every((entry) => entry?.canonicalRow === entry?.row
    && Number.isInteger(entry?.canonicalColumn) && Number(entry?.canonicalColumn) > 0
    && entry?.canonicalCardinality === 'exactly_one')) return null;
  const pairs = new Set(adjacentRowPairs([...sequences], /** @type {{ uid: string, row: number }[]} */ (cells), true)
    .map(([upper, lower]) => `${upper}\0${lower}`));
  return (queryIndex, subjectIndex) => {
    const query = sequences[queryIndex]?.uid;
    const subject = sequences[subjectIndex]?.uid;
    return queryIndex < subjectIndex && (pairs.has(`${query}\0${subject}`) || pairs.has(`${subject}\0${query}`));
  };
};

/**
 * @param {string} code
 * @param {string} message
 * @param {LinearComparisonEdgeInit | null} [edge]
 * @returns {LinearComparisonIssue}
 */
const validationIssue = (code, message, edge = null) => ({
  code,
  message,
  edgeId: cleanUid(edge?.id),
  edgeKey: edge ? linearComparisonEdgeKey(edge.queryUid, edge.subjectUid) : ''
});

/**
 * @param {{
 *   edges?: readonly (LinearComparisonEdgeInit | null | undefined)[] | null,
 *   sequences?: readonly LinearComparisonSequence[] | null,
 *   layout?: readonly LinearComparisonLayoutRow[] | null
 * }} [options]
 */
export const validateLinearComparisonEdges = ({
  edges = [],
  sequences = [],
  layout = []
} = {}) => {
  const sequenceList = Array.isArray(sequences) ? sequences : [];
  const validUids = new Set(sequenceList.map((sequence) => cleanUid(sequence?.uid)).filter(Boolean));
  const { rowByUid, rowPositionByNumber } = linearComparisonRowTopology(sequenceList, layout);
  /** @type {Set<string>} */
  const seen = new Set();
  /** @type {LinearComparisonIssue[]} */
  const issues = [];

  (Array.isArray(edges) ? edges : []).forEach((edge) => {
    const edgeKey = linearComparisonEdgeKey(edge?.queryUid, edge?.subjectUid);
    if (seen.has(edgeKey)) {
      issues.push(validationIssue('duplicate', msg('Selected comparisons cannot contain the same directed pair more than once.'), edge));
    }
    seen.add(edgeKey);

    const queryUid = cleanUid(edge?.queryUid);
    const subjectUid = cleanUid(edge?.subjectUid);
    if (!validUids.has(queryUid) || !validUids.has(subjectUid)) {
      issues.push(validationIssue('missing-uid', msg('A selected comparison references a record that is no longer available.'), edge));
      return;
    }
    if (queryUid === subjectUid) {
      issues.push(validationIssue('self', msg('A comparison cannot connect a record to itself.'), edge));
      return;
    }
    // Both uids are valid here, and every valid uid has a row (linearComparisonRowTopology).
    const queryRow = /** @type {number} */ (rowByUid.get(queryUid));
    const subjectRow = /** @type {number} */ (rowByUid.get(subjectUid));
    if (queryRow === subjectRow) {
      issues.push(validationIssue('same-row', msg('Comparisons must connect records in different rows.'), edge));
    } else if (
      // Every valid uid has a row and every row has a position (linearComparisonRowTopology).
      Math.abs(/** @type {number} */ (rowPositionByNumber.get(queryRow)) - /** @type {number} */ (rowPositionByNumber.get(subjectRow))) !== 1
    ) {
      issues.push(validationIssue('non-adjacent', msg('Comparisons must connect adjacent rows.'), edge));
    }
    if (
      normalizeSource(edge?.source, LINEAR_COMPARISON_SOURCES.UPLOAD) === LINEAR_COMPARISON_SOURCES.UPLOAD &&
      (edge?.fileActive !== true || !edge?.file)
    ) {
      issues.push(validationIssue(
        'missing-upload',
        msg('Every selected uploaded comparison requires an active BLAST TSV file.'),
        edge
      ));
    }
  });
  return issues;
};

/**
 * @param {{
 *   edge: (Omit<LinearComparisonEdgeInit, 'file'> & { file?: InputFileRef | null }) | null | undefined,
 *   queryUid: string, subjectUid: string,
 *   queryIndex: number | undefined, subjectIndex: number | undefined, source: LinearComparisonSource, ordinal: number
 * }} input
 */
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

// A selected pair that names an OFF record is left out with its record and
// comes back with it; a pair that names a record without a card stays an issue.
/**
 * @template {{ uid?: unknown }} T
 * @param {unknown} plan
 * @param {readonly T[]} sequences Every record card.
 * @param {Iterable<unknown> | null | undefined} recordsOff
 */
const drawnComparisonInputs = (plan, sequences, recordsOff) => {
  const normalized = normalizeLinearComparisonPlan(plan);
  const omitted = omittedLinearUids(sequences, recordsOff);
  if (omitted.size > 0) {
    normalized.edges = normalized.edges.filter((edge) => (
      !omitted.has(cleanUid(edge.queryUid)) && !omitted.has(cleanUid(edge.subjectUid))
    ));
  }
  return { normalized, sequenceList: drawnLinearSequences(sequences, recordsOff) };
};

/**
 * The plan a Generate runs: its pairs join the drawn records, and every index
 * is a position in the drawn list (`drawnLinearSequences`).
 * @param {{
 *   plan?: unknown,
 *   sequences?: readonly LinearComparisonSequence[] | null,
 *   recordsOff?: Iterable<unknown>,
 *   layout?: readonly LinearComparisonLayoutRow[] | null,
 *   losatProgram?: string,
 *   blastpMode?: string
 * }} [options]
 * @returns {LinearComparisonResolution}
 */
export const resolveLinearComparisonPlan = ({
  plan = createDefaultLinearComparisonPlan(),
  sequences = [],
  recordsOff = [],
  layout = [],
  losatProgram = 'blastn',
  blastpMode = 'orthogroup'
} = {}) => {
  const { normalized, sequenceList } = drawnComparisonInputs(
    plan, Array.isArray(sequences) ? sequences : [], recordsOff
  );
  const indexByUid = new Map(sequenceList.map((sequence, index) => [cleanUid(sequence?.uid), index]));
  /** @type {Map<string, LinearComparisonEdge>} */
  const draftByEdgeKey = new Map();
  normalized.edges.forEach((edge) => {
    const edgeKey = linearComparisonEdgeKey(edge.queryUid, edge.subjectUid);
    if (!draftByEdgeKey.has(edgeKey)) draftByEdgeKey.set(edgeKey, edge);
  });
  const normalizedProgram = VALID_LOSAT_PROGRAMS.has(String(losatProgram || '').trim().toLowerCase())
    ? String(losatProgram).trim().toLowerCase()
    : 'blastn';
  const normalizedBlastpMode = String(blastpMode || 'orthogroup').trim().toLowerCase();

  /** @type {LinearComparisonEdge[]} */
  let activeDrafts = [];
  if (normalized.mode === LINEAR_COMPARISON_MODES.SELECTED) {
    activeDrafts = normalized.edges.filter((edge) => edge.included === true);
  } else if (normalized.mode === LINEAR_COMPARISON_MODES.ADJACENT) {
    activeDrafts = adjacentRowPairs(
      sequenceList,
      layout,
      true
    ).map(([queryUid, subjectUid]) => {
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
  const selectedHasLosat = normalized.mode === LINEAR_COMPARISON_MODES.SELECTED &&
    activeDrafts.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT);
  if (selectedHasLosat && normalizedProgram === 'blastp' && normalizedBlastpMode !== 'pairwise') {
    issues.push(validationIssue(
      'selected-losat-requires-pairwise',
      msg('Selected LOSAT comparisons require LOSATP Pairwise mode. For Similarity groups or Collinear blocks, choose Run LOSAT for all adjacent pairs.'),
      null
    ));
  }

  const edges = activeDrafts.map((edge, ordinal) => resolvedEdge({
    edge,
    queryUid: cleanUid(edge.queryUid),
    subjectUid: cleanUid(edge.subjectUid),
    // `undefined` for a uid without a card (an issue); consumers of an invalid resolution guard it.
    queryIndex: indexByUid.get(cleanUid(edge.queryUid)),
    subjectIndex: indexByUid.get(cleanUid(edge.subjectUid)),
    source: normalizeSource(edge.source),
    ordinal
  }));
  // A missing upload comes last (D-54): Generate asks about it only when nothing
  // else blocks the run; otherwise the run and the panel report the other issue.
  const frozenIssues = Object.freeze([
    ...issues.filter((issue) => issue.code !== 'missing-upload'),
    ...issues.filter((issue) => issue.code === 'missing-upload')
  ].map((issue) => Object.freeze(issue)));
  const frozenEdges = Object.freeze(edges);
  const hasAnalysisIntent = sequenceList.length > 1
    && normalized.mode === LINEAR_COMPARISON_MODES.ADJACENT
    && normalized.defaultSource === LINEAR_COMPARISON_SOURCES.LOSAT
    && normalizedProgram === 'blastp'
    && ['orthogroup', 'collinear'].includes(normalizedBlastpMode);
  const hasLosatIntent = hasAnalysisIntent
    || edges.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT);
  const hasUploadIntent = edges.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD);

  return Object.freeze({
    mode: normalized.mode,
    defaultSource: normalized.defaultSource,
    edges: frozenEdges,
    errors: frozenIssues,
    error: frozenIssues[0]?.message || '',
    valid: frozenIssues.length === 0,
    hasComparisonIntent: hasAnalysisIntent || edges.length > 0,
    hasLosatIntent,
    hasUploadIntent
  });
};

/**
 * The pairs set to Upload BLAST TSV without an active file, when they are the
 * draft's only issues: Generate asks to choose a file or set them to No
 * comparison (D-54). Any other issue keeps its typed error, since setting the
 * pairs to No comparison would not let the run succeed.
 * @param {LinearComparisonResolution | null | undefined} resolution
 * @returns {{ edgeKey: string, queryIndex: number, subjectIndex: number, label: string }[]}
 */
export const missingUploadPairsToResolve = (resolution) => {
  const issues = resolution?.errors || [];
  if (issues.some((issue) => issue.code !== 'missing-upload')) return [];
  const missing = new Set(issues.map((issue) => issue.edgeKey));
  return (resolution?.edges || []).filter((edge) => missing.has(edge.edgeKey)).map((edge) => ({
    edgeKey: edge.edgeKey,
    // Every issue is `missing-upload` here, so both uids have a card and both indexes are numbers.
    queryIndex: /** @type {number} */ (edge.queryIndex),
    subjectIndex: /** @type {number} */ (edge.subjectIndex),
    label: `#${/** @type {number} */ (edge.queryIndex) + 1} → #${/** @type {number} */ (edge.subjectIndex) + 1}`
  }));
};

/** @param {LinearComparisonResolution | null | undefined} resolution */
export const hasLinearComparisonIntent = (resolution) => Boolean(resolution?.hasComparisonIntent);
/** @param {LinearComparisonResolution | null | undefined} resolution */
export const hasActiveLinearLosatIntent = (resolution) => Boolean(resolution?.hasLosatIntent);
/** @param {LinearComparisonResolution | null | undefined} resolution */
export const hasActiveLinearUploadIntent = (resolution) => Boolean(resolution?.hasUploadIntent);

/** @param {{ file?: unknown, losatFilename?: unknown } | null | undefined} edge */
const hasRetainedPayload = (edge) => Boolean(
  edge?.file || String(edge?.losatFilename || '').trim()
);

/**
 * @param {unknown} plan
 * @param {{ edges?: readonly LinearComparisonResolvedEdge[] } | null | undefined} resolution
 * @returns {LinearComparisonPlan}
 */
export const materializeResolvedEdgesAsSelectedPlan = (plan, resolution) => {
  const normalized = normalizeLinearComparisonPlan(plan);
  const existingById = new Map(normalized.edges.map((edge) => [edge.id, edge]));
  /** @type {Map<string, LinearComparisonEdge>} */
  const existingByKey = new Map();
  normalized.edges.forEach((edge) => {
    const edgeKey = linearComparisonEdgeKey(edge.queryUid, edge.subjectUid);
    if (!existingByKey.has(edgeKey)) existingByKey.set(edgeKey, edge);
  });
  /** @type {Set<string>} */
  const activeKeys = new Set();
  /** @type {Set<string>} */
  const activeIds = new Set();
  const selectedEdges = (Array.isArray(resolution?.edges) ? resolution.edges : []).map((edge) => {
    activeKeys.add(edge.edgeKey);
    const existing = existingById.get(edge.id) || existingByKey.get(edge.edgeKey);
    const selected = createLinearComparisonEdge({
      ...(existing || edge),
      id: existing?.id || '',
      queryUid: edge.queryUid,
      subjectUid: edge.subjectUid,
      included: true,
      source: edge.source
    });
    activeIds.add(selected.id);
    return selected;
  });
  normalized.edges.forEach((edge) => {
    const edgeKey = linearComparisonEdgeKey(edge.queryUid, edge.subjectUid);
    if (activeIds.has(edge.id)) return;
    if (!activeKeys.has(edgeKey) && !hasRetainedPayload(edge)) return;
    selectedEdges.push({ ...edge, included: false });
  });
  return {
    mode: LINEAR_COMPARISON_MODES.SELECTED,
    defaultSource: normalized.defaultSource,
    edges: selectedEdges
  };
};

const STRUCTURAL_PRESENTATION_ISSUES = new Set([
  'duplicate',
  'missing-uid',
  'self',
  'same-row',
  'non-adjacent'
]);

/**
 * @param {number | undefined} upperRow
 * @param {number | undefined} lowerRow
 */
const boundaryKey = (upperRow, lowerRow) => `${upperRow}->${lowerRow}`;

/**
 * A uid without a row looks up `undefined` in `rowPositionByNumber`, which has no such key.
 * @param {{
 *   queryUid: string, subjectUid: string, rowByUid: ReadonlyMap<string, number>,
 *   rowPositionByNumber: ReadonlyMap<number | undefined, number>,
 *   boundariesByKey: ReadonlyMap<string, LinearComparisonBoundary>
 * }} input
 */
const pairBoundary = ({ queryUid, subjectUid, rowByUid, rowPositionByNumber, boundariesByKey }) => {
  const queryRow = rowByUid.get(queryUid);
  const subjectRow = rowByUid.get(subjectUid);
  const queryPosition = rowPositionByNumber.get(queryRow);
  const subjectPosition = rowPositionByNumber.get(subjectRow);
  if (
    queryPosition === undefined ||
    subjectPosition === undefined ||
    Math.abs(queryPosition - subjectPosition) !== 1
  ) {
    return null;
  }
  const upperRow = queryPosition < subjectPosition ? queryRow : subjectRow;
  const lowerRow = queryPosition < subjectPosition ? subjectRow : queryRow;
  return boundariesByKey.get(boundaryKey(upperRow, lowerRow)) || null;
};

/**
 * @param {{
 *   normalizedPlan: LinearComparisonPlan, candidateKind: LinearComparisonCandidateKind,
 *   draft: LinearComparisonEdge | null, resolved: LinearComparisonResolvedEdge | null
 * }} input
 * @returns {LinearComparisonSource | 'none'}
 */
const presentationSource = ({ normalizedPlan, candidateKind, draft, resolved }) => {
  if (resolved) return resolved.source;
  if (
    candidateKind === 'zipped' &&
    normalizedPlan.mode === LINEAR_COMPARISON_MODES.ADJACENT
  ) {
    return normalizedPlan.defaultSource;
  }
  if (normalizedPlan.mode === LINEAR_COMPARISON_MODES.SELECTED && draft?.included) {
    return draft.source;
  }
  return 'none';
};

/** @typedef {'zipped' | 'resolved' | 'draft'} LinearComparisonCandidateKind Where a pair came from. */
/**
 * One pair of the timeline between two adjacent rows; `queryIndex` and `subjectIndex`
 * are positions in the drawn list.
 * @typedef {object} LinearComparisonPair
 * @property {string} edgeKey
 * @property {string} edgeId
 * @property {string} queryUid
 * @property {string} subjectUid
 * @property {number | undefined} queryIndex
 * @property {number | undefined} subjectIndex
 * @property {LinearComparisonSource | 'none'} source
 * @property {boolean} active
 * @property {LinearComparisonCandidateKind} candidateKind
 * @property {LinearComparisonEdge | null} draft
 * @property {LinearComparisonResolvedEdge | null} resolved
 */
/**
 * @typedef {{ upperRow: number, lowerRow: number, pairs: LinearComparisonPair[] }} LinearComparisonBoundary The pairs drawn between two adjacent rows.
 */

/**
 * The rows of drawn records and the pairs between them; a pair that names an
 * OFF record is neither drawn nor listed as unplaced.
 * @template {LinearComparisonSequence} S
 * @param {{ sequences?: readonly S[] | null, recordsOff?: Iterable<unknown>,
 *   layout?: readonly LinearComparisonLayoutRow[] | null, plan?: unknown,
 *   resolution?: LinearComparisonResolution | null }} [options]
 */
export const buildLinearComparisonTimeline = ({
  sequences = [],
  recordsOff = [],
  layout = [],
  plan = createDefaultLinearComparisonPlan(),
  resolution = null
} = {}) => {
  const { normalized: normalizedPlan, sequenceList } = drawnComparisonInputs(
    plan, Array.isArray(sequences) ? sequences : [], recordsOff
  );
  const indexByUid = new Map(sequenceList.map((sequence, index) => [cleanUid(sequence?.uid), index]));
  const {
    rowByUid,
    recordsByRow,
    orderedRows,
    rowPositionByNumber
  } = linearComparisonRowTopology(sequenceList, layout);
  /** @type {Map<string, LinearComparisonBoundary>} */
  const boundariesByKey = new Map();
  const rows = orderedRows.map((row, index) => {
    /** @type {LinearComparisonBoundary | null} */
    let boundaryAfter = null;
    if (index < orderedRows.length - 1) {
      const lowerRow = orderedRows[index + 1];
      boundaryAfter = {
        upperRow: row,
        lowerRow,
        pairs: []
      };
      boundariesByKey.set(boundaryKey(row, lowerRow), boundaryAfter);
    }
    return {
      row,
      // Every ordered row is a key of `recordsByRow`.
      records: /** @type {LinearRowRecord<S>[]} */ (recordsByRow.get(row)),
      boundaryAfter
    };
  });

  /** @type {Map<string, LinearComparisonIssue[]>} */
  const issuesByDraftId = new Map();
  validateLinearComparisonEdges({
    edges: normalizedPlan.edges,
    sequences: sequenceList,
    layout
  }).forEach((issue) => {
    if (!STRUCTURAL_PRESENTATION_ISSUES.has(issue.code)) return;
    if (!issuesByDraftId.has(issue.edgeId)) issuesByDraftId.set(issue.edgeId, []);
    // The line above makes sure the draft has its list.
    /** @type {LinearComparisonIssue[]} */ (issuesByDraftId.get(issue.edgeId)).push(issue);
  });

  const resolutionEdges = Array.isArray(resolution?.edges) ? resolution.edges : [];
  /** @type {Map<string, LinearComparisonResolvedEdge>} */
  const resolvedCandidateByEdgeKey = new Map();
  resolutionEdges.forEach((edge) => {
    if (!resolvedCandidateByEdgeKey.has(edge.edgeKey)) {
      resolvedCandidateByEdgeKey.set(edge.edgeKey, edge);
    }
  });
  /** @type {Map<string, LinearComparisonEdge[]>} */
  const draftsByEdgeKey = new Map();
  normalizedPlan.edges.forEach((draft) => {
    const edgeKey = linearComparisonEdgeKey(draft.queryUid, draft.subjectUid);
    if (!draftsByEdgeKey.has(edgeKey)) draftsByEdgeKey.set(edgeKey, []);
    // The line above makes sure the pair has its list.
    /** @type {LinearComparisonEdge[]} */ (draftsByEdgeKey.get(edgeKey)).push(draft);
  });
  /** @type {Map<string, LinearComparisonEdge>} */
  const ownerDraftByEdgeKey = new Map();
  draftsByEdgeKey.forEach((drafts, edgeKey) => {
    const resolvedId = cleanUid(resolvedCandidateByEdgeKey.get(edgeKey)?.id);
    ownerDraftByEdgeKey.set(
      edgeKey,
      drafts.find((draft) => draft.id === resolvedId) || drafts[0]
    );
  });

  /** @type {{ draft: LinearComparisonEdge, issues: LinearComparisonIssue[] }[]} */
  const unplacedDrafts = [];
  /** @type {Map<string, LinearComparisonEdge>} */
  const draftByEdgeKey = new Map();
  normalizedPlan.edges.forEach((draft) => {
    const edgeKey = linearComparisonEdgeKey(draft.queryUid, draft.subjectUid);
    const owner = ownerDraftByEdgeKey.get(edgeKey);
    const issues = (issuesByDraftId.get(draft.id) || [])
      .filter((issue) => issue.code !== 'duplicate');
    const boundary = pairBoundary({
      queryUid: draft.queryUid,
      subjectUid: draft.subjectUid,
      rowByUid,
      rowPositionByNumber,
      boundariesByKey
    });
    if (draft !== owner) {
      const fallbackIssues = [
        validationIssue(
          'duplicate',
          msg('Selected comparisons cannot contain the same directed pair more than once.'),
          draft
        ),
        ...issues
      ];
      unplacedDrafts.push({ draft, issues: fallbackIssues });
      return;
    }
    if (issues.length > 0 || !boundary) {
      const fallbackIssues = issues.length > 0 ? issues : [validationIssue(
        'non-adjacent',
        msg('Comparisons must connect adjacent rows.'),
        draft
      )];
      unplacedDrafts.push({ draft, issues: fallbackIssues });
      return;
    }
    draftByEdgeKey.set(edgeKey, draft);
  });

  /** @type {Map<string, LinearComparisonResolvedEdge>} */
  const resolvedByEdgeKey = new Map();
  resolutionEdges.forEach((edge) => {
    const edgeDrafts = draftsByEdgeKey.get(edge.edgeKey) || [];
    if (edgeDrafts.length > 0) {
      const owner = draftByEdgeKey.get(edge.edgeKey);
      if (!owner || owner.id !== edge.id) return;
    }
    if (!resolvedByEdgeKey.has(edge.edgeKey)) resolvedByEdgeKey.set(edge.edgeKey, edge);
  });
  /** @type {Map<string, LinearComparisonPair>} */
  const pairByEdgeKey = new Map();
  /** @param {{ queryUid: string, subjectUid: string, candidateKind: LinearComparisonCandidateKind }} input */
  const addPair = ({ queryUid, subjectUid, candidateKind }) => {
    const edgeKey = linearComparisonEdgeKey(queryUid, subjectUid);
    const boundary = pairBoundary({
      queryUid,
      subjectUid,
      rowByUid,
      rowPositionByNumber,
      boundariesByKey
    });
    if (!boundary || pairByEdgeKey.has(edgeKey)) return;
    const draft = draftByEdgeKey.get(edgeKey) || null;
    const resolved = resolvedByEdgeKey.get(edgeKey) || null;
    /** @type {LinearComparisonPair} */
    const pair = {
      edgeKey,
      edgeId: draft?.id || resolved?.id || `linear-derived-${edgeKey}`,
      queryUid,
      subjectUid,
      queryIndex: indexByUid.get(queryUid),
      subjectIndex: indexByUid.get(subjectUid),
      source: presentationSource({ normalizedPlan, candidateKind, draft, resolved }),
      active: Boolean(resolved),
      candidateKind,
      draft,
      resolved
    };
    boundary.pairs.push(pair);
    pairByEdgeKey.set(edgeKey, pair);
  };

  adjacentRowPairs(sequenceList, layout).forEach(([queryUid, subjectUid]) => {
    addPair({ queryUid, subjectUid, candidateKind: 'zipped' });
  });
  resolutionEdges.forEach((edge) => {
    if (resolvedByEdgeKey.get(edge.edgeKey) !== edge) return;
    addPair({
      queryUid: cleanUid(edge?.queryUid),
      subjectUid: cleanUid(edge?.subjectUid),
      candidateKind: 'resolved'
    });
  });
  normalizedPlan.edges.forEach((draft) => {
    if (draftByEdgeKey.get(linearComparisonEdgeKey(draft.queryUid, draft.subjectUid)) !== draft) return;
    addPair({
      queryUid: draft.queryUid,
      subjectUid: draft.subjectUid,
      candidateKind: 'draft'
    });
  });

  return { rows, unplacedDrafts };
};

/**
 * @param {{
 *   resolution?: LinearComparisonResolution,
 *   program?: string,
 *   blastpMode?: string
 * }} [options]
 * @returns {readonly LosatJobSpec[]}
 */
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

/**
 * The LOSAT edge that joins two drawn records: the one in the search's own
 * orientation, else the one in the other orientation, else null. Generate
 * (`buildLosatJobSpecs`) and Session Load (the cache-info rows) read a
 * search's edge here (OV-980, OV-987).
 * @param {readonly LinearComparisonResolvedEdge[]} edges The resolved edges.
 * @param {number} queryIndex
 * @param {number} subjectIndex
 * @returns {LinearComparisonResolvedEdge | null}
 */
export const losatEdgeOfPair = (edges, queryIndex, subjectIndex) => {
  const losatEdges = edges.filter((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT);
  return losatEdges.find((edge) => edge.queryIndex === queryIndex && edge.subjectIndex === subjectIndex)
    || losatEdges.find((edge) => edge.queryIndex === subjectIndex && edge.subjectIndex === queryIndex)
    || null;
};

// The record-pair searches one Generate requests. Generate and the Settings
// job-count estimate both expand them here (N-10). Similarity groups always
// include within-record evidence; Collinear includes it only when inference is
// on (PD-OI-018 item 4).
/**
 * @param {{
 *   resolution?: LinearComparisonResolution,
 *   recordCount?: number,
 *   recordUids?: string[],
 *   program?: string,
 *   blastpMode?: string,
 *   collinearInferOrthogroups?: boolean,
 *   collinearSearchScope?: string
 * }} [options]
 * @returns {readonly LosatJobSpec[]}
 */
export const buildLosatJobSpecs = ({
  resolution,
  recordCount = 0,
  recordUids = [],
  program = 'blastn',
  blastpMode = 'pairwise',
  collinearInferOrthogroups = true,
  collinearSearchScope = 'adjacent'
} = {}) => {
  const normalizedProgram = String(program || '').trim().toLowerCase();
  const mode = String(blastpMode || '').trim().toLowerCase();
  if (normalizedProgram !== 'blastp' || !['orthogroup', 'collinear'].includes(mode)) {
    return buildPairwiseLosatJobSpecs({ resolution, program: normalizedProgram, blastpMode: mode });
  }
  const edges = (resolution?.edges || []).filter((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT);
  const lastOrdinal = Math.max(0, edges.length - 1);
  /** @type {LosatJobSpec[]} */
  const specs = [];
  /**
   * @param {number} queryIndex
   * @param {number} subjectIndex
   * @param {number} ordinal
   */
  const push = (queryIndex, subjectIndex, ordinal) => {
    // The search's own edge is the card that shows it (OV-980) and numbers it
    // (OV-985); a pair no edge joins takes an edge by ordinal.
    const edge = losatEdgeOfPair(edges, queryIndex, subjectIndex)
      || edges.find((candidate) => candidate.ordinal === ordinal) || edges[Math.max(0, Math.min(ordinal, lastOrdinal))];
    specs.push(Object.freeze({
      edgeKey: edge?.edgeKey || '',
      ordinal: edge?.ordinal ?? ordinal,
      queryUid: recordUids[queryIndex],
      subjectUid: recordUids[subjectIndex],
      queryIndex,
      subjectIndex,
      program: normalizedProgram
    }));
  };
  /** @param {boolean} withSelf */
  const pushPairs = (withSelf) => {
    for (let i = 0; i < recordCount; i += 1) {
      if (withSelf) push(i, i, Math.min(i, lastOrdinal));
      for (let j = i + 1; j < recordCount; j += 1) {
        push(i, j, Math.min(i, lastOrdinal));
        push(j, i, Math.min(i, lastOrdinal));
      }
    }
  };
  if (mode === 'orthogroup') {
    pushPairs(true);
    return Object.freeze(specs);
  }
  if (collinearInferOrthogroups !== false) {
    for (let i = 0; i < recordCount; i += 1) push(i, i, Math.min(i, lastOrdinal));
  }
  if (collinearSearchScope === 'all') {
    pushPairs(false);
  } else {
    edges.forEach((edge) => {
      // Both callers pass a valid resolution (Generate requires one, the job estimate returns early on an invalid one), so every uid has a card.
      push(/** @type {number} */ (edge.queryIndex), /** @type {number} */ (edge.subjectIndex), edge.ordinal);
      push(/** @type {number} */ (edge.subjectIndex), /** @type {number} */ (edge.queryIndex), edge.ordinal);
    });
  }
  return Object.freeze(specs);
};

/**
 * Which searches of a LOSAT run are displayed pairs, and the number of each
 * search (a displayed pair's number names its payload, its `/blast_<n>.txt`
 * file and Run info slot). Every search takes its spec ordinal. A displayed
 * pair is its edge's search in the edge's direction or, for a complete row
 * grid, an adjacent-row pair; each such pair is its own ADJACENT edge
 * (`adjacentRowPairs`), so displayed numbers never repeat (OV-1069). The
 * derived LOSATP identity owner refuses a repeated displayed number.
 * @param {readonly LosatJobSpec[]} specs
 * @param {readonly LinearComparisonResolvedEdge[]} edges The edges of the comparison plan.
 * @param {((queryIndex: number, subjectIndex: number) => boolean) | null} gridPair
 *   `canonicalGridDisplayPairs` of a drawing laid out by a complete row grid, else null.
 * @returns {{ displayPair: boolean, pairIndex: number }[]} One per spec, in spec order.
 */
export const numberLosatDisplayPairs = (specs, edges, gridPair) => specs.map((spec) => {
  const edge = edges.find((candidate) => candidate.edgeKey === spec.edgeKey);
  const displayPair = gridPair
    ? gridPair(spec.queryIndex, spec.subjectIndex)
    : Boolean(edge && spec.queryIndex === edge.queryIndex && spec.subjectIndex === edge.subjectIndex);
  return { displayPair, pairIndex: spec.ordinal };
});
