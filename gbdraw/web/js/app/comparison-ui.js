// @ts-check
/** @import { LinearComparisonIssue, LinearComparisonPlan, LinearComparisonResolution, LinearComparisonResolvedEdge } from '../services/linear-comparisons.js' */
import {
  LINEAR_COMPARISON_MODES,
  LINEAR_COMPARISON_SOURCES,
  linearComparisonEdgeKey,
  normalizeLinearComparisonPlan
} from '../services/linear-comparisons.js';
import { comparisonStateForMode } from '../mode-profiles.js';
import { joinSentences, listText, msg, t, tm, tn } from '../utils/i18n.js';
/** @import { ComparisonState } from '../mode-profiles.js' */

const LINEAR_COMPARISON_INTENT_KEYS = Object.freeze({
  NONE: 'none',
  LOSAT: 'losat',
  UPLOAD: 'upload',
  CUSTOM: 'custom'
});

const LINEAR_COMPARISON_DISCLOSURE_KEYS = Object.freeze({
  SETTINGS: 'settings',
  SELECTED_PAIRS: 'selected-pairs'
});

const LINEAR_COMPARISON_LOSAT_MODE_KEYS = Object.freeze({
  LOSATN: 'blastn',
  LOSATP: 'blastp',
  TLOSATX: 'tblastx'
});

const LINEAR_COMPARISON_LOSATP_MODE_KEYS = Object.freeze({
  PAIRWISE: 'pairwise',
  SIMILARITY_GROUPS: 'orthogroup',
  COLLINEAR_BLOCKS: 'collinear'
});

const LINEAR_COMPARISON_SECTION_KEYS = Object.freeze({
  SETTINGS: Object.freeze({
    LOSAT_MODE: 'losat-mode',
    LOSAT_RUNTIME: 'losat-runtime',
    LOSATP_MODE: 'losatp-mode',
    LOSATP_CANDIDATE: 'losatp-candidate',
    BLASTN_TASK: 'blastn-task',
    RECORD_GENETIC_CODES: 'record-genetic-codes',
    BLASTP_MAX_HITS: 'blastp-max-hits',
    MEMBER_HITS: 'member-hits',
    COLLINEAR_PRIMARY: 'collinear-primary',
    COLLINEAR_SEARCH_SCOPE: 'collinear-search-scope',
    RESULT_FILTERS: 'result-filters',
    COMPARISON_APPEARANCE: 'comparison-appearance'
  }),
  SELECTED_PAIRS: Object.freeze({
    PAIR_EDITOR: 'pair-editor',
    RETAINED_DRAFTS: 'retained-drafts'
  }),
  ADVANCED: Object.freeze({
    RECORD_LAYOUT: 'record-layout',
    LOSAT_CACHE: 'losat-cache',
    COLLINEAR_DETAILS: 'collinear-details',
    PARALOG_LINKS: 'paralog-links'
  })
});

/** @typedef {Readonly<{ key: string, label: string, adjacentOnly?: boolean }>} ComparisonModeEntry */

/** @type {readonly ComparisonModeEntry[]} */
const LOSAT_MODES = Object.freeze([
  Object.freeze({
    key: LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATN,
    label: 'LOSATN'
  }),
  Object.freeze({
    key: LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATP,
    label: 'LOSATP'
  }),
  Object.freeze({
    key: LINEAR_COMPARISON_LOSAT_MODE_KEYS.TLOSATX,
    label: 'TLOSATX'
  })
]);

/** @type {readonly ComparisonModeEntry[]} */
const LOSATP_MODES = Object.freeze([
  Object.freeze({
    key: LINEAR_COMPARISON_LOSATP_MODE_KEYS.SIMILARITY_GROUPS,
    label: msg('Similarity groups'),
    adjacentOnly: true
  }),
  Object.freeze({
    key: LINEAR_COMPARISON_LOSATP_MODE_KEYS.COLLINEAR_BLOCKS,
    label: msg('Collinear blocks'),
    adjacentOnly: true
  }),
  Object.freeze({
    key: LINEAR_COMPARISON_LOSATP_MODE_KEYS.PAIRWISE,
    label: msg('Pairwise matches')
  })
]);

const LOSAT_MODE_BY_KEY = new Map(LOSAT_MODES.map((mode) => [mode.key, mode]));
const LOSATP_MODE_BY_KEY = new Map(LOSATP_MODES.map((mode) => [mode.key, mode]));

/** @type {Readonly<Record<string, Readonly<{ disclosureKey: string, focusTargetKey: string }>>>} */
const ISSUE_ROUTES = Object.freeze({
  'missing-upload': Object.freeze({
    disclosureKey: LINEAR_COMPARISON_DISCLOSURE_KEYS.SELECTED_PAIRS,
    focusTargetKey: 'pair-upload'
  }),
  'selected-losat-requires-pairwise': Object.freeze({
    disclosureKey: LINEAR_COMPARISON_DISCLOSURE_KEYS.SETTINGS,
    focusTargetKey: 'losatp-mode'
  }),
  duplicate: Object.freeze({
    disclosureKey: LINEAR_COMPARISON_DISCLOSURE_KEYS.SELECTED_PAIRS,
    focusTargetKey: 'pair-row'
  }),
  'missing-uid': Object.freeze({
    disclosureKey: LINEAR_COMPARISON_DISCLOSURE_KEYS.SELECTED_PAIRS,
    focusTargetKey: 'pair-row'
  }),
  self: Object.freeze({
    disclosureKey: LINEAR_COMPARISON_DISCLOSURE_KEYS.SELECTED_PAIRS,
    focusTargetKey: 'pair-row'
  }),
  'same-row': Object.freeze({
    disclosureKey: LINEAR_COMPARISON_DISCLOSURE_KEYS.SELECTED_PAIRS,
    focusTargetKey: 'pair-row'
  }),
  'non-adjacent': Object.freeze({
    disclosureKey: LINEAR_COMPARISON_DISCLOSURE_KEYS.SELECTED_PAIRS,
    focusTargetKey: 'pair-row'
  })
});

// Why an All adjacent pairs plan that resolves no pair draws no comparison (CI-07a, CI-07b).
/** @type {Readonly<Record<string, string>>} */
const EMPTY_PLAN_NOTICES = Object.freeze({
  [LINEAR_COMPARISON_INTENT_KEYS.LOSAT]: msg('No record pair to compare. Records on the same Linear row are not compared: give each record its own row under Advanced comparison and layout (Record Layout), or add a sequence.'),
  [LINEAR_COMPARISON_INTENT_KEYS.UPLOAD]: msg('No pair has a BLAST TSV yet, so Generate draws no comparison. Choose a file for a pair under Selected pairs.')
});

/**
 * The four result filters of a comparison; each value is unvalidated.
 * @typedef {Partial<Record<keyof ComparisonState, unknown>>} ComparisonFilterValues
 */

/**
 * @param {Record<string, unknown>} value
 * @param {string} key
 */
const own = (value, key) => Object.prototype.hasOwnProperty.call(value, key);

/** @param {unknown} value */
const normalizeLosatProgram = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['blastn', 'tblastx', 'blastp'].includes(normalized) ? normalized : 'blastn';
};

/** @param {unknown} value */
const normalizeBlastpMode = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['pairwise', 'orthogroup', 'collinear'].includes(normalized)
    ? normalized
    : 'orthogroup';
};

// A Selected plan whose included pairs are exactly the adjacent pairs, all
// uploads, is the Upload intent: uploading the table of a pair makes it
// Selected (UJ-05). Selected LOSAT stays Custom, because the Selected
// topology limits the LOSATP modes.
/**
 * @param {ReturnType<typeof normalizeLinearComparisonPlan>} plan
 * @param {readonly string[]} adjacentEdgeKeys
 */
const isAdjacentUploadPlan = (plan, adjacentEdgeKeys) => {
  const included = plan.edges.filter((edge) => edge.included === true);
  const keys = new Set(included.map((edge) => linearComparisonEdgeKey(edge.queryUid, edge.subjectUid)));
  const adjacent = new Set(adjacentEdgeKeys);
  return adjacent.size > 0 && keys.size === included.length && keys.size === adjacent.size
    && [...keys].every((key) => adjacent.has(key))
    && included.every((edge) => edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD);
};

/**
 * @param {ReturnType<typeof normalizeLinearComparisonPlan>} plan
 * @param {readonly string[]} adjacentEdgeKeys
 */
const intentKeyForPlan = (plan, adjacentEdgeKeys) => {
  if (plan.mode === LINEAR_COMPARISON_MODES.NONE) {
    return LINEAR_COMPARISON_INTENT_KEYS.NONE;
  }
  if (plan.mode === LINEAR_COMPARISON_MODES.SELECTED) {
    return isAdjacentUploadPlan(plan, adjacentEdgeKeys)
      ? LINEAR_COMPARISON_INTENT_KEYS.UPLOAD
      : LINEAR_COMPARISON_INTENT_KEYS.CUSTOM;
  }
  return plan.defaultSource === LINEAR_COMPARISON_SOURCES.UPLOAD
    ? LINEAR_COMPARISON_INTENT_KEYS.UPLOAD
    : LINEAR_COMPARISON_INTENT_KEYS.LOSAT;
};

/** @param {LinearComparisonPlan} plan */
const plannedSources = (plan) => {
  if (plan.mode === LINEAR_COMPARISON_MODES.NONE) {
    return Object.freeze({ hasLosat: false, hasUpload: false });
  }
  if (plan.mode === LINEAR_COMPARISON_MODES.ADJACENT) {
    return Object.freeze({
      hasLosat: plan.defaultSource === LINEAR_COMPARISON_SOURCES.LOSAT,
      hasUpload: plan.defaultSource === LINEAR_COMPARISON_SOURCES.UPLOAD
    });
  }
  const activeDrafts = plan.edges.filter((edge) => edge.included === true);
  return Object.freeze({
    hasLosat: activeDrafts.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.LOSAT),
    hasUpload: activeDrafts.some((edge) => edge.source === LINEAR_COMPARISON_SOURCES.UPLOAD)
  });
};

/** @param {readonly LinearComparisonResolvedEdge[]} edges */
const projectSourceBreakdown = (edges) => {
  const losat = edges.filter((edge) => edge?.source === LINEAR_COMPARISON_SOURCES.LOSAT).length;
  const upload = edges.filter((edge) => edge?.source === LINEAR_COMPARISON_SOURCES.UPLOAD).length;
  const key = losat > 0 && upload > 0
    ? 'mixed'
    : losat > 0
      ? LINEAR_COMPARISON_SOURCES.LOSAT
      : upload > 0
        ? LINEAR_COMPARISON_SOURCES.UPLOAD
        : 'none';
  const label = key === 'mixed'
    ? t('{losat} LOSAT, {upload} upload', { losat, upload })
    : key === LINEAR_COMPARISON_SOURCES.LOSAT
      ? t('{losat} LOSAT', { losat })
      : key === LINEAR_COMPARISON_SOURCES.UPLOAD
        ? t('{upload} upload', { upload })
        : t('No active pairs');
  return Object.freeze({ key, losat, upload, label });
};

/**
 * @param {LinearComparisonPlan} plan
 * @param {readonly LinearComparisonResolvedEdge[]} activeEdges
 */
const retainedDormantDraftCount = (plan, activeEdges) => {
  const activeById = new Map(activeEdges.map((edge) => [String(edge?.id || ''), edge]));
  return plan.edges.filter((draft) => {
    const active = activeById.get(String(draft.id || ''));
    const retainedFileIsActive = Boolean(
      active &&
      active.source === LINEAR_COMPARISON_SOURCES.UPLOAD &&
      draft.file &&
      draft.fileActive === true &&
      active.file === draft.file
    );
    const retainedFilenameIsActive = Boolean(
      active &&
      active.source === LINEAR_COMPARISON_SOURCES.LOSAT &&
      String(draft.losatFilename || '').trim() &&
      draft.losatFilenameActive === true &&
      String(active.losatFilename || '') === String(draft.losatFilename || '')
    );
    const hasDormantFile = Boolean(draft.file) && !retainedFileIsActive;
    const hasDormantFilename = Boolean(String(draft.losatFilename || '').trim()) &&
      !retainedFilenameIsActive;
    return hasDormantFile || hasDormantFilename;
  }).length;
};

/**
 * @param {unknown} left
 * @param {unknown} right
 */
const valuesEquivalent = (left, right) => {
  const leftNumber = Number(left);
  const rightNumber = Number(right);
  if (
    String(left ?? '').trim() !== '' &&
    String(right ?? '').trim() !== '' &&
    Number.isFinite(leftNumber) &&
    Number.isFinite(rightNumber)
  ) {
    return leftNumber === rightNumber;
  }
  return String(left ?? '').trim() === String(right ?? '').trim();
};

/** @param {unknown} value */
const displayFilterValue = (value) => {
  if (value === null || value === undefined || String(value).trim() === '') return t('unset');
  return String(value).trim();
};

/** @param {ComparisonFilterValues} filters */
const projectFilterSummary = (filters) => {
  const defaults = comparisonStateForMode('linear');
  const source = filters && typeof filters === 'object' && !Array.isArray(filters)
    ? filters
    : {};
  const values = Object.freeze({
    evalue: own(source, 'evalue') ? source.evalue : defaults.evalue,
    min_bitscore: own(source, 'min_bitscore') ? source.min_bitscore : defaults.min_bitscore,
    identity: own(source, 'identity') ? source.identity : defaults.identity,
    alignment_length: own(source, 'alignment_length')
      ? source.alignment_length
      : defaults.alignment_length
  });
  const text = [
    t('E-value <= {value}', { value: displayFilterValue(values.evalue) }),
    t('Bitscore >= {value}', { value: displayFilterValue(values.min_bitscore) }),
    t('Identity >= {value}%', { value: displayFilterValue(values.identity) }),
    t('Length >= {value}', { value: displayFilterValue(values.alignment_length) })
  ].join(' · ');
  // `values` has exactly the four keys of ComparisonState, built above.
  const isDefault = /** @type {(keyof ComparisonState)[]} */ (Object.keys(values)).every((key) => valuesEquivalent(values[key], defaults[key]));
  return Object.freeze({ text, isDefault });
};

/** @param {string} activeKey */
const projectLosatModes = (activeKey) => Object.freeze(LOSAT_MODES.map((mode) => Object.freeze({
  key: mode.key,
  label: tm(mode.label),
  active: mode.key === activeKey,
  selectable: true,
  disabledReasonKey: null,
  requiredIntentActionKey: null
})));

/**
 * @param {LinearComparisonPlan} plan
 * @param {string} activeKey
 */
const projectLosatpModes = (plan, activeKey) => Object.freeze(LOSATP_MODES.map((mode) => {
  const selectable = !(plan.mode === LINEAR_COMPARISON_MODES.SELECTED && mode.adjacentOnly);
  return Object.freeze({
    key: mode.key,
    label: tm(mode.label),
    active: mode.key === activeKey,
    selectable,
    disabledReasonKey: selectable ? null : 'requires-adjacent-losat',
    requiredIntentActionKey: selectable ? null : LINEAR_COMPARISON_INTENT_KEYS.LOSAT
  });
}));

/**
 * @typedef {object} SettingsSectionInput
 * @property {string} intentKey
 * @property {ReturnType<typeof plannedSources>} planned
 * @property {string} losatModeKey
 * @property {string} losatpModeKey
 * @property {boolean} selectedTopology
 * @property {boolean} inferOrthogroups Collinear builds similarity groups.
 * @property {string} anchorMode The Collinear anchor mode.
 */

/** @param {SettingsSectionInput} input */
const projectSettingsSectionKeys = ({
  intentKey,
  planned,
  losatModeKey,
  losatpModeKey,
  selectedTopology,
  inferOrthogroups,
  anchorMode
}) => {
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.NONE) return Object.freeze([]);

  // In the order the Settings disclosure shows them: search limits and
  // Performance come after the result filters and appearance.
  const keys = [];
  const searchLimits = [];
  if (planned.hasLosat) {
    keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.LOSAT_MODE);
    if (losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATP) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.LOSATP_MODE);
      searchLimits.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.LOSATP_CANDIDATE);
    }
    if (losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATN) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.BLASTN_TASK);
    } else if (losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.TLOSATX) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.RECORD_GENETIC_CODES);
    } else if (losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.PAIRWISE) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.BLASTP_MAX_HITS);
    } else if (
      !selectedTopology &&
      losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.SIMILARITY_GROUPS
    ) {
      searchLimits.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.MEMBER_HITS);
    } else if (!selectedTopology && losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.COLLINEAR_BLOCKS) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.COLLINEAR_PRIMARY);
      // Shown only where they can change the result (OV-2302): the searched
      // pairs feed the similarity groups only, and the member limit reaches
      // the anchors through the groups or All hits. Hidden values stay in
      // the draft.
      if (inferOrthogroups) keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.COLLINEAR_SEARCH_SCOPE);
      if (inferOrthogroups || anchorMode === 'all') searchLimits.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.MEMBER_HITS);
    }
  }
  keys.push(
    LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.RESULT_FILTERS,
    LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.COMPARISON_APPEARANCE,
    ...searchLimits
  );
  if (planned.hasLosat) keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.LOSAT_RUNTIME);
  return Object.freeze(keys);
};

/** @param {SettingsSectionInput & { dormantDraftCount: number }} input */
const projectSectionKeys = ({
  intentKey,
  planned,
  losatModeKey,
  losatpModeKey,
  selectedTopology,
  dormantDraftCount,
  inferOrthogroups,
  anchorMode
}) => {
  const settings = projectSettingsSectionKeys({
    intentKey,
    planned,
    losatModeKey,
    losatpModeKey,
    selectedTopology,
    inferOrthogroups,
    anchorMode
  });
  const selectedPairs = [];
  if (intentKey !== LINEAR_COMPARISON_INTENT_KEYS.NONE || dormantDraftCount > 0) {
    selectedPairs.push(LINEAR_COMPARISON_SECTION_KEYS.SELECTED_PAIRS.PAIR_EDITOR);
  }
  if (dormantDraftCount > 0) {
    selectedPairs.push(LINEAR_COMPARISON_SECTION_KEYS.SELECTED_PAIRS.RETAINED_DRAFTS);
  }

  /** @type {string[]} */
  const advanced = [LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.RECORD_LAYOUT];
  if (planned.hasLosat) {
    advanced.push(LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.LOSAT_CACHE);
  }
  const adjacentLosatp = planned.hasLosat &&
    !selectedTopology &&
    losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATP;
  const collinear = adjacentLosatp && losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.COLLINEAR_BLOCKS;
  if (collinear) advanced.push(LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.COLLINEAR_DETAILS);
  // Paralog links per group shows wherever groups are drawn with it (OV-1128 B).
  if (
    (adjacentLosatp && losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.SIMILARITY_GROUPS)
    || (collinear && inferOrthogroups)
  ) {
    advanced.push(LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.PARALOG_LINKS);
  }
  return Object.freeze({
    settings,
    selectedPairs: Object.freeze(selectedPairs),
    advanced: Object.freeze(advanced)
  });
};

// One entry per issue, worded at its cause, with the actions that fix it.
// The Comparison card lists them once; the pair row of a pair issue
// carries its short text. An action's accessible name adds its pair to the
// visible label.
/** @param {string} code @param {string} pairLabel */
const issueActions = (code, pairLabel) => Object.freeze((code === 'missing-upload'
  ? [{ key: 'choose-file', label: t('Choose file'), ariaLabel: pairLabel ? t('Choose file for {pair}', { pair: pairLabel }) : t('Choose file') }]
  : code === 'selected-losat-requires-pairwise'
    ? [{ key: 'all-losat', label: t('Run LOSAT for all adjacent pairs') }, { key: 'pairwise', label: t('Use Pairwise matches') }]
    : [{ key: 'show', label: t('Show pair'), ariaLabel: pairLabel ? t('Show pair {pair}', { pair: pairLabel }) : t('Show pair') }]
).map((action) => Object.freeze({ ariaLabel: action.label, ...action })));

// Number.isInteger is true only for a number, so a defined edge with both
// indexes is what the casts in the template read.
/** @param {LinearComparisonResolvedEdge | undefined} edge */
const pairLabelOf = (edge) => (Number.isInteger(edge?.queryIndex) && Number.isInteger(edge?.subjectIndex)
  ? `#${/** @type {{ queryIndex: number }} */ (edge).queryIndex + 1} → #${/** @type {{ subjectIndex: number }} */ (edge).subjectIndex + 1}` : '');

// A Selected plan breaks the LOSATP mode: the text names the pairs set to
// Upload BLAST TSV, or says that the pairs are set one by one when every one
// of them runs LOSAT.
/** @param {readonly LinearComparisonResolvedEdge[]} edges */
const selectedPairsText = (edges) => {
  const uploads = edges.filter((edge) => edge?.source === LINEAR_COMPARISON_SOURCES.UPLOAD).map(pairLabelOf).filter(Boolean);
  const cause = uploads.length === 0 ? t('The pairs are set one by one in Selected pairs.')
    : tn(
      uploads.length,
      'Pair {pairs} is set to Upload BLAST TSV, so not every pair runs LOSAT.',
      'Pairs {pairs} are set to Upload BLAST TSV, so not every pair runs LOSAT.',
      { pairs: listText(uploads) }
    );
  return joinSentences(cause, t('Similarity groups and Collinear blocks compare every adjacent pair with LOSAT.'));
};

/** @param {readonly LinearComparisonIssue[] | undefined} errors @param {readonly LinearComparisonResolvedEdge[]} edges */
const projectIssues = (errors, edges) => {
  const edgeList = Array.isArray(edges) ? edges : [];
  const edgeByKey = new Map(edgeList.map((edge) => [String(edge?.edgeKey || ''), edge]));
  return Object.freeze((Array.isArray(errors) ? errors : []).flatMap((issue, index) => {
    const code = String(issue?.code || '').trim();
    const route = ISSUE_ROUTES[code];
    if (!route) return [];
    const edgeKey = String(issue?.edgeKey || '');
    const pairLabel = pairLabelOf(edgeByKey.get(edgeKey));
    const pairText = code === 'missing-upload' ? t('No BLAST file') : tm(String(issue?.message || ''));
    const text = code === 'selected-losat-requires-pairwise' ? selectedPairsText(edgeList)
      : pairLabel ? t('{pair}: {text}.', { pair: pairLabel, text: pairText }).replace(/\.\.$/, '.') : tm(String(issue?.message || ''));
    return [Object.freeze({
      key: `${index}-${code}`,
      issueCode: code,
      disclosureKey: route.disclosureKey,
      focusTargetKey: route.focusTargetKey,
      edgeId: String(issue?.edgeId || ''),
      edgeKey,
      text,
      pairText: pairLabel ? pairText : '',
      actions: issueActions(code, pairLabel)
    })];
  }));
};

/** @typedef {ReturnType<typeof projectIssues>[number]} LinearComparisonUiIssue An entry of the Comparison card's issue list. */

// Each pair row shows every issue of its pair.
/** @param {readonly LinearComparisonUiIssue[]} issues */
const projectPairIssueText = (issues) => {
  /** @type {Record<string, string>} */
  const texts = {};
  for (const issue of issues) {
    if (!issue.edgeKey || !issue.pairText) continue;
    texts[issue.edgeKey] = texts[issue.edgeKey] ? `${texts[issue.edgeKey]}; ${issue.pairText}` : issue.pairText;
  }
  return Object.freeze(texts);
};

/**
 * @param {LinearComparisonPlan} plan
 * @param {number} activePairCount
 */
const projectTopology = (plan, activePairCount) => {
  if (plan.mode === LINEAR_COMPARISON_MODES.NONE) {
    return Object.freeze({
      key: LINEAR_COMPARISON_MODES.NONE,
      label: t('No comparison'),
      pairLabel: t('No comparison')
    });
  }
  if (plan.mode === LINEAR_COMPARISON_MODES.SELECTED) {
    return Object.freeze({
      key: LINEAR_COMPARISON_MODES.SELECTED,
      label: t('Selected pairs'),
      pairLabel: tn(activePairCount, '{count} selected pair', '{count} selected pairs')
    });
  }
  return Object.freeze({
    key: LINEAR_COMPARISON_MODES.ADJACENT,
    label: t('All adjacent pairs'),
    pairLabel: tn(activePairCount, '{count} adjacent pair', '{count} adjacent pairs')
  });
};

export const projectLinearComparisonLosatModeSelection = ({
  modeKey = ''
} = {}) => {
  const mode = LOSAT_MODE_BY_KEY.get(String(modeKey || '').trim());
  if (!mode) {
    return Object.freeze({
      selectable: false,
      reasonKey: 'unknown-losat-mode',
      requiredIntentActionKey: null,
      patch: null
    });
  }
  return Object.freeze({
    selectable: true,
    reasonKey: null,
    requiredIntentActionKey: null,
    patch: Object.freeze({ losatProgram: mode.key })
  });
};

export const projectLinearComparisonLosatpModeSelection = ({
  plan = {},
  modeKey = ''
} = {}) => {
  const normalizedPlan = normalizeLinearComparisonPlan(plan);
  const mode = LOSATP_MODE_BY_KEY.get(String(modeKey || '').trim());
  if (!mode) {
    return Object.freeze({
      selectable: false,
      reasonKey: 'unknown-losatp-mode',
      requiredIntentActionKey: null,
      patch: null
    });
  }
  if (normalizedPlan.mode === LINEAR_COMPARISON_MODES.SELECTED && mode.adjacentOnly) {
    return Object.freeze({
      selectable: false,
      reasonKey: 'requires-adjacent-losat',
      requiredIntentActionKey: LINEAR_COMPARISON_INTENT_KEYS.LOSAT,
      patch: null
    });
  }
  return Object.freeze({
    selectable: true,
    reasonKey: null,
    requiredIntentActionKey: null,
    patch: Object.freeze({ blastpMode: mode.key })
  });
};

/**
 * @param {{
 *   plan?: unknown,
 *   resolution?: Partial<LinearComparisonResolution>,
 *   adjacentEdgeKeys?: readonly string[],
 *   losatProgram?: string,
 *   blastpMode?: string,
 *   inferOrthogroups?: boolean,
 *   anchorMode?: string,
 *   filters?: ComparisonFilterValues
 * }} [options]
 */
export const projectLinearComparisonUi = ({
  plan = {},
  resolution = {},
  adjacentEdgeKeys = [],
  losatProgram = 'blastn',
  blastpMode = 'orthogroup',
  inferOrthogroups = false,
  anchorMode = 'rbh',
  filters = {}
} = {}) => {
  const normalizedPlan = normalizeLinearComparisonPlan(plan);
  const activeEdges = Array.isArray(resolution?.edges) ? resolution.edges : [];
  const intentKey = intentKeyForPlan(normalizedPlan, adjacentEdgeKeys);
  const planned = plannedSources(normalizedPlan);
  const selectedTopology = normalizedPlan.mode === LINEAR_COMPARISON_MODES.SELECTED;
  const losatModeKey = normalizeLosatProgram(losatProgram);
  const losatpModeKey = normalizeBlastpMode(blastpMode);
  // normalizeLosatProgram and normalizeBlastpMode return only keys of these two tables.
  const losatModeLabel = /** @type {ComparisonModeEntry} */ (LOSAT_MODE_BY_KEY.get(losatModeKey)).label;
  const losatpModeLabel = /** @type {ComparisonModeEntry} */ (LOSATP_MODE_BY_KEY.get(losatpModeKey)).label;
  const sourceBreakdown = projectSourceBreakdown(activeEdges);
  const activePairCount = activeEdges.length;
  const dormantDraftCount = retainedDormantDraftCount(normalizedPlan, activeEdges);
  const topology = projectTopology(normalizedPlan, activePairCount);
  const filterSummary = projectFilterSummary(filters);
  const sectionKeys = projectSectionKeys({
    intentKey,
    planned,
    losatModeKey,
    losatpModeKey,
    selectedTopology,
    dormantDraftCount,
    inferOrthogroups,
    anchorMode
  });
  const issues = projectIssues(resolution?.errors, activeEdges);

  return Object.freeze({
    intentKey,
    topologyKey: topology.key,
    topologyLabel: topology.label,
    topologyPairLabel: topology.pairLabel,
    activePairCount,
    sourceBreakdown,
    activeLosatModeKey: losatModeKey,
    activeLosatModeLabel: losatModeLabel,
    losatModes: projectLosatModes(losatModeKey),
    activeLosatpModeKey: losatpModeKey,
    activeLosatpModeLabel: losatpModeLabel,
    losatpModes: projectLosatpModes(normalizedPlan, losatpModeKey),
    emptyPlanNotice: activePairCount === 0 && EMPTY_PLAN_NOTICES[intentKey] ? tm(EMPTY_PLAN_NOTICES[intentKey]) : '',
    filterSummary: filterSummary.text,
    filterSummaryIsDefault: filterSummary.isDefault,
    retainedDormantDraftCount: dormantDraftCount,
    sectionKeys,
    issues,
    pairIssueText: projectPairIssueText(issues)
  });
};
