import {
  LINEAR_COMPARISON_MODES,
  LINEAR_COMPARISON_SOURCES,
  normalizeLinearComparisonPlan
} from './linear-comparisons.js';
import { comparisonStateForMode } from '../mode-profiles.js';

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
    LOSATP_MODE: 'losatp-mode',
    BLASTN_TASK: 'blastn-task',
    RECORD_GENETIC_CODES: 'record-genetic-codes',
    BLASTP_MAX_HITS: 'blastp-max-hits',
    SIMILARITY_GROUPING: 'similarity-grouping',
    COLLINEAR_PRIMARY: 'collinear-primary',
    UPLOAD_READINESS: 'upload-readiness',
    RESULT_FILTERS: 'result-filters',
    COMPARISON_APPEARANCE: 'comparison-appearance'
  }),
  SELECTED_PAIRS: Object.freeze({
    PAIR_EDITOR: 'pair-editor',
    RETAINED_DRAFTS: 'retained-drafts'
  }),
  ADVANCED: Object.freeze({
    RECORD_LAYOUT: 'record-layout',
    LOSAT_RUNTIME: 'losat-runtime',
    LOSAT_CACHE: 'losat-cache',
    RAW_RESULTS: 'raw-results',
    COLLINEAR_DETAILS: 'collinear-details'
  })
});

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

const LOSATP_MODES = Object.freeze([
  Object.freeze({
    key: LINEAR_COMPARISON_LOSATP_MODE_KEYS.SIMILARITY_GROUPS,
    label: 'Similarity groups',
    adjacentOnly: true
  }),
  Object.freeze({
    key: LINEAR_COMPARISON_LOSATP_MODE_KEYS.COLLINEAR_BLOCKS,
    label: 'Collinear blocks',
    adjacentOnly: true
  }),
  Object.freeze({
    key: LINEAR_COMPARISON_LOSATP_MODE_KEYS.PAIRWISE,
    label: 'Pairwise matches'
  })
]);

const LOSAT_MODE_BY_KEY = new Map(LOSAT_MODES.map((mode) => [mode.key, mode]));
const LOSATP_MODE_BY_KEY = new Map(LOSATP_MODES.map((mode) => [mode.key, mode]));

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

const own = (value, key) => Object.prototype.hasOwnProperty.call(value, key);

const pluralizedPairs = (count, qualifier) => (
  `${count} ${qualifier} pair${count === 1 ? '' : 's'}`
);

const normalizeLosatProgram = (value) => {
  const normalized = String(value || '').trim().toLowerCase();
  return ['blastn', 'tblastx', 'blastp'].includes(normalized) ? normalized : 'blastn';
};

const intentKeyForPlan = (plan) => {
  if (plan.mode === LINEAR_COMPARISON_MODES.NONE) {
    return LINEAR_COMPARISON_INTENT_KEYS.NONE;
  }
  if (plan.mode === LINEAR_COMPARISON_MODES.SELECTED) {
    return LINEAR_COMPARISON_INTENT_KEYS.CUSTOM;
  }
  return plan.defaultSource === LINEAR_COMPARISON_SOURCES.UPLOAD
    ? LINEAR_COMPARISON_INTENT_KEYS.UPLOAD
    : LINEAR_COMPARISON_INTENT_KEYS.LOSAT;
};

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
    ? `${losat} LOSAT, ${upload} upload`
    : key === LINEAR_COMPARISON_SOURCES.LOSAT
      ? `${losat} LOSAT`
      : key === LINEAR_COMPARISON_SOURCES.UPLOAD
        ? `${upload} upload`
        : 'No active pairs';
  return Object.freeze({ key, losat, upload, label });
};

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

const retainedRawResultDraftCount = (plan, activeEdges) => {
  const activeById = new Map(activeEdges.map((edge) => [String(edge?.id || ''), edge]));
  return plan.edges.filter((draft) => {
    const retainedName = String(draft.losatFilename || '').trim();
    if (!retainedName) return false;
    const active = activeById.get(String(draft.id || ''));
    return !(
      active?.source === LINEAR_COMPARISON_SOURCES.LOSAT &&
      draft.losatFilenameActive === true &&
      String(active.losatFilename || '') === retainedName
    );
  }).length;
};

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

const displayFilterValue = (value) => {
  if (value === null || value === undefined || String(value).trim() === '') return 'unset';
  return String(value).trim();
};

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
    `E-value <= ${displayFilterValue(values.evalue)}`,
    `Bitscore >= ${displayFilterValue(values.min_bitscore)}`,
    `Identity >= ${displayFilterValue(values.identity)}%`,
    `Length >= ${displayFilterValue(values.alignment_length)}`
  ].join(' · ');
  const isDefault = Object.keys(values).every((key) => valuesEquivalent(values[key], defaults[key]));
  return Object.freeze({ text, isDefault });
};

const projectLosatModes = (activeKey) => Object.freeze(LOSAT_MODES.map((mode) => Object.freeze({
  key: mode.key,
  label: mode.label,
  active: mode.key === activeKey,
  selectable: true,
  disabledReasonKey: null,
  requiredIntentActionKey: null
})));

const projectLosatpModes = (plan, activeKey) => Object.freeze(LOSATP_MODES.map((mode) => {
  const selectable = !(plan.mode === LINEAR_COMPARISON_MODES.SELECTED && mode.adjacentOnly);
  return Object.freeze({
    key: mode.key,
    label: mode.label,
    active: mode.key === activeKey,
    selectable,
    disabledReasonKey: selectable ? null : 'requires-adjacent-losat',
    requiredIntentActionKey: selectable ? null : LINEAR_COMPARISON_INTENT_KEYS.LOSAT
  });
}));

const projectSettingsSectionKeys = ({
  intentKey,
  planned,
  losatModeKey,
  losatpModeKey,
  selectedTopology
}) => {
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.NONE) return Object.freeze([]);

  const keys = [];
  if (planned.hasLosat) {
    keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.LOSAT_MODE);
    if (losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATP) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.LOSATP_MODE);
    }
    if (losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATN) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.BLASTN_TASK);
    } else if (losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.TLOSATX) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.RECORD_GENETIC_CODES);
    } else if (losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.PAIRWISE) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.BLASTP_MAX_HITS);
    } else if (!selectedTopology && losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.SIMILARITY_GROUPS) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.SIMILARITY_GROUPING);
    } else if (!selectedTopology && losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.COLLINEAR_BLOCKS) {
      keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.COLLINEAR_PRIMARY);
    }
  }
  if (planned.hasUpload) {
    keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.UPLOAD_READINESS);
  }
  keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.RESULT_FILTERS);
  if (
    planned.hasUpload ||
    losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATN ||
    losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.TLOSATX ||
    losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.PAIRWISE
  ) {
    keys.push(LINEAR_COMPARISON_SECTION_KEYS.SETTINGS.COMPARISON_APPEARANCE);
  }
  return Object.freeze(keys);
};

const projectSectionKeys = ({
  intentKey,
  planned,
  losatModeKey,
  losatpModeKey,
  selectedTopology,
  dormantDraftCount,
  dormantRawResultDraftCount
}) => {
  const settings = projectSettingsSectionKeys({
    intentKey,
    planned,
    losatModeKey,
    losatpModeKey,
    selectedTopology
  });
  const selectedPairs = [];
  if (intentKey !== LINEAR_COMPARISON_INTENT_KEYS.NONE || dormantDraftCount > 0) {
    selectedPairs.push(LINEAR_COMPARISON_SECTION_KEYS.SELECTED_PAIRS.PAIR_EDITOR);
  }
  if (dormantDraftCount > 0) {
    selectedPairs.push(LINEAR_COMPARISON_SECTION_KEYS.SELECTED_PAIRS.RETAINED_DRAFTS);
  }

  const advanced = [LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.RECORD_LAYOUT];
  if (planned.hasLosat) {
    advanced.push(
      LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.LOSAT_RUNTIME,
      LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.LOSAT_CACHE
    );
  }
  if (planned.hasLosat || dormantRawResultDraftCount > 0) {
    advanced.push(LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.RAW_RESULTS);
  }
  if (
    planned.hasLosat &&
    !selectedTopology &&
    losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATP &&
    losatpModeKey === LINEAR_COMPARISON_LOSATP_MODE_KEYS.COLLINEAR_BLOCKS
  ) {
    advanced.push(LINEAR_COMPARISON_SECTION_KEYS.ADVANCED.COLLINEAR_DETAILS);
  }
  return Object.freeze({
    settings,
    selectedPairs: Object.freeze(selectedPairs),
    advanced: Object.freeze(advanced)
  });
};

const projectErrorRouting = (errors) => {
  const targets = Object.freeze((Array.isArray(errors) ? errors : []).flatMap((issue) => {
    const code = String(issue?.code || '').trim();
    const route = ISSUE_ROUTES[code];
    if (!route) return [];
    return [Object.freeze({
      issueCode: code,
      disclosureKey: route.disclosureKey,
      focusTargetKey: route.focusTargetKey,
      edgeId: String(issue?.edgeId || ''),
      edgeKey: String(issue?.edgeKey || '')
    })];
  }));
  const disclosureKeys = Object.freeze([
    ...new Set(targets.map((target) => target.disclosureKey))
  ]);
  const count = targets.length;
  return Object.freeze({
    badge: count === 0
      ? null
      : Object.freeze({
        count,
        label: `${count} comparison issue${count === 1 ? '' : 's'}`
      }),
    primaryDisclosureKey: targets[0]?.disclosureKey || null,
    disclosureKeys,
    targets
  });
};

const projectTopology = (plan, activePairCount) => {
  if (plan.mode === LINEAR_COMPARISON_MODES.NONE) {
    return Object.freeze({
      key: LINEAR_COMPARISON_MODES.NONE,
      label: 'No comparison',
      pairLabel: 'No comparison'
    });
  }
  if (plan.mode === LINEAR_COMPARISON_MODES.SELECTED) {
    return Object.freeze({
      key: LINEAR_COMPARISON_MODES.SELECTED,
      label: 'Selected pairs',
      pairLabel: pluralizedPairs(activePairCount, 'selected')
    });
  }
  return Object.freeze({
    key: LINEAR_COMPARISON_MODES.ADJACENT,
    label: 'All adjacent pairs',
    pairLabel: pluralizedPairs(activePairCount, 'adjacent')
  });
};

const projectStatusLabel = ({ intentKey, activePairCount, sourceBreakdown }) => {
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.NONE) return 'Current: No comparison';
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.LOSAT) {
    return 'Current: Run LOSAT for all adjacent pairs';
  }
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.UPLOAD) {
    return 'Current: Upload BLAST TSV for all adjacent pairs';
  }
  const sources = sourceBreakdown.key === 'none' ? '' : `; ${sourceBreakdown.label}`;
  return `Current: Selected pairs (${activePairCount}${sources})`;
};

const projectSummaryText = ({
  intentKey,
  topology,
  planned,
  sourceBreakdown,
  losatMethodSummaryLabel,
  filterSummary
}) => {
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.NONE) return 'No comparison';
  const parts = [];
  if (planned.hasLosat) parts.push(losatMethodSummaryLabel);
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.UPLOAD) parts.push('Upload BLAST TSV');
  parts.push(topology.pairLabel);
  if (intentKey === LINEAR_COMPARISON_INTENT_KEYS.CUSTOM) {
    parts.push(sourceBreakdown.label);
  }
  parts.push(filterSummary.text);
  return parts.join(' · ');
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

export const projectLinearComparisonUi = ({
  plan = {},
  resolution = {},
  losatProgram = 'blastn',
  blastpMode = '',
  filters = {}
} = {}) => {
  const normalizedPlan = normalizeLinearComparisonPlan(plan);
  const activeEdges = Array.isArray(resolution?.edges) ? resolution.edges : [];
  const intentKey = intentKeyForPlan(normalizedPlan);
  const planned = plannedSources(normalizedPlan);
  const selectedTopology = normalizedPlan.mode === LINEAR_COMPARISON_MODES.SELECTED;
  const losatModeKey = normalizeLosatProgram(losatProgram);
  const losatpModeKey = String(blastpMode || '');
  const losatModeLabel = LOSAT_MODE_BY_KEY.get(losatModeKey).label;
  const losatpModeLabel = LOSATP_MODE_BY_KEY.get(losatpModeKey)?.label || 'Unrecognized mode';
  const losatMethodSummaryLabel = losatModeKey === LINEAR_COMPARISON_LOSAT_MODE_KEYS.LOSATP
    ? `${losatModeLabel} · ${losatpModeLabel}`
    : losatModeLabel;
  const sourceBreakdown = projectSourceBreakdown(activeEdges);
  const activePairCount = activeEdges.length;
  const dormantDraftCount = retainedDormantDraftCount(normalizedPlan, activeEdges);
  const dormantRawResultDraftCount = retainedRawResultDraftCount(normalizedPlan, activeEdges);
  const topology = projectTopology(normalizedPlan, activePairCount);
  const filterSummary = projectFilterSummary(filters);
  const sectionKeys = projectSectionKeys({
    intentKey,
    planned,
    losatModeKey,
    losatpModeKey,
    selectedTopology,
    dormantDraftCount,
    dormantRawResultDraftCount
  });
  const errorRouting = projectErrorRouting(resolution?.errors);

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
    filterSummary: filterSummary.text,
    filterSummaryIsDefault: filterSummary.isDefault,
    retainedDormantDraftCount: dormantDraftCount,
    retainedRawResultDraftCount: dormantRawResultDraftCount,
    sectionKeys,
    errorBadge: errorRouting.badge,
    errorDisclosureKey: errorRouting.primaryDisclosureKey,
    errorDisclosureKeys: errorRouting.disclosureKeys,
    errorTargets: errorRouting.targets,
    currentStatusLabel: projectStatusLabel({ intentKey, activePairCount, sourceBreakdown }),
    summaryText: projectSummaryText({
      intentKey,
      topology,
      planned,
      sourceBreakdown,
      losatMethodSummaryLabel,
      filterSummary
    })
  });
};
