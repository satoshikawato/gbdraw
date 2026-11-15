// @ts-check
/** @import { AskChoice } from './app-messages.js' */
/** @import { Wording } from '../utils/i18n.js' */
/** @import { ErrorTextPart } from '../utils/error-normalization.js' */
/** @import { WebState } from '../state.js' */
/**
 * @import { AnnotationRecordCatalog } from './annotations/record-catalog.js'
 * @import { AnnotationItem, AnnotationSet, AnnotationStyle, AnnotationTarget } from '../services/annotation-state.js'
 * @import { AnnotationRow, AnnotationSetRow, NumberInputValue } from '../services/drawing/annotation-sets.js'
 * @import { AnnotationSetsDrawing } from './annotations/record-selector.js'
 * @import { ResultRecordBinding } from '../services/sequence-search.js'
 */
import { createAnnotationSet, createDefaultAnnotationStyle, normalizeAnnotationSets, uniqueAnnotationSetId } from '../services/annotation-state.js';
import { nextCopyId } from '../services/track-slot-validation.js';
import { coordinateTarget, featureTarget, featureTargetsFromSelection } from './annotations/target-actions.js';
import { AnnotationTableRefusal, encodeAnnotationTableWithNotice, parseAnnotationTableWithNotice } from './annotations/table-codec.js';
import {
  boundAnnotation,
  createAnnotationRecordSelector,
  reconcileAnnotationRecordBindings
} from './annotations/record-selector.js';
import { getFeatureCaption, getFeatureColorRuleHash } from '../services/feature-utils.js';
import {
  ANNOTATION_SETS_SLICE, addAnnotationSet as addSetToDrawing, addAnnotations, removeAnnotation as removeAnnotationRow,
  removeAnnotationSet as removeSetFromDrawing, renameAnnotationSet as writeSetId, replaceAnnotationTarget, setAnnotationFields,
  setAnnotationSelectorField, setAnnotationSetLegendLabel as writeSetLegendLabel, setAnnotationStyleField, setAnnotationTargetBound
} from '../services/drawing/annotation-sets.js';
import { readFileText } from '../services/file-content-cache.js';
import { downloadTextFile } from '../services/text-download.js';
import { normalizeCaughtError, normalizeUserFacingError } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { toRaw } from '../utils/vue-raw.js';
import { joinSentences, t, tn } from '../utils/i18n.js';
import { hashSelector, recordIndexSelector } from '../utils/selector-syntax.js';

/**
 * The text of an annotation TSV. A failed read rejects as the panel notice tells it: a refusal the diagnostic model
 * knows as an `AnnotationTableRefusal`, in its words, which name the file (OV-2948); any other error as it is, which
 * the notice normalizes (FL-14). A falsy rejection has no words.
 * @param {unknown} file
 * @returns {Promise<string>}
 */
const readAnnotationTableText = (file) => readFileText(file).catch((error) => {
  const shown = normalizeUserFacingError(/** @type {Record<string, unknown> | null} */ (error || null), { stage: 'read' });
  throw shown?.code === 'UNKNOWN' && error instanceof Error ? error : new AnnotationTableRefusal(shown?.summaryParts ?? []);
});

/**
 * The panel notice of a failed import: a refused read or table in its own words, any other failure in its normalized
 * parts, never its raw text; worded where the panel shows it (OV-3683).
 * @param {unknown} error
 * @returns {Wording}
 */
const importFailureNotice = (error) => {
  /** @type {readonly ErrorTextPart[]} */
  const parts = error instanceof AnnotationTableRefusal ? error.parts : normalizeCaughtError(error, { stage: 'parse' }).summaryParts;
  return () => t('Could not import annotations: {reason}', { reason: errorText(parts) });
};

/** The panel notice of Add selected features on a Result without feature identities. @type {Wording} */
const REGENERATE_TO_ANNOTATE_NOTICE = () => t('Generate the diagram again to annotate the selected features.');

/**
 * The panel notice of Download TSV: how the selected-feature rows were written, and the ones it left out.
 * @param {number} placed
 * @param {number} skipped
 * @returns {Wording | null}
 */
const tableExportNotice = (placed, skipped) => (placed > 0 || skipped > 0 ? () => joinSentences(...[
  placed > 0 ? tn(placed,
    '{count} annotation of a selected feature was written as {recordSelector} and {featureSelector}. It names the same feature only while the record order stays as drawn.',
    '{count} annotations of selected features were written as {recordSelector} and {featureSelector}. They name the same features only while the record order stays as drawn.',
    { recordSelector: `record=${recordIndexSelector('<position>')}`, featureSelector: `feature_selector=${hashSelector('<feature hash>')}` }) : '',
  skipped > 0 ? tn(skipped,
    '{count} annotation of a selected feature that the current diagram does not draw was not exported.',
    '{count} annotations of selected features that the current diagram does not draw were not exported.') : ''
].filter(Boolean)) : null);

// The editor keeps every id valid: an empty id takes the default and a used
// one a numbered suffix. The panel notice names the id it stored (FL-14).
/**
 * @param {'set' | 'annotation'} kind @param {any} requested @param {string} stored
 * @returns {Wording | null} Null when the stored id is the requested one: the edit clears the notice.
 */
const idEditNotice = (kind, requested, stored) => {
  const text = String(requested ?? '').trim();
  if (text === stored) return null;
  const params = { text, stored };
  if (kind === 'set') {
    return text
      ? () => t('Set id "{text}" is used by another set, so it was saved as "{stored}".', params)
      : () => t('Set id cannot be empty, so it was saved as "{stored}".', params);
  }
  return text
    ? () => t('Annotation id "{text}" is used by another annotation of this set, so it was saved as "{stored}".', params)
    : () => t('Annotation id cannot be empty, so it was saved as "{stored}".', params);
};

/**
 * @param {readonly AnnotationRow[]} annotations
 * @param {string} prefix
 */
const nextAnnotationId = (annotations, prefix) => {
  const ids = new Set(annotations.map((item) => item.id));
  let index = annotations.length + 1;
  while (ids.has(`${prefix}_${index}`)) index += 1;
  return `${prefix}_${index}`;
};

/**
 * The Web state the editor reads. `sessionOperationAvailability` returns a busy outcome while a Session operation
 * runs.
 * @typedef {Pick<WebState, 'extractedFeatures' | 'biologicalFeatures' | 'selectedFeatures' | 'results'
 *   | 'sessionOperationAvailability'> & {
 *   activeDrawing: () => AnnotationEditorDrawing
 * }} AnnotationEditorState
 * @typedef {AnnotationSetsDrawing & {
 *   adv: { circular_track_slots?: any[], linear_track_slots?: any[] }
 * }} AnnotationEditorDrawing The members of a drawing (`DrawingState` of state.js) the editor reads, and the
 *   handle the set list's writers take.
 * @typedef {object} AnnotationEditorOptions
 * @property {AnnotationEditorState} state
 * @property {() => AnnotationRecordCatalog | null | undefined} getRecordCatalog
 * @property {(notice: Wording | null) => void} onImportNotice Shows the panel notice, worded where it shows, or
 *   clears it (null): the outcome of the last import, export or id edit.
 * @property {<T>(change: () => T) => T} [retireLegendStylesOfUnnamedCaptions] Runs a track data change and
 *   retires the Legend styles of the captions it no longer names (OV-65).
 * @property {AskChoice} [askChoice] The root's choice dialog: an import that removes every set asks first.
 */

/**
 * @param {AnnotationEditorOptions} options
 */
export const createAnnotationEditor = ({
  state, getRecordCatalog, onImportNotice, retireLegendStylesOfUnnamedCaptions = (change) => change(),
  askChoice = async () => 'cancel'
}) => {
  const recordSelector = createAnnotationRecordSelector({ getCatalog: getRecordCatalog });
  // The catalog feature of a selected-feature target in the current Results:
  // a drawn feature, or (unless `drawnOnly`) a listed feature they hide. The
  // shown Results are of the drawing's mode (R2). The lists are replaced on
  // each Generate, so scanning their raw rows tracks only the lists.
  /**
   * @param {Pick<AnnotationTarget, 'recordKey' | 'biologicalFeatureId'> | null | undefined} target
   * @param {{ drawnOnly?: boolean }} [options]
   */
  const catalogFeature = (target, { drawnOnly = false } = {}) => {
    const lists = [state.extractedFeatures.value, drawnOnly ? null : state.biologicalFeatures.value];
    for (const list of lists.map(toRaw)) {
      const feature = (Array.isArray(list) ? list : []).find((item) => item?.record_key === target?.recordKey && item?.biological_feature_id === target?.biologicalFeatureId);
      if (feature) return feature;
    }
    return null;
  };
  /** @param {AnnotationItem | null | undefined} item */
  const featureTargetCaption = (item) => {
    const feature = catalogFeature(item?.target);
    return feature ? getFeatureCaption(feature) : t('{id} (not in the current diagram)', { id: String(item?.target?.biologicalFeatureId) });
  };
  // Design Q4 6.4: the drawn record position and the feature hash of a
  // selected feature, which a `hash` selector matches (OV-401).
  /** @param {AnnotationTarget} target */
  const drawnPlacement = (target) => {
    const feature = catalogFeature(target, { drawnOnly: true });
    const hash = feature ? getFeatureColorRuleHash(feature) : '';
    const recordIndex = Number(feature?.record_idx);
    return hash && Number.isSafeInteger(recordIndex) && recordIndex >= 0 ? { recordIndex, hash } : null;
  };
  // The catalog is read first: when it cannot be read, the sets stay as they were.
  /** @param {Partial<AnnotationSet>[] | null | undefined} sets */
  const replaceSets = (sets) => {
    const drawing = state.activeDrawing();
    const catalog = getRecordCatalog?.();
    retireLegendStylesOfUnnamedCaptions(() => {
      ANNOTATION_SETS_SLICE.install(drawing.slices.annotationSets, { annotationSets: normalizeAnnotationSets(sets) });
      reconcileAnnotationRecordBindings(drawing, catalog);
    });
  };
  const addAnnotationSet = (base = 'annotations') => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const set = createAnnotationSet({ id: uniqueAnnotationSetId(drawing.annotationSets, base) });
    return addSetToDrawing(drawing.slices.annotationSets, set);
  };
  /**
   * @param {AnnotationSetRow} set
   * @param {string} id
   */
  const renameAnnotationSet = (set, id) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const oldId = String(set?.id || '');
    const nextId = uniqueAnnotationSetId(drawing.annotationSets.filter((item) => item !== set), id);
    writeSetId(drawing.slices.annotationSets, set, nextId);
    onImportNotice?.(idEditNotice('set', id, nextId));
    [drawing.adv.circular_track_slots, drawing.adv.linear_track_slots].forEach((slots) => (
      (Array.isArray(slots) ? slots : []).forEach((slot) => {
        if (slot?.renderer === 'annotations' && slot?.params?.set_id === oldId) slot.params.set_id = nextId;
      })
    ));
    return nextId;
  };
  // The legend label of a set names its rows: a new text retires the Legend
  // styles of the old caption in the same step (OV-67).
  /**
   * @param {AnnotationSetRow | null | undefined} set
   * @param {unknown} value
   */
  const setAnnotationSetLegendLabel = (set, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!set) return;
    retireLegendStylesOfUnnamedCaptions(() => writeSetLegendLabel(drawing.slices.annotationSets, set, String(value ?? '').trim()));
  };
  /** @param {AnnotationSetRow} set */
  const duplicateAnnotationSet = (set) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const copy = createAnnotationSet(JSON.parse(JSON.stringify(set)));
    copy.id = nextCopyId(set.id, drawing.annotationSets);
    return addSetToDrawing(drawing.slices.annotationSets, copy);
  };
  /** @param {AnnotationSetRow} set */
  const removeAnnotationSet = (set) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (drawing.annotationSets.includes(set)) {
      retireLegendStylesOfUnnamedCaptions(() => removeSetFromDrawing(drawing.slices.annotationSets, set));
    }
  };
  /**
   * @param {AnnotationSetRow | null | undefined} set
   * @param {{ start?: unknown, end?: unknown, recordId?: string | null, recordIndex?: number | null, coordinateSpace?: string }} [options]
   */
  const addCoordinateAnnotation = (set, options = {}) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!set) return null;
    const id = nextAnnotationId(set.annotations, 'region');
    const item = { id, target: coordinateTarget({ start: 1, end: 1, ...options }), label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
    return addAnnotations(drawing.slices.annotationSets, set, [item])[0] ?? null;
  };
  /**
   * Select region's Add to annotation set: one coordinate-span row labelled "Region N" (the next number no row of the
   * set uses) in the default style, in the drawing's first annotation set, which it adds when there is none. The row
   * is bound to its record of the Result as a search hit is (`record`, `resultRecordBinding`): a record whose card was
   * deleted after Generate leaves its target unavailable until the row is deleted or rebound (OV-2632, D-58).
   * @param {{ record: ResultRecordBinding, start: number, end: number }} region Source coordinates.
   */
  const addRegionAnnotation = ({ record, start, end }) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const set = drawing.annotationSets[0]
      ?? addSetToDrawing(drawing.slices.annotationSets, createAnnotationSet({ id: uniqueAnnotationSetId(drawing.annotationSets, 'annotations') }));
    const used = new Set(set.annotations.map((item) => /^Region (\d+)$/.exec(String(item.label ?? ''))?.[1]).filter(Boolean).map(Number));
    let number = 1;
    while (used.has(number)) number += 1;
    const item = { id: nextAnnotationId(set.annotations, 'region'), target: coordinateTarget({ start, end }), label: `Region ${number}`, mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
    return addAnnotations(drawing.slices.annotationSets, set, [boundAnnotation(item, record)])[0] ?? null;
  };
  /** @param {AnnotationSetRow | null | undefined} set */
  const addSelectedFeatures = (set) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!set) return [];
    const targets = featureTargetsFromSelection(state.selectedFeatures.value);
    if (!targets) {
      onImportNotice?.(REGENERATE_TO_ANNOTATE_NOTICE);
      return [];
    }
    // Each new row's id counts the rows before it.
    /** @type {AnnotationItem[]} */
    const items = [];
    targets.forEach((target) => {
      items.push({ id: nextAnnotationId([...set.annotations, ...items], 'feature'), target, label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} });
    });
    return addAnnotations(drawing.slices.annotationSets, set, items);
  };
  /**
   * @param {AnnotationSetRow} set
   * @param {AnnotationRow} item
   */
  const removeAnnotation = (set, item) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (set?.annotations?.includes(item)) {
      retireLegendStylesOfUnnamedCaptions(() => removeAnnotationRow(drawing.slices.annotationSets, set, item));
    }
  };
  /**
   * @param {AnnotationSetRow} set
   * @param {AnnotationRow} item
   * @param {unknown} value
   */
  const renameAnnotation = (set, item, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const base = String(value || '').trim() || 'region';
    const used = new Set(set.annotations.filter((entry) => entry !== item).map((entry) => entry.id));
    let id = base;
    for (let index = 2; used.has(id); index += 1) id = `${base}_${index}`;
    setAnnotationFields(drawing.slices.annotationSets, set, item, { id });
    onImportNotice?.(idEditNotice('annotation', value, id));
  };
  /**
   * @template {keyof AnnotationStyle} K
   * @param {AnnotationSetRow} set
   * @param {AnnotationRow} item
   * @param {K} field
   * @param {AnnotationStyle[K]} value
   */
  const setAnnotationStyle = (set, item, field, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const defaults = createDefaultAnnotationStyle(set.defaultStyle);
    setAnnotationStyleField(drawing.slices.annotationSets, set, item, field, value, defaults);
  };
  /**
   * @param {AnnotationSetRow | null | undefined} set
   * @param {AnnotationRow | null | undefined} item
   * @param {string} kind
   */
  const setAnnotationTargetKind = (set, item, kind) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    if (!set || !item) return;
    const record = item.target?.record ?? null;
    /** @type {AnnotationTarget} */
    const target = kind === 'featureSpan'
      ? featureTarget({ selector: 'locus_tag=' })
      : coordinateTarget({ start: 1, end: 1 });
    target.record = record && { ...record };
    replaceAnnotationTarget(drawing.slices.annotationSets, set, item, target);
  };
  /**
   * The text of a number input as Vue's `v-model.number` writes it: a number,
   * or the text when it reads as none.
   * @param {unknown} text
   * @returns {NumberInputValue}
   */
  const numberInputValue = (text) => {
    const number = Number.parseFloat(String(text));
    return Number.isNaN(number) ? String(text) : number;
  };
  /**
   * The template's field edits of an annotation row (in place of `v-model`).
   * @param {AnnotationSetRow} set
   * @param {AnnotationRow} item
   * @param {'start' | 'end' | 'selectorKey' | 'selectorValue' | 'mark' | 'label' | 'lane'} field
   * @param {unknown} value The input's value.
   */
  const setAnnotationField = (set, item, field, value) => {
    const drawing = state.activeDrawing();
    const slice = drawing.slices.annotationSets;
    if (field === 'start' || field === 'end') setAnnotationTargetBound(slice, set, item, field, numberInputValue(value));
    else if (field === 'selectorKey') setAnnotationSelectorField(slice, set, item, 'key', String(value).trim());
    else if (field === 'selectorValue') setAnnotationSelectorField(slice, set, item, 'value', String(value).trim());
    else if (field === 'lane') setAnnotationFields(slice, set, item, { lane: numberInputValue(value) });
    else setAnnotationFields(slice, set, item, { [field]: String(value) });
  };
  // An import replaces every set. A table without rows would remove them all,
  // so it asks first (FL-14); the caller's History step stays open until the
  // choice, and Cancel changes nothing.
  /** @param {number} setCount */
  const askReplaceAll = (setCount) => askChoice({
    title: t('Remove all annotation sets?'),
    description: joinSentences(t('The imported table has no annotation rows.'),
      tn(setCount, 'Importing it removes the annotation set.', 'Importing it removes the {count} annotation sets.')),
    actions: [{ id: 'remove', label: t('Remove all sets'), tone: 'danger' }]
  });
  /** @param {string} text */
  const importAnnotationTable = async (text) => {
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    onImportNotice?.(null);
    const { sets, notice } = parseAnnotationTableWithNotice(text);
    const setCount = state.activeDrawing().annotationSets.length;
    if (sets.length === 0 && setCount > 0 && await askReplaceAll(setCount) !== 'remove') return;
    replaceSets(sets);
    onImportNotice?.(notice);
  };
  let importSequence = 0;
  const canDownloadAnnotationTable = () => state.activeDrawing().annotationSets.some((set) => set.annotations.length > 0);
  const downloadAnnotationTable = () => {
    const drawing = state.activeDrawing();
    if (!canDownloadAnnotationTable()) return;
    const { text, placedFeatureIdentityCount: placed, skippedFeatureIdentityCount: skipped } = (
      encodeAnnotationTableWithNotice(drawing.annotationSets, { drawnPlacement })
    );
    downloadTextFile('annotations.tsv', text);
    onImportNotice?.(tableExportNotice(placed, skipped));
  };
  /** @param {{ target?: HTMLInputElement | null } | null | undefined} event The change event of the file input. */
  const importAnnotationTableFile = async (event) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;
    const sequence = ++importSequence;
    onImportNotice?.(null);
    const setsBeforeRead = [...drawing.annotationSets];
    const draftBeforeRead = JSON.stringify(drawing.annotationSets);
    const resultsBeforeRead = [...state.results.value];
    const isCurrent = () => {
      const currentResults = state.results.value;
      return sequence === importSequence && input.files?.[0] === file
      && setsBeforeRead.length === drawing.annotationSets.length
      && setsBeforeRead.every((set, index) => set === drawing.annotationSets[index])
      && draftBeforeRead === JSON.stringify(drawing.annotationSets)
      && resultsBeforeRead.length === currentResults.length
      && resultsBeforeRead.every((result, index) => result === currentResults[index]);
    };
    try {
      const text = await readAnnotationTableText(file);
      if (!isCurrent()) return;
      return await importAnnotationTable(text);
    } catch (error) {
      if (isCurrent()) onImportNotice?.(importFailureNotice(error));
    } finally {
      if (sequence === importSequence && input.files?.[0] === file) input.value = '';
    }
  };
  return {
    addAnnotationSet, renameAnnotationSet, setAnnotationSetLegendLabel, duplicateAnnotationSet, removeAnnotationSet,
    addCoordinateAnnotation, addRegionAnnotation, addSelectedFeatures, removeAnnotation, renameAnnotation, setAnnotationStyle, setAnnotationTargetKind,
    setAnnotationField,
    importAnnotationTable, importAnnotationTableFile, replaceAnnotationSets: replaceSets,
    canDownloadAnnotationTable, downloadAnnotationTable, featureTargetCaption,
    recordOptionsFor: recordSelector.optionsFor,
    recordValueFor: recordSelector.valueFor,
    /**
     * @param {AnnotationRow | null | undefined} annotation
     * @param {unknown} value
     */
    setRecordValue: (annotation, value) => state.sessionOperationAvailability()
      || recordSelector.setValue(state.activeDrawing().slices.annotationSets, annotation, value),
    recordIsRequired: recordSelector.isRequired,
    recordIsMissing: recordSelector.isMissing,
    recordMissingMessage: recordSelector.missingMessage,
    recordShowsSelect: recordSelector.showsSelect,
    recordCaptionFor: recordSelector.captionFor
  };
};
