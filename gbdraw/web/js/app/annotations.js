// @ts-check
/**
 * @import { AnnotationRecordCatalog } from './annotations/record-catalog.js'
 * @import { AnnotationSet } from '../services/annotation-state.js'
 */
import { createAnnotationSet, createDefaultAnnotationStyle, normalizeAnnotationSets, uniqueAnnotationSetId } from '../services/annotation-state.js';
import { coordinateTarget, featureTarget, featureTargetsFromSelection } from './annotations/target-actions.js';
import { encodeAnnotationTableWithNotice, parseAnnotationTableWithNotice } from './annotations/table-codec.js';
import {
  createAnnotationRecordSelector,
  reconcileAnnotationRecordBindings
} from './annotations/record-selector.js';
import { getFeatureCaption, getFeatureColorRuleHash } from '../services/feature-utils.js';
import { readFileText } from '../services/file-content-cache.js';
import { downloadTextFile } from '../services/text-download.js';

// The editor keeps every id valid: an empty id takes the default and a used
// one a numbered suffix. The panel notice names the id it stored (FL-14).
/** @param {string} field @param {any} requested @param {string} stored @param {string} owner */
const idEditNotice = (field, requested, stored, owner) => {
  const text = String(requested ?? '').trim();
  if (text === stored) return '';
  return text
    ? `${field} "${text}" is used by ${owner}, so it was saved as "${stored}".`
    : `${field} cannot be empty, so it was saved as "${stored}".`;
};

const nextAnnotationId = (annotations, prefix) => {
  const ids = new Set(annotations.map((item) => item.id));
  let index = annotations.length + 1;
  while (ids.has(`${prefix}_${index}`)) index += 1;
  return `${prefix}_${index}`;
};

/**
 * The Web state the editor reads; the four feature and result members are refs
 * or plain lists, read only as the Results hold them.
 * @typedef {{
 *   extractedFeatures: any, biologicalFeatures: any, selectedFeatures: any, results: any,
 *   activeDrawing: () => AnnotationEditorDrawing, sessionOperationAvailability?: () => any
 * }} AnnotationEditorState `sessionOperationAvailability` returns a busy outcome while a Session operation runs.
 * @typedef {{
 *   annotationSets: AnnotationSet[], adv: { circular_track_slots?: any[], linear_track_slots?: any[] }
 * }} AnnotationEditorDrawing The members of a drawing (`DrawingState` of state.js) the editor reads.
 * @typedef {object} AnnotationEditorOptions
 * @property {AnnotationEditorState} state
 * @property {() => AnnotationRecordCatalog | null | undefined} getRecordCatalog
 * @property {(notice: string) => void} onImportNotice Shows or clears the panel notice: the outcome of
 *   the last import or id edit.
 * @property {<T>(change: () => T) => T} [retireLegendStylesOfUnnamedCaptions] Runs a track data change and
 *   retires the Legend styles of the captions it no longer names (OV-65).
 * @property {<T extends object>(value: T) => T} [reactive] Makes the dialog state reactive (Vue's `reactive`).
 */

/**
 * @param {AnnotationEditorOptions} options
 */
export const createAnnotationEditor = ({
  state, getRecordCatalog, onImportNotice, retireLegendStylesOfUnnamedCaptions = (change) => change(),
  reactive = (value) => value
}) => {
  const recordSelector = createAnnotationRecordSelector({ getCatalog: getRecordCatalog });
  // The catalog feature of a selected-feature target in the current Results:
  // a drawn feature, or (unless `drawnOnly`) a listed feature they hide. The
  // shown Results are of the drawing's mode (R2). The lists are replaced on
  // each Generate, so scanning their raw rows tracks only the lists.
  const catalogFeature = (target, { drawnOnly = false } = {}) => {
    const raw = (value) => (globalThis.window?.Vue?.toRaw ? window.Vue.toRaw(value) : value);
    const lists = [state.extractedFeatures?.value, drawnOnly ? null : state.biologicalFeatures?.value];
    for (const list of lists.map(raw)) {
      const feature = (Array.isArray(list) ? list : []).find((item) => item?.record_key === target?.recordKey && item?.biological_feature_id === target?.biologicalFeatureId);
      if (feature) return feature;
    }
    return null;
  };
  const featureTargetCaption = (item) => {
    const feature = catalogFeature(item?.target);
    return feature ? getFeatureCaption(feature) : `${item?.target?.biologicalFeatureId} (not in the current diagram)`;
  };
  // Design Q4 6.4: the drawn record position and the feature hash of a
  // selected feature, which a `hash` selector matches (OV-401).
  const drawnPlacement = (target) => {
    const feature = catalogFeature(target, { drawnOnly: true });
    const hash = feature ? getFeatureColorRuleHash(feature) : '';
    const recordIndex = Number(feature?.record_idx);
    return hash && Number.isSafeInteger(recordIndex) && recordIndex >= 0 ? { recordIndex, hash } : null;
  };
  const reconcileRecords = (sets = state.activeDrawing().annotationSets) => {
    reconcileAnnotationRecordBindings(sets, getRecordCatalog?.());
    return sets;
  };
  const replaceSets = (sets) => {
    const drawing = state.activeDrawing();
    const candidate = reconcileRecords(normalizeAnnotationSets(sets));
    retireLegendStylesOfUnnamedCaptions(() => {
      drawing.annotationSets.splice(0, drawing.annotationSets.length, ...candidate);
    });
  };
  const addAnnotationSet = (base = 'annotations') => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const set = createAnnotationSet({ id: uniqueAnnotationSetId(drawing.annotationSets, base) });
    drawing.annotationSets.push(set);
    return set;
  };
  const renameAnnotationSet = (set, id) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const oldId = String(set?.id || '');
    const nextId = uniqueAnnotationSetId(drawing.annotationSets.filter((item) => item !== set), id);
    set.id = nextId;
    onImportNotice?.(idEditNotice('Set id', id, nextId, 'another set'));
    [drawing.adv.circular_track_slots, drawing.adv.linear_track_slots].forEach((slots) => (
      (Array.isArray(slots) ? slots : []).forEach((slot) => {
        if (slot?.renderer === 'annotations' && slot?.params?.set_id === oldId) slot.params.set_id = nextId;
      })
    ));
    return nextId;
  };
  // The legend label of a set names its rows: a new text retires the Legend
  // styles of the old caption in the same step (OV-67).
  const setAnnotationSetLegendLabel = (set, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!set) return;
    retireLegendStylesOfUnnamedCaptions(() => { set.legendLabel = String(value ?? '').trim(); });
  };
  const duplicateAnnotationSet = (set) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const copy = createAnnotationSet(JSON.parse(JSON.stringify(set)));
    copy.id = uniqueAnnotationSetId(drawing.annotationSets, `${set.id}_copy`);
    drawing.annotationSets.push(copy);
    return copy;
  };
  const removeAnnotationSet = (set) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const index = drawing.annotationSets.indexOf(set);
    if (index >= 0) retireLegendStylesOfUnnamedCaptions(() => drawing.annotationSets.splice(index, 1));
  };
  const addCoordinateAnnotation = (set, options = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!set) return null;
    const id = nextAnnotationId(set.annotations, 'region');
    const item = { id, target: coordinateTarget({ start: 1, end: 1, ...options }), label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
    set.annotations.push(item);
    return item;
  };
  const addSelectedFeatures = (set) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!set) return [];
    const targets = featureTargetsFromSelection(state.selectedFeatures?.value ?? state.selectedFeatures ?? []);
    if (!targets) {
      window.alert('Generate the diagram again to annotate the selected features.');
      return [];
    }
    return targets.map((target) => {
      const item = { id: nextAnnotationId(set.annotations, 'feature'), target, label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
      set.annotations.push(item);
      return item;
    });
  };
  const removeAnnotation = (set, item) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const index = set?.annotations?.indexOf(item) ?? -1;
    if (index >= 0) retireLegendStylesOfUnnamedCaptions(() => set.annotations.splice(index, 1));
  };
  const renameAnnotation = (set, item, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const base = String(value || '').trim() || 'region';
    const used = new Set(set.annotations.filter((entry) => entry !== item).map((entry) => entry.id));
    let id = base;
    for (let index = 2; used.has(id); index += 1) id = `${base}_${index}`;
    item.id = id;
    onImportNotice?.(idEditNotice('Annotation id', value, id, 'another annotation of this set'));
  };
  const setAnnotationStyle = (set, item, field, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    item.style ??= createDefaultAnnotationStyle(set.defaultStyle);
    item.style[field] = value;
  };
  const setAnnotationTargetKind = (item, kind) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!item) return;
    const record = item.target?.record ?? null;
    item.target = kind === 'featureSpan'
      ? featureTarget({ selector: 'locus_tag=' })
      : coordinateTarget({ start: 1, end: 1 });
    item.target.record = record;
  };
  // An import replaces every set. A table without rows would remove them all,
  // so it asks first (FL-14); the caller's History step stays open until the
  // choice, and Cancel changes nothing.
  const replaceAllDialog = reactive({ show: false, setCount: 0 });
  /** @type {((choice: string) => void) | null} */
  let answerReplaceAllDialog = null;
  const askReplaceAll = (setCount) => {
    answerReplaceAllDialog?.('cancel');
    replaceAllDialog.setCount = setCount;
    replaceAllDialog.show = true;
    return new Promise((resolve) => { answerReplaceAllDialog = resolve; });
  };
  const handleReplaceAllChoice = (choice) => {
    const answer = answerReplaceAllDialog;
    answerReplaceAllDialog = null;
    replaceAllDialog.show = false;
    answer?.(choice);
  };
  const importAnnotationTable = async (text) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    onImportNotice?.('');
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
    const notice = [
      placed > 0 ? `${placed} annotation(s) of selected features were written as record=#<position> and feature_selector=hash=<feature hash>. They name the same features only while the record order stays as drawn.` : '',
      skipped > 0 ? `${skipped} annotation(s) of selected features that the current diagram does not draw were not exported.` : ''
    ].filter(Boolean).join(' ');
    if (notice) window.alert(notice);
  };
  const importAnnotationTableFile = async (event) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;
    const sequence = ++importSequence;
    onImportNotice?.('');
    const setsBeforeRead = [...drawing.annotationSets];
    const draftBeforeRead = JSON.stringify(drawing.annotationSets);
    const resultsBeforeRead = [...(state.results?.value ?? state.results ?? [])];
    const isCurrent = () => {
      const currentResults = state.results?.value ?? state.results ?? [];
      return sequence === importSequence && input.files?.[0] === file
      && setsBeforeRead.length === drawing.annotationSets.length
      && setsBeforeRead.every((set, index) => set === drawing.annotationSets[index])
      && draftBeforeRead === JSON.stringify(drawing.annotationSets)
      && resultsBeforeRead.length === currentResults.length
      && resultsBeforeRead.every((result, index) => result === currentResults[index]);
    };
    try {
      const text = await readFileText(file);
      if (!isCurrent()) return;
      return await importAnnotationTable(text);
    } catch (error) {
      if (isCurrent()) onImportNotice?.(`Could not import annotations: ${error.message}`);
    } finally {
      if (sequence === importSequence && input.files?.[0] === file) input.value = '';
    }
  };
  return {
    addAnnotationSet, renameAnnotationSet, setAnnotationSetLegendLabel, duplicateAnnotationSet, removeAnnotationSet,
    addCoordinateAnnotation, addSelectedFeatures, removeAnnotation, renameAnnotation, setAnnotationStyle, setAnnotationTargetKind,
    importAnnotationTable, importAnnotationTableFile, replaceAnnotationSets: replaceSets,
    annotationReplaceAllDialog: replaceAllDialog, handleAnnotationReplaceAllChoice: handleReplaceAllChoice,
    canDownloadAnnotationTable, downloadAnnotationTable, featureTargetCaption,
    recordOptionsFor: recordSelector.optionsFor,
    recordValueFor: recordSelector.valueFor,
    setRecordValue: (annotation, value) => state.sessionOperationAvailability?.()
      || recordSelector.setValue(annotation, value),
    recordIsRequired: recordSelector.isRequired,
    recordIsMissing: recordSelector.isMissing,
    recordIsDisabled: recordSelector.isDisabled,
    recordMissingMessage: recordSelector.missingMessage,
    reconcileRecordBindings: reconcileRecords
  };
};
