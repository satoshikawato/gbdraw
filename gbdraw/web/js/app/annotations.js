// @ts-check
/**
 * @import { AnnotationRecordCatalog } from './annotations/record-catalog.js'
 * @import { AnnotationSet } from './annotations/state.js'
 */
import { createAnnotationSet, createDefaultAnnotationStyle, normalizeAnnotationSets, uniqueAnnotationSetId } from './annotations/state.js';
import { coordinateTarget, featureTarget, featureTargetsFromSelection } from './annotations/target-actions.js';
import { encodeAnnotationTableWithNotice, parseAnnotationTableWithNotice } from './annotations/table-codec.js';
import {
  createAnnotationRecordSelector,
  reconcileAnnotationRecordBindings
} from './annotations/record-selector.js';
import { getFeatureCaption, getFeatureColorRuleHash } from '../services/feature-utils.js';
import { readFileText } from '../services/file-content-cache.js';
import { downloadTextFile } from '../services/text-download.js';

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
 *   annotationSets: AnnotationSet[], adv: { circular_track_slots?: any[], linear_track_slots?: any[] },
 *   extractedFeatures: any, biologicalFeatures: any, selectedFeatures: any, results: any,
 *   sessionOperationAvailability?: () => any
 * }} AnnotationEditorState `sessionOperationAvailability` returns a busy outcome while a Session operation runs.
 * @typedef {object} AnnotationEditorOptions
 * @property {AnnotationEditorState} state
 * @property {() => AnnotationRecordCatalog | null | undefined} getRecordCatalog
 * @property {(notice: string) => void} onImportNotice Shows or clears the notice of the last import.
 */

/**
 * @param {AnnotationEditorOptions} options
 */
export const createAnnotationEditor = ({ state, getRecordCatalog, onImportNotice }) => {
  const recordSelector = createAnnotationRecordSelector({ getCatalog: getRecordCatalog });
  // The catalog feature of a selected-feature target in the current Results:
  // a drawn feature, or (unless `drawnOnly`) a listed feature they hide, of the
  // target's mode (R2). The lists are replaced on each Generate, so scanning
  // their raw rows tracks only the lists.
  const catalogFeature = (target, { drawnOnly = false } = {}) => {
    const raw = (value) => (globalThis.window?.Vue?.toRaw ? window.Vue.toRaw(value) : value);
    const lists = [state.extractedFeatures?.value, drawnOnly ? null : state.biologicalFeatures?.value];
    for (const list of lists.map(raw)) {
      const feature = (Array.isArray(list) ? list : []).find((item) => item?.scope === target?.scope
        && item?.record_key === target?.recordKey && item?.biological_feature_id === target?.biologicalFeatureId);
      if (feature) return feature;
    }
    return null;
  };
  const featureTargetCaption = (item) => {
    const feature = catalogFeature(item?.target);
    return feature ? getFeatureCaption(feature) : `${item?.target?.biologicalFeatureId} (not in the current diagram)`;
  };
  // Design Q4 6.4: the drawn record position and drawn hash of a selected
  // feature (catalog 3 and 4 features carry the drawn hash in their rendered ID).
  const drawnPlacement = (target) => {
    const feature = catalogFeature(target, { drawnOnly: true });
    const hash = feature ? feature.drawnSelector?.hash || getFeatureColorRuleHash(feature) : '';
    const recordIndex = Number(feature?.record_idx);
    return hash && Number.isSafeInteger(recordIndex) && recordIndex >= 0 ? { recordIndex, hash } : null;
  };
  const reconcileRecords = (sets = state.annotationSets) => {
    reconcileAnnotationRecordBindings(sets, getRecordCatalog?.());
    return sets;
  };
  const replaceSets = (sets) => {
    const candidate = reconcileRecords(normalizeAnnotationSets(sets));
    state.annotationSets.splice(0, state.annotationSets.length, ...candidate);
  };
  const addAnnotationSet = (base = 'annotations') => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const set = createAnnotationSet({ id: uniqueAnnotationSetId(state.annotationSets, base) });
    state.annotationSets.push(set);
    return set;
  };
  const renameAnnotationSet = (set, id) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const oldId = String(set?.id || '');
    const nextId = uniqueAnnotationSetId(state.annotationSets.filter((item) => item !== set), id);
    set.id = nextId;
    [state.adv.circular_track_slots, state.adv.linear_track_slots].forEach((slots) => (
      (Array.isArray(slots) ? slots : []).forEach((slot) => {
        if (slot?.renderer === 'annotations' && slot?.params?.set_id === oldId) slot.params.set_id = nextId;
      })
    ));
    return nextId;
  };
  const duplicateAnnotationSet = (set) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const copy = createAnnotationSet(JSON.parse(JSON.stringify(set)));
    copy.id = uniqueAnnotationSetId(state.annotationSets, `${set.id}_copy`);
    state.annotationSets.push(copy);
    return copy;
  };
  const removeAnnotationSet = (set) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const index = state.annotationSets.indexOf(set);
    if (index >= 0) state.annotationSets.splice(index, 1);
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
    if (index >= 0) set.annotations.splice(index, 1);
  };
  const renameAnnotation = (set, item, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const base = String(value || '').trim() || 'region';
    const used = new Set(set.annotations.filter((entry) => entry !== item).map((entry) => entry.id));
    let id = base;
    for (let index = 2; used.has(id); index += 1) id = `${base}_${index}`;
    item.id = id;
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
  const importAnnotationTable = (text) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    onImportNotice?.('');
    const { sets, notice } = parseAnnotationTableWithNotice(text);
    replaceSets(sets);
    onImportNotice?.(notice);
  };
  let importSequence = 0;
  const canDownloadAnnotationTable = () => state.annotationSets.some((set) => set.annotations.length > 0);
  const downloadAnnotationTable = () => {
    if (!canDownloadAnnotationTable()) return;
    const { text, placedFeatureIdentityCount: placed, skippedFeatureIdentityCount: skipped } = (
      encodeAnnotationTableWithNotice(state.annotationSets, { drawnPlacement })
    );
    downloadTextFile('annotations.tsv', text);
    const notice = [
      placed > 0 ? `${placed} annotation(s) of selected features were written as record=#<position> and feature_selector=hash=<drawn hash>. They name the same features only while the crop, orientation, and record order stay as drawn.` : '',
      skipped > 0 ? `${skipped} annotation(s) of selected features that the current diagram does not draw were not exported.` : ''
    ].filter(Boolean).join(' ');
    if (notice) window.alert(notice);
  };
  const importAnnotationTableFile = async (event) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;
    const sequence = ++importSequence;
    onImportNotice?.('');
    const setsBeforeRead = [...state.annotationSets];
    const draftBeforeRead = JSON.stringify(state.annotationSets);
    const resultsBeforeRead = [...(state.results?.value ?? state.results ?? [])];
    const isCurrent = () => {
      const currentResults = state.results?.value ?? state.results ?? [];
      return sequence === importSequence && input.files?.[0] === file
      && setsBeforeRead.length === state.annotationSets.length
      && setsBeforeRead.every((set, index) => set === state.annotationSets[index])
      && draftBeforeRead === JSON.stringify(state.annotationSets)
      && resultsBeforeRead.length === currentResults.length
      && resultsBeforeRead.every((result, index) => result === currentResults[index]);
    };
    try {
      const text = await readFileText(file);
      if (!isCurrent()) return;
      return importAnnotationTable(text);
    } catch (error) {
      if (isCurrent()) alert(`Could not import annotations: ${error.message}`);
    } finally {
      if (sequence === importSequence && input.files?.[0] === file) input.value = '';
    }
  };
  return {
    addAnnotationSet, renameAnnotationSet, duplicateAnnotationSet, removeAnnotationSet,
    addCoordinateAnnotation, addSelectedFeatures, removeAnnotation, renameAnnotation, setAnnotationStyle, setAnnotationTargetKind,
    importAnnotationTable, importAnnotationTableFile, replaceAnnotationSets: replaceSets,
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
