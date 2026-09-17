import { createAnnotationSet, createDefaultAnnotationStyle, normalizeAnnotationSets, uniqueAnnotationSetId } from './annotations/state.js';
import { coordinateTarget, featureTarget, featureTargetsFromSelection } from './annotations/target-actions.js';
import { encodeAnnotationTable, parseAnnotationTable } from './annotations/table-codec.js';
import {
  createAnnotationRecordSelector,
  reconcileAnnotationRecordBindings
} from './annotations/record-selector.js';
import { readFileText } from '../services/file-content-cache.js';
import { downloadTextFile } from '../services/text-download.js';

const nextAnnotationId = (annotations, prefix) => {
  const ids = new Set(annotations.map((item) => item.id));
  let index = annotations.length + 1;
  while (ids.has(`${prefix}_${index}`)) index += 1;
  return `${prefix}_${index}`;
};

export const createAnnotationEditor = ({ state, getRecordCatalog }) => {
  const recordSelector = createAnnotationRecordSelector({ getCatalog: getRecordCatalog });
  const reconcileRecords = (sets = state.annotationSets) => {
    reconcileAnnotationRecordBindings(sets, getRecordCatalog?.());
    return sets;
  };
  const replaceSets = (sets) => {
    state.annotationSets.splice(0, state.annotationSets.length, ...normalizeAnnotationSets(sets));
    reconcileRecords();
  };
  const addAnnotationSet = (base = 'annotations') => {
    const set = createAnnotationSet({ id: uniqueAnnotationSetId(state.annotationSets, base) });
    state.annotationSets.push(set);
    return set;
  };
  const renameAnnotationSet = (set, id) => {
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
    const copy = createAnnotationSet(JSON.parse(JSON.stringify(set)));
    copy.id = uniqueAnnotationSetId(state.annotationSets, `${set.id}_copy`);
    state.annotationSets.push(copy);
    return copy;
  };
  const removeAnnotationSet = (set) => {
    const index = state.annotationSets.indexOf(set);
    if (index >= 0) state.annotationSets.splice(index, 1);
  };
  const addCoordinateAnnotation = (set, options = {}) => {
    if (!set) return null;
    const id = nextAnnotationId(set.annotations, 'region');
    const item = { id, target: coordinateTarget({ start: 1, end: 1, ...options }), label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
    set.annotations.push(item);
    return item;
  };
  const addSelectedFeatures = (set) => {
    if (!set) return [];
    const targets = featureTargetsFromSelection(state.selectedFeatures?.value ?? state.selectedFeatures ?? []);
    const items = targets.map((target) => {
      const item = { id: nextAnnotationId(set.annotations, 'feature'), target, label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
      set.annotations.push(item);
      return item;
    });
    reconcileRecords(items.length ? [{ annotations: items }] : []);
    return items;
  };
  const removeAnnotation = (set, item) => {
    const index = set?.annotations?.indexOf(item) ?? -1;
    if (index >= 0) set.annotations.splice(index, 1);
  };
  const renameAnnotation = (set, item, value) => {
    const base = String(value || '').trim() || 'region';
    const used = new Set(set.annotations.filter((entry) => entry !== item).map((entry) => entry.id));
    let id = base;
    for (let index = 2; used.has(id); index += 1) id = `${base}_${index}`;
    item.id = id;
  };
  const setAnnotationStyle = (set, item, field, value) => {
    item.style ??= createDefaultAnnotationStyle(set.defaultStyle);
    item.style[field] = value;
  };
  const setAnnotationTargetKind = (item, kind) => {
    if (!item) return;
    const record = item.target?.record ?? null;
    item.target = kind === 'featureSpan'
      ? featureTarget({ selector: 'locus_tag=' })
      : coordinateTarget({ start: 1, end: 1 });
    item.target.record = record;
  };
  const importAnnotationTable = (text) => replaceSets(parseAnnotationTable(text));
  const canDownloadAnnotationTable = () => state.annotationSets.some((set) => set.annotations.length > 0);
  const downloadAnnotationTable = () => {
    if (!canDownloadAnnotationTable()) return;
    downloadTextFile('annotations.tsv', encodeAnnotationTable(state.annotationSets));
  };
  const importAnnotationTableFile = async (event) => {
    const input = event?.target;
    const file = input?.files?.[0];
    if (!file) return;
    const setsBeforeRead = [...state.annotationSets];
    try {
      const text = await readFileText(file);
      if (input.files?.[0] !== file
        || setsBeforeRead.length !== state.annotationSets.length
        || setsBeforeRead.some((set, index) => set !== state.annotationSets[index])) return;
      importAnnotationTable(text);
    } catch (error) {
      if (input.files?.[0] === file) alert(`Could not import annotations: ${error.message}`);
    } finally {
      if (input.files?.[0] === file) input.value = '';
    }
  };
  return {
    addAnnotationSet, renameAnnotationSet, duplicateAnnotationSet, removeAnnotationSet,
    addCoordinateAnnotation, addSelectedFeatures, removeAnnotation, renameAnnotation, setAnnotationStyle, setAnnotationTargetKind,
    importAnnotationTable, importAnnotationTableFile, replaceAnnotationSets: replaceSets,
    canDownloadAnnotationTable, downloadAnnotationTable,
    recordOptionsFor: recordSelector.optionsFor,
    recordValueFor: recordSelector.valueFor,
    setRecordValue: recordSelector.setValue,
    recordIsRequired: recordSelector.isRequired,
    recordIsMissing: recordSelector.isMissing,
    recordIsDisabled: recordSelector.isDisabled,
    recordMissingMessage: recordSelector.missingMessage,
    reconcileRecordBindings: reconcileRecords
  };
};
