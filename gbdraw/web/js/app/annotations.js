import { createAnnotationSet, normalizeAnnotationSets, uniqueAnnotationSetId } from './annotations/state.js';
import { coordinateTarget, featureTarget, featureTargetsFromSelection } from './annotations/target-actions.js';
import { parseAnnotationTable } from './annotations/table-codec.js';

export const createAnnotationEditor = ({ state }) => {
  const replaceSets = (sets) => state.annotationSets.splice(0, state.annotationSets.length, ...normalizeAnnotationSets(sets));
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
    const id = `region_${set.annotations.length + 1}`;
    const item = { id, target: coordinateTarget({ start: 1, end: 1, ...options }), label: '', mark: 'bracket', lane: null, style: null, legendLabel: null, metadata: {} };
    set.annotations.push(item);
    return item;
  };
  const addSelectedFeatures = (set) => {
    if (!set) return [];
    const targets = featureTargetsFromSelection(state.selectedFeatures?.value ?? state.selectedFeatures ?? []);
    const items = targets.map((target, index) => ({ id: `feature_${set.annotations.length + index + 1}`, target, label: '', mark: 'bracket', lane: null, style: null, legendLabel: null, metadata: {} }));
    set.annotations.push(...items);
    return items;
  };
  const removeAnnotation = (set, item) => {
    const index = set?.annotations?.indexOf(item) ?? -1;
    if (index >= 0) set.annotations.splice(index, 1);
  };
  const setAnnotationTargetKind = (item, kind) => {
    if (!item) return;
    item.target = kind === 'featureSpan'
      ? featureTarget({ selector: 'locus_tag=' })
      : coordinateTarget({ start: 1, end: 1 });
  };
  const importAnnotationTable = (text) => replaceSets(parseAnnotationTable(text));
  const importAnnotationTableFile = async (event) => {
    const file = event?.target?.files?.[0];
    if (!file) return;
    importAnnotationTable(await file.text());
    event.target.value = '';
  };
  return {
    addAnnotationSet, renameAnnotationSet, duplicateAnnotationSet, removeAnnotationSet,
    addCoordinateAnnotation, addSelectedFeatures, removeAnnotation, setAnnotationTargetKind,
    importAnnotationTable, importAnnotationTableFile, replaceAnnotationSets: replaceSets
  };
};
