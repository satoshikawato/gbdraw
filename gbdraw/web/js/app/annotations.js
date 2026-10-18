// @ts-check
/** @import { AskChoice } from './app-messages.js' */
/**
 * @import { AnnotationRecordCatalog } from './annotations/record-catalog.js'
 * @import { AnnotationItem, AnnotationSet, AnnotationStyle, AnnotationTarget } from '../services/annotation-state.js'
 * @import { AnnotationSetsSlice, AnnotationSetsView } from '../services/drawing/annotation-sets.js'
 */
import { createAnnotationSet, createDefaultAnnotationStyle, normalizeAnnotationSets, uniqueAnnotationSetId } from '../services/annotation-state.js';
import { nextCopyId } from '../services/track-slot-validation.js';
import { coordinateTarget, featureTarget, featureTargetsFromSelection } from './annotations/target-actions.js';
import { encodeAnnotationTableWithNotice, parseAnnotationTableWithNotice } from './annotations/table-codec.js';
import {
  createAnnotationRecordSelector,
  reconcileAnnotationRecordBindings
} from './annotations/record-selector.js';
import { getFeatureCaption, getFeatureColorRuleHash } from '../services/feature-utils.js';
import { ANNOTATION_SETS_SLICE, addAnnotationSet as addSetToDrawing, removeAnnotationSet as removeSetFromDrawing } from '../services/drawing/annotation-sets.js';
import { readFileText } from '../services/file-content-cache.js';
import { downloadTextFile } from '../services/text-download.js';
import { toRaw } from '../utils/vue-raw.js';

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

/**
 * @param {readonly AnnotationItem[]} annotations
 * @param {string} prefix
 */
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
 * @typedef {AnnotationSetsView & {
 *   adv: { circular_track_slots?: any[], linear_track_slots?: any[] },
 *   readonly slices: { readonly annotationSets: AnnotationSetsSlice }
 * }} AnnotationEditorDrawing The members of a drawing (`DrawingState` of state.js) the editor reads, and the
 *   handle the set list's writers take.
 * @typedef {object} AnnotationEditorOptions
 * @property {AnnotationEditorState} state
 * @property {() => AnnotationRecordCatalog | null | undefined} getRecordCatalog
 * @property {(notice: string) => void} onImportNotice Shows or clears the panel notice: the outcome of
 *   the last import or id edit.
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
    const lists = [state.extractedFeatures?.value, drawnOnly ? null : state.biologicalFeatures?.value];
    for (const list of lists.map(toRaw)) {
      const feature = (Array.isArray(list) ? list : []).find((item) => item?.record_key === target?.recordKey && item?.biological_feature_id === target?.biologicalFeatureId);
      if (feature) return feature;
    }
    return null;
  };
  /** @param {AnnotationItem | null | undefined} item */
  const featureTargetCaption = (item) => {
    const feature = catalogFeature(item?.target);
    return feature ? getFeatureCaption(feature) : `${item?.target?.biologicalFeatureId} (not in the current diagram)`;
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
  const reconcileRecords = (sets = state.activeDrawing().annotationSets) => {
    reconcileAnnotationRecordBindings(sets, getRecordCatalog?.());
    return sets;
  };
  /** @param {Partial<AnnotationSet>[] | null | undefined} sets */
  const replaceSets = (sets) => {
    const drawing = state.activeDrawing();
    const candidate = reconcileRecords(normalizeAnnotationSets(sets));
    retireLegendStylesOfUnnamedCaptions(() => {
      ANNOTATION_SETS_SLICE.install(drawing.slices.annotationSets, { annotationSets: candidate });
    });
  };
  const addAnnotationSet = (base = 'annotations') => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const set = createAnnotationSet({ id: uniqueAnnotationSetId(drawing.annotationSets, base) });
    addSetToDrawing(drawing.slices.annotationSets, set);
    return set;
  };
  /**
   * @param {AnnotationSet} set
   * @param {string} id
   */
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
  /**
   * @param {AnnotationSet | null | undefined} set
   * @param {unknown} value
   */
  const setAnnotationSetLegendLabel = (set, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!set) return;
    retireLegendStylesOfUnnamedCaptions(() => { set.legendLabel = String(value ?? '').trim(); });
  };
  /** @param {AnnotationSet} set */
  const duplicateAnnotationSet = (set) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const copy = createAnnotationSet(JSON.parse(JSON.stringify(set)));
    copy.id = nextCopyId(set.id, drawing.annotationSets);
    addSetToDrawing(drawing.slices.annotationSets, copy);
    return copy;
  };
  /** @param {AnnotationSet} set */
  const removeAnnotationSet = (set) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (drawing.annotationSets.includes(set)) {
      retireLegendStylesOfUnnamedCaptions(() => removeSetFromDrawing(drawing.slices.annotationSets, set));
    }
  };
  /**
   * @param {AnnotationSet | null | undefined} set
   * @param {{ start?: unknown, end?: unknown, recordId?: string | null, recordIndex?: number | null, coordinateSpace?: string }} [options]
   */
  const addCoordinateAnnotation = (set, options = {}) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!set) return null;
    const id = nextAnnotationId(set.annotations, 'region');
    const item = { id, target: coordinateTarget({ start: 1, end: 1, ...options }), label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
    set.annotations.push(item);
    return item;
  };
  /** @param {AnnotationSet | null | undefined} set */
  const addSelectedFeatures = (set) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!set) return [];
    const targets = featureTargetsFromSelection(state.selectedFeatures?.value ?? state.selectedFeatures ?? []);
    if (!targets) {
      onImportNotice?.('Generate the diagram again to annotate the selected features.');
      return [];
    }
    return targets.map((target) => {
      const item = { id: nextAnnotationId(set.annotations, 'feature'), target, label: '', mark: 'highlight', lane: null, style: null, legendLabel: null, metadata: {} };
      set.annotations.push(item);
      return item;
    });
  };
  /**
   * @param {AnnotationSet} set
   * @param {AnnotationItem} item
   */
  const removeAnnotation = (set, item) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const index = set?.annotations?.indexOf(item) ?? -1;
    if (index >= 0) retireLegendStylesOfUnnamedCaptions(() => set.annotations.splice(index, 1));
  };
  /**
   * @param {AnnotationSet} set
   * @param {AnnotationItem} item
   * @param {unknown} value
   */
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
  /**
   * @template {keyof AnnotationStyle} K
   * @param {AnnotationSet} set
   * @param {AnnotationItem} item
   * @param {K} field
   * @param {AnnotationStyle[K]} value
   */
  const setAnnotationStyle = (set, item, field, value) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    item.style ??= createDefaultAnnotationStyle(set.defaultStyle);
    item.style[field] = value;
  };
  /**
   * @param {AnnotationItem | null | undefined} item
   * @param {string} kind
   */
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
  /** @param {number} setCount */
  const askReplaceAll = (setCount) => askChoice({
    title: 'Remove all annotation sets?',
    description: `The imported table has no annotation rows. Importing it removes all ${setCount} annotation ${setCount === 1 ? 'set' : 'sets'}.`,
    actions: [{ id: 'remove', label: 'Remove all sets', tone: 'danger' }]
  });
  /** @param {string} text */
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
    onImportNotice?.(notice);
  };
  /** @param {{ target?: HTMLInputElement | null } | null | undefined} event The change event of the file input. */
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
    canDownloadAnnotationTable, downloadAnnotationTable, featureTargetCaption,
    recordOptionsFor: recordSelector.optionsFor,
    recordValueFor: recordSelector.valueFor,
    /**
     * @param {AnnotationItem | null | undefined} annotation
     * @param {unknown} value
     */
    setRecordValue: (annotation, value) => state.sessionOperationAvailability?.()
      || recordSelector.setValue(annotation, value),
    recordIsRequired: recordSelector.isRequired,
    recordIsMissing: recordSelector.isMissing,
    recordMissingMessage: recordSelector.missingMessage,
    recordShowsSelect: recordSelector.showsSelect,
    recordCaptionFor: recordSelector.captionFor,
    reconcileRecordBindings: reconcileRecords
  };
};
