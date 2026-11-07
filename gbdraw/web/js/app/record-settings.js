// @ts-check
// The record list's records and the transitions it asks the composition root
// for: each record's settings owner (a Linear card, or a record of the
// Circular input) and the edits that name it, the record sources the list
// shows, and Delete settings, Remove File, and the read of the shown drawing's
// records after a Session load. The root creates the record selection owner
// (`createRecordSelection`) with these as its ports. The root passes the
// state, the shown-drawing readers, the Linear File cards, the record
// selector's and the record catalog's readers, the rule rows' take-back, the
// palette and rule projection, the Linear File and Circular source edits, the
// dialog and History pending readers, and the Circular source inspection as
// ports.
import { annotationRecordKeysOf } from './annotations/record-catalog.js';
import { ANNOTATION_RECORD_BINDING_KEY } from '../services/annotation-state.js';
import { removeAnnotations } from '../services/drawing/annotation-sets.js';
import { namedFeatureInstance } from '../services/feature-selector.js';
import { CIRCULAR_INPUT_FIELDS, circularInputSlots, inputSourceFiles, inputTypeSpec } from '../services/input-types.js';
import { groupLinearSourceRecords, planLinearSourceRemoval } from '../services/linear-sources.js';
import {
  circularRecordDrawKey,
  circularRecordRequestKeys,
  linearCardHasSettings,
  RECORD_SETTINGS_DEFAULTS,
  recordsWithEdits,
  removeRecordEdits
} from '../services/record-draw-selection.js';
import { buildDisambiguatedRecordEntries, resolveDisambiguatedRecordSelection } from '../services/record-options.js';
import { msg } from '../utils/i18n.js';
/** @import { DiagramMode, DrawingState, LinearSeq, SessionBusy, WebState } from '../state.js' */
/** @import { CircularRecordEntry } from '../services/drawing/inputs.js' */
/** @import { LinearSourceGroup } from '../services/linear-sources.js' */
/** @import { CircularRecordDiscovery, SequenceRecord } from '../services/record-discovery.js' */
/** @import { RecordDrawKey, RecordEditData, RecordSettingsOwner } from '../services/record-draw-selection.js' */
/** @import { SpecificColorRule } from '../services/specific-color-rules.js' */
/** @import { AnnotationCatalogSource } from './annotations/record-catalog.js' */
/** @import { ApplyLinearSeqMutation, AskLinearFileRemoval, RemoveLinearSource } from './linear-source-edits.js' */
/** @import { RecordSelectionMode, RecordSource } from './record-selection.js' */

/** @typedef {(drawingMode: DiagramMode) => DrawingState | null} ShownDrawingOfMode The shown drawing when it is of that mode. */
/**
 * What a Circular drawing's record discovery holds for its current input.
 * @typedef {{ hasInput: boolean, status: CircularRecordDiscovery['status'], records: CircularRecordEntry[] }} CircularDiscoveryState
 */
/** @typedef {(drawing: DrawingState) => CircularDiscoveryState} CircularDiscoveryStateOf */

/**
 * @typedef {object} RecordSettingsSourcePorts
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {{ readonly value: LinearSourceGroup[] }} linearSourceGroups The shown Linear drawing's Files.
 * @property {(source: LinearSourceGroup) => string} linearSourceFileName The name a File's card shows.
 * @property {(drawing: DrawingState, seq: LinearSeq) => SequenceRecord[]} linearRecordsFor
 *   The records a Linear card's selector holds (the Linear record selector's `recordsFor`).
 * @property {(drawing: DrawingState, seq: LinearSeq, drawn: boolean) => AnnotationCatalogSource} linearAnnotationSource
 *   A Linear card as the annotation record catalog reads it.
 * @property {CircularDiscoveryStateOf} getCircularRecordDiscoveryState
 * @property {(drawing: DrawingState) => string} circularAnnotationSourceKey The Circular input's annotation source key.
 * @property {(drawing: DrawingState, keep: (rules: SpecificColorRule[]) => SpecificColorRule[]) => unknown} takeBackRuleRows
 *   Takes back the drawing's color rule rows that `keep` leaves out, in the open History step.
 */

/**
 * Registers nothing. The root creates it with the transitions below, before
 * the record selection owner that takes `recordSourceOf`.
 * @param {RecordSettingsSourcePorts} ports
 */
export const createRecordSettingsSources = ({
  shownDrawingOfMode, linearSourceGroups, linearSourceFileName, linearRecordsFor, linearAnnotationSource,
  getCircularRecordDiscoveryState, circularAnnotationSourceKey, takeBackRuleRows
}) => {
  // Record selection (D-01..D-10): the record sources and edits its owner
  // reads, and the transitions it asks the root for.
  /** @param {DrawingState} drawing the drawing that holds the card @param {LinearSeq} seq The discovered record a Linear card draws. */
  const linearCardRecord = (drawing, seq) => {
    const selection = resolveDisambiguatedRecordSelection(linearRecordsFor(drawing, seq), seq.region_record_id);
    return selection.record || (selection.entries.length === 1 ? selection.entries[0] : null);
  };
  /** @param {DrawingState} drawing the drawing that holds the card @param {LinearSeq} seq @returns {RecordSettingsOwner} */
  const linearRecordSettingsOwner = (drawing, seq) => ({
    key: seq.uid, requestKeys: [seq.uid], ownsExpansions: true,
    bindingKeys: annotationRecordKeysOf(linearAnnotationSource(drawing, seq, false), seq.uid),
    displaySource: seq.uid, displaySelector: null
  });
  /** @param {DrawingState} drawing */
  const circularRecordEntries = (drawing) => buildDisambiguatedRecordEntries(
    getCircularRecordDiscoveryState(drawing).records.map((record) => ({
      ...record, recordId: record?.record_id ?? record?.recordId, recordLength: record?.record_length ?? record?.recordLength
    }))
  );
  /** @returns {RecordSettingsOwner[]} */
  /** @param {DrawingState} drawing */
  const circularRecordSettingsOwners = (drawing) => {
    const discovery = getCircularRecordDiscoveryState(drawing);
    const bindingKeys = annotationRecordKeysOf(
      { sourceKey: circularAnnotationSourceKey(drawing), hasInput: discovery.hasInput, status: discovery.status, records: discovery.records },
      'circular-source'
    );
    return circularRecordEntries(drawing).map((entry, index) => ({
      key: circularRecordDrawKey(entry),
      requestKeys: circularRecordRequestKeys(entry),
      bindingKeys: bindingKeys[index] ? [bindingKeys[index]] : [],
      displaySource: 'circular', displaySelector: circularRecordDrawKey(entry)
    }));
  };
  /** The record settings owners of a drawing's records. @param {DrawingState} drawing */
  const recordSettingsOwners = (drawing) => (drawing.mode === 'circular' ? circularRecordSettingsOwners(drawing) : drawing.linearSeqs.map((seq) => linearRecordSettingsOwner(drawing, seq)));
  /** @param {DrawingState} drawing @returns {RecordEditData} */
  const recordEditData = (drawing) => ({
    featureEdits: drawing.slices.featureEdits,
    paint: drawing.slices.paint,
    annotationSets: drawing.annotationSets,
    removeAnnotation: (/** @type {object} */ annotation) => removeAnnotations(drawing.slices.annotationSets, [annotation]),
    colorRuleEdits: drawing.manualSpecificRules.flatMap((row) => {
      const instance = namedFeatureInstance(row);
      return instance ? [{
        recordKey: instance.recordKey,
        remove: () => takeBackRuleRows(drawing, (rules) => rules.filter((rule) => rule !== row))
      }] : [];
    }),
    recordDisplayDrafts: drawing.recordDisplayDrafts,
    comparisonEdges: drawing.mode === 'linear' ? drawing.linearComparisonPlan.edges : [],
    annotationBindingField: ANNOTATION_RECORD_BINDING_KEY
  });
  /** @param {RecordSelectionMode} listMode @param {string} sourceKey @returns {RecordSource | null} */
  const recordSourceOf = (listMode, sourceKey) => {
    const drawing = shownDrawingOfMode(listMode);
    if (!drawing) return null;
    if (listMode === 'circular') {
      if (sourceKey !== 'circular' || getCircularRecordDiscoveryState(drawing).status !== 'ready') return null;
      const owners = circularRecordSettingsOwners(drawing);
      if (owners.length === 0) return null;
      const edited = recordsWithEdits(owners, recordEditData(drawing));
      const entries = circularRecordEntries(drawing);
      return {
        key: 'circular',
        name: (({ primaryFile, pairedFile }) => [primaryFile, pairedFile])(inputSourceFiles(drawing.inputType.value, circularInputSlots(drawing.files)))
          .filter(Boolean).map((file) => String(file?.name || 'Unnamed file')).join(' + '),
        records: owners.map((owner, index) => ({
          key: owner.key, recordId: entries[index].recordId, length: entries[index].recordLength ?? null,
          hasSettings: edited.has(owner.key)
        }))
      };
    }
    const group = linearSourceGroups.value.find((entry) => entry.uid === sourceKey);
    if (!group) return null;
    const cards = group.records.map(({ sequence }) => sequence);
    const edited = recordsWithEdits(cards.map((seq) => linearRecordSettingsOwner(drawing, seq)), recordEditData(drawing));
    return {
      key: group.uid,
      name: linearSourceFileName(group),
      records: cards.map((seq) => {
        const record = linearCardRecord(drawing, seq);
        return {
          key: seq.uid, recordId: String(record?.recordId || seq.region_record_id || ''),
          length: record?.recordLength ?? null, hasSettings: edited.has(seq.uid) || linearCardHasSettings(seq)
        };
      })
    };
  };
  return { circularRecordEntries, recordSettingsOwners, recordEditData, recordSourceOf };
};

/**
 * @typedef {object} RecordSettingsTransitionPorts
 * @property {WebState} state
 * @property {ShownDrawingOfMode} shownDrawingOfMode
 * @property {(drawing: DrawingState) => RecordSettingsOwner[]} recordSettingsOwners
 *   The record settings owners of a drawing's records (`createRecordSettingsSources`).
 * @property {(drawing: DrawingState) => RecordEditData} recordEditData The edits that name a drawing's records.
 * @property {<T>(drawing: DrawingState, change: () => T) => T} retireLegendStylesOfUnnamedCaptions
 *   Runs `change` and retires the Legend styles of the captions it leaves unnamed.
 * @property {() => unknown} projectPaletteAndRules Projects the palette and the color rules' matches again.
 * @property {ApplyLinearSeqMutation} applyLinearSeqMutation
 * @property {(field: 'c_gb' | 'c_gff' | 'c_fasta', value: File | null) => unknown} setCircularSourceFile
 *   Writes a Circular source slot of the shown Circular drawing.
 * @property {RemoveLinearSource} removeLinearSource
 * @property {AskLinearFileRemoval} askLinearFileRemoval
 * @property {{ readonly value: boolean }} dialogChoicePending Whether a dialog choice is open.
 * @property {() => boolean} mutationPending Whether a History step is open.
 * @property {(drawing: DrawingState) => SessionBusy | Promise<unknown>} refreshLinearRecords
 *   Reads the records of the drawing's Linear cards (the Linear record selector's refresh).
 * @property {CircularDiscoveryStateOf} getCircularRecordDiscoveryState
 * @property {() => unknown} inspectCircularSourceRecords Reads the records of the shown Circular drawing's input.
 */

/**
 * Registers nothing. The root creates it after the record sources, before
 * the record selection owner that takes these transitions.
 * @param {RecordSettingsTransitionPorts} ports
 */
export const createRecordSettingsTransitions = ({
  state, shownDrawingOfMode, recordSettingsOwners, recordEditData, retireLegendStylesOfUnnamedCaptions,
  projectPaletteAndRules, applyLinearSeqMutation, setCircularSourceFile, removeLinearSource, askLinearFileRemoval,
  dialogChoicePending, mutationPending, refreshLinearRecords, getCircularRecordDiscoveryState, inspectCircularSourceRecords
}) => {
  /**
   * Delete settings (D-08), inside the owner's History step: the records' card
   * settings and rotation, feature edits and placements, annotations, and
   * comparison pairs.
   * @param {RecordSelectionMode} listMode
   * @param {RecordDrawKey[]} keys
   */
  const deleteRecordSettings = (listMode, keys) => {
    const drawing = shownDrawingOfMode(listMode);
    if (!drawing) return false;
    const deleted = new Set(keys);
    const owners = recordSettingsOwners(drawing).filter((owner) => deleted.has(owner.key));
    const ruleCount = drawing.manualSpecificRules.length;
    retireLegendStylesOfUnnamedCaptions(drawing, () => removeRecordEdits(owners, recordEditData(drawing)));
    // The color rows it took back color nothing drawn; the rule matches follow the table.
    if (drawing.manualSpecificRules.length !== ruleCount) void projectPaletteAndRules();
    if (listMode === 'linear') {
      applyLinearSeqMutation(drawing, drawing.linearSeqs.map((seq) => (
        deleted.has(seq.uid) ? { ...seq, ...RECORD_SETTINGS_DEFAULTS } : seq
      )), { alignmentMutation: msg('record settings deleted.') });
    }
    return true;
  };
  /** @param {DrawingState} drawing @param {string} sourceKey The Linear File's key in the record list. */
  const recordFileRemovalPlan = (drawing, sourceKey) => {
    const intent = /** @type {'clear' | 'delete'} */ (groupLinearSourceRecords(drawing.linearSeqs).length > 1 ? 'delete' : 'clear');
    return { intent, plan: planLinearSourceRemoval({ sequences: drawing.linearSeqs, sourceUid: sourceKey, intent }) };
  };
  /**
   * D-06 Remove File, inside the dialog's History step: a Linear File is
   * deleted, or cleared when it is the only File; the Circular input loses its
   * GenBank file or its GFF3 and FASTA pair.
   * @param {RecordSelectionMode} listMode
   * @param {string} sourceKey
   */
  const removeRecordSource = (listMode, sourceKey) => {
    if (listMode === 'circular') {
      // The shown Circular drawing's input, which its File card removes from.
      const drawing = shownDrawingOfMode('circular');
      if (drawing === null) return undefined;
      const { primary, paired } = inputTypeSpec(drawing.inputType.value);
      if (paired) setCircularSourceFile(CIRCULAR_INPUT_FIELDS[primary], null);
      return setCircularSourceFile(CIRCULAR_INPUT_FIELDS[paired ?? primary], null);
    }
    const drawing = shownDrawingOfMode('linear');
    if (!drawing) return false;
    const { intent, plan } = recordFileRemovalPlan(drawing, sourceKey);
    return plan.allowed ? removeLinearSource(drawing, plan, intent) : false;
  };
  /** @param {RecordSelectionMode} listMode @param {string} sourceKey @param {string} label */
  const askRecordFileRemoval = (listMode, sourceKey, label) => {
    const drawing = shownDrawingOfMode('linear');
    if (listMode !== 'linear' || !drawing || dialogChoicePending.value || mutationPending()) return null;
    const { intent, plan } = recordFileRemovalPlan(drawing, sourceKey);
    return plan.allowed
      ? askLinearFileRemoval(drawing, plan.retainedSequences, label, (turnOffReaching) => removeLinearSource(drawing, plan, intent, turnOffReaching))
      : null;
  };
  // After a Session load the records of the shown drawing's inputs are read on demand, as Generate reads them: the
  // record list and Select region's Add to annotation set ask for them.
  /** @param {RecordSelectionMode} listMode */
  const readShownRecords = async (listMode) => {
    if (listMode === 'linear') return refreshLinearRecords(state.activeDrawing());
    const drawing = shownDrawingOfMode('circular');
    if (!drawing) return null;
    return getCircularRecordDiscoveryState(drawing).status === 'ready' ? null : inspectCircularSourceRecords();
  };
  return { deleteRecordSettings, removeRecordSource, askRecordFileRemoval, readShownRecords };
};
