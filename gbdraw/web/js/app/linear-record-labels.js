// @ts-check
// The Record Labels of a Linear drawing: the definition line rows and their
// styles (size, weight, color), which fields Auto shows or hides and why, and
// the focus of a field's visibility control. The composition root passes the
// state, the mode, the session admission, the committed visibility reader, the
// settings pane's reveal, and the optional number input helpers as ports.
import {
  describeLinearLabelVisibility,
  requireLinearLabelVisibilityMode,
  resolveLinearLabelVisibility
} from '../services/linear-label-visibility.js';
import { linearRecordLayoutHasSharedRow } from '../services/linear-record-layout.js';
import { colorValueMode, toNativeColorInputValue } from '../utils/color-utils.js';
import { joinSentences, listText, t } from '../utils/i18n.js';
/** @import { DiagramMode, DrawingState, SessionBusy, WebState } from '../state.js' */

const { computed } = window.Vue;

export const definitionLineStyleRows = Object.freeze([
  { key: 'name', get label() { return t('Name / Species'); } },
  { key: 'subtitle', get label() { return t('Subtitle'); } },
  { key: 'replicon', get label() { return t('Replicon'); }, visibilityKey: 'linear_show_replicon', visibilityType: 'boolean' },
  { key: 'accession', get label() { return t('Accession'); }, visibilityKey: 'linear_accession_visibility', visibilityType: 'mode' },
  { key: 'length', get label() { return t('Length / Coordinates'); }, visibilityKey: 'linear_length_visibility', visibilityType: 'mode' }
]);
/** @param {DrawingState} drawing @param {string} kind */
const ensureDefinitionLineStyle = (drawing, kind) => {
  const key = String(kind || '');
  if (
    !drawing.adv.linear_definition_line_styles ||
    typeof drawing.adv.linear_definition_line_styles !== 'object' ||
    Array.isArray(drawing.adv.linear_definition_line_styles)
  ) {
    drawing.adv.linear_definition_line_styles = {};
  }
  const existing = drawing.adv.linear_definition_line_styles[key];
  if (!existing || typeof existing !== 'object' || Array.isArray(existing)) {
    drawing.adv.linear_definition_line_styles[key] = {
      font_size: null,
      font_weight: null,
      fill: null
    };
  }
  return drawing.adv.linear_definition_line_styles[key];
};

/**
 * @typedef {object} LinearRecordLabelsPorts
 * @property {WebState} state
 * @property {{ readonly value: DiagramMode }} mode The shown drawing's mode.
 * @property {() => SessionBusy | null} sessionOperationAvailability
 * @property {() => Record<string, boolean | null>} getCommittedLinearDefinitionVisibility
 *   The Linear Accession and Length visibility the shown Result was drawn with.
 * @property {(target: HTMLElement) => Promise<void>} revealSetting The settings pane shows and focuses a control.
 * @property {(value: unknown) => unknown} optionalNumberInputValue An optional number's input value.
 * @property {(target: Record<string, unknown> | null | undefined, key: string, value: unknown) => SessionBusy | undefined} setOptionalNumberInputValue
 *   An optional number's write from its input.
 */

/**
 * The Record Labels controls of the shown Linear drawing.
 * @param {LinearRecordLabelsPorts} ports
 */
export const createLinearRecordLabels = ({
  state, mode, sessionOperationAvailability, getCommittedLinearDefinitionVisibility, revealSetting,
  optionalNumberInputValue, setOptionalNumberInputValue
}) => {
  const linearLabelHasSharedRow = computed(() => {
    const drawing = state.activeDrawing();
    return linearRecordLayoutHasSharedRow(
      drawing.drawnLinearSeqs.value,
      drawing.linearRecordRows,
      { enabled: Boolean(drawing.linearRecordLayoutEnabled.value) }
    );
  });
  /** @param {unknown} mode */
  const linearLabelVisibilitySummary = (mode) => describeLinearLabelVisibility(mode, {
    hasSharedRow: linearLabelHasSharedRow.value
  });
  const linearLabelAutoFields = computed(() => {
    const drawing = state.activeDrawing();
    return definitionLineStyleRows.filter((row) => (
      row.visibilityType === 'mode'
      && requireLinearLabelVisibilityMode(drawing.adv[row.visibilityKey]) === 'auto'
    ));
  });
  const linearLabelAutoDisclosure = computed(() => {
    if (!linearLabelAutoFields.value.length) return '';
    const fields = listText(linearLabelAutoFields.value.map((/** @type {{ label: string }} */ row) => row.label));
    const shown = resolveLinearLabelVisibility('auto', {
      hasSharedRow: linearLabelHasSharedRow.value
    });
    // The shown Result already has Auto's outcome: say what it shows.
    const applied = getCommittedLinearDefinitionVisibility();
    const effect = linearLabelAutoFields.value
      .every((/** @type {{ visibilityKey: string }} */ row) => applied[row.visibilityKey] === shown)
      ? (shown ? t('Auto shows these fields throughout the diagram') : t('Auto hides these fields throughout the diagram'))
      : (shown
        ? t('Auto will show these fields throughout the diagram on the next successful Generate')
        : t('Auto will hide these fields throughout the diagram on the next successful Generate'));
    return shown
      ? t('{fields}: {effect} because no rendered row contains multiple records.', { fields, effect })
      : joinSentences(t('{fields}: {effect} because at least one rendered row contains multiple records.', { fields, effect }),
        t('Choose Show in Record Labels to keep a field visible.'));
  });
  /** @param {string} key */
  const focusLinearLabelVisibility = async (key) => {
    if (mode.value !== 'linear') return;
    const select = document.getElementById(`linear-label-visibility-${key}`);
    if (select) await revealSetting(select);
  };
  /** @param {string} kind */
  const getDefinitionLineStyleSize = (kind) => {
    const drawing = state.activeDrawing();
    return optionalNumberInputValue(
      ensureDefinitionLineStyle(drawing, kind).font_size
    );
  };
  /** @param {string} kind @param {unknown} value */
  const setDefinitionLineStyleSize = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    setOptionalNumberInputValue(ensureDefinitionLineStyle(drawing, kind), 'font_size', value);
  };
  /** @param {string} kind */
  const getDefinitionLineStyleWeight = (kind) => {
    const drawing = state.activeDrawing();
    return ensureDefinitionLineStyle(drawing, kind).font_weight ?? '';
  };
  /** @param {string} kind @param {unknown} value */
  const setDefinitionLineStyleWeight = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = String(value || '').trim().toLowerCase();
    ensureDefinitionLineStyle(drawing, kind).font_weight = normalized === 'bold' ? 'bold' : null;
  };
  /** @param {string} kind */
  const getDefinitionLineStyleFill = (kind) => {
    const drawing = state.activeDrawing();
    return ensureDefinitionLineStyle(drawing, kind).fill ?? '';
  };
  /** @param {string} kind @param {unknown} value */
  const setDefinitionLineStyleColor = (kind, value) => {
    const drawing = state.activeDrawing();
    const sessionBusy = sessionOperationAvailability();
    if (sessionBusy) return sessionBusy;
    const normalized = String(value || '').trim();
    ensureDefinitionLineStyle(drawing, kind).fill = normalized || null;
  };
  /** @param {string} kind */
  const getDefinitionLineStyleColorMode = (kind) => (
    colorValueMode(getDefinitionLineStyleFill(kind))
  );
  /** @param {string} kind */
  const getDefinitionLineStyleSwatchValue = (kind) => {
    return toNativeColorInputValue(getDefinitionLineStyleFill(kind));
  };
  /** @param {{ visibilityKey?: string, visibilityType?: string }} row */
  const isDefinitionLineStyleMuted = (row) => {
    const drawing = state.activeDrawing();
    const key = row?.visibilityKey;
    if (!key) return false;
    if (row.visibilityType === 'mode') {
      return !resolveLinearLabelVisibility(drawing.adv[key], {
        hasSharedRow: linearLabelHasSharedRow.value
      });
    }
    return drawing.adv[key] === false;
  };
  return {
    linearLabelHasSharedRow, linearLabelVisibilitySummary, linearLabelAutoFields, linearLabelAutoDisclosure,
    focusLinearLabelVisibility, getDefinitionLineStyleSize, setDefinitionLineStyleSize, getDefinitionLineStyleWeight,
    setDefinitionLineStyleWeight, getDefinitionLineStyleFill, setDefinitionLineStyleColor, getDefinitionLineStyleColorMode,
    getDefinitionLineStyleSwatchValue, isDefinitionLineStyleMuted
  };
};
