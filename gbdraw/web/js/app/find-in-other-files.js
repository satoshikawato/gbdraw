// @ts-check
// Find in other files (PD-OI-123): the dialog, its search and proposal, and Create. Python owns every rule (the
// reference genes, the search plan, the windows the search ticks, the kept windows, the new drawing); this owner holds
// the dialog's own values (the reference rows, the files to search, the search settings, the ticks and edits) and
// sends them through the Worker operation `findInOtherFiles`. Nothing it holds reaches the Session or the History; only
// Create adds a project step, through the drawing list.
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { errorText } from '../utils/error-text.js';
import { t, tm, tn } from '../utils/i18n.js';
import { reactive } from '../utils/reactivity.js';
import { resolveComparisonThresholds } from '../mode-profiles.js';
import { formatRecordLength } from '../services/record-options.js';
import {
  LOCATE_OPERATION, fileNamesOf, locateLosatJobs, locateOutcomes, readDrawingReply, readFilesReply, readPlanReply,
  readReferenceReply, readSearchPlanReply, referenceWindowOf, wholeNumber
} from '../services/find-in-other-files.js';
import { derivedDrawingNotes } from './region-drawing.js';
import { proposalCountText, proposalView, windowChoices, windowDrafts, windowRefusalText } from './find-in-other-files/proposal.js';

/** @import { DrawingId } from '../services/drawing-view.js' */
/** @import { ProjectedDraftRequest, RegionDerivation } from '../services/config.js' */
/** @import { LosatOverlayRun } from '../services/losat.js' */
/** @import { LocateHelperPayload, LocateStep } from '../services/find-in-other-files.js' */
/** @import { DrawingListOutcome } from './drawing-list.js' */
/** @import { LinearSeq } from '../services/drawing/inputs.js' */
/** @import { ProposalFileView, WindowDraft } from './find-in-other-files/proposal.js' */
/**
 * @import { LocateDrawingReady, LocateFileRecord, LocateOptionsInput, LocatePlanReply, LocateReferenceInput,
 *   LocateSearchInput } from '../contracts/locate.generated.js'
 */

/**
 * @typedef {object} FindInOtherFilesPorts
 * @property {() => DrawingId} shownDrawingId
 * @property {() => Promise<ProjectedDraftRequest>} projectShownRequest The shown drawing's draft request.
 * @property {(payload: LocateHelperPayload) => Promise<{ result?: unknown }>} runLocateStep Runs the Worker operation
 *   `findInOtherFiles`; the reply is read here.
 * @property {(task: (run: LosatOverlayRun) => Promise<LocatePlanReply>) => Promise<{ status: 'ok', value: LocatePlanReply }
 *   | { status: 'canceled' } | { status: 'busy', reason: string }>} runInOverlay The search under the full-window overlay
 *   from its first step to its last (Python reads the Files, LOSAT by File, Python plans), with its status and Cancel.
 * @property {() => ShownDrawingInputs} shownDrawing The shown drawing's Pairwise Comparisons thresholds as its draft
 *   holds them (`adv`) and its File cards.
 * @property {(sourceId: DrawingId, derivation: RegionDerivation, kind: 'search') => Promise<DrawingListOutcome>}
 *   newDrawingFromSearch The drawing list's new drawing, right after the source drawing's tab.
 */

/** @typedef {{ adv: Readonly<Record<string, unknown>>, linearSeqs: readonly Pick<LinearSeq, 'gb' | 'gff' | 'fasta'>[] }} ShownDrawingInputs */
/** @typedef {{ id: number, recordKey: string, start: string, end: string, acrossOrigin: boolean }} ReferenceRow */
/** @typedef {{ fileKey: string, name: string, records: LocateFileRecord[], size: string }} SearchFileView */
/** @typedef {{ program: LocateSearchInput['program'], bitscore: string, evalue: string, identity: string, alignmentLength: string }} SearchSettings */
/** @typedef {{ key: string, start: string, end: string, reverseComplement: boolean, message: string }} WindowEdit */

/**
 * The dialog's view: display values and the user's input only. `screen`: the reference and the search, or the
 * proposal. `status`: 'loading', 'ready', 'checking' (the genes of the windows are counted), 'searching', 'deriving',
 * 'proposal', 'creating' or 'error'.
 * @typedef {object} FindDialogView
 * @property {boolean} show
 * @property {'search' | 'proposal'} screen
 * @property {'loading' | 'ready' | 'checking' | 'searching' | 'deriving' | 'proposal' | 'creating' | 'error'} status
 * @property {string} message A failure, or why a tick or an edit was refused.
 * @property {SearchFileView[]} files Every File of the drawing, in its order.
 * @property {ReferenceRow[]} rows
 * @property {string} referenceReason Why the reference windows cannot be searched, or ''.
 * @property {number | null} genes The number of query genes, once the windows are read.
 * @property {Record<string, boolean>} targets Whether each File is searched, by File key.
 * @property {SearchSettings} settings
 * @property {string} reference The proposal's line about the reference.
 * @property {string} cutGenes The reference CDS a window edge cuts, or ''.
 * @property {ProposalFileView[]} proposal
 * @property {string} count "N windows in M files".
 * @property {string} createReason Why Create is disabled, or ''.
 * @property {{ setting: string, label: string, explicit: string, auto: string }[]} returnedToAuto
 * @property {string[]} notCarried
 * @property {WindowEdit | null} edit
 */

/** How long typing in Start or End pauses before the genes are counted. */
const TYPING_PAUSE_MS = 250;

/** @returns {FindDialogView} */
const closedDialog = () => ({
  show: false, screen: 'search', status: 'loading', message: '', files: [], rows: [], referenceReason: '', genes: null,
  targets: {}, settings: { program: 'tlosatn', bitscore: '', evalue: '', identity: '', alignmentLength: '' }, reference: '',
  cutGenes: '', proposal: [], count: '', createReason: '', returnedToAuto: [], notCarried: [], edit: null
});

/** The record of a key among the Files of the drawing, with its File. @param {FindDialogView} dialog @param {string} key */
export const recordOf = (dialog, key) => dialog.files.flatMap((file) => file.records.map((record) => ({ file, record })))
  .find(({ record }) => record.recordKey === key) ?? null;

/** The Files that can be searched: every File but the reference's. @param {FindDialogView} dialog */
export const searchFiles = (dialog) => {
  const reference = dialog.rows.length ? recordOf(dialog, dialog.rows[0].recordKey)?.file.fileKey : undefined;
  return dialog.files.filter((file) => file.fileKey !== reference);
};

/** @type {Readonly<Record<'EMPTY' | 'NOT_WHOLE' | 'ORIGIN_NEEDS_START_AFTER_END', () => string>>} */
const ROW_REASONS = Object.freeze({
  EMPTY: () => t('Enter the Start and End of each reference window.'),
  NOT_WHOLE: () => t('Enter Start and End as whole numbers.'),
  ORIGIN_NEEDS_START_AFTER_END: () => t('For a window across the origin, enter a Start greater than the End.')
});

/**
 * The reference the rows give, or why they do not give one yet.
 * @param {FindDialogView} dialog
 * @returns {{ reference: LocateReferenceInput } | { reason: string }}
 */
const referenceInput = (dialog) => {
  const keys = new Set(dialog.rows.map((row) => row.recordKey));
  if (!dialog.rows.length) return { reason: t('Add a reference window.') };
  if (keys.has('')) return { reason: t('Choose the record of each reference window.') };
  if (keys.size > 1) return { reason: t('Choose the reference windows in one record.') };
  const windows = [];
  for (const row of dialog.rows) {
    const read = referenceWindowOf(row);
    if ('reason' in read) return { reason: ROW_REASONS[read.reason]() };
    windows.push(read.window);
  }
  return { reference: { recordKey: dialog.rows[0].recordKey, windows } };
};

/**
 * Why Search is disabled now, or ''.
 * @param {FindDialogView} dialog
 */
export const searchReason = (dialog) => {
  if (dialog.status === 'loading') return t('Reading the files of the drawing…');
  if (dialog.files.length < 2) return t('Add another file to search in');
  if (dialog.referenceReason) return dialog.referenceReason;
  if (dialog.status === 'checking' || dialog.genes === null) return t('Counting the genes of the windows…');
  if (!searchFiles(dialog).some((file) => dialog.targets[file.fileKey])) return t('Tick at least one file to search in.');
  return '';
};

/**
 * The dialog's search settings, from the Pairwise Comparisons thresholds of the drawing, each field on its own (its
 * mode's default for a blank or an invalid value).
 * @param {Readonly<Record<string, unknown>>} adv
 * @returns {SearchSettings}
 */
const settingsFrom = (adv) => {
  /** @param {string} draftKey @param {string} field */
  const one = (draftKey, field) => {
    try {
      return String(resolveComparisonThresholds({ [draftKey]: adv[draftKey] }, 'linear')[field]);
    } catch {
      return String(resolveComparisonThresholds({}, 'linear')[field]);
    }
  };
  return {
    program: 'tlosatn', bitscore: one('min_bitscore', 'bitscore'), evalue: one('evalue', 'evalue'),
    identity: one('identity', 'identity'), alignmentLength: one('alignment_length', 'alignmentLength')
  };
};

/**
 * The thresholds of a search, checked on Python's domains as Generate checks them (a typed INPUT_INVALID).
 * @param {SearchSettings} settings
 * @returns {LocateOptionsInput}
 */
const searchOptions = (settings) => {
  const resolved = resolveComparisonThresholds({
    min_bitscore: settings.bitscore, evalue: settings.evalue, identity: settings.identity, alignment_length: settings.alignmentLength
  }, 'linear');
  return {
    evalue: Number(resolved.evalue), bitscore: Number(resolved.bitscore), identity: Number(resolved.identity),
    alignmentLength: Number(resolved.alignmentLength)
  };
};

/** @param {Error | Record<string, unknown> | string | null} error */
const failureText = (error) => {
  const normalized = normalizeUserFacingError(error);
  return normalized ? errorText(normalized.summaryParts) : String(error);
};

/**
 * The Files of the drawing as Search in lists them.
 * @param {readonly { fileKey: string, records: LocateFileRecord[] }[]} files @param {Readonly<Record<string, string>>} names
 * @returns {SearchFileView[]}
 */
const fileViews = (files, names) => files.map(({ fileKey, records }) => ({
  fileKey, name: names[fileKey] ?? fileKey, records,
  size: t('{records} · {length}', {
    records: tn(records.length, '{count} record', '{count} records'),
    length: formatRecordLength(records.reduce((sum, record) => sum + record.length, 0))
  })
}));

/**
 * The proposal's line about the reference: its file, record, windows and genes.
 * @param {FindDialogView} dialog @param {LocatePlanReply} reply
 */
const referenceLine = (dialog, reply) => {
  const { reference } = reply.plan;
  return t('Reference: {record} in {file}, {windows}, {genes}.', {
    file: dialog.files.find((file) => file.fileKey === reference.file_key)?.name ?? reference.file_key,
    record: reference.record_id,
    windows: tn(reference.windows.length, '{count} window', '{count} windows'),
    genes: tn(reference.genes.length, '{count} gene', '{count} genes')
  });
};

/**
 * One search and its proposal: the source drawing and its request, the plan of the search, the window drafts, and the
 * derived drawing of the kept windows.
 * @typedef {{ sourceId: DrawingId, sent: ProjectedDraftRequest, names: Record<string, string>, search: LocatePlanReply | null,
 *   drafts: WindowDraft[], derived: LocateDrawingReady | null }} FindSession
 */

/**
 * What the dialog's parts share: its view, its search, one step of the Worker operation, and the request counter by
 * which a reply that a newer action superseded is dropped.
 * @typedef {object} FindContext
 * @property {FindDialogView} dialog
 * @property {() => FindSession | null} session
 * @property {() => number} begin Starts a request; its number.
 * @property {(mine: number) => boolean} current Whether no request started after `mine`.
 * @property {(sent: ProjectedDraftRequest, locate: LocateStep) => Promise<unknown>} step
 * @property {(error: Error | Record<string, unknown> | string | null) => void} fail
 */

/**
 * The reference rows: their edits, and the count of the genes they hold (Python's).
 * @param {FindContext} context
 */
const referenceRows = ({ dialog, session, begin, current, step, fail }) => {
  let rowId = 0;
  /** Counts the genes of the reference windows. */
  const checkReference = async () => {
    const mine = begin();
    const input = referenceInput(dialog);
    const found = session();
    if (!('reference' in input) || !found) {
      Object.assign(dialog, { status: 'ready', genes: null, referenceReason: 'reason' in input ? input.reason : '' });
      return;
    }
    Object.assign(dialog, { status: 'checking', genes: null, referenceReason: '', message: '' });
    try {
      const reply = readReferenceReply(await step(found.sent, { step: 'reference', reference: input.reference }));
      if (!current(mine)) return;
      Object.assign(dialog, reply.status === 'refused'
        ? { status: 'ready', referenceReason: windowRefusalText(reply.reason, reply.start, reply.end) }
        : { status: 'ready', genes: reply.genes });
    } catch (error) {
      if (current(mine)) fail(error);
    }
  };
  /**
   * The row of the picked span: the Region popup names its record by key and gives the span in source coordinates on
   * the source strand, so on a circular record a Start after the End is a span across the origin.
   * @param {{ recordKey: string, start: number, end: number } | null} region
   */
  const firstRow = (region) => {
    const found = region === null ? null : recordOf(dialog, region.recordKey);
    const picked = found ? region : null;
    dialog.rows = [{
      id: ++rowId, recordKey: picked ? picked.recordKey : '', start: picked ? String(picked.start) : '',
      end: picked ? String(picked.end) : '', acrossOrigin: Boolean(picked && found?.record.circular && picked.start > picked.end)
    }];
    return checkReference();
  };
  /** @type {ReturnType<typeof setTimeout> | undefined} */
  let typing;
  /** @type {((checked: Promise<void>) => void)[]} */
  let waiting = [];
  /** Drops a count that waits for a pause in typing (the dialog closes). */
  const stopWaiting = () => {
    clearTimeout(typing);
    for (const resolve of waiting) resolve(Promise.resolve());
    waiting = [];
  };
  /** Counts the genes once typing pauses: each keystroke in Start or End would read the reference file again. */
  const checkAfterTyping = () => {
    clearTimeout(typing);
    begin(); // a count still running is for the text before this keystroke
    Object.assign(dialog, { status: 'checking', genes: null, referenceReason: '' });
    return new Promise((resolve) => {
      waiting.push(resolve);
      typing = setTimeout(() => {
        const checked = checkReference();
        for (const done of waiting) done(checked);
        waiting = [];
      }, TYPING_PAUSE_MS);
    });
  };
  /**
   * A row edit, found by the row's id; a record that is not circular has no window across the origin. Start and End
   * are counted once typing pauses.
   * @param {number} id @param {Partial<Omit<ReferenceRow, 'id'>>} patch
   */
  const editRow = (id, patch) => {
    const row = dialog.rows.find((candidate) => candidate.id === id);
    if (!row) return Promise.resolve();
    Object.assign(row, patch);
    if (!recordOf(dialog, row.recordKey)?.record.circular) row.acrossOrigin = false;
    return 'start' in patch || 'end' in patch ? checkAfterTyping() : checkReference();
  };
  /** "+ Add window": a row of the last row's record. */
  const addRow = () => {
    dialog.rows.push({ id: ++rowId, recordKey: dialog.rows.at(-1)?.recordKey ?? '', start: '', end: '', acrossOrigin: false });
    return checkReference();
  };
  /** @param {number} id */
  const removeRow = (id) => {
    dialog.rows = dialog.rows.filter((row) => row.id !== id);
    return checkReference();
  };
  return { firstRow, editRow, addRow, removeRow, stopWaiting };
};

/**
 * The proposal: the ticks and edits of the windows, Back, and Create.
 * @param {FindContext} context
 * @param {FindInOtherFilesPorts['newDrawingFromSearch']} newDrawingFromSearch
 * @param {() => void} close
 */
const proposalActions = ({ dialog, session, begin, current, step, fail }, newDrawingFromSearch, close) => {
  /**
   * The kept windows' plan and drawing, from Python. A refused edit or tick returns its reason and keeps the proposal
   * shown before it.
   * @returns {Promise<string>}
   */
  const derive = async () => {
    const found = session();
    if (!found?.search) return '';
    const mine = begin();
    const { search, drafts, names } = found;
    dialog.status = 'deriving';
    try {
      const reply = readDrawingReply(await step(found.sent, {
        step: 'drawing', drawing: { plan: search.plan, choices: windowChoices(drafts), labels: names, adaptSizes: true }
      }));
      if (!current(mine)) return '';
      if (reply.status === 'refused') {
        dialog.status = 'proposal';
        return windowRefusalText(reply.reason, reply.start, reply.end);
      }
      found.derived = reply;
      Object.assign(dialog, {
        status: 'proposal', message: '', proposal: proposalView(search, reply, drafts, names), count: proposalCountText(drafts),
        createReason: drafts.some((draft) => draft.ticked) ? '' : t('Tick at least one window to create the drawing.'),
        ...derivedDrawingNotes(reply.drawing)
      });
    } catch (error) {
      if (current(mine)) fail(error);
    }
    return '';
  };
  /** @param {string} key */
  const draftOf = (key) => session()?.drafts.find((draft) => draft.key === key) ?? null;
  /**
   * A tick, or an edit of Start, End and Reverse complement: changes that window only. A change Python refuses is
   * undone, with the reason shown. While Python answers the last change the proposal takes none (the dialog is busy).
   * @param {string} key @param {Partial<Pick<WindowDraft, 'ticked' | 'start' | 'end' | 'reverseComplement'>>} change
   * @returns {Promise<string>} The refusal, or ''.
   */
  const changeWindow = async (key, change) => {
    const draft = draftOf(key);
    if (!draft || dialog.status !== 'proposal') return '';
    const before = { ...draft };
    Object.assign(draft, change, { edited: draft.edited || 'start' in change });
    const refused = await derive();
    if (refused) Object.assign(draft, before);
    dialog.message = refused;
    return refused;
  };
  /** @param {string} key */
  const startEdit = (key) => {
    const draft = draftOf(key);
    if (draft) dialog.edit = { key, start: String(draft.start), end: String(draft.end), reverseComplement: draft.reverseComplement, message: '' };
  };
  const applyEdit = async () => {
    const { edit } = dialog;
    if (!edit) return;
    const start = wholeNumber(edit.start);
    const end = wholeNumber(edit.end);
    if (start === null || end === null) {
      edit.message = t('Enter Start and End as whole numbers.');
      return;
    }
    const refused = await changeWindow(edit.key, { start, end, reverseComplement: edit.reverseComplement });
    if (refused) edit.message = refused;
    else dialog.edit = null;
  };
  /** Back: the search settings again; the proposal is dropped and nothing is added. */
  const back = () => {
    begin();
    const found = session();
    if (found) Object.assign(found, { search: null, drafts: [], derived: null });
    Object.assign(dialog, { screen: 'search', status: 'ready', message: '', proposal: [], edit: null, returnedToAuto: [], notCarried: [] });
  };
  /** Create: the new drawing as one project step, right after the source drawing; the dialog closes. */
  const confirm = async () => {
    const found = session();
    if (dialog.status !== 'proposal' || dialog.createReason || !found?.derived) return null;
    const { sourceId, sent, derived } = found;
    dialog.status = 'creating';
    try {
      const outcome = await newDrawingFromSearch(sourceId, { derived: derived.drawing, sent }, 'search');
      close();
      return outcome;
    } catch (error) {
      fail(error);
      return null;
    }
  };
  return {
    derive, startEdit, applyEdit, back, confirm,
    cancelEdit: () => { dialog.edit = null; },
    toggleWindow: (/** @type {string} */ key, /** @type {boolean} */ ticked) => changeWindow(key, { ticked })
  };
};

/**
 * Find in other files: the dialog's one owner (open, close, every input edit, Search, the proposal, Create); the
 * template assigns none of its state.
 * @param {FindInOtherFilesPorts} ports
 */
export const createFindInOtherFiles = ({
  shownDrawingId, projectShownRequest, runLocateStep, runInOverlay, shownDrawing, newDrawingFromSearch
}) => {
  /** @type {FindDialogView} */
  const dialog = reactive(closedDialog());
  /** @type {FindSession | null} */
  let found = null;
  let request = 0;
  /** @type {FindContext} */
  const context = {
    dialog,
    session: () => found,
    begin: () => ++request,
    current: (mine) => mine === request,
    step: async (sent, locate) => (await runLocateStep({
      canonicalRequest: sent.envelope.renderRequest, resources: sent.envelope.resources, locate
    })).result,
    fail: (error) => { Object.assign(dialog, { status: 'error', message: failureText(error) }); }
  };
  const { step, fail } = context;
  const rows = referenceRows(context);
  const close = () => {
    request += 1;
    found = null;
    rows.stopWaiting();
    Object.assign(dialog, closedDialog());
  };
  const proposal = proposalActions(context, newDrawingFromSearch, close);

  /**
   * Opens on the Files of the shown drawing, with one reference row: the span picked on the figure (Region popup),
   * or an empty row (Drawing menu). A close before the Files are read keeps it closed.
   * @param {{ recordKey: string, start: number, end: number } | null} [region]
   */
  const open = async (region = null) => {
    close();
    const mine = request;
    const sourceId = shownDrawingId();
    dialog.show = true;
    try {
      const sent = await projectShownRequest();
      const files = readFilesReply(await step(sent, { step: 'files' })).files;
      if (mine !== request) return;
      const names = fileNamesOf(sent.envelope.renderRequest.records, sent.filesByResource);
      found = { sourceId, sent, names, search: null, drafts: [], derived: null };
      Object.assign(dialog, {
        files: fileViews(files, names), targets: Object.fromEntries(files.map((file) => [file.fileKey, true])),
        settings: settingsFrom(shownDrawing().adv)
      });
      await rows.firstRow(region);
    } catch (error) {
      if (mine === request) fail(error);
    }
  };

  /**
   * The search's three steps under the overlay: Python reads the Files and gives one LOSAT job per File, LOSAT runs
   * them by File, and Python plans the windows. The dialog steps aside while the overlay is up ('searching').
   * @param {FindSession} session @param {LocateSearchInput} searchInput @param {LocateOptionsInput} options
   * @param {LosatOverlayRun} overlay
   */
  const searchSteps = async (session, searchInput, options, overlay) => {
    dialog.status = 'searching';
    overlay.status(t('Reading the files of the drawing…'));
    const { jobs, sequences } = locateLosatJobs(readSearchPlanReply(
      await overlay.until(step(session.sent, { step: 'search', search: searchInput }))));
    const results = await overlay.losat({ jobs, sequences, byFile: { operation: LOCATE_OPERATION } });
    overlay.status(t('Finding the windows in each file…'));
    return readPlanReply(await overlay.until(step(session.sent, {
      step: 'plan', plan: { search: searchInput, outcomes: locateOutcomes(jobs, results), options }
    })));
  };

  /** Search, under the overlay from its first step to its last; Cancel adds nothing and keeps the dialog on the search. */
  const search = async () => {
    const input = referenceInput(dialog);
    if (searchReason(dialog) || !found || !('reference' in input)) return;
    const mine = ++request;
    const session = found;
    dialog.message = '';
    try {
      const options = searchOptions(dialog.settings);
      /** @type {LocateSearchInput} */
      const searchInput = {
        reference: input.reference, program: dialog.settings.program,
        targets: searchFiles(dialog).filter((file) => dialog.targets[file.fileKey]).map((file) => file.fileKey)
      };
      const run = await runInOverlay((overlay) => searchSteps(session, searchInput, options, overlay));
      if (mine !== request) return;
      if (run.status !== 'ok') {
        Object.assign(dialog, { status: 'ready', message: run.status === 'busy' ? tm(run.reason) : '' });
        return;
      }
      const reply = run.value;
      Object.assign(session, { search: reply, drafts: windowDrafts(reply.plan) });
      const cut = reply.plan.reference.cut_genes.map((cds) => reply.geneNames[cds.protein_id] ?? cds.protein_id);
      Object.assign(dialog, {
        screen: 'proposal', edit: null, reference: referenceLine(dialog, reply),
        cutGenes: cut.length ? t('Cut by a window edge and not searched: {genes}.', { genes: cut.join(', ') }) : ''
      });
      await proposal.derive();
    } catch (error) {
      if (mine === request) fail(error);
    }
  };

  return {
    dialog, open, close, search, editRow: rows.editRow, addRow: rows.addRow, removeRow: rows.removeRow,
    toggleWindow: proposal.toggleWindow, startEdit: proposal.startEdit, applyEdit: proposal.applyEdit,
    cancelEdit: proposal.cancelEdit, back: proposal.back, confirm: proposal.confirm,
    setTarget: (/** @type {string} */ fileKey, /** @type {boolean} */ ticked) => { dialog.targets[fileKey] = ticked; },
    setSetting: (/** @type {Partial<SearchSettings>} */ patch) => { Object.assign(dialog.settings, patch); },
    searchReason: () => searchReason(dialog),
    searchFiles: () => searchFiles(dialog),
    editWindowField: (/** @type {Partial<Omit<WindowEdit, 'key' | 'message'>>} */ patch) => { if (dialog.edit) Object.assign(dialog.edit, patch); },
    /** Why the Region popup's Find in other files is disabled: a drawing of fewer than two files. */
    entryReason: () => (shownDrawing().linearSeqs.filter((seq) => seq.gb || seq.gff || seq.fasta).length < 2
      ? t('Add another file to search in') : '')
  };
};
