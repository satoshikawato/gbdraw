// @ts-check
/** @import { DrawingState, LegendEntry } from '../../state.js' */
/** @import { PythonLegendKey } from '../../services/legend-svg.js' */
/** @import { SpecificColorRule, SpecificRuleBuild } from '../../services/specific-color-rules.js' */
/** @import { AskChoice, ShowNotice } from '../app-messages.js' */
import { createSpecificRulePatternDrafts } from './pattern-drafts.js';
import { normalizeUserFacingError } from '../../utils/error-normalization.js';
import { runWhenPrepared } from '../rule-matching.js';
import { ruleMatcher, ruleMatchesFeature } from '../../services/rule-matchers.js';
import { appliedFeatureColorOf, resolveColorToHex } from '../../utils/color-utils.js';
import { parseSpecificRules, serializeSpecificRules } from '../../services/file-imports.js';
import { formatFeatureRange, getFeatureColorRuleHash } from '../../services/feature-utils.js';
import {
  buildLegendIntents, createRuleLegendCaptions, legendRowRules, prepareSpecificColorImport, rekeyedRuleRows,
  ruleCommitLegendRows, ruleLegendCaptions
} from '../../services/specific-color-rules.js';
import { resolveFeatureLabelSelector } from '../../services/feature-selector.js';
import { downloadTextFile } from '../../services/text-download.js';
import {
  featureOverrideKey,
  getFeatureOverride,
  migrateLegacyFeatureOverrides
} from '../../services/feature-override-identity.js';
import {
  defaultFeatureRendering,
  normalizeFeatureRendering
} from '../../utils/feature-rendering.js';
import { featureOverrideValue } from '../../services/feature-placement.js';
import { hasOwnEntry, ownEntry } from '../../utils/own-entry.js';
import {
  displayedLegendRowContext, featureDrawnContext, legendRowsShowable, legendRuleOrder, resultLegendRowKeys
} from '../../services/feature-visibility.js';
import { legendEntryKey, legendRowTakesPlace, shownLegendKeys, shownPythonLegendRow } from '../../services/legend-svg.js';
/** @import { PatternDraft } from './pattern-drafts.js' */
/** @import { ColorFeature } from './color-actions.js' */

// R13: `projectPaletteAndRules` is the composition root's projection of the
// palette and the specific-color rules (R3); this owner calls it after a rule
// commit, whose candidate rules it has prepared. The composition root
// registers the label owner's `requestAutomaticRerender` in `ports` once that
// owner exists; this owner only calls it.
/**
 * The label owner's reaction the composition root registers once that owner exists.
 * @typedef {object} RuleActionsPorts
 * @property {() => boolean} requestAutomaticRerender
 *   Asks for the automatic rerender of an edit Python must draw again (OV-43).
 * @property {(drawing: DrawingState) => boolean} isRedrawOwed
 *   Whether a rerender asked for in `drawing` is not yet shown drawn from its current paint (D-36).
 */

/**
 * The root's way to bring a setting of the settings pane into view and focus
 * it: the pattern field a failed rule write offers for editing.
 * @typedef {(setting: HTMLElement) => unknown} RevealSetting
 */

/**
 * A Legend row a rule commit draws, named by its name and color (the commit
 * allocates its key), and where it goes: with `replaces`, in the place of the
 * row of key `at`, which it replaces with every record stored under that key
 * (a whole row renamed, OV-158); else right after the row of key `at` (a row
 * split off, D-29).
 * @typedef {{ name: string, color: string, at: PythonLegendKey, replaces?: boolean }} LegendPlacement
 */

/**
 * The Legend rows a rule commit draws, prepared and not yet applied.
 * @typedef {object} PreparedFileLegend
 * @property {{ add: any[], update: any[], remove: any[], unchanged: any[] }} diff
 * @property {() => boolean} isCurrent
 * @property {() => import('../candidate-render.js').RuleLegendRows | null} apply Writes the rows into the
 *   intent; returns the rows the displayed Result shows at once.
 */

/**
 * @typedef {object} FeatureRuleActionsOptions
 * @property {Record<string, any>} state App state (state.js; not yet typed).
 * @property {(intents: Record<string, any>[], options?: { drawing?: DrawingState, previousFileIntents?: Record<string, any>[], isCurrent?: () => boolean, placement?: { caption: PythonLegendKey, after: PythonLegendKey } | null, rekeyed?: ReadonlyMap<PythonLegendKey, PythonLegendKey> }) => Promise<PreparedFileLegend | false>} prepareFileLegendEntries
 *   The Legend owner's preparation of the rows the rules draw.
 * @property {import('../rule-matching.js').RulePreparation} rulePreparation
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoable History's undoable step.
 * @property {(label: string, fn: () => any, options?: Record<string, any>) => any} runUndoableCheckpoint
 *   History's undoable step that stores a checkpoint of the Result.
 * @property {(options?: { recolor?: Record<string, any>, prepareRules?: boolean, legendRows?: import('../candidate-render.js').RuleLegendRows | null }) => boolean | Promise<boolean>} projectPaletteAndRules
 *   The root's projection of the palette and the rules (R3).
 * @property {RuleActionsPorts} ports
 * @property {() => ({ diagramOptions?: Record<string, any> } | null)} [getCommittedRequest]
 *   The committed canonical request (Python owns the option fields, R7).
 * @property {(value?: any) => { value: any }} ref Vue `ref`
 * @property {<T>(getter: () => T) => { value: T }} computed Vue `computed`
 * @property {(source: { value: any }, callback: (value: any) => void) => () => void} watch
 *   Vue `watch`; a rule commit waits on it for an automatic rerender in flight.
 * @property {() => boolean} [isPatternEditAvailable] False while a Session import is pending.
 * @property {AskChoice} [askChoice] The root's choice dialog: a pattern that may match slowly asks first.
 * @property {ShowNotice} [showNotice] The root's notice banner.
 * @property {RevealSetting} [revealSetting] Without a settings pane, the field takes the focus where it is.
 */

/** @param {FeatureRuleActionsOptions} options */
export const createFeatureRuleActions = ({ state, prepareFileLegendEntries, rulePreparation, runUndoable, runUndoableCheckpoint, projectPaletteAndRules, ports, getCommittedRequest = () => null, ref, computed, watch, isPatternEditAvailable = () => true, askChoice = async () => 'cancel', showNotice = () => {}, revealSetting = (setting) => setting.focus() }) => {
  const {
    newColorFeat,
    newColorVal,
    newSpecRule,
    specificRulePresets,
    selectedSpecificPreset,
    specificRulePresetLoading,
    newPriorityRule,
    newFeatureToAdd,
    extractedFeatures,
    editableLabels
  } = state;

  /** @param {unknown} value */
  const normalizeCaption = (value) => String(value || '').trim();
  /** @param {unknown} value */
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  /** @param {unknown} value */
  const normalizeFeatureIdKey = (value) => String(value || '').trim().toLowerCase();
  /** @param {unknown} value @param {unknown} target */
  const captionMatches = (value, target) => normalizeCaptionKey(value) === normalizeCaptionKey(target);
  const specificRuleFields = new Set(['feat', 'qual', 'val', 'color', 'cap']);

  const patternDrafts = createSpecificRulePatternDrafts({
    rules: () => state.activeDrawing().manualSpecificRules, ref, changed: () => releaseCancelledWrites()
  });
  const ruleFailure = ref(null);
  const canRetrySpecificRuleFailure = computed(() => Boolean(ruleFailure.value
    && state.errorLog?.value === ruleFailure.value.error
    && rulePreparation.isCurrent(ruleFailure.value.snapshot)
    && (!ruleFailure.value.presetId || selectedSpecificPreset.value === ruleFailure.value.presetId)));
  const canEditSpecificRuleFailure = computed(() => canRetrySpecificRuleFailure.value
    && Boolean(ruleFailure.value.input?.isConnected));
  const retrySpecificRuleFailure = () => {
    const busy = state.sessionOperationAvailability?.();
    return busy || (canRetrySpecificRuleFailure.value ? ruleFailure.value.retry() : false);
  };
  const editSpecificRuleFailure = () => {
    if (canEditSpecificRuleFailure.value) revealSetting(ruleFailure.value.input);
  };
  // The legend rows the candidate rules draw on the features the displayed
  // Result draws (OV-344: not on another batch Result's features or a hidden
  // feature, where Generate draws no row), and the rows
  // the commit retires: those the current rules draw (so it retires the row
  // Generate drew), the listed rows of the Python rows the rules take (by key,
  // `takenKeys`, with the caption and the swatch they show, OV-310), and
  // `retiredLegendIntents`, rows this commit replaces, which are no
  // renderer rows for the N-06 caption allocation. `removedRuleRows` are
  // the rows of current rules that no candidate rule names, drawn or not;
  // `rekeyed` are the rows both rule tables hold under another key (OV-645).
  /**
   * @param {DrawingState} drawing
   * @param {Partial<SpecificColorRule>[]} candidateRules
   * @param {Array<{ caption?: string, color?: string }>} retiredLegendIntents
   */
  const candidateLegendIntents = (drawing, candidateRules, retiredLegendIntents) => {
    const displayed = displayedLegendRowContext(state, drawing);
    const candidateMatches = ruleMatcher(candidateRules);
    const used = new Set([...displayed.features].map(feature => candidateMatches.first(feature)).filter(Boolean));
    const { rendererRows, takenKeys } = ruleCommitLegendRows({
      ...displayed,
      rules: [...drawing.manualSpecificRules, ...candidateRules,
        ...retiredLegendIntents.map(intent => ({ cap: intent?.caption, color: intent?.color }))]
    });
    const currentCaption = createRuleLegendCaptions(drawing.manualSpecificRules, rendererRows);
    /** @type {Array<{ caption: PythonLegendKey, color: string }>} */
    const ruleRows = drawing.manualSpecificRules.filter((/** @type {SpecificColorRule} */ rule) => rule.cap)
      .map((/** @type {SpecificColorRule} */ rule) => ({ caption: currentCaption(rule), color: rule.color }));
    /** @type {LegendEntry[]} */
    const listed = drawing.legendEntries.value || [];
    const takenRows = listed.filter(entry => takenKeys.has(legendEntryKey(entry)))
      .map(entry => ({ caption: legendEntryKey(entry), color: entry.color }));
    const candidateCaptions = new Set(buildLegendIntents(candidateRules, rendererRows).intents.map(intent => intent.caption));
    return {
      intents: buildLegendIntents(candidateRules, rendererRows, rule => used.has(rule)).intents,
      previousIntents: [...ruleRows, ...takenRows, ...retiredLegendIntents],
      removedRuleRows: ruleRows.filter(row => !candidateCaptions.has(row.caption)),
      rekeyed: rekeyedRuleRows(drawing.manualSpecificRules, candidateRules, rendererRows)
    };
  };
  // OV-43 (Owner decision 2026-10-06, option A), narrowed by R14-8 (U3a 1d):
  // a rule table change asks for the automatic rerender when the displayed
  // Result cannot show the Legend rows it regroups (`legendRowsShowable`): a
  // row that loses part of its features, a row built from parts of others, a
  // merge into a new row, an unknown match, or a row of another batch Result.
  // Python then redraws the rows, their order, and the "other <type>s" rows.
  // The rows are read with the current Feature visibility before and after.
  /** @param {DrawingState} drawing @param {Partial<SpecificColorRule>[]} rules */
  const legendSourceContext = (drawing, rules) => ({
    ...featureDrawnContext(drawing, { diagramOptions: getCommittedRequest()?.diagramOptions }),
    colorRules: rules
  });
  // OV-152: rules that recolor a whole type keep its caption, so their Legend
  // row is that of the type's default row and the live Legend follows them.
  // A change that removes the last such rule of a type makes Python draw the
  // default row again, which the live removal of the rule row does not; the
  // automatic rerender draws it.
  /** @param {Record<string, any>[]} rules */
  const typeRowCaptions = (rules) => new Set(rules
    .filter(rule => String(rule?.cap ?? '').trim() && rule.cap === rule.feat)
    .map(rule => rule.cap));
  /**
   * @param {DrawingState} drawing
   * @param {Record<string, any>[]} before
   * @param {Record<string, any>[]} after
   * @param {{ rows?: import('../candidate-render.js').RuleLegendRows | null, replaced?: ReadonlyMap<string, string>, holdsRow?: (caption: string) => boolean }} [options]
   *   `rows`: the rows the commit's own show adds and retires on the displayed
   *   Result (else the Result shows what History restored); `replaced`: the
   *   rows the commit's show draws in the place of the rows they replace (old
   *   key to new key); `holdsRow`: whether the Legend already holds the default
   *   row of a type.
   */
  const redrawsLegendFor = (drawing, before, after, { rows = null, replaced = new Map(), holdsRow = () => false } = {}) => {
    const kept = typeRowCaptions(after);
    const svg = state.svgContainer?.value?.querySelector?.('svg');
    const drawnKeys = shownLegendKeys(svg);
    const added = new Set((rows?.add || []).map(({ caption }) => caption));
    const retired = new Set(rows?.retire || []);
    const shows = (/** @type {string} */ key) => added.has(key) || (!retired.has(key) && drawnKeys.has(key));
    // The commit's show appends a new row unless it takes the place of the row
    // it relabels (OV-158), where Python keeps it. A restore shows the bytes
    // History kept, which may be a commit's before its rerender (M1): its new
    // row takes the old row's place only as the commit's placement left it.
    // The rows a restore splits are exact when they are Python's (L1).
    const placed = rows
      ? (/** @type {string} */ key, /** @type {string} */ nextKey) => replaced.get(key) === nextKey
      : (/** @type {string} */ key, /** @type {string} */ nextKey) => legendRowTakesPlace(svg, key, nextKey);
    // A commit keys rule rows as Python does (OV-1235); a restore by caption (D-15-38).
    const byColor = Boolean(rows);
    return !legendRowsShowable(
      resultLegendRowKeys(state, legendSourceContext(drawing, before), { byColor }),
      resultLegendRowKeys(state, legendSourceContext(drawing, after), { byColor }),
      {
        displayed: Number(state.selectedResultIndex?.value) || 0,
        shows,
        placed,
        firstRule: legendRuleOrder(before, { byColor }),
        restored: rows ? null : (/** @type {string} */ key) => shownPythonLegendRow(svg, key)
      }
    ) || [...typeRowCaptions(before)].some(caption => !kept.has(caption) && !holdsRow(caption));
  };
  // The Legend rows a History restore returned, shown or deleted in the editor.
  /** @param {DrawingState} drawing @param {string} caption */
  const legendHoldsRow = (drawing, caption) => [
    ...(drawing.legendEntries.value || []), ...(drawing.deletedLegendEntries.value || [])
  ].some((entry) => captionMatches(entry?.caption, caption) || captionMatches(entry?.originalCaption, caption));
  // Undo and Redo of a rule edit restore the rules and, for a step that changed
  // Legend rows, the Result the step left; the composition root passes the
  // rules they replaced, and rows the restored Result cannot show ask for the
  // rerender, as the edit did. A removed whole-type rule asks for it only
  // when the restored Legend lacks the type's default row (OV-169): Undo of
  // an Apply to all returns the row it recolored, so it needs no Python run.
  // The composition root does not call this for a figure History recorded as
  // drawn (D-24 (OV-731)).
  /** @param {Partial<SpecificColorRule>[]} previousRules */
  const followRestoredRules = (previousRules) => {
    const drawing = state.activeDrawing();
    return (
      JSON.stringify(previousRules) !== JSON.stringify(drawing.manualSpecificRules)
      && redrawsLegendFor(drawing, previousRules, drawing.manualSpecificRules,
        { holdsRow: (caption) => legendHoldsRow(drawing, caption) })
      && ports.requestAutomaticRerender()
    );
  };
  // The automatic rerender replaces the Results the candidate was prepared
  // against, which makes the candidate stale (#857). An edit made while one runs
  // waits for it and, when it replaced the Results under an unchanged rule
  // table, prepares again instead of dropping the edit. The rerender is idle
  // before the binder has read the Result it wrote, so a preparation can also
  // go stale on that binding with the Results unchanged (OV-166): any change
  // of the inputs the candidate was prepared from counts, unless the rule
  // table changed, or the caller or a Session operation superseded this one.
  /** @returns {Promise<void>} */
  const rerenderIdle = () => (state.labelReflowProcessing?.value
    ? new Promise((resolve) => {
      const stop = watch(state.labelReflowProcessing, (busy) => { if (!busy) { stop(); resolve(); } });
    })
    : Promise.resolve());
  // OV-61: the rule table is the latest explicit color of the rows its rules
  // draw. A Legend-only color an earlier popup edit left on such a row (the
  // popup stores the caption's color when it recolors a whole group) would
  // win at Generate over the recolored rule, so a commit that recolors the row
  // retires it; `afterCommit` may store a new one.
  // OV-152: a commit that removes a rule row (Clear All, a deleted or
  // recaptioned rule) retires that copy too; left behind, it would paint the
  // row Generate draws for the features' palette color. A copy has the removed
  // rule's caption and color. A color set in the Legend editor is no copy: the
  // editor recolors the rules of a row a rule draws (legend.js), so its colors
  // sit on rows no rule draws and stay.
  /**
   * @param {DrawingState} drawing
   * @param {{ caption: string, color: string }[]} intents
   * @param {{ caption: string, color: string }[]} removedRuleRows
   */
  const retireSupersededLegendColors = (drawing, intents, removedRuleRows) => {
    const overrides = drawing.legendColorOverrides;
    for (const { caption, color } of intents) {
      if (Object.hasOwn(overrides, caption) && resolveColorToHex(overrides[caption]) !== resolveColorToHex(color)) {
        delete overrides[caption];
      }
    }
    for (const { caption, color } of removedRuleRows) {
      if (Object.hasOwn(overrides, caption) && resolveColorToHex(overrides[caption]) === resolveColorToHex(color)) {
        delete overrides[caption];
      }
    }
  };
  // OV-444 (OV-443, OV-445): this owner is the one writer of the rule table,
  // and every write goes through one line (`queueRuleWrite`), in the order it
  // was made. A write builds its next table from the table at its turn
  // (`SpecificRuleBuild`), so no write is lost to another. Only a pending
  // pattern edit can be cancelled: by any newer write, a keystroke, or a
  // drawer close. Its draft stays, and it leaves the line at once
  // (`releaseCancelledWrites`), so no write waits for a cancelled one (OV-448).
  /** @type {Promise<unknown>} */
  let writeLine = Promise.resolve();
  // The rule writes applied so far: a run whose preparation one of them made
  // stale prepares again (`runWithRuleMatches`).
  let appliedWrites = 0;
  /** @type {Set<{ cancelled: () => boolean, leave: () => void }>} */
  const cancellableWrites = new Set();
  /**
   * @template T
   * @param {() => Promise<T>} write Runs at the write's turn.
   * @param {() => boolean} [cancelled] Of a pattern edit: whether it was cancelled before it ended.
   * @returns {Promise<T>}
   */
  const queueRuleWrite = (write, cancelled) => {
    if (!cancelled) patternDrafts.suspend();
    const turn = writeLine.then(write);
    const entry = { cancelled: cancelled || (() => false), leave: () => {} };
    writeLine = Promise.all([writeLine, new Promise((resolve) => { entry.leave = () => resolve(undefined); })]);
    cancellableWrites.add(entry);
    const settle = () => { cancellableWrites.delete(entry); entry.leave(); };
    turn.then(settle, settle);
    return turn;
  };
  const releaseCancelledWrites = () => cancellableWrites.forEach(entry => entry.cancelled() && entry.leave());
  /**
   * `drawing` is the drawing of the action (default: the shown mode's at the
   * call); `isCurrent` is read from the write's turn on. `joinsOpenStep`: the
   * write is part of its caller's open History step and keeps that step's name.
   * @typedef {{ drawing?: DrawingState, isCurrent?: () => boolean, afterCommit?: (intents: Record<string, any>[]) => void, previousLegendIntents?: Record<string, any>[], legendPlacement?: LegendPlacement | null, joinsOpenStep?: boolean }} CommitSpecificRulesOptions
   */
  /**
   * The write port of the rule table (R13): other owners pass how the next
   * table follows from the current one, never a table.
   * @param {SpecificRuleBuild} build
   * @param {string} [label]
   * @param {CommitSpecificRulesOptions} [options]
   */
  const commitSpecificRules = (build, label = 'Change specific color rules', options = {}) => {
    const drawing = options.drawing || state.activeDrawing();
    return queueRuleWrite(() => runSpecificRuleWrite(drawing, build, label, options));
  };
  /**
   * One write at its turn: its table, built from the current one, then its commit.
   * @param {DrawingState} drawing
   * @param {SpecificRuleBuild} build
   * @param {string} label
   * @param {CommitSpecificRulesOptions} commitOptions
   */
  const runSpecificRuleWrite = async (drawing, build, label, commitOptions) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (commitOptions.isCurrent?.() === false) return false;
    await rerenderIdle();
    const current = drawing.manualSpecificRules;
    /** @type {Map<Partial<SpecificColorRule>, SpecificColorRule>} */
    const revised = new Map();
    const rules = build(current, (row, patch) => {
      const next = { ...row, ...patch };
      revised.set(next, row);
      return next;
    });
    if (!rules) return false;
    // The row each rule is written into, so a kept or revised row keeps its identity.
    const sourceRows = rules.map(rule => revised.get(rule) || (current.includes(rule) ? rule : null));
    for (let attempt = 0; ; attempt += 1) {
      const results = state.results.value;
      const table = JSON.stringify(drawing.manualSpecificRules);
      const inputs = rulePreparation.snapshot();
      const outcome = await commitOnce(drawing, rules, sourceRows, label, commitOptions);
      const stale = state.results.value !== results || !rulePreparation.isCurrent(inputs);
      if (outcome || attempt === 2 || !stale
        || JSON.stringify(drawing.manualSpecificRules) !== table
        || commitOptions.isCurrent?.() === false
        || state.sessionOperationAvailability?.()) return outcome;
      await rerenderIdle();
    }
  };
  // R13: this owner builds the candidate rules, so it is the one place that asks
  // the rule preparation for their matches: a commit prepares the candidate
  // rows, a color action prepares the rules it reads as they are.
  /**
   * @param {Partial<SpecificColorRule>[]} rules
   * @param {Parameters<import('../rule-matching.js').RulePreparation['prepareCandidate']>[1]} [options]
   */
  const prepareRules = (rules, options) => rulePreparation.prepareCandidate(rules, options);
  // Runs `commit` once the color rule matches of `rules` are prepared; a
  // prepared table answers at once, so the commit runs in the caller's tick.
  // Rules Python must still match wait for the automatic rerender in flight:
  // it replaces the Result, the catalog, and the Legend rows the preparation
  // must stay current with, so a preparation that spanned it went stale and
  // the run did nothing (OV-280). A rule write that commits meanwhile makes it
  // stale the same way: the run prepares again, so it is not dropped (OV-444).
  /**
   * @param {Record<string, any>[]} rules
   * @param {() => any} commit
   * @returns {unknown}
   */
  const runWithRuleMatches = (rules, commit) => {
    const writes = appliedWrites;
    return runWhenPrepared(state, () => {
      const prepare = () => prepareRules(rules, { asGiven: true });
      return [rulePreparation.isPrepared(rules)
        || (state.labelReflowProcessing?.value ? rerenderIdle().then(prepare) : prepare())];
    }, commit, () => (appliedWrites !== writes ? runWithRuleMatches(rules, commit) : undefined));
  };
  /**
   * @param {DrawingState} drawing
   * @param {Partial<SpecificColorRule>[]} rules
   * @param {(Record<string, any> | null)[]} sourceRows The row each rule is written into, or null for a new row.
   * @param {string} label
   * @param {CommitSpecificRulesOptions} [options]
   */
  const commitOnce = async (drawing, rules, sourceRows, label, { isCurrent = () => true, afterCommit = () => {}, previousLegendIntents = [], legendPlacement = null, joinsOpenStep = false } = {}) => {
    const candidate = await prepareRules(rules);
    if (!candidate) return false;
    const current = () => !state.sessionOperationAvailability?.() && isCurrent()
      && rulePreparation.isCurrent(candidate.snapshot);
    if (!current()) return false;
    const { intents, previousIntents, removedRuleRows, rekeyed } = candidateLegendIntents(drawing, candidate.rules, previousLegendIntents);
    const previousCaptions = new Set(previousIntents.map(intent => intent.caption));
    // The placed row is the intent of its name and color, by the key the
    // commit allocates to it (`X [<fill>]` while another row shows `X`).
    /** @type {PythonLegendKey | undefined} */
    const placedKey = legendPlacement ? intents.find(intent => captionMatches(intent.name ?? intent.caption, legendPlacement.name)
      && resolveColorToHex(intent.color) === resolveColorToHex(legendPlacement.color))?.caption : undefined;
    /** @type {Map<PythonLegendKey, PythonLegendKey>} */
    const replaced = new Map(placedKey && legendPlacement?.replaces && placedKey !== legendPlacement.at ? [[legendPlacement.at, placedKey]] : []);
    // The rows the Legend owner keeps in their place under a new key (OV-645, OV-1235).
    const inPlace = new Map([...rekeyed, ...replaced]);
    const legend = await prepareFileLegendEntries(intents.filter(intent => !(drawing.deletedLegendEntries?.value || [])
      .some((/** @type {LegendEntry} */ entry) => (entry.originalCaption || entry.caption) === intent.caption)), {
      drawing,
      previousFileIntents: previousIntents,
      isCurrent: current,
      placement: placedKey && legendPlacement && !legendPlacement.replaces ? { caption: placedKey, after: legendPlacement.at } : null,
      rekeyed: inPlace
    });
    if (!legend) return false;
    const redrawsLegend = redrawsLegendFor(drawing, [...drawing.manualSpecificRules], candidate.rules, {
      rows: { add: legend.diff.add, retire: legend.diff.remove.map(({ caption }) => caption) },
      replaced: inPlace
    });
    let applied = false;
    // One History step: the rule transition first, then the legend rows it
    // draws (R13); a checkpoint when the legend gains or loses a row.
    const transact = legend.diff.add.length || legend.diff.remove.length
      ? runUndoableCheckpoint : runUndoable;
    // D-24 (OV-731): a checkpoint step records whether the figure of each
    // side is the one its state draws, so Undo and Redo of a drawn side ask
    // for no Python run. The figure the commit starts from is drawn unless a
    // rerender asked for is still owed when History captures it (D-36). The
    // commit's own figure is drawn when it also shows every row; otherwise
    // History keeps the figure of the rerender the commit asks for once that
    // is shown, while the step is the last one (`keepRedrawnResult`).
    await transact(label, () => {
      if (!current() || !legend.isCurrent()) return false;
      drawing.manualSpecificRules.splice(0, drawing.manualSpecificRules.length, ...candidate.rules.map((rule, index) => {
        const row = sourceRows[index];
        if (!row) return rule;
        Object.assign(row, rule);
        if (!Object.hasOwn(rule, 'fromFile')) delete row.fromFile;
        return row;
      }));
      patternDrafts.reconcile();
      drawing.fileLegendCaptions.value = new Set(candidate.rules.filter(rule => rule.fromFile && rule.cap).map(rule => rule.cap));
      drawing.addedLegendCaptions.value = new Set([
        ...[...drawing.addedLegendCaptions.value].filter(caption => !previousCaptions.has(caption)),
        ...intents.map(intent => intent.caption)
      ]);
      refreshFeatureOverrides(extractedFeatures.value);
      retireSupersededLegendColors(drawing, intents, removedRuleRows);
      afterCommit(intents);
      applied = true;
      const legendRows = legend.apply();
      // The step shows its fills and Legend rows once its rules, feature
      // fills, Legend colors, and rows are written: every feature and row fill
      // is an operation of the one compile (OV-146).
      projectPaletteAndRules({ prepareRules: false, legendRows });
      return legend.diff;
    }, { joinsOpenStep, drawn: { before: () => !ports.isRedrawOwed(drawing), change: !redrawsLegend } });
    if (applied) appliedWrites += 1;
    if (applied && redrawsLegend) ports.requestAutomaticRerender();
    return applied;
  };
  /**
   * A write of the rule table's own controls, whose failure the table shows
   * with Retry and Edit.
   * @param {DrawingState} drawing
   * @param {SpecificRuleBuild} build
   * @param {string} label
   * @param {Element | null} [input] The input a failed commit offers for editing.
   */
  const commitPrepared = (drawing, build, label, afterCommit = () => {}, input = null) => queueRuleWrite(async () => {
    const snapshot = rulePreparation.snapshot();
    const previousError = state.errorLog?.value;
    try {
      const applied = await runSpecificRuleWrite(drawing, build, label, { afterCommit });
      if (applied && state.errorLog?.value === ruleFailure.value?.error) state.errorLog.value = null;
      if (applied) ruleFailure.value = null;
      return applied;
    } catch (cause) {
      if (!rulePreparation.isCurrent(snapshot) || state.errorLog?.value !== previousError) return false;
      const error = normalizeUserFacingError(cause, { operation: 'evaluateRules', stage: 'resource-staging' });
      if (state.errorLog) state.errorLog.value = error;
      ruleFailure.value = { error, snapshot, input, retry: () => commitPrepared(drawing, build, label, afterCommit, input) };
      return false;
    }
  });
  /** @type {(row: SpecificColorRule, patch: Partial<SpecificColorRule>) => SpecificRuleBuild} The table with `row` revised, or null when it is gone. */
  const reviseRow = (row, patch) => (current, revise) => (current.includes(row)
    ? current.map(rule => (rule === row ? revise(rule, patch) : rule)) : null);
  /** @param {SpecificColorRule} row @param {string} value */
  const editSpecificRulePattern = (row, value) => {
    const busy = state.sessionOperationAvailability?.();
    return busy || patternDrafts.edit(row, value);
  };
  /** @param {SpecificColorRule} row @param {HTMLInputElement | null} [input] */
  const revertSpecificRulePattern = (row, input = null) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    patternDrafts.revert(row);
    if (input?.isConnected) {
      input.value = row.val;
      input.focus();
    }
  };
  /** @param {DrawingState} drawing @param {SpecificColorRule} row @param {string} value */
  const applySpecificRulePattern = async (drawing, row, value) => {
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    if (!drawing.manualSpecificRules.includes(row) || !isPatternEditAvailable()
      || (state.generatedMode?.value && state.mode?.value !== state.generatedMode.value)) return false;
    patternDrafts.suspend();
    const token = patternDrafts.begin(row, value);
    // OV-1133 (S4-Q5): the edit stays current while its draft does. A Result
    // or catalog that an automatic rerender replaces meanwhile makes its
    // preparation stale, and the write prepares again (`runSpecificRuleWrite`).
    return queueRuleWrite(async () => {
      const current = () => patternDrafts.isCurrent(row, token) && isPatternEditAvailable();
      try {
        const applied = await runSpecificRuleWrite(drawing, reviseRow(row, { val: String(value ?? ''), fromFile: false }),
          'Edit specific color rule', { isCurrent: current });
        if (applied && patternDrafts.isCurrent(row, token)) patternDrafts.revert(row);
        return applied;
      } catch (cause) {
        if (!current() || cause?.canceled || cause?.name === 'AbortError' || ['canceled', 'stale', 'superseded'].includes(cause?.status)) return false;
        // `current()` above holds, so the row has a draft of this token.
        /** @type {PatternDraft} */ (patternDrafts.get(row)).error = normalizeUserFacingError(cause, { operation: 'evaluateRules', stage: 'resource-staging' });
        return false;
      } finally {
        // `isCurrent` holds, so the row has a draft of this token.
        if (patternDrafts.isCurrent(row, token)) /** @type {PatternDraft} */ (patternDrafts.get(row)).pending = false;
      }
    }, () => !patternDrafts.isCurrent(row, token));
  };
  /** @param {SpecificColorRule} row */
  const retrySpecificRulePattern = (row) => {
    const drawing = state.activeDrawing();
    const draft = patternDrafts.get(row);
    return draft && !draft.pending ? applySpecificRulePattern(drawing, row, draft.text) : false;
  };
  /**
   * @param {number} index
   * @param {'feat' | 'qual' | 'val' | 'color' | 'cap'} field
   * @param {string} value
   * @param {HTMLInputElement | HTMLSelectElement | null} [input]
   */
  const setSpecificRuleField = (index, field, value, input = null) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!specificRuleFields.has(field)) return;
    const current = drawing.manualSpecificRules[index];
    if (!current) return;
    if (field === 'val') return applySpecificRulePattern(drawing, current, value);
    const nextValue = field === 'color' ? resolveColorToHex(String(value || '#000000')) : String(value ?? '');
    return commitPrepared(drawing, reviseRow(current, { [field]: nextValue, fromFile: false }), 'Edit specific color rule', () => {}, input).finally(() => {
      if (input?.isConnected && drawing.manualSpecificRules.includes(current)) input.value = current[field] ?? '';
    });
  };

  /** @param {DrawingState} drawing @param {number} index @param {number} offset */
  const moveSpecificRule = (drawing, index, offset) => {
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const row = drawing.manualSpecificRules[index];
    return commitPrepared(drawing, (current) => {
      const rows = [...current];
      const from = rows.indexOf(row);
      const target = from + offset;
      if (from < 0 || target < 0 || target >= rows.length) return null;
      rows.splice(from, 1);
      rows.splice(target, 0, row);
      return rows;
    }, 'Move specific color rule');
  };

  const downloadSpecificRulesTsv = () => {
    const drawing = state.activeDrawing();
    const text = serializeSpecificRules(drawing.manualSpecificRules);
    if (!text.trim()) {
      showNotice({ text: 'There are no color rules to export.' });
      return;
    }
    downloadTextFile('gbdraw_specific_rules.tsv', text);
  };

  /** @param {ColorFeature} feat */
  const getIndividualFeatureLabel = (feat) => {
    return feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${formatFeatureRange(feat)}`;
  };

  // First entry per normalized feature ID (the first-match rule of a linear scan).
  // Entries are edited in place, so callers build it per operation, not once.
  const indexEditableLabels = () => {
    const index = new Map();
    if (!Array.isArray(editableLabels.value)) return index;
    for (const entry of editableLabels.value) {
      const key = normalizeFeatureIdKey(entry?.featureId);
      if (key && !index.has(key)) index.set(key, entry);
    }
    return index;
  };

  /** @param {ColorFeature | null | undefined} feat @param {Map<string, any>} [entryIndex] */
  const getDisplayedFeatureLabel = (feat, entryIndex = indexEditableLabels()) => {
    const drawing = state.activeDrawing();
    if (!feat) return '';

    const featureIdKey = normalizeFeatureIdKey(feat.svg_id || feat.id);
    const editableEntry = featureIdKey ? entryIndex.get(featureIdKey) || null : null;
    const editableText = normalizeCaption(editableEntry?.text);
    if (editableText) return editableText;

    const normalizedOverride = normalizeCaption(featureOverrideValue(drawing.featureOverrides, feat, 'labelText'));
    if (normalizedOverride) return normalizedOverride;

    const sourceText = normalizeCaption(editableEntry?.sourceText);
    if (sourceText) {
      // The source text is label text: only an own key of the table is an override.
      const normalizedBulk = normalizeCaption(ownEntry(drawing.labelTextBulkOverrides, sourceText));
      if (normalizedBulk) return normalizedBulk;
    }

    return normalizeCaption(getIndividualFeatureLabel(feat));
  };

  // The rules a legend row draws, by its Python key (`legendEntryKey`);
  // editing the row edits them (N-06).
  /** @param {PythonLegendKey} key */
  const getLegendRowRules = (key) => legendRowRules(key, displayedLegendRowContext(state, state.activeDrawing()));

  // Resolve the effective legend item label used by current SVG coloring
  // priority: a rule's feature belongs to the row Generate draws for it (N-06).
  // One reader serves a pass over many features (OV-193). It answers "which
  // row" with the row's key; a caption written from it is the row's name
  // (the Legend owner's `legendRowName`, OV-640).
  /** @returns {(feat: Record<string, any> | null | undefined) => PythonLegendKey} */
  const effectiveLegendCaptions = () => {
    const drawing = state.activeDrawing();
    const ruleMatches = ruleMatcher(drawing.manualSpecificRules);
    const legendCaption = ruleLegendCaptions(displayedLegendRowContext(state, drawing));
    /** @param {unknown} caption A row whose key is its caption. */
    const rowOf = (caption) => /** @type {PythonLegendKey} */ (normalizeCaption(caption));
    return (feat) => {
      if (!feat) return rowOf('');
      const rule = ruleMatches.first(feat);
      if (rule && normalizeCaption(rule.cap)) return legendCaption(rule);

      // A color edit carries a name (its rule's caption): its row is the one
      // the rules of that caption and color draw, else the row it names.
      const override = getFeatureOverride(drawing.featureColorOverrides, feat);
      const overrideCaption = normalizeCaption(override?.caption);
      if (overrideCaption) return legendCaption({ cap: overrideCaption, color: override.color });

      return rowOf(normalizeCaption(feat.type) || getIndividualFeatureLabel(feat));
    };
  };
  /** @param {ColorFeature | null | undefined} feat */
  const getEffectiveLegendCaption = (feat) => effectiveLegendCaptions()(feat);

  const addCustomColor = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newColorFeat.value) return;
    drawing.currentColors.value = {
      ...drawing.currentColors.value,
      [newColorFeat.value]: newColorVal.value
    };
  };

  /** @param {string} value */
  const setLabelFilterMode = (value) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.filterMode.value = value;
  };
  const addWhitelistRule = () => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.manualWhitelist.push({ feat: 'CDS', qual: 'product', key: '' });
  };
  /** @param {number} index */
  const removeWhitelistRule = (index) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.manualWhitelist.splice(index, 1);
  };
  /** @param {number} index */
  const removePriorityRule = (index) => {
    const drawing = state.activeDrawing();
    const busy = state.sessionOperationAvailability?.();
    if (busy) return busy;
    drawing.manualPriorityRules.splice(index, 1);
  };

  const addPriorityRule = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newPriorityRule.order) return;
    const idx = drawing.manualPriorityRules.findIndex((/** @type {{ feat?: string, order?: string }} */ r) => r.feat === newPriorityRule.feat);
    if (idx >= 0) {
      drawing.manualPriorityRules[idx].order = newPriorityRule.order;
    } else {
      drawing.manualPriorityRules.push({ feat: newPriorityRule.feat, order: newPriorityRule.order });
    }
  };

  const addSpecificRule = async () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!newSpecRule.val) return;

    // One question names every reason the pattern may match slowly.
    const pattern = String(newSpecRule.val);
    const risks = [
      pattern.length > 50 ? 'It is longer than 50 characters.' : '',
      /\(.+[\+\*]\)[\+\*]/.test(pattern) || /\(.*\)\+/.test(pattern)
        ? 'It repeats a group, as in (a+)+, which can make the page stop responding.'
        : ''
    ].filter(Boolean);
    if (risks.length > 0) {
      const choice = await askChoice({
        title: 'Add this color rule?',
        description: 'This pattern may make color matching slow:',
        items: risks,
        actions: [{ id: 'add', label: 'Add rule', tone: 'primary' }],
        cancelLabel: 'Edit pattern'
      });
      if (choice !== 'add' || newSpecRule.val !== pattern) return;
      const busy = state.sessionOperationAvailability?.();
      if (busy) return busy;
    }

    const rule = {
      feat: String(newSpecRule.feat || ''),
      qual: String(newSpecRule.qual || ''),
      val: String(newSpecRule.val),
      color: String(newSpecRule.color || '#000000'),
      cap: String(newSpecRule.cap || '')
    };
    return commitPrepared(drawing, (current) => [...current, rule], 'Add specific color rule', () => {
      if (newSpecRule.val === rule.val) newSpecRule.val = '';
    }, document.querySelector?.('[data-new-specific-rule-pattern]'));
  };

  // A rules TSV replaces the rules an earlier file added and keeps the others.
  // The file is read at once, so a bad file fails before the write waits.
  /**
   * @param {string} text
   * @param {{ drawing?: DrawingState, isCurrent?: () => boolean }} [options]
   * @returns {Promise<number | false>} The count of imported rules, or false.
   */
  const importSpecificRules = async (text, { drawing = state.activeDrawing(), isCurrent = () => true } = {}) => {
    const { importedCount } = prepareSpecificColorImport(text);
    const applied = await commitSpecificRules((current) => prepareSpecificColorImport(text, [...current]).nextRules,
      'Import specific color rules', { drawing, isCurrent });
    return applied === true ? importedCount : false;
  };

  const applySpecificRulePreset = async () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (specificRulePresetLoading.value) return;
    const presetId = selectedSpecificPreset.value;
    if (!presetId) return;
    const preset = specificRulePresets.find((/** @type {{ id: string, path: string }} */ entry) => entry.id === presetId);
    if (!preset) {
      showNotice({ text: 'That rule preset is not available.', tone: 'warning' });
      return;
    }

    specificRulePresetLoading.value = true;
    // The preset takes its place in the write line when it is chosen.
    return queueRuleWrite(async () => {
      const presetContext = rulePreparation.snapshot();
      const previousError = state.errorLog?.value;
      try {
        const response = await fetch(preset.path, { cache: 'no-store' });
        if (!response.ok) {
          throw { code: 'INPUT_UNREADABLE', operation: 'evaluateRules', stage: 'resource-staging' };
        }
        const text = await response.text();
        const { rules } = parseSpecificRules(text);
        if (selectedSpecificPreset.value !== presetId || !rulePreparation.isCurrent(presetContext)) return;
        return await runSpecificRuleWrite(drawing, () => rules, 'Apply specific color preset', {
          isCurrent: () => selectedSpecificPreset.value === presetId,
          afterCommit: () => {
            if (presetId === 'bakta') {
              drawing.currentColors.value = { ...drawing.currentColors.value, CDS: '#cccccc' };
              drawing.adv.legend_box_size = 12;
              drawing.adv.legend_font_size = 12;
            }
          }
        });
      } catch (e) {
        if (selectedSpecificPreset.value !== presetId || !rulePreparation.isCurrent(presetContext)
          || state.errorLog?.value !== previousError) return;
        const error = normalizeUserFacingError(e, { operation: 'evaluateRules', stage: 'resource-staging' });
        if (state.errorLog) state.errorLog.value = error;
        ruleFailure.value = { error, snapshot: presetContext, presetId, retry: applySpecificRulePreset };
      } finally {
        specificRulePresetLoading.value = false;
      }
    });
  };

  const addFeature = () => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (newFeatureToAdd.value && !drawing.adv.features.includes(newFeatureToAdd.value)) {
      drawing.adv.features.push(newFeatureToAdd.value);
      if (!drawing.adv.feature_shapes || typeof drawing.adv.feature_shapes !== 'object') {
        drawing.adv.feature_shapes = {};
      }
      if (!Object.prototype.hasOwnProperty.call(drawing.adv.feature_shapes, newFeatureToAdd.value)) {
        drawing.adv.feature_shapes[newFeatureToAdd.value] = defaultFeatureRendering(newFeatureToAdd.value);
      }
    }
  };

  /** @param {string} featureType */
  const removeFeature = (featureType) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    const idx = drawing.adv.features.indexOf(featureType);
    if (idx >= 0) {
      drawing.adv.features.splice(idx, 1);
    }
  };

  /** @param {string} featureType */
  const getFeatureShape = (featureType) => {
    const drawing = state.activeDrawing();
    if (!drawing.adv.feature_shapes || typeof drawing.adv.feature_shapes !== 'object') {
      return defaultFeatureRendering(featureType);
    }
    return hasOwnEntry(drawing.adv.feature_shapes, featureType)
      ? normalizeFeatureRendering(drawing.adv.feature_shapes[featureType])
      : defaultFeatureRendering(featureType);
  };

  /** @param {string} featureType @param {string} shape */
  const setFeatureShape = (featureType, shape) => {
    const drawing = state.activeDrawing();
    const sessionBusy = state.sessionOperationAvailability?.();
    if (sessionBusy) return sessionBusy;
    if (!drawing.adv.feature_shapes || typeof drawing.adv.feature_shapes !== 'object') {
      drawing.adv.feature_shapes = {};
    }
    drawing.adv.feature_shapes[featureType] = normalizeFeatureRendering(shape);
  };

  /** @param {ColorFeature} feat */
  const getFeatureColor = (feat) => {
    const drawing = state.activeDrawing();
    const override = getFeatureOverride(drawing.featureColorOverrides, feat);
    if (override) {
      return resolveColorToHex(override.color || override);
    }
    return resolveColorToHex(appliedFeatureColorOf(state, feat.type)) || '#cccccc';
  };

  /** @param {ColorFeature} feat */
  const getFeatureColorValue = (feat) => {
    const drawing = state.activeDrawing();
    const override = getFeatureOverride(drawing.featureColorOverrides, feat);
    if (!override) return null;
    return Object.prototype.hasOwnProperty.call(override, 'color')
      ? override.color
      : override;
  };

  const canEditFeatureColor = () => true;

  // Python matches a single-feature color rule only by the stable generation
  // hash, so duplicate record instances that share one hash share the rule.
  /** @param {ColorFeature} feat */
  const getFeatureQualifier = (feat) => {
    const colorRuleHash = getFeatureColorRuleHash(feat);
    return colorRuleHash ? { qual: 'hash', val: colorRuleHash } : null;
  };

  /** @param {ColorFeature | null | undefined} feat @param {string} label */
  const getLabelSpecificRule = (feat, label) => {
    const drawing = state.activeDrawing();
    if (!feat) return null;
    const priorityRule = drawing.manualPriorityRules.find((/** @type {{ feat?: string, order?: string }} */ rule) => rule.feat === feat.type);
    const priority = String(priorityRule?.order || '')
      .split(',')
      .map((qualifier) => qualifier.trim())
      .filter(Boolean);
    const selector = resolveFeatureLabelSelector(feat, label, { priority });
    if (!selector) return null;
    return {
      feat: feat.type,
      qual: selector.qualifier,
      val: selector.pattern
    };
  };

  /** @param {ColorFeature[] | null | undefined} features */
  const refreshFeatureOverrides = (features) => {
    const drawing = state.activeDrawing();
    const ruleMatches = ruleMatcher(drawing.manualSpecificRules);
    if (!features || features.length === 0 || !ruleMatches.ready(features)) return;
    migrateLegacyFeatureOverrides(drawing.featureColorOverrides, features);

    for (const feat of features) {
      const rule = ruleMatches.first(feat);
      const key = featureOverrideKey(feat);
      if (key) {
        if (rule) drawing.featureColorOverrides[key] = { color: rule.color, caption: rule.cap };
        else delete drawing.featureColorOverrides[key];
      }
    }
  };

  /** @param {ColorFeature} feat */
  const findMatchingRegexRule = (feat) => {
    const drawing = state.activeDrawing();
    return ruleMatcher(drawing.manualSpecificRules.filter((/** @type {SpecificColorRule} */ rule) => rule.qual !== 'hash')).first(feat);
  };

  /** @param {Partial<SpecificColorRule> | null | undefined} rule */
  const countFeaturesMatchingRule = (rule) => {
    if (!rule || String(rule.qual || '').toLowerCase() === 'hash') return 0;

    let count = 0;
    for (const feat of extractedFeatures.value) {
      if (rule.feat !== '*' && feat.type !== rule.feat) continue;

      if (ruleMatchesFeature(feat, rule)) count++;
    }
    return count;
  };

  /** @param {ColorFeature} currentFeat @param {string | null} [caption] */
  const findFeaturesWithSameLegendItem = (currentFeat, caption = null) => {
    const captionOf = effectiveLegendCaptions();
    const targetCaption = normalizeCaption(caption || captionOf(currentFeat));
    if (!targetCaption) return [];
    return /** @type {ColorFeature[]} */ (extractedFeatures.value).filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(captionOf(f), targetCaption);
    });
  };

  /** @param {ColorFeature} currentFeat @param {string | null} [label] */
  const findFeaturesWithSameIndividualLabel = (currentFeat, label = null) => {
    const targetLabel = normalizeCaption(label || getIndividualFeatureLabel(currentFeat));
    if (!targetLabel) return [];

    return /** @type {ColorFeature[]} */ (extractedFeatures.value).filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getIndividualFeatureLabel(f), targetLabel);
    });
  };

  /** @param {ColorFeature} currentFeat @param {string | null} [label] */
  const findFeaturesWithSameDisplayedLabel = (currentFeat, label = null) => {
    const entryIndex = indexEditableLabels();
    const targetLabel = normalizeCaption(label || getDisplayedFeatureLabel(currentFeat, entryIndex));
    if (!targetLabel) return [];

    return /** @type {ColorFeature[]} */ (extractedFeatures.value).filter((f) => {
      if (f.svg_id === currentFeat.svg_id) return false;
      return captionMatches(getDisplayedFeatureLabel(f, entryIndex), targetLabel);
    });
  };

  /** @param {ColorFeature} _currentFeat @param {string} caption */
  const findExistingColorForCaption = (_currentFeat, caption) => {
    const drawing = state.activeDrawing();
    const targetCaption = normalizeCaption(caption);
    if (!targetCaption) return null;

    for (const rule of drawing.manualSpecificRules) {
      if (captionMatches(rule.cap, targetCaption) && String(rule.qual || '').toLowerCase() === 'hash') {
        return { rule, color: rule.color };
      }
    }

    for (const rule of drawing.manualSpecificRules) {
      if (captionMatches(rule.cap, targetCaption) && String(rule.qual || '').toLowerCase() !== 'hash') {
        return { rule, color: rule.color };
      }
    }

    return null;
  };

  return {
    specificRulePattern: patternDrafts.text,
    specificRulePatternDraft: patternDrafts.get,
    specificRulePatternFieldId: patternDrafts.fieldId,
    editSpecificRulePattern, retrySpecificRulePattern, revertSpecificRulePattern,
    suspendSpecificRulePatternDrafts: patternDrafts.suspend,
    clearSpecificRulePatternDrafts: patternDrafts.clear,
    captureSpecificRulePatternDrafts: patternDrafts.capture,
    restoreSpecificRulePatternDrafts: patternDrafts.restore,
    canRetrySpecificRuleFailure, canEditSpecificRuleFailure, retrySpecificRuleFailure, editSpecificRuleFailure,
    addCustomColor,
    addFeature,
    removeFeature,
    getFeatureShape,
    setFeatureShape,
    addPriorityRule,
    setLabelFilterMode,
    addWhitelistRule,
    removeWhitelistRule,
    removePriorityRule,
    addSpecificRule,
    applySpecificRulePreset,
    canEditFeatureColor,
    clearAllSpecificRules: () => {
      const drawing = state.activeDrawing();
      return commitPrepared(drawing, () => [], 'Clear specific color rules');
    },
    commitSpecificRules,
    importSpecificRules,
    followRestoredRules,
    countFeaturesMatchingRule,
    downloadSpecificRulesTsv,
    findExistingColorForCaption,
    findFeaturesWithSameCaption: findFeaturesWithSameLegendItem,
    findFeaturesWithSameDisplayedLabel,
    findFeaturesWithSameIndividualLabel,
    findFeaturesWithSameLegendItem,
    findMatchingRegexRule,
    getFeatureColor,
    getFeatureColorValue,
    getDisplayedFeatureLabel,
    effectiveLegendCaptions,
    getEffectiveLegendCaption,
    getIndividualFeatureLabel,
    getFeatureQualifier,
    getLabelSpecificRule,
    getLegendRowRules,
    /** @param {number} index */
    moveSpecificRuleDown: (index) => {
      const drawing = state.activeDrawing();
      return moveSpecificRule(drawing, index, 1);
    },
    /** @param {number} index */
    moveSpecificRuleUp: (index) => {
      const drawing = state.activeDrawing();
      return moveSpecificRule(drawing, index, -1);
    },
    runWithRuleMatches,
    /** @param {number} index */
    removeSpecificRule: (index) => {
      const drawing = state.activeDrawing();
      const row = drawing.manualSpecificRules[index];
      return commitPrepared(drawing, (current) => (current.includes(row) ? current.filter(rule => rule !== row) : null),
        'Remove specific color rule');
    },
    setSpecificRuleField
  };
};
