// @ts-check
// Copy to another drawing into an existing drawing (PD-OI-090 〔作成〕(4)): the drawings the dialog offers, the payload
// that sends the target drawing's draft request with the source's, the lines of a paste's reply, and the last target of
// the open Session.
import { alreadyPresentLine, leftOutLine, pasteNotCarriedLine } from './copy-text.js';

/** @import { DrawingTab, DrawingState } from '../../state.js' */
/** @import { DrawingId } from '../../services/drawing-view.js' */
/** @import { PasteDerivation, ProjectedDraftRequest } from '../../services/config.js' */
/** @import { PastedRecords } from '../../contracts/locate.generated.js' */
/** @import { PageExistingDrawingCopy, RegionDrawingHelperPayload } from '../../services/diagram-worker-protocol.js' */

/**
 * What the dialog reads and does of the existing drawings.
 * @typedef {object} PasteTargetPorts
 * @property {() => readonly DrawingTab[]} tabs The drawings of the Session, in tab order.
 * @property {(drawingId: DrawingId) => DrawingState} drawingOf A drawing of the Session, by its ID.
 * @property {(drawingId: DrawingId) => Promise<ProjectedDraftRequest>} projectDrawing A drawing's draft request, shown
 *   or not, its records read first.
 * @property {(drawingId: DrawingId, paste: PasteDerivation) => Promise<unknown>} pasteInto The paste, as one History
 *   step of that drawing, which is shown.
 */

/** The target of a copy: a new drawing, or an existing drawing by its ID. @typedef {'new' | DrawingId} CopyTarget */

/**
 * The drawings a copy can paste into: every drawing but the shown one (D-18), in tab order.
 * @param {readonly DrawingTab[]} tabs @param {DrawingId} shown
 * @returns {{ id: DrawingId, name: string }[]}
 */
export const pasteTargetsOf = (tabs, shown) => tabs.filter(({ id }) => id !== shown).map(({ id, name }) => ({ id, name }));

/**
 * The payload parts of a paste: the one resource table (the source's resources, then the target's under staging IDs
 * of their own, `paste-target-<n>`) and the copy that names them (`PageExistingDrawingCopy`).
 * @param {RegionDrawingHelperPayload['resources']} sourceResources
 * @param {ProjectedDraftRequest} targetSent
 * @param {Pick<PageExistingDrawingCopy, 'edits' | 'duplicates'>} choices
 * @returns {{ resources: RegionDrawingHelperPayload['resources'], copy: PageExistingDrawingCopy }}
 */
export const pastePayloadOf = (sourceResources, targetSent, choices) => {
  /** @type {Record<string, unknown>} */
  const resources = { ...sourceResources };
  /** @type {Record<string, string>} */
  const resourceIds = {};
  let next = 0;
  for (const [resourceId, descriptor] of Object.entries(targetSent.envelope.resources)) {
    let stagingId;
    do { next += 1; stagingId = `paste-target-${next}`; } while (Object.hasOwn(resources, stagingId));
    resources[stagingId] = descriptor;
    resourceIds[resourceId] = stagingId;
  }
  return {
    resources: /** @type {RegionDrawingHelperPayload['resources']} */ (resources),
    copy: { into: 'existing', ...choices, target: { renderRequest: targetSent.envelope.renderRequest, resourceIds } }
  };
};

/**
 * The dialog's lines of a paste: what is already in the drawing, what it leaves out (with the reason), and what it
 * does not bring. Records are named by the captions of the dialog's rows (by record key).
 * @param {PastedRecords} pasted @param {(recordKey: string) => string} captionOf @param {string} drawingName
 */
export const pasteLinesOf = ({ alreadyPresent, leftOut, notCarried }, captionOf, drawingName) => ({
  alreadyPresent: alreadyPresent.map((item) => alreadyPresentLine(item, captionOf, drawingName)),
  leftOut: leftOut.map((item) => leftOutLine(item, captionOf)),
  notCarried: notCarried.map(pasteNotCarriedLine)
});

/**
 * The last target of the open Session (not saved in it): a new drawing at first; after a copy, the target chosen.
 * A remembered drawing that is no longer the same drawing (removed, or replaced by a Load or New Session) gives a new
 * drawing.
 * @param {PasteTargetPorts['drawingOf']} drawingOf
 */
export const createLastTarget = (drawingOf) => {
  /** @type {{ id: DrawingId, drawing: DrawingState } | null} */
  let last = null;
  return {
    /** @param {readonly { id: DrawingId }[]} targets The drawings offered. @returns {CopyTarget} */
    read: (targets) => {
      const kept = last;
      return kept && targets.some(({ id }) => id === kept.id) && drawingOf(kept.id) === kept.drawing ? kept.id : 'new';
    },
    /** @param {CopyTarget} target */
    remember: (target) => {
      last = target === 'new' ? null : { id: target, drawing: drawingOf(target) };
    }
  };
};
