// @ts-check
// What a Linear comparison-plan edge holds of the user's own input, read by the plan's edits
// (services/linear-comparisons.js, app/linear-comparison-edits.js) and by Reset Settings (the settings slice,
// services/drawing/settings.js). It imports nothing, so the slice may import it.
/** @import { LinearComparisonEdge, WindowLosatFilenames } from './linear-comparisons.js' */

/**
 * The window-pair names of an edge (OV-2742); an edge without any holds none.
 * @param {{ windowLosatFilenames?: WindowLosatFilenames }} edge
 * @returns {WindowLosatFilenames}
 */
export const windowLosatFilenamesOf = (edge) => edge.windowLosatFilenames ?? {};

/**
 * Whether an edge the plan leaves out still holds what the user gave it (its BLAST TSV, a typed LOSAT name of the
 * record pair or of one of its window pairs), so the plan keeps it: the one rule the plan's edits and Reset Settings
 * keep edges by.
 * @param {Pick<LinearComparisonEdge, 'file' | 'losatFilename' | 'windowLosatFilenames'>} edge
 */
export const edgeRetainsUserData = (edge) => Boolean(
  edge.file || edge.losatFilename.trim() || Object.keys(windowLosatFilenamesOf(edge)).length > 0
);
