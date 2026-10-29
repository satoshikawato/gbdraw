// @ts-check

import { defaultDepthSeriesLabel } from '../services/depth-track-state.js';

/** @import { DepthSlotFile } from '../services/depth-track-state.js' */

/** @param {DepthSlotFile | null | undefined} file */
export const getDepthTrackFileBaseName = (file) => {
  // split always returns at least one element, so pop is not undefined.
  const rawName = /** @type {string} */ (String(file?.name || '').split(/[\\/]/).pop()).trim();
  if (!rawName) return '';
  const withoutExtension = rawName.replace(/\.[^.]+$/, '').trim();
  return withoutExtension || rawName;
};

/**
 * @param {unknown} label
 * @param {number} index
 * @param {DepthSlotFile | null} [file]
 */
export const isDepthTrackAutoLabel = (label, index, file = null) => {
  const normalized = String(label ?? '').trim();
  if (!normalized) return true;
  // A default name the series has for some series count: one a Session saved
  // before OV-1508 stored (Depth for the first of several series included).
  if ([index + 1, index + 2].some((count) => normalized === defaultDepthSeriesLabel(index, count))) return true;
  const fileBaseName = getDepthTrackFileBaseName(file);
  return Boolean(fileBaseName && normalized === fileBaseName);
};

/**
 * Whether a series label is still the drawing's auto label: a default name, or
 * the base name of a file the series had in that drawing before a change. The
 * drawing's own files decide it, so a label stays with its drawing through an
 * Undo, a Reset, or a Load (OV-2202).
 * @param {unknown} label
 * @param {number} index
 * @param {readonly (DepthSlotFile | null | undefined)[]} files The series' files before the change
 */
export const isDepthSeriesAutoLabel = (label, index, files) => (
  isDepthTrackAutoLabel(label, index) || files.some((file) => Boolean(file) && isDepthTrackAutoLabel(label, index, file))
);
