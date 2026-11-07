// @ts-check
// The preserved settings of the shown drawing: the config overrides a loaded
// Session or request set that no control of the form manages, as the
// Preserved settings list shows them, and their Reset, one History step each.
// The composition root passes the state and History's step runner as ports.
import { msg } from '../utils/i18n.js';
/** @import { WebState } from '../state.js' */

const { computed } = window.Vue;

/**
 * @typedef {object} PreservedConfigOverridePorts
 * @property {WebState} state
 * @property {(label: string, fn: () => unknown) => Promise<unknown>} runUndoable History's undoable step.
 */

/**
 * Registers the `unmanagedConfigOverrideEntries` computed; the root calls it where it stood.
 * @param {PreservedConfigOverridePorts} ports
 */
export const createPreservedConfigOverrides = ({ state, runUndoable }) => {
  /** @param {unknown} value */
  const displayPreservedConfigValue = (value) => {
    const serialized = JSON.stringify(value) ?? String(value);
    return serialized.length <= 240
      ? serialized
      : `${serialized.slice(0, 237)}...`;
  };
  const unmanagedConfigOverrideEntries = computed(() => (
    Object.entries(state.activeDrawing().unmanagedConfigOverrides)
      .sort(([left], [right]) => left.localeCompare(right))
      .map(([path, value]) => ({
        path,
        value,
        displayValue: displayPreservedConfigValue(value)
      }))
  ));
  /** @param {string} path */
  const resetUnmanagedConfigOverride = (path) => {
    const drawing = state.activeDrawing();
    return runUndoable(
      msg('Reset preserved setting'),
      () => {
        if (!Object.prototype.hasOwnProperty.call(drawing.unmanagedConfigOverrides, path)) {
          return false;
        }
        delete drawing.unmanagedConfigOverrides[path];
        return true;
      }
    );
  };
  return { unmanagedConfigOverrideEntries, resetUnmanagedConfigOverride };
};
