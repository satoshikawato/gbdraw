// @ts-check
// What the header and the error notice show: the Session title as typed and
// as the header names it, the error as the notice shows it, whether Save
// needs one Generate first, and the reload the notice offers. The
// composition root passes the state as a port.
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { t } from '../utils/i18n.js';
/** @import { WebState } from '../state.js' */

const { computed } = window.Vue;

/**
 * Registers the `errorDisplay`, `sessionSaveNeedsGenerate` and
 * `sessionTitleLabel` computeds; the root calls it where they stood.
 * @param {{ state: WebState }} ports
 */
export const createSessionStatus = ({ state }) => {
  const { errorLog, sessionTitle, results } = state;
  /** @param {unknown} value */
  const normalizeSessionTitle = (value) => {
    if (value === null || value === undefined) return '';
    return String(value).trim();
  };

  const errorDisplay = computed(() => normalizeUserFacingError(errorLog.value));
  // D-25: a Result without current feature metadata (a Session older than 40)
  // is saved only after one Generate; the notice comes from the same state.
  const sessionSaveNeedsGenerate = computed(() => results.value.length > 0 && !state.featureCatalog.value);
  const reloadAfterOperationError = () => window.location.reload();

  const sessionTitleLabel = computed(() => {
    const title = normalizeSessionTitle(sessionTitle.value);
    return title || t('Untitled session');
  });
  return { normalizeSessionTitle, errorDisplay, sessionSaveNeedsGenerate, reloadAfterOperationError, sessionTitleLabel };
};
