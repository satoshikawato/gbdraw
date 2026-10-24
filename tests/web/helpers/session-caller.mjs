// Save (`exportSession`) and Load (`importSession`) take the shown drawing
// from their caller, as app-setup passes it (`shownDrawing`); a test calls
// them through this helper, which passes `state.activeDrawing`. Node tests
// and Playwright page code (`import('/tests/web/helpers/session-caller.mjs')`)
// both use it, so it imports nothing. It returns the module's exports with
// those two wrapped.
/**
 * @param {{ exportSession: Function, importSession: Function }} config The
 *   services/config.js module the test loaded.
 * @param {{ activeDrawing: () => unknown }} state The state.js `state` of
 *   that module.
 */
export const callerSession = (config, state) => ({
  ...config,
  exportSession: (title = null, options = {}) => config.exportSession(title, {
    shownDrawing: state.activeDrawing, ...options
  }),
  importSession: (event, options) => config.importSession(event, {
    shownDrawing: state.activeDrawing, ...options
  })
});
