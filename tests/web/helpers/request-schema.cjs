// Promote a canonical render request to the current schema with the app's own
// promoter, so Save/Load comparisons hold against requests captured before the
// session writer upgraded them.
const promoteRequest = (page, request) => page.evaluate(async value => {
  const { promoteCanonicalRenderRequestToCurrent } = await import('./js/services/session-request.js');
  return promoteCanonicalRenderRequestToCurrent(value);
}, request);

module.exports = { promoteRequest };
