// Reads session-authority.js after the writer flip (the generated
// SESSION_VERSION is the drawings layout's, 50) and prints the band answers
// the flip must leave as they are, as JSON. A test runs this file in a child
// process (tests/web/session-authority.test.mjs).
const real = await import(new URL('../../../gbdraw/web/js/contracts/schema.generated.js', import.meta.url).href);
const { isSettingsOnlySessionDocument, projectsWithoutRequest, validateSessionAuthorityInventory } = await import(
  new URL('../../../gbdraw/web/js/services/session-authority.js', import.meta.url).href
);

// A settings-only Session of `version`, or a Session 50 drawing's view
// (`drawingViewOf`) of a drawing that holds no request.
/** @param {number} version */
const requestLess = (version) => ({
  format: 'gbdraw-session', version, resources: {}, renderRequest: null, results: [],
  editorState: { featureCatalog: null }, ui: { mode: 'circular' }, modes: { circular: { config: { form: {}, adv: {} } } },
  webFiles: { bindings: { schema: 2, c_gb: null, c_gff: null, c_fasta: null, linearSeqs: [] } }
});
/** @param {number} version */
const bindingsRead = (version) => {
  try {
    validateSessionAuthorityInventory(requestLess(version), version, {}, { load: true });
    return 'read';
  } catch (error) {
    return String(/** @type {Error} */ (error).message);
  }
};
console.log(JSON.stringify({
  sessionVersion: real.SESSION_DRAWINGS_LAYOUT.version,
  settingsOnly: Object.fromEntries([46, 47, 48, 49, 50].map((version) => [version, isSettingsOnlySessionDocument(requestLess(version))])),
  bindings: Object.fromEntries([46, 47, 48, 49, 50].map((version) => [version, bindingsRead(version)])),
  // Load's migration of a Session 27-39 names the writer's version (50), and a
  // Session 27-30 holds no request: Load reads it by the file's band.
  projectsWithoutRequest: {
    30: projectsWithoutRequest({ format: 'gbdraw-session', version: 50 }, 30),
    48: projectsWithoutRequest(requestLess(48), 48),
    49: projectsWithoutRequest(requestLess(49), 49),
    50: projectsWithoutRequest(requestLess(50), 50)
  }
}));
