// Ambient globals for the R14 typed-boundaries check (gbdraw/web/CLAUDE.md).
// This is the only TypeScript-syntax file of the check and lives outside the
// shipped tree. It declares only what vendored scripts install and the test
// hooks; every other type is JSDoc in the module that owns it.
export {};

declare global {
  interface Window {
    // Vendored scripts loaded by index.html.
    Vue: any;
    jspdf: any;
    DOMPurify: any;
    // Test and debugging hooks (`__GBDRAW_HISTORY__`, `__GBDRAW_APP__`).
    [hook: `__GBDRAW_${string}`]: any;
    // The pairwise-match colors of the last Generate (`app/run-analysis.js`
    // writes them, `app/svg-styles.js` reads them).
    _origPairwiseMin: string;
    _origPairwiseMax: string;
  }
  var DOMPurify: any;
  var __GBDRAW_TEST_HOOKS__: any;
  var __GBDRAW_LAST_LOSAT_TELEMETRY__: any;
  var __GBDRAW_LOSAT_EXECUTOR__: any;
}
