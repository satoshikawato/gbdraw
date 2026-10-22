// Module hook (tests only): records each verdict of R4's one owner,
// `restoresNothing` in services/session-authority.js, in
// `globalThis.restoresNothingVerdicts`, so a test can tell an R4 refusal from
// another refusal with the same diagnostic (E-OD22). Register it with
// `register('./helpers/restores-nothing-spy.mjs', import.meta.url)` before the
// app modules are imported.
const NEEDLE = 'export const restoresNothing = (stored, repairs) => {';
export async function load(url, context, nextLoad) {
  const result = await nextLoad(url, context);
  if (!url.endsWith('/gbdraw/web/js/services/session-authority.js')) return result;
  const source = String(result.source);
  if (!source.includes(NEEDLE)) throw new Error('restores-nothing spy: restoresNothing not found');
  return {
    ...result,
    shortCircuit: true,
    source: source.replace(NEEDLE, 'export const restoresNothing = (stored, repairs) => {\n'
      + '  const verdict = unspiedRestoresNothing(stored, repairs);\n'
      + '  globalThis.restoresNothingVerdicts?.push(verdict);\n'
      + '  return verdict;\n'
      + '};\n'
      + 'const unspiedRestoresNothing = (stored, repairs) => {')
  };
}
