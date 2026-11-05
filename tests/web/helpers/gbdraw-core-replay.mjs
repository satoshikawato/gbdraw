// Replays the engine calls that tests/test_native_ribbons_differential.py recorded
// through the Wasm engine, and compares each result with the native engine's.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { createGbdrawCore } from '../../../gbdraw/web/js/workers/gbdraw-core.js';

const [wasmPath, casesPath, contract] = process.argv.slice(2);
const core = await createGbdrawCore(readFileSync(wasmPath));
if (core.CONTRACT_VERSION !== Number(contract)) {
  throw new Error(`The Wasm engine has contract ${core.CONTRACT_VERSION}, not ${contract}. Rebuild it.`);
}
const cases = JSON.parse(readFileSync(casesPath, 'utf8'));
const identical = cases.filter(({ params, columns, sha256 }) => {
  const bytes = Buffer.from(columns, 'hex');
  const text = core.pairwise_ribbons(params, { toJs: () => bytes, destroy() {} });
  return createHash('sha256').update(text).digest('hex') === sha256;
}).length;
console.log(`${identical} of ${cases.length} identical`);
