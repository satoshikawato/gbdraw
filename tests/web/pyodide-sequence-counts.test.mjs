// gbdraw's sequence counts run in the vendored Pyodide, where numpy's default
// integer and intp are 32 bits, as they do in the Worker. CPython on a 64-bit
// machine cannot show a value that wraps or refuses at 2^31 there: two such
// defects reached a branch (OV-3356: reduceat refused int64 row marks; R4-1:
// a run's chunk offsets past 2 GiB, the T2T genome's sequence, refused or
// wrapped). The checks are in helpers/pyodide_sequence_counts.py.
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { loadVendoredGbdraw } from './helpers/pyodide-high-heap.mjs';

const pyodide = await loadVendoredGbdraw();

test('sequence counts and chunk offsets are exact in Pyodide past 2 GiB', () => {
  pyodide.runPython(readFileSync(new URL('./helpers/pyodide_sequence_counts.py', import.meta.url), 'utf8'));
});
