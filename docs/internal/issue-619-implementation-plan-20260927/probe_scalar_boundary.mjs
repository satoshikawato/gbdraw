import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { execFileSync } from 'node:child_process';

const repoRoot = resolve(process.cwd());
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-issue-619-boundary-'));
try {
  await cp(join(repoRoot, 'gbdraw/web/js'), join(tempRoot, 'js'), { recursive: true });
  await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}');
  const slotsModule = await import(pathToFileURL(join(tempRoot, 'js/app/circular-track-slots.js')));
  const contract = await import(pathToFileURL(join(tempRoot, 'js/services/session-active-config-contract.js')));
  const cases = [
    ['typed_px', { value: 20, unit: 'px' }, { value: 20, unit: 'px' }],
    ['typed_factor', { value: 0.08, unit: 'factor' }, { value: 0.08, unit: 'factor' }],
    ['typed_numeric_text', { value: '1.', unit: 'px' }, { value: 1, unit: 'px' }],
    ['typed_exponent_text', { value: '1e-3', unit: 'factor' }, { value: 0.001, unit: 'factor' }],
    ['legacy_px', '20px', { value: 20, unit: 'px' }],
    ['legacy_percent', '65%', { value: 0.65, unit: 'factor' }],
    ['legacy_bare', '1.5', { value: 1.5, unit: 'factor' }],
    ['auto', null, null]
  ];
  const observations = [];
  for (const [name, input, expected] of cases) {
    const config = { form: contract.createDefaultForm(), adv: contract.createDefaultAdv('circular') };
    const slot = config.adv.circular_track_slots.find(row => row.renderer === 'dinucleotide_content');
    slot.width = structuredClone(input);
    const before = structuredClone(config);
    contract.validateCurrentWriterActiveConfig({ mode: 'circular', storedConfig: config });
    contract.validateImportedCircularTrackSlots(config);
    const normalized = slotsModule.normalizeCircularTrackSlot(slot);
    assert.deepEqual(normalized.width, input);
    const payload = slotsModule.buildCircularTrackSlotPayload(normalized);
    assert.deepEqual(payload.width, expected);
    assert.deepEqual(config, before);
    observations.push({ name, input, acceptedByActiveConfigAndSlotValidation: true,
      draftRetained: true, canonical: payload.width });
  }
  for (const input of [{ value: 'bad', unit: 'px' }, { value: '1e', unit: 'factor' },
    { value: 0, unit: 'px' }, { value: -1, unit: 'factor' }, { value: 1, unit: 'em' }]) {
    const config = { form: contract.createDefaultForm(), adv: contract.createDefaultAdv('circular') };
    const slot = config.adv.circular_track_slots.find(row => row.renderer === 'dinucleotide_content');
    slot.width = structuredClone(input);
    const normalized = slotsModule.normalizeCircularTrackSlot(slot);
    assert.deepEqual(normalized.width, input);
    assert.throws(() => contract.validateImportedCircularTrackSlots(config));
    assert.throws(() => slotsModule.buildCircularTrackSlotPayload(normalized));
    observations.push({ input, rejectedBySlotValidationAndPayload: true, draftRetained: true });
  }
  const sourcePaths = ['gbdraw/web/js/app/circular-track-slots.js',
    'gbdraw/web/js/app/track-slot-validation.js', 'gbdraw/web/js/services/session-active-config-contract.js'];
  const hashes = {};
  for (const path of sourcePaths) hashes[path] = createHash('sha256').update(await readFile(join(repoRoot, path))).digest('hex');
  process.stdout.write(JSON.stringify({ sourceSha: execFileSync('git', ['rev-parse', 'HEAD'], { cwd: repoRoot, encoding: 'utf8' }).trim(),
    sourceSha256: hashes, observations,
    limits: ['DOM-free normalization, active-config validation, slot validation and payload only.',
      'No actual Save/Load, native Session reader, History or new UI was executed.'] }, null, 2) + '\n');
} finally {
  await rm(tempRoot, { recursive: true, force: true });
}
