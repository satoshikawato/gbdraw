import assert from 'node:assert/strict';
import { cp, mkdtemp, readFile, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { createHash } from 'node:crypto';
import { execFileSync } from 'node:child_process';

const root = resolve(process.cwd());
const cases = JSON.parse(await readFile(new URL('./scalar-fixtures.json', import.meta.url)));
const decode = value => {
  if (value && typeof value === 'object' && '$number' in value) return Number(value.$number);
  if (value && typeof value === 'object') return Object.fromEntries(Object.entries(value).map(([k,v]) => [k,decode(v)]));
  return value;
};
const capture = fn => {
  try { return { accepted: true, output: fn() ?? null }; }
  catch (e) { return { accepted: false, error: e.message }; }
};
const temp = await mkdtemp(join(tmpdir(), 'gbdraw-619-s00-node-'));
try {
  await cp(join(root, 'gbdraw/web/js'), join(temp, 'js'), { recursive: true });
  await writeFile(join(temp, 'package.json'), '{"type":"module"}');
  const slots = await import(pathToFileURL(join(temp, 'js/app/circular-track-slots.js')));
  const contract = await import(pathToFileURL(join(temp, 'js/services/session-active-config-contract.js')));
  const observations = [];
  for (const sample of cases) for (const field of ['width', 'radius']) {
    const config = { form: contract.createDefaultForm(), adv: contract.createDefaultAdv('circular') };
    const row = config.adv.circular_track_slots.find(s => s.renderer === 'dinucleotide_content');
    row[field] = decode(sample.input);
    const before = structuredClone(config);
    const normalized = slots.normalizeCircularTrackSlot(row);
    assert.deepEqual(normalized[field], row[field]);
    const validation = capture(() => contract.validateImportedCircularTrackSlots(config));
    const activeConfig = capture(() => contract.validateCurrentWriterActiveConfig({mode:'circular',storedConfig:config}));
    const payload = capture(() => slots.buildCircularTrackSlotPayload(normalized)[field]);
    assert.deepEqual(config, before);
    if (sample.valid) {
      assert.equal(validation.accepted, true, sample.name);
      assert.equal(activeConfig.accepted, true, sample.name);
      assert.deepEqual(payload.output, sample.canonical, sample.name);
    }
    if (payload.accepted && payload.output !== null) {
      assert.equal(typeof payload.output.value, 'number');
      assert.ok(Number.isFinite(payload.output.value) && payload.output.value > 0);
      assert.ok(['px','factor'].includes(payload.output.unit));
    }
    observations.push({name:sample.name,field,input:sample.input,normalizationRetained:true,validation,activeConfig,payload,sourceUnchanged:true});
  }
  const sources = ['gbdraw/web/js/app/circular-track-slots.js','gbdraw/web/js/app/track-slot-validation.js','gbdraw/web/js/services/session-active-config-contract.js'];
  const sourceSha256 = Object.fromEntries(await Promise.all(sources.map(async path => [path,createHash('sha256').update(await readFile(join(root,path))).digest('hex')])));
  process.stdout.write(JSON.stringify({sourceSha:execFileSync('git',['rev-parse','HEAD'],{encoding:'utf8'}).trim(),nodeVersion:process.version,sourceSha256,observations,limits:['DOM-free only. Actual browser Save/Load, native reader and History are recorded separately.','Observed coercions are not Product authority or permission to expand admission.']},null,2)+'\n');
} finally { await rm(temp,{recursive:true,force:true}); }
