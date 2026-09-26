// Read-only probes. Import a disposable ESM copy, as the repository's Node tests do.
import { cp, mkdtemp, writeFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const root = await mkdtemp(join(tmpdir(), 'gbdraw-issue-600-'));
try {
  await cp(join(process.cwd(), 'gbdraw/web/js'), join(root, 'js'), { recursive: true });
  await writeFile(join(root, 'package.json'), '{"type":"module"}');
  const load = (path) => import(pathToFileURL(join(root, 'js', path)));
  const { parseAnnotationTable } = await load('app/annotations/table-codec.js');
  const { prepareSpecificColorImport } = await load('app/specific-color-rules.js');
  const { validateCustomTrackPlan, customTrackPlanIssues } = await load('app/track-slot-validation.js');
  const { buildCircularTrackSlotPayload } = await load('app/circular-track-slots.js');
  const slot = (gap) => ({ id: 'gc', renderer: 'dinucleotide_content', enabled: true, side: 'inside', inner_gap_px: gap, params: { nt: 'GC' } });
  const observations = [];
  const observe = (name, operation) => {
    try { observations.push({ case: name, value: operation() }); }
    catch (error) { observations.push({ case: name, error: error.message }); }
  };
  observe('BUG-08 auxiliary column', () => parseAnnotationTable('set_id\tid\tmark\tstart\tend\tnotes\ns\ta\tband\t1\t5\tauxiliary\n').length);
  observe('BUG-10 same caption two colors', () => prepareSpecificColorImport('CDS\tgene\tabc\t#112233\tTransporter\nCDS\tgene\tmfs\t#445566\tTransporter\n').intents);
  for (const gap of ['10', '10px', '10PX', 'oops', 'px']) {
    observe(`BUG-14 raw validation ${gap}`, () => customTrackPlanIssues(validateCustomTrackPlan({ mode: 'circular', slots: [slot(gap)], axisIndex: 0, trackType: 'tuckin' })).map((issue) => issue.message));
    observe(`BUG-14 payload ${gap}`, () => buildCircularTrackSlotPayload(slot(gap)).innerGapPx);
  }
  observe('BUG-14 linear control 10px', () => customTrackPlanIssues(validateCustomTrackPlan({ mode: 'linear', slots: [{ id: 'gc', renderer: 'dinucleotide_content', enabled: true, side: 'below', height: '10px', params: { nt: 'GC' } }], axisIndex: 0 })).map((issue) => issue.message));
  console.log(JSON.stringify(observations, null, 2));
} finally {
  await rm(root, { recursive: true, force: true });
}
