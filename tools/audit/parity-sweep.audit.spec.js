// GUI Generate vs. Source recipe sweep: a baseline plus 31 Circular and 25 Linear option probes.
// Each probe resets the form, applies one option family, Generates, and saves the Result SVG,
// the Source recipe command, and its helper files. Replay the commands with parity_replay.py.
// Filter probes with AUDIT_ONLY=name1,name2; choose the Circular input with AUDIT_INPUT.
const { test } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const A = require('./helpers/audit-common.cjs');

const ONLY = process.env.AUDIT_ONLY ? process.env.AUDIT_ONLY.split(',') : null;
const select = (probes) => (ONLY ? probes.filter(([name]) => ONLY.includes(name)) : probes);

const CIRCULAR_PROBES = [
  ['labels_both_font20', 'a.form.labels_mode="both"; a.adv.label_font_size=20;'],
  ['track_middle', 'a.form.track_type="middle";'],
  ['track_spreadout_labels_out', 'a.form.track_type="spreadout"; a.form.labels_mode="out";'],
  ['legend_upper_right', 'a.form.legend="upper_right";'],
  ['legend_none', 'a.form.legend="none";'],
  ['title_top', 'a.adv.plot_title_position="top"; a.form.plot_title="My Title"; a.adv.plot_title_font_size=40;'],
  ['species_strain', 'a.form.species="Homo sapiens"; a.form.strain="Ref";'],
  ['species_title_bottom', 'a.form.species="Homo sapiens"; a.form.strain="Ref"; a.adv.plot_title_position="bottom"; a.adv.keep_full_definition_with_plot_title=true;'],
  ['nt_AT_window', 'a.adv.nt="AT"; a.adv.window_size=200; a.adv.step_size=50;'],
  ['scale_interval', 'a.adv.scale_interval=2000;'],
  ['scale_interval_neg', 'a.adv.scale_interval=-1000;'],
  ['strokes', 'a.adv.block_stroke_width=1.5; a.adv.block_stroke_color="#ff0000"; a.adv.axis_stroke_width=4; a.adv.axis_stroke_color="#0000ff";'],
  ['def_font', 'a.adv.def_font_size=30;'],
  ['legend_sizes', 'a.adv.legend_font_size=10; a.adv.legend_box_size=12;'],
  ['gc_mode_absolute', 'a.adv.gc_content_mode="absolute"; a.adv.gc_content_min_percent=20; a.adv.gc_content_max_percent=60;'],
  ['no_gc_no_skew', 'a.form.suppress_gc=true; a.form.suppress_skew=true;'],
  ['single_strand', 'a.form.separate_strands=false;'],
  ['reverse', 'a.form.circular_reverse=true;'],
  ['region', 'a.form.circular_region_start=1000; a.form.circular_region_end=9000;'],
  ['no_multi_canvas', 'a.form.multi_record_canvas=false;'],
  ['arrow_ratio', 'a.adv.arrow_head_length_ratio=0.5; a.adv.arrow_shaft_width_ratio=0.5;'],
  ['label_rendering_curved', 'a.form.labels_mode="both"; a.adv.label_rendering="curved";'],
  ['label_spacing_offsets', 'a.form.labels_mode="both"; a.adv.circular_label_spacing=8; a.adv.outer_label_x_offset=1.1; a.adv.inner_label_y_offset=0.9;'],
  ['tick_font', 'a.adv.tick_label_font_size=20;'],
  ['hide_scale', 'a.form.show_scale=false;'],
  ['palette_forest', 'a.selectedPalette="forest";'],
  ['center_radius', 'a.adv.center_reserved_radius=100;'],
  ['feature_width', 'a.adv.feature_width_circular=40;'],
  ['gc_width_radius', 'a.adv.gc_content_width_circular=30; a.adv.gc_content_radius_circular=0.5;'],
  ['labels_placement_radial', 'a.form.labels_mode="out"; a.adv.circular_label_placement="radial";'],
  ['resolve_overlaps', 'a.adv.resolve_overlaps=true; a.adv.feature_overlap_tolerance_bp=10;']
];

const LINEAR_PROBES = [
  ['labels_all_font', 'a.form.show_labels_linear="all"; a.adv.label_font_size=12;'],
  ['layout_above', 'a.form.linear_track_layout="above";'],
  ['layout_below_axisgap', 'a.form.linear_track_layout="below"; a.adv.track_axis_gap=10;'],
  ['feature_height', 'a.adv.feature_height=40;'],
  ['legend_top', 'a.form.legend="top";'],
  ['legend_none', 'a.form.legend="none";'],
  ['plot_title', 'a.form.plot_title="Linear Title"; a.adv.plot_title_position="top"; a.adv.plot_title_font_size=30;'],
  ['def_font', 'a.adv.def_font_size=14;'],
  ['gc_skew_on', 'a.form.show_gc=true; a.form.show_skew=true; a.adv.gc_height=30;'],
  ['scale_ruler', 'a.form.scale_style="ruler"; a.adv.scale_interval=10000; a.adv.scale_font_size=10;'],
  ['ruler_on_axis', 'a.form.scale_style="ruler"; a.form.linear_ruler_on_axis=true;'],
  ['normalize', 'a.form.normalize_length=true;'],
  ['align_center', 'a.form.align_center=true;'],
  ['separate_strands_off', 'a.form.separate_strands=false;'],
  ['blast_filters', 'a.adv.evalue="1e-10"; a.adv.min_bitscore=100; a.adv.identity=50; a.adv.alignment_length=100;'],
  ['comparison_height', 'a.adv.comparison_height=100; a.adv.pairwise_match_style="ribbon";'],
  ['strokes', 'a.adv.block_stroke_width=1; a.adv.line_stroke_width=2; a.adv.axis_stroke_width=3;'],
  ['label_rotation', 'a.form.show_labels_linear="all"; a.adv.label_rotation=45;'],
  ['record_gap', 'a.linearRecordGap=50;'],
  ['show_replicon', 'a.adv.linear_show_replicon=true;'],
  ['hide_scale', 'a.form.show_scale=false;'],
  ['keep_def_left_off', 'a.form.keep_definition_left_aligned=false;'],
  ['legend_sizes', 'a.adv.legend_font_size=10; a.adv.legend_box_size=10;'],
  ['label_placement_above', 'a.form.show_labels_linear="all"; a.adv.label_placement="above_feature";'],
  ['region_rc', 'a.linearSeqs[0].region_start=1000; a.linearSeqs[0].region_end=60000; a.linearSeqs[0].region_reverse=true;']
];

// Every probe starts from a freshly opened app with the inputs loaded, so no option leaks
// from one probe into the next.
const sweep = async (page, outdir, probes, setup) => {
  const log = A.collect(page);
  await setup();
  const base = await A.runGenerate(page);
  A.writeEvidence(outdir, 'baseline.svg', base.svg);
  const baseHelpers = base.status === 'ok' ? await A.captureHelpers(page, join(outdir, 'helpers_baseline')) : [];
  const manifest = [{
    name: 'baseline', status: base.status, command: base.command, exact: base.exactReplay,
    err: base.errorSummary, details: base.errorDetails, helpers: baseHelpers
  }];
  for (const [name, body] of probes) {
    await setup();
    await page.evaluate(new Function(`const a = window.__GBDRAW_APP__; ${body}`));
    await page.waitForTimeout(300);
    const out = await A.runGenerate(page);
    A.writeEvidence(outdir, `${name}.svg`, out.svg);
    const helperFiles = out.status === 'ok' ? await A.captureHelpers(page, join(outdir, `helpers_${name}`)) : [];
    manifest.push({
      name, status: out.status, command: out.command, exact: out.exactReplay, recipe: out.sourceRecipe,
      err: out.errorSummary, details: out.errorDetails, sameAsBaseline: out.svg === base.svg, helpers: helperFiles
    });
  }
  manifest.push({ log });
  A.writeEvidence(outdir, 'manifest.json', manifest);
};

test('Circular GUI Generate vs. Source recipe sweep', async ({ page }) => {
  test.setTimeout(1_800_000);
  const input = process.env.AUDIT_INPUT || 'HmmtDNA.gbk';
  const outdir = A.outDir('parity', `circular-${input.replace(/\W+/g, '_')}`);
  page.on('dialog', (d) => d.dismiss());
  await sweep(page, outdir, select(CIRCULAR_PROBES), async () => {
    await A.helpers().openApp(page);
    await A.setFile(page, 'c_gb', join(A.REPO_ROOT, 'tests/test_inputs', input));
  });
});

test('Linear GUI Generate vs. Source recipe sweep', async ({ page }) => {
  test.setTimeout(1_800_000);
  const outdir = A.outDir('parity', 'linear-MJNV');
  const examples = join(A.REPO_ROOT, 'examples');
  const inputs = {
    a1: readFileSync(join(examples, 'MjeNMV.gb'), 'utf8'),
    a2: readFileSync(join(examples, 'MelaMJNV.gb'), 'utf8'),
    b: readFileSync(join(examples, 'MjeNMV.MelaMJNV.tblastx.out'), 'utf8')
  };
  page.on('dialog', (d) => d.dismiss());
  await sweep(page, outdir, select(LINEAR_PROBES), async () => {
    await A.helpers().openApp(page);
    await page.evaluate(async ({ a1, a2, b }) => {
      const app = window.__GBDRAW_APP__;
      app.setDiagramMode('linear');
      await window.Vue.nextTick();
      if (app.linearSeqs.length < 2) app.addLinearSeq();
      app.setLinearSeqPrimaryFile(0, 'gb', new File([a1], 'MjeNMV.gb', { type: 'text/plain' }));
      app.setLinearSeqPrimaryFile(1, 'gb', new File([a2], 'MelaMJNV.gb', { type: 'text/plain' }));
      app.losatProgram = 'blastn';
      app.linearComparisonPlan.mode = 'adjacent';
      app.linearComparisonPlan.defaultSource = 'upload';
      app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
        id: 'audit-edge', queryUid: app.linearSeqs[0].uid, subjectUid: app.linearSeqs[1].uid,
        included: true, fileActive: true, losatFilenameActive: false, source: 'upload',
        file: new File([b], 'MjeNMV.MelaMJNV.tblastx.out', { type: 'text/plain' }), losatFilename: ''
      });
      await window.Vue.nextTick();
    }, inputs);
  });
});
