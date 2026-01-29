import { runLosatPair } from '../services/losat.js';

const downloadTextFile = (filename, text) => {
  const safeName = filename || 'losat.tsv';
  const blob = new Blob([text], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = safeName;
  link.click();
  URL.revokeObjectURL(url);
};

const hashText = async (text) => {
  if (globalThis.crypto?.subtle) {
    const buffer = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text));
    return Array.from(new Uint8Array(buffer))
      .map((b) => b.toString(16).padStart(2, '0'))
      .join('');
  }
  let hash = 2166136261;
  for (let i = 0; i < text.length; i++) {
    hash ^= text.charCodeAt(i);
    hash = Math.imul(hash, 16777619);
  }
  return `fnv1a-${(hash >>> 0).toString(16)}`;
};

const makeSafeFilename = (name) => {
  const cleaned = String(name || '').replace(/[^\w.-]+/g, '_').replace(/^_+|_+$/g, '');
  return cleaned || 'losat';
};

export const createRunAnalysis = ({ state, getPyodide, writeFileToFs, refreshFeatureOverrides }) => {
  const {
    pyodideReady,
    processing,
    results,
    selectedResultIndex,
    errorLog,
    zoom,
    skipCaptureBaseConfig,
    skipPositionReapply,
    pairwiseMatchFactors,
    addedLegendCaptions,
    fileLegendCaptions,
    featureColorOverrides,
    legendEntries,
    deletedLegendEntries,
    legendColorOverrides,
    originalLegendOrder,
    originalLegendColors,
    currentColors,
    filterMode,
    manualSpecificRules,
    manualWhitelist,
    manualBlacklist,
    manualPriorityRules,
    form,
    adv,
    mode,
    cInputType,
    lInputType,
    blastSource,
    losatProgram,
    losat,
    losatCacheInfo,
    losatCache,
    files,
    linearSeqs,
    generatedLegendPosition,
    extractedFeatures,
    featureRecordIds,
    selectedFeatureRecordIdx
  } = state;

  const runAnalysis = async () => {
    if (!pyodideReady.value) return;
    const pyodide = getPyodide();
    if (!pyodide) return;
    processing.value = true;
    results.value = [];
    selectedResultIndex.value = 0;
    errorLog.value = null;
    zoom.value = 1.0;
    skipCaptureBaseConfig.value = false;
    skipPositionReapply.value = false;
    pairwiseMatchFactors.value = {};
    addedLegendCaptions.value = new Set();
    fileLegendCaptions.value = new Set();
    Object.keys(featureColorOverrides).forEach((k) => delete featureColorOverrides[k]);
    legendEntries.value = [];
    deletedLegendEntries.value = [];
    Object.keys(legendColorOverrides).forEach((k) => delete legendColorOverrides[k]);
    originalLegendOrder.value = [];
    originalLegendColors.value = {};
    window._origPairwiseMin = currentColors.value.pairwise_match_min || '#FFE7E7';
    window._origPairwiseMax = currentColors.value.pairwise_match_max || '#FF7272';

    try {
      let args = [];
      let regionSpecs = [];

      if (form.prefix && form.prefix.trim() !== '') args.push('-o', form.prefix.trim());
      if (form.species) args.push('--species', form.species);
      if (form.strain) args.push('--strain', form.strain);
      if (form.separate_strands) args.push('--separate_strands');

      if (adv.features.length) args.push('-k', adv.features.join(','));
      if (adv.window_size) args.push('--window', adv.window_size);
      if (adv.step_size) args.push('--step', adv.step_size);
      if (adv.nt && adv.nt !== 'GC') args.push('--nt', adv.nt);

      if (adv.def_font_size) args.push('--definition_font_size', adv.def_font_size);
      if (adv.label_font_size) args.push('--label_font_size', adv.label_font_size);

      if (adv.block_stroke_width !== null) args.push('--block_stroke_width', adv.block_stroke_width);
      if (adv.block_stroke_color) args.push('--block_stroke_color', adv.block_stroke_color);
      if (adv.line_stroke_width !== null) args.push('--line_stroke_width', adv.line_stroke_width);
      if (adv.line_stroke_color) args.push('--line_stroke_color', adv.line_stroke_color);
      if (adv.axis_stroke_width !== null) args.push('--axis_stroke_width', adv.axis_stroke_width);
      if (adv.axis_stroke_color) args.push('--axis_stroke_color', adv.axis_stroke_color);

      if (adv.legend_box_size) args.push('--legend_box_size', adv.legend_box_size);
      if (adv.legend_font_size) args.push('--legend_font_size', adv.legend_font_size);

      let dContent = '';
      for (const [k, v] of Object.entries(currentColors.value)) dContent += `${k}\t${v}\n`;
      pyodide.FS.writeFile('/combined_d.tsv', dContent);
      args.push('-d', '/combined_d.tsv');

      let tContent = '';
      manualSpecificRules.forEach((r) => {
        tContent += `${r.feat}\t${r.qual}\t${r.val}\t${r.color}\t${r.cap}\n`;
      });
      if (tContent.trim() !== '') {
        pyodide.FS.writeFile('/combined_t.tsv', tContent);
        args.push('-t', '/combined_t.tsv');
      }

      if (filterMode.value === 'Blacklist') {
        if (manualBlacklist.value) {
          args.push('--label_blacklist', manualBlacklist.value.replace(/\n/g, ','));
        }
      } else if (filterMode.value === 'Whitelist') {
        if (manualWhitelist.length > 0) {
          let wlContent = '';
          manualWhitelist.forEach((r) => {
            if (r.feat && r.qual) wlContent += `${r.feat}\t${r.qual}\t${r.key}\n`;
          });
          pyodide.FS.writeFile('/manual_wl.tsv', wlContent);
          args.push('--label_whitelist', '/manual_wl.tsv');
        }
      }

      let pContent = '';
      manualPriorityRules.forEach((r) => {
        pContent += `${r.feat}\t${r.order}\n`;
      });
      if (pContent.trim() !== '') {
        pyodide.FS.writeFile('/priority.tsv', pContent);
        args.push('--qualifier_priority', '/priority.tsv');
      }

      if (mode.value === 'circular') {
        args.push('--track_type', form.track_type, '-l', form.legend);
        if (form.show_labels) args.push('--show_labels');
        if (form.allow_inner_labels) args.push('--allow_inner_labels');
        if (form.suppress_gc) args.push('--suppress_gc');
        if (form.suppress_skew) args.push('--suppress_skew');

        if (adv.outer_label_x_offset) args.push('--outer_label_x_radius_offset', adv.outer_label_x_offset);
        if (adv.outer_label_y_offset) args.push('--outer_label_y_radius_offset', adv.outer_label_y_offset);
        if (adv.inner_label_x_offset) args.push('--inner_label_x_radius_offset', adv.inner_label_x_offset);
        if (adv.inner_label_y_offset) args.push('--inner_label_y_radius_offset', adv.inner_label_y_offset);
        if (adv.scale_interval) args.push('--scale_interval', adv.scale_interval);

        if (cInputType.value === 'gb') {
          if (!files.c_gb) throw new Error('Please upload a GenBank file.');
          await writeFileToFs(files.c_gb, '/input.gb');
          args.push('--gbk', '/input.gb');
        } else {
          if (!files.c_gff || !files.c_fasta) throw new Error('GFF3 and FASTA are required.');
          await writeFileToFs(files.c_gff, '/input.gff');
          await writeFileToFs(files.c_fasta, '/input.fasta');
          args.push('--gff', '/input.gff', '--fasta', '/input.fasta');
        }
      } else {
        args.push('--scale_style', form.scale_style);
        if (form.align_center) args.push('--align_center');
        if (form.show_gc) args.push('--show_gc');
        if (form.show_skew) args.push('--show_skew');
        if (form.normalize_length) args.push('--normalize_length');
        if (form.legend !== 'right') args.push('-l', form.legend);
        args.push('--bitscore', adv.min_bitscore, '--evalue', adv.evalue, '--identity', adv.identity);

        if (form.show_labels_linear !== 'none') {
          args.push('--show_labels');
          if (form.show_labels_linear === 'first') args.push('first');
        }

        if (adv.resolve_overlaps) args.push('--resolve_overlaps');
        if (adv.feature_height) args.push('--feature_height', adv.feature_height);
        if (adv.gc_height) args.push('--gc_height', adv.gc_height);
        if (adv.comparison_height) args.push('--comparison_height', adv.comparison_height);

        if (adv.scale_interval) args.push('--scale_interval', adv.scale_interval);
        if (adv.scale_font_size) args.push('--scale_font_size', adv.scale_font_size);
        if (adv.scale_stroke_width) args.push('--scale_stroke_width', adv.scale_stroke_width);
        if (adv.scale_stroke_color) args.push('--scale_stroke_color', adv.scale_stroke_color);

        const recordLabels = linearSeqs.map((seq) => (seq.definition ?? '').toString());
        const hasRecordLabels = recordLabels.some((label) => label.trim() !== '');
        if (hasRecordLabels) {
          recordLabels.forEach((label) => {
            args.push('--record_label', label);
          });
        }

        const buildRegionSpec = (seq, idx) => {
          const hasStart = seq.region_start !== null && seq.region_start !== undefined && seq.region_start !== '';
          const hasEnd = seq.region_end !== null && seq.region_end !== undefined && seq.region_end !== '';
          const recordIdRaw = seq.region_record_id ? String(seq.region_record_id).trim() : '';
          const hasAny =
            hasStart ||
            hasEnd ||
            recordIdRaw !== '' ||
            Boolean(seq.region_reverse);
          if (!hasAny) return null;
          if (!hasStart || !hasEnd) {
            throw new Error(`Sequence #${idx + 1}: Region start and end are required.`);
          }
          const start = Number(seq.region_start);
          const end = Number(seq.region_end);
          if (!Number.isFinite(start) || !Number.isFinite(end)) {
            throw new Error(`Sequence #${idx + 1}: Region start/end must be numbers.`);
          }
          if (!Number.isInteger(start) || !Number.isInteger(end)) {
            throw new Error(`Sequence #${idx + 1}: Region start/end must be integers.`);
          }
          if (start < 1 || end < 1) {
            throw new Error(`Sequence #${idx + 1}: Region start/end must be >= 1.`);
          }
          const specBody = `${start}-${end}${seq.region_reverse ? ':rc' : ''}`;
          const cliSpec = recordIdRaw ? `${recordIdRaw}:${specBody}` : `#${idx + 1}:${specBody}`;
          const fileSpec = recordIdRaw ? `${recordIdRaw}:${specBody}` : specBody;
          return { cli: cliSpec, file: fileSpec };
        };

        regionSpecs = linearSeqs.map((seq, idx) => buildRegionSpec(seq, idx));
        regionSpecs.forEach((spec) => {
          if (spec?.cli) args.push('--region', spec.cli);
        });

        let inputArgs = [];
        let blastArgs = [];
        const useLosat = blastSource.value === 'losat';
        const fastaCache = new Map();
        const fastaHashCache = new Map();
        const textEncoder = new TextEncoder();
        let extractFirstFasta = null;
        let cacheInfo = [];
        const cacheMap = losatCache.value || new Map();

        if (useLosat) {
          extractFirstFasta = pyodide.globals.get('extract_first_fasta');
        } else {
          losatCacheInfo.value = [];
        }

        const getSeqEntry = (idx) => {
          if (fastaCache.has(idx)) return fastaCache.get(idx);
          const path = lInputType.value === 'gb' ? `/seq_${idx}.gb` : `/seq_${idx}.fasta`;
          const fmt = lInputType.value === 'gb' ? 'genbank' : 'fasta';
          const regionSpec = regionSpecs[idx]?.file || null;
          const res = JSON.parse(extractFirstFasta(path, fmt, regionSpec));
          if (res.error) throw new Error(res.error);
          const entry = {
            fasta: res.fasta,
            recordId: res.record_id || `seq_${idx + 1}`
          };
          fastaCache.set(idx, entry);
          return entry;
        };

        const getSeqHash = async (idx) => {
          if (fastaHashCache.has(idx)) return fastaHashCache.get(idx);
          const entry = getSeqEntry(idx);
          const hash = await hashText(entry.fasta);
          fastaHashCache.set(idx, hash);
          return hash;
        };

        const buildCacheKey = async (argsKey, queryIdx, subjectIdx) => {
          const queryHash = await getSeqHash(queryIdx);
          const subjectHash = await getSeqHash(subjectIdx);
          const payload = JSON.stringify({
            program: losatProgram.value,
            outfmt: String(losat.outfmt || '6'),
            args: argsKey,
            queryHash,
            subjectHash
          });
          return hashText(payload);
        };

        const buildCacheFilename = (pairIndex, queryEntry, subjectEntry) => {
          const left = makeSafeFilename(queryEntry.recordId || `seq_${pairIndex + 1}`);
          const right = makeSafeFilename(subjectEntry.recordId || `seq_${pairIndex + 2}`);
          return `${left}_vs_${right}.losat.tsv`;
        };

        const pushArg = (arr, flag, value) => {
          if (value === null || value === undefined || value === '') return;
          if (typeof value === 'number' && !Number.isFinite(value)) return;
          const valueStr = String(value);
          if (valueStr.startsWith('-')) {
            arr.push(`${flag}=${valueStr}`);
          } else {
            arr.push(flag, valueStr);
          }
        };

        const getGencode = (idx) => {
          const raw = linearSeqs[idx]?.losat_gencode;
          if (raw === null || raw === undefined || raw === '') return null;
          const num = Number(raw);
          if (!Number.isFinite(num)) return null;
          return num;
        };

        const buildLosatArgs = (queryIdx, subjectIdx) => {
          const args = [];
          if (losatProgram.value === 'blastn') {
            pushArg(args, '--task', losat.blastn.task);
          } else {
            pushArg(args, '--query-gencode', getGencode(queryIdx));
            pushArg(args, '--db-gencode', getGencode(subjectIdx));
          }
          return args;
        };

        for (let i = 0; i < linearSeqs.length; i++) {
          const seq = linearSeqs[i];
          if (lInputType.value === 'gb') {
            if (!seq.gb) throw new Error(`Sequence #${i + 1}: Missing GenBank file.`);
            await writeFileToFs(seq.gb, `/seq_${i}.gb`);
            inputArgs.push(`/seq_${i}.gb`);
          } else {
            if (!seq.gff || !seq.fasta) throw new Error(`Sequence #${i + 1}: GFF3 and FASTA are required.`);
            await writeFileToFs(seq.gff, `/seq_${i}.gff`);
            await writeFileToFs(seq.fasta, `/seq_${i}.fasta`);
          }
        }

        for (let i = 0; i < linearSeqs.length - 1; i++) {
          const seq = linearSeqs[i];
          if (useLosat) {
            const queryEntry = getSeqEntry(i);
            const subjectEntry = getSeqEntry(i + 1);
            const losatArgs = buildLosatArgs(i, i + 1);
            const cacheKey = await buildCacheKey(losatArgs, i, i + 1);
            const cached = cacheMap.get(cacheKey);
            const blastText = cached
              ? cached.text
              : await runLosatPair({
                  program: losatProgram.value,
                  queryFasta: queryEntry.fasta,
                  subjectFasta: subjectEntry.fasta,
                  outfmt: losat.outfmt || '6',
                  extraArgs: losatArgs
                });
            if (!cached) {
              cacheMap.set(cacheKey, { text: blastText });
            }
            cacheInfo.push({
              key: cacheKey,
              filename: buildCacheFilename(i, queryEntry, subjectEntry)
            });
            pyodide.FS.writeFile(`/blast_${i}.txt`, textEncoder.encode(blastText));
            blastArgs.push(`/blast_${i}.txt`);
          } else if (seq.blast) {
            await writeFileToFs(seq.blast, `/blast_${i}.txt`);
            blastArgs.push(`/blast_${i}.txt`);
          }
        }
        if (extractFirstFasta) {
          extractFirstFasta.destroy();
        }
        if (useLosat) {
          losatCacheInfo.value = cacheInfo;
          losatCache.value = cacheMap;
        }
        if (lInputType.value === 'gb') args.push('--gbk', ...inputArgs);
        else {
          let gffs = [];
          let fastas = [];
          for (let i = 0; i < linearSeqs.length; i++) {
            gffs.push(`/seq_${i}.gff`);
            fastas.push(`/seq_${i}.fasta`);
          }
          args.push('--gff', ...gffs, '--fasta', ...fastas);
        }
        if (blastArgs.length) args.push('-b', ...blastArgs);
      }

      console.log('CMD:', args.join(' '));
      const jsonResult = pyodide
        .globals
        .get('run_gbdraw_wrapper')(mode.value, pyodide.toPy(args.map(String)));
      const res = JSON.parse(jsonResult);
      if (res.error) throw new Error(res.error);
      results.value = res;

      generatedLegendPosition.value = form.legend;

      extractedFeatures.value = [];

      if (mode.value === 'circular' && cInputType.value === 'gb') {
        try {
          const featJson = pyodide.globals.get('extract_features_from_genbank')('/input.gb');
          const featData = JSON.parse(featJson);
          if (!featData.error && featData.features) {
            extractedFeatures.value = featData.features;
            featureRecordIds.value = featData.record_ids || [];
            selectedFeatureRecordIdx.value = 0;
            refreshFeatureOverrides(featData.features);
            console.log(
              `Extracted ${featData.features.length} features from ${featData.record_ids.length} record(s) for color editor.`
            );
          }
        } catch (e) {
          console.log('Could not extract features:', e);
        }
      } else if (mode.value === 'linear' && lInputType.value === 'gb' && linearSeqs.length > 0) {
        try {
          let allFeatures = [];
          let allRecordLabels = [];
          for (let i = 0; i < linearSeqs.length; i++) {
            const regionSpec = regionSpecs[i]?.file || null;
            const featJson = pyodide.globals.get('extract_features_from_genbank')(`/seq_${i}.gb`, regionSpec);
            const featData = JSON.parse(featJson);
            if (!featData.error && featData.features) {
              featData.features.forEach((f) => {
                f.fileIdx = i;
                f.displayRecordId = `File ${i + 1}: ${f.record_id}`;
                f.id = `file${i}_${f.id}`;
              });
              allFeatures = allFeatures.concat(featData.features);
              featData.record_ids.forEach((rid, ridx) => {
                allRecordLabels.push({ label: `File ${i + 1}: ${rid}`, fileIdx: i, recordIdx: ridx });
              });
            }
          }
          extractedFeatures.value = allFeatures;
          featureRecordIds.value = allRecordLabels.map((r) => r.label);
          selectedFeatureRecordIdx.value = 0;
          refreshFeatureOverrides(allFeatures);
          console.log(
            `Extracted ${allFeatures.length} features from ${linearSeqs.length} file(s) for color editor.`
          );
        } catch (e) {
          console.log('Could not extract features:', e);
        }
      }
    } catch (e) {
      errorLog.value = e.message;
    } finally {
      processing.value = false;
    }
  };

  const downloadLosatCache = async () => {
    if (!losatCacheInfo.value || losatCacheInfo.value.length === 0) return;
    const cacheMap = losatCache.value;
    if (!cacheMap || cacheMap.size === 0) return;

    const totalChars = losatCacheInfo.value.reduce((sum, entry) => {
      const cached = cacheMap.get(entry.key);
      return sum + (cached?.text ? cached.text.length : 0);
    }, 0);

    if (totalChars > 50 * 1024 * 1024) {
      const proceed = confirm(
        `LOSAT TSV export will download about ${(totalChars / (1024 * 1024)).toFixed(1)} MB. Continue?`
      );
      if (!proceed) return;
    }

    for (let idx = 0; idx < losatCacheInfo.value.length; idx += 1) {
      const entry = losatCacheInfo.value[idx];
      const cached = cacheMap.get(entry.key);
      if (!cached) continue;
      const filename = entry.filename || `losat_pair_${idx + 1}.tsv`;
      downloadTextFile(filename, cached.text);
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
  };

  const clearLosatCache = () => {
    if (losatCache.value) {
      losatCache.value.clear();
    }
    losatCacheInfo.value = [];
  };

  return { runAnalysis, downloadLosatCache, clearLosatCache };
};
