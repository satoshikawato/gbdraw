import { getFeatureCaption } from '../feature-utils.js';

export const createFeatureColorActions = ({
  state,
  nextTick,
  legendActions,
  svgActions,
  ruleActions,
  featureSvgActions
}) => {
  const {
    pyodideReady,
    results,
    selectedResultIndex,
    currentColors,
    manualSpecificRules,
    extractedFeatures,
    featureColorOverrides,
    svgContainer,
    clickedFeature,
    colorScopeDialog,
    resetColorDialog,
    legendEntries,
    legendStrokeOverrides,
    originalSvgStroke,
    skipCaptureBaseConfig,
    skipExtractOnSvgChange,
    addedLegendCaptions
  } = state;

  const {
    addLegendEntry,
    removeLegendEntry,
    updateLegendEntryColorByCaption,
    extractLegendEntries,
    getAllFeatureLegendGroups
  } = legendActions;
  const { applySpecificRulesToSvg } = svgActions;
  const {
    countFeaturesMatchingRule,
    findExistingColorForCaption,
    findFeaturesWithSameIndividualLabel,
    findFeaturesWithSameLegendItem,
    findMatchingRegexRule,
    getEffectiveLegendCaption,
    getIndividualFeatureLabel,
    getFeatureQualifier
  } = ruleActions;
  const { applyInstantPreview } = featureSvgActions;
  const normalizeCaption = (value) => String(value || '').trim();
  const normalizeCaptionKey = (value) => normalizeCaption(value).toLowerCase();
  const captionsMatch = (left, right) => normalizeCaptionKey(left) === normalizeCaptionKey(right);
  const SUFFIXED_CAPTION_PATTERN = /^(.*?)\s*\((\d+)\)$/;

  const findLegendEntryByCaption = (caption) => {
    const normalizedCaption = normalizeCaptionKey(caption);
    if (!normalizedCaption) return null;

    return (
      legendEntries.value.find(
        (entry) => normalizeCaptionKey(entry?.caption) === normalizedCaption
      ) || null
    );
  };

  const findExistingCaptionColor = (feat, caption) => {
    const existingCaption = findExistingColorForCaption(feat, caption);
    if (existingCaption?.color) {
      return {
        caption: existingCaption.rule?.cap || caption,
        color: existingCaption.color,
        rule: existingCaption.rule || null
      };
    }

    const legendEntry = findLegendEntryByCaption(caption);
    if (legendEntry?.color) {
      return {
        caption: legendEntry.caption,
        color: legendEntry.color,
        rule: null
      };
    }

    return null;
  };

  const clearColorScopeDialog = () => {
    colorScopeDialog.show = false;
    colorScopeDialog.feat = null;
    colorScopeDialog.color = null;
    colorScopeDialog.matchingRule = null;
    colorScopeDialog.ruleMatchCount = 0;
    colorScopeDialog.legendName = null;
    colorScopeDialog.siblingCount = 0;
    colorScopeDialog.individualLabel = null;
    colorScopeDialog.individualLabelSiblingCount = 0;
    colorScopeDialog.existingCaptionRule = null;
    colorScopeDialog.existingCaptionColor = null;
  };

  const reclaimOrphanedBaseCaptions = async () => {
    if (extractedFeatures.value.length === 0) return false;

    const catalog = new Map();
    const rememberCaption = (caption) => {
      const normalized = normalizeCaption(caption);
      if (!normalized) return;
      const key = normalizeCaptionKey(normalized);
      if (!catalog.has(key)) {
        catalog.set(key, normalized);
      }
    };

    legendEntries.value.forEach((entry) => rememberCaption(entry?.caption));
    manualSpecificRules.forEach((rule) => rememberCaption(rule?.cap));
    Object.values(featureColorOverrides).forEach((override) => rememberCaption(override?.caption));

    const usageByCaptionKey = new Map();
    extractedFeatures.value.forEach((feat) => {
      const effectiveCaption = normalizeCaption(getEffectiveLegendCaption(feat));
      if (!effectiveCaption) return;
      rememberCaption(effectiveCaption);
      const key = normalizeCaptionKey(effectiveCaption);
      usageByCaptionKey.set(key, (usageByCaptionKey.get(key) || 0) + 1);
    });

    const suffixCandidatesByBase = new Map();
    for (const [captionKey, captionRaw] of catalog.entries()) {
      const match = captionRaw.match(SUFFIXED_CAPTION_PATTERN);
      if (!match) continue;

      const baseRaw = normalizeCaption(match[1]);
      if (!baseRaw) continue;

      const baseKey = normalizeCaptionKey(baseRaw);
      rememberCaption(baseRaw);
      const parsedIndex = Number.parseInt(match[2], 10);
      const index = Number.isFinite(parsedIndex) ? parsedIndex : Number.MAX_SAFE_INTEGER;

      if (!suffixCandidatesByBase.has(baseKey)) {
        suffixCandidatesByBase.set(baseKey, []);
      }
      suffixCandidatesByBase.get(baseKey).push({
        captionKey,
        captionRaw,
        baseRaw,
        index
      });
    }

    let changed = false;

    for (const [baseKey, candidates] of suffixCandidatesByBase.entries()) {
      const baseUsage = usageByCaptionKey.get(baseKey) || 0;
      if (baseUsage > 0) continue;

      const ordered = [...candidates].sort((a, b) => {
        if (a.index !== b.index) return a.index - b.index;
        return a.captionRaw.localeCompare(b.captionRaw);
      });

      const activeCandidate =
        ordered.find((candidate) => (usageByCaptionKey.get(candidate.captionKey) || 0) > 0) || ordered[0];
      if (!activeCandidate) continue;
      if (activeCandidate.captionKey === baseKey) continue;

      const baseCaption = catalog.get(baseKey) || activeCandidate.baseRaw;
      const suffixCaption = activeCandidate.captionRaw;

      let sourceColor = null;
      const suffixLegendEntry = findLegendEntryByCaption(suffixCaption);
      if (suffixLegendEntry?.color) {
        sourceColor = suffixLegendEntry.color;
      }

      manualSpecificRules.forEach((rule) => {
        if (!captionsMatch(rule.cap, suffixCaption)) return;
        if (!sourceColor && rule.color) sourceColor = rule.color;
        rule.cap = baseCaption;
        changed = true;
      });

      Object.values(featureColorOverrides).forEach((override) => {
        if (!override || !captionsMatch(override.caption, suffixCaption)) return;
        if (!sourceColor && override.color) sourceColor = override.color;
        override.caption = baseCaption;
        changed = true;
      });

      const baseLegendEntry = findLegendEntryByCaption(baseCaption);
      if (baseLegendEntry) {
        if (sourceColor) {
          updateLegendEntryColorByCaption(baseLegendEntry.caption, sourceColor);
          changed = true;
        }
      } else {
        const colorToUse = sourceColor || '#cccccc';
        const addedCaption = await addLegendEntry(baseCaption, colorToUse);
        if (addedCaption && typeof addedCaption === 'string') {
          addedLegendCaptions.value.add(addedCaption);
          const addedKey = normalizeCaptionKey(addedCaption);
          if (addedKey !== baseKey) {
            manualSpecificRules.forEach((rule) => {
              if (captionsMatch(rule.cap, baseCaption)) {
                rule.cap = addedCaption;
              }
            });
            Object.values(featureColorOverrides).forEach((override) => {
              if (override && captionsMatch(override.caption, baseCaption)) {
                override.caption = addedCaption;
              }
            });
          }
          changed = true;
        }
      }

      removeLegendEntry(suffixLegendEntry?.caption || suffixCaption);
      changed = true;
    }

    if (changed) {
      applySpecificRulesToSvg();
      await nextTick();
      extractLegendEntries();
    }

    return changed;
  };

  const requestFeatureColorChange = async (feat, color, requestedLegendName = null, options = {}) => {
    if (!feat) return;
    const requestedCaption = normalizeCaption(requestedLegendName);
    const effectiveCaption = normalizeCaption(getEffectiveLegendCaption(feat));
    const fallbackCaption = normalizeCaption(getFeatureCaption(feat));
    const legendName = requestedCaption || effectiveCaption || fallbackCaption;
    if (!legendName) return;

    const matchingRule = findMatchingRegexRule(feat);
    const ruleMatchCount = matchingRule ? countFeaturesMatchingRule(matchingRule) : 0;

    const siblings = findFeaturesWithSameLegendItem(feat, legendName);
    const siblingCount = siblings.length;

    const existingCaption = findExistingCaptionColor(feat, legendName);
    const individualLabel = getIndividualFeatureLabel(feat);
    let individualLabelSiblingCount = 0;
    if (matchingRule && ruleMatchCount > 1) {
      const ruleCaption = matchingRule.cap || matchingRule.val;
      if (individualLabel !== ruleCaption) {
        individualLabelSiblingCount = findFeaturesWithSameIndividualLabel(feat, individualLabel).length;
      }
    }

    const needsDialog = Boolean(matchingRule) || siblingCount > 0;
    if (needsDialog) {
      colorScopeDialog.show = true;
      colorScopeDialog.feat = feat;
      colorScopeDialog.color = color;
      colorScopeDialog.matchingRule = matchingRule;
      colorScopeDialog.ruleMatchCount = ruleMatchCount;
      colorScopeDialog.legendName = legendName;
      colorScopeDialog.siblingCount = siblingCount;
      colorScopeDialog.individualLabel = individualLabel;
      colorScopeDialog.individualLabelSiblingCount = individualLabelSiblingCount;
      colorScopeDialog.existingCaptionRule = existingCaption?.rule || null;
      colorScopeDialog.existingCaptionColor = existingCaption?.color || null;
      if (options.closePopupOnDialog) {
        clickedFeature.value = null;
      }
      return;
    }

    if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id) {
      clickedFeature.value.color = color;
      if (requestedCaption) {
        clickedFeature.value.legendName = requestedCaption;
      }
    }
    await setFeatureColor(feat, color, legendName);
  };

  const updateClickedFeatureColor = async (color) => {
    if (!clickedFeature.value) return;
    const feat = clickedFeature.value.feat;
    if (!feat) return;
    const customName = normalizeCaption(clickedFeature.value.legendName);
    await requestFeatureColorChange(feat, color, customName, { closePopupOnDialog: true });
  };

  const handleLegendNameCommit = async () => {
    if (!clickedFeature.value) return;
    const feat = clickedFeature.value.feat;
    if (!feat) return;

    const requestedCaption = clickedFeature.value.legendName?.trim();
    if (!requestedCaption) return;

    const existingCaption = findExistingCaptionColor(feat, requestedCaption);
    if (!existingCaption?.color) return;

    const currentColor = String(clickedFeature.value.color || '').toLowerCase();
    const existingColor = String(existingCaption.color).toLowerCase();
    if (!currentColor || currentColor === existingColor) return;

    const matchedCaption = existingCaption.caption || requestedCaption;
    const siblings = findFeaturesWithSameLegendItem(feat, matchedCaption);
    colorScopeDialog.show = true;
    colorScopeDialog.feat = feat;
    colorScopeDialog.color = clickedFeature.value.color;
    colorScopeDialog.matchingRule = null;
    colorScopeDialog.ruleMatchCount = 0;
    colorScopeDialog.legendName = matchedCaption;
    colorScopeDialog.siblingCount = siblings.length;
    colorScopeDialog.individualLabel = null;
    colorScopeDialog.individualLabelSiblingCount = 0;
    colorScopeDialog.existingCaptionRule = existingCaption.rule || null;
    colorScopeDialog.existingCaptionColor = existingCaption.color;
  };

  const selectLegendNameOption = async (caption) => {
    if (!clickedFeature.value) return;
    const selectedCaption = String(caption || '').trim();
    if (!selectedCaption) return;
    clickedFeature.value.legendName = selectedCaption;
    await handleLegendNameCommit();
  };

  const handleColorScopeChoice = async (choice) => {
    const { feat, color, matchingRule, legendName, existingCaptionColor } = colorScopeDialog;
    if (choice === 'cancel' || !feat || !color) {
      clearColorScopeDialog();
      return;
    }

    if (choice === 'rule') {
      if (matchingRule) {
        matchingRule.color = color;
        if (matchingRule.cap) {
          updateLegendEntryColorByCaption(matchingRule.cap, color);
          extractLegendEntries();
        }
      }
      applySpecificRulesToSvg();
    } else if (choice === 'caption') {
      const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
      if (!targetLegendName) {
        clearColorScopeDialog();
        return;
      }
      const siblings = findFeaturesWithSameLegendItem(feat, targetLegendName);
      const allFeatures = [feat, ...siblings];
      const existingLegendEntry = findLegendEntryByCaption(targetLegendName);
      let finalCaption = existingLegendEntry?.caption || targetLegendName;

      if (existingLegendEntry) {
        updateLegendEntryColorByCaption(existingLegendEntry.caption, color);
      } else {
        const addedCaption = await addLegendEntry(targetLegendName, color);
        if (addedCaption && typeof addedCaption === 'string') {
          finalCaption = addedCaption;
          addedLegendCaptions.value.add(addedCaption);
        }
      }

      for (const f of allFeatures) {
        const existingIdx = manualSpecificRules.findIndex(
          (r) => String(r.qual || '').toLowerCase() === 'hash' && r.val === f.svg_id
        );
        if (existingIdx >= 0) {
          manualSpecificRules[existingIdx].color = color;
          manualSpecificRules[existingIdx].cap = finalCaption;
        } else {
          manualSpecificRules.push({
            feat: f.type,
            qual: 'hash',
            val: f.svg_id,
            color: color,
            cap: finalCaption
          });
        }
      }
      applySpecificRulesToSvg();
      await reclaimOrphanedBaseCaptions();
      extractLegendEntries();
    } else if (choice === 'single') {
      let singleCaption = legendName;
      if (matchingRule && colorScopeDialog.ruleMatchCount > 1) {
        const ruleCaption = matchingRule.cap || matchingRule.val;
        if (legendName === ruleCaption) {
          singleCaption = getIndividualFeatureLabel(feat);
        }
      }
      await setFeatureColor(feat, color, singleCaption);
    } else if (choice === 'individualLabel') {
      const individualLabel =
        normalizeCaption(colorScopeDialog.individualLabel) || normalizeCaption(getIndividualFeatureLabel(feat));
      if (!individualLabel) {
        clearColorScopeDialog();
        return;
      }
      const individualSiblings = findFeaturesWithSameIndividualLabel(feat, individualLabel);
      const allFeatures = [feat, ...individualSiblings];
      const existingLegendEntry = findLegendEntryByCaption(individualLabel);
      let finalCaption = existingLegendEntry?.caption || individualLabel;

      if (existingLegendEntry) {
        updateLegendEntryColorByCaption(existingLegendEntry.caption, color);
      } else {
        const addedCaption = await addLegendEntry(individualLabel, color);
        if (addedCaption && typeof addedCaption === 'string') {
          finalCaption = addedCaption;
          addedLegendCaptions.value.add(addedCaption);
        }
      }

      for (const f of allFeatures) {
        const existingIdx = manualSpecificRules.findIndex(
          (r) => String(r.qual || '').toLowerCase() === 'hash' && r.val === f.svg_id
        );
        if (existingIdx >= 0) {
          manualSpecificRules[existingIdx].color = color;
          manualSpecificRules[existingIdx].cap = finalCaption;
        } else {
          manualSpecificRules.push({
            feat: f.type,
            qual: 'hash',
            val: f.svg_id,
            color: color,
            cap: finalCaption
          });
        }
      }
      applySpecificRulesToSvg();
      await reclaimOrphanedBaseCaptions();
      extractLegendEntries();
    } else if (choice === 'useExisting') {
      if (existingCaptionColor) {
        const targetLegendName = normalizeCaption(legendName) || normalizeCaption(getEffectiveLegendCaption(feat));
        const existingRule = manualSpecificRules.find(
          (r) => String(r.qual || '').toLowerCase() === 'hash' && r.val === feat.svg_id
        );
        if (existingRule) {
          existingRule.color = existingCaptionColor;
        } else {
          manualSpecificRules.push({
            feat: feat.type,
            qual: 'hash',
            val: feat.svg_id,
            color: existingCaptionColor,
            cap: targetLegendName
          });
        }
        if (clickedFeature.value && clickedFeature.value.svg_id === feat.svg_id) {
          clickedFeature.value.color = existingCaptionColor;
          clickedFeature.value.legendName = targetLegendName;
        }
        applyInstantPreview(feat, existingCaptionColor, targetLegendName);
        applySpecificRulesToSvg();
      }
    }

    clearColorScopeDialog();
  };

  const updateClickedFeatureStroke = (strokeColor, strokeWidth) => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const svgId = clickedFeature.value.svg_id;
    const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);

    elements.forEach((el) => {
      if (strokeColor !== null) {
        el.setAttribute('stroke', strokeColor);
        clickedFeature.value.strokeColor = strokeColor;
      }
      if (strokeWidth !== null) {
        const widthVal = parseFloat(strokeWidth);
        if (!isNaN(widthVal)) {
          el.setAttribute('stroke-width', widthVal);
          clickedFeature.value.strokeWidth = widthVal;
        }
      }
    });

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }
  };

  const resetClickedFeatureStroke = () => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const svgId = clickedFeature.value.svg_id;
    const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);

    const originalColor = originalSvgStroke.value.color;
    const originalWidth = originalSvgStroke.value.width;

    elements.forEach((el) => {
      if (originalColor === null) {
        el.removeAttribute('stroke');
      } else {
        el.setAttribute('stroke', originalColor);
      }
      if (originalWidth === null) {
        el.removeAttribute('stroke-width');
      } else {
        el.setAttribute('stroke-width', originalWidth);
      }
    });

    clickedFeature.value.strokeColor = originalColor || '';
    clickedFeature.value.strokeWidth = originalWidth ?? '';

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }
  };

  const resetClickedFeatureFillColor = () => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const feat = clickedFeature.value.feat;
    if (!feat) return;

    const defaultColor = currentColors.value[feat.type];
    if (!defaultColor) {
      console.warn('No default color found for feature type:', feat.type);
      return;
    }

    const caption = getFeatureCaption(feat);

    const siblings = extractedFeatures.value.filter(
      (f) => getFeatureCaption(f) === caption && f.svg_id !== clickedFeature.value.svg_id
    );

    if (siblings.length > 0) {
      resetColorDialog.show = true;
      resetColorDialog.caption = caption;
      resetColorDialog.defaultColor = defaultColor;
      resetColorDialog.siblingCount = siblings.length;
    } else {
      doResetFillColor('this');
    }
  };

  const handleResetColorChoice = async (choice) => {
    resetColorDialog.show = false;
    await doResetFillColor(choice);
  };

  const doResetFillColor = async (choice) => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const feat = clickedFeature.value.feat;
    if (!feat) return;

    const defaultColor = resetColorDialog.defaultColor || currentColors.value[feat.type];
    const caption =
      resetColorDialog.caption ||
      feat.product ||
      feat.gene ||
      feat.locus_tag ||
      feat.note ||
      `${feat.type} at ${feat.start}..${feat.end}`;
    const svgId = clickedFeature.value.svg_id;

    if (choice === 'this' || choice === 'this_with_legend') {
      const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
      elements.forEach((el) => {
        el.setAttribute('fill', defaultColor);
      });

      clickedFeature.value.color = defaultColor;

      const ruleIdx = manualSpecificRules.findIndex((r) => r.qual === 'hash' && r.val === svgId);
      if (ruleIdx !== -1) {
        manualSpecificRules.splice(ruleIdx, 1);
      }

      if (choice === 'this_with_legend') {
        console.log(`Attempting to add legend entry: caption="${caption}", color="${defaultColor}"`);
        const addedCaption = await addLegendEntry(caption, defaultColor);
        console.log(`addLegendEntry returned: ${addedCaption}`);
        if (addedCaption) {
          extractLegendEntries();
          console.log(`Added legend entry: ${addedCaption} with color: ${defaultColor}`);
        } else {
          console.error(`Failed to add legend entry for caption="${caption}"`);
        }
      }

      console.log(`Reset fill color to default (${defaultColor}) for feature: ${svgId}`);
    } else if (choice === 'all') {
      const matchingFeatures = extractedFeatures.value.filter((f) => getFeatureCaption(f) === caption);

      for (const matchFeat of matchingFeatures) {
        const elements = svg.querySelectorAll(`#${CSS.escape(matchFeat.svg_id)}`);
        elements.forEach((el) => {
          el.setAttribute('fill', defaultColor);
        });
      }

      clickedFeature.value.color = defaultColor;

      for (let i = manualSpecificRules.length - 1; i >= 0; i--) {
        const rule = manualSpecificRules[i];
        if (rule.cap === caption) {
          manualSpecificRules.splice(i, 1);
        }
      }

      const legendIdx = legendEntries.value.findIndex((e) => e.caption === caption);
      if (legendIdx !== -1) {
        await removeLegendEntry(caption);
        extractLegendEntries();
      }

      console.log(
        `Reset fill color to default (${defaultColor}) for all ${matchingFeatures.length} features with caption: ${caption}`
      );
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    clickedFeature.value = null;
  };

  const applyStrokeToAllSiblings = async () => {
    if (!clickedFeature.value) return;
    if (!svgContainer.value) return;

    const svg = svgContainer.value.querySelector('svg');
    if (!svg) return;

    const currentStrokeColor = clickedFeature.value.strokeColor;
    const currentStrokeWidth = clickedFeature.value.strokeWidth;
    const currentFillColor = clickedFeature.value.color;
    const clickedSvgId = clickedFeature.value.svg_id;

    const feat = clickedFeature.value.feat;
    if (!feat) {
      console.warn('No feature data found in clickedFeature');
      return;
    }
    const caption = getFeatureCaption(feat);

    let targetLegendEntry = legendEntries.value.find((e) => e.featureIds && e.featureIds.includes(clickedSvgId));

    if (!targetLegendEntry) {
      targetLegendEntry = legendEntries.value.find(
        (e) => e.color && e.color.toLowerCase() === currentFillColor?.toLowerCase()
      );
    }

    let siblingFeatureIds = [];
    if (targetLegendEntry && targetLegendEntry.featureIds && targetLegendEntry.featureIds.length > 0) {
      siblingFeatureIds = targetLegendEntry.featureIds;
      console.log(`Found legend entry "${targetLegendEntry.caption}" with ${siblingFeatureIds.length} features`);
    } else {
      const siblingFeatures = extractedFeatures.value.filter((f) => {
        if (getFeatureCaption(f) !== caption) return false;
        const el = svg.querySelector(`#${CSS.escape(f.svg_id)}`);
        if (!el) return false;
        const fillColor = el.getAttribute('fill');
        return fillColor && fillColor.toLowerCase() === currentFillColor?.toLowerCase();
      });
      siblingFeatureIds = siblingFeatures.map((f) => f.svg_id);
      console.log(
        `Fallback: Found ${siblingFeatureIds.length} features by caption="${caption}" and color="${currentFillColor}"`
      );
    }

    for (const svgId of siblingFeatureIds) {
      const elements = svg.querySelectorAll(`#${CSS.escape(svgId)}`);
      elements.forEach((el) => {
        if (currentStrokeColor) {
          el.setAttribute('stroke', currentStrokeColor);
        }
        if (currentStrokeWidth !== null && currentStrokeWidth !== '') {
          el.setAttribute('stroke-width', parseFloat(currentStrokeWidth));
        }
      });
    }

    console.log(`Applied stroke to ${siblingFeatureIds.length} features`);

    if (targetLegendEntry) {
      const legendGroups = getAllFeatureLegendGroups(svg);
      for (const targetGroup of legendGroups) {
        const entryGroup = targetGroup.querySelector(
          `g[data-legend-key="${CSS.escape(targetLegendEntry.caption)}"]`
        );
        if (entryGroup) {
          const paths = entryGroup.querySelectorAll('path');
          for (const path of paths) {
            const fill = path.getAttribute('fill');
            if (fill && fill !== 'none' && !fill.startsWith('url(')) {
              if (currentStrokeColor) {
                path.setAttribute('stroke', currentStrokeColor);
              }
              if (currentStrokeWidth !== null && currentStrokeWidth !== '') {
                path.setAttribute('stroke-width', parseFloat(currentStrokeWidth));
              }
              console.log(`Updated legend rect stroke for "${targetLegendEntry.caption}"`);
              break;
            }
          }
        }
      }
    }

    const overrideKey = targetLegendEntry ? targetLegendEntry.caption : caption;

    if (!targetLegendEntry) {
      const fillColor = clickedFeature.value.color || '#cccccc';
      const addedCaption = await addLegendEntry(caption, fillColor);
      if (addedCaption) {
        extractLegendEntries();
        targetLegendEntry = legendEntries.value.find((e) => e.caption === addedCaption);
        if (targetLegendEntry) {
          targetLegendEntry.showStroke = true;
          if (!targetLegendEntry.featureIds) targetLegendEntry.featureIds = [];
          if (!targetLegendEntry.featureIds.includes(clickedSvgId)) {
            targetLegendEntry.featureIds.push(clickedSvgId);
          }
        }
        console.log(`Created new legend entry for caption: ${addedCaption}`);
      } else {
        targetLegendEntry = {
          caption: caption,
          originalCaption: caption,
          color: fillColor,
          visible: true,
          showStroke: true,
          featureIds: [clickedSvgId]
        };
        legendEntries.value.push(targetLegendEntry);
      }
    }

    if (!legendStrokeOverrides[overrideKey]) {
      legendStrokeOverrides[overrideKey] = {};
    }
    if (currentStrokeColor) {
      legendStrokeOverrides[overrideKey].strokeColor = currentStrokeColor;
    }
    if (currentStrokeWidth !== null && currentStrokeWidth !== '') {
      legendStrokeOverrides[overrideKey].strokeWidth = parseFloat(currentStrokeWidth);
    }

    skipCaptureBaseConfig.value = true;
    const resultIdx = selectedResultIndex.value;
    if (resultIdx >= 0 && results.value.length > resultIdx) {
      const serializer = new XMLSerializer();
      results.value[resultIdx] = { ...results.value[resultIdx], content: serializer.serializeToString(svg) };
    }

    console.log(
      `Applied stroke (color: ${currentStrokeColor}, width: ${currentStrokeWidth}) to ${siblingFeatureIds.length} features`
    );
  };

  const setFeatureColor = async (feat, color, customCaption = null) => {
    const qualInfo = getFeatureQualifier(feat);
    if (!qualInfo) {
      console.warn(
        `Cannot identify feature: ${feat.type} at ${feat.start}..${feat.end} (no locus_tag, gene, or product)`
      );
      return;
    }
    const { qual, val } = qualInfo;

    const featureKey = feat.id;

    const caption = normalizeCaption(
      customCaption || feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${feat.start}..${feat.end}`
    );
    if (!caption) return;

    const oldOverride = featureColorOverrides[featureKey];
    const oldCaption = normalizeCaption(oldOverride?.caption);
    featureColorOverrides[featureKey] = { color, caption };

    applyInstantPreview(feat, color, caption);

    await nextTick();
    let actualCaption = caption;
    if (pyodideReady.value && caption) {
      if (oldCaption) {
        if (captionsMatch(oldCaption, caption)) {
          const hasNonHashRule = manualSpecificRules.some(
            (rule) => captionsMatch(rule.cap, caption) && String(rule.qual || '').toLowerCase() !== 'hash'
          );
          const hasOtherUses = extractedFeatures.value.some((f) => {
            if (f.svg_id === feat.svg_id) return false;
            return captionsMatch(getEffectiveLegendCaption(f), caption);
          });

          if (hasNonHashRule || hasOtherUses) {
            actualCaption = await addLegendEntry(caption, color);
            if (actualCaption && typeof actualCaption === 'string') {
              addedLegendCaptions.value.add(actualCaption);
            }
          } else {
            updateLegendEntryColorByCaption(oldCaption, color);
          }
        } else {
          removeLegendEntry(oldCaption);
          actualCaption = await addLegendEntry(caption, color);
          if (actualCaption && typeof actualCaption === 'string') {
            addedLegendCaptions.value.add(actualCaption);
          }
        }
      } else {
        actualCaption = await addLegendEntry(caption, color);
        if (actualCaption && typeof actualCaption === 'string') {
          addedLegendCaptions.value.add(actualCaption);
        }
      }

      if (actualCaption && typeof actualCaption === 'string' && actualCaption !== caption) {
        featureColorOverrides[featureKey] = { color, caption: actualCaption };
      }
    }

    skipExtractOnSvgChange.value = true;

    const finalCaption = actualCaption && typeof actualCaption === 'string' ? actualCaption : caption;
    const existingIdx = manualSpecificRules.findIndex(
      (r) => r.feat === feat.type && String(r.qual || '').toLowerCase() === String(qual || '').toLowerCase() && r.val === val
    );

    if (existingIdx >= 0) {
      manualSpecificRules[existingIdx].color = color;
      manualSpecificRules[existingIdx].cap = finalCaption;
    } else {
      manualSpecificRules.push({
        feat: feat.type,
        qual: qual,
        val: val,
        color: color,
        cap: finalCaption
      });
    }

    await reclaimOrphanedBaseCaptions();

    await nextTick();
    await nextTick();
    skipExtractOnSvgChange.value = false;
    extractLegendEntries();
  };

  return {
    applyStrokeToAllSiblings,
    handleColorScopeChoice,
    handleLegendNameCommit,
    requestFeatureColorChange,
    selectLegendNameOption,
    handleResetColorChoice,
    resetClickedFeatureFillColor,
    resetClickedFeatureStroke,
    setFeatureColor,
    updateClickedFeatureColor,
    updateClickedFeatureStroke
  };
};
