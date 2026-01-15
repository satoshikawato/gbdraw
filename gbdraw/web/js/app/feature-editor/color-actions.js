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
    findFeaturesWithSameCaption,
    findMatchingRegexRule,
    getFeatureQualifier
  } = ruleActions;
  const { applyInstantPreview } = featureSvgActions;

  const updateClickedFeatureColor = async (color) => {
    if (!clickedFeature.value) return;
    const feat = clickedFeature.value.feat;
    const customName = clickedFeature.value.legendName?.trim();
    const legendName = customName ? customName : clickedFeature.value.label;

    const matchingRule = findMatchingRegexRule(feat);
    const ruleMatchCount = matchingRule ? countFeaturesMatchingRule(matchingRule) : 0;

    const siblings = findFeaturesWithSameCaption(feat, legendName);
    const siblingCount = siblings.length;

    const existingCaption = findExistingColorForCaption(feat, legendName);

    const individualLabel =
      feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${feat.start}..${feat.end}`;
    let individualLabelSiblingCount = 0;
    if (matchingRule && ruleMatchCount > 1) {
      const ruleCaption = matchingRule.cap || matchingRule.val;
      if (individualLabel !== ruleCaption) {
        const individualSiblings = findFeaturesWithSameCaption(feat, individualLabel);
        individualLabelSiblingCount = individualSiblings.length;
      }
    }

    const needsDialog = matchingRule || siblingCount > 0 || existingCaption;

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
      clickedFeature.value = null;
    } else {
      clickedFeature.value.color = color;
      setFeatureColor(feat, color, legendName);
    }
  };

  const handleColorScopeChoice = async (choice) => {
    const { feat, color, matchingRule, legendName, existingCaptionColor } = colorScopeDialog;

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
      const siblings = findFeaturesWithSameCaption(feat, legendName);
      const allFeatures = [feat, ...siblings];

      for (const f of allFeatures) {
        const existingIdx = manualSpecificRules.findIndex((r) => r.qual === 'hash' && r.val === f.svg_id);
        if (existingIdx >= 0) {
          manualSpecificRules[existingIdx].color = color;
          manualSpecificRules[existingIdx].cap = legendName;
        } else {
          manualSpecificRules.push({
            feat: f.type,
            qual: 'hash',
            val: f.svg_id,
            color: color,
            cap: legendName
          });
        }
      }
      addLegendEntry(legendName, color);
      applySpecificRulesToSvg();
    } else if (choice === 'single') {
      let singleCaption = legendName;
      if (matchingRule && colorScopeDialog.ruleMatchCount > 1) {
        const ruleCaption = matchingRule.cap || matchingRule.val;
        if (legendName === ruleCaption) {
          singleCaption =
            feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${feat.start}..${feat.end}`;
        }
      }
      await setFeatureColor(feat, color, singleCaption);
    } else if (choice === 'individualLabel') {
      const individualLabel = colorScopeDialog.individualLabel;
      const individualSiblings = findFeaturesWithSameCaption(feat, individualLabel);
      const allFeatures = [feat, ...individualSiblings];

      for (const f of allFeatures) {
        const existingIdx = manualSpecificRules.findIndex((r) => r.qual === 'hash' && r.val === f.svg_id);
        if (existingIdx >= 0) {
          manualSpecificRules[existingIdx].color = color;
          manualSpecificRules[existingIdx].cap = individualLabel;
        } else {
          manualSpecificRules.push({
            feat: f.type,
            qual: 'hash',
            val: f.svg_id,
            color: color,
            cap: individualLabel
          });
        }
      }
      addLegendEntry(individualLabel, color);
      applySpecificRulesToSvg();
    } else if (choice === 'useExisting') {
      if (existingCaptionColor) {
        const existingRule = manualSpecificRules.find((r) => r.qual === 'hash' && r.val === feat.svg_id);
        if (existingRule) {
          existingRule.color = existingCaptionColor;
        } else {
          manualSpecificRules.push({
            feat: feat.type,
            qual: 'hash',
            val: feat.svg_id,
            color: existingCaptionColor,
            cap: legendName
          });
        }
        applyInstantPreview(feat, existingCaptionColor, legendName);
        applySpecificRulesToSvg();
      }
    }

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

    const caption =
      customCaption || feat.product || feat.gene || feat.locus_tag || `${feat.type} at ${feat.start}..${feat.end}`;

    const oldOverride = featureColorOverrides[featureKey];
    const oldCaption = oldOverride?.caption;
    featureColorOverrides[featureKey] = { color, caption };

    applyInstantPreview(feat, color, caption);

    await nextTick();
    let actualCaption = caption;
    if (pyodideReady.value && caption) {
      if (oldCaption) {
        if (oldCaption === caption) {
          const hasNonHashRule = manualSpecificRules.some(
            (rule) => rule.cap === caption && String(rule.qual || '').toLowerCase() !== 'hash'
          );
          const hasOtherUses = extractedFeatures.value.some((f) => {
            if (f.svg_id === feat.svg_id) return false;
            const override = featureColorOverrides[f.id];
            const featCaption = override?.caption || f.product || f.gene || f.locus_tag || f.type;
            return featCaption === caption;
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
    const existingIdx = manualSpecificRules.findIndex((r) => r.feat === feat.type && r.qual === qual && r.val === val);

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

    await nextTick();
    await nextTick();
    skipExtractOnSvgChange.value = false;
    extractLegendEntries();
  };

  return {
    applyStrokeToAllSiblings,
    handleColorScopeChoice,
    handleResetColorChoice,
    resetClickedFeatureFillColor,
    resetClickedFeatureStroke,
    setFeatureColor,
    updateClickedFeatureColor,
    updateClickedFeatureStroke
  };
};
