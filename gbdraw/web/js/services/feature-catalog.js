const FEATURE_CATALOG_RELOAD_MESSAGE =
  'The diagram engine returned incompatible feature metadata. Reload the page and Generate again.';

export const FEATURE_CATALOG_SCHEMA = 3;

const isObject = (value) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
);

const text = (value) => String(value ?? '').trim();

const cloneJson = (value) => {
  if (typeof structuredClone === 'function') return structuredClone(value);
  return JSON.parse(JSON.stringify(value));
};

const catalogError = () => new Error(FEATURE_CATALOG_RELOAD_MESSAGE);

const requireArray = (value) => {
  if (!Array.isArray(value)) throw catalogError();
  return value;
};

export const biologicalFeatureKey = (recordKey, biologicalFeatureId) => {
  const normalizedRecordKey = text(recordKey);
  const normalizedFeatureId = text(biologicalFeatureId);
  return normalizedRecordKey && normalizedFeatureId
    ? `${normalizedRecordKey}\u0000${normalizedFeatureId}`
    : '';
};

export const stableFeatureOverrideKey = (feature) => biologicalFeatureKey(
  feature?.record_key ?? feature?.recordKey,
  feature?.biological_feature_id ?? feature?.biologicalFeatureId
);

const validateReferences = (entries, knownFeatures, fields) => {
  entries.forEach((entry) => {
    if (!isObject(entry)) throw catalogError();
    fields.forEach(([recordField, featureField]) => {
      const recordKey = text(entry[recordField]);
      const featureId = text(entry[featureField]);
      if (!recordKey && !featureId) return;
      if (!knownFeatures.has(biologicalFeatureKey(recordKey, featureId))) {
        throw catalogError();
      }
    });
  });
};

const canonicalSequenceSourceStrings = (sequenceSources) => (
  sequenceSources.map((source) => {
    if (
      !isObject(source)
      || typeof source.sequence !== 'string'
      || !source.sequence
      || /\s/.test(source.sequence)
    ) {
      throw catalogError();
    }
    return source.sequence;
  })
);

const validateCatalogItem = (item, result, resultIndex) => {
  if (
    !isObject(item)
    || item.resultIndex !== resultIndex
    || text(item.resultName) !== text(result?.name)
  ) {
    throw catalogError();
  }

  const recordKeys = requireArray(item.recordKeys).map(text);
  if (recordKeys.some((recordKey) => !recordKey) || new Set(recordKeys).size !== recordKeys.length) {
    throw catalogError();
  }
  const knownRecordKeys = new Set(recordKeys);
  const sequenceSources = item.sequenceSources === undefined
    ? []
    : requireArray(item.sequenceSources);
  const sourceSequences = canonicalSequenceSourceStrings(sequenceSources);
  const biologicalFeatures = requireArray(item.biologicalFeatures);
  const knownFeatures = new Set();
  biologicalFeatures.forEach((feature) => {
    if (!isObject(feature)) throw catalogError();
    const recordKey = text(feature.recordKey);
    const featureId = text(feature.biologicalFeatureId);
    const key = biologicalFeatureKey(recordKey, featureId);
    if (
      !key
      || !knownRecordKeys.has(recordKey)
      || knownFeatures.has(key)
    ) {
      throw catalogError();
    }
    if (feature.translationFromAminoAcidSequence !== undefined) {
      const qualifiers = isObject(feature.qualifiers) ? feature.qualifiers : {};
      const aminoAcidSequence = Object.prototype.hasOwnProperty.call(
        feature,
        'amino_acid_sequence'
      )
        ? feature.amino_acid_sequence
        : feature.aminoAcidSequence;
      if (
        feature.translationFromAminoAcidSequence !== true
        || typeof aminoAcidSequence !== 'string'
        || !aminoAcidSequence.trim()
        || Object.prototype.hasOwnProperty.call(qualifiers, 'translation')
      ) {
        throw catalogError();
      }
    }
    if (feature.sequenceSourceIndex !== undefined) {
      const sourceIndex = feature.sequenceSourceIndex;
      const source = sequenceSources[sourceIndex];
      const expectedRecordIndex = recordKeys.indexOf(recordKey);
      if (
        !Number.isInteger(sourceIndex)
        || sourceIndex < 0
        || Object.prototype.hasOwnProperty.call(feature, 'nucleotide_sequence')
        || Object.prototype.hasOwnProperty.call(feature, 'nucleotideSequence')
        || !isObject(source)
        || !reconstructCatalogNucleotideSequence(
          feature,
          sequenceSources,
          sourceSequences,
          expectedRecordIndex
        )
      ) {
        throw catalogError();
      }
    }
    knownFeatures.add(key);
  });

  const svgIds = new Set();
  const features = requireArray(item.features);
  features.forEach((feature) => {
    if (!isObject(feature)) throw catalogError();
    const svgId = text(feature.svgId);
    const key = biologicalFeatureKey(
      feature.recordKey,
      feature.biologicalFeatureId
    );
    if (!svgId || svgIds.has(svgId) || !knownFeatures.has(key)) {
      throw catalogError();
    }
    svgIds.add(svgId);
  });

  const orthogroups = requireArray(item.orthogroups);
  orthogroups.forEach((group) => {
    if (!isObject(group)) throw catalogError();
    validateReferences(
      requireArray(group.members),
      knownFeatures,
      [['recordKey', 'biologicalFeatureId']]
    );
  });
  const comparisonMatches = requireArray(item.comparisonMatches);
  validateReferences(comparisonMatches, knownFeatures, [
    ['queryRecordKey', 'queryBiologicalFeatureId'],
    ['subjectRecordKey', 'subjectBiologicalFeatureId']
  ]);
  requireArray(item.annotations);
};

export const validateFeatureCatalog = (catalog, results) => {
  const logicalResults = requireArray(results);
  if (!isObject(catalog) || catalog.schema !== FEATURE_CATALOG_SCHEMA) {
    throw catalogError();
  }
  const validated = cloneJson(catalog);
  const items = requireArray(validated.items);
  if (items.length !== logicalResults.length) throw catalogError();
  items.forEach((item, resultIndex) => {
    validateCatalogItem(item, logicalResults[resultIndex], resultIndex);
  });
  return validated;
};

const selectorForFeature = (feature, stableId) => {
  if (isObject(feature.selector)) return cloneJson(feature.selector);
  const selector = {
    type: text(feature.type),
    start: feature.start ?? null,
    end: feature.end ?? null,
    strand: feature.strand ?? null,
    qualifiers: isObject(feature.qualifiers) ? cloneJson(feature.qualifiers) : {},
    hash: stableId
  };
  return selector;
};

const reverseComplementSequence = (sequence) => {
  const complements = {
    A: 'T', C: 'G', G: 'C', T: 'A', U: 'A', R: 'Y', Y: 'R',
    S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', D: 'H', H: 'D',
    V: 'B', N: 'N', '-': '-'
  };
  const normalized = String(sequence || '').replace(/\s+/g, '').toUpperCase();
  return [...normalized]
    .reverse()
    .map((base) => complements[base] || 'N')
    .join('');
};

const reconstructCatalogNucleotideSequence = (
  feature,
  sequenceSources,
  sourceSequences,
  expectedRecordIndex
) => {
  const sourceIndex = feature?.sequenceSourceIndex;
  const source = sequenceSources?.[sourceIndex];
  if (
    !Number.isInteger(sourceIndex)
    || !isObject(source)
    || typeof sourceSequences?.[sourceIndex] !== 'string'
    || !['linear-record', 'circular-reference'].includes(text(source.origin))
    || !Number.isInteger(source.recordIndex)
    || source.recordIndex !== expectedRecordIndex
  ) {
    return '';
  }
  const sourceSequence = sourceSequences[sourceIndex];
  const rawParts = feature?.location_parts ?? feature?.locationParts;
  const parts = Array.isArray(rawParts) && rawParts.length
    ? rawParts
    : [feature];
  const sequences = [];
  for (const part of parts) {
    const start = Number(part?.start);
    const end = Number(part?.end);
    if (
      !Number.isInteger(start)
      || !Number.isInteger(end)
      || start < 0
      || end < start
      || end > sourceSequence.length
    ) {
      return '';
    }
    const sequence = sourceSequence.slice(start, end);
    const strand = part?.strand ?? feature?.strand;
    sequences.push(
      strand === '-' || Number(strand) < 0
        ? reverseComplementSequence(sequence)
        : sequence
    );
  }
  return sequences.join('');
};

const deriveCatalogNucleotideSequence = (
  feature,
  sequenceSources,
  sourceSequences,
  expectedRecordIndex
) => (
  text(feature?.nucleotide_sequence ?? feature?.nucleotideSequence)
  || reconstructCatalogNucleotideSequence(
    feature,
    sequenceSources,
    sourceSequences,
    expectedRecordIndex
  )
);

const firstQualifierValue = (value) => {
  const values = Array.isArray(value) ? value : [value];
  const selected = values.find((entry) => text(entry));
  return selected === undefined ? '' : String(selected);
};

const restoreCompactBiologicalDefaults = (feature) => {
  const restored = feature;
  const qualifiers = isObject(restored.qualifiers)
    ? restored.qualifiers
    : {};
  const hasSnakeAminoAcidSequence = Object.prototype.hasOwnProperty.call(
    restored,
    'amino_acid_sequence'
  );
  const rawAminoAcidSequence = hasSnakeAminoAcidSequence
    ? restored.amino_acid_sequence
    : restored.aminoAcidSequence;
  const aminoAcidSequence = text(rawAminoAcidSequence);
  if (restored.translationFromAminoAcidSequence === true) {
    if (
      typeof rawAminoAcidSequence !== 'string'
      || !rawAminoAcidSequence.trim()
    ) {
      throw catalogError();
    }
    qualifiers.translation = [aminoAcidSequence];
    restored.qualifiers = qualifiers;
    delete restored.translationFromAminoAcidSequence;
  }
  [
    'protein_id',
    'locus_tag',
    'gene_id',
    'old_locus_tag',
    'gene',
    'product'
  ].forEach((field) => {
    if (!text(restored[field])) {
      const value = firstQualifierValue(qualifiers[field]);
      if (value) restored[field] = value;
    }
  });
  if (!text(restored.note)) {
    const note = firstQualifierValue(qualifiers.note);
    if (note) restored.note = Array.from(note).slice(0, 50).join('');
  }
  const parts = restored.location_parts ?? restored.locationParts;
  const start = Number(restored.start);
  const end = Number(restored.end);
  if (
    !Array.isArray(parts)
    && Number.isInteger(start)
    && Number.isInteger(end)
  ) {
    restored.location_parts = [{
      start,
      end,
      strand: restored.strand ?? '',
      display: `${start + 1}..${end}`
    }];
  }
  return restored;
};

const expandBiologicalFeature = (
  feature,
  recordIndex,
  displayRecordId,
  sequenceSources,
  sourceSequences,
  sequenceRecordIndex
) => {
  const recordKey = text(feature.recordKey);
  const biologicalFeatureId = text(feature.biologicalFeatureId);
  const stableId = text(feature.stableFeatureId) || biologicalFeatureId;
  const overrideKey = biologicalFeatureKey(recordKey, biologicalFeatureId);
  const expanded = restoreCompactBiologicalDefaults({
    ...cloneJson(feature),
    id: overrideKey,
    record_key: recordKey,
    recordKey,
    biological_feature_id: biologicalFeatureId,
    biologicalFeatureId,
    stable_override_key: overrideKey,
    stable_feature_id: stableId,
    stable_svg_id: stableId,
    svg_id: stableId,
    fileIdx: recordIndex,
    record_idx: recordIndex,
    displayRecordId,
    selector: selectorForFeature(feature, stableId)
  });
  const nucleotideSequence = deriveCatalogNucleotideSequence(
    expanded,
    sequenceSources,
    sourceSequences,
    sequenceRecordIndex
  );
  if (
    expanded.sequenceSourceIndex !== undefined
    && (
      Object.prototype.hasOwnProperty.call(expanded, 'nucleotide_sequence')
      || Object.prototype.hasOwnProperty.call(expanded, 'nucleotideSequence')
    )
  ) {
    throw catalogError();
  }
  if (expanded.sequenceSourceIndex !== undefined && !nucleotideSequence) {
    throw catalogError();
  }
  if (nucleotideSequence && !text(expanded.nucleotide_sequence)) {
    expanded.nucleotide_sequence = nucleotideSequence;
  }
  return expanded;
};

const renderedFeaturesByBiologicalKey = (features) => {
  const index = new Map();
  features.forEach((feature) => {
    const key = biologicalFeatureKey(
      feature.recordKey,
      feature.biologicalFeatureId
    );
    if (!index.has(key)) index.set(key, []);
    index.get(key).push(feature);
  });
  return index;
};

const expandOrthogroup = (group, biologicalByKey, renderedByKey) => ({
  ...cloneJson(group),
  members: group.members.map((member) => {
    const key = biologicalFeatureKey(
      member.recordKey,
      member.biologicalFeatureId
    );
    const biological = biologicalByKey.get(key) || {};
    const rendered = renderedByKey.get(key)?.[0] || null;
    const recordIndex = Number(biological.record_idx);
    const stableId = text(biological.stable_feature_id);
    return {
      ...cloneJson(member),
      record_key: text(member.recordKey),
      biological_feature_id: text(member.biologicalFeatureId),
      recordKey: text(member.recordKey),
      biologicalFeatureId: text(member.biologicalFeatureId),
      recordIndex: Number.isInteger(recordIndex) ? recordIndex : null,
      record_id: text(biological.record_id),
      featureIndex: Number.isInteger(Number(biological.feature_index))
        ? Number(biological.feature_index)
        : null,
      stableFeatureSvgId: stableId,
      featureSvgId: stableId,
      renderedFeatureSvgId: text(rendered?.svgId),
      proteinId: text(biological.protein_id),
      sourceProteinId: text(
        biological.source_protein_id || biological.protein_id
      ),
      start: biological.start ?? null,
      end: biological.end ?? null,
      strand: biological.strand ?? null,
      gene: text(biological.gene),
      locusTag: text(biological.locus_tag),
      geneId: text(biological.gene_id),
      oldLocusTag: text(biological.old_locus_tag),
      product: text(biological.product),
      note: text(biological.note)
    };
  })
});

export const featureStateFromCatalog = (catalog, { mode = '' } = {}) => {
  if (!isObject(catalog) || catalog.schema !== FEATURE_CATALOG_SCHEMA) {
    throw catalogError();
  }

  const recordKeys = [];
  const recordIndexByKey = new Map();
  catalog.items.forEach((item) => {
    item.recordKeys.forEach((rawRecordKey) => {
      const recordKey = text(rawRecordKey);
      if (recordIndexByKey.has(recordKey)) return;
      recordIndexByKey.set(recordKey, recordKeys.length);
      recordKeys.push(recordKey);
    });
  });
  const publicRecordIdByKey = new Map();
  catalog.items.forEach((item) => {
    item.biologicalFeatures.forEach((feature) => {
      const recordKey = text(feature.recordKey);
      const recordId = text(feature.record_id ?? feature.recordId);
      if (recordKey && recordId && !publicRecordIdByKey.has(recordKey)) {
        publicRecordIdByKey.set(recordKey, recordId);
      }
    });
  });
  const normalizedMode = text(mode).toLowerCase();
  const featureRecordIds = recordKeys.map((recordKey, index) => {
    const recordId = publicRecordIdByKey.get(recordKey) || recordKey;
    return normalizedMode === 'linear'
      ? `File ${index + 1}: ${recordId}`
      : recordId;
  });
  const displayRecordIdByKey = new Map(
    recordKeys.map((recordKey, index) => [recordKey, featureRecordIds[index]])
  );

  const extractedFeatures = [];
  const biologicalFeatures = [];
  const orthogroups = [];
  const annotations = [];
  const comparisonMatches = [];
  const sequenceSources = [];

  catalog.items.forEach((item) => {
    const itemSequenceSources = cloneJson(item.sequenceSources || []);
    const itemSourceSequences = canonicalSequenceSourceStrings(
      itemSequenceSources
    );
    const itemRecordIndexByKey = new Map(
      item.recordKeys.map((recordKey, recordIndex) => [
        text(recordKey),
        recordIndex
      ])
    );
    const expandedBiological = item.biologicalFeatures.map((feature) => (
      expandBiologicalFeature(
        feature,
        recordIndexByKey.get(text(feature.recordKey)) ?? 0,
        displayRecordIdByKey.get(text(feature.recordKey)) || text(feature.recordKey),
        itemSequenceSources,
        itemSourceSequences,
        itemRecordIndexByKey.get(text(feature.recordKey)) ?? -1
      )
    ));
    const biologicalByKey = new Map(
      expandedBiological.map((feature) => [stableFeatureOverrideKey(feature), feature])
    );
    const renderedByKey = renderedFeaturesByBiologicalKey(item.features);
    biologicalFeatures.push(...expandedBiological);
    item.features.forEach((rendered) => {
      const key = biologicalFeatureKey(
        rendered.recordKey,
        rendered.biologicalFeatureId
      );
      const biological = biologicalByKey.get(key);
      if (!biological) throw catalogError();
      extractedFeatures.push({
        ...cloneJson(biological),
        ...cloneJson(rendered),
        id: key,
        record_key: text(rendered.recordKey),
        biological_feature_id: text(rendered.biologicalFeatureId),
        stable_override_key: key,
        stable_feature_id: biological.stable_feature_id,
        stable_svg_id: biological.stable_svg_id,
        rendered_feature_svg_id: text(rendered.svgId),
        svg_id: text(rendered.svgId),
        fill_color: text(rendered.fillColor)
      });
    });
    orthogroups.push(
      ...item.orthogroups.map((group) => (
        expandOrthogroup(group, biologicalByKey, renderedByKey)
      ))
    );
    annotations.push(...cloneJson(item.annotations));
    comparisonMatches.push(...cloneJson(item.comparisonMatches));
    sequenceSources.push(...itemSequenceSources);
  });

  return {
    extractedFeatures,
    biologicalFeatures,
    featureRecordIds,
    featureSelectorSafetyScope: [],
    selectedFeatureRecordIdx: 0,
    orthogroups,
    annotations,
    comparisonMatches,
    sequenceSources
  };
};

export const featureCatalogReloadMessage = FEATURE_CATALOG_RELOAD_MESSAGE;
