const { expect } = require('@playwright/test');

const expectBiologicalOwnersUnchanged = (snapshot) => {
  expect(snapshot.structural).toMatchObject({
    directEditFeatureCatalogReplacementCount: 0,
    directEditExtractedFeatureReplacementCount: 0,
    directEditBiologicalFeatureReplacementCount: 0,
    directEditOrthogroupReplacementCount: 0
  });
};

const expectDirectEditFlushed = (before, after) => {
  expect(after.result.sha256).not.toBe(before.result.sha256);
  expect(
    after.structural.directEditResultFlushCount
      - before.structural.directEditResultFlushCount
  ).toBe(1);
  expect(
    after.structural.directEditSvgSerializationCount
      - before.structural.directEditSvgSerializationCount
  ).toBe(1);
  expect(after.history.artifactReplacementHistoryEntryCount).toBe(
    before.history.artifactReplacementHistoryEntryCount
  );
  expect(after.history.artifactCheckpointBuilds).toBe(before.history.artifactCheckpointBuilds);
  expect(after.history.generatedArtifactFullCloneCount).toBe(
    before.history.generatedArtifactFullCloneCount
  );
  expect(after.history.generatedArtifactFullSerializationCount).toBe(
    before.history.generatedArtifactFullSerializationCount
  );
  expectBiologicalOwnersUnchanged(after);
};

module.exports = {
  expectBiologicalOwnersUnchanged,
  expectDirectEditFlushed
};

