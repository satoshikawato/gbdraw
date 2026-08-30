#!/usr/bin/env node

import { webcrypto } from 'node:crypto';
import { readFile, writeFile } from 'node:fs/promises';
import { gunzipSync, gzipSync } from 'node:zlib';
import { resolveLinearComparisonPlan } from '../gbdraw/web/js/app/linear-comparisons.js';
import {
  createGallerySessionPublication
} from '../gbdraw/web/js/services/gallery-session-publication.js';
import { promoteGallerySessionToCurrent } from '../gbdraw/web/js/services/gallery-session-migration.js';
import {
  assertCanonicalRenderRequestsEquivalent,
  buildCanonicalRenderRequest,
  buildCanonicalRequestState,
  promoteCanonicalRenderRequestToCurrent,
  projectCanonicalSessionRequest
} from '../gbdraw/web/js/services/session-request.js';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const {
  finalizeGallerySessionPublication,
  prepareGallerySessionForPublication,
  validateGalleryPublicationReadiness
} = createGallerySessionPublication({
  promoteSession: promoteGallerySessionToCurrent,
  assertRequestsEquivalent: assertCanonicalRenderRequestsEquivalent,
  buildRequest: buildCanonicalRenderRequest,
  buildRequestState: buildCanonicalRequestState,
  promoteRequest: promoteCanonicalRenderRequestToCurrent,
  projectRequest: projectCanonicalSessionRequest,
  resolveComparisonPlan: resolveLinearComparisonPlan
});

const readSession = async (path) => {
  const payload = await readFile(path);
  const decoded = payload[0] === 0x1f && payload[1] === 0x8b
    ? gunzipSync(payload)
    : payload;
  return JSON.parse(decoded.toString('utf8'));
};

const writeSession = async (path, session) => {
  const payload = Buffer.from(JSON.stringify(session), 'utf8');
  await writeFile(
    path,
    path.toLowerCase().endsWith('.gz')
      ? gzipSync(payload, { level: 6, mtime: 0 })
      : payload
  );
};

const safeReport = (session, equivalence) => ({
  version: session.version,
  requestSchema: session.renderRequest?.schema,
  mode: session.renderRequest?.mode,
  grouping: session.renderRequest?.grouping,
  recordCount: session.renderRequest?.records?.length || 0,
  committedDigest: equivalence?.expected?.digest || null,
  rebuiltDigest: equivalence?.actual?.digest || null,
  resourceBindings: equivalence?.expected?.resourceBindings || []
});

const writeReport = async (path, report) => {
  if (!path) return;
  await writeFile(path, `${JSON.stringify(report, null, 2)}\n`, 'utf8');
};

const main = async () => {
  const [command, ...args] = process.argv.slice(2);
  if (command === 'prepare') {
    const [inputPath, outputPath, reportPath] = args;
    if (!inputPath || !outputPath) {
      throw new Error('Usage: publish_gallery_session.mjs prepare INPUT OUTPUT [REPORT]');
    }
    const result = await prepareGallerySessionForPublication(await readSession(inputPath));
    await writeSession(outputPath, result.session);
    await writeReport(reportPath, safeReport(result.session, result.equivalence));
    return;
  }
  if (command === 'finalize') {
    const [preparedPath, replayedPath, outputPath, reportPath] = args;
    if (!preparedPath || !replayedPath || !outputPath) {
      throw new Error(
        'Usage: publish_gallery_session.mjs finalize PREPARED REPLAYED OUTPUT [REPORT]'
      );
    }
    const result = await finalizeGallerySessionPublication({
      prepared: await readSession(preparedPath),
      replayed: await readSession(replayedPath)
    });
    await writeSession(outputPath, result.session);
    await writeReport(reportPath, safeReport(result.session, result.equivalence));
    return;
  }
  if (command === 'verify') {
    const [inputPath, reportPath] = args;
    if (!inputPath) {
      throw new Error('Usage: publish_gallery_session.mjs verify INPUT [REPORT]');
    }
    const session = await readSession(inputPath);
    const result = await validateGalleryPublicationReadiness(session);
    await writeReport(reportPath, safeReport(session, result.equivalence));
    return;
  }
  throw new Error('Expected command: prepare, finalize, or verify.');
};

main().catch((error) => {
  console.error(error?.message || error);
  process.exitCode = 1;
});
