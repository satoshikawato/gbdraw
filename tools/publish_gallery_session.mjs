#!/usr/bin/env node

import { spawn } from 'node:child_process';
import { webcrypto } from 'node:crypto';
import { readFile, writeFile } from 'node:fs/promises';
import { resolve as resolvePath } from 'node:path';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { resolveLinearComparisonPlan } from '../gbdraw/web/js/services/linear-comparisons.js';
import {
  createGallerySessionPublication, galleryPublicationDocument, galleryPublicationView
} from '../gbdraw/web/js/services/gallery-session-publication.js';
import {
  assertCanonicalRenderRequestsEquivalent,
  buildCanonicalRenderRequest,
  buildCanonicalRequestState,
  readPublishedRenderRequest,
  projectCanonicalSessionRequest
} from '../gbdraw/web/js/services/session-request.js';

if (!globalThis.crypto) globalThis.crypto = webcrypto;

const REPO_ROOT = fileURLToPath(new URL('..', import.meta.url));
// The publication's draft split is Python's `split_draft_into_modes`, and its
// Session files are written by Python's `write_serialized_session` (both in
// gbdraw/session_io.py), so a `.gz` Session is compressed as the CLI's.
// `tools/refresh_gallery_sessions.py` names its own interpreter in PYTHON.
const SPLIT_DRAFT = `
import json, sys
from gbdraw.session_io import mode_split_palette_colors, split_draft_into_modes
payload = json.load(sys.stdin)
context = payload['context']
json.dump(split_draft_into_modes(
    payload['draft'],
    committed_mode=context['committedMode'],
    palette_colors=mode_split_palette_colors(payload['draft'].get('config')),
), sys.stdout)
`;
const WRITE_SESSION = `
import sys
from gbdraw.session_io import write_serialized_session
write_serialized_session(sys.argv[1], sys.stdin.buffer.read())
`;

/**
 * Run Python `source` with `args` and `input` on its standard input.
 * @param {string} task what the code does, for the error
 * @param {string} source
 * @param {string[]} args
 * @param {string} input
 * @returns {Promise<string>} its standard output
 */
const runPython = (task, source, args, input) => new Promise((resolve, reject) => {
  const child = spawn(process.env.PYTHON || 'python', ['-c', source, ...args], {
    cwd: REPO_ROOT,
    stdio: ['pipe', 'pipe', 'pipe']
  });
  let stdout = '';
  let stderr = '';
  child.stdout.setEncoding('utf8').on('data', (chunk) => { stdout += chunk; });
  child.stderr.setEncoding('utf8').on('data', (chunk) => { stderr += chunk; });
  child.on('error', reject);
  // A Python that exits before it reads its input (EPIPE) is reported by `close` with its stderr.
  child.stdin.on('error', () => {});
  child.on('close', (code) => {
    if (code === 0) resolve(stdout);
    else reject(new Error(`Python's ${task} failed (exit ${code}):\n${stderr}`));
  });
  child.stdin.end(input);
});

/** @type {import('../gbdraw/web/js/services/gallery-session-publication.js').PublicationDraftSplit} */
const pythonDraftSplit = async (draft, context) => JSON.parse(
  await runPython('draft split', SPLIT_DRAFT, [], JSON.stringify({ draft, context }))
);

const {
  finalizeGallerySessionPublication,
  prepareGallerySessionForPublication,
  validateGalleryPublicationReadiness
} = createGallerySessionPublication({
  assertRequestsEquivalent: assertCanonicalRenderRequestsEquivalent,
  buildRequest: buildCanonicalRenderRequest,
  buildRequestState: buildCanonicalRequestState,
  readRequest: readPublishedRenderRequest,
  projectRequest: projectCanonicalSessionRequest,
  resolveComparisonPlan: resolveLinearComparisonPlan,
  splitDraft: pythonDraftSplit
});

// A Session file as the publication reads it: the one-drawing view of the
// drawing it shows, and that drawing's ID and name (Session 51).
const readSession = async (path) => {
  const payload = await readFile(path);
  const decoded = payload[0] === 0x1f && payload[1] === 0x8b
    ? gunzipSync(payload)
    : payload;
  return galleryPublicationView(JSON.parse(decoded.toString('utf8')));
};

const writeSession = async (path, session) => {
  await runPython('Session write', WRITE_SESSION, [resolvePath(path)], JSON.stringify(session));
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
    const input = await readSession(inputPath);
    const result = await prepareGallerySessionForPublication(input.view);
    await writeSession(outputPath, galleryPublicationDocument(result.session, input.drawing));
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
    const prepared = await readSession(preparedPath);
    const result = await finalizeGallerySessionPublication({
      prepared: prepared.view,
      replayed: (await readSession(replayedPath)).view
    });
    await writeSession(outputPath, galleryPublicationDocument(result.session, prepared.drawing));
    await writeReport(reportPath, safeReport(result.session, result.equivalence));
    return;
  }
  if (command === 'verify') {
    const [inputPath, reportPath] = args;
    if (!inputPath) {
      throw new Error('Usage: publish_gallery_session.mjs verify INPUT [REPORT]');
    }
    const session = (await readSession(inputPath)).view;
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
