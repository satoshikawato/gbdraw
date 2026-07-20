#!/usr/bin/env node

import { readFile, writeFile } from 'node:fs/promises';
import { gunzipSync, gzipSync } from 'node:zlib';
import { promoteGallerySessionToCanonicalV3 } from '../gbdraw/web/js/services/gallery-session-migration.js';

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error(
    'Usage: node --experimental-default-type=module tools/promote_gallery_session.mjs INPUT OUTPUT'
  );
}

const input = await readFile(inputPath);
const decoded = input[0] === 0x1f && input[1] === 0x8b
  ? gunzipSync(input)
  : input;
const session = JSON.parse(decoded.toString('utf8'));
const promoted = promoteGallerySessionToCanonicalV3(session);
const payload = Buffer.from(JSON.stringify(promoted), 'utf8');
const output = outputPath.toLowerCase().endsWith('.gz')
  ? gzipSync(payload, { level: 6, mtime: 0 })
  : payload;
await writeFile(outputPath, output);

