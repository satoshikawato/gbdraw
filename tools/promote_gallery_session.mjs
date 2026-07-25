#!/usr/bin/env node

import fs from 'fs';
import zlib from 'zlib';
import { promoteGallerySessionToCanonicalV3 } from '../gbdraw/web/js/services/gallery-session-migration.js';

const { readFile, writeFile } = fs.promises;
const { gunzipSync, gzipSync } = zlib;

const main = async () => {
  const [inputPath, outputPath] = process.argv.slice(2);
  if (!inputPath || !outputPath) {
    throw new Error(
      'Usage: node tools/promote_gallery_session.mjs INPUT OUTPUT'
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
};

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});

