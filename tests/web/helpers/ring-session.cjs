// The synthetic three-ring Session, with its ring rows as a command line search recorded them, so the rule that
// reuses rows by the runtime that searched them (`isReusableLosatEntry`) lets a test replay them without a search.
// The file itself keeps Session 40's rows, which have no runtime record and search again once.
const { readFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');

const RING_SESSION_FILE = resolve(__dirname, '..', '..', 'fixtures', 'sessions', 'synthetic_conservation.gbdraw-session.json.gz');
const CLI_RUNTIME = Object.freeze({ kind: 'losat', version: '0.1.0', source: 'managed', path: 'LOSAT', program: 'blastn' });

const ringSessionText = () => {
  const session = JSON.parse(gunzipSync(readFileSync(RING_SESSION_FILE)).toString('utf8'));
  session.losatCache.entries.forEach((entry) => { entry.runtime = { ...CLI_RUNTIME }; });
  return JSON.stringify(session);
};

// For `setInputFiles`.
const ringSessionUpload = () => ({
  name: 'synthetic_conservation.gbdraw-session.json',
  mimeType: 'application/json',
  buffer: Buffer.from(ringSessionText())
});

module.exports = { RING_SESSION_FILE, ringSessionText, ringSessionUpload };
