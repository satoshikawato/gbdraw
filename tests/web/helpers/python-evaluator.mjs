// One Python process per test file evaluates every snippet the file sends.
// Starting `python` and importing gbdraw costs about a second, so a process
// per evaluation made the Web JavaScript tests spend minutes on start-up.
import { spawn } from 'node:child_process';
import { fileURLToPath } from 'node:url';

const REPO_ROOT = fileURLToPath(new URL('../../../', import.meta.url));
const STDERR_TAIL_BYTES = 16 * 1024;

// Reads one JSON request per line, runs its `code` with `payload` bound, and
// writes `{result}` (the snippet's `result` variable) or `{error}` (the
// traceback) as one JSON line. Snippets print to stderr, never into the replies.
const SERVER = `
import json, sys, traceback
replies = sys.stdout
sys.stdout = sys.stderr
for line in sys.stdin:
    request = json.loads(line)
    scope = {'payload': request['payload']}
    try:
        exec(request['code'], scope)
        reply = {'result': scope.get('result')}
    except BaseException:
        reply = {'error': traceback.format_exc()}
    replies.write(json.dumps(reply) + '\\n')
    replies.flush()
`;

let server = null;
let queue = Promise.resolve();

const setHeld = (child, held) => {
  for (const handle of [child, child.stdin, child.stdout, child.stderr]) {
    if (held) handle.ref();
    else handle.unref();
  }
};

const startServer = () => {
  const child = spawn(process.env.PYTHON || 'python', ['-c', SERVER], {
    cwd: REPO_ROOT,
    stdio: ['pipe', 'pipe', 'pipe']
  });
  const current = { child, buffer: '', stderr: '', pending: null };
  child.stdout.setEncoding('utf8');
  child.stderr.setEncoding('utf8');
  child.stdout.on('data', (chunk) => {
    current.buffer += chunk;
    const end = current.buffer.indexOf('\n');
    if (end === -1 || !current.pending) return;
    const line = current.buffer.slice(0, end);
    current.buffer = current.buffer.slice(end + 1);
    const { resolve } = current.pending;
    current.pending = null;
    resolve(JSON.parse(line));
  });
  child.stderr.on('data', (chunk) => {
    current.stderr = (current.stderr + chunk).slice(-STDERR_TAIL_BYTES);
  });
  const fail = (error) => {
    if (server === current) server = null;
    if (current.pending) {
      const { reject } = current.pending;
      current.pending = null;
      reject(new Error(`Python evaluator stopped: ${error}\n${current.stderr}`));
    }
  };
  child.on('error', (error) => fail(error.message));
  child.on('exit', (code, signal) => fail(`exit ${code ?? signal}`));
  child.stdin.on('error', (error) => fail(error.message));
  process.on('exit', () => child.kill());
  setHeld(child, false);
  return current;
};

const send = (code, payload) => {
  server ??= startServer();
  const current = server;
  return new Promise((resolve, reject) => {
    current.pending = { resolve, reject };
    setHeld(current.child, true);
    current.child.stdin.write(`${JSON.stringify({ code, payload })}\n`);
  }).finally(() => setHeld(current.child, false));
};

/**
 * Runs `code` in the shared Python process with `payload` bound and resolves
 * with the JSON value the code assigns to `result`. One request runs at a time;
 * a Python exception rejects with its traceback.
 * @param {string} code
 * @param {unknown} [payload]
 * @returns {Promise<unknown>}
 */
export const evaluatePython = (code, payload = null) => {
  const reply = queue.then(() => send(code, payload));
  queue = reply.catch(() => {});
  return reply.then((response) => {
    if (response.error !== undefined) throw new Error(`Python evaluation failed:\n${response.error}`);
    return response.result;
  });
};

/** @param {{ features: unknown, rules: unknown, kind: string }} payload */
export const evaluatePythonRules = async (payload) => JSON.parse(await evaluatePython(`
import json
from gbdraw.web_support.rule_matching import evaluate_rules_json
result = evaluate_rules_json(json.dumps(payload['features']), json.dumps(payload['rules']), payload['kind'])
`, payload));

/**
 * The Worker's typed override check, for `setUnmanagedConfigOverrideValidator`.
 * @param {{ mode: string, config: unknown, configOverrides: unknown, managedPaths: unknown, source: string }} payload
 */
export const validateWebConfigOverrides = async (payload) => JSON.parse(await evaluatePython(`
import json
from gbdraw.web_support.config_overrides import validate_web_config_overrides_json
result = validate_web_config_overrides_json(payload['mode'], json.dumps(payload['config']),
    json.dumps(payload['configOverrides']), json.dumps(payload['managedPaths']), payload['source'])
`, payload));

/**
 * The Worker's `synthesizeColorRules` helper (OV-244), for a rule preparation's `synthesize`.
 * @param {unknown} request
 */
export const synthesizePythonColorRules = async (request) => JSON.parse(await evaluatePython(`
import json
from gbdraw.web_support.rule_matching import synthesize_color_rules_json
result = synthesize_color_rules_json(json.dumps(payload))
`, request));

/**
 * The rule owner's `groupColorRules` (app/feature-editor/rule-actions.js) over
 * `preparation`, for a test that stands in for the rule owner.
 * @param {{ synthesizeGroup: Function }} preparation
 * @param {(feature: any) => string[]} [labelsOf]
 */
export const groupColorRulesOver = (preparation, labelsOf = () => []) => (features, color, caption, scope) => async (current) => {
  const { table } = await preparation.synthesizeGroup({ scope, features, labels: features.map(labelsOf), color, caption }, current);
  return table.map((row) => (typeof row === 'number' ? current[row] : row));
};

// Load of a Session 27-44 (DESIGN 4.3): the Worker port as the composition
// root sets it, with Python's upgrade evaluated here; a helper error is
// thrown as the Worker client throws it. `calls` collects each document sent.
export const pythonLegacySessionUpgrader = (calls = []) => async ({ document }) => {
  calls.push(document);
  const reply = JSON.parse(await evaluatePython(`
from gbdraw.web_support.legacy_session import legacy_session_to_drawings_json
result = legacy_session_to_drawings_json(payload)
`, JSON.stringify(document)));
  if (reply.error) throw reply.error;
  return { result: reply };
};
