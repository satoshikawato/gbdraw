import { execFileSync } from 'node:child_process';

export const evaluatePythonRules = async (payload) => JSON.parse(execFileSync(
  'python', ['-c', `import json,sys
from gbdraw.web_support.rule_matching import evaluate_rules_json
p=json.load(sys.stdin)
print(evaluate_rules_json(json.dumps(p['features']),json.dumps(p['rules']),p['kind']))`],
  { input: JSON.stringify(payload), encoding: 'utf8', maxBuffer: 32 * 1024 * 1024, stdio: ['pipe', 'pipe', 'pipe'] }
));
