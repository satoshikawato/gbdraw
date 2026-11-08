import { execFileSync } from 'node:child_process';

/**
 * Python's reply to the Worker helper `planLinearUnits` for its payload (`{ records, pairs }`): the real
 * `plan_linear_units_json`, so a test never holds a second description of the unit plan.
 * @param {{ records: unknown, pairs: unknown }} payload
 */
export const pythonLinearUnitPlan = (payload) => JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json, sys
from gbdraw.web_support.linear_units import plan_linear_units_json
payload = json.load(sys.stdin)
print(plan_linear_units_json(json.dumps(payload["records"]), json.dumps(payload["pairs"])))
`], { input: JSON.stringify(payload), encoding: 'utf8' }));
